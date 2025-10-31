"""
"""

import sys
import numpy as np
from path import Path
import MDAnalysis as mda
import subprocess
from tempfile import TemporaryDirectory
from .pdbqt import PDBQT
from gdee.misc import DataContainer
import warnings
from oddt.virtualscreening import virtualscreening as vs
import math

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


def external_command(arguments, name):
    proc = subprocess.run(
        arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if proc.returncode:
        raise RuntimeError("Error processing job '{}':\n{}\n".format(name, proc.stderr.decode("UTF-8")))


class RescoringDocking:
    def __init__(self, parameters):
        self.parameters = parameters
        self.ligand = parameters["ligand"]
        self.prepare_receptor = Path(parameters["mgltools"]) / "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        self.name = "vina"
        self.program = parameters["vina"]
        self.n_cpu = -1
        self.pickle_path = parameters["function"] # It can be a path to a pickle file or the name of the scoring function to train
        self.ligands_type = "pdbqt"


    def run(self, job_data):
        job_dir = job_data.job_dir
        temp_dir = TemporaryDirectory(prefix="gdee_docking")
        temp_path = Path(temp_dir.name)
        Path(self.ligand.filename).copy(temp_path / "ligand.pdbqt")

        with temp_path:
            for idx, model in enumerate(job_data.modeling.models):
                if "evals" not in model:
                    model.evals = {}

                if model.rejected:
                    continue

                protein = mda.Universe(str(job_dir / model.pdb))
                pos = protein.atoms.positions
                size = np.array(self.parameters["box_size"], np.float32) + 1
                center = np.array(self.parameters["box_center"], np.float32)
                upper = center + size
                lower = center - size
                atoms = np.all((pos <= upper) & (pos >= lower), axis=1)
                smaller = protein.atoms[atoms].residues.atoms
                if not len(smaller):
                    raise RuntimeError("No protein atoms inside the docking search box")
                smaller.write(str(temp_path / "model.pdb"))

                try:
                    self.run_docking(job_data)
                    pdbqt = PDBQT("results.pdbqt")

                except Exception as error:
                    print(error)

                else:
                    # Sanity check
                    if not pdbqt.size():
                        continue

                    # Process and save results
                    results_pdb = "docking_{}_{:04d}.pdb".format(self.ligand.name, idx)
                    pdbqt.write_pdb(job_dir / results_pdb)

                    docking = DataContainer()
                    docking.ligand_name = self.ligand.name
                    docking.ligand_file = self.ligand.filename
                    docking.method = self.name
                    docking.pdb = results_pdb
                    docking.energies = [model.energy for model in pdbqt]

                    model.evals[self.ligand.name] = [docking]

                    # run rescoring

                    results, rescoring_method = self.run_rescoring("results.pdbqt", str(job_dir / model.pdb), self.pickle_path)

                    rescoring = DataContainer()
                    rescoring.ligand_name = self.ligand.name
                    rescoring.ligand_file = self.ligand.filename
                    rescoring.method = rescoring_method
                    rescoring.pdb = results_pdb
                    rescoring.energies = [self.kd_to_energy(float(results[index])) for index in range(0,len(results))]
                       

                    # Adding the results to job_data                    
                    model.evals[self.ligand.name].append(rescoring)

                    
        return job_data

    def run_docking(self, job_data):
        # Generate model's PDBQT
        command = [
            self.prepare_receptor,
            "-r", "model.pdb",
            "-o", "model.pdbqt",
            "-A", "checkhydrogens",
        ]

        external_command(command, job_data.variant.name)

        # Run docking
        box_center = "--center_x {:.2f} --center_y {:.2f} --center_z {:.2f}".format(*self.parameters["box_center"])
        box_size = "--size_x {:.2f} --size_y {:.2f} --size_z {:.2f}".format(*self.parameters["box_size"])

        command = [
            self.program,
            "--exhaustiveness", self.parameters["exhaustiveness"],
            "--cpu", "1",
            "--num_modes", "500",   # Exaggerated
            "--energy_range", "30", # numbers
            "--receptor", "model.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--out", "results.pdbqt"
        ] + box_center.split(" ") + box_size.split(" ")
        command = list(map(str, command))

        external_command(command, job_data.variant.name)


    
    def run_rescoring(self, ligand, protein, pickle_path):
        vs_rescore = vs(n_cpu=self.n_cpu)
        vs_rescore.load_ligands(self.ligands_type, ligand)
        vs_rescore.score(function=pickle_path, protein=protein)

        # Save results
        results =  []
        rescoring_method = ''
        for mol in vs_rescore.fetch():
            data = mol.data.to_dict()
            
            if len(data) > 0:
                data['name'] = mol.title
            else:
                print('There is no data', file=sys.stderr)
                return False

            for key in data:
                if key.startswith('rf'):
                    results.append(data[key])
                    rescoring_method = key
                elif key.startswith('nn'):
                    results.append(data[key])
                    rescoring_method = key
                elif key.startswith('PLEC'):
                    results.append(data[key])
                    rescoring_method = key

        return results, rescoring_method




    def kd_to_energy(self, pkd):
        R = 8.314
        T = 298.15 
        kd = pow(10, -(float(pkd)))
        j = R * T * (math.log(float(kd)))
        kj = j / 1000
        kcal = kj / 4.18
        return kcal

    


