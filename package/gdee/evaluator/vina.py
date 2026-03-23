"""AutoDock Vina/Vinardo molecular docking implementation."""

import os
import numpy as np
from path import Path
import MDAnalysis as mda
import subprocess
from tempfile import TemporaryDirectory
from .pdbqt import PDBQT
from gdee.misc import DataContainer
import warnings

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


def external_command(arguments, name):
    """Execute external command and check for errors.

    Args:
        arguments: List of command arguments
        name: Job name for error reporting

    Returns:
        subprocess.CompletedProcess: Command result

    Raises:
        RuntimeError: If command fails
    """
    proc = subprocess.run(
        arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if proc.returncode:
        raise RuntimeError("Error processing job '{}':\n{}\n".format(name, proc.stderr.decode("UTF-8")))


class BaseVina:
    """Base class for Vina-based docking engines.

    Handles protein preparation, search box validation, and result processing.
    """
    def __init__(self, parameters):
        """Initialize Vina docking base.

        Args:
            parameters: Docking configuration dictionary
        """
        self.parameters = parameters
        self.name = ""
        self.program = ""
        self.ligand = parameters["ligand"]
        self.extra_arguments = []
        self.prepare_receptor = Path(parameters["mgltools"]) / "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
        self.atom_type = parameters['atom_type']

    def run(self, job_data):
        """Execute docking for all models in job.

        Args:
            job_data: Job data with models to dock

        Returns:
            DataContainer: Job data with docking results
        """
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
                    # Process and save results
                    if pdbqt.size():
                        results_pdb = "docking_{}_{:04d}.pdb".format(self.ligand.name,
                                                                     idx)
                        pdbqt.write_pdb(job_dir / results_pdb)

                        docking = DataContainer()
                        docking.ligand_name = self.ligand.name
                        docking.ligand_file = self.ligand.filename
                        docking.method = self.name
                        docking.pdb = results_pdb
                        docking.energies = [model.energy for model in pdbqt]

                        model.evals[self.ligand.name] = docking

        return job_data

    def run_docking(self, job_data):
        """Prepare receptor PDBQT from PDB.
           Execute docking command.

        Args:
            job_data: Job data
        """
        # Generate model's PDBQT
        command = [
            self.prepare_receptor,
            "-r", "model.pdb",
            "-o", "model.pdbqt",
            "-A", "checkhydrogens",
            "-U", "nphs_lps_waters"
        ]

        external_command(command, job_data.variant.name)

        # Run atom type patch
        if self.atom_type:
            self.atom_type_patch('model.pdbqt')


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
        command += self.extra_arguments
        command = list(map(str, command))

        external_command(command, job_data.variant.name)


    def atom_type_patch(self, pdbqt):
        """
        Patch PDBQT file to update atom types based on provided mapping.
        """
        new_pdbqt = 'model_patched.pdbqt'

        with open(pdbqt, 'r') as f, open(new_pdbqt, 'w') as nf:
            for line in f:
                write_line = True
                if line.startswith('ATOM'):
                    # "{name}:{chain}:{id}"
                    atom = "{}:{}:{}".format(line[12:16].strip(), line[21:22], line[22:26].strip())
                    if atom in self.atom_type:
                        if self.atom_type[atom] == "r":
                            write_line = False
                        else:
                            line = line[:77] + "{}\n".format(self.atom_type[atom])
                if write_line:
                    nf.write(line)

        os.remove(pdbqt)
        os.rename(new_pdbqt, pdbqt)


class VinaDocking(BaseVina):
    """AutoDock Vina docking implementation."""
    def __init__(self, parameters, *args, **kwargs):
        """Initialize Vina docking.

        Args:
            parameters: Docking configuration
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vina"
        self.program = parameters["vina"]


class VinardoDocking(BaseVina):
    """Vinardo docking implementation.

    Uses Vinardo scoring function via Smina.
    """
    def __init__(self, parameters, *args, **kwargs):
        """Initialize Vinardo docking.

        Args:
            parameters: Docking configuration
            *args: Additional positional arguments
            **kwargs: Additional keyword arguments
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vinardo"
        self.program = parameters["vinardo"]
        self.extra_arguments = ["--scoring", "vinardo"]
