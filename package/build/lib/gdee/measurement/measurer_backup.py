"""
"""


from path import Path
import MDAnalysis as mda
from gdee.misc import DataContainer
import warnings

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


class Task:
    def __init__(self, metric, name, prot_text, lig_text):
        self.name = name
        self.identifier = "{}|{}|{}".format(metric.name(), prot_text, lig_text)
        self.metric = metric
        self.prot_text = prot_text
        self.lig_text = lig_text
        self.data = []
        self.enabled = False

    def clear(self):
        self.data = []

    def set_system(self, protein, ligand):
        self.prot_sel = protein.select_atoms(self.prot_text)
        self.lig_sel = ligand.select_atoms(self.lig_text)
        self.enabled = bool(len(self.prot_sel) and len(self.lig_sel))

    def compute(self):
        self.data.append(self.metric.compute(self.prot_sel.positions, self.lig_sel.positions))

    def to_container(self):
        data = DataContainer()
        data.name = self.name
        data.identifier = self.identifier
        data.data = self.data
        return data


class Measurer:
    def __init__(self, ligand_name):
        self.ligand_name = ligand_name
        self.task_list = []
        # Avoid repeated measurements
        self.task_names = set()

    def add(self, metric, name, prot_sel, lig_sel):
        if name in self.task_names:
            raise RuntimeError("Repeat of measurement '{}'".format(name))

        self.task_names.add(name)
        self.task_list.append(Task(metric, name, prot_sel, lig_sel))

    def run(self, job_data):
        job_dir = Path(job_data.job_dir)
        modeling = job_data.modeling
        lig_name = self.ligand_name

        # Check if there is at least one evaluation
        ligand_pdb = None
        protein_pdb = None
        for model in modeling.models:
            if lig_name in model.evals:
                ligand_pdb = model.evals[lig_name].pdb
                protein_pdb = model.pdb
                break

        if ligand_pdb is None or protein_pdb is None:
            return job_data

        protein = mda.Universe(str(job_dir / protein_pdb))
        ligand = mda.Universe(str(job_dir / ligand_pdb))

        for task in self.task_list:
            task.set_system(protein, ligand)
            if not task.enabled:
                print("Could not apply measurement '{}' to variant '{}'".format(
                          task.name, job_data.variant.name))

        for model in job_data.modeling.models:
            if lig_name not in model.evals:
                continue

            evaluation = model.evals[lig_name]
            protein.load_new(str(job_dir / model.pdb))
            ligand.load_new(str(job_dir / evaluation.pdb))

            for ts in ligand.trajectory:
                for task in self.task_list:
                    if task.enabled:
                        task.compute()

            if "measurements" not in evaluation:
                evaluation.measurements = []

            for task in self.task_list:
                if task.enabled:
                    evaluation.measurements.append(task.to_container())
                task.clear()

        return job_data
