"""
"""


from gdee.pipeline import PipelineFactory
from gdee.platform import PlatformFactory
import os
import signal
import socket


__all__ = ["ProteinEngineering"]


class Ligand:
    def __init__(self, name, filename):
        self.name = name
        self.filename = filename
        self.measurements = []

    def add_measurement(self, name, metric, protein_sel, ligand_sel):
        self.measurements.append((name, metric, protein_sel, ligand_sel))


class ProteinEngineering:
    def __init__(self, protein_name, database):
        self.work_dir = os.getcwd()
        self.protein_name = protein_name
        self.db_file = database
        self.pdb = None
        self.ligands = {}
        self.io = {"output": "files", "output_format": ".{:06d}", "output_freq": 1000}
        self.platform = {"name": "simple", "local_cpu": 1}
        self.programs = {"mgltools": "mgltools", "vina": "vina", "vinardo": "smina", "voromqa": "voronota-voromqa"}
        self.variant = {"name": "mutation", "matrix": "blosum62", "selection": "", "fixed": "", "conservative": True, "max_iterations": 1000, "combinations": -1, "msa": "", "excluded": {}, "excluded_all": ""}
        self.model = {"name": "modeller", "optimize_radius": 0, "num_models": 5, "optimize_level": 0}
        self.model_quality = {"norm_dope": -1, "voromqa": 0.4}
        self.evaluator = {"name": "vina", "exhaustiveness": 50, "atm_type": None}
        self._pipeline = None
        self._terminate = False
        signal.signal(signal.SIGUSR1, self.catch_signals)

    def add_ligand(self, name, filename):
        if name in self.ligands:
            raise RuntimeError("Ligand '{}' already exists".format(name))

        ligand = Ligand(name, filename)
        self.ligands[name] = ligand
        return ligand

    def run(self):
        pipeline_factory = PipelineFactory()
        pipeline_factory.protein_name = self.protein_name
        pipeline_factory.programs = self.programs
        pipeline_factory.work_dir = self.work_dir
        pipeline_factory.pdb = self.pdb
        pipeline_factory.ligands = tuple(self.ligands.values())
        pipeline_factory.db_file = self.db_file
        pipeline_factory.io = self.io
        pipeline_factory.variant_parameters = self.variant
        pipeline_factory.model_parameters = self.model
        pipeline_factory.model_quality_parameters = self.model_quality
        pipeline_factory.evaluator_parameters = self.evaluator
        self.pipeline = pipeline_factory.make()

        platform_factory = PlatformFactory()
        platform_factory.parameters = self.platform
        platform_factory.pipeline = self.pipeline
        platform = platform_factory.make()

        platform.run()

        if self._terminate:
            raise RuntimeError("Processing interrupted by a system signal. Everything should be fine")

    def catch_signals(self, signal, frame):
        self._terminate = True
        self.pipeline.terminate()
        hostname = socket.gethostname()
        pid = os.getpid()
        print("\nProcess {} on {} caught a user termination signal. Finalizing all workers. This may take some time\n".format(pid, hostname))
