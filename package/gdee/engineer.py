"""Main interface for protein engineering workflows."""

from gdee.pipeline import PipelineFactory, RescoreFactory
from gdee.platform import PlatformFactory
import os
import signal
import socket


__all__ = ["ProteinEngineering", "RescoreVariants"]


class Ligand:
    """Represents a ligand for molecular docking.

    Attributes:
        name: Ligand identifier
        filename: Path to PDBQT ligand file
        measurements: List of measurement specifications
    """

    def __init__(self, name, filename):
        """Initialize a ligand.

        Args:
            name: Unique ligand identifier
            filename: Path to PDBQT format file
        """
        self.name = name
        self.filename = filename
        self.measurements = []

    def add_measurement(self, name, metric, protein_sel, ligand_sel):
        """Add a measurement to compute during docking.

        Args:
            name: Measurement identifier
            metric: Metric type (e.g., "distance")
            protein_sel: MDAnalysis selection string for protein atoms
            ligand_sel: MDAnalysis selection string for ligand atoms
        """
        self.measurements.append((name, metric, protein_sel, ligand_sel))


class ProteinEngineering:
    """Main interface for the GDEE protein engineering platform.

    Orchestrates variant generation, 3D modeling, quality assessment,
    molecular docking, and result storage in SQLite database.
    """

    def __init__(self, protein_name, database):
        """Initialize protein engineering workflow.

        Args:
            protein_name: Name identifier for target protein
            database: Path to SQLite database file
        """
        self.work_dir = os.getcwd()
        self.protein_name = protein_name
        self.db_file = database
        self.pdb = None
        self.ligands = {}
        self.io = {"output": "files", "output_format": ".{:06d}", "output_freq": 1000}
        self.platform = {"name": "simple", "local_cpu": 1}
        self.programs = {
            "mgltools": "mgltools",
            "vina": "vina",
            "vinardo": "smina",
            "voromqa": "voronota-voromqa"
        }
        self.variant = {
            "name": "mutation",
            "matrix": "blosum62",
            "selection": "",
            "fixed": "",
            "conservative": True,
            "max_iterations": 1000,
            "combinations": -1,
            "msa": "",
            "excluded": {},
            "excluded_all": ""
        }
        self.model = {
            "name": "modeller",
            "optimize_radius": 0,
            "num_models": 5,
            "optimize_level": 0
        }
        self.model_quality = {
            "norm_dope": -1,
            "voromqa": 0.4
        }
        self.evaluator = {
            "name": "vina",
            "exhaustiveness": 50,
            "atom_type": None
        }
        self._pipeline = None
        self._terminate = False
        signal.signal(signal.SIGUSR1, self.catch_signals)

    def add_ligand(self, name, filename):
        """Add a ligand for docking calculations.

        Args:
            name: Unique ligand identifier
            filename: Path to PDBQT format ligand file

        Returns:
            Ligand: Created ligand object for adding measurements

        Raises:
            RuntimeError: If ligand name already exists
        """
        if name in self.ligands:
            raise RuntimeError("Ligand '{}' already exists".format(name))

        ligand = Ligand(name, filename)
        self.ligands[name] = ligand
        return ligand

    def run(self):
        """Execute the complete protein engineering workflow.

        Runs variant generation, 3D modeling, quality assessment,
        docking, and saves results to database.

        Raises:
            RuntimeError: If processing fails or is interrupted
        """
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
        self._pipeline = pipeline_factory.make()

        platform_factory = PlatformFactory()
        platform_factory.parameters = self.platform
        platform_factory.pipeline = self._pipeline
        platform = platform_factory.make()

        platform.run()

        if self._terminate:
            raise RuntimeError("Processing interrupted by a system signal. Everything should be fine")

    def catch_signals(self, signal, frame):
        """Handle termination signals gracefully.

        Args:
            signal: Signal number
            frame: Stack frame
        """
        self._terminate = True
        self._pipeline.terminate()
        hostname = socket.gethostname()
        pid = os.getpid()
        print("\nProcess {} on {} caught a user termination signal. "
              "Finalizing all workers. This may take some time\n".format(pid, hostname))


class RescoreVariants:
    """Re-evaluate existing docking results with trained metamodel."""

    def __init__(self, in_db, out_db, table, functions, files_path):
        """Initialize rescoring workflow.

        Args:
            in_db: Input database path with docking results
            out_db: Output database path for rescored results
            table: Table name for rescored poses
            functions: Path to pickle file with scoring functions
            files_path: Path to structure files directory
        """
        self.work_dir = os.getcwd()
        self.db_file = in_db
        self.out_db = out_db
        self.table = table
        self.functions = functions
        self.files_path = files_path
        self.platform = {"name": "simple", "local_cpu": 1}
        self._pipeline = None
        self._terminate = False
        signal.signal(signal.SIGUSR1, self.catch_signals)

    def run(self):
        """Execute the rescoring workflow.

        Re-scores poses from input database and saves to output database.

        Raises:
            RuntimeError: If processing fails or is interrupted
        """
        rescore_factory = RescoreFactory()
        rescore_factory.table = self.table
        rescore_factory.pickle_path = self.functions
        rescore_factory.input_db = self.db_file
        rescore_factory.output_db = self.out_db
        rescore_factory.files_path = self.files_path
        self._pipeline = rescore_factory.make()

        platform_factory = PlatformFactory()
        platform_factory.parameters = self.platform
        platform_factory.pipeline = self._pipeline
        platform = platform_factory.make()

        platform.run()

        if self._terminate:
            raise RuntimeError("Processing interrupted by a system signal. Everything should be fine")

    def catch_signals(self, signal, frame):
        """Handle termination signals gracefully.

        Args:
            signal: Signal number
            frame: Stack frame
        """
        self._terminate = True
        self._pipeline.terminate()
        hostname = socket.gethostname()
        pid = os.getpid()
        print("\nProcess {} on {} caught a user termination signal. "
              "Finalizing all workers. This may take some time\n".format(pid, hostname))
