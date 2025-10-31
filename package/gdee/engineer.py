"""
Main entry point for the Gene Discovery and Enzyme Engineering (GDEE) platform.

This module provides the primary interface for setting up and running protein engineering
workflows, including variant generation, structure modeling, quality assessment and 
molecular docking.
"""


from gdee.pipeline import PipelineFactory
from gdee.platform import PlatformFactory
import os
import signal
import socket


__all__ = ["ProteinEngineering"]


class Ligand:
    """
    Represents a ligand molecule for molecular docking calculations.
    
    A ligand contains information about the molecule file and associated
    measurements to be computed during the pipeline execution.

    Attributes:
        name (str): Unique identifier for the ligand
        filename (str): Path to the ligand structure file (PDBQT format)
        measurements (list): List of tuples containing measurement specifications
    """
    
    def __init__(self, name, filename):
        """
        Initialize a new ligand.
        
        Args:
            name (str): Unique name identifier for the ligand
            filename (str): Path to the ligand structure file
        """
        self.name = name
        self.filename = filename
        self.measurements = []

    def add_measurement(self, name, metric, protein_sel, ligand_sel):
        """
        Add a measurement to be computed for this ligand.
        
        Args:
            name (str): Name identifier for this metric
            metric (str): Type of metric to compute (At this moment only "distance" is supported)
            protein_sel (str): MDAnalysis selection string for protein atom
            ligand_sel (str): MDAnalysis selection string for ligand atom
        """
        self.measurements.append((name, metric, protein_sel, ligand_sel))


class ProteinEngineering:
    """
    Main interface for the Gene Discovery and Enzyme Engineering platform.
    
    This class orchestrates the entire GDEE workflow including:
    - Protein variant generation (MSA-based, mutation-based)
    - 3D structure modeling using MODELLER
    - Model quality assessment with VoroMQA and Normalized DOPE
    - Molecular docking with AutoDock Vina/Smina (With Vina/Vinardo scoring, respectively)
    - Distance measurements
    - Results storage in SQLite database
    
    The workflow can be executed on single machines or distributed across
    multiple nodes using MPI.
    """
    
    def __init__(self, protein_name, database):
        """
        Initialize a new protein engineering workflow.
        
        Args:
            protein_name (str): Name identifier for the target protein
            database (str): Path to SQLite database file for storing results
        """
        self.work_dir = os.getcwd()
        self.protein_name = protein_name
        self.db_file = database
        self.pdb = None
        self.ligands = {}
        
        # I/O configuration for result archiving
        self.io = {
            "output": "files",           # Base name for output archives
            "output_format": ".{:06d}",  # Format string for archive numbering
            "output_freq": 1000          # Number of files per archive
        }
        
        # Platform configuration for execution
        self.platform = {
            "name": "simple",    # Execution platform: "simple" or "mpi"
            "local_cpu": 1       # Number of local CPU cores to use
        }
        
        # External program paths
        self.programs = {
            "mgltools": "mgltools",           # MGLTools path for receptor preparation
            "vina": "vina",                   # AutoDock Vina path
            "vinardo": "smina",               # Smina path for Vinardo scoring
            "voromqa": "voronota-voromqa"     # VoroMQA quality assessment tool path
        }
        
        # Variant generation parameters
        self.variant = {
            "name": "mutation",         # Variant builder: "mutation", "exhaustive", or "msa"
            "matrix": "blosum62",       # Substitution matrix for mutations
            "selection": "",            # Residue selection string for mutations
            "fixed": "",                # Residues to keep fixed during optimization
            "conservative": True,       # Use conservative mutation weights
            "max_iterations": 1000,     # Maximum variants to generate
            "combinations": -1,         # Number of simultaneous mutations (-1 = all)
            "msa": "",                  # FASTA file path
            "excluded": {},             # Per-residue amino acid exclusions
            "excluded_all": ""          # Amino acids to exclude globally
        }
        
        # 3D modeling parameters using MODELLER
        self.model = {
            "name": "modeller",      # Modeling method (currently only MODELLER supported)
            "optimize_radius": 0,    # Radius for local optimization (0 = only mutated)
            "num_models": 5,         # Number of models per variant to be passed to the docking step
            "optimize_level": 0      # Optimization level: 0=fast, 1=normal, 2=slow
        }
        
        # Model quality assessment thresholds
        self.model_quality = {
            "norm_dope": -1,    # Normalized DOPE score threshold (lower is better)
            "voromqa": 0.4      # VoroMQA score threshold (higher is better)
        }
        
        # Molecular docking parameters
        self.evaluator = {
            "name": "vina",         # Scoring method: "vina" or "vinardo"
            "exhaustiveness": 50    # Vina exhaustiveness parameter
        }
        
        self._pipeline = None
        self._terminate = False
        signal.signal(signal.SIGUSR1, self.catch_signals)

    def add_ligand(self, name, filename):
        """
        Add a ligand for molecular docking calculations.
        
        Args:
            name (str): Unique identifier for the ligand
            filename (str): Path to ligand structure file (PDBQT format)
            
        Returns:
            Ligand: The created ligand object for adding measurements
            
        Raises:
            RuntimeError: If a ligand with the same name already exists
        """
        if name in self.ligands:
            raise RuntimeError("Ligand '{}' already exists".format(name))

        ligand = Ligand(name, filename)
        self.ligands[name] = ligand
        return ligand

    def run(self):
        """
        Execute the complete protein engineering workflow.
        
        This method orchestrates the entire pipeline by:
        1. Creating a PipelineFactory and configuring all components
        2. Setting up modeling (MODELLER) and quality assessment (VoroMQA/DOPE)
        3. Configuring molecular docking (Vina/Vinardo) for each ligand
        4. Setting up measurement systems for distance calculations
        5. Creating the execution platform (SimplePlatform or MPIPlatform)
        6. Running the workflow and storing results in the database
        
        Raises:
            RuntimeError: If processing is interrupted by a signal
        """
        # Create and configure the processing pipeline
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

        # Create and configure the execution platform
        platform_factory = PlatformFactory()
        platform_factory.parameters = self.platform
        platform_factory.pipeline = self.pipeline
        platform = platform_factory.make()

        # Execute the workflow
        platform.run()

        if self._terminate:
            raise RuntimeError("Processing interrupted by a system signal. Everything should be fine")

    def catch_signals(self, signal, frame):
        """
        Signal handler for graceful termination.
        
        Handles SIGUSR1 signals to allow clean shutdown of the workflow.
        This is particularly important for MPI-based distributed execution.
        
        Args:
            signal: The received signal
            frame: Current stack frame
        """
        self._terminate = True
        self.pipeline.terminate()
        hostname = socket.gethostname()
        pid = os.getpid()
        print("\nProcess {} on {} caught a user termination signal. Finalizing all workers. This may take some time\n".format(pid, hostname))
