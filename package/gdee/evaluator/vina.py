"""
AutoDock Vina and Vinardo molecular docking implementations for GDEE platform.

This module provides molecular docking functionality using AutoDock Vina and
Vinardo engines for protein-ligand interaction evaluation and pose generation.
"""


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
    """
    Execute external docking programs with error handling.
    
    Runs AutoDock Vina/Vinardo or MGLTools programs and captures
    output for error reporting and debugging.
    
    Args:
        arguments (list): Command-line arguments for external program
        name (str): Job name for error reporting
        
    Raises:
        RuntimeError: If external program fails or returns non-zero exit code
    """
    proc = subprocess.run(
        arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if proc.returncode:
        raise RuntimeError("Error processing job '{}':\n{}\n".format(
            name, proc.stderr.decode("UTF-8")))


class BaseVina:
    """
    Base class for AutoDock Vina-based molecular docking engines.
    
    This class provides common functionality for Vina and Vinardo docking
    including protein preparation, search box validation, and result processing.
    It handles the complete docking workflow from structure preparation to
    pose generation and energy evaluation.
    
    Attributes:
        parameters (dict): Docking configuration including box geometry,
                          exhaustiveness, and program paths
        name (str): Docking engine identifier
        program (str): Path to docking executable
        ligand (Ligand): Target ligand with structure file
        extra_arguments (list): Engine-specific command-line arguments
        prepare_receptor (Path): Path to MGLTools receptor preparation script
    """
    
    def __init__(self, parameters):
        """
        Initialize base docking engine.
        
        Args:
            parameters (dict): Docking configuration including program paths,
                             search box parameters, and ligand specification
        """
        self.parameters = parameters
        self.name = ""
        self.program = ""
        self.ligand = parameters["ligand"]
        self.extra_arguments = []
        
        # Configure MGLTools receptor preparation script
        self.prepare_receptor = Path(parameters["mgltools"]) / \
            "MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"

    def run(self, job_data):
        """
        Execute molecular docking for all protein models.
        
        This method:
        1. Creates temporary workspace with ligand structure
        2. Processes each non-rejected protein model
        3. Validates protein atoms within search box
        4. Runs docking calculation for each model
        5. Converts PDBQT results to PDB format
        6. Stores docking results in job data
        
        Args:
            job_data (DataContainer): Job data with modeling results
            
        Returns:
            DataContainer: Updated job data with docking evaluations
        """
        job_dir = job_data.job_dir
        temp_dir = TemporaryDirectory(prefix="gdee_docking")
        temp_path = Path(temp_dir.name)
        
        # Copy ligand structure to temporary workspace
        Path(self.ligand.filename).copy(temp_path / "ligand.pdbqt")

        with temp_path:
            for idx, model in enumerate(job_data.modeling.models):
                # Initialize evaluation storage
                if "evals" not in model:
                    model.evals = {}

                # Skip models rejected by quality assessment
                if model.rejected:
                    continue

                # Validate search box contains protein atoms
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
                
                # Write reduced protein structure for docking
                smaller.write(str(temp_path / "model.pdb"))

                try:
                    # Execute docking calculation
                    self.run_docking(job_data)
                    pdbqt = PDBQT("results.pdbqt")

                except Exception as error:
                    print(error)

                else:
                    # Process and save docking results
                    if pdbqt.size():
                        results_pdb = "docking_{}_{:04d}.pdb".format(
                            self.ligand.name, idx)
                        pdbqt.write_pdb(job_dir / results_pdb)

                        # Create docking evaluation data container
                        docking = DataContainer()
                        docking.ligand_name = self.ligand.name
                        docking.ligand_file = self.ligand.filename
                        docking.method = self.name
                        docking.pdb = results_pdb
                        docking.energies = [model.energy for model in pdbqt]

                        model.evals[self.ligand.name] = docking

        return job_data

    def run_docking(self, job_data):
        """
        Execute the docking calculation workflow.
        
        This method:
        1. Converts protein PDB to PDBQT format using MGLTools
        2. Configures docking parameters and search box
        3. Runs AutoDock Vina/Vinardo docking calculation
        4. Generates multiple poses with binding energies
        
        Args:
            job_data (DataContainer): Job data with variant information
            
        Raises:
            RuntimeError: If receptor preparation or docking fails
        """
        # Convert protein PDB to PDBQT format
        command = [
            self.prepare_receptor,
            "-r", "model.pdb",      # Input protein structure
            "-o", "model.pdbqt",    # Output PDBQT format
            "-A", "checkhydrogens", # Add hydrogens if missing
        ]

        external_command(command, job_data.variant.name)

        # Configure docking search box
        box_center = "--center_x {:.2f} --center_y {:.2f} --center_z {:.2f}".format(
            *self.parameters["box_center"])
        box_size = "--size_x {:.2f} --size_y {:.2f} --size_z {:.2f}".format(
            *self.parameters["box_size"])

        # Build docking command
        command = [
            self.program,
            "--exhaustiveness", self.parameters["exhaustiveness"],  # Search thoroughness
            "--cpu", "1",           # Single-threaded execution
            "--num_modes", "500",   # Maximum number of poses (will be filtered)
            "--energy_range", "30", # Energy range for pose selection
            "--receptor", "model.pdbqt",
            "--ligand", "ligand.pdbqt",
            "--out", "results.pdbqt"
        ] + box_center.split(" ") + box_size.split(" ")
        
        # Add engine-specific arguments
        command += self.extra_arguments
        command = list(map(str, command))  # Ensure all arguments are strings

        external_command(command, job_data.variant.name)


class VinaDocking(BaseVina):
    """
    AutoDock Vina molecular docking implementation.
    
    This class implements the standard AutoDock Vina docking engine
    with default scoring function and search algorithms. Vina is widely
    used for virtual screening and structure-based drug design.
    
    Attributes:
        name (str): Engine identifier ("vina")
        program (str): Path to vina executable
    """
    
    def __init__(self, parameters, *args, **kwargs):
        """
        Initialize AutoDock Vina docking engine.
        
        Args:
            parameters (dict): Configuration including vina executable path
            *args, **kwargs: Additional arguments passed to BaseVina
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vina"
        self.program = parameters["vina"]


class VinardoDocking(BaseVina):
    """
    Vinardo molecular docking implementation.
    
    This class implements the Vinardo docking engine, which uses
    an improved scoring function compared to standard Vina. Vinardo
    often provides better accuracy for certain protein-ligand systems.
    
    Attributes:
        name (str): Engine identifier ("vinardo")
        program (str): Path to vinardo executable
        extra_arguments (list): Vinardo-specific command-line flags
    """
    
    def __init__(self, parameters, *args, **kwargs):
        """
        Initialize Vinardo docking engine.
        
        Args:
            parameters (dict): Configuration including vinardo executable path
            *args, **kwargs: Additional arguments passed to BaseVina
        """
        super().__init__(parameters, *args, **kwargs)
        self.name = "vinardo"
        self.program = parameters["vinardo"]
        self.extra_arguments = ["--scoring", "vinardo"]  # Use Vinardo scoring function
