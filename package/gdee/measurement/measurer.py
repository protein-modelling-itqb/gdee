"""
Molecular measurement orchestration for the GDEE platform.

This module provides the main measurement system that coordinates distance
calculations and geometric analysis between protein-ligand complexes from
molecular docking results.
"""

from path import Path
import MDAnalysis as mda
from gdee.misc import DataContainer
import warnings

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


class Task:
    """
    Individual measurement task for protein-ligand analysis.
    
    This class encapsulates a single measurement calculation including
    atom selections, metric computation, and result storage. Each task
    corresponds to one measurement specification from a ligand configuration.
    
    Attributes:
        name (str): Human-readable measurement name
        identifier (str): Database identifier with metric and selection info
        metric: Measurement calculation instance (e.g., EuclideanDistance)
        prot_text (str): MDAnalysis selection string for protein atoms
        lig_text (str): MDAnalysis selection string for ligand atoms
        data (list): Accumulated measurement values across poses
        enabled (bool): Whether selections are valid for computation
        prot_sel: MDAnalysis atom selection for protein
        lig_sel: MDAnalysis atom selection for ligand
    """
    
    def __init__(self, metric, name, prot_text, lig_text):
        """
        Initialize measurement task.
        
        Args:
            metric: Measurement calculation instance
            name (str): Human-readable measurement name
            prot_text (str): MDAnalysis selection for protein atoms
            lig_text (str): MDAnalysis selection for ligand atoms
        """
        self.name = name
        self.identifier = "{}|{}|{}".format(metric.name(), prot_text, lig_text)
        self.metric = metric
        self.prot_text = prot_text
        self.lig_text = lig_text
        self.data = []
        self.enabled = False
        self.prot_sel = None
        self.lig_sel = None

    def clear(self):
        """Clear accumulated measurement data for next model."""
        self.data = []

    def set_system(self, protein, ligand):
        """
        Configure atom selections for measurement calculation.
        
        Validates that both protein and ligand selections return atoms.
        If either selection is empty, the task is disabled for this variant.
        
        Args:
            protein (mda.Universe): Protein structure universe
            ligand (mda.Universe): Ligand docking poses universe
        """
        try:
            self.prot_sel = protein.select_atoms(self.prot_text)
            self.lig_sel = ligand.select_atoms(self.lig_text)
            self.enabled = bool(len(self.prot_sel) and len(self.lig_sel))
        except Exception:
            # Invalid selection strings disable the task
            self.enabled = False
            self.prot_sel = None
            self.lig_sel = None

    def compute(self):
        """
        Calculate measurement for current protein-ligand configuration.
        
        Computes the metric value using current atomic coordinates
        and appends the result to the data list.
        """
        if self.enabled and self.prot_sel is not None and self.lig_sel is not None:
            value = self.metric.compute(self.prot_sel.positions, self.lig_sel.positions)
            self.data.append(value)

    def to_container(self):
        """
        Package measurement results for database storage.
        
        Creates a data container with measurement metadata and results
        suitable for storage in the database.
        
        Returns:
            DataContainer: Container with name, identifier, and measurement data
        """
        container = DataContainer()
        container.name = self.name
        container.identifier = self.identifier
        container.data = self.data.copy()
        return container


class Measurer:
    """
    Orchestrates molecular measurements for protein-ligand complexes.
    
    This class coordinates multiple measurement tasks across all docking
    poses and models for a specific ligand. It handles structure loading,
    atom selection validation, and measurement computation for pose analysis.
    
    Attributes:
        ligand_name (str): Name of target ligand
        task_list (list): List of measurement tasks to execute
        task_names (set): Cache of task names to prevent duplicates
    """
    
    def __init__(self, ligand_name):
        """
        Initialize measurement orchestrator for a specific ligand.
        
        Args:
            ligand_name (str): Name of the ligand to measure
        """
        self.ligand_name = ligand_name
        self.task_list = []
        self.task_names = set()  # Prevent duplicate measurements

    def add(self, metric, name, prot_sel, lig_sel):
        """
        Add a measurement task to the execution list.
        
        Args:
            metric: Measurement calculation instance
            name (str): Human-readable measurement name
            prot_sel (str): MDAnalysis selection for protein atoms
            lig_sel (str): MDAnalysis selection for ligand atoms
            
        Raises:
            RuntimeError: If measurement name is already used
        """
        if name in self.task_names:
            raise RuntimeError("Repeat of measurement '{}'".format(name))

        self.task_names.add(name)
        self.task_list.append(Task(metric, name, prot_sel, lig_sel))

    def run(self, job_data):
        """
        Execute all measurements for protein-ligand docking results.
        
        This method:
        1. Loads protein and ligand structures from docking results
        2. Validates atom selections for all measurement tasks
        3. Computes measurements across all docking poses
        4. Stores results in job data for database integration
        5. Handles cases where no valid evaluations exist
        
        Args:
            job_data (DataContainer): Job data with modeling and docking results
            
        Returns:
            DataContainer: Updated job data with measurement results
        """
        job_dir = Path(job_data.job_dir)
        modeling = job_data.modeling
        lig_name = self.ligand_name

        # Find representative structures for atom selection validation
        ligand_pdb = None
        protein_pdb = None
        for model in modeling.models:
            if hasattr(model, 'evals') and lig_name in model.evals:
                ligand_pdb = model.evals[lig_name].pdb
                protein_pdb = model.pdb
                break

        # Skip if no docking evaluations found
        if ligand_pdb is None or protein_pdb is None:
            return job_data

        # Load structures and validate atom selections
        protein = mda.Universe(str(job_dir / protein_pdb))
        ligand = mda.Universe(str(job_dir / ligand_pdb))

        for task in self.task_list:
            task.set_system(protein, ligand)
            if not task.enabled:
                print("Could not apply measurement '{}' to variant '{}' - invalid selections".format(
                    task.name, job_data.variant.name))

        # Compute measurements for each model and docking evaluation
        for model in modeling.models:
            if not hasattr(model, 'evals') or lig_name not in model.evals:
                continue

            evaluation = model.evals[lig_name]
            
            # Load model-specific structures
            protein.load_new(str(job_dir / model.pdb))
            ligand.load_new(str(job_dir / evaluation.pdb))

            # Compute measurements for each docking pose
            for pose_frame in ligand.trajectory:
                for task in self.task_list:
                    if task.enabled:
                        task.compute()

            # Store measurement results in evaluation data
            if not hasattr(evaluation, 'measurements'):
                evaluation.measurements = []

            for task in self.task_list:
                if task.enabled:
                    evaluation.measurements.append(task.to_container())
                task.clear()  # Reset for next model

        return job_data
