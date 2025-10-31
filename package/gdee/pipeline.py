"""
Pipeline orchestration for GDEE workflows.

This module provides factory and pipeline classes for coordinating the complete
GDEE workflow, including variant generation, structure modeling,
quality assessment, molecular docking, and measurement analysis.
"""


from gdee import files
from gdee.variant import VariantBuilderFactory
from gdee.modeling import ModelBuilderFactory, ModelQualityBuilderFactory
from gdee.evaluator import EvaluatorFactory
from gdee.measurement import MeasurerFactory
from path import Path
import tarfile
import os
import traceback


__all__ = ["PipelineFactory"]


class PipelineFactory:
    """
    Factory class for creating and configuring GDEE pipelines.
    
    This factory orchestrates the creation of all pipeline components including
    variant builders, modeling tools, quality checkers, molecular docking
    evaluators, and measurement systems.
    
    Attributes:
        protein_name (str): Name identifier for the target protein
        programs (dict): Paths to external programs (vina, voromqa, etc.)
        work_dir (str): Working directory path for pipeline execution
        pdb (str): Path to the template PDB structure file
        ligands (list): List of Ligand objects for docking calculations
        db_file (str): Path to SQLite database for results storage
        io (dict): I/O configuration for result archiving
        variant_parameters (dict): Configuration for variant generation
        model_parameters (dict): Configuration for 3D modeling
        model_quality_parameters (dict): Configuration for quality assessment
        evaluator_parameters (dict): Configuration for molecular docking
    """
    
    def __init__(self):
        """Initialize factory with default empty configuration."""
        self.protein_name = None
        self.programs = {}
        self.work_dir = None
        self.pdb = None
        self.ligands = []
        self.db_file = None
        self.io = {}
        self.variant_parameters = {}
        self.model_parameters = {}
        self.model_quality_parameters = {}
        self.evaluator_parameters = {}

    def make(self):
        """
        Create and configure a complete protein engineering pipeline.
        
        This method:
        1. Sets up the working directory and archiving system
        2. Creates variant builder using VariantBuilderFactory
        3. Configures 3D modeling with ModellerBuilder
        4. Sets up quality assessment with VoroMQA/DOPE scoring
        5. Creates docking evaluators for each ligand
        6. Sets up measurement systems for pose filtering
        
        Returns:
            Pipeline: Fully configured pipeline ready for execution
            
        Raises:
            RuntimeError: If required configuration parameters are missing
        """
        self.pdb = Path(self.pdb).abspath()
        pipeline = Pipeline()
        base_dir = Path(self.work_dir).abspath()
        pipeline.work_dir = base_dir / "files"
        pipeline.work_dir.makedirs_p()

        # Configure result archiving system
        pipeline.archiver = files.Archiver(self.io["output"],
                                           self.io["output_format"],
                                           self.io["output_freq"])

        # Create variant builder (MSA, mutation, or exhaustive)
        variant_factory = VariantBuilderFactory()
        self.variant_parameters["db_file"] = self.db_file
        self.variant_parameters["pdb_file"] = self.pdb
        self.variant_parameters["protein_name"] = self.protein_name
        variant_factory.parameters = self.variant_parameters
        pipeline.variant_builder = variant_factory.make()

        # Create 3D modeling component (MODELLER)
        model_factory = ModelBuilderFactory()
        self.model_parameters["pdb_file"] = self.pdb
        model_factory.parameters = self.model_parameters
        pipeline.add_task(model_factory.make())

        # Create quality assessment component (VoroMQA, DOPE)
        quality_factory = ModelQualityBuilderFactory()
        self.model_quality_parameters["programs"] = self.programs
        quality_factory.parameters = self.model_quality_parameters
        pipeline.add_task(quality_factory.make())

        # Create docking and measurement components for each ligand
        self.evaluator_parameters.update(self.programs)
        evaluator_factory = EvaluatorFactory()
        measurer_factory = MeasurerFactory()
        for ligand in self.ligands:
            # Configure docking evaluator (Vina/Vinardo)
            evaluator_parameters = self.evaluator_parameters.copy()
            evaluator_parameters["ligand"] = ligand
            evaluator_factory.parameters = evaluator_parameters
            pipeline.add_task(evaluator_factory.make())

            # Configure measurement system for distance calculations
            measurer_factory.ligand = ligand
            pipeline.add_task(measurer_factory.make())

        return pipeline


class Pipeline:
    """
    Main pipeline executor for GDEE workflows.
    
    This class coordinates the execution of all pipeline tasks in sequence:
    variant generation → 3D modeling → quality assessment → docking → measurements.
    It handles job scheduling, error management, result archiving, and database storage.
    
    Attributes:
        database: Database file path 
        work_dir (Path): Working directory
        archiver (Archiver): File archiving system for results
        task_list (list): Ordered list of pipeline tasks to execute
        _variant_builder: Variant generation component
        _terminate (bool): Flag for graceful termination
    """
    
    def __init__(self):
        """Initialize pipeline with default configuration."""
        self.database = None
        self.work_dir = Path().abspath()
        self.archiver = files.Archiver("files", ".{:06d}", 1000)
        self._variant_builder = None
        self.task_list = []
        self._terminate = False

    @property
    def variant_builder(self):
        """Get the variant builder component."""
        return self._variant_builder

    @variant_builder.setter
    def variant_builder(self, obj):
        """Set the variant builder component."""
        self._variant_builder = obj

    def add_task(self, task):
        """
        Add a processing task to the pipeline.
        
        Tasks are executed in the order they are added:
        1. ModellerBuilder (3D modeling)
        2. ModelQualityChecker (quality assessment)  
        3. VinaDocking/VinardoDocking (molecular docking)
        4. Measurer (distance measurements)
        
        Args:
            task: Pipeline task object with run() method
        """
        self.task_list.append(task)

    def next_job(self, size):
        """
        Get the next batch of variants to process.
        
        Args:
            size (int): Maximum number of variants to retrieve
            
        Returns:
            list: List of job data objects for processing
        """
        job_list = []

        if self._terminate:
            return job_list

        for i in range(size):
            job = self.variant_builder.next_job()
            if job is None:
                break

            job_list.append(job)

        return job_list

    def run_pipeline(self, job_data):
        """
        Execute the complete pipeline for a single variant.
        
        This method:
        1. Creates working directory for the variant
        2. Executes all pipeline tasks in sequence
        3. Handles errors and sets fatal_error flag
        4. Returns processed job data with results
        
        Args:
            job_data: Job data container with variant information
            
        Returns:
            job_data: Updated container with processing results
        """
        try:
            job_dir = self.work_dir / job_data.variant_dir
            job_data.job_dir = job_dir
            job_dir.makedirs_p()
            job_data.fatal_error = False

            # Execute all pipeline tasks sequentially
            for step in self.task_list:
                job_data = step.run(job_data)
                if job_data.fatal_error:
                    return job_data

        except Exception as error:
            job_data.fatal_error = True
            print("Exception caught:", traceback.print_exc(), error, sep="\n")

        return job_data

    def save_results(self, data):
        """
        Save pipeline results to database and archive files.
        
        This method:
        1. Calls variant builder to save results to database
        2. Archives successful job directories for storage
        3. Handles errors during save operations
        4. Cleans up failed results from database
        
        Args:
            data (list): List of processed job data containers
        """
        for result in data:
            try:
                # Save results to SQLite database
                self._variant_builder.save_results(result)

                # Archive job directory if processing was successful
                if not result.fatal_error:
                    self.archiver.add(result.job_dir, result.variant_dir)

            except Exception as error:
                print("Exception caught:", traceback.print_exc(), error, sep="\n")
                # Clean up partial database entries on error
                self._variant_builder.unsave_results(data)

    def terminate(self):
        """Set termination flag for graceful shutdown."""
        self._terminate = True

    def finalize(self):
        """Finalize pipeline execution and close archive files."""
        self.archiver.finalize()
