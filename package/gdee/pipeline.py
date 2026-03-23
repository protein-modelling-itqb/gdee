"""Pipeline orchestration for workflow execution."""

from gdee import files
from gdee.variant import VariantBuilderFactory
from gdee.modeling import ModelBuilderFactory, ModelQualityBuilderFactory
from gdee.evaluator import EvaluatorFactory
from gdee.measurement import MeasurerFactory
from path import Path
import tarfile
import os
import traceback
from gdee.analysis.rescoring import Rescore


__all__ = ["PipelineFactory", "RescoreFactory"]


class PipelineFactory:
    """Factory for creating pipeline instances."""
    def __init__(self):
        """Initialize pipeline factory."""
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
        """Create and configure pipeline instance.

        Returns:
            Pipeline: Configured pipeline object
        """
        self.pdb = Path(self.pdb).absolute()
        pipeline = Pipeline()
        base_dir = Path(self.work_dir).absolute()
        pipeline.work_dir = base_dir / "files"
        pipeline.work_dir.makedirs_p()

        pipeline.archiver = files.Archiver(self.io["output"],
                                           self.io["output_format"],
                                           self.io["output_freq"])

        variant_factory = VariantBuilderFactory()
        self.variant_parameters["db_file"] = self.db_file
        self.variant_parameters["pdb_file"] = self.pdb
        self.variant_parameters["protein_name"] = self.protein_name
        variant_factory.parameters = self.variant_parameters
        pipeline.variant_builder = variant_factory.make()

        model_factory = ModelBuilderFactory()
        self.model_parameters["pdb_file"] = self.pdb
        model_factory.parameters = self.model_parameters
        pipeline.add_task(model_factory.make())

        quality_factory = ModelQualityBuilderFactory()
        self.model_quality_parameters["programs"] = self.programs
        quality_factory.parameters = self.model_quality_parameters
        pipeline.add_task(quality_factory.make())

        self.evaluator_parameters.update(self.programs)
        evaluator_factory = EvaluatorFactory()
        measurer_factory = MeasurerFactory()
        for ligand in self.ligands:
            evaluator_parameters = self.evaluator_parameters.copy()
            evaluator_parameters["ligand"] = ligand
            evaluator_factory.parameters = evaluator_parameters
            pipeline.add_task(evaluator_factory.make())

            measurer_factory.ligand = ligand
            pipeline.add_task(measurer_factory.make())

        return pipeline


class Pipeline:
    """Main pipeline for processing protein variants."""
    def __init__(self):
        """Initialize pipeline."""
        self.database = None
        self.work_dir = Path().absolute()
        self.archiver = files.Archiver("files", ".{:06d}", 1000)
        self._variant_builder = None
        self.task_list = []
        self._terminate = False

    @property
    def variant_builder(self):
        """Get variant builder.

        Returns:
            BaseBuilder: Variant generation strategy
        """
        return self._variant_builder

    @variant_builder.setter
    def variant_builder(self, obj):
        """Set variant builder.

        Args:
            obj: Variant builder instance
        """
        self._variant_builder = obj

    def add_task(self, task):
        """Add processing task to pipeline.

        Args:
            task: Task object with run() method
        """
        self.task_list.append(task)

    def next_job(self, size):
        """Get next batch of jobs.

        Args:
            size: Number of jobs to retrieve

        Returns:
            list: List of job data containers
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
        """Execute all pipeline tasks on job.

        Args:
            job_data: Job data to process

        Returns:
            DataContainer: Processed job data
        """
        try:
            job_dir = self.work_dir / job_data.variant_dir
            job_data.job_dir = job_dir
            job_dir.makedirs_p()
            job_data.fatal_error = False

            for step in self.task_list:
                job_data = step.run(job_data)
                if job_data.fatal_error:
                    return job_data

        except Exception as error:
            job_data.fatal_error = True
            print("Exception caught:", traceback.print_exc(), error, sep="\n")

        return job_data

    def save_results(self, data):
        """Save processing results.

        Args:
            data: List of job results
        """
        for result in data:
            try:
                self._variant_builder.save_results(result)

                if not result.fatal_error:
                    self.archiver.add(result.job_dir, result.variant_dir)

            except Exception as error:
                print("Exception caught:", traceback.print_exc(), error, sep="\n")
                self._variant_builder.unsave_results(data)

    def terminate(self):
        """Request pipeline termination."""
        self._terminate = True

    def finalize(self):
        """Finalize pipeline execution."""
        self.archiver.finalize()


class RescoreFactory:
    """Factory for creating rescoring pipeline."""
    def __init__(self):
        """Initialize rescoring factory."""
        self.table = None
        self.pickle_path = None
        self.input_db = None
        self.output_db = None
        self.files_path = None

    def make(self):
        """Create rescoring pipeline.

        Returns:
            RescorePipeline: Configured rescoring pipeline
        """
        pipeline = RescorePipeline()

        rescore = Rescore(self.output_db, self.table)
        rescore.pickle_path = self.pickle_path
        rescore.path = Path(self.files_path).absolute()
        rescore.get_variants(self.input_db)

        pipeline.rescore = rescore

        return pipeline

class RescorePipeline:
    """Pipeline for re-scoring docking poses."""
    def __init__(self):
        """Initialize rescoring pipeline."""
        self.work_dir = Path().absolute()
        self.rescore = None
        self._terminate = False
        self._created = False

    def next_job(self, size):
        """Get next batch of rescoring jobs.

        Args:
            size: Number of jobs to retrieve

        Returns:
            list: List of job data containers
        """
        job_list = []

        if self._terminate:
            return job_list

        for i in range(size):
            job = self.rescore.fetch_next_job()
            if job is None:
                break

            job_list.append(job)

        return job_list

    def run_pipeline(self, job_data):
        """Execute rescoring on job.

        Args:
            job_data: Job data with poses to rescore

        Returns:
            DataContainer: Job data with rescored results
        """
        try:
            job_data.fatal_error = False

            job_data = self.rescore.run(job_data)
            if job_data.fatal_error:
                return job_data

        except Exception as error:
            job_data.fatal_error = True
            print("Exception caught:", traceback.print_exc(), error, sep="\n")

        return job_data

    def save_results(self, data):
        """Save rescoring results.

        Args:
            data: List of job results
        """
        for result in data:
            try:
                self.rescore.save_results(result)

            except Exception as error:
                print("Exception caught:", traceback.print_exc(), error, sep="\n")

    def finalize(self):
        """Finalize pipeline execution."""
        self._terminate = True
