"""
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
    def __init__(self):
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
    def __init__(self):
        self.database = None
        self.work_dir = Path().absolute()
        self.archiver = files.Archiver("files", ".{:06d}", 1000)
        self._variant_builder = None
        self.task_list = []
        self._terminate = False

    @property
    def variant_builder(self):
        return self._variant_builder

    @variant_builder.setter
    def variant_builder(self, obj):
        self._variant_builder = obj

    def add_task(self, task):
        self.task_list.append(task)

    def next_job(self, size):
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
        for result in data:
            try:
                self._variant_builder.save_results(result)

                if not result.fatal_error:
                    self.archiver.add(result.job_dir, result.variant_dir)

            except Exception as error:
                print("Exception caught:", traceback.print_exc(), error, sep="\n")
                self._variant_builder.unsave_results(data)

    def terminate(self):
        self._terminate = True

    def finalize(self):
        self.archiver.finalize()
