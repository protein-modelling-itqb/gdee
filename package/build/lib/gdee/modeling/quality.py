"""
"""
import subprocess
from path import Path


def external_command(arguments, name):
    proc = subprocess.run(
        arguments,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if proc.returncode:
        raise RuntimeError("Error processing job '{}':\n{}\n".format(name, proc.stderr.decode("UTF-8")))

    return proc


class BaseComparison:
    def __init__(self, coff):
        self.coff = coff

    def __call__(self, value):
        raise NotImplementedError()


class HigherThan(BaseComparison):
    def __call__(self, value):
        return value > self.coff


class LowerThan(BaseComparison):
    def __call__(self, value):
        return value < self.coff


class ModelQualityChecker:
    def __init__(self):
        self.cutoff = {}
        self.scorers = []

    def enable_scorer(self, name, program):
        if name == "voromqa":
            self.scorers.append(VoroMQA(program))

        else:
            raise ValueError("Invalid scorer: {}".format(name))

    def add_cutoff(self, name, comparison):
        self.cutoff[name] = comparison

    def add_lower_cutoff(self, name, cutoff):
        self.add_cutoff(name, LowerThan(cutoff))

    def add_higher_cutoff(self, name, cutoff):
        self.add_cutoff(name, HigherThan(cutoff))

    def run(self, job_data):
        for model in job_data.modeling.models:
            for scorer in self.scorers:
                value = scorer(job_data.job_dir / model.pdb)
                model.scores[scorer.name] = value

            accepted = True
            for name, comparison in self.cutoff.items():
                accepted = accepted and comparison(model.scores[name])

            model.rejected = not accepted

        return job_data


class VoroMQA:
    name = "voromqa"

    def __init__(self, program):
        self.program = program

    def __call__(self, model_pdb):
        command = [
            self.program,
            "--processors", "1",
            "-i", model_pdb
        ]

        process = external_command(command, "voromqa")
        return float(process.stdout.decode().split()[1])
