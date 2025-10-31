"""
"""


from .modeller_builder import ModellerBuilder
from .quality import ModelQualityChecker


__all__ = ["ModelBuilderFactory"]


class ModelBuilderFactory:
    def __init__(self):
        self.parameters = {}

    def make(self):
        if self.parameters["pdb_file"] is None:
            raise RuntimeError("Invalid template PDB file.")

        name = self.parameters["name"]
        if name == "modeller":
            return ModellerBuilder(self.parameters)

        # if name == "rosetta":
        #     return RosettaBuilder(self.parameters)

        else:
            raise RuntimeError("Modeler '{}' does not exists.".format(name))


class ModelQualityBuilderFactory:
    def __init__(self):
        self.parameters = {}

    def make(self):
        checker = ModelQualityChecker()

        if self.parameters["norm_dope"] is not None:
            # Normalized DOPE is computed by Modeller for now
            # checker.enable_scorer("norm_dope")
            checker.add_lower_cutoff("norm_dope", self.parameters["norm_dope"])

        if self.parameters["voromqa"] is not None:
            checker.enable_scorer("voromqa", self.parameters["programs"]["voromqa"])
            checker.add_higher_cutoff("voromqa", self.parameters["voromqa"])

        return checker
