"""
"""


from ..database import Database
from .msa_builder import MSABuilder
from .mutation_builder import MutationBuilder
from .exhaustive_builder import ExhaustiveBuilder


__all__ = ["VariantBuilderFactory"]


class VariantBuilderFactory:
    def __init__(self):
        self.parameters = {}

    def make(self):
        if self.parameters["pdb_file"] is None:
            raise RuntimeError("Invalid template PDB file.")

        excludes = self.parameters["excluded"]
        all_excludes = self.parameters["excluded_all"]
        for res in self.parameters["selection"].split():
            excludes[res] = excludes.get(res, "") + all_excludes

        database = Database(self.parameters["db_file"])

        name = self.parameters["name"]
        if name == "msa":
            return MSABuilder(self.parameters, database)

        elif name == "mutation":
            return MutationBuilder(self.parameters, database)

        elif name == "exhaustive":
            return ExhaustiveBuilder(self.parameters, database)

        else:
            raise RuntimeError("Variant builder '{}' does not exists.".format(name))
