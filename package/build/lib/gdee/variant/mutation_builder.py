"""
"""


import itertools
from .base_builder import BaseBuilder
from .sequence import ResidueIndex, MatrixMutation, Blosum62Mutation
from gdee.misc import DataContainer, get_valid_filename


class MutationBuilder(BaseBuilder):
    def __init__(self, parameters, database):
        super().__init__(parameters, database)
        self.variant = None
        self.mut_sel = []
        self.invert_weights = self.parameters["conservative"]
        self.max_iter = self.parameters["max_iterations"]
        self.combinations = None
        self.combinations_size = self.parameters["combinations"]
        self.iterations = 0

        if self.parameters["matrix"] == "blosum62":
            self.matrix = Blosum62Mutation()

        else:
            self.matrix = MatrixMutation(
                self.parameters["matrix_aa"],
                self.parameters["matrix_weights"]
            )

    def next_sel(self):
        # Clear previous mutations that won't be selected
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            mut_res.code = wt_res.code

        if self.combinations is None:
            return zip(self.wildtype_sel, self.variant_sel)

        indices = next(self.combinations)
        wildtype = tuple(self.wildtype_sel[i] for i in indices)
        variant = tuple(self.variant_sel[i] for i in indices)
        return zip(wildtype, variant)

    def special_initialize(self):
        self.variant = self.protein.copy()

        if not self.parameters["selection"]:
            raise RuntimeError("No residues were selected to be mutated")

        selection = ResidueIndex(self.protein, self.parameters["selection"])
        self.wildtype_sel = selection.apply(self.protein)
        self.variant_sel = selection.apply(self.variant)
        for wt, mut in zip(self.wildtype_sel, self.variant_sel):
            assert wt == mut # Order check

        size = len(self.variant_sel)
        if self.combinations_size > 0 and self.combinations_size < size:
            self.combinations = itertools.cycle(itertools.combinations(range(size), self.combinations_size))

        else:
            self.combinations = None

        fixed = ResidueIndex(self.protein, self.parameters["fixed"])
        fixed_sel = fixed.apply(self.protein)
        for res in fixed_sel:
            if res in self.wildtype_sel:
                raise RuntimeError("Residue {}:{} marked as fixed and mutable".format(res.chain, res.resid))

        self.fixed_index = [res.index for res in fixed_sel]

    def mutations(self):
        mutations = []
        mut_index = []
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            if wt_res.code != mut_res.code:
                mutations.append("{}:{}{}{}".format(wt_res.chain, wt_res.code, wt_res.resid, mut_res.code))
                mut_index.append(mut_res.index)
        mut_index.sort()

        if mutations:
            return "|".join(mutations), tuple(mut_index)

        return self.protein.name, tuple(res.index for res in self.wildtype_sel)

    def fetch_next_job(self):
        if self.iterations >= self.max_iter:
            return None

        max_tries = 20 * len(self.wildtype_sel) # We can only go linear
        while True:
            mut_name, mut_index = self.mutations()
            if not self.variant_exists(mut_name):
                break

            if max_tries == 0:
                return None
            max_tries -= 1

            for wt_res, mut_res in self.next_sel():
                while True:
                    code = self.matrix.mutate(wt_res.code, self.invert_weights)
                    if not self.is_excluded(mut_res, code):
                        mut_res.code = code
                        break

        variant_dir = get_valid_filename(mut_name.replace("|", "_"))
        variant = self.variant.copy()
        variant.name = mut_name
        self.iterations += 1

        job = DataContainer()
        job.variant_dir = variant_dir
        job.wildtype = self.protein.copy()
        job.variant = variant
        job.is_wildtype = mut_name == self.protein.name
        job.mut_index = mut_index
        job.fixed_index = self.fixed_index.copy()

        return job
