"""
"""


import itertools
from .base_builder import BaseBuilder
from .sequence import Blosum, ResidueIndex
from gdee.misc import DataContainer, get_valid_filename


class CombinatorialMutation:
    def __init__(self, size, data_size):
        self.k = size
        self.data_size = data_size
        self.index_iter = itertools.combinations(range(data_size), self.k)
        self.indices = []
        self.mutations = None
        self.stop = False
        self.change_group()

    def change_group(self):
        self.mutations = itertools.product(Blosum()[62][0], repeat=self.k)
        try:
            self.indices = next(self.index_iter)
        except StopIteration:
            self.stop = True

    def __iter__(self):
        return self

    def __next__(self):
        try:
            mutations = next(self.mutations)
        except StopIteration:
            self.change_group()
            mutations = next(self.mutations)

        if self.stop:
            raise StopIteration()

        return zip(self.indices, mutations)


class ExhaustiveBuilder(BaseBuilder):
    def __init__(self, parameters, database):
        super().__init__(parameters, database)
        self.combinations = None

    def special_initialize(self):
        self.variant = self.protein.copy()

        if not self.parameters["selection"]:
            raise RuntimeError("No residues were selected to be mutated")

        selection = ResidueIndex(self.protein, self.parameters["selection"])
        self.wildtype_sel = selection.apply(self.protein)
        self.variant_sel = selection.apply(self.variant)
        for wt, mut in zip(self.wildtype_sel, self.variant_sel):
            assert wt == mut # Order check

        fixed = ResidueIndex(self.protein, self.parameters["fixed"])
        fixed_sel = fixed.apply(self.protein)
        for res in fixed_sel:
            if res in self.wildtype_sel:
                raise RuntimeError("Residue {}:{} marked as fixed and mutable".format(res.chain, res.resid))

        self.fixed_index = [res.index for res in fixed_sel]

        size = len(self.wildtype_sel)
        k = self.parameters["combinations"]
        if k > 0 and k < size:
            self.combinations = CombinatorialMutation(k, size)
        else:
            self.combinations = CombinatorialMutation(size, size)

    def mutations(self):
        mutations = []
        mut_index = []
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            if wt_res.code != mut_res.code:
                mutations.append("{}:{}{}{}".format(wt_res.chain, wt_res.code, wt_res.resid, mut_res.code))
                mut_index.append(mut_res.index)
        mut_index.sort()

        if mutations:
            # return "|".join(mutations), tuple(res.index for res in self.wildtype_sel)
            return "|".join(mutations), tuple(mut_index)

        return self.protein.name, tuple(res.index for res in self.wildtype_sel)

    def apply_mutations(self, rules):
        # Clear previous mutations that won't be selected
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            mut_res.code = wt_res.code

        for pos, new_code in rules:
            res = self.variant_sel[pos]

            if not self.is_excluded(res, new_code):
                res.code = new_code

    def fetch_next_job(self):
        while True:
            mut_name, mut_index = self.mutations()
            if not self.variant_exists(mut_name):
                break

            try:
                self.apply_mutations(next(self.combinations))
            except StopIteration:
                return None

        variant_dir = get_valid_filename(mut_name.replace("|", "_"))
        variant = self.variant.copy()
        variant.name = mut_name

        job = DataContainer()
        job.variant_dir = variant_dir
        job.wildtype = self.protein.copy()
        job.variant = variant
        job.is_wildtype = mut_name == self.protein.name
        job.mut_index = mut_index
        job.fixed_index = self.fixed_index.copy()

        return job
