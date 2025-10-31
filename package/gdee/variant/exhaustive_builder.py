"""
Exhaustive variant enumeration for GDEE platform.

This module provides a variant builder that systematically generates all possible
amino acid combinations at selected positions using combinatorial enumeration.
"""


import itertools
from .base_builder import BaseBuilder
from .sequence import Blosum, ResidueIndex
from gdee.misc import DataContainer, get_valid_filename


class CombinatorialMutation:
    """
    Iterator for exhaustive amino acid combinations.
    
    This class generates all possible combinations of amino acids at
    selected positions. It handles both position combinations and
    amino acid permutations for systematic variant enumeration.
    
    Attributes:
        k (int): Number of positions to mutate simultaneously
        data_size (int): Total number of mutable positions
        index_iter (iterator): Position combination iterator
        indices (tuple): Current position combination
        mutations (iterator): Amino acid permutation iterator
        stop (bool): Enumeration completion flag
    """
    
    def __init__(self, size, data_size):
        """
        Initialize combinatorial mutation generator.
        
        Args:
            size (int): Number of simultaneous mutations
            data_size (int): Number of available mutation positions
        """
        self.k = size
        self.data_size = data_size
        self.index_iter = itertools.combinations(range(data_size), self.k)
        self.indices = []
        self.mutations = None
        self.stop = False
        self.change_group()

    def change_group(self):
        """
        Move to next position combination and reset amino acid iterator.
        
        When all amino acid combinations for current positions are exhausted,
        advance to the next set of positions and restart amino acid enumeration.
        """
        self.mutations = itertools.product(Blosum()[62][0], repeat=self.k)
        try:
            self.indices = next(self.index_iter)
        except StopIteration:
            self.stop = True

    def __iter__(self):
        return self

    def __next__(self):
        """
        Get next mutation combination.
        
        Returns:
            zip: Pairs of (position_index, amino_acid) for current combination
            
        Raises:
            StopIteration: When all combinations are exhausted
        """
        try:
            mutations = next(self.mutations)
        except StopIteration:
            self.change_group()
            mutations = next(self.mutations)

        if self.stop:
            raise StopIteration()

        return zip(self.indices, mutations)


class ExhaustiveBuilder(BaseBuilder):
    """
    Generates all possible amino acid variants through exhaustive enumeration.
    
    This builder systematically creates every possible combination of amino
    acids at selected positions. It's useful for comprehensive exploration
    of small mutation spaces but computationally intensive for large selections.
    
    Attributes:
        variant (ProtSeq): Working copy for mutation application
        wildtype_sel (list): Selected wildtype positions
        variant_sel (list): Corresponding variant positions
        fixed_index (list): Indices of residues kept fixed
        combinations (CombinatorialMutation): Systematic enumeration iterator
    """
    
    def __init__(self, parameters, database):
        """
        Initialize exhaustive variant builder.
        
        Args:
            parameters (dict): Configuration including selection and combination size
            database (Database): Database connection for result storage
        """
        super().__init__(parameters, database)
        self.variant = None
        self.wildtype_sel = []
        self.variant_sel = []
        self.fixed_index = []
        self.combinations = None

    def special_initialize(self):
        """
        Set up exhaustive enumeration data structures.
        
        Parses residue selections, validates fixed vs mutable positions,
        and configures the combinatorial mutation iterator.
        
        Raises:
            RuntimeError: If no residues selected or fixed/mutable conflicts
        """
        self.variant = self.protein.copy()

        if not self.parameters["selection"]:
            raise RuntimeError("No residues were selected to be mutated")

        # Parse mutable residue selection
        selection = ResidueIndex(self.protein, self.parameters["selection"])
        self.wildtype_sel = selection.apply(self.protein)
        self.variant_sel = selection.apply(self.variant)
        for wt, mut in zip(self.wildtype_sel, self.variant_sel):
            assert wt == mut  # Order consistency check

        # Parse fixed residue selection and validate compatibility
        fixed = ResidueIndex(self.protein, self.parameters["fixed"])
        fixed_sel = fixed.apply(self.protein)
        for res in fixed_sel:
            if res in self.wildtype_sel:
                raise RuntimeError("Residue {}:{} marked as fixed and mutable".format(res.chain, res.resid))

        self.fixed_index = [res.index for res in fixed_sel]

        # Configure combinatorial enumeration
        size = len(self.wildtype_sel)
        k = self.parameters["combinations"]
        if k > 0 and k < size:
            self.combinations = CombinatorialMutation(k, size)
        else:
            # Mutate all selected positions simultaneously
            self.combinations = CombinatorialMutation(size, size)

    def mutations(self):
        """
        Generate mutation description for current variant state.
        
        Returns:
            tuple: (mutation_string, mutation_indices) describing current mutations
        """
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

    def apply_mutations(self, rules):
        """
        Apply specific amino acid substitutions to variant sequence.
        
        Args:
            rules (iterable): Pairs of (position_index, amino_acid_code)
                             from CombinatorialMutation iterator
        """
        # Reset to wildtype sequence
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            mut_res.code = wt_res.code

        # Apply new mutations respecting exclusion rules
        for pos, new_code in rules:
            res = self.variant_sel[pos]

            if not self.is_excluded(res, new_code):
                res.code = new_code

    def fetch_next_job(self):
        """
        Generate next variant in exhaustive enumeration.
        
        Applies next combination from iterator, ensures uniqueness,
        and creates job data for pipeline processing.
        
        Returns:
            DataContainer or None: Job data with variant information,
                                   or None when enumeration is complete
        """
        while True:
            mut_name, mut_index = self.mutations()
            if not self.variant_exists(mut_name):
                break

            try:
                self.apply_mutations(next(self.combinations))
            except StopIteration:
                return None

        # Create job data container
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
