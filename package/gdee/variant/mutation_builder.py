"""
Random mutation-based variant generator for GDEE platform.

This module provides a variant builder that generates protein variants through
random amino acid substitutions using substitution matrices like BLOSUM62.
"""


import itertools
from .base_builder import BaseBuilder
from .sequence import ResidueIndex, MatrixMutation, Blosum62Mutation
from gdee.misc import DataContainer, get_valid_filename


class MutationBuilder(BaseBuilder):
    """
    Generates protein variants through random amino acid mutations.
    
    This builder uses substitution matrices (typically BLOSUM62) to generate
    biologically reasonable amino acid substitutions. It supports both single
    and multiple simultaneous mutations with configurable conservative/
    non-conservative mutation bias.
    
    Attributes:
        variant (ProtSeq): Working copy of protein sequence for mutations
        wildtype_sel (list): Selected wildtype residues for mutation
        variant_sel (list): Corresponding positions in variant sequence
        matrix (MatrixMutation): Substitution matrix for mutation probabilities
        invert_weights (bool): Use conservative vs non-conservative mutations
        max_iter (int): Maximum number of variants to generate
        combinations (iterator): For multi-residue simultaneous mutations
        iterations (int): Current number of generated variants
        fixed_index (list): Residue indices to keep fixed during optimization
    """
    
    def __init__(self, parameters, database):
        """
        Initialize mutation-based variant builder.
        
        Args:
            parameters (dict): Configuration including matrix type, selection,
                             iteration limits, and mutation strategy
            database (Database): Database connection for result storage
        """
        super().__init__(parameters, database)
        self.variant = None
        self.wildtype_sel = []
        self.variant_sel = []
        self.invert_weights = not self.parameters["conservative"]
        self.max_iter = self.parameters["max_iterations"]
        self.combinations = None
        self.combinations_size = self.parameters["combinations"]
        self.iterations = 0
        self.fixed_index = []

        # Set up mutation matrix
        if self.parameters["matrix"] == "blosum62":
            self.matrix = Blosum62Mutation()
        else:
            self.matrix = MatrixMutation(
                self.parameters["matrix_aa"],
                self.parameters["matrix_weights"]
            )

    def next_sel(self):
        """
        Generate next set of residues for simultaneous mutation.
        
        Handles both single residue mutations and combinations of multiple
        residues for simultaneous substitution.
        
        Yields:
            tuple: Pairs of (wildtype_residue, variant_residue) for mutation
        """
        # Clear previous mutations that won't be selected
        for wt_res, mut_res in zip(self.wildtype_sel, self.variant_sel):
            mut_res.code = wt_res.code

        if self.combinations is None:
            # Single residue mutations
            return zip(self.wildtype_sel, self.variant_sel)

        # Multi-residue combinations
        indices = next(self.combinations)
        wildtype = tuple(self.wildtype_sel[i] for i in indices)
        variant = tuple(self.variant_sel[i] for i in indices)
        return zip(wildtype, variant)

    def special_initialize(self):
        """
        Set up mutation-specific data structures.
        
        Parses residue selection strings, validates fixed vs mutable residues,
        and configures multi-residue combination strategies.
        
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

        # Set up multi-residue combinations if requested
        size = len(self.variant_sel)
        if self.combinations_size > 0 and self.combinations_size < size:
            self.combinations = itertools.cycle(itertools.combinations(range(size), self.combinations_size))
        else:
            self.combinations = None

        # Parse fixed residue selection and validate compatibility
        fixed = ResidueIndex(self.protein, self.parameters["fixed"])
        fixed_sel = fixed.apply(self.protein)
        for res in fixed_sel:
            if res in self.wildtype_sel:
                raise RuntimeError("Residue {}:{} marked as fixed and mutable".format(res.chain, res.resid))

        self.fixed_index = [res.index for res in fixed_sel]

    def mutations(self):
        """
        Generate mutation string and index tuple for current variant.
        
        Creates human-readable mutation names like "A:A123V|B:T456K"
        and corresponding residue indices for MODELLER optimization.
        
        Returns:
            tuple: (mutation_string, mutation_indices) or (protein_name, all_indices)
                   for wildtype
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

    def fetch_next_job(self):
        """
        Generate next random mutation variant.
        
        Creates random mutations using the substitution matrix, ensures
        uniqueness, and respects amino acid exclusion rules.
        
        Returns:
            DataContainer or None: Job data with variant information,
                                   or None if maximum iterations reached
        """
        if self.iterations >= self.max_iter:
            return None

        max_tries = 20 * len(self.wildtype_sel)  # Avoid infinite loops
        while True:
            mut_name, mut_index = self.mutations()
            if not self.variant_exists(mut_name):
                break

            if max_tries == 0:
                return None
            max_tries -= 1

            # Generate new random mutations
            for wt_res, mut_res in self.next_sel():
                while True:
                    code = self.matrix.mutate(wt_res.code, not self.invert_weights)
                    if not self.is_excluded(mut_res, code):
                        mut_res.code = code
                        break

        # Create job data container
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
