"""
Protein variant generation for the GDEE platform.

This module provides factory classes and builders for generating protein variants
through different strategies: MSA-based (FASTA file resulting from the BLAST search), random mutations, and exhaustive mutations.
"""

from ..database import Database
from .msa_builder import MSABuilder
from .mutation_builder import MutationBuilder
from .exhaustive_builder import ExhaustiveBuilder

__all__ = ["VariantBuilderFactory"]


class VariantBuilderFactory:
    """
    Factory for creating variant generation strategies.
    
    This factory creates different types of variant builders based on
    configuration parameters. It handles database setup and amino acid
    exclusion rules for all builder types.
    
    Attributes:
        parameters (dict): Configuration parameters including builder type,
                          PDB file, exclusion rules, and strategy-specific settings
    """
    
    def __init__(self):
        """Initialize factory with empty parameters."""
        self.parameters = {}

    def make(self):
        """
        Create and configure a variant builder.
        
        Validates configuration, sets up amino acid exclusions, creates
        database connection, and instantiates the appropriate builder type.
        
        Returns:
            BaseBuilder: Configured variant builder (MSABuilder, MutationBuilder, 
                        or ExhaustiveBuilder)
            
        Raises:
            RuntimeError: If PDB file is invalid or builder type is unknown
        """
        if self.parameters["pdb_file"] is None:
            raise RuntimeError("Invalid template PDB file.")

        # Combine per-residue and global amino acid exclusions
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
            raise RuntimeError("Variant builder '{}' does not exist.".format(name))
