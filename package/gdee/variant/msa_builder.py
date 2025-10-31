"""
FASTA-based variant generator for the GEE platform.

This module provides a variant builder that generates protein variants based
on sequences from FASTA files, using information from the BLAST search to guide variant selection.
"""


import warnings
from .base_builder import BaseBuilder
from gdee.misc import DataContainer, get_valid_filename

from Bio import SeqIO, Align, BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonWarning)
    from Bio.Align import substitution_matrices


class MSABuilder(BaseBuilder):
    """
    Generates protein variants from multiple sequence alignment data.
    
    This builder uses sequences from FASTA files to create variants. It performs pairwise alignment
    between the target sequence and FASTA entries to identify mutations.
    
    Attributes:
        wt_seq (str): Wildtype sequence in flat format
        msa (tuple): Parsed MSA sequences from FASTA file
        _iter (iterator): Iterator over MSA sequences
        _is_wildtype (bool): Flag for generating wildtype first
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialize FASTA-based variant builder.
        
        Args:
            *args, **kwargs: Arguments passed to BaseBuilder
        """
        super().__init__(*args, **kwargs)
        self.wt_seq = ""
        self.msa = tuple()
        self._iter = iter(self.msa)
        self._is_wildtype = True

    def special_initialize(self):
        """
        Load and parse FASTA file.

        Converts wildtype sequence to flat format and loads sequences
        from FASTA file for variant generation.
        """
        self.wt_seq = self.protein.to_modeller().replace("/", "")
        self.msa = tuple(SeqIO.parse(self.parameters["msa"], "fasta"))
        self._iter = iter(self.msa)

    def variant_from_alignment(self, name, other_seq):
        """
        Generate variant from FASTA sequence through pairwise alignment.

        Performs BLOSUM62-based pairwise alignment between wildtype and
        FASTA sequence to identify mutation positions and create variant.

        Args:
            name (str): Sequence identifier from FASTA
            other_seq (str): Amino acid sequence from FASTA

        Returns:
            tuple: (variant_protein, is_wildtype_flag)
            
        Raises:
            RuntimeError: If no suitable alignment can be found
        """
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        
        # Align wildtype sequence with MSA sequence
        alignment = aligner.align(self.wt_seq.replace(".", ""), other_seq)
        if not alignment:
            raise RuntimeError("No suitable alignment found for sequence '{}'".format(other_seq))

        # Create variant by applying mutations from alignment
        variant = self.protein.copy()
        variant.name = name
        variant_iter = iter(variant.flatten())

        is_wildtype = True
        for query, match, target in zip(*alignment[0].format().split()):
            if query == "-":  # Gap in wildtype sequence
                continue

            # Find corresponding position in variant sequence
            seq_pos = next(variant_iter)
            while seq_pos.is_blk:  # Skip unknown positions
                seq_pos = next(variant_iter)
            
            code = seq_pos.code
            assert code == query  # Alignment consistency check

            # Apply mutation if sequences differ
            if target != code:
                is_wildtype = False

                if target != "-":  # Not a gap in target
                    seq_pos.code = target

        return variant, is_wildtype

    def fetch_next_job(self):
        """
        Generate next variant from FASTA data.

        First generates the wildtype sequence, then processes FASTA entries
        sequentially to create variants based on evolutionary sequences.
        
        Returns:
            DataContainer or None: Job data with variant information,
                                   or None when all FASTA sequences processed
        """
        is_wildtype = False
        
        if self._is_wildtype:
            # Generate wildtype first
            self._is_wildtype = False
            is_wildtype = True
            variant = self.protein.copy()

        else:
            # Process next sequence from MSA
            try:
                next_seq = next(self._iter)
                variant, is_wildtype = self.variant_from_alignment(next_seq.name, str(next_seq.seq))

            except StopIteration:
                return None

        # Create job data container
        job = DataContainer()
        job.variant_dir = get_valid_filename(variant.name)
        job.wildtype = self.protein.copy()
        job.variant = variant
        job.is_wildtype = is_wildtype
        job.mut_index = tuple()  # MSA variants don't specify mutation indices
        job.fixed_index = tuple()  # No fixed residues in MSA mode

        return job
