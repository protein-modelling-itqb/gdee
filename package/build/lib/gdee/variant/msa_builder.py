"""
"""


import warnings
from .base_builder import BaseBuilder
from gdee.misc import DataContainer, get_valid_filename

from Bio import SeqIO, Align, BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonWarning)
    from Bio.Align import substitution_matrices


class MSABuilder(BaseBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.wt_seq = ""
        self.msa = tuple()
        self._iter = iter(self.msa)
        self._is_wildtype = True

    def special_initialize(self):
        self.wt_seq = self.protein.to_modeller().replace("/", "")
        self.msa = tuple(SeqIO.parse(self.parameters["msa"], "fasta"))
        self._iter = iter(self.msa)

    def variant_from_alignment(self, name, other_seq):
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -10
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        alignment = aligner.align(self.wt_seq.replace(".", ""), other_seq)
        if not alignment:
            raise RuntimeError("No suitable alignment found for sequence '{}'".format(other_seq))

        variant = self.protein.copy()
        variant.name = name
        variant_iter = iter(variant.flatten())

        is_wildtype = True
        for query, target in zip(alignment[0][0], alignment[0][1]):
            if query == "-":
                continue

            seq_pos = next(variant_iter)
            while seq_pos.is_blk:
                seq_pos = next(variant_iter)
            code = seq_pos.code
            assert code == query

            if target != code:
                is_wildtype = False

                if target != "-":
                    seq_pos.code = target

        return variant, is_wildtype

    def fetch_next_job(self):
        is_wildtype = False
        if self._is_wildtype:
            self._is_wildtype = False
            is_wildtype = True
            variant = self.protein.copy()

        else:
            try:
                next_seq = next(self._iter)
                variant, is_wildtype = self.variant_from_alignment(next_seq.name, str(next_seq.seq))

            except StopIteration:
                return None

        job = DataContainer()
        job.variant_dir = get_valid_filename(variant.name)
        job.wildtype = self.protein.copy()
        job.variant = variant
        job.is_wildtype = is_wildtype
        job.mut_index = tuple()
        job.fixed_index = tuple()

        return job
