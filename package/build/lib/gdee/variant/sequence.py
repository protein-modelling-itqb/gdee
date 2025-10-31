"""
"""


import numpy as np
import MDAnalysis as mda
import io
import os
import copy
import random
import string
import pkgutil
import itertools
from collections import defaultdict
import warnings

warnings.filterwarnings("ignore", module=r"MDAnalysis.*")


__all__ = ["three_to_one", "SEQ_3_1", "ProtSeq", "ChainSeq", "SeqPos", "ResidueIndex", "MatrixMutation", "Blosum62Mutation"]


SEQ_1_3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}


SEQ_3_1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def three_to_one(aa):
    if len(aa) == 1:
        return aa

    return SEQ_3_1.get(aa, ".")


def one_to_three(aa):
    if len(aa) == 3:
        return aa

    elif len(aa) == 1:
        return SEQ_1_3.get(aa, "BLK")

    else:
        raise ValueError("Unknown aminoacids notation")


def non_negative(array):
    return array + 1 - np.min(array)


class SeqPos:
    def __init__(self, index, chain, resid, resname=None):
        self.index = index
        self.chain = chain
        self.resid = int(resid)
        self._resname = "gap"
        self._gap = True
        self._blk = False
        self._code = "-"

        if resname is not None:
            self.resname = resname

    def __repr__(self):
        return "SeqPos(index={}, chain='{}', resid={}, resname='{}')".format(self.index, self.chain, self.resid, self.resname)

    def __str__(self):
        return self._code

    def __eq__(self, other):
        return self.index == other.index and self.chain == other.chain and self.resid == other.resid and self.resname == other.resname

    @property
    def code(self):
        return self._code

    @code.setter
    def code(self, value):
        self._code = str(value)

        if value == "-":
            self._resname = "gap"
            self._gap = True

        elif value == ".":
            self._resname = "BLK"
            self._gap = False
            self._blk = True

        else:
            self._resname = one_to_three(value)
            self._gap = False

    @property
    def is_blk(self):
        return self._blk

    @property
    def is_gap(self):
        return self._gap

    @property
    def resname(self):
        return self._resname

    @resname.setter
    def resname(self, value):
        self.code = three_to_one(str(value))


class ChainSeq(list):
    def __init__(self, code, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.code = code

    def __repr__(self):
        if len(self) > 13:
            sequence = "".join(map(str, self[:5]))
            sequence += "..."
            sequence += "".join(map(str, self[-5:]))
        else:
            sequence = "".join(map(str, self))

        return "<ChainSeq '{}' with sequence '{}'>".format(self.code, sequence)

    def __str__(self):
        return "".join(map(str, self))


class ProtSeq:
    def __init__(self, name, input_file=None, sequence=None):
        self.name = name
        self.input_file = input_file
        self._chains = []
        self._chain_ids = {}

        if input_file is not None:
            ext = os.path.splitext(input_file)[1].lower()
            if ext == ".pdb":
                self._init_from_pdb()

            elif ext == ".fasta":
                self._init_from_fasta()

            else:
                raise RuntimeError("Unknown file type")

        elif sequence is not None:
            self._init_from_sequence(sequence)

        else:
            raise RuntimeError("Input file or sequence must be provided")

    def __repr__(self):
        return "<ProtSeq object '{}' with {} chains>".format(self.name, len(self._chains))

    def __getitem__(self, key):
        try:
            return self._chains[key]
        except TypeError:
            pass

        return self._chain_ids[key]

    def __iter__(self):
        return iter(self._chains)

    def _init_from_fasta(self):
        raise NotImplementedError("This method should be implemented")

    def _init_from_pdb(self):
        pdb = mda.Universe(str(self.input_file))

        seq_idx = 0
        for segment in pdb.segments:
            chain = ChainSeq(segment.segid)
            self._chains.append(chain)
            self._chain_ids[segment.segid] = chain

            for res in segment.atoms.residues:
                chain.append(SeqPos(seq_idx, chain.code, res.resid, res.resname))
                seq_idx += 1

    def _init_from_sequence(self, sequence):
        seq_idx = 0
        resid = 1
        for idx, seq in enumerate(sequence.split("/")):
            chain = ChainSeq(string.ascii_uppercase[idx])
            self._chains.append(chain)
            self._chain_ids[chain.code] = chain

            for code in seq:
                chain.append(SeqPos(seq_idx, chain.code, resid, one_to_three(code)))
                seq_idx += 1
                resid += 1

    def flatten(self):
        return tuple(itertools.chain.from_iterable(self._chains))

    def copy(self):
        return copy.deepcopy(self)

    def keys(self):
        return self._chain_ids.keys()

    def to_modeller(self):
        seq = ""
        for chain in self._chains:
            seq += "/" + str(chain)

        return seq[1:]


class ResidueIndex:
    def __init__(self, prot_seq, selection):
        self.protein = prot_seq
        self.sel_text = selection
        self._sel_index = defaultdict(list)
        self.update()

    def apply(self, prot_seq):
        selected = []
        for chain, indexes in self._sel_index.items():
            for index in indexes:
                selected.append(prot_seq[chain][index])

        return selected

    def update(self):
        input_sel = defaultdict(set)

        if not self.sel_text:
            return

        for residue in self.sel_text.split(" "):
            chain, resid = residue.split(":")
            input_sel[chain].add(int(resid))

        for chain, sel_list in input_sel.items():
            for index, residue in enumerate(self.protein[chain]):
                if residue.resid in sel_list:
                    self._sel_index[chain].append(index)
                    sel_list.remove(residue.resid)

            if sel_list:
                not_found = ", ".join(map(str, sorted(sel_list)))
                raise RuntimeError("Residues not found in chain '{}': {}".format(chain, not_found))


class MatrixMutation:
    def __init__(self, aminoacids, weights, inverted=False):
        self._aa = []
        self._weights = {}
        self._set_data(aminoacids, weights, inverted)

    def __getitem__(self, key):
        return self._weights[key]

    def _set_data(self, aminoacids, weights, inverted):
        self._aa = np.array([three_to_one(code) for code in aminoacids])

        weights = np.array(weights, "float64")
        shape = weights.shape

        if len(self._aa) != shape[0] or len(self._aa) != shape[1]:
            raise ValueError("Incompatible shapes between weights matrix and amino acids specification")

        # Store normalized individual weights
        for idx, code in enumerate(self._aa):
            aa_weights = weights[idx, :]
            aa_weights /= aa_weights.sum()
            inverted_weights = np.reciprocal(aa_weights)
            inverted_weights /= inverted_weights.sum()

            if inverted:
                # When non-conservative has higher probabilities
                self._weights[code] = [inverted_weights, aa_weights]

            else:
                # When conservative has higher probabilities
                self._weights[code] = [aa_weights, inverted_weights]

    def mutate(self, aa, conservative=True):
        wild_aa = three_to_one(aa)
        if conservative:
            weights = self._weights[wild_aa][0]
        else:
            weights = self._weights[wild_aa][1]

        return random.choices(self._aa, weights, k=1)[0]


class Blosum62Mutation(MatrixMutation):
    def __init__(self):
        # Include only usual amino acids
        aminoacids, blosum = Blosum()[62]
        super().__init__(aminoacids, blosum, False)


class Blosum:
    __data = None

    def __init__(self):
        if Blosum.__data is None:
            # Load only once
            stream = io.BytesIO()
            stream.write(pkgutil.get_data("gdee", "data/blosum_all.npz"))
            stream.seek(0)
            file = np.load(stream)
            self.__data = {}

            for name in file.files:
                self.__data[name] = file[name]

    def __getitem__(self, key):
        aa = self.__data["aminoacids_{}".format(key)]
        matrix = self.__data["prob_{}".format(key)]
        return aa, matrix
