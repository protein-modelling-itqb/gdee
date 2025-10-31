"""
Protein sequence handling and mutation management for the GDEE platform.

This module provides classes for representing protein sequences, managing mutations,
and handling amino acid conversions. It supports loading sequences from PDB files
or FASTA format, and provides mutation matrices for variant generation.
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

# Amino acid conversion dictionaries
SEQ_1_3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}

SEQ_3_1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def three_to_one(aa):
    """
    Convert three-letter amino acid code to single-letter code.
    
    Args:
        aa (str): Three-letter amino acid code (e.g., "ALA") or single-letter code
        
    Returns:
        str: Single-letter amino acid code (e.g., "A") or "." for unknown
    """
    if len(aa) == 1:
        return aa

    return SEQ_3_1.get(aa, ".")


def one_to_three(aa):
    """
    Convert single-letter amino acid code to three-letter code.
    
    Args:
        aa (str): Single-letter amino acid code (e.g., "A") or three-letter code
        
    Returns:
        str: Three-letter amino acid code (e.g., "ALA") or "BLK" for unknown
        
    Raises:
        ValueError: If amino acid notation is not 1 or 3 characters
    """
    if len(aa) == 3:
        return aa

    elif len(aa) == 1:
        return SEQ_1_3.get(aa, "BLK")

    else:
        raise ValueError("Unknown aminoacids notation")


def non_negative(array):
    """
    Convert array to non-negative values by shifting minimum to 1.
    
    Args:
        array (numpy.ndarray): Input array
        
    Returns:
        numpy.ndarray: Array with all non-negative values
    """
    return array + 1 - np.min(array)


class SeqPos:
    """
    Represents a single position in a protein sequence.
    
    This class stores information about amino acid residues including their
    position, chain, residue number, and amino acid type. It handles gaps
    and unknown residues in sequences.
    
    Attributes:
        index (int): Sequential position in flattened sequence
        chain (str): Chain identifier
        resid (int): Residue number from PDB
        code (str): Single-letter amino acid code
        resname (str): Three-letter amino acid code
        is_gap (bool): Whether position is a gap ("-")
        is_blk (bool): Whether position is unknown (".")
    """
    
    def __init__(self, index, chain, resid, resname=None):
        """
        Initialize sequence position.
        
        Args:
            index (int): Sequential position index
            chain (str): Chain identifier
            resid (int): Residue number
            resname (str, optional): Three-letter amino acid code
        """
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
        """Single-letter amino acid code."""
        return self._code

    @code.setter
    def code(self, value):
        """Set amino acid code and update derived properties."""
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
        """Whether position represents unknown amino acid."""
        return self._blk

    @property
    def is_gap(self):
        """Whether position represents a gap."""
        return self._gap

    @property
    def resname(self):
        """Three-letter amino acid code."""
        return self._resname

    @resname.setter
    def resname(self, value):
        """Set three-letter code and update single-letter code."""
        self.code = three_to_one(str(value))


class ChainSeq(list):
    """
    Represents a single protein chain as a list of sequence positions.
    
    Extends list to store SeqPos objects with an associated chain identifier.
    Provides string representation of the amino acid sequence.
    
    Attributes:
        code (str): Chain identifier (e.g., "A", "B")
    """
    
    def __init__(self, code, *args, **kwargs):
        """
        Initialize protein chain.
        
        Args:
            code (str): Chain identifier
            *args, **kwargs: Additional arguments for list initialization
        """
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
        """Return amino acid sequence as string."""
        return "".join(map(str, self))


class ProtSeq:
    """
    Represents a complete protein sequence with multiple chains.
    
    This class manages protein sequences loaded from PDB files, FASTA files,
    or provided as strings. It supports multi-chain proteins and provides
    methods for sequence manipulation and format conversion.
    
    Attributes:
        name (str): Protein name identifier
        input_file (str): Source file path
        _chains (list): List of ChainSeq objects
        _chain_ids (dict): Mapping of chain IDs to ChainSeq objects
    """
    
    def __init__(self, name, input_file=None, sequence=None):
        """
        Initialize protein sequence.
        
        Args:
            name (str): Protein identifier
            input_file (str, optional): Path to PDB or FASTA file
            sequence (str, optional): Amino acid sequence string
            
        Raises:
            RuntimeError: If neither input_file nor sequence is provided,
                         or if file type is not supported
        """
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
        """Access chains by index or chain ID."""
        try:
            return self._chains[key]
        except TypeError:
            pass

        return self._chain_ids[key]

    def __iter__(self):
        """Iterate over chains."""
        return iter(self._chains)

    def _init_from_fasta(self):
        """Initialize from FASTA file (not implemented)."""
        raise NotImplementedError("This method should be implemented")

    def _init_from_pdb(self):
        """
        Initialize sequence from PDB file.
        
        Uses MDAnalysis to parse PDB structure and extract sequence
        information including chain IDs and residue numbers.
        """
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
        """
        Initialize from amino acid sequence string.
        
        Multi-chain sequences are separated by "/" characters.
        Chains are labeled with uppercase letters (A, B, C, ...).
        
        Args:
            sequence (str): Amino acid sequence with "/" separating chains
        """
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
        """
        Return all sequence positions as a flat tuple.
        
        Returns:
            tuple: All SeqPos objects from all chains in sequential order
        """
        return tuple(itertools.chain.from_iterable(self._chains))

    def copy(self):
        """
        Create a deep copy of the protein sequence.
        
        Returns:
            ProtSeq: Independent copy of the protein sequence
        """
        return copy.deepcopy(self)

    def keys(self):
        """
        Get chain identifiers.
        
        Returns:
            dict_keys: Chain IDs present in the protein
        """
        return self._chain_ids.keys()

    def to_modeller(self):
        """
        Convert sequence to MODELLER format.
        
        Returns slash-separated sequence string suitable for MODELLER
        alignment files and template specification.
        
        Returns:
            str: Sequence in MODELLER format (e.g., "ACGT/DEFG")
        """
        seq = ""
        for chain in self._chains:
            seq += "/" + str(chain)

        return seq[1:]


class ResidueIndex:
    """
    Manages residue selection and indexing for protein sequences.
    
    This class parses residue selection strings and maps them to sequence
    positions. It's used by variant builders to identify mutable positions
    and fixed residues during optimization.
    
    Attributes:
        protein (ProtSeq): Target protein sequence
        sel_text (str): Selection string (e.g., "A:123 A:456 B:789")
        _sel_index (dict): Mapping of chain IDs to residue indices
    """
    
    def __init__(self, prot_seq, selection):
        """
        Initialize residue selection.
        
        Args:
            prot_seq (ProtSeq): Protein sequence to select from
            selection (str): Space-separated residue identifiers (chain:resid)
        """
        self.protein = prot_seq
        self.sel_text = selection
        self._sel_index = defaultdict(list)
        self.update()

    def apply(self, prot_seq):
        """
        Apply selection to a protein sequence.
        
        Args:
            prot_seq (ProtSeq): Protein sequence to select residues from
            
        Returns:
            list: Selected SeqPos objects
        """
        selected = []
        for chain, indexes in self._sel_index.items():
            for index in indexes:
                selected.append(prot_seq[chain][index])

        return selected

    def update(self):
        """
        Parse selection string and update internal indices.
        
        Parses selection strings like "A:123 B:456" and maps them to
        sequence positions within each chain.
        
        Raises:
            RuntimeError: If specified residues are not found in the protein
        """
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
    """
    Manages amino acid mutations based on substitution matrices.
    
    This class provides weighted mutation selection using substitution matrices
    like BLOSUM. It supports both conservative and non-conservative mutation
    strategies based on matrix probabilities.
    
    Attributes:
        _aa (numpy.ndarray): Array of amino acid single-letter codes
        _weights (dict): Normalized probability weights for each amino acid
    """
    
    def __init__(self, aminoacids, weights, inverted=False):
        """
        Initialize mutation matrix.
        
        Args:
            aminoacids (array-like): Amino acid codes
            weights (array-like): Substitution matrix weights
            inverted (bool): Whether to invert weight interpretation
        """
        self._aa = []
        self._weights = {}
        self._set_data(aminoacids, weights, inverted)

    def __getitem__(self, key):
        """Get mutation weights for amino acid."""
        return self._weights[key]

    def _set_data(self, aminoacids, weights, inverted):
        """
        Process and normalize substitution matrix data.
        
        Args:
            aminoacids (array-like): Amino acid identifiers
            weights (array-like): 2D substitution matrix
            inverted (bool): Weight interpretation mode
            
        Raises:
            ValueError: If matrix dimensions don't match amino acid count
        """
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
        """
        Generate a mutation for the given amino acid.
        
        Args:
            aa (str): Wildtype amino acid code
            conservative (bool): Use conservative mutation weights
            
        Returns:
            str: Mutated amino acid single-letter code
        """
        wild_aa = three_to_one(aa)
        if conservative:
            weights = self._weights[wild_aa][0]
        else:
            weights = self._weights[wild_aa][1]

        return random.choices(self._aa, weights, k=1)[0]


class Blosum62Mutation(MatrixMutation):
    """
    BLOSUM62-based mutation matrix for conservative amino acid substitutions.
    
    This class implements the BLOSUM62 substitution matrix commonly used
    in protein sequence analysis. It provides biologically reasonable
    mutation probabilities based on observed substitution frequencies.
    """
    
    def __init__(self):
        """Initialize with BLOSUM62 matrix data."""
        # Include only usual amino acids
        aminoacids, blosum = Blosum()[62]
        super().__init__(aminoacids, blosum, False)


class Blosum:
    """
    BLOSUM substitution matrix data loader.
    
    This class provides access to BLOSUM matrices stored in compressed
    NumPy format. It implements singleton pattern to load matrix data
    only once during program execution.
    
    Class Attributes:
        __data (dict): Cached matrix data (shared across instances)
    """
    
    __data = None

    def __init__(self):
        """Initialize and load BLOSUM matrix data if not already cached."""
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
        """
        Get BLOSUM matrix by number.
        
        Args:
            key (int): BLOSUM matrix number (e.g., 62 for BLOSUM62)
            
        Returns:
            tuple: (amino_acids, probability_matrix)
        """
        aa = self.__data["aminoacids_{}".format(key)]
        matrix = self.__data["prob_{}".format(key)]
        return aa, matrix
