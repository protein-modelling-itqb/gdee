"""
PDBQT file parser and PDB converter for molecular docking results.

This module provides classes for parsing AutoDock PDBQT files containing
docking poses and converting them to standard PDB format for analysis
and visualization.
"""


import numpy as np


class PDBQT:
    """
    Parser and converter for AutoDock PDBQT format files.
    
    This class reads PDBQT files produced by AutoDock Vina/Vinardo
    containing multiple docking poses with binding energies. It extracts
    atomic coordinates, binding energies, and converts results to PDB format.
    
    Attributes:
        atom_data (list): Atom record lines (first 27 characters)
        atom_qt (list): AutoDock atom type information (after column 70)
        models (list): List of DockingModel objects with poses and energies
    """
    
    def __init__(self, filename):
        """
        Initialize PDBQT parser and process file.
        
        Args:
            filename (str): Path to PDBQT file from docking calculation
        """
        self.atom_data = []
        self.atom_qt = []
        self.models = []
        self.parse(filename)

    def size(self):
        """
        Get number of docking poses.
        
        Returns:
            int: Number of poses in the PDBQT file
        """
        return len(self.models)

    def __getitem__(self, key):
        """Access docking pose by index."""
        return self.models[key]

    def __iter__(self):
        """Iterate over docking poses."""
        return iter(self.models)

    def parse(self, filename):
        """
        Parse PDBQT file and extract docking poses.
        
        Processes PDBQT format with support for multiple MODEL/ENDMDL blocks,
        extracts binding energies from various comment formats (USER, REMARK),
        and builds coordinate arrays for each pose.
        
        Args:
            filename (str): Path to PDBQT file to parse
        """
        coords = []
        charges = []
        types = []
        energy = 0
        has_models = False
        first_pass = True

        with open(filename) as file:
            for line in file:
                trimmed = line.strip()

                if trimmed.startswith("MODEL"):
                    # Start of new docking pose
                    has_models = True
                    coords.clear()
                    charges.clear()
                    types.clear()

                elif trimmed.startswith("USER") and \
                    "Estimated Free Energy of Binding" in trimmed:
                    # Extract binding energy from USER record
                    columns = trimmed.split()
                    energy = float(columns[7])

                elif trimmed.startswith("REMARK"):
                    # Extract binding energy from REMARK records
                    columns = trimmed.split()
                    if "VINA" in columns and "RESULT:" in columns:
                        energy = float(columns[3])

                    elif "minimizedAffinity" in columns:
                        energy = float(columns[2])

                elif trimmed.startswith("ATOM") or trimmed.startswith("HETATM"):
                    # Process atom records
                    if first_pass:
                        # Store atom metadata on first pass
                        self.atom_data.append(trimmed[:27])  # PDB atom record prefix
                        self.atom_qt.append(trimmed[70:])    # AutoDock type info

                    # Extract coordinates for current pose
                    x = float(trimmed[30:38])
                    y = float(trimmed[38:46])
                    z = float(trimmed[46:54])
                    coords.append((x, y, z))

                elif trimmed.startswith("ENDMDL"):
                    # End of current pose - create model
                    first_pass = False
                    if coords:
                        self.models.append(DockingModel(coords, energy))

            # Handle single-pose files without MODEL/ENDMDL blocks
            if not has_models:
                if coords:
                    self.models.append(DockingModel(coords, energy))

    def write_pdb(self, filename):
        """
        Convert PDBQT poses to multi-model PDB format.
        
        Writes all docking poses to a single PDB file with MODEL/ENDMDL
        blocks. Binding energies are stored in REMARK records and
        B-factor columns for visualization.
        
        Args:
            filename (str): Output PDB file path
        """
        with open(filename, "w") as fd:
            for model_idx, model in enumerate(self.models):
                # Write model header with binding energy
                fd.write("MODEL {:d}\nREMARK ENERGY: {:.2f}\n".format(
                    model_idx, model.energy))

                # Write atom coordinates with binding energy in B-factor column
                for atom_idx in range(len(self.atom_data)):
                    fd.write("{:<30s}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
                        self.atom_data[atom_idx],        # Atom record prefix
                        model.coords[atom_idx, 0],       # X coordinate
                        model.coords[atom_idx, 1],       # Y coordinate
                        model.coords[atom_idx, 2],       # Z coordinate
                        0,                               # Occupancy (unused)
                        model.energy,                    # B-factor = binding energy
                    ))

                fd.write("ENDMDL\n")


class DockingModel:
    """
    Individual docking pose with coordinates and binding energy.
    
    This class represents a single docking pose including atomic coordinates
    and associated binding energy. It provides methods for pose comparison
    and geometric analysis.
    
    Attributes:
        coords (numpy.ndarray): Atomic coordinates array (N_atoms x 3)
        energy (float): Binding energy in kcal/mol (lower is better)
    """
    
    def __init__(self, coords, energy):
        """
        Initialize docking pose.
        
        Args:
            coords (list): List of (x, y, z) coordinate tuples
            energy (float): Binding energy in kcal/mol
        """
        self.coords = np.array(coords, "float64")
        self.energy = energy

    def rmsd(self, other):
        """
        Calculate root mean square deviation between poses.
        
        Computes RMSD between atomic positions for pose comparison
        and clustering analysis.
        
        Args:
            other (DockingModel): Reference pose for comparison
            
        Returns:
            float: RMSD value in Angstroms
            
        Raises:
            ValueError: If poses have different numbers of atoms
        """
        if self.coords.shape != other.coords.shape:
            raise ValueError("Models must have the same number of atoms")

        return np.sqrt(np.sum((self.coords - other.coords) ** 2) / self.coords.shape[0])

    def centroid(self):
        """
        Calculate geometric center of the pose.
        
        Returns:
            numpy.ndarray: Centroid coordinates [x, y, z]
        """
        return self.coords.mean(0)

    def box(self):
        """
        Calculate bounding box dimensions.
        
        Returns:
            numpy.ndarray: Box dimensions [width, height, depth] in Angstroms
        """
        return self.coords.max(0) - self.coords.min(0)
