"""
Distance measurement metrics for molecular analysis.

This module provides distance calculation functions for measuring geometric
relationships between protein and ligand atoms in docking poses.
"""

import numpy as np


class EuclideanDistance:
    """
    Euclidean distance calculator for protein-ligand measurements.
    
    This class computes center-of-mass distances between selected protein
    and ligand atom groups. It's commonly used for measuring binding site
    proximity and pose filtering in molecular docking analysis.
    """
    
    @staticmethod
    def name():
        """
        Get metric identifier for database storage.
        
        Returns:
            str: Metric name used in measurement specifications
        """
        return "distance"

    def compute(self, mat_pos1, mat_pos2):
        """
        Calculate Euclidean distance between two atom groups.
        
        Computes the distance between the centers of mass of two sets
        of atomic coordinates. This is useful for measuring binding site
        proximity and validating docking poses.
        
        Args:
            mat_pos1 (numpy.ndarray): Positions of first atom group (protein)
                                     Shape: (n_atoms, 3)
            mat_pos2 (numpy.ndarray): Positions of second atom group (ligand)
                                     Shape: (m_atoms, 3)
                                     
        Returns:
            float: Euclidean distance in Angstroms between group centroids
            
        Examples:
            >>> metric = EuclideanDistance()
            >>> protein_coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
            >>> ligand_coords = np.array([[7.0, 8.0, 9.0]])
            >>> distance = metric.compute(protein_coords, ligand_coords)
        """
        protein_center = mat_pos1.mean(axis=0)
        ligand_center = mat_pos2.mean(axis=0)
        return np.linalg.norm(protein_center - ligand_center)
