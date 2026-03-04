"""Distance measurement calculations."""

import numpy as np


class EuclideanDistance:
    """Computes Euclidean distance between atom selections."""

    @staticmethod
    def name():
        """Get measurement metric name.
        
        Returns:
            str: Metric identifier "distance"
        """
        return "distance"

    def compute(self, mat_pos1, mat_pos2):
        """Compute Euclidean distance between two atom selections.
        
        Args:
            mat_pos1: N x 3 array of first atom positions
            mat_pos2: M x 3 array of second atom positions
            
        Returns:
            float: Distance between centers of mass in Angstroms
        """
        return np.linalg.norm(mat_pos1.mean(0) - mat_pos2.mean(0))
