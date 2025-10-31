"""
"""


import numpy as np


class EuclideanDistance:
    @staticmethod
    def name():
        return "distance"

    def compute(self, mat_pos1, mat_pos2):
        return np.linalg.norm(mat_pos1.mean(0) - mat_pos2.mean(0))
