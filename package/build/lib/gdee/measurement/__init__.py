"""
"""


from .measurer import Measurer
from .distance import EuclideanDistance


__all__ = ["MeasurerFactory"]


class MeasurerFactory:
    def __init__(self):
        self.ligand = None
        self.names = set()
        self.metrics = {
            EuclideanDistance.name(): EuclideanDistance,
        }

    def make(self):
        measurer = Measurer(self.ligand.name)

        for name, metric, prot_sel, lig_sel in self.ligand.measurements:
            # Sanity check
            if name in self.names:
                raise RuntimeError("Duplicate metric '{}'".format(name))
            self.names.add(name)

            if metric in self.metrics:
                measurer.add(self.metrics[metric](), name, prot_sel, lig_sel)

            else:
                raise RuntimeError("Measurement category '{}' does not exists.".format(name))
        
        return measurer
