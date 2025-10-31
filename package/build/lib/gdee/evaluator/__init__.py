"""
"""


from .vina import VinaDocking, VinardoDocking


__all__ = ["EvaluatorFactory"]


class EvaluatorFactory:
    def __init__(self):
        self.parameters = {}

    def make(self):
        sizes = self.parameters["box_size"]
        if sizes is None or len(sizes) != 3:
            raise RuntimeError("Invalid evaluator box sizes.")

        center = self.parameters["box_center"]
        if center is None or len(center) != 3:
            raise RuntimeError("Invalid evaluator box center.")

        name = self.parameters["name"]
        if name == "vina":
            return VinaDocking(self.parameters)

        elif name == "vinardo":
            return VinardoDocking(self.parameters)

        else:
            raise RuntimeError("Evaluator '{}' does not exists.".format(name))
