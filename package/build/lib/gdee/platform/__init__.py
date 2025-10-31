"""
"""


from .simple_platform import SimplePlatform
from .mpi_platform import MPIPlatform


__all__ = ["PlatformFactory"]


class PlatformFactory:
    def __init__(self):
        self.parameters = {}
        self.pipeline = None

    def make(self):
        name = self.parameters["name"]

        if name == "simple":
            return SimplePlatform(self.parameters, self.pipeline)

        elif name == "mpi":
            return MPIPlatform(self.parameters, self.pipeline)

        else:
            raise RuntimeError("Platform '{}' does not exists.".format(name))
