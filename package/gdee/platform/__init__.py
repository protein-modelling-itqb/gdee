"""
Execution platforms for the GDEE platform.

This module provides factory and platform classes for executing GDEE workflows
on different computational environments, including single-machine and
distributed MPI-based execution.
"""

from .simple_platform import SimplePlatform
from .mpi_platform import MPIPlatform

__all__ = ["PlatformFactory"]


class PlatformFactory:
    """
    Factory for creating execution platforms.

    This factory creates platform objects based on configuration parameters
    to execute GDEE workflows on different computational
    infrastructures.

    Attributes:
        parameters (dict): Platform configuration including execution type
                          and resource allocation
        pipeline (Pipeline): Configured pipeline for variant processing
    """

    def __init__(self):
        """Initialize factory with empty configuration."""
        self.parameters = {}
        self.pipeline = None

    def make(self):
        """
        Create and configure an execution platform.

        Creates either a SimplePlatform for single-machine execution or
        an MPIPlatform for distributed execution based on configuration.

        Returns:
            Platform: Configured execution platform (SimplePlatform or MPIPlatform)

        Raises:
            RuntimeError: If platform type is not recognized
        """
        name = self.parameters["name"]

        if name == "simple":
            return SimplePlatform(self.parameters, self.pipeline)

        elif name == "mpi":
            return MPIPlatform(self.parameters, self.pipeline)

        else:
            raise RuntimeError("Platform '{}' does not exist.".format(name))
