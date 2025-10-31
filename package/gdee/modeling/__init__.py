"""
3D protein structure modeling for the GDEE platform.

This module provides factory classes and builders for generating 3D protein structures
from sequence variants using MODELLER, and for assessing model quality using various
scoring functions including VoroMQA and Normalized DOPE.
"""


from .modeller_builder import ModellerBuilder
from .quality import ModelQualityChecker

__all__ = ["ModelBuilderFactory", "ModelQualityBuilderFactory"]


class ModelBuilderFactory:
    """
    Factory for creating 3D structure modeling components.

    This factory creates modeling tools based on configuration parameters.
    Currently supports MODELLER for homology modeling with mutation-specific
    optimization strategies.

    Attributes:
        parameters (dict): Modeling configuration including method name,
                          optimization settings, and template PDB file
    """

    def __init__(self):
        """Initialize factory with empty parameters."""
        self.parameters = {}

    def make(self):
        """
        Create and configure a structure modeling component.

        Validates template PDB file and creates the appropriate modeling
        tool based on the specified method.

        Returns:
            ModellerBuilder: Configured 3D structure modeling component

        Raises:
            RuntimeError: If PDB file is invalid or modeling method is unknown
        """
        if self.parameters["pdb_file"] is None:
            raise RuntimeError("Invalid template PDB file.")

        name = self.parameters["name"]
        if name == "modeller":
            return ModellerBuilder(self.parameters)

        # Future expansion for other modeling tools
        # if name == "rosetta":
        #     return RosettaBuilder(self.parameters)

        else:
            raise RuntimeError("Modeler '{}' does not exist.".format(name))


class ModelQualityBuilderFactory:
    """
    Factory for creating model quality assessment components.

    This factory creates quality checkers that evaluate 3D models using
    various scoring functions and apply cutoff thresholds for model
    acceptance or rejection.

    Attributes:
        parameters (dict): Quality assessment configuration including
                          scoring methods, cutoff values, and program paths
    """

    def __init__(self):
        """Initialize factory with empty parameters."""
        self.parameters = {}

    def make(self):
        """
        Create and configure a model quality checker.

        Sets up scoring functions (VoroMQA, Normalized DOPE) and applies
        cutoff thresholds based on configuration parameters.

        Returns:
            ModelQualityChecker: Configured quality assessment component
        """
        checker = ModelQualityChecker()

        # Configure Normalized DOPE scoring (computed by MODELLER)
        if self.parameters["norm_dope"] is not None:
            # Normalized DOPE is computed by MODELLER during model building
            checker.add_lower_cutoff("norm_dope", self.parameters["norm_dope"])

        # Configure VoroMQA scoring (external program)
        if self.parameters["voromqa"] is not None:
            checker.enable_scorer("voromqa", self.parameters["programs"]["voromqa"])
            checker.add_higher_cutoff("voromqa", self.parameters["voromqa"])

        return checker
