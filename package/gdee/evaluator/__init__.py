"""
Molecular docking evaluation system for the GDEE platform.

This module provides factory classes and docking tools for performing
protein-ligand molecular docking using AutoDock Vina and Vinardo functions.
It supports automated pose generation, energy scoring, and result analysis.
"""

from .vina import VinaDocking, VinardoDocking

__all__ = ["EvaluatorFactory"]


class EvaluatorFactory:
    """
    Factory for creating molecular docking evaluation components.
    
    This factory creates docking evaluators based on configuration parameters.
    It validates docking box dimensions and center coordinates before
    instantiating the appropriate docking engine.
    
    Attributes:
        parameters (dict): Docking configuration including engine type,
                          search box parameters, and program paths
    """
    
    def __init__(self):
        """Initialize factory with empty parameters."""
        self.parameters = {}

    def make(self):
        """
        Create and configure a molecular docking evaluator.
        
        Validates search box geometry and creates the appropriate docking
        engine (Vina or Vinardo) based on configuration parameters.
        
        Returns:
            BaseVina: Configured docking evaluator (VinaDocking or VinardoDocking)
            
        Raises:
            RuntimeError: If box dimensions are invalid, center coordinates
                         are missing, or docking engine is not recognized
        """
        # Validate docking box dimensions (must be 3D coordinates)
        sizes = self.parameters["box_size"]
        if sizes is None or len(sizes) != 3:
            raise RuntimeError("Invalid evaluator box sizes.")

        # Validate docking box center coordinates
        center = self.parameters["box_center"]
        if center is None or len(center) != 3:
            raise RuntimeError("Invalid evaluator box center.")

        # Create appropriate docking engine
        name = self.parameters["name"]
        if name == "vina":
            return VinaDocking(self.parameters)

        elif name == "vinardo":
            return VinardoDocking(self.parameters)

        else:
            raise RuntimeError("Evaluator '{}' does not exist.".format(name))
