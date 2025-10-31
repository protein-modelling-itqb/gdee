"""
Molecular measurement system for GDEE platform.

This module provides factory classes and measurement tools for calculating
geometric properties and distances between protein-ligand complexes from
molecular docking results.
"""

from .measurer import Measurer
from .distance import EuclideanDistance


__all__ = ["MeasurerFactory"]


class MeasurerFactory:
    """
    Factory for creating molecular measurement systems.
    
    This factory creates measurement components that calculate geometric
    properties between protein and ligand atoms during pose analysis.
    It supports extensible metric types and validates measurement configurations.
    
    Attributes:
        ligand (Ligand): Target ligand with measurement specifications
        names (set): Cache of measurement names to prevent duplicates
        metrics (dict): Registry of available measurement metrics
    """
    
    def __init__(self):
        """Initialize factory with default measurement metrics."""
        self.ligand = None
        self.names = set()
        self.metrics = {
            EuclideanDistance.name(): EuclideanDistance,
            # Future expansion for other metrics:
            # "angle": AngleMeasurement,
            # "dihedral": DihedralMeasurement,
        }

    def make(self):
        """
        Create and configure a molecular measurement system.
        
        Processes all measurement specifications from the ligand configuration
        and creates corresponding measurement tasks with proper validation.
        
        Returns:
            Measurer: Configured measurement system for ligand analysis
            
        Raises:
            RuntimeError: If duplicate measurement names exist or 
                         if measurement category is not supported
        """
        measurer = Measurer(self.ligand.name)

        for name, metric, prot_sel, lig_sel in self.ligand.measurements:
            # Validate measurement name uniqueness
            if name in self.names:
                raise RuntimeError("Duplicate metric '{}'".format(name))
            self.names.add(name)

            # Create measurement task if metric is supported
            if metric in self.metrics:
                metric_instance = self.metrics[metric]()
                measurer.add(metric_instance, name, prot_sel, lig_sel)

            else:
                raise RuntimeError("Measurement category '{}' does not exist.".format(metric))
        
        return measurer
