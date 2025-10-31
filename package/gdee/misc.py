"""
Miscellaneous utility functions and classes for the GDEE platform.

This module provides general-purpose utilities including data containers,
JSON serialization helpers, and file naming functions used throughout
the GDEE workflow.
"""

import re
import json
import numpy as np


def _jsonfy(obj):
    """
    JSON serialization helper for NumPy arrays.
    
    Converts NumPy arrays to Python lists for JSON serialization.
    This function is used as a default serializer in DataContainer.jsonfy().
    
    Args:
        obj: Object to serialize (typically numpy.ndarray)
        
    Returns:
        list: Python list representation of NumPy array
        
    Raises:
        TypeError: If object type is not supported for serialization
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()

    raise TypeError()


class DataContainer(dict):
    """
    Dictionary-like container with attribute-style access.
    
    This class extends dict to provide both dictionary and attribute-style
    access to data. It's used throughout the GDEE pipeline to pass job data
    between components like variant builders, modeling tools, docking evaluators,
    and measurement systems.
    
    The container supports JSON serialization with automatic handling of
    NumPy arrays and other complex data types commonly used in the workflow.
    
    Examples:
        >>> container = DataContainer()
        >>> container.variant_name = "A123V"
        >>> container['energy_scores'] = [1.2, 3.4, 5.6]
        >>> print(container.variant_name)  # "A123V"
        >>> print(container['variant_name'])  # "A123V"
    """
    
    def __setattr__(self, name, value):
        """
        Set attribute value using dictionary assignment.
        
        Args:
            name (str): Attribute name
            value: Value to assign
        """
        self[name] = value

    def __getattr__(self, name):
        """
        Get attribute value from dictionary.
        
        Args:
            name (str): Attribute name
            
        Returns:
            Value associated with the attribute name
            
        Raises:
            AttributeError: If attribute name is not found
        """
        if name in self:
            return self[name]
        raise AttributeError(name)

    def jsonfy(self):
        """
        Serialize container contents to JSON string.
        
        Automatically handles NumPy arrays by converting them to Python lists.
        This method is useful for storing job data in databases or log files.
        
        Returns:
            str: JSON representation of container contents
        """
        return json.dumps(self, default=_jsonfy)


def get_valid_filename(name):
    """
    Generate a valid filename from a string.
    
    Sanitizes input strings to create safe filenames for job directories
    and result files. Removes special characters that might cause issues
    with filesystem operations while preserving readability.
    
    This function is used by variant builders to create directory names
    for individual variant processing jobs.
    
    Args:
        name (str): Input string to convert to filename
        
    Returns:
        str: Sanitized filename string
        
    Examples:
        >>> get_valid_filename("A123V|T456K")
        "A123VT456K"
        >>> get_valid_filename("wild type")
        "wild_type"
    """
    name = str(name).strip().replace(" ", "_")
    return re.sub(r"(?u)[^-\w.]", "", name)
