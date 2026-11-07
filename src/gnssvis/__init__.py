"""GNSSVis - Radio Telescope Observation Analysis and Visualization

A Python library for radio astronomy data analysis, focusing on detection
and visualization of GNSS satellites and other celestial objects."""

__version__ = "0.1.0"
__author__ = "Erik Limeira"

# Main imports
from .core import GNSSVis, DataAnalyzer
from .io import IOManager
from .sky import SkyAnalyzer
from .plot import Visualizer

# Simplified API
__all__ = [
    'GNSSVis',
    'DataAnalyzer', 
    'IOManager',
    'SkyAnalyzer',
    'Visualizer'
]
 