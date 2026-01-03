"""
ORCA Parser Package - Modular parsers refactored from orca_praser.py

Usage:
    from src.parser import ParserFactory, BatchParser
    
    # Single file
    factory = ParserFactory()
    result = factory.parse("molecule.out")
    
    # Batch files
    batch = BatchParser("data/**/*.out")
    df = batch.parse_all()
"""

from .factory import ParserFactory
from .batch import BatchParser
from .geometry import GeometryParser
from .energy import EnergyParser
from .orbitals import OrbitalParser
from .spectroscopy import SpectroscopyParser
from .tddft import TDDFTParser
from .spectrum_file import SpectrumFileParser

__all__ = [
    "ParserFactory",
    "BatchParser",
    "GeometryParser",
    "EnergyParser",
    "OrbitalParser",
    "SpectroscopyParser",
    "TDDFTParser",
    "SpectrumFileParser",
]
