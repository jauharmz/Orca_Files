"""Parser factory for orchestrating all parsers."""

from typing import Optional
from pathlib import Path

from ..core.data_models import ParseResult
from ..core.exceptions import ORCAFileError
from ..logger import get_logger
from .geometry import GeometryParser
from .energy import EnergyParser
from .orbitals import OrbitalParser
from .spectroscopy import SpectroscopyParser
from .tddft import TDDFTParser


class ParserFactory:
    """Factory for parsing ORCA output files."""
    
    def __init__(self):
        self.logger = get_logger("ParserFactory")
    
    def parse(self, filepath: str) -> ParseResult:
        """
        Parse an ORCA output file.
        
        Args:
            filepath: Path to .out file
            
        Returns:
            ParseResult with all extracted data
        """
        self.logger.info(f"Parsing: {filepath}")
        
        # Load file
        path = Path(filepath)
        if not path.exists():
            raise ORCAFileError(f"File not found: {filepath}")
        
        text = path.read_text(encoding='utf-8', errors='ignore')
        
        # Create result
        result = ParseResult(filename=str(path))
        
        # Run all parsers
        result.geometry = GeometryParser(text, filepath).parse()
        result.energy = EnergyParser(text, filepath).parse()
        result.orbitals = OrbitalParser(text, filepath).parse()
        result.spectra = SpectroscopyParser(text, filepath).parse()
        result.tddft = TDDFTParser(text, filepath).parse()
        
        # Detect calculation type
        result.is_optimization = "GEOMETRY OPTIMIZATION" in text
        result.has_tddft = result.tddft.states is not None
        
        # Detect optimized state
        if "OPTIMIZE THE S0 STATE" in text.upper():
            result.optimized_state = "S0"
        elif "OPTIMIZE THE S1 STATE" in text.upper():
            result.optimized_state = "S1"
        elif "OPTIMIZE THE T1 STATE" in text.upper():
            result.optimized_state = "T1"
        
        self.logger.info(f"Parsed: {path.name} (OPT={result.is_optimization}, TDDFT={result.has_tddft})")
        
        return result
    
    def parse_text(self, text: str, filename: str = "unknown") -> ParseResult:
        """
        Parse ORCA output from text string.
        
        Args:
            text: ORCA output text
            filename: Optional filename for reference
            
        Returns:
            ParseResult with all extracted data
        """
        result = ParseResult(filename=filename)
        
        result.geometry = GeometryParser(text, filename).parse()
        result.energy = EnergyParser(text, filename).parse()
        result.orbitals = OrbitalParser(text, filename).parse()
        result.spectra = SpectroscopyParser(text, filename).parse()
        result.tddft = TDDFTParser(text, filename).parse()
        
        result.is_optimization = "GEOMETRY OPTIMIZATION" in text
        result.has_tddft = result.tddft.states is not None
        
        return result
