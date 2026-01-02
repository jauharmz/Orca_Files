"""Energy parser for ORCA output."""

from typing import Optional

from ..core.base_parser import BaseParser
from ..core.data_models import EnergyData
from . import regex_patterns as rx


class EnergyParser(BaseParser):
    """Parse energy data from ORCA output."""
    
    def parse(self) -> EnergyData:
        """Parse energy data."""
        self.logger.debug("Parsing energy...")
        
        data = EnergyData()
        
        # Parse Gibbs energy
        data.gibbs_Eh = self._parse_gibbs()
        if data.gibbs_Eh:
            self._log_found(f"Gibbs energy: {data.gibbs_Eh} Eh")
        
        # Parse single point energy
        data.single_point_Eh = self._parse_single_point()
        if data.single_point_Eh:
            self._log_found(f"Single point energy: {data.single_point_Eh} Eh")
        
        return data
    
    def _parse_gibbs(self) -> Optional[float]:
        """Parse Gibbs free energy."""
        match = rx.GIBBS_ENERGY.search(self.text)
        if match:
            return self._parse_float(match.group(1))
        self._log_not_found("Gibbs energy")
        return None
    
    def _parse_single_point(self) -> Optional[float]:
        """Parse final single point energy."""
        matches = rx.SINGLE_POINT_ENERGY.findall(self.text)
        if matches:
            # Return last match (optimized)
            return self._parse_float(matches[-1])
        self._log_not_found("Single point energy")
        return None
