"""TD-DFT parser for ORCA output."""

from typing import Optional, Dict
import pandas as pd
import re

from ..core.base_parser import BaseParser
from ..core.data_models import TDDFTData
from . import regex_patterns as rx


class TDDFTParser(BaseParser):
    """Parse TD-DFT excited state data."""
    
    def parse(self) -> TDDFTData:
        """Parse TD-DFT data."""
        self.logger.debug("Parsing TD-DFT...")
        
        data = TDDFTData()
        
        # Parse excited states
        data.states = self._parse_states()
        
        # Parse electric dipole spectrum
        data.electric_dipole = self._parse_electric_dipole()
        
        # Parse velocity dipole spectrum
        data.velocity_dipole = self._parse_velocity_dipole()
        
        return data
    
    def _parse_states(self) -> Optional[pd.DataFrame]:
        """Parse TD-DFT excited states with transitions."""
        states = []
        
        # Find HOMO index from transitions
        from_orbs = set()
        to_orbs = set()
        
        for match in rx.TDDFT_STATE.finditer(self.text):
            transitions_block = match.group(5)
            for trans in rx.TDDFT_TRANSITION.finditer(transitions_block):
                from_orbs.add(int(trans.group(1)))
                to_orbs.add(int(trans.group(2)))
        
        homo_idx = max(from_orbs) if from_orbs else None
        lumo_idx = min(to_orbs) if to_orbs else None
        
        # Parse states
        for match in rx.TDDFT_STATE.finditer(self.text):
            state_num = int(match.group(1))
            energy_au = self._parse_float(match.group(2))
            energy_ev = self._parse_float(match.group(3))
            energy_cm1 = self._parse_float(match.group(4))
            transitions_block = match.group(5)
            
            # Parse transitions
            for trans in rx.TDDFT_TRANSITION.finditer(transitions_block):
                from_orb = int(trans.group(1))
                to_orb = int(trans.group(2))
                weight = self._parse_float(trans.group(3))
                coeff = self._parse_float(trans.group(4))
                
                # Convert to HOMO/LUMO notation
                homo = from_orb - homo_idx if homo_idx and from_orb <= homo_idx else None
                lumo = to_orb - lumo_idx if lumo_idx and to_orb >= lumo_idx else None
                
                states.append({
                    "state": state_num,
                    "energy_au": energy_au,
                    "energy_eV": energy_ev,
                    "energy_cm-1": energy_cm1,
                    "from_orb": from_orb,
                    "to_orb": to_orb,
                    "homo": homo,
                    "lumo": lumo,
                    "weight": weight,
                    "coeff": coeff
                })
        
        if states:
            self._log_found("TD-DFT states", len(set(s["state"] for s in states)))
            return pd.DataFrame(states)
        
        self._log_not_found("TD-DFT states")
        return None
    
    def _parse_dipole_spectrum(self, pattern: re.Pattern) -> Dict[str, pd.DataFrame]:
        """Parse absorption spectrum (electric or velocity dipole)."""
        result = {"abs": [], "soc": []}
        
        match = pattern.search(self.text)
        if not match:
            return {"abs": pd.DataFrame(), "soc": pd.DataFrame()}
        
        for line in match.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) < 9:
                continue
            
            try:
                transition = parts[0]
                if "->" in " ".join(parts[:3]):
                    transition = " ".join(parts[:3]).replace(" ", "")
                    parts = parts[2:]
                
                result["abs"].append({
                    "transition": transition,
                    "energy_eV": self._parse_float(parts[1]),
                    "energy_cm-1": self._parse_float(parts[2]),
                    "wavelength_nm": self._parse_float(parts[3]),
                    "fosc": self._parse_float(parts[4]),
                    "D2": self._parse_float(parts[5]),
                    "dx": self._parse_float(parts[6]),
                    "dy": self._parse_float(parts[7]),
                    "dz": self._parse_float(parts[8])
                })
            except (ValueError, IndexError):
                continue
        
        if result["abs"]:
            self._log_found("dipole transitions", len(result["abs"]))
        
        return {
            "abs": pd.DataFrame(result["abs"]) if result["abs"] else pd.DataFrame(),
            "soc": pd.DataFrame(result["soc"]) if result["soc"] else pd.DataFrame()
        }
    
    def _parse_electric_dipole(self) -> Dict[str, pd.DataFrame]:
        """Parse electric dipole absorption spectrum."""
        return self._parse_dipole_spectrum(rx.ELECTRIC_DIPOLE)
    
    def _parse_velocity_dipole(self) -> Dict[str, pd.DataFrame]:
        """Parse velocity dipole absorption spectrum."""
        return self._parse_dipole_spectrum(rx.VELOCITY_DIPOLE)
