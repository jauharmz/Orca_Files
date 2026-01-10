"""
Orbital Parser - Refactored from orca_praser.py

Parses: orbital energies, HOMO/LUMO, spin-up/down
"""

import re
from typing import Optional, List, Dict, Any, Union
import pandas as pd

from .regex_patterns import RE_ORBITAL_BLOCK, RE_SPIN_UP, RE_SPIN_DOWN
from ..core.base_parser import BaseParser
from ..core.data_models import OrbitalData
from ..logger import get_logger


def _to_float(x: Any) -> Optional[float]:
    """Safe float conversion."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


class OrbitalParser(BaseParser):
    """Parser for orbital energies."""
    
    def __init__(self, text: str):
        super().__init__(text)
        self.logger = get_logger("OrbitalParser")
    
    def parse(self, as_df: bool = True) -> OrbitalData:
        """Parse orbital energies."""
        self.logger.debug("Starting orbital parsing...")
        
        orbitals_df = self.parse_orbitals(as_df=as_df)
        
        homo_energy = None
        lumo_energy = None
        
        if orbitals_df is not None and not orbitals_df.empty:
            occupied = orbitals_df[orbitals_df["OCC"] > 0]
            virtual = orbitals_df[orbitals_df["OCC"] == 0]
            
            if not occupied.empty:
                homo_energy = occupied["eV"].max()
            if not virtual.empty:
                lumo_energy = virtual["eV"].min()
            
            self.logger.info(f"  Orbitals: {len(occupied)} occupied, {len(virtual)} virtual")
            if homo_energy:
                self.logger.info(f"  HOMO: {homo_energy:.3f} eV")
            if lumo_energy:
                self.logger.info(f"  LUMO: {lumo_energy:.3f} eV")
            if homo_energy and lumo_energy:
                gap = lumo_energy - homo_energy
                self.logger.info(f"  HOMO-LUMO gap: {gap:.3f} eV")
        else:
            self.logger.debug("  No orbital data found")
        
        return OrbitalData(
            orbitals=orbitals_df,
            homo_energy=homo_energy,
            lumo_energy=lumo_energy,
            homo_lumo_gap=lumo_energy - homo_energy if (homo_energy and lumo_energy) else None
        )
    
    def parse_orbitals(self, as_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, Any]], None]:
        """Return DataFrame of last orbital block. Handles open and closed shell."""
        if not self.text:
            return pd.DataFrame() if as_df else None
        
        # Try spin-up/down blocks first (open-shell)
        up_blocks = RE_SPIN_UP.findall(self.text)
        down_blocks = RE_SPIN_DOWN.findall(self.text)
        
        spins = []
        if up_blocks and down_blocks:
            spins = [("up", up_blocks[-1]), ("down", down_blocks[-1])]
            self.logger.debug("  Found open-shell (spin-up/down) orbitals")
        else:
            # fallback: single ORBITAL ENERGIES block (closed-shell)
            blocks = RE_ORBITAL_BLOCK.findall(self.text)
            if not blocks:
                self.logger.debug("  No orbital blocks found")
                return pd.DataFrame() if as_df else None
            spins = [("na", blocks[-1])]
            self.logger.debug("  Found closed-shell orbitals")
        
        orbitals: List[Dict[str, Any]] = []
        for spin_label, block in spins:
            for L in block.strip().splitlines():
                parts = L.split()
                if not parts or not parts[0].lstrip().isdigit():
                    continue
                if len(parts) < 4:
                    continue
                occ = _to_float(parts[1])
                eh = _to_float(parts[2])
                ev = _to_float(parts[3])
                if None in (occ, eh, ev):
                    continue
                if occ not in (0.0, 1.0, 2.0):
                    continue
                if eh == occ + 1.0 and ev == occ + 2.0:
                    continue
                orbitals.append({"OCC": occ, "Eh": eh, "eV": ev, "spin": spin_label})
        
        if not orbitals:
            self.logger.debug("  No valid orbitals parsed")
            return pd.DataFrame() if as_df else None
        
        # separate occupied & virtual, sort, assign levels
        occupied = [o for o in orbitals if o["OCC"] > 0.0]
        virtual = [o for o in orbitals if o["OCC"] == 0.0]
        
        occupied.sort(key=lambda x: x["Eh"], reverse=True)
        virtual.sort(key=lambda x: x["Eh"])
        
        for i, o in enumerate(occupied):
            o["lvl"] = i  # 0 = HOMO, 1 = HOMO-1
        for i, o in enumerate(virtual):
            o["lvl"] = i  # 0 = LUMO, 1 = LUMO+1
        
        all_orbs = occupied + virtual
        all_orbs.sort(key=lambda x: x["Eh"])
        
        self.logger.debug(f"  Parsed {len(all_orbs)} total orbitals")
        
        if as_df:
            return pd.DataFrame(all_orbs)
        return all_orbs
