"""Orbital parser for ORCA output."""

from typing import Optional, List, Dict
import pandas as pd

from ..core.base_parser import BaseParser
from ..core.data_models import OrbitalData
from . import regex_patterns as rx


class OrbitalParser(BaseParser):
    """Parse orbital energy data from ORCA output."""
    
    def parse(self) -> OrbitalData:
        """Parse orbital data."""
        self.logger.debug("Parsing orbitals...")
        
        data = OrbitalData()
        data.orbitals = self._parse_orbitals()
        
        if data.orbitals is not None and not data.orbitals.empty:
            # Calculate HOMO/LUMO
            occupied = data.orbitals[data.orbitals["OCC"] > 0]
            virtual = data.orbitals[data.orbitals["OCC"] == 0]
            
            if not occupied.empty:
                data.homo_energy = occupied["eV"].max()
            if not virtual.empty:
                data.lumo_energy = virtual["eV"].min()
            if data.homo_energy and data.lumo_energy:
                data.homo_lumo_gap = data.lumo_energy - data.homo_energy
            
            self._log_found(f"orbitals", len(data.orbitals))
            if data.homo_lumo_gap:
                self.logger.debug(f"HOMO-LUMO gap: {data.homo_lumo_gap:.2f} eV")
        
        return data
    
    def _parse_orbitals(self) -> Optional[pd.DataFrame]:
        """Parse orbital energies with spin handling."""
        # Check for open-shell (spin polarized)
        spin_up = rx.SPIN_UP_BLOCK.findall(self.text)
        spin_down = rx.SPIN_DOWN_BLOCK.findall(self.text)
        
        if spin_up and spin_down:
            spins = [("up", spin_up[-1]), ("down", spin_down[-1])]
        else:
            # Closed-shell
            blocks = rx.ORBITAL_BLOCK.findall(self.text)
            if not blocks:
                self._log_not_found("orbital energies")
                return None
            spins = [("na", blocks[-1])]
        
        orbitals = []
        for spin, block in spins:
            lines = block.strip().splitlines()
            for line in lines:
                parts = line.split()
                if not parts or not parts[0].isdigit() or len(parts) < 4:
                    continue
                
                occ = self._parse_float(parts[1])
                energy_h = self._parse_float(parts[2])
                energy_ev = self._parse_float(parts[3])
                
                if None in (occ, energy_h, energy_ev):
                    continue
                
                # Filter fake data
                if occ not in [0.0, 1.0, 2.0]:
                    continue
                
                orbitals.append({
                    "OCC": occ,
                    "Eh": energy_h,
                    "eV": energy_ev,
                    "spin": spin
                })
        
        if not orbitals:
            return None
        
        # Assign HOMO/LUMO levels
        occupied = [o for o in orbitals if o["OCC"] > 0]
        virtual = [o for o in orbitals if o["OCC"] == 0]
        
        occupied.sort(key=lambda x: x["Eh"], reverse=True)
        virtual.sort(key=lambda x: x["Eh"])
        
        for i, orb in enumerate(occupied):
            orb["lvl"] = -i  # HOMO=0, HOMO-1=-1, etc.
            orb["label"] = "HOMO" if i == 0 else f"HOMO-{i}"
        
        for i, orb in enumerate(virtual):
            orb["lvl"] = i + 1  # LUMO=1, LUMO+1=2, etc.
            orb["label"] = "LUMO" if i == 0 else f"LUMO+{i}"
        
        all_orbitals = occupied + virtual
        all_orbitals.sort(key=lambda x: x["Eh"])
        
        return pd.DataFrame(all_orbitals)
