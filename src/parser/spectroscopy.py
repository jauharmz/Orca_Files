"""Spectroscopy parser for ORCA output."""

from typing import Optional, Dict
import pandas as pd

from ..core.base_parser import BaseParser
from ..core.data_models import SpectraData
from . import regex_patterns as rx


class SpectroscopyParser(BaseParser):
    """Parse spectroscopy data (IR, Raman, NMR, vibrations)."""
    
    def parse(self) -> SpectraData:
        """Parse spectroscopy data."""
        self.logger.debug("Parsing spectroscopy...")
        
        data = SpectraData()
        
        # Parse vibrations
        data.vibrations = self._parse_vibrations()
        
        # Parse IR spectrum
        data.ir = self._parse_ir()
        
        # Parse Raman spectrum
        data.raman = self._parse_raman()
        
        # Parse NMR
        data.nmr_shielding, data.nmr_coupling = self._parse_nmr()
        
        return data
    
    def _parse_vibrations(self) -> Optional[pd.DataFrame]:
        """Parse vibrational frequencies."""
        blocks = rx.VIBRATIONS.findall(self.text)
        if not blocks:
            return None
        
        last_block = blocks[-1].strip().splitlines()
        vib_data = []
        
        for line in last_block:
            parts = line.split()
            if not parts or not parts[0].rstrip(":").isdigit():
                continue
            try:
                freq = float(parts[1])
                img = 1 if freq < 0 else 0
                vib_data.append({"freq_cm-1": freq, "imaginary": img})
            except (ValueError, IndexError):
                continue
        
        if vib_data:
            self._log_found("vibrations", len(vib_data))
            return pd.DataFrame(vib_data)
        return None
    
    def _parse_ir(self) -> Optional[pd.DataFrame]:
        """Parse IR spectrum."""
        match = rx.IR_SPECTRUM.search(self.text)
        if not match:
            return None
        
        block = match.group(1).strip().splitlines()
        ir_data = []
        
        for line in block:
            parts = line.split()
            if not parts or not parts[0].rstrip(":").isdigit():
                continue
            try:
                freq = float(parts[1])
                eps = float(parts[2])
                intensity = float(parts[3])
                ir_data.append({
                    "freq_cm-1": freq,
                    "eps": eps,
                    "intensity_km/mol": intensity
                })
            except (ValueError, IndexError):
                continue
        
        if ir_data:
            self._log_found("IR peaks", len(ir_data))
            return pd.DataFrame(ir_data)
        return None
    
    def _parse_raman(self) -> Optional[pd.DataFrame]:
        """Parse Raman spectrum."""
        match = rx.RAMAN_SPECTRUM.search(self.text)
        if not match:
            return None
        
        block = match.group(1).strip().splitlines()
        raman_data = []
        
        for line in block:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                freq = float(parts[1])
                activity = float(parts[2])
                depol = float(parts[3])
                raman_data.append({
                    "freq_cm-1": freq,
                    "activity": activity,
                    "depolarization": depol
                })
            except (ValueError, IndexError):
                continue
        
        if raman_data:
            self._log_found("Raman peaks", len(raman_data))
            return pd.DataFrame(raman_data)
        return None
    
    def _parse_nmr(self):
        """Parse NMR shielding and coupling."""
        shielding = None
        coupling = None
        
        # Parse shielding
        match = rx.NMR_SHIELDING.search(self.text)
        if match:
            lines = match.group(1).strip().split('\n')
            data = []
            for line in lines:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        data.append({
                            "Nucleus": parts[0],
                            "Element": parts[1],
                            "Isotropic": float(parts[2]),
                            "Anisotropy": float(parts[3])
                        })
                    except ValueError:
                        continue
            if data:
                shielding = pd.DataFrame(data)
                self._log_found("NMR shielding", len(data))
        
        # Parse coupling
        match = rx.NMR_COUPLING.search(self.text)
        if match:
            lines = match.group(1).strip().split('\n')
            data = []
            if lines:
                header = None
                for i, line in enumerate(lines):
                    parts = line.split()
                    if parts and ('H' in parts or 'C' in parts or 'N' in parts):
                        header = parts
                        break
                
                if header:
                    for line in lines[i+1:]:
                        parts = line.split()
                        if len(parts) >= len(header) + 2:
                            nucleus1 = f"{parts[0]} {parts[1]}"
                            for j, nucleus2 in enumerate(header):
                                try:
                                    j_hz = float(parts[j + 2])
                                    if abs(j_hz) > 1e-6:
                                        data.append({
                                            "Nucleus1": nucleus1,
                                            "Nucleus2": nucleus2,
                                            "J_Hz": j_hz
                                        })
                                except (ValueError, IndexError):
                                    continue
            if data:
                coupling = pd.DataFrame(data)
                self._log_found("NMR couplings", len(data))
        
        return shielding, coupling
