"""
Spectrum File Parser - Parse external spectrum files (VG/AH/AHAS/FLUOR/PHOSP)

Refactored from orca_praser.py
"""

import os
import re
from typing import Optional, List, Dict, Any, Union
import pandas as pd

from ..logger import get_logger


def _to_float(x: Any) -> Optional[float]:
    """Safe float conversion."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


class SpectrumFileParser:
    """Parser for external spectrum files."""
    
    def __init__(self):
        self.logger = get_logger("SpectrumFileParser")
    
    def parse_spectrum_file(
        self, 
        filename: str, 
        spectrum_type: Optional[str] = None, 
        as_df: bool = True
    ) -> Union[pd.DataFrame, List[Dict[str, Any]]]:
        """
        Parse external spectrum file.
        
        Args:
            filename: Path to spectrum file
            spectrum_type: Type of spectrum (VG, AH, AHAS, FLUOR, PHOSP)
            as_df: Return as DataFrame if True
            
        Returns:
            DataFrame or list of dicts with spectrum data
        """
        self.logger.debug(f"Parsing spectrum file: {filename}")
        
        try:
            with open(filename, "r") as fh:
                lines = fh.readlines()
        except Exception as e:
            self.logger.warning(f"Could not read {filename}: {e}")
            return pd.DataFrame() if as_df else []
        
        data = []
        for L in lines[1:]:  # Skip header
            parts = re.split(r"\s+", L.strip())
            if not parts:
                continue
            
            energy = _to_float(parts[0])
            if energy is None:
                continue
            
            if spectrum_type and spectrum_type.upper() in ("VG", "AH", "AHAS"):
                # VG/AH/AHAS format: energy, total, FC, HT
                if len(parts) < 4:
                    continue
                tot = _to_float(parts[1])
                fc = _to_float(parts[2])
                ht = _to_float(parts[3])
                if None in (tot, fc, ht):
                    continue
                data.append({
                    "energy_cm1": energy,
                    "total_spectrum": tot,
                    "intensity_fc": fc,
                    "intensity_ht": ht
                })
            elif spectrum_type and spectrum_type.upper() in ("FLUOR", "PHOSP"):
                # Fluorescence/Phosphorescence format: energy, intensity
                if len(parts) < 2:
                    continue
                tot = _to_float(parts[1])
                data.append({
                    "energy_cm1": energy,
                    "total_spectrum": tot
                })
            else:
                # Generic format: energy + columns
                nums = [_to_float(x) for x in parts[1:]]
                entry = {"energy_cm1": energy}
                for i, v in enumerate(nums):
                    if v is not None:
                        entry[f"col{i+1}"] = v
                data.append(entry)
        
        if data:
            self.logger.info(f"  Parsed {len(data)} spectrum points")
        
        return pd.DataFrame(data) if as_df else data
    
    def find_spectrum_files(
        self, 
        base_filename: str, 
        esd_flag: Optional[str] = None,
        geometry_filename: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Find spectrum files associated with an ORCA output.
        
        Args:
            base_filename: Path to ORCA .out file
            esd_flag: ESD type (VG, AH, AHAS, etc.)
            geometry_filename: Optional molecule filename from geometry
            
        Returns:
            Dict mapping spectrum type to file path
        """
        base_name = os.path.splitext(base_filename)[0]
        
        candidates = [f"{base_name}.spectrum"]
        if esd_flag:
            candidates.append(f"{base_name}_{esd_flag.lower()}.spectrum")
        if geometry_filename:
            geom_base = os.path.splitext(geometry_filename)[0]
            candidates.append(f"{geom_base}.spectrum")
            if esd_flag:
                candidates.append(f"{geom_base}_{esd_flag.lower()}.spectrum")
        
        found = {}
        for cand in candidates:
            if os.path.exists(cand):
                spec_type = esd_flag or "UNKNOWN"
                found[spec_type] = cand
                self.logger.debug(f"  Found spectrum file: {cand}")
                break
        
        return found
