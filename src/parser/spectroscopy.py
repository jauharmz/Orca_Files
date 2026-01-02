"""
Spectroscopy Parser - Refactored from orca_praser.py

Parses: vibrations, IR, Raman, Mulliken, NMR
"""

import re
from typing import Optional, List, Dict, Any, Union
import pandas as pd

from .regex_patterns import RE_VIBS, RE_IR, RE_RAMAN, RE_MULLIKEN, RE_NMR_SHIELD, RE_NMR_COUPLING
from ..core.base_parser import BaseParser
from ..core.data_models import SpectraData, MullikenData
from ..logger import get_logger


def _to_float(x: Any) -> Optional[float]:
    """Safe float conversion."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


class SpectroscopyParser(BaseParser):
    """Parser for spectroscopic data."""
    
    def __init__(self, text: str):
        super().__init__(text)
        self.logger = get_logger("SpectroscopyParser")
    
    def parse(self, as_df: bool = True) -> SpectraData:
        """Parse all spectroscopic data."""
        self.logger.debug("Starting spectroscopy parsing...")
        
        vib = self.parse_vibrations()
        ir = self.parse_ir_spectrum()
        raman = self.parse_raman_spectrum(as_df=as_df)
        nmr = self.parse_nmr(as_df=as_df)
        
        return SpectraData(
            vibrations=pd.DataFrame(vib) if vib and as_df else vib,
            ir=pd.DataFrame(ir) if ir and as_df else ir,
            raman=raman,
            nmr_shielding=nmr.get("shielding"),
            nmr_coupling=nmr.get("coupling")
        )
    
    def parse_vibrations(self) -> Optional[List[Dict[str, Any]]]:
        """Parse vibrational frequencies."""
        self.logger.debug("Parsing vibrations...")
        
        if not self.text:
            return None
        
        blocks = RE_VIBS.findall(self.text)
        if not blocks:
            self.logger.debug("  No vibration blocks found")
            return None
        
        last = blocks[-1].strip().splitlines()
        out: List[Dict[str, Any]] = []
        
        for L in last:
            parts = L.split()
            if not parts:
                continue
            if not parts[0].rstrip(":").isdigit():
                continue
            freq = _to_float(parts[1])
            if freq is None:
                continue
            out.append({"freq_cm-1": freq, "img": 1 if freq < 0 else 0})
        
        if out:
            img_count = sum(1 for v in out if v["img"])
            self.logger.info(f"  Vibrations: {len(out)} modes ({img_count} imaginary)")
        else:
            self.logger.debug("  No vibrations parsed")
        
        return out if out else None
    
    def parse_ir_spectrum(self) -> Optional[List[Dict[str, Any]]]:
        """Parse IR spectrum."""
        self.logger.debug("Parsing IR spectrum...")
        
        if not self.text:
            return None
        
        m = RE_IR.search(self.text)
        if not m:
            self.logger.debug("  No IR spectrum found")
            return None
        
        block = m.group(1).strip().splitlines()
        out = []
        
        for L in block:
            parts = L.split()
            if len(parts) < 4:
                continue
            if not parts[0].rstrip(":").isdigit():
                continue
            freq = _to_float(parts[1])
            eps = _to_float(parts[2])
            intensity = _to_float(parts[3])
            if None in (freq, eps, intensity):
                continue
            out.append({"freq_cm-1": freq, "eps": eps, "intensity_km/mol": intensity})
        
        if out:
            self.logger.info(f"  IR spectrum: {len(out)} peaks")
            max_int = max(p["intensity_km/mol"] for p in out)
            self.logger.debug(f"    Max intensity: {max_int:.2f} km/mol")
        else:
            self.logger.debug("  No IR peaks parsed")
        
        return out if out else None
    
    def parse_raman_spectrum(self, as_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, Any]], None]:
        """Parse Raman spectrum."""
        self.logger.debug("Parsing Raman spectrum...")
        
        if not self.text:
            return pd.DataFrame() if as_df else None
        
        m = RE_RAMAN.search(self.text)
        if not m:
            self.logger.debug("  No Raman spectrum found")
            return pd.DataFrame() if as_df else None
        
        block = m.group(1).strip().splitlines()
        out: List[Dict[str, Any]] = []
        
        for L in block:
            parts = re.split(r"\s+", L.strip())
            if not parts or not re.match(r"^\d+:?$", parts[0]):
                if not parts[0].isdigit():
                    continue
            nums = [p for p in parts if re.match(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$", p)]
            if len(nums) < 3:
                continue
            freq = _to_float(nums[0])
            activity = _to_float(nums[1])
            depol = _to_float(nums[2])
            if None in (freq, activity, depol):
                continue
            out.append({"freq_cm-1": freq, "activity": activity, "depolarization": depol})
        
        if out:
            self.logger.info(f"  Raman spectrum: {len(out)} peaks")
        else:
            self.logger.debug("  No Raman peaks parsed")
        
        return pd.DataFrame(out) if as_df else out
    
    def parse_mulliken(self, as_df: bool = True) -> MullikenData:
        """Parse Mulliken population analysis."""
        self.logger.debug("Parsing Mulliken charges...")
        
        if not self.text:
            return MullikenData(charges=pd.DataFrame() if as_df else [])
        
        m = RE_MULLIKEN.search(self.text)
        if not m:
            self.logger.debug("  No Mulliken data found")
            return MullikenData(charges=pd.DataFrame() if as_df else [])
        
        lines = m.group(1).strip().splitlines()
        out = []
        
        for L in lines:
            parts = re.split(r"\s+", L.strip())
            if len(parts) < 4:
                continue
            pop = _to_float(parts[2])
            charge = _to_float(parts[3])
            if None in (pop, charge):
                continue
            out.append({"Nucleus": parts[0], "Element": parts[1], "Population": pop, "Charge": charge})
        
        if out:
            self.logger.info(f"  Mulliken charges: {len(out)} atoms")
            # Show most charged atoms
            sorted_by_charge = sorted(out, key=lambda x: abs(x["Charge"]), reverse=True)
            for atom in sorted_by_charge[:3]:
                self.logger.debug(f"    {atom['Element']}{atom['Nucleus']}: {atom['Charge']:+.3f}")
        else:
            self.logger.debug("  No Mulliken charges parsed")
        
        return MullikenData(charges=pd.DataFrame(out) if as_df else out)
    
    def parse_nmr(self, as_df: bool = True) -> Dict[str, Any]:
        """Parse NMR data (shielding and coupling)."""
        self.logger.debug("Parsing NMR data...")
        
        out = {
            "shielding": pd.DataFrame() if as_df else [],
            "coupling": pd.DataFrame() if as_df else []
        }
        
        if not self.text:
            return out
        
        # Shielding
        m = RE_NMR_SHIELD.search(self.text)
        if m:
            rows = []
            for L in m.group(1).strip().splitlines():
                parts = re.split(r"\s+", L.strip())
                if len(parts) < 4:
                    continue
                iso = _to_float(parts[2])
                aniso = _to_float(parts[3])
                if None in (iso, aniso):
                    continue
                rows.append({"Nucleus": parts[0], "Element": parts[1], "Isotropic": iso, "Anisotropy": aniso})
            
            if rows:
                self.logger.info(f"  NMR shielding: {len(rows)} nuclei")
            out["shielding"] = pd.DataFrame(rows) if as_df else rows
        
        # Coupling
        m2 = RE_NMR_COUPLING.search(self.text)
        if m2:
            rows = []
            lines = m2.group(1).strip().splitlines()
            header = None
            data_start = 0
            
            for i, L in enumerate(lines):
                parts = re.split(r"\s+", L.strip())
                if parts and any(re.match(r"^[A-Z][a-z]?$", p) for p in parts):
                    header = parts
                    data_start = i + 1
                    break
            
            if header:
                for L in lines[data_start:]:
                    parts = re.split(r"\s+", L.strip())
                    if len(parts) < len(header) + 2:
                        continue
                    nucleus1 = f"{parts[0]} {parts[1]}"
                    for i_h, nuc2 in enumerate(header):
                        j_hz = _to_float(parts[i_h + 2])
                        if j_hz and abs(j_hz) > 1e-6:
                            rows.append({"Nucleus1": nucleus1, "Nucleus2": nuc2, "J_Hz": j_hz})
            
            if rows:
                self.logger.info(f"  NMR coupling: {len(rows)} pairs")
            out["coupling"] = pd.DataFrame(rows) if as_df else rows
        
        if out["shielding"].empty if as_df else not out["shielding"]:
            if out["coupling"].empty if as_df else not out["coupling"]:
                self.logger.debug("  No NMR data found")
        
        return out
