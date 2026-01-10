"""
Energy Parser - Refactored from orca_praser.py

Parses: Gibbs free energy, single point energy, ESD flags, calc type
"""

import re
from typing import Optional, Dict, Any

from .regex_patterns import RE_FINAL_GIBBS, RE_FINAL_SPE, RE_INPUT_FILE
from ..core.base_parser import BaseParser
from ..core.data_models import EnergyData
from ..logger import get_logger


class EnergyParser(BaseParser):
    """Parser for energy data."""
    
    def __init__(self, text: str):
        super().__init__(text)
        self.logger = get_logger("EnergyParser")
    
    def parse(self) -> EnergyData:
        """Parse energy values."""
        self.logger.debug("Starting energy parsing...")
        
        energies = self.parse_energies()
        
        gibbs = energies.get("gibbs_Eh")
        sp = energies.get("single_point_Eh")
        
        if gibbs:
            self.logger.info(f"  Gibbs energy: {gibbs:.6f} Eh")
        if sp:
            self.logger.info(f"  Single-point energy: {sp:.6f} Eh")
        if not gibbs and not sp:
            self.logger.debug("  No energy values found")
        
        return EnergyData(
            gibbs_Eh=gibbs,
            single_point_Eh=sp
        )
    
    def parse_energies(self) -> Dict[str, Optional[float]]:
        """Parse Gibbs and single point energies."""
        if not self.text:
            return {"gibbs_Eh": None, "single_point_Eh": None}
        
        gib = RE_FINAL_GIBBS.search(self.text)
        spe = RE_FINAL_SPE.search(self.text)
        
        return {
            "gibbs_Eh": float(gib.group(1)) if gib else None,
            "single_point_Eh": float(spe.group(1)) if spe else None
        }
    
    def parse_esd_flag(self) -> Optional[str]:
        """Parse ESD flag (VG/AH/AHAS/FLUOR/PHOSP)."""
        self.logger.debug("Parsing ESD flag...")
        
        if not self.text:
            return None
        
        m = RE_INPUT_FILE.search(self.text)
        if m:
            s = m.group(2)
            m1 = re.search(r"ESD\s*\(\s*([A-Za-z]+)\s*\)", s, flags=re.I)
            if m1:
                flag = m1.group(1).upper()
                self.logger.info(f"  ESD flag: {flag}")
                return flag
            m2 = re.search(r"HESSFLAG\s+(\w+)", s, flags=re.I)
            if m2:
                flag = m2.group(1).upper()
                self.logger.info(f"  HESSFLAG: {flag}")
                return flag
            if re.search(r"HESSFLAG\s+VG", s, flags=re.I):
                self.logger.info("  ESD flag: VG")
                return "VG"
            if re.search(r"HESSFLAG\s+AHAS", s, flags=re.I):
                self.logger.info("  ESD flag: AHAS")
                return "AHAS"
            if re.search(r"ESD\s*\(\s*ABS", s, flags=re.I):
                self.logger.info("  ESD flag: ABS")
                return "ABS"
        
        # fallback
        if re.search(r"HESSFLAG\s+VG", self.text, flags=re.I):
            self.logger.info("  ESD flag (fallback): VG")
            return "VG"
        if re.search(r"HESSFLAG\s+AHAS", self.text, flags=re.I):
            self.logger.info("  ESD flag (fallback): AHAS")
            return "AHAS"
        
        self.logger.debug("  No ESD flag found")
        return None
    
    def detect_calc_type(self) -> Dict[str, Any]:
        """Detect calculation type."""
        self.logger.debug("Detecting calculation type...")
        
        out = {
            "is_optimization": False,
            "has_tddft": False,
            "has_esd": False,
            "esd_type": None,
            "optimized_state": "S0",
            "calc_class": "single_point"
        }
        
        if not self.text:
            return out
        
        # Check for optimization
        if re.search(r"GEOMETRY OPTIMIZATION", self.text, flags=re.I) or re.search(r"OPTIMIZATION CONVERGED", self.text, flags=re.I):
            out["is_optimization"] = True
            out["calc_class"] = "optimization"
            self.logger.info("  Calc type: Geometry Optimization")
        
        # Check for TD-DFT
        if re.search(r"TD-DFT", self.text, flags=re.I):
            out["has_tddft"] = True
            if out["calc_class"] == "single_point":
                out["calc_class"] = "tddft"
            self.logger.info("  TD-DFT: Yes")
        
        # Check ESD
        esd = self.parse_esd_flag()
        if esd:
            out["has_esd"] = True
            out["esd_type"] = esd
            if out["calc_class"] == "single_point":
                out["calc_class"] = "spectrum"
        
        # find IROOT in %TDDFT block
        tblock = re.search(r"%TDDFT\s*(.*?)END", self.text, flags=re.S|re.I)
        if tblock:
            mroot = re.search(r"IROOT\s+(\d+)", tblock.group(1), flags=re.I)
            if mroot:
                try:
                    root = int(mroot.group(1))
                    out["optimized_state"] = f"S{root}"
                    self.logger.info(f"  Optimized state: S{root}")
                except Exception:
                    pass
        
        # multiplicity inference
        in_m = RE_INPUT_FILE.search(self.text)
        if in_m:
            mxyz = re.search(r"\*\s*xyzfile\s+(\d+)\s+(\d+)", in_m.group(2), flags=re.I)
            if mxyz:
                try:
                    mult = int(mxyz.group(2))
                    if mult == 3 and out["optimized_state"] == "S0":
                        out["optimized_state"] = "T1"
                        self.logger.info("  Optimized state: T1 (triplet)")
                except Exception:
                    pass
        
        self.logger.info(f"  Calc class: {out['calc_class']}")
        
        return out
