"""
Method Parser - Parse ORCA method descriptor from input header

Extracts: formalism, functional, basis_set, dispersion, relativistic, SOC, solvent

See ARCHITECTURE.md for the conceptual model.
"""

import re
from typing import Optional
from ..core.data_models import MethodData
from ..core.base_parser import BaseParser
from ..logger import get_logger


# Known ORCA keywords organized by category
FORMALISMS = {
    "DFT", "HF", "RHF", "UHF", "ROHF", "RKS", "UKS", "ROKS",
    "MP2", "RI-MP2", "DLPNO-MP2", 
    "CCSD", "CCSD(T)", "DLPNO-CCSD", "DLPNO-CCSD(T)",
    "CASSCF", "NEVPT2", "CASPT2",
    "TDDFT", "TDA", "CIS", "CISD"
}

FUNCTIONALS = {
    # GGA
    "PBE", "BP86", "BLYP", "OLYP", "PWLDA",
    # Meta-GGA
    "TPSS", "M06L", "SCAN",
    # Hybrid
    "B3LYP", "PBE0", "BHANDHLYP", "M06", "M062X", "M06-2X", "TPSSh", 
    # Range-separated hybrid
    "CAM-B3LYP", "WB97", "WB97X", "WB97X-D3", "WB97X-D4", "LC-BLYP", "LCBLYP",
    "ωB97X", "ωB97X-D3", "ωB97X-D4", "wB97X", "wB97X-D3", "wB97X-D4",
    # Double hybrid
    "B2PLYP", "B2PLYP-D3", "PWPB95", "DSD-BLYP"
}

BASIS_SETS = {
    # Karlsruhe
    "DEF2-SVP", "DEF2-TZVP", "DEF2-TZVPP", "DEF2-QZVP", "DEF2-QZVPP",
    "DEF2-SV(P)", "DEF2-TZV(P)", 
    # Auxiliary
    "DEF2/J", "DEF2-TZVP/C", "DEF2/JK",
    # Pople
    "6-31G", "6-31G*", "6-31G**", "6-311G", "6-311G*", "6-311+G**", "6-311++G**",
    # Dunning
    "CC-PVDZ", "CC-PVTZ", "CC-PVQZ", "AUG-CC-PVDZ", "AUG-CC-PVTZ",
    # Other
    "STO-3G", "MINI", "LANL2DZ", "LANL2TZ"
}

DISPERSION = {
    "D3", "D3BJ", "D3(BJ)", "D4", "D3(0)", "D2", "VV10", "NL"
}

RELATIVISTIC = {
    "ZORA", "DKH", "DKH2", "X2C", "NORI"
}

SOLVENTS = {
    "WATER", "ETHANOL", "METHANOL", "ACETONITRILE", "DMSO", "THF",
    "CHLOROFORM", "DCM", "DICHLOROMETHANE", "ACETONE", "TOLUENE", "BENZENE",
    "HEXANE", "DIETHYLETHER", "DMF"
}


class MethodParser(BaseParser):
    """Parser for ORCA method descriptor."""
    
    def __init__(self, text: str):
        super().__init__(text)
        self.logger = get_logger("MethodParser")
    
    def parse(self) -> MethodData:
        """Parse method descriptor from ORCA output."""
        self.logger.debug("Parsing method descriptor...")
        
        method = MethodData()
        
        # Find the input line (! ... )
        input_line = self._find_input_line()
        if input_line:
            method.input_line = input_line
            self._parse_input_line(input_line, method)
        
        # Parse %cpcm block for solvent
        self._parse_cpcm_block(method)
        
        # Parse %rel block for relativistic
        self._parse_rel_block(method)
        
        # Log result
        method_id = method.to_id()
        self.logger.info(f"  Method: {method_id}")
        if method.functional:
            self.logger.debug(f"    Functional: {method.functional}")
        if method.basis_set:
            self.logger.debug(f"    Basis: {method.basis_set}")
        if method.dispersion:
            self.logger.debug(f"    Dispersion: {method.dispersion}")
        if method.solvent:
            self.logger.debug(f"    Solvent: {method.solvent}")
        
        return method
    
    def _find_input_line(self) -> Optional[str]:
        """Find the main input line (! keywords...)."""
        # Look for lines starting with !
        pattern = re.compile(r'^[#\s]*!\s*(.+)$', re.MULTILINE)
        
        for match in pattern.finditer(self.text[:5000]):  # Only check first 5000 chars
            line = match.group(1).strip()
            if line:
                self.logger.debug(f"  Found input line: {line}")
                return line
        
        # Alternative: look in INPUT FILE section
        input_section = re.search(
            r'INPUT FILE\s*[-=]+\s*\n(.*?)\n\s*END OF INPUT',
            self.text[:10000], 
            re.DOTALL | re.IGNORECASE
        )
        if input_section:
            for line in input_section.group(1).splitlines():
                if line.strip().startswith('!'):
                    return line.strip()[1:].strip()
        
        return None
    
    def _parse_input_line(self, line: str, method: MethodData):
        """Parse keywords from input line."""
        tokens = line.upper().replace("!", "").split()
        
        for token in tokens:
            # Check formalism
            if token in FORMALISMS:
                if method.formalism is None:
                    method.formalism = token
                # Also set wavefunction type
                if token in {"RKS", "UKS", "ROKS", "RHF", "UHF", "ROHF"}:
                    method.wavefunction = token
                    if token.startswith("R"):
                        method.formalism = "DFT" if "KS" in token else "HF"
                    elif token.startswith("U"):
                        method.formalism = "DFT" if "KS" in token else "HF"
            
            # Check functional
            elif token in FUNCTIONALS or self._normalize_functional(token) in FUNCTIONALS:
                method.functional = self._normalize_functional(token)
                if method.formalism is None:
                    method.formalism = "DFT"
            
            # Check basis set
            elif token in BASIS_SETS or self._is_basis(token):
                if "/" not in token:  # Not auxiliary
                    method.basis_set = token
                else:
                    method.aux_basis = token
            
            # Check dispersion
            elif token in DISPERSION or token.startswith("D3") or token.startswith("D4"):
                method.dispersion = token
            
            # Check relativistic
            elif token in RELATIVISTIC:
                method.relativistic = token
            
            # Check for CPCM/SMD keywords
            elif token in {"CPCM", "SMD"}:
                method.environment = token
            
            # Check for RI/RIJCOSX
            elif token in {"RI", "RIJCOSX", "RIJK", "RIJONX"}:
                pass  # Could track this if needed
            
            # Check for optimization/task keywords (not method)
            elif token in {"OPT", "FREQ", "NUMFREQ", "ENGRAD", "SP", "TIGHTSCF", "NORMALPRINT"}:
                pass  # Task, not method
    
    def _normalize_functional(self, token: str) -> str:
        """Normalize functional names."""
        # Handle ω vs w
        if token.startswith("W") or token.startswith("Ω"):
            return token.replace("W", "ω").replace("Ω", "ω")
        return token
    
    def _is_basis(self, token: str) -> bool:
        """Check if token looks like a basis set."""
        patterns = [
            r"DEF2-\w+",
            r"\d-\d+\+*G\**",
            r"CC-PV[DTQ5]Z",
            r"AUG-CC-PV[DTQ5]Z",
            r"STO-\dG",
            r"LANL\d\w+"
        ]
        return any(re.match(p, token, re.IGNORECASE) for p in patterns)
    
    def _parse_cpcm_block(self, method: MethodData):
        """Parse %cpcm block for solvent."""
        # Look for %cpcm ... end block
        cpcm_pattern = re.compile(
            r'%cpcm\s*(.*?)\s*end',
            re.DOTALL | re.IGNORECASE
        )
        
        match = cpcm_pattern.search(self.text[:10000])
        if match:
            block = match.group(1)
            method.environment = "CPCM"
            
            # Look for SMD
            if re.search(r'smd\s+true', block, re.IGNORECASE):
                method.environment = "SMD"
            
            # Look for solvent
            solvent_match = re.search(r'solvent\s+"?(\w+)"?', block, re.IGNORECASE)
            if solvent_match:
                method.solvent = solvent_match.group(1).upper()
    
    def _parse_rel_block(self, method: MethodData):
        """Parse %rel block for relativistic treatment."""
        rel_pattern = re.compile(
            r'%rel\s*(.*?)\s*end',
            re.DOTALL | re.IGNORECASE
        )
        
        match = rel_pattern.search(self.text[:10000])
        if match:
            block = match.group(1)
            
            # Look for method
            method_match = re.search(r'method\s+(\w+)', block, re.IGNORECASE)
            if method_match:
                method.relativistic = method_match.group(1).upper()
            
            # Look for SOC
            if re.search(r'dosoc\s+true', block, re.IGNORECASE):
                method.soc = "perturbative"
            elif re.search(r'soctype\s+(\w+)', block, re.IGNORECASE):
                soc_match = re.search(r'soctype\s+(\w+)', block, re.IGNORECASE)
                if soc_match:
                    method.soc = soc_match.group(1).lower()
