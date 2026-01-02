"""
Geometry Parser - Refactored from orca_praser.py

Parses: cartesian coordinates, geometry info, internal coordinates, SMILES
"""

import re
import os
from typing import Optional, List, Tuple, Dict, Any, Union
from tempfile import NamedTemporaryFile
import pandas as pd

from .regex_patterns import (
    RE_CARTESIAN_BLOCKS, RE_INTERNAL_BLOCKS, RE_INPUT_FILE, 
    RE_OPT_PARAMS, FLOAT_RE
)
from ..core.base_parser import BaseParser
from ..core.data_models import GeometryData, InternalCoordsData
from ..logger import get_logger

# Try to import pybel for SMILES
try:
    from openbabel import pybel
    _HAS_PYBEL = True
except Exception:
    _HAS_PYBEL = False


def _to_float(x: Any) -> Optional[float]:
    """Safe float conversion."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def _safe_unlink(path: str):
    """Safe file removal."""
    try:
        os.unlink(path)
    except Exception:
        pass


class GeometryParser(BaseParser):
    """Parser for geometry data."""
    
    def __init__(self, text: str):
        super().__init__(text)
        self.logger = get_logger("GeometryParser")
    
    def parse(self) -> GeometryData:
        """Parse geometry info and coordinates."""
        self.logger.debug("Starting geometry parsing...")
        
        if not self.text:
            self.logger.warning("No text to parse")
            return GeometryData()
        
        filename = self._parse_filename()
        charge = self._parse_charge()
        multiplicity = self._parse_multiplicity()
        cart_coords = self.parse_last_cartesian()
        smiles = self._generate_smiles()
        
        self.logger.info(f"Geometry parsed: filename={filename}, charge={charge}, mult={multiplicity}")
        if cart_coords is not None:
            self.logger.info(f"  Atoms: {len(cart_coords)} coordinates found")
        if smiles:
            self.logger.info(f"  SMILES: {smiles}")
        
        return GeometryData(
            filename=filename,
            charge=charge,
            multiplicity=multiplicity,
            cart_coords=cart_coords,
            smiles=smiles
        )
    
    def parse_geometry_info(self) -> Optional[Dict[str, Any]]:
        """Parse geometry info from input section."""
        self.logger.debug("Parsing geometry info...")
        
        if not self.text:
            return None
        
        info = {"filename": None, "charge": None, "multiplicity": None}
        
        # Try to find filename from coordinates declaration
        m = re.search(r"The coordinates will be read from file:\s*(\S+)", self.text)
        if m:
            info["filename"] = os.path.splitext(m.group(1))[0]
            self.logger.debug(f"  Found filename from coords declaration: {info['filename']}")
        
        # Parse from INPUT FILE section
        in_m = RE_INPUT_FILE.search(self.text)
        if in_m:
            inp = in_m.group(2)
            m2 = re.search(r"\*\s*xyzfile\s+(\d+)\s+(\d+)\s+(\S+)", inp, flags=re.I)
            if m2:
                try:
                    info["charge"] = int(m2.group(1))
                    info["multiplicity"] = int(m2.group(2))
                    info["filename"] = os.path.splitext(m2.group(3))[0]
                    self.logger.debug(f"  Found from xyzfile: charge={info['charge']}, mult={info['multiplicity']}, file={info['filename']}")
                except Exception:
                    pass
        
        return info if info.get("filename") else None
    
    def _parse_filename(self) -> Optional[str]:
        """Extract filename."""
        info = self.parse_geometry_info()
        return info.get("filename") if info else None
    
    def _parse_charge(self) -> Optional[int]:
        """Extract charge."""
        info = self.parse_geometry_info()
        return info.get("charge") if info else None
    
    def _parse_multiplicity(self) -> Optional[int]:
        """Extract multiplicity."""
        info = self.parse_geometry_info()
        return info.get("multiplicity") if info else None
    
    def parse_last_cartesian(self) -> Optional[pd.DataFrame]:
        """Parse the last Cartesian coordinates block."""
        self.logger.debug("Parsing Cartesian coordinates...")
        
        if not self.text:
            return None
        
        blocks = RE_CARTESIAN_BLOCKS.findall(self.text)
        if not blocks:
            self.logger.debug("  No Cartesian blocks found")
            return None
        
        self.logger.debug(f"  Found {len(blocks)} Cartesian blocks, using last one")
        
        last = blocks[-1].strip().splitlines()
        coords: List[Tuple[str, float, float, float]] = []
        
        for L in last:
            parts = L.split()
            if len(parts) < 4:
                continue
            atom = parts[0]
            x = _to_float(parts[1])
            y = _to_float(parts[2])
            z = _to_float(parts[3])
            if None in (x, y, z):
                continue
            coords.append((atom, x, y, z))
        
        if not coords:
            self.logger.debug("  No valid coordinates parsed")
            return None
        
        self.logger.info(f"  Parsed {len(coords)} atoms from Cartesian coordinates")
        
        # Log first few atoms
        for atom, x, y, z in coords[:3]:
            self.logger.debug(f"    {atom}: ({x:.4f}, {y:.4f}, {z:.4f})")
        if len(coords) > 3:
            self.logger.debug(f"    ... and {len(coords)-3} more atoms")
        
        return pd.DataFrame(coords, columns=["atom", "x", "y", "z"])
    
    def coords_to_smiles(self, coords: List[Tuple[str, float, float, float]]) -> Optional[str]:
        """Convert coordinates to SMILES string."""
        if not coords:
            return None
        if not _HAS_PYBEL:
            self.logger.debug("OpenBabel not available, skipping SMILES generation")
            return None
        
        self.logger.debug(f"Generating SMILES from {len(coords)} atoms...")
        
        with NamedTemporaryFile("w", suffix=".xyz", delete=False) as tmp:
            tmp.write(f"{len(coords)}\n")
            tmp.write("orca coords\n")
            for atom, x, y, z in coords:
                tmp.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")
            tmp_path = tmp.name
        
        try:
            mol = next(pybel.readfile("xyz", tmp_path))
            smi = mol.write("smi").strip()
            result = smi.split()[0] if smi else None
            if result:
                self.logger.debug(f"  Generated SMILES: {result}")
            return result
        except Exception as e:
            self.logger.debug(f"  SMILES generation failed: {e}")
            return None
        finally:
            _safe_unlink(tmp_path)
    
    def _generate_smiles(self) -> Optional[str]:
        """Generate SMILES from parsed coordinates."""
        cart_df = self.parse_last_cartesian()
        if cart_df is None or cart_df.empty:
            return None
        
        coords = list(cart_df.itertuples(index=False, name=None))
        return self.coords_to_smiles(coords)
    
    def parse_internal(self, as_df: bool = True) -> InternalCoordsData:
        """
        Parse internal coordinates with priority:
        1. Redundant Internal Coordinates (from optimization) - most detailed
        2. Optimized Parameters block (fallback)
        3. Basic Internal Coordinates block (last resort)
        """
        self.logger.debug("Parsing internal coordinates...")
        
        if not self.text:
            self.logger.debug("  No text to parse")
            return InternalCoordsData(
                bonds=pd.DataFrame() if as_df else [],
                angles=pd.DataFrame() if as_df else [],
                dihedrals=pd.DataFrame() if as_df else []
            )
        
        out = {"bonds": [], "angles": [], "dihedrals": []}
        
        # PRIORITY 1: Try Redundant Internal Coordinates
        pattern = re.compile(
            r"Redundant Internal Coordinates.*?Definition.*?(?:OldVal|Value).*?dE/dq.*?Step.*?(?:FinalVal|New-Value).*?"
            r"\n\s*-+\s*\n(.*?)(?:\n\s*-{40,}|\Z)",
            flags=re.S | re.I
        )
        
        all_matches = list(pattern.finditer(self.text))
        
        if all_matches:
            self.logger.debug(f"  Found {len(all_matches)} Redundant Internal Coords blocks, using last")
            redundant_match = all_matches[-1]
            block = redundant_match.group(1).strip()
            lines = block.split('\n')
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('-') or len(line) < 10:
                    continue
                
                match_line = re.match(
                    r'^\s*(\d+)\.\s+([BAD])\(([^\)]+)\)\s+(' + FLOAT_RE + r')\s+(' + 
                    FLOAT_RE + r')\s+(' + FLOAT_RE + r')\s+(' + FLOAT_RE + r')\s*$',
                    line
                )
                
                if not match_line:
                    continue
                
                try:
                    index = int(match_line.group(1))
                    coord_type = match_line.group(2)
                    definition = match_line.group(3).strip()
                    old_val = _to_float(match_line.group(4))
                    dedq = _to_float(match_line.group(5))
                    step = _to_float(match_line.group(6))
                    final_val = _to_float(match_line.group(7))
                    
                    if None in (old_val, dedq, step, final_val):
                        continue
                    
                    atoms_clean = re.sub(r'\s+', '', definition).replace(',', '-')
                    
                    entry = {
                        "Index": index,
                        "Type": coord_type,
                        "Definition": f"{coord_type}({definition})",
                        "Atoms": atoms_clean,
                        "OldVal": old_val,
                        "dE_dq": dedq,
                        "Step": step,
                        "FinalVal": final_val,
                        "Unit": "Angstrom" if coord_type == "B" else "Degrees"
                    }
                    
                    if coord_type == "B":
                        out["bonds"].append(entry)
                    elif coord_type == "A":
                        out["angles"].append(entry)
                    elif coord_type == "D":
                        out["dihedrals"].append(entry)
                        
                except (ValueError, IndexError, AttributeError):
                    continue
            
            if out["bonds"] or out["angles"] or out["dihedrals"]:
                self.logger.info(f"  Internal coords (Redundant): {len(out['bonds'])} bonds, {len(out['angles'])} angles, {len(out['dihedrals'])} dihedrals")
                return InternalCoordsData(
                    bonds=pd.DataFrame(out["bonds"]) if as_df else out["bonds"],
                    angles=pd.DataFrame(out["angles"]) if as_df else out["angles"],
                    dihedrals=pd.DataFrame(out["dihedrals"]) if as_df else out["dihedrals"]
                )
        
        # PRIORITY 2: Optimized Parameters block
        m = RE_OPT_PARAMS.search(self.text)
        if m:
            self.logger.debug("  Trying Optimized Parameters block")
            for L in m.group(1).strip().splitlines():
                dm = re.search(r"([BAD])\((.*?)\).*?([-+]?\d*\.\d+)$", L.strip())
                if not dm:
                    continue
                ctype = dm.group(1).upper()
                defn = dm.group(2).strip()
                val = _to_float(dm.group(3))
                
                entry = {
                    "Type": ctype,
                    "Definition": f"{ctype}({defn})",
                    "Atoms": re.sub(r'\s+', '', defn).replace(',', '-'),
                    "FinalVal": val,
                    "Unit": "Angstrom" if ctype == "B" else "Degrees"
                }
                
                if ctype == "B":
                    out["bonds"].append(entry)
                elif ctype == "A":
                    out["angles"].append(entry)
                elif ctype == "D":
                    out["dihedrals"].append(entry)
            
            if out["bonds"] or out["angles"] or out["dihedrals"]:
                self.logger.info(f"  Internal coords (OptParams): {len(out['bonds'])} bonds, {len(out['angles'])} angles, {len(out['dihedrals'])} dihedrals")
                return InternalCoordsData(
                    bonds=pd.DataFrame(out["bonds"]) if as_df else out["bonds"],
                    angles=pd.DataFrame(out["angles"]) if as_df else out["angles"],
                    dihedrals=pd.DataFrame(out["dihedrals"]) if as_df else out["dihedrals"]
                )
        
        # PRIORITY 3: Basic Internal Coordinates block
        b = RE_INTERNAL_BLOCKS.search(self.text)
        if not b:
            self.logger.debug("  No internal coordinates found")
            return InternalCoordsData(
                bonds=pd.DataFrame() if as_df else [],
                angles=pd.DataFrame() if as_df else [],
                dihedrals=pd.DataFrame() if as_df else []
            )
        
        self.logger.debug("  Using basic Internal Coordinates block")
        lines = b.group(1).strip().splitlines()
        atom_list = []
        for idx, L in enumerate(lines):
            parts = L.split()
            if len(parts) < 5:
                continue
            atom = parts[0]
            i = int(parts[1]) if parts[1].isdigit() else 0
            j = int(parts[2]) if parts[2].isdigit() else 0
            k = int(parts[3]) if parts[3].isdigit() else 0
            bond = _to_float(parts[4])
            angle = _to_float(parts[5]) if len(parts) > 5 else None
            dihedral = _to_float(parts[6]) if len(parts) > 6 else None
            cur = f"{atom}{idx}"
            
            if i > 0 and (i-1) < len(atom_list) and bond is not None:
                out["bonds"].append({
                    "Type": "B",
                    "Definition": f"B({cur},{atom_list[i-1]})",
                    "Atoms": f"{cur}-{atom_list[i-1]}",
                    "FinalVal": bond,
                    "Unit": "Angstrom"
                })
            if i > 0 and j > 0 and (i-1) < len(atom_list) and (j-1) < len(atom_list) and angle is not None:
                out["angles"].append({
                    "Type": "A",
                    "Definition": f"A({cur},{atom_list[i-1]},{atom_list[j-1]})",
                    "Atoms": f"{cur}-{atom_list[i-1]}-{atom_list[j-1]}",
                    "FinalVal": angle,
                    "Unit": "Degrees"
                })
            if i > 0 and j > 0 and k > 0 and (i-1) < len(atom_list) and (j-1) < len(atom_list) and (k-1) < len(atom_list) and dihedral is not None:
                out["dihedrals"].append({
                    "Type": "D",
                    "Definition": f"D({cur},{atom_list[i-1]},{atom_list[j-1]},{atom_list[k-1]})",
                    "Atoms": f"{cur}-{atom_list[i-1]}-{atom_list[j-1]}-{atom_list[k-1]}",
                    "FinalVal": dihedral,
                    "Unit": "Degrees"
                })
            atom_list.append(cur)
        
        self.logger.info(f"  Internal coords (Basic): {len(out['bonds'])} bonds, {len(out['angles'])} angles, {len(out['dihedrals'])} dihedrals")
        
        return InternalCoordsData(
            bonds=pd.DataFrame(out["bonds"]) if as_df else out["bonds"],
            angles=pd.DataFrame(out["angles"]) if as_df else out["angles"],
            dihedrals=pd.DataFrame(out["dihedrals"]) if as_df else out["dihedrals"]
        )
