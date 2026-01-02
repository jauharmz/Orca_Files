"""
ORCA Output Parser (single-file)
Option A: orbitals stored as dict of DataFrames keyed by state (e.g. "S0_OPT", "S0_SP")

Usage:
    from orca_parser import ORCAParser, ORCABatchParser, parse_orca_file, parse_orca_batch

    p = ORCAParser("p1xs0.out")
    out = p.parse(as_df=True)

    df, states = parse_orca_batch("p1x*.out", "OH*.out")
"""

from __future__ import annotations
import os
import re
import glob
from tempfile import NamedTemporaryFile
from typing import Any, Dict, List, Optional, Tuple, Union, Iterable
import pandas as pd

# Optional OpenBabel for smiles (graceful fallback)
try:
    from openbabel import pybel  # type: ignore
    _HAS_PYBEL = True
except Exception:
    _HAS_PYBEL = False

# -----------------------
# Helpers & regexes
# -----------------------

_FLOAT_RE = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

def _to_float(x: Any) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None

def _safe_unlink(path: str) -> None:
    try:
        os.unlink(path)
    except Exception:
        pass

def _is_empty_obj(x: Any) -> bool:
    if x is None:
        return True
    if isinstance(x, pd.DataFrame):
        return x.empty
    if isinstance(x, dict):
        return all(_is_empty_obj(v) for v in x.values()) if x else True
    if isinstance(x, list):
        return len(x) == 0
    return False

# Common regex blocks
_RE_CARTESIAN_BLOCKS = re.compile(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S)
_RE_INTERNAL_BLOCKS = re.compile(r"INTERNAL COORDINATES \(ANGSTROEM\)\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S)
_RE_OPT_PARAMS = re.compile(r"Optimized Parameters.*?Definition.*?FinalVal\n\s*-+\n(.*?)(?:\n\s*-+|\Z)", flags=re.S|re.I)
_RE_VIBS = re.compile(r"VIBRATIONAL FREQUENCIES\n-+\n.*?\n\n(.*?)(?:\n\n|\Z)", flags=re.S)
_RE_IR = re.compile(r"IR SPECTRUM\n-+\n.*?\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S)
# For RAMAN: read the table after header lines; tolerant to different spacing
_RE_RAMAN = re.compile(r"RAMAN SPECTRUM\s*\n-+\s*\n.*?Mode.*?Depolarization.*?\n-+\s*\n(.*?)(?:\n\n|\Z)", flags=re.S|re.I)
_RE_MULLIKEN = re.compile(r"MULLIKEN POPULATION ANALYSIS\s*\n-+\n.*?\n-+\n(.*?)(?:\n\n|\Z)", flags=re.S|re.I)
_RE_NMR_SHIELD = re.compile(r"CHEMICAL SHIELDING SUMMARY \(ppm\)\s*-+\s*Nucleus\s+Element\s+Isotropic\s+Anisotropy\s*-+\s*([\s\S]*?)(?=\n{2,}|\Z)", flags=re.S|re.I)
_RE_NMR_COUPLING = re.compile(r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS\s*J\s*\(Hz\)\s*-+\s*([\s\S]*?)(?=\n{2,}|\Z)", flags=re.S|re.I)
_RE_INPUT_FILE = re.compile(r"={10,}\s*(INPUT FILE|INPUT)\s*\n={10,}\n(.*?)\*{4}END OF INPUT\*{4}", flags=re.S|re.I)

_RE_FINAL_GIBBS = re.compile(r"Final Gibbs free energy\s*[:\.]*\s*(-?\d+\.\d+)\s+Eh", flags=re.I)
_RE_FINAL_SPE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", flags=re.I)
_RE_ORBITAL_BLOCK = re.compile(r"ORBITAL ENERGIES\n-+\n\s*NO.*?\n(.*?)(?:\n\n|\Z)", flags=re.S)
_RE_SPIN_UP = re.compile(r"SPIN UP ORBITALS\n.*?\n(.*?)(?=\n\s*SPIN DOWN ORBITALS|\Z)", flags=re.S)
_RE_SPIN_DOWN = re.compile(r"SPIN DOWN ORBITALS\n.*?\n(.*?)(?:\n\n|\Z)", flags=re.S)
_RE_TDDFT_STATE = re.compile(
    r"STATE\s+(\d+):\s+E=\s+(" + _FLOAT_RE + r")\s+au\s+(" + _FLOAT_RE + r")\s+eV\s+(" + _FLOAT_RE + r")\s+cm\*{1,2}-1([\s\S]*?)(?=STATE|\Z)",
    flags=re.S | re.I
)
_RE_TD_TRANS = re.compile(r"(\d+)[a-zA-Z]?\s*->\s*(\d+)[a-zA-Z]?\s*:\s*(" + _FLOAT_RE + r")\s*\(c=\s*(" + _FLOAT_RE + r")\)")

# absorption / velocity generic
_RE_ABS_BLOCK = re.compile(r"ABSORPTION SPECTRUM.*?Transition.*?\n-+\n(.*?)(?:\n{2,}|\Z)", flags=re.S|re.I)

# -----------------------
# Parser
# -----------------------

class ORCAParser:
    def __init__(self, filename: str):
        self.filename = filename
        self.text: Optional[str] = None
        self._load()

    def _load(self) -> None:
        try:
            with open(self.filename, "r", encoding="utf-8", errors="ignore") as fh:
                self.text = fh.read()
        except FileNotFoundError:
            self.text = None

    # ---------- coords -> smiles (optional) ----------
    def coords_to_smiles(self, coords: List[Tuple[str, float, float, float]]) -> Optional[str]:
        if not coords or not _HAS_PYBEL:
            return None
        with NamedTemporaryFile("w", suffix=".xyz", delete=False) as tmp:
            tmp.write(f"{len(coords)}\n")
            tmp.write("orca coords\n")
            for atom, x, y, z in coords:
                tmp.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")
            tmp_path = tmp.name
        try:
            mol = next(pybel.readfile("xyz", tmp_path))
            smi = mol.write("smi").strip()
            return smi.split()[0] if smi else None
        except Exception:
            return None
        finally:
            _safe_unlink(tmp_path)

    # ---------- geometry / coords ----------
    def parse_geometry_info(self) -> Optional[Dict[str, Any]]:
        if not self.text:
            return None
        info = {"filename": None, "charge": None, "multiplicity": None}
        m = re.search(r"The coordinates will be read from file:\s*(\S+)", self.text)
        if m:
            info["filename"] = os.path.splitext(m.group(1))[0]
        in_m = _RE_INPUT_FILE.search(self.text)
        if in_m:
            inp = in_m.group(2)
            m2 = re.search(r"\*\s*xyzfile\s+(\d+)\s+(\d+)\s+(\S+)", inp, flags=re.I)
            if m2:
                try:
                    info["charge"] = int(m2.group(1))
                    info["multiplicity"] = int(m2.group(2))
                    info["filename"] = os.path.splitext(m2.group(3))[0]
                except Exception:
                    pass
        return info if info.get("filename") else None

    def parse_last_cartesian(self) -> Optional[List[Tuple[str, float, float, float]]]:
        if not self.text:
            return None
        blocks = _RE_CARTESIAN_BLOCKS.findall(self.text)
        if not blocks:
            return None
        last = blocks[-1].strip().splitlines()
        coords: List[Tuple[str, float, float, float]] = []
        for L in last:
            parts = L.split()
            if len(parts) < 4:
                continue
            atom = parts[0]
            x = _to_float(parts[1]); y = _to_float(parts[2]); z = _to_float(parts[3])
            if None in (x, y, z):
                continue
            coords.append((atom, x, y, z))
        return coords if coords else None

    # ---------- internal coords ----------
    def parse_internal(self, as_df: bool = True) -> Union[Dict[str, Any], Dict[str, pd.DataFrame]]:
        """
        Parse internal coordinates with priority:
        1. Redundant Internal Coordinates (from optimization) - most detailed
        2. Optimized Parameters block (fallback)
        3. Basic Internal Coordinates block (last resort)
        
        Returns dict with 'bonds', 'angles', 'dihedrals' containing detailed structural data.
        """
        if not self.text:
            return {"bonds": pd.DataFrame(), "angles": pd.DataFrame(), "dihedrals": pd.DataFrame()} if as_df else {"bonds": [], "angles": [], "dihedrals": []}
        
        out = {"bonds": [], "angles": [], "dihedrals": []}
        
        # PRIORITY 1: Try Redundant Internal Coordinates (most detailed)
        # Find ALL matches and use the LAST one (the optimized parameters at end)
        pattern = re.compile(
            r"Redundant Internal Coordinates.*?Definition.*?(?:OldVal|Value).*?dE/dq.*?Step.*?(?:FinalVal|New-Value).*?"
            r"\n\s*-+\s*\n(.*?)(?:\n\s*-{40,}|\Z)",
            flags=re.S | re.I
        )
        
        all_matches = list(pattern.finditer(self.text))
        
        if all_matches:
            # Use the LAST match (optimized parameters at end of optimization)
            redundant_match = all_matches[-1]
            block = redundant_match.group(1).strip()
            lines = block.split('\n')
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('-') or len(line) < 10:
                    continue
                
                # Match pattern: "1. B(H   1,C   0)                1.0834 -0.000029  0.0000    1.0835"
                match_line = re.match(
                    r'^\s*(\d+)\.\s+([BAD])\(([^\)]+)\)\s+(' + _FLOAT_RE + r')\s+(' + 
                    _FLOAT_RE + r')\s+(' + _FLOAT_RE + r')\s+(' + _FLOAT_RE + r')\s*$',
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
                    
                    # Clean definition: "H   1,C   0" -> "H1-C0"
                    atoms_clean = re.sub(r'\s+', '', definition)
                    atoms_clean = atoms_clean.replace(',', '-')
                    
                    # Build entry
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
                        
                except (ValueError, IndexError, AttributeError) as e:
                    continue
            
            # If we got data from Redundant section, return it
            if out["bonds"] or out["angles"] or out["dihedrals"]:
                return {k: pd.DataFrame(v) for k, v in out.items()} if as_df else out
        
        # PRIORITY 2: Try Optimized Parameters block (less detailed fallback)
        m = _RE_OPT_PARAMS.search(self.text)
        if m:
            for L in m.group(1).strip().splitlines():
                dm = re.search(r"([BAD])\((.*?)\).*?([-+]?\d*\.\d+)$", L.strip())
                if not dm:
                    continue
                ctype = dm.group(1).upper()
                defn = dm.group(2).strip()
                val = _to_float(dm.group(3))
                
                # Simplified entry (no gradient info)
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
                return {k: pd.DataFrame(v) for k, v in out.items()} if as_df else out
        
        # PRIORITY 3: Basic Internal Coordinates block (last resort)
        b = _RE_INTERNAL_BLOCKS.search(self.text)
        if not b:
            return {k: pd.DataFrame() for k in out} if as_df else out
        
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
        
        return {k: pd.DataFrame(v) for k, v in out.items()} if as_df else out

    # ---------- orbitals ----------
    def parse_orbitals(self, as_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, Any]], None]:
        """Return DataFrame of last orbital block. Handles open and closed shell.
        If as_df=True returns pd.DataFrame, else list of dicts."""
        if not self.text:
            return pd.DataFrame() if as_df else None

        # Try spin-up/down blocks first (open-shell)
        up_blocks = _RE_SPIN_UP.findall(self.text)
        down_blocks = _RE_SPIN_DOWN.findall(self.text)

        spins = []
        if up_blocks and down_blocks:
            spins = [("up", up_blocks[-1]), ("down", down_blocks[-1])]
        else:
            # fallback: single ORBITAL ENERGIES block (closed-shell)
            blocks = _RE_ORBITAL_BLOCK.findall(self.text)
            if not blocks:
                return pd.DataFrame() if as_df else None
            spins = [("na", blocks[-1])]

        orbitals: List[Dict[str, Any]] = []
        for spin_label, block in spins:
            for L in block.strip().splitlines():
                parts = L.split()
                if not parts or not parts[0].lstrip().isdigit():
                    continue
                # Expect: NO OCC E(Eh) E(eV)
                if len(parts) < 4:
                    continue
                occ = _to_float(parts[1])
                eh = _to_float(parts[2])
                ev = _to_float(parts[3])
                if None in (occ, eh, ev):
                    continue
                # only valid occ values
                if occ not in (0.0, 1.0, 2.0):
                    continue
                # drop silly fake rows
                if eh == occ + 1.0 and ev == occ + 2.0:
                    continue
                orbitals.append({"OCC": occ, "Eh": eh, "eV": ev, "spin": spin_label})

        if not orbitals:
            return pd.DataFrame() if as_df else None

        # separate occupied & virtual, sort, assign levels
        occupied = [o for o in orbitals if o["OCC"] > 0.0]
        virtual = [o for o in orbitals if o["OCC"] == 0.0]

        # HOMO = highest Eh among occupied
        occupied.sort(key=lambda x: x["Eh"], reverse=True)
        virtual.sort(key=lambda x: x["Eh"])

        for i, o in enumerate(occupied):
            o["lvl"] = i  # 0 = HOMO, 1 = HOMO-1, ...
        for i, o in enumerate(virtual):
            o["lvl"] = i  # 0 = LUMO, 1 = LUMO+1, ...

        all_orbs = occupied + virtual
        all_orbs.sort(key=lambda x: x["Eh"])

        if as_df:
            return pd.DataFrame(all_orbs)
        return all_orbs

    # ---------- Build HOMO/LUMO label mapping ----------
    def _build_state_to_orbital_map(self) -> Dict[int, Tuple[int, int]]:
        """Map TDDFT state number to the dominant orbital transition (from_orb, to_orb).
        
        Returns dict like {1: (52, 54), 2: (53, 54), ...}
        where keys are state numbers and values are (from_orbital_idx, to_orbital_idx).
        """
        tddft = self.parse_tddft_states()
        if tddft.empty:
            return {}
        
        # For each state, find the transition with the largest weight
        state_orbs = {}
        for state_num in sorted(tddft['state'].unique()):
            state_data = tddft[tddft['state'] == state_num]
            # Sort by abs(weight) to find dominant contribution
            state_data = state_data.sort_values('weight', key=abs, ascending=False)
            if not state_data.empty:
                top = state_data.iloc[0]
                state_orbs[int(state_num)] = (int(top['from_orb']), int(top['to_orb']))
        
        return state_orbs

    def _build_orbital_homo_lumo_map(self) -> Dict[int, str]:
        """Build mapping from orbital index to HOMO/LUMO label using the existing 'lvl' field.
        
        Returns dict like {53: 'H', 54: 'L', 55: 'L+1', 52: 'H-1', ...}
        Uses the 'lvl' field from parse_orbitals which already has proper HOMO/LUMO level assignments.
        """
        orbs_df = self.parse_orbitals(as_df=True)
        if orbs_df.empty:
            return {}
        
        label_map = {}
        
        for idx, row in orbs_df.iterrows():
            lvl = int(row['lvl'])
            occ = row['OCC']
            
            if occ > 0.0:
                # Occupied orbital
                if lvl == 0:
                    label_map[idx] = 'H'
                else:
                    label_map[idx] = f'H-{lvl}'
            else:
                # Virtual orbital
                if lvl == 0:
                    label_map[idx] = 'L'
                else:
                    label_map[idx] = f'L+{lvl}'
        
        return label_map

    def _convert_state_label(self, state_label: str) -> str:
        """Convert ORCA state label to HOMO/LUMO notation.
        
        ORCA format: 'N-MS' where N = state number (0=ground), MS = multiplicity suffix (1A, 1B, etc.)
        
        Ground state (0-1A, 0-3A) → kept as is (or can be labeled "GS")
        Excited states (1-1A, 2-1A, etc.) → converted using dominant orbital transition
        
        Returns label in format 'H→L', 'H-1→L+1', etc.
        Or returns original label if conversion fails.
        """
        # Extract numeric part (state number)
        match = re.match(r'(\d+)([a-zA-Z0-9]*)', state_label.strip())
        if not match:
            return state_label
        
        state_num_str = match.group(1)
        try:
            state_num = int(state_num_str)
        except ValueError:
            return state_label
        
        # Ground state (0) - keep as is
        if state_num == 0:
            return state_label
        
        # Build maps on first call
        if not hasattr(self, '_state_orbital_cache'):
            self._state_orbital_cache = self._build_state_to_orbital_map()
            self._orbital_label_cache = self._build_orbital_homo_lumo_map()
        
        # Get the dominant orbital transition for this excited state
        orb_pair = self._state_orbital_cache.get(state_num)
        if orb_pair is None:
            return state_label
        
        from_orb_idx, to_orb_idx = orb_pair
        from_label = self._orbital_label_cache.get(from_orb_idx)
        to_label = self._orbital_label_cache.get(to_orb_idx)
        
        if from_label is None or to_label is None:
            return state_label
        
        # Return in format "H→L"
        return f"{from_label}→{to_label}"

    # ---------- vibrations ----------
    def parse_vibrations(self) -> Optional[List[Dict[str, Any]]]:
        if not self.text:
            return None
        blocks = _RE_VIBS.findall(self.text)
        if not blocks:
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
        return out if out else None

    # ---------- IR ----------
    def parse_ir_spectrum(self) -> Optional[List[Dict[str, Any]]]:
        if not self.text:
            return None
        m = _RE_IR.search(self.text)
        if not m:
            return None
        block = m.group(1).strip().splitlines()
        out = []
        for L in block:
            parts = L.split()
            if len(parts) < 4:
                continue
            if not parts[0].rstrip(":").isdigit():
                continue
            freq = _to_float(parts[1]); eps = _to_float(parts[2]); intensity = _to_float(parts[3])
            if None in (freq, eps, intensity):
                continue
            out.append({"freq_cm-1": freq, "eps": eps, "intensity_km/mol": intensity})
        return out if out else None

    # ---------- RAMAN ----------
    def parse_raman_spectrum(self, as_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, Any]], None]:
        if not self.text:
            return pd.DataFrame() if as_df else None
        m = _RE_RAMAN.search(self.text)
        if not m:
            return pd.DataFrame() if as_df else None
        block = m.group(1).strip().splitlines()
        out: List[Dict[str, Any]] = []
        for L in block:
            parts = re.split(r"\s+", L.strip())
            # Accept lines like: "  6:    44.92   4.325850   0.749978"
            # sometimes Mode index with colon or not
            if not parts or not re.match(r"^\d+:?$", parts[0]):
                # allow if first token is number without colon
                if not parts[0].isdigit():
                    continue
            # find numeric tokens after the label
            nums = [p for p in parts if re.match(r"^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$", p)]
            if len(nums) < 3:
                continue
            freq = _to_float(nums[0]); activity = _to_float(nums[1]); depol = _to_float(nums[2])
            if None in (freq, activity, depol):
                continue
            out.append({"freq_cm-1": freq, "activity": activity, "depolarization": depol})
        return pd.DataFrame(out) if as_df else out

    # ---------- Mulliken ----------
    def parse_mulliken(self, as_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, Any]]]:
        if not self.text:
            return pd.DataFrame() if as_df else []
        m = _RE_MULLIKEN.search(self.text)
        if not m:
            return pd.DataFrame() if as_df else []
        lines = m.group(1).strip().splitlines()
        out = []
        for L in lines:
            parts = re.split(r"\s+", L.strip())
            if len(parts) < 4:
                continue
            pop = _to_float(parts[2]); charge = _to_float(parts[3])
            if None in (pop, charge):
                continue
            out.append({"Nucleus": parts[0], "Element": parts[1], "Population": pop, "Charge": charge})
        return pd.DataFrame(out) if as_df else out

    # ---------- NMR ----------
    def parse_nmr(self, as_df: bool = True) -> Dict[str, Any]:
        out = {"shielding": pd.DataFrame(), "coupling": pd.DataFrame()} if as_df else {"shielding": [], "coupling": []}
        if not self.text:
            return out
        m = _RE_NMR_SHIELD.search(self.text)
        if m:
            rows = []
            for L in m.group(1).strip().splitlines():
                parts = re.split(r"\s+", L.strip())
                if len(parts) < 4:
                    continue
                iso = _to_float(parts[2]); aniso = _to_float(parts[3])
                if None in (iso, aniso):
                    continue
                rows.append({"Nucleus": parts[0], "Element": parts[1], "Isotropic": iso, "Anisotropy": aniso})
            if as_df:
                out["shielding"] = pd.DataFrame(rows) if rows else pd.DataFrame()
            else:
                out["shielding"] = rows
        m2 = _RE_NMR_COUPLING.search(self.text)
        if m2:
            rows = []
            lines = m2.group(1).strip().splitlines()
            header = None; data_start = 0
            for i, L in enumerate(lines):
                parts = re.split(r"\s+", L.strip())
                # crude heuristic: header contains uppercase nucleus labels like 'H', 'C'
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
            if as_df:
                out["coupling"] = pd.DataFrame(rows) if rows else pd.DataFrame()
            else:
                out["coupling"] = rows
        return out

    # ---------- energies ----------
    def parse_energies(self) -> Dict[str, Optional[float]]:
        if not self.text:
            return {"gibbs_Eh": None, "single_point_Eh": None}
        gib = _RE_FINAL_GIBBS.search(self.text)
        spe = _RE_FINAL_SPE.search(self.text)
        return {
            "gibbs_Eh": float(gib.group(1)) if gib else None,
            "single_point_Eh": float(spe.group(1)) if spe else None
        }

    # ---------- TDDFT states ----------
    def parse_tddft_states(self) -> pd.DataFrame:
        """Parse TD-DFT states (singlets & triplets). Return DataFrame with transitions."""
        if not self.text:
            return pd.DataFrame()
        states = []
        
        # Split text into sections to identify singlet vs triplet blocks
        # Look for section headers that indicate multiplicity
        singlet_sections = []
        triplet_sections = []
        
        # Find all TD-DFT section headers
        singlet_pattern = r'TD-DFT(?:/TDA)?\s+EXCITED\s+STATES\s+\(SINGLETS\)'
        triplet_pattern = r'TD-DFT(?:/TDA)?\s+EXCITED\s+STATES\s+\(TRIPLETS\)'
        
        # Find positions of singlet and triplet sections
        for m in re.finditer(singlet_pattern, self.text, flags=re.IGNORECASE):
            singlet_sections.append(m.start())
        
        for m in re.finditer(triplet_pattern, self.text, flags=re.IGNORECASE):
            triplet_sections.append(m.start())
        
        # Combine and sort all section positions
        all_sections = [(pos, 1) for pos in singlet_sections] + [(pos, 3) for pos in triplet_sections]
        all_sections.sort(key=lambda x: x[0])
        
        # collect from regex
        for st_m in _RE_TDDFT_STATE.finditer(self.text):
            stnum = int(st_m.group(1))
            energy_au = _to_float(st_m.group(2))
            energy_ev = _to_float(st_m.group(3))
            energy_cm1 = _to_float(st_m.group(4))
            block = st_m.group(5)
            
            # Determine multiplicity based on section position
            mult = 1  # default to singlet
            match_pos = st_m.start()
            
            # Find which section this state belongs to
            for i, (sec_pos, sec_mult) in enumerate(all_sections):
                if match_pos > sec_pos:
                    # Check if there's a next section
                    if i + 1 < len(all_sections):
                        next_pos = all_sections[i + 1][0]
                        if match_pos < next_pos:
                            mult = sec_mult
                            break
                    else:
                        # This is the last section
                        mult = sec_mult
                        break
            
            # Also check for explicit Mult declaration in header line (newer format)
            header_line = st_m.group(0).splitlines()[0]
            m_mult = re.search(r"Mult\s+(\d+)", header_line, flags=re.I)
            if m_mult:
                try:
                    mult = int(m_mult.group(1))
                except Exception:
                    pass
            
            for tr in _RE_TD_TRANS.finditer(block):
                frm = int(tr.group(1)); to = int(tr.group(2))
                weight = _to_float(tr.group(3)); coeff = _to_float(tr.group(4))
                states.append({
                    "state": stnum,
                    "energy_au": energy_au,
                    "energy_ev": energy_ev,
                    "energy_cm1": energy_cm1,
                    "from_orb": frm,
                    "to_orb": to,
                    "weight": weight,
                    "coeff": coeff,
                    "mult": mult
                })
        
        # Calculate HOMO/LUMO levels for each transition
        # HOMO = max(from_orbs), LUMO = min(to_orbs)
        if states:
            from_orbs = set(s["from_orb"] for s in states)
            to_orbs = set(s["to_orb"] for s in states)
            homo_idx = max(from_orbs) if from_orbs else None
            lumo_idx = min(to_orbs) if to_orbs else None
            
            if homo_idx is not None and lumo_idx is not None:
                for state in states:
                    from_orb = state["from_orb"]
                    to_orb = state["to_orb"]
                    # Calculate HOMO/LUMO level numbers (for display purposes)
                    # HOMO = 0, HOMO-1 = -1, HOMO-2 = -2, etc.
                    # LUMO = 0, LUMO+1 = 1, LUMO+2 = 2, etc.
                    state["homo_level"] = from_orb - homo_idx
                    state["lumo_level"] = to_orb - lumo_idx
        
        return pd.DataFrame(states)
    
    # ---------- electric dipole spectrum ----------
    def parse_electric_dipole_spectrum(self) -> Dict[str, pd.DataFrame]:
        if not self.text:
            return {"abs": pd.DataFrame(), "soc": pd.DataFrame()}

        pat = re.compile(
            r"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*"
            r"-+\s*(.*?)\n\s*-{10,}",
            re.S | re.I
        )

        abs_df = self._parse_abs_block(pat)

        soc_pattern = re.compile(
            r"SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS.*?"
            r"State\s+Energy.*?\n-+\n(.*?)(?:\n{2,}|\Z)",
            flags=re.S | re.I
        )
        soc_df = self._parse_abs_block(soc_pattern)

        return {"abs": abs_df, "soc": soc_df}

    # ---------- velocity dipole spectrum ----------
    def parse_velocity_dipole_spectrum(self) -> Dict[str, pd.DataFrame]:
        if not self.text:
            return {"abs": pd.DataFrame(), "soc": pd.DataFrame()}

        pat = re.compile(
            r"ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*"
            r"-+\s*(.*?)\n\s*-{10,}",
            re.S | re.I
        )

        abs_df = self._parse_abs_block(pat)

        soc_pattern = re.compile(
            r"SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS.*?"
            r"State\s+Energy.*?\n-+\n(.*?)(?:\n{2,}|\Z)",
            flags=re.S | re.I
        )
        soc_df = self._parse_abs_block(soc_pattern)

        return {"abs": abs_df, "soc": soc_df}

    # ---------- absorption / velocity ----------
    def _parse_abs_block(self, pattern: re.Pattern) -> pd.DataFrame:
        """Parse absorption spectrum block (electric or velocity dipole).

        Corrected implementation: locate the section header and extract data
        lines after the second dashed separator. Returns a DataFrame of rows
        with fields like from_state, to_state, energy_ev, wavelength_nm, fosc,
        dipole_strength and dx/dy/dz.
        """
        if not self.text:
            return pd.DataFrame()

        # Determine which section to find based on the pattern
        pat_text = pattern.pattern.upper()
        if "ELECTRIC DIPOLE" in pat_text:
            section_start = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
        elif "VELOCITY DIPOLE" in pat_text:
            section_start = "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS"
        else:
            return pd.DataFrame()

        # Locate the section start
        start_idx = self.text.find(section_start)
        if start_idx == -1:
            return pd.DataFrame()

        # Extract a chunk of the section (header + data)
        section_text = self.text[start_idx:start_idx + 10000]
        lines = section_text.split('\n')

        rows = []
        dash_count = 0

        for line in lines:
            s = line.strip()

            # count dashed separator lines (---)
            if s.startswith('-') and set(s) <= {'-', ' '}:
                dash_count += 1
                continue

            # Skip empty lines (don't stop parsing)
            if not s:
                continue

            # Data begins after the second dashed separator
            if dash_count >= 2:
                # stop if another major section appears (only truly different section)
                if 'CD SPECTRUM' in s.upper():
                    break

                if '->' in s:
                    parts = s.split()
                    if len(parts) >= 10:
                        try:
                            from_state_raw = parts[0]
                            to_state_raw = parts[2]
                            
                            # Extract state numbers from raw labels (e.g., "1" from "1-1A")
                            # These are for TDDFT state matching
                            from_state_match = re.match(r'(\d+)', from_state_raw)
                            to_state_match = re.match(r'(\d+)', to_state_raw)
                            from_state_num = int(from_state_match.group(1)) if from_state_match else None
                            to_state_num = int(to_state_match.group(1)) if to_state_match else None
                            
                            # Convert to HOMO/LUMO notation for display
                            from_state_homo_lumo = self._convert_state_label(from_state_raw)
                            to_state_homo_lumo = self._convert_state_label(to_state_raw)
                            
                            rows.append({
                                "from_state": from_state_raw,  # Keep original for TDDFT matching
                                "to_state": to_state_raw,      # Keep original for TDDFT matching
                                "from_state_homo_lumo": from_state_homo_lumo,  # HOMO/LUMO display
                                "to_state_homo_lumo": to_state_homo_lumo,      # HOMO/LUMO display
                                "energy_ev": float(parts[3]),
                                "energy_cm": float(parts[4]),
                                "wavelength_nm": float(parts[5]),
                                "fosc": float(parts[6]),
                                "dipole_strength": float(parts[7]),
                                "dx": float(parts[8]),
                                "dy": float(parts[9]),
                                "dz": float(parts[10]) if len(parts) > 10 else None
                            })
                        except (ValueError, IndexError):
                            # skip malformed lines
                            continue

        return pd.DataFrame(rows)

    # ---------- spectrum file parser ----------
    def parse_spectrum_file(self, filename: str, spectrum_type: Optional[str] = None, as_df: bool = True) -> Union[pd.DataFrame, List[Dict[str, Any]]]:
        try:
            with open(filename, "r") as fh:
                lines = fh.readlines()
        except Exception:
            return pd.DataFrame() if as_df else []
        data = []
        for L in lines[1:]:
            parts = re.split(r"\s+", L.strip())
            if not parts:
                continue
            energy = _to_float(parts[0])
            if energy is None:
                continue
            if spectrum_type and spectrum_type.upper() in ("VG", "AH", "AHAS"):
                if len(parts) < 4:
                    continue
                tot = _to_float(parts[1]); fc = _to_float(parts[2]); ht = _to_float(parts[3])
                if None in (tot, fc, ht):
                    continue
                data.append({"energy_cm1": energy, "total_spectrum": tot, "intensity_fc": fc, "intensity_ht": ht})
            elif spectrum_type and spectrum_type.upper() in ("FLUOR", "PHOSP"):
                if len(parts) < 2:
                    continue
                tot = _to_float(parts[1])
                data.append({"energy_cm1": energy, "total_spectrum": tot})
            else:
                nums = [ _to_float(x) for x in parts[1:] ]
                entry = {"energy_cm1": energy}
                for i, v in enumerate(nums):
                    entry[f"col{i+1}"] = v
                data.append(entry)
        return pd.DataFrame(data) if as_df else data

    # ---------- parse esd flag (VG/AH/AHAS/FLUOR/PHOSP/...) ----------
    def parse_esd_flag(self) -> Optional[str]:
        if not self.text:
            return None
        m = _RE_INPUT_FILE.search(self.text)
        if m:
            s = m.group(2)
            # look for ESD or HESSFLAG keywords
            m1 = re.search(r"ESD\s*\(\s*([A-Za-z]+)\s*\)", s, flags=re.I)
            if m1:
                return m1.group(1).upper()
            m2 = re.search(r"HESSFLAG\s+(\w+)", s, flags=re.I)
            if m2:
                return m2.group(1).upper()
            # fallback: DOHT/HESSFLAG/ESD presence
            if re.search(r"HESSFLAG\s+VG", s, flags=re.I):
                return "VG"
            if re.search(r"HESSFLAG\s+AHAS", s, flags=re.I):
                return "AHAS"
            if re.search(r"ESD\s*\(\s*ABS", s, flags=re.I):
                return "ABS"
        # fallback to scanning file text for HESSFLAG lines
        if re.search(r"HESSFLAG\s+VG", self.text, flags=re.I):
            return "VG"
        if re.search(r"HESSFLAG\s+AHAS", self.text, flags=re.I):
            return "AHAS"
        return None

    # ---------- detect calc type ----------
    def detect_calc_type(self) -> Dict[str, Any]:
        out = {"is_optimization": False, "has_tddft": False, "has_esd": False, "esd_type": None, "optimized_state": "S0", "calc_class": "single_point"}
        if not self.text:
            return out
        if re.search(r"GEOMETRY OPTIMIZATION", self.text, flags=re.I) or re.search(r"OPTIMIZATION CONVERGED", self.text, flags=re.I):
            out["is_optimization"] = True
            out["calc_class"] = "optimization"
        if re.search(r"TD-DFT", self.text, flags=re.I):
            out["has_tddft"] = True
            if out["calc_class"] == "single_point":
                out["calc_class"] = "tddft"
        esd = self.parse_esd_flag()
        if esd:
            out["has_esd"] = True
            out["esd_type"] = esd
            if out["calc_class"] == "single_point":
                out["calc_class"] = "spectrum"
        # find IROOT in %TDDFT block to infer optimized_state root
        tblock = re.search(r"%TDDFT\s*(.*?)END", self.text, flags=re.S|re.I)
        if tblock:
            mroot = re.search(r"IROOT\s+(\d+)", tblock.group(1), flags=re.I)
            if mroot:
                try:
                    root = int(mroot.group(1))
                    out["optimized_state"] = f"S{root}"
                except Exception:
                    pass
        # multiplicity inference from input
        in_m = _RE_INPUT_FILE.search(self.text)
        if in_m:
            mxyz = re.search(r"\*\s*xyzfile\s+(\d+)\s+(\d+)", in_m.group(2), flags=re.I)
            if mxyz:
                try:
                    mult = int(mxyz.group(2))
                    if mult == 3 and out["optimized_state"] == "S0":
                        out["optimized_state"] = "T1"
                except Exception:
                    pass
        return out

    # ---------- master parse ----------
    def parse(self, as_df: bool = True) -> Dict[str, Any]:
        """Return parsed dictionary. If as_df True common fields are DataFrames where possible."""
        if not self.text:
            # empty structure
            return {
                "filename": self.filename,
                "geometry": None,
                "cart_coords": pd.DataFrame() if as_df else None,
                "smiles": None,
                "internal": {"bonds": pd.DataFrame(), "angles": pd.DataFrame(), "dihedrals": pd.DataFrame()} if as_df else None,
                "orbitals": {},
                "vibrations": pd.DataFrame() if as_df else None,
                "ir_spectrum": pd.DataFrame() if as_df else None,
                "raman_spectrum": pd.DataFrame() if as_df else None,
                "gibbs_energy_Eh": None,
                "single_point_energy_Eh": None,
                "tddft_states": pd.DataFrame() if as_df else None,
                "electric_dipole_spectrum": {"abs": pd.DataFrame(), "soc": pd.DataFrame()} if as_df else None,
                "velocity_dipole_spectrum": {"abs": pd.DataFrame(), "soc": pd.DataFrame()} if as_df else None,
                "nmr_data": {"shielding": pd.DataFrame(), "coupling": pd.DataFrame()} if as_df else None,
                "mulliken_charges": pd.DataFrame() if as_df else None,
                "is_optimization": False,
                "optimized_state": "S0",
                "esd_type": None,
                "calc_class": "single_point",
                "spectrum_file": {}
            }

        geometry = self.parse_geometry_info()
        cart_coords = self.parse_last_cartesian()
        smiles = self.coords_to_smiles(cart_coords) if cart_coords else None
        cart_out = pd.DataFrame(cart_coords, columns=["atom","x","y","z"]) if as_df and cart_coords else (cart_coords if cart_coords else (pd.DataFrame() if as_df else None))

        energies = self.parse_energies()
        calc_info = self.detect_calc_type()
        tddft = self.parse_tddft_states()
        electric = self.parse_electric_dipole_spectrum()
        velocity = self.parse_velocity_dipole_spectrum()
        raman = self.parse_raman_spectrum(as_df=as_df)
        nmr = self.parse_nmr(as_df=as_df)
        mull = self.parse_mulliken(as_df=as_df)
        internal = self.parse_internal(as_df=as_df)
        orbitals_df = self.parse_orbitals(as_df=as_df)
        vibrations = self.parse_vibrations()
        ir = self.parse_ir_spectrum()

        # Spectrum file heuristics
        spec_files = {}
        base_name = os.path.splitext(self.filename)[0]
        esd_flag = self.parse_esd_flag() or calc_info.get("esd_type")
        candidates = [f"{base_name}.spectrum"]
        if esd_flag:
            candidates.append(f"{base_name}_{esd_flag.lower()}.spectrum")
        if geometry and geometry.get("filename"):
            candidates.append(f"{geometry['filename']}.spectrum")
            if esd_flag:
                candidates.append(f"{geometry['filename']}_{esd_flag.lower()}.spectrum")
        for cand in candidates:
            if os.path.exists(cand):
                sp = self.parse_spectrum_file(cand, spectrum_type=esd_flag, as_df=as_df)
                if (isinstance(sp, pd.DataFrame) and not sp.empty) or (isinstance(sp, dict) and sp):
                    spec_files[esd_flag or "UNKNOWN"] = sp
                    break

        # Pack results
        result = {
            "filename": self.filename,
            "geometry": geometry,
            "cart_coords": cart_out,
            "smiles": smiles,
            "internal": internal,
            "orbitals": {},  # Option A: dict of DataFrames keyed by state
            "vibrations": vibrations if vibrations else (pd.DataFrame() if as_df else None),
            "ir_spectrum": ir if ir else (pd.DataFrame() if as_df else None),
            "raman_spectrum": raman if raman is not None else (pd.DataFrame() if as_df else None),
            "gibbs_energy_Eh": energies.get("gibbs_Eh"),
            "single_point_energy_Eh": energies.get("single_point_Eh"),
            "tddft_states": tddft if not tddft.empty else (pd.DataFrame() if as_df else None),
            "electric_dipole_spectrum": electric,
            "velocity_dipole_spectrum": velocity,
            "nmr_data": nmr,
            "mulliken_charges": mull,
            "spectrum_file": spec_files,
            "is_optimization": calc_info.get("is_optimization"),
            "optimized_state": calc_info.get("optimized_state"),
            "esd_type": calc_info.get("esd_type"),
            "calc_class": calc_info.get("calc_class")
        }

        # Decide state key: optimized_state + OPT/SP marker
        state_key = f"{result['optimized_state']}_{'OPT' if result['is_optimization'] else 'SP'}"
        # store orbitals under this key if present
        if orbitals_df is not None and (not isinstance(orbitals_df, pd.DataFrame) or not orbitals_df.empty):
            result["orbitals"][state_key] = orbitals_df if as_df else orbitals_df.to_dict("records")

        # store energies under energies dict (for easier grouping later)
        result["energies"] = {}
        if result["gibbs_energy_Eh"] is not None or result["single_point_energy_Eh"] is not None:
            e_dict = {}
            if result["gibbs_energy_Eh"] is not None:
                e_dict["gibbs_Eh"] = result["gibbs_energy_Eh"]
            if result["single_point_energy_Eh"] is not None:
                e_dict["single_point_Eh"] = result["single_point_energy_Eh"]
            result["energies"][state_key] = e_dict

        return result

# -----------------------
# Batch parser
# -----------------------

class ORCABatchParser:
    """Batch parser that accepts multiple glob patterns.
    Produces a DataFrame where each row is a parsed calculation dict.
    Then groups by molecule_id to produce a standardized record per molecule.
    """

    def __init__(self, *patterns: str):
        # accept multiple patterns
        if not patterns:
            patterns = ("*.out",)
        self.patterns = patterns
        self.files = sorted({f for pat in patterns for f in glob.glob(pat)})
        self.states: Dict[str, int] = {}

    @staticmethod
    def _extract_molecule_id(filename: str, geometry_info: Optional[Dict[str, Any]]) -> str:
        if geometry_info and geometry_info.get("filename"):
            geom_file = geometry_info["filename"]
            mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt)$', '', geom_file, flags=re.I)
            return mol_id.strip("_-")
        base = os.path.splitext(os.path.basename(filename))[0]
        mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt)$', '', base, flags=re.I)
        return mol_id.strip("_-")

    def _is_empty(self, val: Any) -> bool:
        return _is_empty_obj(val)

    def parse_all(self, verbose: bool = True, limit_files: Optional[int] = None) -> pd.DataFrame:
        parsed = []
        files = self.files[:limit_files] if limit_files else self.files
        if verbose:
            print(f"Found files: {len(files)}")
        for f in files:
            try:
                p = ORCAParser(f)
                res = p.parse(as_df=True)
                res["molecule_id"] = self._extract_molecule_id(f, res.get("geometry"))
                parsed.append(res)
                if verbose:
                    print(f"Parsed: {os.path.basename(f)} -> {res['molecule_id']}")
            except Exception as ex:
                if verbose:
                    print(f"Error parsing {f}: {ex}")

        if not parsed:
            return pd.DataFrame()

        raw = pd.DataFrame(parsed)

        # Propagate geometry-dependent fields from OPT to SP of same molecule/state where missing
        geometry_fields = [
            "cart_coords", "smiles", "geometry", "vibrations", "ir_spectrum",
            "raman_spectrum", "nmr_data", "mulliken_charges", "internal"
        ]
        # energies & orbitals are state-specific and normally not propagated
        filled = raw.copy()

        for mid in raw["molecule_id"].unique():
            mask_mid = raw["molecule_id"] == mid
            states = raw[mask_mid]["optimized_state"].unique()
            for st in states:
                opt_mask = mask_mid & (raw["is_optimization"] == True) & (raw["optimized_state"] == st)
                if not raw[opt_mask].empty:
                    opt_row = raw[opt_mask].iloc[0]
                    non_opt_mask = mask_mid & (raw["is_optimization"] == False)
                    for idx in raw[non_opt_mask].index:
                        for col in geometry_fields:
                            if col not in filled.columns:
                                continue
                            try:
                                if (self._is_empty(filled.at[idx, col]) or pd.isna(filled.at[idx, col])) and not self._is_empty(opt_row[col]):
                                    filled.at[idx, col] = opt_row[col]
                            except Exception:
                                # cope with different types
                                try:
                                    if (filled.at[idx, col] is None or (isinstance(filled.at[idx, col], float) and pd.isna(filled.at[idx, col]))) and not self._is_empty(opt_row[col]):
                                        filled.at[idx, col] = opt_row[col]
                                except Exception:
                                    pass

        # Build standardized molecule records grouped by molecule_id
        standardized = []
        for mid in filled["molecule_id"].unique():
            mol_rows = filled[filled["molecule_id"] == mid]
            row = {
                "molecule_id": mid,
                "energies": {},
                "orbitals": {},
                "geometries": {},
                "internal_coords": {},
                # REMOVE: "absorption": {},
                "smiles": None,
                "cart_coords": pd.DataFrame(),
                "geometry": None,
                "vibrations": None,
                "ir_spectrum": None,
                "raman_spectrum": None,
                "nmr_data": None,
                "mulliken_charges": None,
                "tddft_states": pd.DataFrame(),
                # ADD THESE:
                "electric_dipole_spectrum": {"abs": pd.DataFrame(), "soc": pd.DataFrame()},
                "velocity_dipole_spectrum": {"abs": pd.DataFrame(), "soc": pd.DataFrame()},
                "spectrum_file": {}
            }

            # Prefer S0 OPT, but fallback to S0 SP if OPT not available
            s0_opt = mol_rows[(mol_rows["is_optimization"] == True) & (mol_rows["optimized_state"] == "S0")]
            s0_sp = mol_rows[(mol_rows["is_optimization"] == False) & (mol_rows["optimized_state"] == "S0")]
            s0_all = mol_rows[mol_rows["optimized_state"] == "S0"]

            # Use OPT if available, otherwise SP
            s0_source = s0_opt if not s0_opt.empty else s0_sp

            if not s0_source.empty:
                s0 = s0_source.iloc[0]
                row["smiles"] = s0.get("smiles")
                row["cart_coords"] = s0.get("cart_coords")
                row["geometry"] = s0.get("geometry")
                row["vibrations"] = s0.get("vibrations")
                row["ir_spectrum"] = s0.get("ir_spectrum")
                row["nmr_data"] = s0.get("nmr_data")
                row["mulliken_charges"] = s0.get("mulliken_charges")
                row["internal"] = s0.get("internal")
                row["raman_spectrum"] = s0.get("raman_spectrum")
            
            # Special handling for Raman: search ALL S0 calculations if not found yet
            # This is needed because Raman might only be in SP calculations
            if self._is_empty(row.get("raman_spectrum")):
                for _, calc in s0_all.iterrows():
                    raman = calc.get("raman_spectrum")
                    if not self._is_empty(raman):
                        row["raman_spectrum"] = raman
                        break

            # iterate all calcs for this molecule
            for idx, calc in mol_rows.iterrows():
                state = calc["optimized_state"]
                key_type = "OPT" if calc["is_optimization"] else "SP"
                key = f"{state}_{key_type}"

                # energies
                gibbs = calc.get("gibbs_energy_Eh")
                spe = calc.get("single_point_energy_Eh")
                if pd.notna(gibbs) or pd.notna(spe):
                    d = {}
                    if pd.notna(gibbs):
                        d["gibbs_Eh"] = float(gibbs)
                    if pd.notna(spe):
                        d["single_point_Eh"] = float(spe)
                    if d:
                        row["energies"][key] = d

                # orbitals
                orbs = calc.get("orbitals", {})
                if isinstance(orbs, dict) and orbs:
                    # orbs is a dict keyed by state (usually only this state)
                    for k_orb, v_orb in orbs.items():
                        # if DataFrame -> keep as is
                        row["orbitals"][k_orb] = v_orb

                # geometries
                geom = calc.get("cart_coords")
                if not self._is_empty(geom):
                    row["geometries"][key] = geom

                # internal coords
                internal = calc.get("internal")
                if not self._is_empty(internal):
                    row["internal_coords"][key] = internal

                # tddft
                tdd = calc.get("tddft_states")
                if isinstance(tdd, pd.DataFrame) and not tdd.empty:
                    # if already have tddft states, append
                    if isinstance(row.get("tddft_states"), pd.DataFrame) and not row["tddft_states"].empty:
                        try:
                            row["tddft_states"] = pd.concat([row["tddft_states"], tdd], ignore_index=True)
                        except Exception:
                            row["tddft_states"] = tdd
                    else:
                        row["tddft_states"] = tdd

                # electric dipole spectrum
                elec = calc.get("electric_dipole_spectrum")
                if elec and not self._is_empty(elec):
                    for key in ["abs", "soc"]:
                        if key in elec and not self._is_empty(elec[key]):
                            if self._is_empty(row["electric_dipole_spectrum"][key]):
                                row["electric_dipole_spectrum"][key] = elec[key]
                            else:
                                try:
                                    row["electric_dipole_spectrum"][key] = pd.concat(
                                        [row["electric_dipole_spectrum"][key], elec[key]], 
                                        ignore_index=True
                                    )
                                except Exception:
                                    row["electric_dipole_spectrum"][key] = elec[key]

                # velocity dipole spectrum
                vel = calc.get("velocity_dipole_spectrum")
                if vel and not self._is_empty(vel):
                    for key in ["abs", "soc"]:
                        if key in vel and not self._is_empty(vel[key]):
                            if self._is_empty(row["velocity_dipole_spectrum"][key]):
                                row["velocity_dipole_spectrum"][key] = vel[key]
                            else:
                                try:
                                    row["velocity_dipole_spectrum"][key] = pd.concat(
                                        [row["velocity_dipole_spectrum"][key], vel[key]], 
                                        ignore_index=True
                                    )
                                except Exception:
                                    row["velocity_dipole_spectrum"][key] = vel[key]

                # absorption / spectrum file
                spf = calc.get("spectrum_file")
                if spf and not self._is_empty(spf) and not row["spectrum_file"]:
                    row["spectrum_file"] = spf

            standardized.append(row)

        df = pd.DataFrame(standardized)
        self.states = {row["molecule_id"]: idx for idx, row in df.iterrows()}
        if verbose:
            print(f"Created {len(df)} molecule records from {len(files)} files")
        return df

# -----------------------
# Convenience functions
# -----------------------

def parse_orca_file(filename: str, as_df: bool = True) -> Dict[str, Any]:
    p = ORCAParser(filename)
    return p.parse(as_df=as_df)

def parse_orca_batch(*patterns: str, verbose: bool = True, limit_files: Optional[int] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    batch = ORCABatchParser(*patterns)
    df = batch.parse_all(verbose=verbose, limit_files=limit_files)
    return df, batch.states

# End of module