"""
TD-DFT Parser - Refactored from orca_praser.py

Parses: TD-DFT states, electric/velocity dipole spectra, state->orbital mapping
"""

import re
from typing import Optional, List, Dict, Any, Tuple
import pandas as pd

from .regex_patterns import RE_TDDFT_STATE, RE_TD_TRANS
from ..core.base_parser import BaseParser
from ..core.data_models import TDDFTData
from ..logger import get_logger


def _to_float(x: Any) -> Optional[float]:
    """Safe float conversion."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


class TDDFTParser(BaseParser):
    """Parser for TD-DFT data."""
    
    def __init__(self, text: str):
        super().__init__(text)
        self.logger = get_logger("TDDFTParser")
        self._state_orbital_cache = None
        self._orbital_label_cache = None
        self._orbitals_df = None
    
    def set_orbitals(self, orbitals_df: pd.DataFrame):
        """Set orbital data for HOMO/LUMO label conversion."""
        self._orbitals_df = orbitals_df
    
    def parse(self) -> TDDFTData:
        """Parse all TD-DFT data."""
        self.logger.debug("Starting TD-DFT parsing...")
        
        states = self.parse_tddft_states()
        elec = self.parse_electric_dipole_spectrum()
        vel = self.parse_velocity_dipole_spectrum()
        
        return TDDFTData(
            states=states,
            electric_dipole_abs=elec.get("abs"),
            electric_dipole_soc=elec.get("soc"),
            velocity_dipole_abs=vel.get("abs"),
            velocity_dipole_soc=vel.get("soc")
        )
    
    def parse_tddft_states(self) -> pd.DataFrame:
        """Parse TD-DFT states (singlets & triplets)."""
        self.logger.debug("Parsing TD-DFT excited states...")
        
        if not self.text:
            return pd.DataFrame()
        
        states = []
        
        # Find singlet/triplet section positions
        singlet_pattern = r'TD-DFT(?:/TDA)?\s+EXCITED\s+STATES\s+\(SINGLETS\)'
        triplet_pattern = r'TD-DFT(?:/TDA)?\s+EXCITED\s+STATES\s+\(TRIPLETS\)'
        
        singlet_sections = [m.start() for m in re.finditer(singlet_pattern, self.text, flags=re.IGNORECASE)]
        triplet_sections = [m.start() for m in re.finditer(triplet_pattern, self.text, flags=re.IGNORECASE)]
        
        if singlet_sections:
            self.logger.debug(f"  Found {len(singlet_sections)} singlet section(s)")
        if triplet_sections:
            self.logger.debug(f"  Found {len(triplet_sections)} triplet section(s)")
        
        all_sections = [(pos, 1) for pos in singlet_sections] + [(pos, 3) for pos in triplet_sections]
        all_sections.sort(key=lambda x: x[0])
        
        # Parse states
        for st_m in RE_TDDFT_STATE.finditer(self.text):
            stnum = int(st_m.group(1))
            energy_au = _to_float(st_m.group(2))
            energy_ev = _to_float(st_m.group(3))
            energy_cm1 = _to_float(st_m.group(4))
            block = st_m.group(5)
            
            # Determine multiplicity
            mult = 1
            match_pos = st_m.start()
            
            for i, (sec_pos, sec_mult) in enumerate(all_sections):
                if match_pos > sec_pos:
                    if i + 1 < len(all_sections):
                        next_pos = all_sections[i + 1][0]
                        if match_pos < next_pos:
                            mult = sec_mult
                            break
                    else:
                        mult = sec_mult
                        break
            
            # Check for explicit Mult declaration
            header_line = st_m.group(0).splitlines()[0]
            m_mult = re.search(r"Mult\s+(\d+)", header_line, flags=re.I)
            if m_mult:
                try:
                    mult = int(m_mult.group(1))
                except Exception:
                    pass
            
            # Parse transitions
            for tr in RE_TD_TRANS.finditer(block):
                frm = int(tr.group(1))
                to = int(tr.group(2))
                weight = _to_float(tr.group(3))
                coeff = _to_float(tr.group(4))
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
        
        # Calculate HOMO/LUMO levels
        if states:
            from_orbs = set(s["from_orb"] for s in states)
            to_orbs = set(s["to_orb"] for s in states)
            homo_idx = max(from_orbs) if from_orbs else None
            lumo_idx = min(to_orbs) if to_orbs else None
            
            if homo_idx is not None and lumo_idx is not None:
                for state in states:
                    state["homo_level"] = state["from_orb"] - homo_idx
                    state["lumo_level"] = state["to_orb"] - lumo_idx
            
            # Log summary
            unique_states = len(set(s["state"] for s in states))
            singlet_count = sum(1 for s in states if s["mult"] == 1)
            triplet_count = sum(1 for s in states if s["mult"] == 3)
            self.logger.info(f"  TD-DFT states: {unique_states} states, {len(states)} transitions")
            self.logger.info(f"    Singlets: {singlet_count}, Triplets: {triplet_count}")
        else:
            self.logger.debug("  No TD-DFT states found")
        
        return pd.DataFrame(states)
    
    def _build_state_to_orbital_map(self) -> Dict[int, Tuple[int, int]]:
        """Map TDDFT state number to the dominant orbital transition (from_orb, to_orb)."""
        tddft = self.parse_tddft_states()
        if tddft.empty:
            return {}
        
        state_orbs = {}
        for state_num in sorted(tddft['state'].unique()):
            state_data = tddft[tddft['state'] == state_num]
            state_data = state_data.sort_values('weight', key=abs, ascending=False)
            if not state_data.empty:
                top = state_data.iloc[0]
                state_orbs[int(state_num)] = (int(top['from_orb']), int(top['to_orb']))
        
        return state_orbs
    
    def _build_orbital_homo_lumo_map(self) -> Dict[int, str]:
        """Build mapping from orbital index to HOMO/LUMO label."""
        if self._orbitals_df is None or self._orbitals_df.empty:
            return {}
        
        label_map = {}
        for idx, row in self._orbitals_df.iterrows():
            lvl = int(row['lvl'])
            occ = row['OCC']
            
            if occ > 0.0:
                if lvl == 0:
                    label_map[idx] = 'H'
                else:
                    label_map[idx] = f'H-{lvl}'
            else:
                if lvl == 0:
                    label_map[idx] = 'L'
                else:
                    label_map[idx] = f'L+{lvl}'
        
        return label_map
    
    def _convert_state_label(self, state_label: str) -> str:
        """Convert ORCA state label to HOMO/LUMO notation."""
        match = re.match(r'(\d+)([a-zA-Z0-9]*)', state_label.strip())
        if not match:
            return state_label
        
        try:
            state_num = int(match.group(1))
        except ValueError:
            return state_label
        
        # Ground state (0) - keep as is
        if state_num == 0:
            return state_label
        
        # Build maps on first call
        if self._state_orbital_cache is None:
            self._state_orbital_cache = self._build_state_to_orbital_map()
            self._orbital_label_cache = self._build_orbital_homo_lumo_map()
        
        orb_pair = self._state_orbital_cache.get(state_num)
        if orb_pair is None:
            return state_label
        
        from_orb_idx, to_orb_idx = orb_pair
        from_label = self._orbital_label_cache.get(from_orb_idx)
        to_label = self._orbital_label_cache.get(to_orb_idx)
        
        if from_label is None or to_label is None:
            return state_label
        
        return f"{from_label}→{to_label}"
    
    def parse_electric_dipole_spectrum(self) -> Dict[str, pd.DataFrame]:
        """Parse electric dipole absorption spectrum."""
        self.logger.debug("Parsing electric dipole spectrum...")
        
        if not self.text:
            return {"abs": pd.DataFrame(), "soc": pd.DataFrame()}
        
        abs_df = self._parse_abs_block("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS")
        
        if not abs_df.empty:
            self.logger.info(f"  Electric dipole absorption: {len(abs_df)} transitions")
        
        soc_pattern = re.compile(
            r"SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS.*?"
            r"State\s+Energy.*?\n-+\n(.*?)(?:\n{2,}|\Z)",
            flags=re.S | re.I
        )
        soc_df = self._parse_soc_block(soc_pattern)
        
        if not soc_df.empty:
            self.logger.info(f"  Electric dipole SOC: {len(soc_df)} transitions")
        
        return {"abs": abs_df, "soc": soc_df}
    
    def parse_velocity_dipole_spectrum(self) -> Dict[str, pd.DataFrame]:
        """Parse velocity dipole absorption spectrum."""
        self.logger.debug("Parsing velocity dipole spectrum...")
        
        if not self.text:
            return {"abs": pd.DataFrame(), "soc": pd.DataFrame()}
        
        abs_df = self._parse_abs_block("ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS")
        
        if not abs_df.empty:
            self.logger.info(f"  Velocity dipole absorption: {len(abs_df)} transitions")
        
        soc_pattern = re.compile(
            r"SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS.*?"
            r"State\s+Energy.*?\n-+\n(.*?)(?:\n{2,}|\Z)",
            flags=re.S | re.I
        )
        soc_df = self._parse_soc_block(soc_pattern)
        
        if not soc_df.empty:
            self.logger.info(f"  Velocity dipole SOC: {len(soc_df)} transitions")
        
        return {"abs": abs_df, "soc": soc_df}
    
    def _parse_abs_block(self, section_start: str) -> pd.DataFrame:
        """Parse absorption spectrum block with HOMO/LUMO conversion."""
        if not self.text:
            return pd.DataFrame()
        
        start_idx = self.text.find(section_start)
        if start_idx == -1:
            return pd.DataFrame()
        
        section_text = self.text[start_idx:start_idx + 10000]
        lines = section_text.split('\n')
        
        rows = []
        dash_count = 0
        
        for line in lines:
            s = line.strip()
            
            if s.startswith('-') and set(s) <= {'-', ' '}:
                dash_count += 1
                continue
            
            if not s:
                continue
            
            if dash_count >= 2:
                if 'CD SPECTRUM' in s.upper():
                    break
                
                if '->' in s:
                    parts = s.split()
                    if len(parts) >= 10:
                        try:
                            from_state_raw = parts[0]
                            to_state_raw = parts[2]
                            
                            # Convert to HOMO/LUMO notation
                            from_state_homo_lumo = self._convert_state_label(from_state_raw)
                            to_state_homo_lumo = self._convert_state_label(to_state_raw)
                            
                            rows.append({
                                "from_state": from_state_raw,
                                "to_state": to_state_raw,
                                "from_state_homo_lumo": from_state_homo_lumo,
                                "to_state_homo_lumo": to_state_homo_lumo,
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
                            continue
        
        return pd.DataFrame(rows)
    
    def _parse_soc_block(self, pattern: re.Pattern) -> pd.DataFrame:
        """Parse SOC block."""
        m = pattern.search(self.text)
        if not m:
            return pd.DataFrame()
        
        rows = []
        for line in m.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 6:
                try:
                    rows.append({
                        "state": parts[0],
                        "energy_ev": float(parts[1]),
                        "wavelength_nm": float(parts[2]),
                        "fosc": float(parts[3]),
                        "dipole_strength": float(parts[4])
                    })
                except (ValueError, IndexError):
                    continue
        
        return pd.DataFrame(rows)
