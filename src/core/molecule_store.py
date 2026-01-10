"""
MoleculeStore - Hierarchical storage with canonical projection

Implements the ORCA data architecture from ARCHITECTURE.md:
- Storage: Molecule → Method → State → Task → Properties
- Access: Canonical projection for simple use cases
- Advanced: Explicit method/state selection for comparison

See ARCHITECTURE.md for the conceptual model.
"""

import re
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
import pandas as pd

from .data_models import ParseResult, MethodData
from ..logger import get_logger


@dataclass
class CalculationEntry:
    """Single calculation result in the hierarchy."""
    result: ParseResult
    task: str  # OPT, SP, TDDFT, ESD
    

class MoleculeStore:
    """
    Hierarchical storage for ORCA calculation results.
    
    Structure:
        molecule_id → method_id → state → task → ParseResult
    
    Usage:
        store = MoleculeStore()
        store.add(result)  # Auto-extracts molecule_id, method_id, state, task
        
        # Simple access (canonical projection)
        df = store.to_dataframe()
        mol = store.get("p1x")
        
        # Advanced access
        methods = store.get_methods("p1x")
        result = store.get("p1x", method_id="DFT/B3LYP/def2-TZVP/D3BJ", state="S0")
    """
    
    # Basis set quality ranking (for canonical selection)
    BASIS_RANK = {
        "DEF2-QZVPP": 7, "DEF2-QZVP": 6,
        "DEF2-TZVPP": 5, "DEF2-TZVP": 4,
        "DEF2-SVP": 2, "DEF2-SV(P)": 1,
        "6-311++G**": 4, "6-311+G**": 3, "6-31G**": 2, "6-31G*": 1,
        "CC-PVQZ": 6, "CC-PVTZ": 4, "CC-PVDZ": 2,
    }
    
    # State priority (for canonical selection)
    STATE_PRIORITY = {"S0": 0, "S1": 1, "T1": 2}
    
    # Task priority (for canonical selection)
    TASK_PRIORITY = {"OPT": 0, "SP": 1, "TDDFT": 2, "ESD": 3}
    
    def __init__(self):
        self.logger = get_logger("MoleculeStore")
        # Deep storage: molecule_id → method_id → state → task → CalculationEntry
        self._data: Dict[str, Dict[str, Dict[str, Dict[str, CalculationEntry]]]] = {}
        self._files: List[str] = []
    
    def add(self, result: ParseResult, 
            molecule_id: Optional[str] = None,
            task: Optional[str] = None) -> None:
        """
        Add a ParseResult to the store.
        
        Auto-extracts:
        - molecule_id from geometry.filename (or filename)
        - method_id from result.method
        - state from result.optimized_state
        - task from result.calc_class
        """
        # Extract identifiers
        if molecule_id is None:
            molecule_id = self._extract_molecule_id(result)
        
        method_id = result.method.to_id()
        state = result.optimized_state or "S0"
        
        if task is None:
            task = self._infer_task(result)
        
        # Store in hierarchy
        if molecule_id not in self._data:
            self._data[molecule_id] = {}
        
        if method_id not in self._data[molecule_id]:
            self._data[molecule_id][method_id] = {}
        
        if state not in self._data[molecule_id][method_id]:
            self._data[molecule_id][method_id][state] = {}
        
        # Store entry (overwrites if same path exists)
        entry = CalculationEntry(result=result, task=task)
        self._data[molecule_id][method_id][state][task] = entry
        
        self._files.append(result.filename)
        
        self.logger.debug(f"Added: {molecule_id}/{method_id}/{state}/{task}")
    
    def _extract_molecule_id(self, result: ParseResult) -> str:
        """Extract molecule_id from result."""
        # Try geometry.filename first
        if result.geometry.filename:
            return self._strip_suffixes(result.geometry.filename)
        
        # Fallback to filename
        if result.filename:
            import os
            base = os.path.splitext(os.path.basename(result.filename))[0]
            return self._strip_suffixes(base)
        
        return "unknown"
    
    def _strip_suffixes(self, name: str) -> str:
        """Remove state/task suffixes from name."""
        pattern = r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt|_?sp)$'
        mol_id = re.sub(pattern, '', name, flags=re.I)
        return mol_id.strip("_-") or name
    
    def _infer_task(self, result: ParseResult) -> str:
        """Infer task from result."""
        if result.is_optimization:
            return "OPT"
        if result.has_tddft:
            return "TDDFT"
        if result.esd_type:
            return "ESD"
        return "SP"
    
    # =========================================================================
    # SIMPLE ACCESS (Canonical Projection)
    # =========================================================================
    
    def get(self, molecule_id: str, 
            method_id: Optional[str] = None,
            state: Optional[str] = None,
            task: Optional[str] = None) -> Optional[ParseResult]:
        """
        Get a calculation result.
        
        With no arguments: returns canonical (best) result.
        With arguments: returns specific result.
        """
        if molecule_id not in self._data:
            return None
        
        mol_data = self._data[molecule_id]
        
        # Select method
        if method_id is None:
            method_id = self._select_canonical_method(mol_data)
        
        if method_id not in mol_data:
            return None
        
        method_data = mol_data[method_id]
        
        # Select state
        if state is None:
            state = self._select_canonical_state(method_data)
        
        if state not in method_data:
            return None
        
        state_data = method_data[state]
        
        # Select task
        if task is None:
            task = self._select_canonical_task(state_data)
        
        if task not in state_data:
            return None
        
        return state_data[task].result
    
    def get_canonical(self, molecule_id: str) -> Dict[str, Any]:
        """
        Get canonical projected view of molecule.
        
        Returns combined data from best available sources:
        - geometry from OPT
        - orbitals from SP
        - spectra from SP
        - tddft from TDDFT
        """
        if molecule_id not in self._data:
            return {}
        
        mol_data = self._data[molecule_id]
        method_id = self._select_canonical_method(mol_data)
        
        if method_id not in mol_data:
            return {}
        
        method_data = mol_data[method_id]
        
        # Build canonical view
        canonical = {
            "molecule_id": molecule_id,
            "method_id": method_id,
        }
        
        # Get geometry (prefer OPT > SP)
        for state in ["S0", "S1", "T1"]:
            if state not in method_data:
                continue
            state_data = method_data[state]
            
            if "OPT" in state_data:
                opt = state_data["OPT"].result
                canonical["geometry"] = opt.geometry.cart_coords
                canonical["smiles"] = opt.geometry.smiles
                canonical["charge"] = opt.geometry.charge
                canonical["multiplicity"] = opt.geometry.multiplicity
                canonical["vibrations"] = opt.spectra.vibrations
                break
            elif "SP" in state_data and "geometry" not in canonical:
                sp = state_data["SP"].result
                canonical["geometry"] = sp.geometry.cart_coords
        
        # Get spectra (IR, Raman from SP or OPT)
        for state in ["S0", "S1", "T1"]:
            if state not in method_data:
                continue
            for task in ["SP", "OPT"]:
                if task not in method_data[state]:
                    continue
                result = method_data[state][task].result
                if result.spectra.ir is not None and "ir" not in canonical:
                    canonical["ir"] = result.spectra.ir
                if result.spectra.raman is not None and "raman" not in canonical:
                    canonical["raman"] = result.spectra.raman
        
        # Get orbitals
        for state in ["S0", "S1", "T1"]:
            if state not in method_data:
                continue
            for task in ["SP", "OPT"]:
                if task not in method_data[state]:
                    continue
                result = method_data[state][task].result
                if result.orbitals.orbitals is not None:
                    canonical["orbitals"] = result.orbitals.orbitals
                    canonical["homo_energy"] = result.orbitals.homo_energy
                    canonical["lumo_energy"] = result.orbitals.lumo_energy
                    break
            if "orbitals" in canonical:
                break
        
        # Get TDDFT
        for state in ["S0", "S1", "T1"]:
            if state not in method_data:
                continue
            if "TDDFT" in method_data[state]:
                tddft = method_data[state]["TDDFT"].result
                canonical["tddft_states"] = tddft.tddft.states
                canonical["electric_dipole"] = tddft.tddft.electric_dipole_abs
                break
        
        # Get energies
        for state in ["S0", "S1", "T1"]:
            if state not in method_data:
                continue
            for task in ["OPT", "SP"]:
                if task not in method_data[state]:
                    continue
                result = method_data[state][task].result
                if result.energy.gibbs_Eh and "gibbs_Eh" not in canonical:
                    canonical["gibbs_Eh"] = result.energy.gibbs_Eh
                if result.energy.single_point_Eh and "single_point_Eh" not in canonical:
                    canonical["single_point_Eh"] = result.energy.single_point_Eh
        
        return canonical
    
    def _select_canonical_method(self, mol_data: Dict) -> str:
        """Select best method by basis set quality."""
        if not mol_data:
            return "unknown"
        
        methods = list(mol_data.keys())
        
        # Rank by basis set
        def rank_method(m: str) -> int:
            for basis, rank in self.BASIS_RANK.items():
                if basis in m.upper():
                    return rank
            return 0
        
        methods.sort(key=rank_method, reverse=True)
        return methods[0]
    
    def _select_canonical_state(self, method_data: Dict) -> str:
        """Select best state (S0 > S1 > T1)."""
        states = list(method_data.keys())
        states.sort(key=lambda s: self.STATE_PRIORITY.get(s, 99))
        return states[0] if states else "S0"
    
    def _select_canonical_task(self, state_data: Dict) -> str:
        """Select best task (OPT > SP > TDDFT)."""
        tasks = list(state_data.keys())
        tasks.sort(key=lambda t: self.TASK_PRIORITY.get(t, 99))
        return tasks[0] if tasks else "SP"
    
    # =========================================================================
    # DISCOVERY METHODS
    # =========================================================================
    
    def molecules(self) -> List[str]:
        """List all molecule IDs."""
        return list(self._data.keys())
    
    def get_methods(self, molecule_id: str) -> List[str]:
        """List all methods for a molecule."""
        if molecule_id not in self._data:
            return []
        return list(self._data[molecule_id].keys())
    
    def get_states(self, molecule_id: str, method_id: str) -> List[str]:
        """List all states for a molecule/method."""
        if molecule_id not in self._data:
            return []
        if method_id not in self._data[molecule_id]:
            return []
        return list(self._data[molecule_id][method_id].keys())
    
    def get_tasks(self, molecule_id: str, method_id: str, state: str) -> List[str]:
        """List all tasks for a molecule/method/state."""
        if molecule_id not in self._data:
            return []
        if method_id not in self._data[molecule_id]:
            return []
        if state not in self._data[molecule_id][method_id]:
            return []
        return list(self._data[molecule_id][method_id][state].keys())
    
    # =========================================================================
    # COMPARISON METHODS
    # =========================================================================
    
    def compare_methods(self, molecule_id: str, 
                       property_name: str = "single_point_Eh") -> pd.DataFrame:
        """Compare a property across methods for a molecule."""
        if molecule_id not in self._data:
            return pd.DataFrame()
        
        rows = []
        for method_id, method_data in self._data[molecule_id].items():
            for state, state_data in method_data.items():
                for task, entry in state_data.items():
                    result = entry.result
                    value = None
                    
                    if property_name == "single_point_Eh":
                        value = result.energy.single_point_Eh
                    elif property_name == "gibbs_Eh":
                        value = result.energy.gibbs_Eh
                    elif property_name == "homo_energy":
                        value = result.orbitals.homo_energy
                    elif property_name == "lumo_energy":
                        value = result.orbitals.lumo_energy
                    
                    rows.append({
                        "method_id": method_id,
                        "state": state,
                        "task": task,
                        property_name: value
                    })
        
        return pd.DataFrame(rows)
    
    # =========================================================================
    # EXPORT METHODS
    # =========================================================================
    
    def to_dataframe(self, canonical: bool = True) -> pd.DataFrame:
        """
        Export to DataFrame.
        
        Args:
            canonical: If True, returns one row per molecule (projected view).
                      If False, returns one row per calculation (full data).
        """
        if canonical:
            return self._to_dataframe_canonical()
        else:
            return self._to_dataframe_full()
    
    def _to_dataframe_canonical(self) -> pd.DataFrame:
        """Export canonical projected view."""
        rows = []
        for mol_id in self._data:
            row = self.get_canonical(mol_id)
            rows.append(row)
        return pd.DataFrame(rows)
    
    def _to_dataframe_full(self) -> pd.DataFrame:
        """Export all calculations."""
        rows = []
        for mol_id, mol_data in self._data.items():
            for method_id, method_data in mol_data.items():
                for state, state_data in method_data.items():
                    for task, entry in state_data.items():
                        row = entry.result.to_dict()
                        row["molecule_id"] = mol_id
                        row["_method_id"] = method_id
                        row["_state"] = state
                        row["_task"] = task
                        rows.append(row)
        return pd.DataFrame(rows)
    
    def summary(self) -> str:
        """Return summary string."""
        lines = [
            f"MoleculeStore: {len(self._data)} molecules, {len(self._files)} files",
            ""
        ]
        
        for mol_id in sorted(self._data.keys())[:10]:  # First 10
            mol_data = self._data[mol_id]
            n_methods = len(mol_data)
            n_calcs = sum(
                len(state_data)
                for method_data in mol_data.values()
                for state_data in method_data.values()
            )
            lines.append(f"  {mol_id}: {n_methods} method(s), {n_calcs} calc(s)")
        
        if len(self._data) > 10:
            lines.append(f"  ... and {len(self._data) - 10} more")
        
        return "\n".join(lines)
    
    def __repr__(self) -> str:
        return f"MoleculeStore({len(self._data)} molecules)"
    
    def __len__(self) -> int:
        return len(self._data)
