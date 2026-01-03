"""
Data Models - UPDATED to match original orca_praser.py output
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import pandas as pd


@dataclass
class GeometryData:
    """Geometry data from ORCA output."""
    filename: Optional[str] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    cart_coords: Optional[pd.DataFrame] = None  # columns: atom, x, y, z
    smiles: Optional[str] = None


@dataclass
class EnergyData:
    """Energy data from ORCA output."""
    gibbs_Eh: Optional[float] = None
    single_point_Eh: Optional[float] = None


@dataclass
class OrbitalData:
    """Orbital data from ORCA output."""
    orbitals: Optional[pd.DataFrame] = None  # columns: OCC, Eh, eV, spin, lvl
    homo_energy: Optional[float] = None  # eV
    lumo_energy: Optional[float] = None  # eV
    homo_lumo_gap: Optional[float] = None  # eV


@dataclass
class SpectraData:
    """Spectroscopy data from ORCA output."""
    vibrations: Optional[pd.DataFrame] = None  # columns: freq_cm-1, img
    ir: Optional[pd.DataFrame] = None  # columns: freq_cm-1, eps, intensity_km/mol
    raman: Optional[pd.DataFrame] = None  # columns: freq_cm-1, activity, depolarization
    nmr_shielding: Optional[pd.DataFrame] = None  # columns: Nucleus, Element, Isotropic, Anisotropy
    nmr_coupling: Optional[pd.DataFrame] = None  # columns: Nucleus1, Nucleus2, J_Hz


@dataclass
class TDDFTData:
    """TD-DFT data from ORCA output."""
    states: Optional[pd.DataFrame] = None  # columns: state, energy_au, energy_ev, from_orb, to_orb, weight, mult
    electric_dipole_abs: Optional[pd.DataFrame] = None
    electric_dipole_soc: Optional[pd.DataFrame] = None
    velocity_dipole_abs: Optional[pd.DataFrame] = None
    velocity_dipole_soc: Optional[pd.DataFrame] = None


@dataclass
class MullikenData:
    """Mulliken population analysis."""
    charges: Optional[pd.DataFrame] = None  # columns: Nucleus, Element, Population, Charge


@dataclass
class InternalCoordsData:
    """Internal coordinates."""
    bonds: Optional[pd.DataFrame] = None  # columns: Index, Type, Definition, Atoms, FinalVal, Unit
    angles: Optional[pd.DataFrame] = None
    dihedrals: Optional[pd.DataFrame] = None


@dataclass
class ParseResult:
    """Complete parse result from ORCA output."""
    filename: str = ""
    geometry: GeometryData = field(default_factory=GeometryData)
    energy: EnergyData = field(default_factory=EnergyData)
    orbitals: OrbitalData = field(default_factory=OrbitalData)
    spectra: SpectraData = field(default_factory=SpectraData)
    tddft: TDDFTData = field(default_factory=TDDFTData)
    mulliken: MullikenData = field(default_factory=MullikenData)
    internal_coords: InternalCoordsData = field(default_factory=InternalCoordsData)
    is_optimization: bool = False
    has_tddft: bool = False
    optimized_state: str = "S0"
    esd_type: Optional[str] = None
    calc_class: str = "single_point"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame."""
        return {
            "filename": self.filename,
            "molecule_id": self.geometry.filename,
            "smiles": self.geometry.smiles,
            "charge": self.geometry.charge,
            "multiplicity": self.geometry.multiplicity,
            "gibbs_Eh": self.energy.gibbs_Eh,
            "single_point_Eh": self.energy.single_point_Eh,
            "homo_energy": self.orbitals.homo_energy,
            "lumo_energy": self.orbitals.lumo_energy,
            "homo_lumo_gap": self.orbitals.homo_lumo_gap,
            "is_optimization": self.is_optimization,
            "has_tddft": self.has_tddft,
            "optimized_state": self.optimized_state,
            "esd_type": self.esd_type,
            "calc_class": self.calc_class,
            # Nested data
            "cart_coords": self.geometry.cart_coords,
            "orbitals": self.orbitals.orbitals,
            "vibrations": self.spectra.vibrations,
            "ir": self.spectra.ir,
            "raman": self.spectra.raman,
            "nmr_shielding": self.spectra.nmr_shielding,
            "nmr_coupling": self.spectra.nmr_coupling,
            "tddft_states": self.tddft.states,
            "electric_dipole_abs": self.tddft.electric_dipole_abs,
            "electric_dipole_soc": self.tddft.electric_dipole_soc,
            "velocity_dipole_abs": self.tddft.velocity_dipole_abs,
            "velocity_dipole_soc": self.tddft.velocity_dipole_soc,
            "mulliken": self.mulliken.charges,
            "bonds": self.internal_coords.bonds,
            "angles": self.internal_coords.angles,
            "dihedrals": self.internal_coords.dihedrals,
        }


# Analysis models
@dataclass
class Hierarchy:
    """Molecule hierarchy detection result."""
    roots: List[str] = field(default_factory=list)
    variants: Dict[str, List[str]] = field(default_factory=dict)
    
    def to_tree(self) -> str:
        """Return tree representation."""
        lines = []
        for root in sorted(self.roots):
            lines.append(f"├── {root}")
            if root in self.variants:
                for var in sorted(self.variants[root]):
                    lines.append(f"│   └── {var}")
        return "\n".join(lines) if lines else "(empty)"


@dataclass
class Partitions:
    """Data partitions by state/type."""
    by_state: Dict[str, List[str]] = field(default_factory=dict)  # S0, S1, T1 -> molecule_ids
    by_calc_type: Dict[str, List[str]] = field(default_factory=dict)  # OPT, SP -> molecule_ids
    by_esd_type: Dict[str, List[str]] = field(default_factory=dict)  # VG, AH, AHAS -> molecule_ids


@dataclass
class Pathway:
    """Degradation pathway."""
    nodes: List[str] = field(default_factory=list)
    edges: List[tuple] = field(default_factory=list)
    color: str = "blue"
    reactions: Dict[tuple, Dict] = field(default_factory=dict)
