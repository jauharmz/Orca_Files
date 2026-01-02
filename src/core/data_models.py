"""Data models for parsed ORCA data."""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import pandas as pd


@dataclass
class GeometryData:
    """Geometry data from ORCA output."""
    cart_coords: Optional[pd.DataFrame] = None  # atom, x, y, z
    internal_coords: Optional[pd.DataFrame] = None  # atom, bond, angle, dihedral
    smiles: Optional[str] = None
    filename: Optional[str] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None


@dataclass
class EnergyData:
    """Energy data from ORCA output."""
    gibbs_Eh: Optional[float] = None
    single_point_Eh: Optional[float] = None


@dataclass
class OrbitalData:
    """Orbital energy data."""
    orbitals: Optional[pd.DataFrame] = None  # OCC, Eh, eV, spin, lvl
    homo_energy: Optional[float] = None
    lumo_energy: Optional[float] = None
    homo_lumo_gap: Optional[float] = None


@dataclass
class SpectraData:
    """Spectroscopy data."""
    ir: Optional[pd.DataFrame] = None
    raman: Optional[pd.DataFrame] = None
    nmr_shielding: Optional[pd.DataFrame] = None
    nmr_coupling: Optional[pd.DataFrame] = None
    vibrations: Optional[pd.DataFrame] = None


@dataclass
class TDDFTData:
    """TD-DFT excited state data."""
    states: Optional[pd.DataFrame] = None
    electric_dipole: Optional[Dict[str, pd.DataFrame]] = None
    velocity_dipole: Optional[Dict[str, pd.DataFrame]] = None


@dataclass
class MullikenData:
    """Mulliken population analysis."""
    charges: Optional[pd.DataFrame] = None


@dataclass
class InternalCoordsData:
    """Internal coordinates (bonds, angles, dihedrals)."""
    bonds: Optional[pd.DataFrame] = None
    angles: Optional[pd.DataFrame] = None
    dihedrals: Optional[pd.DataFrame] = None


@dataclass
class ParseResult:
    """Complete parse result containing all data types."""
    filename: str
    geometry: GeometryData = field(default_factory=GeometryData)
    energy: EnergyData = field(default_factory=EnergyData)
    orbitals: OrbitalData = field(default_factory=OrbitalData)
    spectra: SpectraData = field(default_factory=SpectraData)
    tddft: TDDFTData = field(default_factory=TDDFTData)
    mulliken: MullikenData = field(default_factory=MullikenData)
    internal_coords: InternalCoordsData = field(default_factory=InternalCoordsData)
    spectrum_files: Dict[str, pd.DataFrame] = field(default_factory=dict)
    
    # Metadata
    is_optimization: bool = False
    has_tddft: bool = False
    optimized_state: Optional[str] = None
    esd_flag: Optional[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame creation."""
        return {
            "filename": self.filename,
            "cart_coords": self.geometry.cart_coords,
            "smiles": self.geometry.smiles,
            "charge": self.geometry.charge,
            "multiplicity": self.geometry.multiplicity,
            "gibbs_Eh": self.energy.gibbs_Eh,
            "single_point_Eh": self.energy.single_point_Eh,
            "orbitals": self.orbitals.orbitals,
            "homo_energy": self.orbitals.homo_energy,
            "lumo_energy": self.orbitals.lumo_energy,
            "ir": self.spectra.ir,
            "raman": self.spectra.raman,
            "tddft_states": self.tddft.states,
            "is_optimization": self.is_optimization,
            "has_tddft": self.has_tddft,
        }


@dataclass
class Pathway:
    """Represents a degradation pathway."""
    nodes: List[str]
    edges: List[Tuple[str, str]]
    reactions: Dict[Tuple[str, str], Dict] = field(default_factory=dict)
    color: Optional[str] = None


@dataclass
class Hierarchy:
    """Represents molecule hierarchy."""
    roots: List[str]
    variants: Dict[str, List[str]]  # root -> [variants]
    
    def to_tree(self) -> str:
        """Generate tree string representation."""
        lines = []
        for root in self.roots:
            lines.append(f"├── {root}")
            if root in self.variants:
                for variant in self.variants[root]:
                    lines.append(f"│   └── {variant}")
        return "\n".join(lines)


@dataclass
class Partitions:
    """Represents data partitions."""
    by_state: Dict[str, List[str]] = field(default_factory=dict)  # S0, S1, T1
    by_calc_type: Dict[str, List[str]] = field(default_factory=dict)  # OPT, SP
    by_esd_type: Dict[str, List[str]] = field(default_factory=dict)  # VG, AH
