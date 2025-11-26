"""
ORCA Main Output File Parser (.out)

Parses the main ORCA output file to extract:
- Job information (method, basis set, charge, multiplicity)
- Final energy and SCF energies
- Orbital energies (occupied and virtual)
- Dipole moment and polarizability
- Mulliken and Loewdin charges
- Mayer bond orders
- Vibrational frequencies and IR spectrum
- Thermochemistry (ZPE, enthalpy, entropy, Gibbs)
- NMR chemical shifts and J-couplings
"""

import re
import logging
from dataclasses import dataclass, field
from typing import Optional

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)


@dataclass
class JobInfo:
    """ORCA job information."""
    method: str = ""
    basis_set: str = ""
    charge: int = 0
    multiplicity: int = 1
    num_electrons: int = 0
    job_types: list[str] = field(default_factory=list)


@dataclass
class DipoleMoment:
    """Dipole moment data."""
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    magnitude_au: float = 0.0
    magnitude_debye: float = 0.0


@dataclass
class Polarizability:
    """Static polarizability tensor."""
    tensor: list[list[float]] = field(default_factory=list)  # 3x3 matrix
    isotropic: float = 0.0
    eigenvalues: list[float] = field(default_factory=list)


@dataclass
class OrbitalEnergy:
    """Single orbital energy."""
    index: int
    occupation: float
    energy_eh: float
    energy_ev: float


@dataclass
class ChemicalShieldingTensor:
    """Full chemical shielding tensor for a nucleus."""
    atom_index: int
    element: str
    # 3x3 tensors (stored as list of 9 values: [xx, xy, xz, yx, yy, yz, zx, zy, zz])
    diamagnetic_tensor: list[float] = field(default_factory=list)  # 3x3 matrix
    paramagnetic_tensor: list[float] = field(default_factory=list)  # 3x3 matrix
    total_tensor: list[float] = field(default_factory=list)  # 3x3 matrix
    # Diagonalized components
    sdso_components: list[float] = field(default_factory=list)  # [comp1, comp2, comp3]
    spso_components: list[float] = field(default_factory=list)  # [comp1, comp2, comp3]
    total_components: list[float] = field(default_factory=list)  # [comp1, comp2, comp3]
    # Isotropic values
    sdso_iso: float = 0.0
    spso_iso: float = 0.0
    total_iso: float = 0.0
    # Orientation eigenvectors (3 vectors of 3 components each)
    orientation_x: list[float] = field(default_factory=list)  # [x, y, z]
    orientation_y: list[float] = field(default_factory=list)  # [x, y, z]
    orientation_z: list[float] = field(default_factory=list)  # [x, y, z]


@dataclass
class NMRShift:
    """NMR chemical shift for a nucleus (summary data)."""
    atom_index: int
    element: str
    isotropic: float  # ppm
    anisotropy: float  # ppm


@dataclass
class JCoupling:
    """NMR J-coupling between two nuclei."""
    atom1_index: int
    atom2_index: int
    value: float  # Hz


@dataclass
class NMRData:
    """NMR spectroscopy data."""
    chemical_shifts: list[NMRShift] = field(default_factory=list)
    j_couplings: list[JCoupling] = field(default_factory=list)


@dataclass
class IRMode:
    """IR spectrum mode."""
    index: int
    frequency: float  # cm^-1
    intensity: float  # km/mol
    tx: float = 0.0
    ty: float = 0.0
    tz: float = 0.0


@dataclass
class RamanMode:
    """Raman spectrum mode."""
    index: int
    frequency: float  # cm^-1
    activity: float
    depolarization: float


@dataclass
class DispersionCorrection:
    """DFT dispersion correction data."""
    method: str = ""  # e.g., "DFTD3"
    total_correction: float = 0.0  # Eh
    total_correction_kcal: float = 0.0  # kcal/mol
    e6_kcal: float = 0.0
    e8_kcal: float = 0.0


@dataclass
class Thermochemistry:
    """Thermochemistry data."""
    temperature: float = 298.15
    pressure: float = 1.0
    total_mass: float = 0.0
    zero_point_energy: Optional[float] = None
    inner_energy: Optional[float] = None
    enthalpy: Optional[float] = None
    entropy: Optional[float] = None
    entropy_electronic: Optional[float] = None
    entropy_vibrational: Optional[float] = None
    entropy_rotational: Optional[float] = None
    entropy_translational: Optional[float] = None
    gibbs_free_energy: Optional[float] = None


@dataclass
class NormalMode:
    """Vibrational normal mode with displacement vectors."""
    mode_index: int
    frequency: float  # cm^-1
    displacements: list[tuple[float, float, float]] = field(default_factory=list)  # (dx, dy, dz) per atom


@dataclass
class SCFIteration:
    """Single SCF iteration data."""
    iteration: int
    energy: float  # Eh
    delta_e: float  # Eh
    rmsdp: float  # RMS density change
    maxdp: float  # Max density change
    diis_error: Optional[float] = None


@dataclass
class TimingData:
    """Computational timing breakdown."""
    total_time: float = 0.0  # seconds
    total_run_time: float = 0.0  # seconds (from TOTAL RUN TIME at end)
    scf_time: float = 0.0
    fock_matrix_time: float = 0.0
    diagonalization_time: float = 0.0
    scf_iterations_time: float = 0.0
    property_time: float = 0.0
    timings: dict[str, float] = field(default_factory=dict)  # All timing entries


@dataclass
class DFTGridInfo:
    """DFT numerical integration grid parameters."""
    integration_accuracy: float = 0.0
    radial_grid_type: str = ""
    angular_grid: str = ""
    pruning_method: str = ""
    total_grid_points: int = 0
    avg_points_per_atom: int = 0


@dataclass
class BasisSetInfo:
    """Basis set information."""
    name: str = ""
    num_basis_functions: int = 0
    num_primitive_gaussians: int = 0
    element_basis: dict[str, str] = field(default_factory=dict)  # element -> basis description


@dataclass
class EnergyComponents:
    """Detailed energy breakdown."""
    nuclear_repulsion: float = 0.0  # Eh
    electronic_energy: float = 0.0  # Eh
    one_electron_energy: float = 0.0  # Eh
    two_electron_energy: float = 0.0  # Eh
    kinetic_energy: float = 0.0  # Eh
    potential_energy: float = 0.0  # Eh
    virial_ratio: float = 0.0
    exchange_energy: float = 0.0  # Eh (DFT)
    correlation_energy: float = 0.0  # Eh (DFT)
    xc_energy: float = 0.0  # Eh (DFT)


@dataclass
class CPCMSolvation:
    """CPCM solvation model properties."""
    surface_charge: float = 0.0
    corrected_charge: float = 0.0
    outlying_charge_correction: float = 0.0  # Eh
    outlying_charge_correction_ev: float = 0.0  # eV
    dielectric_energy: float = 0.0  # Eh (CPCM Dielectric term)


@dataclass
class GeometricPerturbations:
    """Geometric perturbations data for numerical derivatives."""
    num_nuclei: int = 0
    num_perturbations: int = 0
    num_batches: int = 0
    total_time: float = 0.0  # seconds
    max_memory_mb: float = 0.0


@dataclass
class SharkIntegrals:
    """SHARK integral package configuration."""
    num_basis_functions: int = 0
    num_shells: int = 0
    max_angular_momentum: int = 0
    num_aux_j_functions: int = 0
    num_aux_j_shells: int = 0
    max_angular_momentum_aux_j: int = 0
    total_shell_pairs: int = 0
    shell_pairs_after_screening: int = 0
    schwartz_threshold: float = 0.0
    tcut_threshold: float = 0.0


@dataclass
class COSXGrid:
    """COSX grid generation data for a single grid."""
    grid_number: int = 0
    integration_accuracy: float = 0.0
    radial_grid_type: str = ""
    angular_grid: str = ""
    total_grid_points: int = 0
    num_batches: int = 0
    avg_points_per_batch: int = 0
    avg_points_per_atom: int = 0


@dataclass
class SCFConvergence:
    """SCF convergence criteria values."""
    energy_change: float = 0.0
    max_density_change: float = 0.0
    rms_density_change: float = 0.0
    diis_error: float = 0.0
    orbital_gradient: float = 0.0
    orbital_rotation: float = 0.0


@dataclass
class OrbitalPopulation:
    """Orbital population breakdown for an atom."""
    atom_index: int
    element: str
    s_total: float = 0.0
    p_total: float = 0.0
    d_total: float = 0.0
    f_total: float = 0.0
    g_total: float = 0.0
    orbital_details: dict[str, float] = field(default_factory=dict)  # e.g., {'px': 0.885, 'py': 1.012, ...}


@dataclass
class InternalCoordinate:
    """Internal coordinate (Z-matrix) entry."""
    atom_index: int
    element: str
    bond_to: int  # Atom index for bond
    angle_to: int  # Atom index for angle
    dihedral_to: int  # Atom index for dihedral
    bond_length: float  # Angstrom or a.u.
    bond_angle: float  # degrees
    dihedral_angle: float  # degrees


@dataclass
class MOCharge:
    """Molecular orbital charge entry."""
    mo_index: int
    atom_index: int
    element: str
    orbital: str  # e.g., '1s', '2px', '3dz2'
    charge: float
    uncorrected_charge: Optional[float] = None  # Only for Mulliken


@dataclass
class TDDFTState:
    """TD-DFT excited state data."""
    state: int
    energy_au: float
    energy_ev: float
    energy_cm1: float
    homo: Optional[int] = None  # HOMO index relative (e.g., -1 for HOMO-1)
    lumo: Optional[int] = None  # LUMO index relative (e.g., 0 for LUMO, 1 for LUMO+1)
    weight: float = 0.0
    coeff: float = 0.0
    from_orb: int = 0
    to_orb: int = 0


@dataclass
class ElectricDipoleTransition:
    """Electric dipole absorption spectrum transition."""
    transition: str  # e.g., "0 -> 1"
    energy_ev: float
    energy_cm1: float
    wavelength_nm: float
    fosc_d2: float  # Oscillator strength
    d2: float  # Dipole strength
    dx: float
    dy: float
    dz: float


@dataclass
class VelocityDipoleTransition:
    """Velocity dipole absorption spectrum transition."""
    transition: str
    energy_ev: float
    energy_cm1: float
    wavelength_nm: float
    fosc_p2: float  # Oscillator strength
    p2: float  # Velocity dipole strength
    px: float
    py: float
    pz: float


@dataclass
class AbsorptionSpectra:
    """Absorption spectra data (regular and SOC-corrected)."""
    regular: list[ElectricDipoleTransition] = field(default_factory=list)
    soc_corrected: list[ElectricDipoleTransition] = field(default_factory=list)


@dataclass
class VelocityAbsorptionSpectra:
    """Velocity dipole absorption spectra data."""
    regular: list[VelocityDipoleTransition] = field(default_factory=list)
    soc_corrected: list[VelocityDipoleTransition] = field(default_factory=list)


@dataclass
class OrcaOutput:
    """Complete parsed ORCA output."""
    job_info: JobInfo
    final_energy: Optional[float] = None
    scf_energies: list[float] = field(default_factory=list)
    optimization_energies: list[float] = field(default_factory=list)
    coordinates: list[tuple[str, float, float, float]] = field(default_factory=list)  # (element, x, y, z) in Angstrom
    coordinates_au: list[tuple[str, float, float, float, float, float]] = field(default_factory=list)  # (element, x, y, z, atomic_number, mass) in a.u.
    internal_coords: list[InternalCoordinate] = field(default_factory=list)  # Z-matrix in Angstrom
    internal_coords_au: list[InternalCoordinate] = field(default_factory=list)  # Z-matrix in a.u.
    dipole_moment: Optional[DipoleMoment] = None
    polarizability: Optional[Polarizability] = None
    orbital_energies: list[OrbitalEnergy] = field(default_factory=list)
    frequencies: list[float] = field(default_factory=list)
    ir_spectrum: list[IRMode] = field(default_factory=list)
    raman_spectrum: list[RamanMode] = field(default_factory=list)
    normal_modes: list[NormalMode] = field(default_factory=list)
    scf_iterations: list[SCFIteration] = field(default_factory=list)
    dispersion_correction: Optional[DispersionCorrection] = None
    timing_data: Optional[TimingData] = None
    dft_grid_info: Optional[DFTGridInfo] = None
    basis_set_info: Optional[BasisSetInfo] = None
    energy_components: Optional[EnergyComponents] = None
    cpcm_solvation: Optional[CPCMSolvation] = None
    scf_convergence: Optional[SCFConvergence] = None
    mulliken_charges: dict[int, tuple[str, float]] = field(default_factory=dict)
    mulliken_orbital_populations: list[OrbitalPopulation] = field(default_factory=list)
    mulliken_orbital_charges: list[MOCharge] = field(default_factory=list)  # Per-MO charges
    loewdin_charges: dict[int, tuple[str, float]] = field(default_factory=dict)
    loewdin_orbital_charges: list[MOCharge] = field(default_factory=list)  # Per-MO charges
    mayer_bond_orders: list[tuple[int, int, float]] = field(default_factory=list)
    loewdin_bond_orders: list[tuple[int, int, float]] = field(default_factory=list)
    mulliken_overlap_charges: list[tuple[int, int, float]] = field(default_factory=list)  # (atom1, atom2, overlap_charge)
    thermochemistry: Optional[Thermochemistry] = None
    nmr_data: Optional[NMRData] = None
    chemical_shielding_tensors: list[ChemicalShieldingTensor] = field(default_factory=list)  # Full tensor data
    geometric_perturbations: Optional[GeometricPerturbations] = None
    shark_integrals: Optional[SharkIntegrals] = None
    cosx_grids: list[COSXGrid] = field(default_factory=list)
    # TD-DFT and UV-Vis spectroscopy
    tddft_states: list[TDDFTState] = field(default_factory=list)
    electric_dipole_spectrum: Optional[AbsorptionSpectra] = None
    velocity_dipole_spectrum: Optional[VelocityAbsorptionSpectra] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'job_info': {
                'method': self.job_info.method,
                'basis_set': self.job_info.basis_set,
                'charge': self.job_info.charge,
                'multiplicity': self.job_info.multiplicity,
                'num_electrons': self.job_info.num_electrons,
                'job_types': self.job_info.job_types
            },
            'final_energy': self.final_energy,
            'scf_energies': self.scf_energies,
            'optimization_energies': self.optimization_energies,
            'coordinates': [
                {'element': c[0], 'x': c[1], 'y': c[2], 'z': c[3]}
                for c in self.coordinates
            ],
            'coordinates_au': [
                {'element': c[0], 'x': c[1], 'y': c[2], 'z': c[3], 'atomic_number': c[4], 'mass': c[5]}
                for c in self.coordinates_au
            ],
            'internal_coords': [
                {
                    'atom_index': ic.atom_index,
                    'element': ic.element,
                    'bond_to': ic.bond_to,
                    'angle_to': ic.angle_to,
                    'dihedral_to': ic.dihedral_to,
                    'bond_length': ic.bond_length,
                    'bond_angle': ic.bond_angle,
                    'dihedral_angle': ic.dihedral_angle
                }
                for ic in self.internal_coords
            ],
            'internal_coords_au': [
                {
                    'atom_index': ic.atom_index,
                    'element': ic.element,
                    'bond_to': ic.bond_to,
                    'angle_to': ic.angle_to,
                    'dihedral_to': ic.dihedral_to,
                    'bond_length': ic.bond_length,
                    'bond_angle': ic.bond_angle,
                    'dihedral_angle': ic.dihedral_angle
                }
                for ic in self.internal_coords_au
            ],
            'dipole_moment': {
                'x': self.dipole_moment.x,
                'y': self.dipole_moment.y,
                'z': self.dipole_moment.z,
                'magnitude_au': self.dipole_moment.magnitude_au,
                'magnitude_debye': self.dipole_moment.magnitude_debye
            } if self.dipole_moment else None,
            'polarizability': {
                'tensor': self.polarizability.tensor,
                'isotropic': self.polarizability.isotropic,
                'eigenvalues': self.polarizability.eigenvalues
            } if self.polarizability else None,
            'orbital_energies': [
                {'index': o.index, 'occupation': o.occupation,
                 'energy_eh': o.energy_eh, 'energy_ev': o.energy_ev}
                for o in self.orbital_energies
            ],
            'frequencies': self.frequencies,
            'ir_spectrum': [
                {'index': m.index, 'frequency': m.frequency, 'intensity': m.intensity}
                for m in self.ir_spectrum
            ],
            'raman_spectrum': [
                {'index': m.index, 'frequency': m.frequency, 'activity': m.activity, 'depolarization': m.depolarization}
                for m in self.raman_spectrum
            ],
            'normal_modes': [
                {'mode_index': nm.mode_index, 'frequency': nm.frequency, 'num_atoms': len(nm.displacements)}
                for nm in self.normal_modes
            ],
            'scf_iterations': [
                {'iteration': it.iteration, 'energy': it.energy, 'delta_e': it.delta_e,
                 'rmsdp': it.rmsdp, 'maxdp': it.maxdp, 'diis_error': it.diis_error}
                for it in self.scf_iterations
            ],
            'dispersion_correction': {
                'method': self.dispersion_correction.method,
                'total_correction': self.dispersion_correction.total_correction,
                'total_correction_kcal': self.dispersion_correction.total_correction_kcal,
                'e6_kcal': self.dispersion_correction.e6_kcal,
                'e8_kcal': self.dispersion_correction.e8_kcal
            } if self.dispersion_correction else None,
            'timing_data': {
                'total_time': self.timing_data.total_time,
                'total_run_time': self.timing_data.total_run_time,
                'scf_time': self.timing_data.scf_time,
                'fock_matrix_time': self.timing_data.fock_matrix_time,
                'timings': self.timing_data.timings
            } if self.timing_data else None,
            'dft_grid_info': {
                'integration_accuracy': self.dft_grid_info.integration_accuracy,
                'radial_grid_type': self.dft_grid_info.radial_grid_type,
                'angular_grid': self.dft_grid_info.angular_grid,
                'total_grid_points': self.dft_grid_info.total_grid_points,
                'avg_points_per_atom': self.dft_grid_info.avg_points_per_atom
            } if self.dft_grid_info else None,
            'basis_set_info': {
                'name': self.basis_set_info.name,
                'num_basis_functions': self.basis_set_info.num_basis_functions,
                'num_primitive_gaussians': self.basis_set_info.num_primitive_gaussians,
                'element_basis': self.basis_set_info.element_basis
            } if self.basis_set_info else None,
            'energy_components': {
                'nuclear_repulsion': self.energy_components.nuclear_repulsion,
                'electronic_energy': self.energy_components.electronic_energy,
                'one_electron_energy': self.energy_components.one_electron_energy,
                'two_electron_energy': self.energy_components.two_electron_energy,
                'kinetic_energy': self.energy_components.kinetic_energy,
                'potential_energy': self.energy_components.potential_energy,
                'virial_ratio': self.energy_components.virial_ratio,
                'exchange_energy': self.energy_components.exchange_energy,
                'correlation_energy': self.energy_components.correlation_energy,
                'xc_energy': self.energy_components.xc_energy
            } if self.energy_components else None,
            'cpcm_solvation': {
                'surface_charge': self.cpcm_solvation.surface_charge,
                'corrected_charge': self.cpcm_solvation.corrected_charge,
                'outlying_charge_correction': self.cpcm_solvation.outlying_charge_correction,
                'outlying_charge_correction_ev': self.cpcm_solvation.outlying_charge_correction_ev,
                'dielectric_energy': self.cpcm_solvation.dielectric_energy
            } if self.cpcm_solvation else None,
            'scf_convergence': {
                'energy_change': self.scf_convergence.energy_change,
                'max_density_change': self.scf_convergence.max_density_change,
                'rms_density_change': self.scf_convergence.rms_density_change,
                'diis_error': self.scf_convergence.diis_error,
                'orbital_gradient': self.scf_convergence.orbital_gradient,
                'orbital_rotation': self.scf_convergence.orbital_rotation
            } if self.scf_convergence else None,
            'mulliken_charges': {
                str(k): {'element': v[0], 'charge': v[1]}
                for k, v in self.mulliken_charges.items()
            },
            'mulliken_orbital_populations': [
                {
                    'atom_index': pop.atom_index,
                    'element': pop.element,
                    's': pop.s_total,
                    'p': pop.p_total,
                    'd': pop.d_total,
                    'f': pop.f_total,
                    'g': pop.g_total,
                    'orbital_details': pop.orbital_details
                }
                for pop in self.mulliken_orbital_populations
            ],
            'mulliken_orbital_charges': [
                {
                    'mo_index': charge.mo_index,
                    'atom_index': charge.atom_index,
                    'element': charge.element,
                    'orbital': charge.orbital,
                    'charge': charge.charge,
                    'uncorrected_charge': charge.uncorrected_charge
                }
                for charge in self.mulliken_orbital_charges
            ],
            'loewdin_charges': {
                str(k): {'element': v[0], 'charge': v[1]}
                for k, v in self.loewdin_charges.items()
            },
            'loewdin_orbital_charges': [
                {
                    'mo_index': charge.mo_index,
                    'atom_index': charge.atom_index,
                    'element': charge.element,
                    'orbital': charge.orbital,
                    'charge': charge.charge
                }
                for charge in self.loewdin_orbital_charges
            ],
            'mayer_bond_orders': [
                {'atom1': b[0], 'atom2': b[1], 'order': b[2]}
                for b in self.mayer_bond_orders
            ],
            'loewdin_bond_orders': [
                {'atom1': b[0], 'atom2': b[1], 'order': b[2]}
                for b in self.loewdin_bond_orders
            ],
            'mulliken_overlap_charges': [
                {'atom1': oc[0], 'atom2': oc[1], 'overlap_charge': oc[2]}
                for oc in self.mulliken_overlap_charges
            ],
            'thermochemistry': {
                'temperature': self.thermochemistry.temperature,
                'pressure': self.thermochemistry.pressure,
                'total_mass': self.thermochemistry.total_mass,
                'zero_point_energy': self.thermochemistry.zero_point_energy,
                'inner_energy': self.thermochemistry.inner_energy,
                'enthalpy': self.thermochemistry.enthalpy,
                'entropy': self.thermochemistry.entropy,
                'entropy_electronic': self.thermochemistry.entropy_electronic,
                'entropy_vibrational': self.thermochemistry.entropy_vibrational,
                'entropy_rotational': self.thermochemistry.entropy_rotational,
                'entropy_translational': self.thermochemistry.entropy_translational,
                'gibbs_free_energy': self.thermochemistry.gibbs_free_energy
            } if self.thermochemistry else None,
            'nmr_data': {
                'chemical_shifts': [
                    {'atom_index': s.atom_index, 'element': s.element,
                     'isotropic': s.isotropic, 'anisotropy': s.anisotropy}
                    for s in self.nmr_data.chemical_shifts
                ],
                'j_couplings': [
                    {'atom1': j.atom1_index, 'atom2': j.atom2_index, 'value': j.value}
                    for j in self.nmr_data.j_couplings
                ]
            } if self.nmr_data else None,
            'chemical_shielding_tensors': [
                {
                    'atom_index': t.atom_index,
                    'element': t.element,
                    'diamagnetic_tensor': t.diamagnetic_tensor,
                    'paramagnetic_tensor': t.paramagnetic_tensor,
                    'total_tensor': t.total_tensor,
                    'sdso_components': t.sdso_components,
                    'spso_components': t.spso_components,
                    'total_components': t.total_components,
                    'sdso_iso': t.sdso_iso,
                    'spso_iso': t.spso_iso,
                    'total_iso': t.total_iso,
                    'orientation_x': t.orientation_x,
                    'orientation_y': t.orientation_y,
                    'orientation_z': t.orientation_z
                }
                for t in self.chemical_shielding_tensors
            ],
            'geometric_perturbations': {
                'num_nuclei': self.geometric_perturbations.num_nuclei,
                'num_perturbations': self.geometric_perturbations.num_perturbations,
                'num_batches': self.geometric_perturbations.num_batches,
                'total_time': self.geometric_perturbations.total_time,
                'max_memory_mb': self.geometric_perturbations.max_memory_mb
            } if self.geometric_perturbations else None,
            'shark_integrals': {
                'num_basis_functions': self.shark_integrals.num_basis_functions,
                'num_shells': self.shark_integrals.num_shells,
                'max_angular_momentum': self.shark_integrals.max_angular_momentum,
                'num_aux_j_functions': self.shark_integrals.num_aux_j_functions,
                'num_aux_j_shells': self.shark_integrals.num_aux_j_shells,
                'max_angular_momentum_aux_j': self.shark_integrals.max_angular_momentum_aux_j,
                'total_shell_pairs': self.shark_integrals.total_shell_pairs,
                'shell_pairs_after_screening': self.shark_integrals.shell_pairs_after_screening,
                'schwartz_threshold': self.shark_integrals.schwartz_threshold,
                'tcut_threshold': self.shark_integrals.tcut_threshold
            } if self.shark_integrals else None,
            'cosx_grids': [
                {
                    'grid_number': g.grid_number,
                    'integration_accuracy': g.integration_accuracy,
                    'radial_grid_type': g.radial_grid_type,
                    'angular_grid': g.angular_grid,
                    'total_grid_points': g.total_grid_points,
                    'num_batches': g.num_batches,
                    'avg_points_per_batch': g.avg_points_per_batch,
                    'avg_points_per_atom': g.avg_points_per_atom
                }
                for g in self.cosx_grids
            ]
        }


def parse_out_file(filepath: str) -> OrcaOutput:
    """Parse ORCA output file."""
    logger.info(f"Starting to parse ORCA output file: {filepath}")
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        logger.debug(f"Successfully read file {filepath}, size: {len(content)} bytes")
        result = parse_out_content(content)
        logger.info(f"Successfully parsed {filepath}")
        return result
    except FileNotFoundError:
        logger.error(f"File not found: {filepath}")
        raise
    except IOError as e:
        logger.error(f"Error reading file {filepath}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error parsing {filepath}: {e}", exc_info=True)
        raise


def parse_out_content(content: str) -> OrcaOutput:
    """Parse ORCA output content from string."""
    logger.debug("Starting content parsing...")
    logger.debug("Parsing job information...")
    job_info = parse_job_info(content)
    final_energy = parse_final_energy(content)
    scf_energies = parse_scf_energies(content)
    opt_energies = parse_optimization_energies(content)
    coordinates = parse_coordinates(content)
    coordinates_au = parse_coordinates_au(content)
    internal_coords = parse_internal_coordinates(content, "ANGSTROEM")
    internal_coords_au = parse_internal_coordinates(content, "A.U.")
    dipole = parse_dipole_moment(content)
    polarizability = parse_polarizability(content)
    orbitals = parse_orbital_energies(content)
    frequencies = parse_frequencies(content)
    ir_spectrum = parse_ir_spectrum(content)
    raman_spectrum = parse_raman_spectrum(content)
    dispersion = parse_dispersion_correction(content)
    mulliken = parse_mulliken_charges(content)
    loewdin = parse_loewdin_charges(content)
    bond_orders = parse_mayer_bond_orders(content)
    loewdin_bond_orders = parse_loewdin_bond_orders(content)
    mulliken_overlap = parse_mulliken_overlap_charges(content)
    thermo = parse_thermochemistry(content)
    nmr = parse_nmr_data(content)
    shielding_tensors = parse_chemical_shielding_tensors(content)
    geom_pert = parse_geometric_perturbations(content)
    shark = parse_shark_integrals(content)
    cosx = parse_cosx_grids(content)
    scf_iterations = parse_scf_iterations(content)
    timing = parse_timing_data(content)
    dft_grid = parse_dft_grid_info(content)
    basis_set = parse_basis_set_info(content)
    energy_comp = parse_energy_components(content)
    cpcm_solv = parse_cpcm_solvation(content)
    scf_conv = parse_scf_convergence(content)
    mulliken_orb_pop = parse_mulliken_orbital_populations(content)
    mulliken_orb_charges = parse_mulliken_orbital_charges(content)
    loewdin_orb_charges = parse_loewdin_orbital_charges(content)

    # Get number of atoms for normal modes parsing
    num_atoms = len(mulliken) if mulliken else 0
    normal_modes = parse_normal_modes(content, num_atoms) if num_atoms > 0 else []

    return OrcaOutput(
        job_info=job_info,
        final_energy=final_energy,
        scf_energies=scf_energies,
        optimization_energies=opt_energies,
        coordinates=coordinates,
        coordinates_au=coordinates_au,
        internal_coords=internal_coords,
        internal_coords_au=internal_coords_au,
        dipole_moment=dipole,
        polarizability=polarizability,
        orbital_energies=orbitals,
        frequencies=frequencies,
        ir_spectrum=ir_spectrum,
        raman_spectrum=raman_spectrum,
        normal_modes=normal_modes,
        scf_iterations=scf_iterations,
        dispersion_correction=dispersion,
        timing_data=timing,
        dft_grid_info=dft_grid,
        basis_set_info=basis_set,
        energy_components=energy_comp,
        cpcm_solvation=cpcm_solv,
        scf_convergence=scf_conv,
        mulliken_charges=mulliken,
        mulliken_orbital_populations=mulliken_orb_pop,
        mulliken_orbital_charges=mulliken_orb_charges,
        loewdin_charges=loewdin,
        loewdin_orbital_charges=loewdin_orb_charges,
        mayer_bond_orders=bond_orders,
        loewdin_bond_orders=loewdin_bond_orders,
        mulliken_overlap_charges=mulliken_overlap,
        thermochemistry=thermo,
        nmr_data=nmr,
        chemical_shielding_tensors=shielding_tensors,
        geometric_perturbations=geom_pert,
        shark_integrals=shark,
        cosx_grids=cosx
    )


def parse_coordinates(content: str) -> list[tuple[str, float, float, float]]:
    """Extract Cartesian coordinates in Angstrom."""
    coordinates = []

    # Find the first CARTESIAN COORDINATES (ANGSTROEM) section
    coord_section = re.search(
        r'CARTESIAN COORDINATES \(ANGSTROEM\)\s*-+\s*(.*?)(?:\n\n|$)',
        content, re.DOTALL
    )

    if not coord_section:
        return coordinates

    section_text = coord_section.group(1)

    # Parse each coordinate line: "  C      0.976427    4.012579   -0.068330"
    coord_lines = re.findall(
        r'^\s*([A-Z][a-z]?)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)',
        section_text, re.MULTILINE
    )

    for element, x, y, z in coord_lines:
        coordinates.append((element, float(x), float(y), float(z)))

    return coordinates


def parse_coordinates_au(content: str) -> list[tuple[str, float, float, float, float, float]]:
    """Extract Cartesian coordinates in atomic units with atomic number and mass."""
    coordinates = []

    # Find the first CARTESIAN COORDINATES (A.U.) section
    coord_section = re.search(
        r'CARTESIAN COORDINATES \(A\.U\.\)\s*-+\s*NO\s+LB\s+ZA\s+FRAG\s+MASS\s+X\s+Y\s+Z\s+(.*?)(?:\n\n|$)',
        content, re.DOTALL
    )

    if not coord_section:
        return coordinates

    section_text = coord_section.group(1)

    # Parse each coordinate line: "   0 C     6.0000    0    12.011    1.845179    7.582676   -0.129125"
    coord_lines = re.findall(
        r'^\s*\d+\s+([A-Z][a-z]?)\s+(\d+\.?\d*)\s+\d+\s+(\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)',
        section_text, re.MULTILINE
    )

    for element, za, mass, x, y, z in coord_lines:
        coordinates.append((element, float(x), float(y), float(z), float(za), float(mass)))

    return coordinates


def parse_internal_coordinates(content: str, units: str = "ANGSTROEM") -> list[InternalCoordinate]:
    """Extract internal coordinates (Z-matrix format)."""
    coords = []

    # Find the INTERNAL COORDINATES section
    pattern = rf'INTERNAL COORDINATES \({units}\)\s*-+\s*(.*?)(?:\n-+|$)'
    coord_section = re.search(pattern, content, re.DOTALL)

    if not coord_section:
        return coords

    section_text = coord_section.group(1)

    # Parse each line: " C      0   0   0     0.000000000000     0.00000000     0.00000000"
    # Format: Element, bond_to, angle_to, dihedral_to, bond_length, bond_angle, dihedral_angle
    coord_lines = re.findall(
        r'^\s*([A-Z][a-z]?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)',
        section_text, re.MULTILINE
    )

    for idx, (element, bond_to, angle_to, dihedral_to, bond_length, bond_angle, dihedral_angle) in enumerate(coord_lines):
        coords.append(InternalCoordinate(
            atom_index=idx,
            element=element,
            bond_to=int(bond_to),
            angle_to=int(angle_to),
            dihedral_to=int(dihedral_to),
            bond_length=float(bond_length),
            bond_angle=float(bond_angle),
            dihedral_angle=float(dihedral_angle)
        ))

    return coords


def parse_job_info(content: str) -> JobInfo:
    """Extract job information from output."""
    info = JobInfo()

    # Extract method and job types from input line
    input_match = re.search(r'\|\s*>\s*!\s*(.+)', content)
    if input_match:
        keywords = input_match.group(1).split()
        methods = ['B3LYP', 'HF', 'MP2', 'CCSD', 'PBE', 'PBE0', 'M06', 'wB97X']
        job_types = ['OPT', 'FREQ', 'SP', 'TD', 'NMR', 'TDDFT']

        for kw in keywords:
            kw_upper = kw.upper()
            if kw_upper in methods or any(m in kw_upper for m in methods):
                info.method = kw
            if kw_upper in job_types:
                info.job_types.append(kw_upper)

    # Extract basis set
    basis_match = re.search(r'Your calculation utilizes the basis:\s*(\S+)', content)
    if basis_match:
        info.basis_set = basis_match.group(1)

    # Extract charge
    charge_match = re.search(r'Total Charge\s+Charge\s+\.\.\.\.\s+(-?\d+)', content)
    if charge_match:
        info.charge = int(charge_match.group(1))

    # Extract multiplicity
    mult_match = re.search(r'Multiplicity\s+Mult\s+\.\.\.\.\s+(\d+)', content)
    if mult_match:
        info.multiplicity = int(mult_match.group(1))

    # Extract number of electrons
    nel_match = re.search(r'Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)', content)
    if nel_match:
        info.num_electrons = int(nel_match.group(1))

    return info


def parse_final_energy(content: str) -> Optional[float]:
    """Extract final single point energy."""
    matches = re.findall(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.?\d*)', content)
    if matches:
        return float(matches[-1])
    return None


def parse_scf_energies(content: str) -> list[float]:
    """Extract all SCF energies."""
    energies = []
    matches = re.findall(r'Total Energy\s+:\s+(-?\d+\.?\d*)\s+Eh', content)
    for match in matches:
        energies.append(float(match))
    return energies


def parse_optimization_energies(content: str) -> list[float]:
    """Extract energies from optimization cycles."""
    energies = []
    matches = re.findall(r'FINAL SINGLE POINT ENERGY\s+(-?\d+\.?\d*)', content)
    for match in matches:
        energies.append(float(match))
    return energies


def parse_dipole_moment(content: str) -> Optional[DipoleMoment]:
    """Extract dipole moment."""
    dipole_section = re.search(
        r'DIPOLE MOMENT\s*-+\s*.*?Total Dipole Moment\s+:\s+'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*).*?'
        r'Magnitude \(a\.u\.\)\s+:\s+(-?\d+\.?\d*).*?'
        r'Magnitude \(Debye\)\s+:\s+(-?\d+\.?\d*)',
        content, re.DOTALL
    )

    if dipole_section:
        return DipoleMoment(
            x=float(dipole_section.group(1)),
            y=float(dipole_section.group(2)),
            z=float(dipole_section.group(3)),
            magnitude_au=float(dipole_section.group(4)),
            magnitude_debye=float(dipole_section.group(5))
        )
    return None


def parse_polarizability(content: str) -> Optional[Polarizability]:
    """Extract static polarizability tensor with eigenvalues."""
    pol_section = re.search(
        r'(?:THE|STATIC) POLARIZABILITY TENSOR.*?'
        r'The raw cartesian tensor \(atomic units\):\s*'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*).*?'
        r'diagonalized tensor:\s*'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*).*?'
        r'Isotropic polarizability\s*:\s*(-?\d+\.?\d*)',
        content, re.DOTALL
    )

    if pol_section:
        tensor = [
            [float(pol_section.group(1)), float(pol_section.group(2)), float(pol_section.group(3))],
            [float(pol_section.group(4)), float(pol_section.group(5)), float(pol_section.group(6))],
            [float(pol_section.group(7)), float(pol_section.group(8)), float(pol_section.group(9))]
        ]
        eigenvalues = [
            float(pol_section.group(10)),
            float(pol_section.group(11)),
            float(pol_section.group(12))
        ]
        return Polarizability(
            tensor=tensor,
            isotropic=float(pol_section.group(13)),
            eigenvalues=eigenvalues
        )
    return None


def parse_orbital_energies(content: str) -> list[OrbitalEnergy]:
    """Extract orbital energies (occupied and virtual)."""
    orbitals = []

    # Find the last orbital energies section
    sections = re.findall(
        r'ORBITAL ENERGIES\s*-+\s*(.*?)(?=-{20,}|\Z)',
        content, re.DOTALL
    )

    if sections:
        last_section = sections[-1]
        # Match orbital lines: NO, OCC, E(Eh), E(eV)
        matches = re.findall(
            r'^\s*(\d+)\s+(\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)',
            last_section, re.MULTILINE
        )
        for match in matches:
            orbitals.append(OrbitalEnergy(
                index=int(match[0]),
                occupation=float(match[1]),
                energy_eh=float(match[2]),
                energy_ev=float(match[3])
            ))

    return orbitals


def parse_frequencies(content: str) -> list[float]:
    """Extract vibrational frequencies in cm^-1."""
    frequencies = []

    freq_section = re.search(
        r'VIBRATIONAL FREQUENCIES\s*-+\s*(.*?)(?=\s*-{5,}|NORMAL MODES)',
        content, re.DOTALL
    )

    if freq_section:
        freq_matches = re.findall(
            r'^\s*\d+:\s+(-?\d+\.?\d*)\s+cm\*\*-1',
            freq_section.group(1), re.MULTILINE
        )
        for match in freq_matches:
            freq = float(match)
            if freq > 0:
                frequencies.append(freq)

    return frequencies


def parse_ir_spectrum(content: str) -> list[IRMode]:
    """Extract IR spectrum data."""
    modes = []

    ir_section = re.search(
        r'IR SPECTRUM\s*-+\s*Mode.*?-+\s*(.*?)(?:\n\s*\n|\Z)',
        content, re.DOTALL
    )

    if ir_section:
        # Match: Mode, freq (cm^-1), eps, Int (km/mol), T^2, TX, TY, TZ
        matches = re.findall(
            r'^\s*(\d+):\s+([0-9.]+)\s+[0-9.]+\s+([0-9.]+)\s+[0-9.]+\s*'
            r'\(\s*(-?[0-9.]+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)\)',
            ir_section.group(1), re.MULTILINE
        )
        for match in matches:
            freq = float(match[1])
            if freq > 0:
                modes.append(IRMode(
                    index=int(match[0]),
                    frequency=freq,
                    intensity=float(match[2]),
                    tx=float(match[3]),
                    ty=float(match[4]),
                    tz=float(match[5])
                ))

    return modes


def parse_mulliken_charges(content: str) -> dict[int, tuple[str, float]]:
    """Extract Mulliken atomic charges."""
    charges = {}

    sections = re.findall(
        r'MULLIKEN ATOMIC CHARGES\s*-+\s*(.*?)Sum of atomic charges',
        content, re.DOTALL
    )

    if sections:
        last_section = sections[-1]
        matches = re.findall(
            r'(\d+)\s+(\w+)\s+:\s+(-?\d+\.?\d*)',
            last_section
        )
        for match in matches:
            idx = int(match[0])
            element = match[1]
            charge = float(match[2])
            charges[idx] = (element, charge)

    return charges


def parse_loewdin_charges(content: str) -> dict[int, tuple[str, float]]:
    """Extract Loewdin atomic charges."""
    charges = {}

    sections = re.findall(
        r'LOEWDIN ATOMIC CHARGES\s*-+\s*(.*?)(?=Total charge|LOEWDIN REDUCED)',
        content, re.DOTALL
    )

    if sections:
        last_section = sections[-1]
        matches = re.findall(
            r'(\d+)\s+(\w+)\s+:\s+(-?\d+\.?\d*)',
            last_section
        )
        for match in matches:
            idx = int(match[0])
            element = match[1]
            charge = float(match[2])
            charges[idx] = (element, charge)

    return charges


def parse_mayer_bond_orders(content: str) -> list[tuple[int, int, float]]:
    """Extract Mayer bond orders."""
    bond_orders = []

    # Find all bond order entries
    matches = re.findall(
        r'B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s+(\d+\.?\d*)',
        content
    )

    # Use a set to avoid duplicates from multiple sections
    seen = set()
    for match in matches:
        atom1 = int(match[0])
        atom2 = int(match[1])
        order = float(match[2])
        key = (min(atom1, atom2), max(atom1, atom2))
        if key not in seen and order > 0.1:  # Filter weak bonds
            seen.add(key)
            bond_orders.append((atom1, atom2, order))

    return bond_orders


def parse_loewdin_bond_orders(content: str) -> list[tuple[int, int, float]]:
    """Extract Loewdin bond orders."""
    bond_orders = []

    # Find the Loewdin bond orders section
    loewdin_section = re.search(
        r'LOEWDIN BOND ORDERS.*?-+\s*(.*?)(?:\n\n|--+)',
        content, re.DOTALL
    )

    if loewdin_section:
        # Find all bond order entries in the Loewdin section only
        matches = re.findall(
            r'B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s+(\d+\.?\d*)',
            loewdin_section.group(1)
        )

        for match in matches:
            atom1 = int(match[0])
            atom2 = int(match[1])
            order = float(match[2])
            if order > 0.05:  # Loewdin uses 0.05 threshold
                bond_orders.append((atom1, atom2, order))

    return bond_orders


def parse_mulliken_overlap_charges(content: str) -> list[tuple[int, int, float]]:
    """Extract Mulliken overlap charges (bonding analysis)."""
    overlap_charges = []

    # Find the Mulliken overlap charges section
    overlap_section = re.search(
        r'MULLIKEN OVERLAP CHARGES\s*-+\s*(.*?)(?:\n\n|$)',
        content, re.DOTALL
    )

    if not overlap_section:
        return overlap_charges

    section_text = overlap_section.group(1)

    # Find all overlap charge entries: B(  0-C ,  1-H ) :   1.1736
    matches = re.findall(
        r'B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s+(-?\d+\.?\d*)',
        section_text
    )

    for match in matches:
        atom1 = int(match[0])
        atom2 = int(match[1])
        overlap = float(match[2])
        overlap_charges.append((atom1, atom2, overlap))

    return overlap_charges


def parse_mulliken_orbital_charges(content: str) -> list[MOCharge]:
    """Extract Mulliken orbital charges (per-MO charge distribution)."""
    charges = []

    # Find MULLIKEN ORBITAL CHARGES section
    section = re.search(
        r'MULLIKEN ORBITAL CHARGES\s*-+\s*.*?\n(.*?)(?:\n\n|LOEWDIN)',
        content, re.DOTALL
    )

    if not section:
        return charges

    # Parse format: "   0:   0C   1s           1.103241 (  0.658373)"
    # MO_index: atom_index+element orbital charge (uncorrected_charge)
    matches = re.findall(
        r'^\s*(\d+):\s+(\d+)([A-Z][a-z]?)\s+(\S+)\s+(-?\d+\.?\d+)\s+\(\s*(-?\d+\.?\d+)\)',
        section.group(1), re.MULTILINE
    )

    for mo_idx, atom_idx, element, orbital, charge, uncorrected in matches:
        charges.append(MOCharge(
            mo_index=int(mo_idx),
            atom_index=int(atom_idx),
            element=element,
            orbital=orbital,
            charge=float(charge),
            uncorrected_charge=float(uncorrected)
        ))

    return charges


def parse_loewdin_orbital_charges(content: str) -> list[MOCharge]:
    """Extract Loewdin orbital charges (per-MO charge distribution)."""
    charges = []

    # Find LOEWDIN ORBITAL CHARGES section
    section = re.search(
        r'LOEWDIN ORBITAL CHARGES\s*-+\s*(.*?)(?:\n\n|$)',
        content, re.DOTALL
    )

    if not section:
        return charges

    # Parse format: "   0:   0C   1s           1.046092"
    # MO_index: atom_index+element orbital charge
    matches = re.findall(
        r'^\s*(\d+):\s+(\d+)([A-Z][a-z]?)\s+(\S+)\s+(-?\d+\.?\d+)',
        section.group(1), re.MULTILINE
    )

    for mo_idx, atom_idx, element, orbital, charge in matches:
        charges.append(MOCharge(
            mo_index=int(mo_idx),
            atom_index=int(atom_idx),
            element=element,
            orbital=orbital,
            charge=float(charge),
            uncorrected_charge=None  # Loewdin doesn't have uncorrected charges
        ))

    return charges


def parse_thermochemistry(content: str) -> Optional[Thermochemistry]:
    """Extract thermochemistry data."""
    thermo = Thermochemistry()

    # Temperature and pressure
    temp_match = re.search(r'Temperature\s+\.\.\.\s+(\d+\.?\d*)\s+K', content)
    if temp_match:
        thermo.temperature = float(temp_match.group(1))

    press_match = re.search(r'Pressure\s+\.\.\.\s+(\d+\.?\d*)\s+atm', content)
    if press_match:
        thermo.pressure = float(press_match.group(1))

    mass_match = re.search(r'Total Mass\s+\.\.\.\s+(\d+\.?\d*)\s+AMU', content)
    if mass_match:
        thermo.total_mass = float(mass_match.group(1))

    # Zero point energy
    zpe_match = re.search(r'Zero point energy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if zpe_match:
        thermo.zero_point_energy = float(zpe_match.group(1))

    # Inner energy
    inner_match = re.search(r'Total thermal energy\s+(-?\d+\.?\d*)\s+Eh', content)
    if inner_match:
        thermo.inner_energy = float(inner_match.group(1))

    # Total enthalpy
    enthalpy_match = re.search(r'Total enthalpy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if enthalpy_match:
        thermo.enthalpy = float(enthalpy_match.group(1))

    # Entropy components
    elec_s = re.search(r'Electronic entropy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if elec_s:
        thermo.entropy_electronic = float(elec_s.group(1))

    vib_s = re.search(r'Vibrational entropy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if vib_s:
        thermo.entropy_vibrational = float(vib_s.group(1))

    rot_s = re.search(r'Rotational entropy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if rot_s:
        thermo.entropy_rotational = float(rot_s.group(1))

    trans_s = re.search(r'Translational entropy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if trans_s:
        thermo.entropy_translational = float(trans_s.group(1))

    # Total entropy
    entropy_match = re.search(r'Final entropy term\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if entropy_match:
        thermo.entropy = float(entropy_match.group(1))

    # Final Gibbs free energy
    gibbs_match = re.search(r'Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if gibbs_match:
        thermo.gibbs_free_energy = float(gibbs_match.group(1))

    # Only return if we found thermochemistry data
    if thermo.total_mass > 0:
        return thermo
    return None


def parse_nmr_data(content: str) -> Optional[NMRData]:
    """Extract NMR chemical shifts and J-couplings."""
    nmr = NMRData()

    # Parse chemical shift summary
    shift_section = re.search(
        r'CHEMICAL SHIELDING SUMMARY \(ppm\)\s*-+\s*(.*?)(?=-{20,}|\Z)',
        content, re.DOTALL
    )

    if shift_section:
        # Match: Nucleus, Element, Isotropic, Anisotropy
        matches = re.findall(
            r'^\s*(\d+)\s+(\w+)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)',
            shift_section.group(1), re.MULTILINE
        )
        for match in matches:
            nmr.chemical_shifts.append(NMRShift(
                atom_index=int(match[0]),
                element=match[1],
                isotropic=float(match[2]),
                anisotropy=float(match[3])
            ))

    # Parse J-coupling summary
    jcoupling_section = re.search(
        r'SUMMARY OF ISOTROPIC COUPLING CONSTANTS.*?-+\s*(.*?)(?=-{20,}|\Z)',
        content, re.DOTALL
    )

    if jcoupling_section:
        # Match pairs of nuclei and J values
        matches = re.findall(
            r'(\d+)\s+(\d+)\s+(-?\d+\.?\d*)',
            jcoupling_section.group(1)
        )
        for match in matches:
            nmr.j_couplings.append(JCoupling(
                atom1_index=int(match[0]),
                atom2_index=int(match[1]),
                value=float(match[2])
            ))

    # Only return if we found NMR data
    if nmr.chemical_shifts or nmr.j_couplings:
        return nmr
    return None


def parse_chemical_shielding_tensors(content: str) -> list[ChemicalShieldingTensor]:
    """Extract full chemical shielding tensors (diamagnetic, paramagnetic, total)."""
    tensors = []

    # Find all nucleus sections in CHEMICAL SHIELDINGS (ppm) block
    # Pattern: "Nucleus   1H :" or "Nucleus   0C :"
    nucleus_pattern = r'-{10,}\s*Nucleus\s+(\d+)([A-Z][a-z]?)\s*:\s*-{10,}\s*(.*?)(?=-{10,}\s*Nucleus|\Z)'

    nucleus_matches = re.findall(nucleus_pattern, content, re.DOTALL)

    for atom_idx, element, nucleus_content in nucleus_matches:
        tensor = ChemicalShieldingTensor(
            atom_index=int(atom_idx),
            element=element
        )

        # Parse diamagnetic tensor (3x3 matrix)
        dia_match = re.search(
            r'Diamagnetic contribution.*?:\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)',
            nucleus_content
        )
        if dia_match:
            tensor.diamagnetic_tensor = [float(dia_match.group(i)) for i in range(1, 10)]

        # Parse paramagnetic tensor (3x3 matrix)
        para_match = re.search(
            r'Paramagnetic contribution.*?:\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)',
            nucleus_content
        )
        if para_match:
            tensor.paramagnetic_tensor = [float(para_match.group(i)) for i in range(1, 10)]

        # Parse total shielding tensor (3x3 matrix)
        total_match = re.search(
            r'Total shielding tensor.*?:\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)',
            nucleus_content
        )
        if total_match:
            tensor.total_tensor = [float(total_match.group(i)) for i in range(1, 10)]

        # Parse diagonalized components
        # Format: sDSO  11.490  44.520  27.669  iso=  27.893
        diag_match = re.search(
            r'sDSO\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+iso=\s*(-?\d+\.?\d+)\s*\n\s*sPSO\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+iso=\s*(-?\d+\.?\d+)\s*\n.*?Total\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+iso=\s*(-?\d+\.?\d+)',
            nucleus_content, re.DOTALL
        )
        if diag_match:
            tensor.sdso_components = [float(diag_match.group(i)) for i in [1, 2, 3]]
            tensor.sdso_iso = float(diag_match.group(4))
            tensor.spso_components = [float(diag_match.group(i)) for i in [5, 6, 7]]
            tensor.spso_iso = float(diag_match.group(8))
            tensor.total_components = [float(diag_match.group(i)) for i in [9, 10, 11]]
            tensor.total_iso = float(diag_match.group(12))

        # Parse orientation eigenvectors
        # Format: X  -0.0461294  0.0468023  -0.9978385
        orient_match = re.search(
            r'Orientation:\s*\n\s*X\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*Y\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s*\n\s*Z\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.?\d+)',
            nucleus_content
        )
        if orient_match:
            tensor.orientation_x = [float(orient_match.group(i)) for i in [1, 2, 3]]
            tensor.orientation_y = [float(orient_match.group(i)) for i in [4, 5, 6]]
            tensor.orientation_z = [float(orient_match.group(i)) for i in [7, 8, 9]]

        tensors.append(tensor)

    return tensors


def parse_geometric_perturbations(content: str) -> Optional[GeometricPerturbations]:
    """Extract geometric perturbations data."""
    geom_pert = GeometricPerturbations()

    section = re.search(
        r'GEOMETRIC PERTURBATIONS \((\d+) nuclei\)\s*-+(.*?)(?:Maximum memory|Property integrals)',
        content, re.DOTALL
    )

    if not section:
        return None

    geom_pert.num_nuclei = int(section.group(1))

    # Extract number of perturbations from BATCH line
    batch_match = re.search(r'BATCH.*?\(\s*(\d+)\s+perturbations\)', section.group(2))
    if batch_match:
        geom_pert.num_perturbations = int(batch_match.group(1))

    # Extract total time
    time_match = re.search(r'geometrical perturbations done \(\s*(\d+\.?\d*)\s+sec\)', section.group(2))
    if time_match:
        geom_pert.total_time = float(time_match.group(1))

    # Extract max memory
    mem_match = re.search(r'MaxCore\s+\.\.\.\s+(\d+)\s+MB', section.group(2))
    if mem_match:
        geom_pert.max_memory_mb = float(mem_match.group(1))

    # Extract number of batches
    batch_count = re.search(r'Number of batches\s+\.\.\.\s+(\d+)', section.group(2))
    if batch_count:
        geom_pert.num_batches = int(batch_count.group(1))

    return geom_pert


def parse_shark_integrals(content: str) -> Optional[SharkIntegrals]:
    """Extract SHARK integral package configuration."""
    shark = SharkIntegrals()

    section = re.search(
        r'SHARK INTEGRAL PACKAGE\s*-+(.*?)(?:Checking pre-screening|$)',
        content, re.DOTALL
    )

    if not section:
        return None

    section_text = section.group(1)

    # Extract basic info
    patterns = {
        'num_basis_functions': r'Number of basis functions\s+\.\.\.\s+(\d+)',
        'num_shells': r'Number of shells\s+\.\.\.\s+(\d+)',
        'max_angular_momentum': r'Maximum angular momentum\s+\.\.\.\s+(\d+)',
        'num_aux_j_functions': r'# of basis functions in Aux-J\s+\.\.\.\s+(\d+)',
        'num_aux_j_shells': r'# of shells in Aux-J\s+\.\.\.\s+(\d+)',
        'max_angular_momentum_aux_j': r'Maximum angular momentum in Aux-J\s+\.\.\.\s+(\d+)',
        'total_shell_pairs': r'Total number of shell pairs\s+\.\.\.\s+(\d+)',
        'shell_pairs_after_screening': r'Shell pairs after pre-screening\s+\.\.\.\s+(\d+)',
        'schwartz_threshold': r'Thresh\s+\.\.\.\s+(\d+\.\d+e[+-]?\d+)',
        'tcut_threshold': r'Tcut\s+\.\.\.\s+(\d+\.\d+e[+-]?\d+)'
    }

    for field, pattern in patterns.items():
        match = re.search(pattern, section_text)
        if match:
            value = match.group(1)
            if field in ['schwartz_threshold', 'tcut_threshold']:
                setattr(shark, field, float(value))
            else:
                setattr(shark, field, int(value))

    return shark


def parse_cosx_grids(content: str) -> list[COSXGrid]:
    """Extract COSX grid generation data."""
    grids = []

    # Find all GRIDX sections
    grid_sections = re.findall(
        r'GRIDX\s+(\d+)\s*-+(.*?)(?=GRIDX\s+\d+|$)',
        content, re.DOTALL
    )

    for grid_num, grid_content in grid_sections:
        grid = COSXGrid(grid_number=int(grid_num))

        # Extract grid properties
        acc_match = re.search(r'General Integration Accuracy\s+IntAcc\s+\.\.\.\s+(\d+\.?\d*)', grid_content)
        if acc_match:
            grid.integration_accuracy = float(acc_match.group(1))

        rad_match = re.search(r'Radial Grid Type\s+RadialGrid\s+\.\.\.\s+(.*?)(?:\n|$)', grid_content)
        if rad_match:
            grid.radial_grid_type = rad_match.group(1).strip()

        ang_match = re.search(r'Angular Grid.*?AngularGrid\s+\.\.\.\s+(.*?)(?:\n|$)', grid_content)
        if ang_match:
            grid.angular_grid = ang_match.group(1).strip()

        points_match = re.search(r'Total number of grid points\s+\.\.\.\s+(\d+)', grid_content)
        if points_match:
            grid.total_grid_points = int(points_match.group(1))

        batch_match = re.search(r'Total number of batches\s+\.\.\.\s+(\d+)', grid_content)
        if batch_match:
            grid.num_batches = int(batch_match.group(1))

        avg_batch_match = re.search(r'Average number of points per batch\s+\.\.\.\s+(\d+)', grid_content)
        if avg_batch_match:
            grid.avg_points_per_batch = int(avg_batch_match.group(1))

        avg_atom_match = re.search(r'Average number of grid points per atom\s+\.\.\.\s+(\d+)', grid_content)
        if avg_atom_match:
            grid.avg_points_per_atom = int(avg_atom_match.group(1))

        grids.append(grid)

    return grids


def parse_raman_spectrum(content: str) -> list[RamanMode]:
    """Extract Raman spectrum data."""
    modes = []

    # Find the Raman section and extract all mode data
    raman_section = re.search(
        r'RAMAN SPECTRUM.*?-{50,}\s*(.*?)(?:\n\s*\n|\Z)',
        content, re.DOTALL
    )

    if raman_section:
        # Match: Mode, freq (cm^-1), Activity, Depolarization
        matches = re.findall(
            r'^\s*(\d+):\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)',
            raman_section.group(1), re.MULTILINE
        )
        for match in matches:
            freq = float(match[1])
            if freq > 0:
                modes.append(RamanMode(
                    index=int(match[0]),
                    frequency=freq,
                    activity=float(match[2]),
                    depolarization=float(match[3])
                ))

    return modes


def parse_dispersion_correction(content: str) -> Optional[DispersionCorrection]:
    """Extract DFT dispersion correction data."""
    disp = DispersionCorrection()

    # Look for DFTD3 section with format: Edisp/kcal,au: value1 value2
    disp_section = re.search(
        r'DFT DISPERSION CORRECTION.*?DFTD(\d+).*?'
        r'Edisp/kcal,au:\s*(-?\d+\.?\d*)\s+(-?\d+\.?\d*).*?'
        r'E6\s+/kcal\s+:\s*(-?\d+\.?\d*).*?'
        r'E8\s+/kcal\s+:\s*(-?\d+\.?\d*)',
        content, re.DOTALL | re.IGNORECASE
    )

    if disp_section:
        disp.method = f"DFTD{disp_section.group(1)}"
        disp.total_correction_kcal = float(disp_section.group(2))
        disp.total_correction = float(disp_section.group(3))
        disp.e6_kcal = float(disp_section.group(4))
        disp.e8_kcal = float(disp_section.group(5))
        return disp

    return None


def parse_normal_modes(content: str, num_atoms: int) -> list[NormalMode]:
    """Extract normal mode displacement vectors."""
    modes = []

    # Find NORMAL MODES section
    normal_section = re.search(
        r'NORMAL MODES.*?weighted by.*?\n\s*(.*?)(?:\n\s*\n|\Z)',
        content, re.DOTALL | re.IGNORECASE
    )

    if not normal_section:
        return modes

    section_text = normal_section.group(1)
    lines = section_text.strip().split('\n')

    # Parse mode indices from header line
    if not lines:
        return modes

    # First line should have mode indices
    header_match = re.findall(r'\d+', lines[0])
    if not header_match:
        return modes

    mode_indices = [int(x) for x in header_match]
    num_modes = len(mode_indices)

    # Initialize modes
    mode_data = {idx: [] for idx in mode_indices}

    # Parse displacement data (3 lines per atom: X, Y, Z)
    line_idx = 1
    for atom_idx in range(num_atoms):
        if line_idx + 2 >= len(lines):
            break

        for coord_idx in range(3):  # X, Y, Z
            if line_idx >= len(lines):
                break

            parts = lines[line_idx].split()
            if len(parts) >= num_modes + 1:
                for mode_offset, mode_id in enumerate(mode_indices):
                    try:
                        value = float(parts[mode_offset + 1])
                        if coord_idx == 0:  # First coordinate for this atom
                            if len(mode_data[mode_id]) <= atom_idx:
                                mode_data[mode_id].append([value, 0.0, 0.0])
                            else:
                                mode_data[mode_id][atom_idx][0] = value
                        elif coord_idx == 1:  # Y
                            mode_data[mode_id][atom_idx][1] = value
                        else:  # Z
                            mode_data[mode_id][atom_idx][2] = value
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Failed to parse normal mode displacement for mode {mode_id}, "
                                     f"atom {atom_idx}, coord {coord_idx}: {e}")

            line_idx += 1

    # Convert to NormalMode objects
    # Need to get frequencies for each mode
    freq_match = re.search(r'VIBRATIONAL FREQUENCIES.*?\n-+\s*(.*?)(?:\n\s*\n|\Z)', content, re.DOTALL)
    frequencies = {}
    if freq_match:
        freq_lines = freq_match.group(1).strip().split('\n')
        for line in freq_lines:
            parts = line.split(':')
            if len(parts) == 2:
                try:
                    mode_num = int(parts[0].strip())
                    freq_str = parts[1].strip().split()[0]
                    frequencies[mode_num] = float(freq_str)
                except (ValueError, IndexError) as e:
                    logger.warning(f"Failed to parse frequency from line '{line.strip()}': {e}")

    for mode_id in mode_indices:
        if mode_id in mode_data:
            modes.append(NormalMode(
                mode_index=mode_id,
                frequency=frequencies.get(mode_id, 0.0),
                displacements=[tuple(disp) for disp in mode_data[mode_id]]
            ))

    return modes


def parse_scf_iterations(content: str) -> list[SCFIteration]:
    """Extract SCF convergence data."""
    iterations = []

    # Find all SCF iteration tables
    scf_sections = re.finditer(
        r'Iteration\s+Energy.*?Delta-E.*?RMSDP.*?MaxDP.*?DIISErr.*?\n(.*?)(?:\n\*+|Total SCF|SCF CONVERGED|\Z)',
        content, re.DOTALL | re.IGNORECASE
    )

    for section in scf_sections:
        section_text = section.group(1)
        # Match iteration lines: number, energy, delta-E, RMSDP, MaxDP, DIISErr
        matches = re.findall(
            r'^\s*(\d+)\s+(-?\d+\.?\d+)\s+(-?\d+\.\d+[eE][-+]?\d+)\s+'
            r'(\d+\.\d+[eE][-+]?\d+)\s+(\d+\.\d+[eE][-+]?\d+)\s+(\d+\.\d+[eE][-+]?\d+)',
            section_text, re.MULTILINE
        )

        for match in matches:
            try:
                iterations.append(SCFIteration(
                    iteration=int(match[0]),
                    energy=float(match[1]),
                    delta_e=float(match[2]),
                    rmsdp=float(match[3]),
                    maxdp=float(match[4]),
                    diis_error=float(match[5])
                ))
            except ValueError as e:
                logger.warning(f"Failed to parse SCF iteration data from match {match}: {e}")

    return iterations


def parse_dft_grid_info(content: str) -> Optional[DFTGridInfo]:
    """Extract DFT grid parameters."""
    grid_section = re.search(
        r'DFT GRID GENERATION\s*-+(.*?)(?:\n-{20,}|\Z)',
        content, re.DOTALL | re.IGNORECASE
    )

    if not grid_section:
        return None

    section_text = grid_section.group(1)
    grid = DFTGridInfo()

    # Parse integration accuracy
    acc_match = re.search(r'General Integration Accuracy\s+IntAcc\s+\.\.\.\s+(\d+\.?\d*)', section_text)
    if acc_match:
        grid.integration_accuracy = float(acc_match.group(1))

    # Parse radial grid type
    radial_match = re.search(r'Radial Grid Type\s+RadialGrid\s+\.\.\.\s+(.+)', section_text)
    if radial_match:
        grid.radial_grid_type = radial_match.group(1).strip()

    # Parse angular grid
    angular_match = re.search(r'Angular Grid \(max\. ang\.\)\s+AngularGrid\s+\.\.\.\s+(.+)', section_text)
    if angular_match:
        grid.angular_grid = angular_match.group(1).strip()

    # Parse pruning method
    pruning_match = re.search(r'Angular grid pruning method\s+GridPruning\s+\.\.\.\s+(.+)', section_text)
    if pruning_match:
        grid.pruning_method = pruning_match.group(1).strip()

    # Parse total grid points
    points_match = re.search(r'Total number of grid points\s+\.\.\.\s+(\d+)', section_text)
    if points_match:
        grid.total_grid_points = int(points_match.group(1))

    # Parse average points per atom
    avg_match = re.search(r'Average number of grid points per atom\s+\.\.\.\s+(\d+)', section_text)
    if avg_match:
        grid.avg_points_per_atom = int(avg_match.group(1))

    if grid.total_grid_points > 0:
        return grid
    return None


def parse_basis_set_info(content: str) -> Optional[BasisSetInfo]:
    """Extract basis set information."""
    basis = BasisSetInfo()

    # Parse basis set name
    basis_match = re.search(r'Your calculation utilizes the basis:\s*(\S+)', content)
    if basis_match:
        basis.name = basis_match.group(1)

    # Parse number of basis functions
    nbf_match = re.search(r'Number of basis functions\s+\.+\s+(\d+)', content, re.IGNORECASE)
    if nbf_match:
        basis.num_basis_functions = int(nbf_match.group(1))

    # Parse number of primitive gaussians
    nprim_match = re.search(r'Number of primitive gaussian shells\s+\.+\s+(\d+)', content, re.IGNORECASE)
    if nprim_match:
        basis.num_primitive_gaussians = int(nprim_match.group(1))

    # Parse basis set groups per element
    basis_section = re.search(
        r'BASIS SET INFORMATION\s*-+\s*There are \d+ groups of distinct atoms\s*(.*?)(?:\n\n|Atom\s+\d+)',
        content, re.DOTALL
    )

    if basis_section:
        group_text = basis_section.group(1)
        # Match: Group   1 Type C   : description
        groups = re.findall(r'Group\s+\d+\s+Type\s+(\w+)\s*:\s*(.+)', group_text)
        for element, description in groups:
            basis.element_basis[element] = description.strip()

    if basis.num_basis_functions > 0:
        return basis
    return None


def parse_energy_components(content: str) -> Optional[EnergyComponents]:
    """Extract detailed energy breakdown."""
    energy = EnergyComponents()

    # Find energy components section
    energy_section = re.search(
        r'Nuclear Repulsion\s*:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'Electronic Energy\s*:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'One Electron Energy:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'Two Electron Energy:\s*(-?\d+\.?\d+)\s*Eh',
        content, re.DOTALL
    )

    if energy_section:
        energy.nuclear_repulsion = float(energy_section.group(1))
        energy.electronic_energy = float(energy_section.group(2))
        energy.one_electron_energy = float(energy_section.group(3))
        energy.two_electron_energy = float(energy_section.group(4))

    # Parse virial components
    virial_section = re.search(
        r'Potential Energy\s*:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'Kinetic Energy\s*:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'Virial Ratio\s*:\s*(-?\d+\.?\d+)',
        content, re.DOTALL
    )

    if virial_section:
        energy.potential_energy = float(virial_section.group(1))
        energy.kinetic_energy = float(virial_section.group(2))
        energy.virial_ratio = float(virial_section.group(3))

    # Parse DFT components
    xc_section = re.search(
        r'E\(X\)\s*:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'E\(C\)\s*:\s*(-?\d+\.?\d+)\s*Eh.*?'
        r'E\(XC\)\s*:\s*(-?\d+\.?\d+)\s*Eh',
        content, re.DOTALL
    )

    if xc_section:
        energy.exchange_energy = float(xc_section.group(1))
        energy.correlation_energy = float(xc_section.group(2))
        energy.xc_energy = float(xc_section.group(3))

    # Return if we found at least some data
    if energy.nuclear_repulsion != 0.0 or energy.electronic_energy != 0.0:
        return energy
    return None


def parse_cpcm_solvation(content: str) -> Optional[CPCMSolvation]:
    """Extract CPCM solvation model properties."""
    cpcm = CPCMSolvation()

    # Find CPCM Solvation Model Properties section
    cpcm_section = re.search(
        r'CPCM Solvation Model Properties:\s*'
        r'Surface-charge\s*:\s*(-?\d+\.?\d*)\s*'
        r'Corrected charge\s*:\s*(-?\d+\.?\d*)\s*'
        r'Outlying charge corr\.\s*:\s*(-?\d+\.?\d+)\s*Eh\s*(-?\d+\.?\d+)\s*eV',
        content, re.DOTALL
    )

    if cpcm_section:
        cpcm.surface_charge = float(cpcm_section.group(1))
        cpcm.corrected_charge = float(cpcm_section.group(2))
        cpcm.outlying_charge_correction = float(cpcm_section.group(3))
        cpcm.outlying_charge_correction_ev = float(cpcm_section.group(4))

    # Also extract CPCM Dielectric energy from components section
    dielectric_match = re.search(r'CPCM Dielectric\s*:\s*(-?\d+\.?\d+)\s*Eh', content)
    if dielectric_match:
        cpcm.dielectric_energy = float(dielectric_match.group(1))

    # Return if we found at least some data
    if cpcm.surface_charge != 0.0 or cpcm.dielectric_energy != 0.0:
        return cpcm
    return None


def parse_scf_convergence(content: str) -> Optional[SCFConvergence]:
    """Extract SCF convergence criteria values."""
    conv = SCFConvergence()

    # Find SCF CONVERGENCE section
    conv_section = re.search(
        r'SCF CONVERGENCE\s*-+\s*'
        r'Last Energy change\s*\.\.\.\s*(-?\d+\.?\d+e[+-]?\d+).*?'
        r'Last MAX-Density change\s*\.\.\.\s*(-?\d+\.?\d+e[+-]?\d+).*?'
        r'Last RMS-Density change\s*\.\.\.\s*(-?\d+\.?\d+e[+-]?\d+).*?'
        r'Last DIIS Error\s*\.\.\.\s*(-?\d+\.?\d+e[+-]?\d+).*?'
        r'Last Orbital Gradient\s*\.\.\.\s*(-?\d+\.?\d+e[+-]?\d+).*?'
        r'Last Orbital Rotation\s*\.\.\.\s*(-?\d+\.?\d+e[+-]?\d+)',
        content, re.DOTALL
    )

    if conv_section:
        conv.energy_change = float(conv_section.group(1))
        conv.max_density_change = float(conv_section.group(2))
        conv.rms_density_change = float(conv_section.group(3))
        conv.diis_error = float(conv_section.group(4))
        conv.orbital_gradient = float(conv_section.group(5))
        conv.orbital_rotation = float(conv_section.group(6))
        return conv

    return None


def parse_mulliken_orbital_populations(content: str) -> list[OrbitalPopulation]:
    """Extract Mulliken reduced orbital charges (orbital population analysis)."""
    populations = []

    # Find the MULLIKEN REDUCED ORBITAL CHARGES section
    section_match = re.search(
        r'MULLIKEN REDUCED ORBITAL CHARGES\s*-+\s*(.*?)(?:\n\n\n)',
        content, re.DOTALL
    )

    if not section_match:
        return populations

    section_text = section_match.group(1)
    lines = section_text.strip().split('\n')

    current_pop = None

    for line in lines:
        # Check if this line starts a new atom: "  0 C s       :     3.171400  s :     3.171400"
        atom_match = re.match(r'^\s*(\d+)\s+(\w+)\s+', line)
        if atom_match:
            # Save previous atom if any
            if current_pop:
                populations.append(current_pop)

            # Start new atom
            atom_idx = int(atom_match.group(1))
            element = atom_match.group(2)
            current_pop = OrbitalPopulation(atom_index=atom_idx, element=element)

            # Also process the rest of this line (contains first orbital data)
            # Extract totals from the same line
            for orbital_type in ['s', 'p', 'd', 'f', 'g']:
                match = re.search(rf'{orbital_type}\s*:\s*(\d+\.?\d+)', line)
                if match:
                    setattr(current_pop, f'{orbital_type}_total', float(match.group(1)))

        elif current_pop:
            # This is a continuation line for the current atom
            # Extract orbital details and totals
            # Format: "      pz      :     0.990710  p :     2.888037"
            # or just: "      px      :     0.885174"

            # Check if line has both individual orbital and total
            # Pattern: "orbital_name : value  orbital_type : total"
            combined_match = re.match(r'([a-z][a-z0-9+\-]*)\s*:\s*(\d+\.?\d+)\s+([a-z])\s*:\s*(\d+\.?\d+)', line.strip())
            if combined_match:
                # Has both individual orbital and total
                orb_name = combined_match.group(1)
                orb_value = float(combined_match.group(2))
                orb_type = combined_match.group(3)
                orb_total = float(combined_match.group(4))

                if orb_name not in ['s', 'p', 'd', 'f', 'g']:
                    current_pop.orbital_details[orb_name] = orb_value
                setattr(current_pop, f'{orb_type}_total', orb_total)
            else:
                # Just has individual orbital
                orb_match = re.match(r'([a-z][a-z0-9+\-]*)\s*:\s*(\d+\.?\d+)', line.strip())
                if orb_match:
                    orb_name = orb_match.group(1)
                    orb_value = float(orb_match.group(2))
                    if orb_name not in ['s', 'p', 'd', 'f', 'g']:
                        current_pop.orbital_details[orb_name] = orb_value

    # Don't forget the last atom
    if current_pop:
        populations.append(current_pop)

    return populations


def parse_timing_data(content: str) -> Optional[TimingData]:
    """Extract computational timing breakdown."""
    timing = TimingData()

    # Find TIMINGS section (usually the last one)
    timing_sections = list(re.finditer(r'TIMINGS\s*-+\s*(.*?)(?:\n\*+|\Z)', content, re.DOTALL | re.IGNORECASE))
    if not timing_sections:
        return None

    # Use the last timing section (final summary)
    section_text = timing_sections[-1].group(1)

    # Parse total time
    total_match = re.search(r'Total time\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec', section_text)
    if total_match:
        timing.total_time = float(total_match.group(1))

    # Parse TOTAL RUN TIME at the very end of the file
    runtime_match = re.search(
        r'TOTAL RUN TIME:\s+(\d+)\s+days?\s+(\d+)\s+hours?\s+(\d+)\s+minutes?\s+(\d+)\s+seconds?\s+(\d+)\s+msec',
        content
    )
    if runtime_match:
        days = int(runtime_match.group(1))
        hours = int(runtime_match.group(2))
        minutes = int(runtime_match.group(3))
        seconds = int(runtime_match.group(4))
        msec = int(runtime_match.group(5))
        # Convert to total seconds
        timing.total_run_time = days * 86400 + hours * 3600 + minutes * 60 + seconds + msec / 1000.0

    # Parse SCF time
    scf_match = re.search(r'Total SCF time:\s+\d+\s+days?\s+\d+\s+hours?\s+\d+\s+min\s+(\d+)', section_text)
    if scf_match:
        # Convert to seconds - already in section as total time
        pass

    # Parse individual timing components
    timing_patterns = {
        'scf_preparation': r'SCF preparation\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec',
        'fock_matrix': r'Fock matrix formation\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec',
        'diagonalization': r'Diagonalization\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec',
        'density_formation': r'Density matrix formation\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec',
        'population_analysis': r'Population analysis\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec',
        'diis_solution': r'DIIS solution\s+\.\.\.\.\s+(\d+\.?\d*)\s+sec',
    }

    for key, pattern in timing_patterns.items():
        match = re.search(pattern, section_text)
        if match:
            timing.timings[key] = float(match.group(1))
            if key == 'fock_matrix':
                timing.fock_matrix_time = float(match.group(1))

    if timing.total_time > 0:
        return timing
    return None


if __name__ == '__main__':
    import sys
    import json

    if len(sys.argv) > 1:
        result = parse_out_file(sys.argv[1])
        print(f"Method: {result.job_info.method}")
        print(f"Basis: {result.job_info.basis_set}")
        print(f"Charge: {result.job_info.charge}, Mult: {result.job_info.multiplicity}")
        print(f"Final Energy: {result.final_energy:.6f} Eh")

        if result.coordinates:
            print(f"Coordinates: {len(result.coordinates)} atoms")

        if result.coordinates_au:
            print(f"Coordinates (a.u.): {len(result.coordinates_au)} atoms with atomic numbers and masses")

        if result.internal_coords:
            print(f"Internal Coordinates: {len(result.internal_coords)} atoms (Z-matrix)")
            if len(result.internal_coords) > 1:
                # Show second atom as example (first atom has no bonds)
                ic = result.internal_coords[1]
                print(f"  Example: {ic.element} bonded to atom {ic.bond_to}, r={ic.bond_length:.4f} Å")

        if result.dipole_moment:
            print(f"Dipole: {result.dipole_moment.magnitude_debye:.3f} Debye")

        if result.polarizability:
            print(f"Isotropic Polarizability: {result.polarizability.isotropic:.2f} a.u.")
            if result.polarizability.eigenvalues:
                eigs = result.polarizability.eigenvalues
                print(f"  Eigenvalues: {eigs[0]:.2f}, {eigs[1]:.2f}, {eigs[2]:.2f} a.u.")

        if result.orbital_energies:
            occupied = [o for o in result.orbital_energies if o.occupation > 0]
            virtual = [o for o in result.orbital_energies if o.occupation == 0]
            print(f"Orbitals: {len(occupied)} occupied, {len(virtual)} virtual")
            if occupied:
                homo = occupied[-1]
                print(f"  HOMO: {homo.energy_ev:.3f} eV")
            if virtual:
                lumo = virtual[0]
                print(f"  LUMO: {lumo.energy_ev:.3f} eV")

        if result.frequencies:
            print(f"Frequencies: {len(result.frequencies)} modes")
            print(f"  Range: {min(result.frequencies):.1f} - {max(result.frequencies):.1f} cm^-1")

        if result.ir_spectrum:
            strongest = max(result.ir_spectrum, key=lambda x: x.intensity)
            print(f"IR Spectrum: {len(result.ir_spectrum)} modes")
            print(f"  Strongest: {strongest.frequency:.1f} cm^-1 ({strongest.intensity:.1f} km/mol)")

        if result.mayer_bond_orders:
            print(f"Mayer Bond Orders: {len(result.mayer_bond_orders)} bonds")

        if result.loewdin_bond_orders:
            print(f"Loewdin Bond Orders: {len(result.loewdin_bond_orders)} bonds")

        if result.mulliken_overlap_charges:
            print(f"Mulliken Overlap Charges: {len(result.mulliken_overlap_charges)} atom pairs")
            # Show a strong bonding example
            strong = max(result.mulliken_overlap_charges, key=lambda x: abs(x[2]))
            print(f"  Strongest overlap: atoms {strong[0]}-{strong[1]}, charge={strong[2]:.4f}")

        if result.mulliken_orbital_populations:
            print(f"Mulliken Orbital Populations: {len(result.mulliken_orbital_populations)} atoms")
            # Show first few atoms as example
            for pop in result.mulliken_orbital_populations[:3]:
                orb_str = f"s={pop.s_total:.2f}, p={pop.p_total:.2f}"
                if pop.d_total > 0.01:
                    orb_str += f", d={pop.d_total:.2f}"
                if pop.f_total > 0.001:
                    orb_str += f", f={pop.f_total:.3f}"
                print(f"  {pop.atom_index} {pop.element}: {orb_str}")
            if len(result.mulliken_orbital_populations) > 3:
                print(f"  ... and {len(result.mulliken_orbital_populations) - 3} more")

        if result.mulliken_orbital_charges:
            print(f"Mulliken Orbital Charges: {len(result.mulliken_orbital_charges)} MO entries")
            # Show first few entries
            for charge in result.mulliken_orbital_charges[:3]:
                print(f"  MO {charge.mo_index}: {charge.atom_index}{charge.element} {charge.orbital} = {charge.charge:.6f} ({charge.uncorrected_charge:.6f})")
            if len(result.mulliken_orbital_charges) > 3:
                print(f"  ... and {len(result.mulliken_orbital_charges) - 3} more")

        if result.loewdin_orbital_charges:
            print(f"Loewdin Orbital Charges: {len(result.loewdin_orbital_charges)} MO entries")
            # Show first few entries
            for charge in result.loewdin_orbital_charges[:3]:
                print(f"  MO {charge.mo_index}: {charge.atom_index}{charge.element} {charge.orbital} = {charge.charge:.6f}")
            if len(result.loewdin_orbital_charges) > 3:
                print(f"  ... and {len(result.loewdin_orbital_charges) - 3} more")

        if result.thermochemistry:
            print(f"Gibbs Free Energy: {result.thermochemistry.gibbs_free_energy:.6f} Eh")

        if result.nmr_data:
            print(f"NMR Shifts: {len(result.nmr_data.chemical_shifts)} nuclei")
            print(f"J-Couplings: {len(result.nmr_data.j_couplings)} pairs")

        if result.chemical_shielding_tensors:
            print(f"Chemical Shielding Tensors: {len(result.chemical_shielding_tensors)} nuclei")
            # Show first tensor as example
            if result.chemical_shielding_tensors:
                t = result.chemical_shielding_tensors[0]
                print(f"  Example: {t.atom_index}{t.element}, total_iso={t.total_iso:.2f} ppm")

        if result.geometric_perturbations:
            print(f"Geometric Perturbations: {result.geometric_perturbations.num_perturbations} perturbations")
            print(f"  Time: {result.geometric_perturbations.total_time:.1f}s, Memory: {result.geometric_perturbations.max_memory_mb:.0f} MB")

        if result.shark_integrals:
            print(f"SHARK Integrals: {result.shark_integrals.num_basis_functions} basis functions, {result.shark_integrals.num_shells} shells")
            print(f"  Shell pairs: {result.shark_integrals.shell_pairs_after_screening:,} (after screening)")

        if result.cosx_grids:
            print(f"COSX Grids: {len(result.cosx_grids)} grids")
            total_points = sum(g.total_grid_points for g in result.cosx_grids)
            print(f"  Total grid points: {total_points:,}")

        if result.raman_spectrum:
            strongest = max(result.raman_spectrum, key=lambda x: x.activity)
            print(f"Raman Spectrum: {len(result.raman_spectrum)} modes")
            print(f"  Strongest: {strongest.frequency:.1f} cm^-1 (activity: {strongest.activity:.1f})")

        if result.dispersion_correction:
            print(f"Dispersion: {result.dispersion_correction.method}")
            print(f"  Energy: {result.dispersion_correction.total_correction_kcal:.2f} kcal/mol")

        if result.normal_modes:
            print(f"Normal Modes: {len(result.normal_modes)} modes with displacement vectors")
            if result.normal_modes:
                print(f"  First mode: {result.normal_modes[0].frequency:.1f} cm^-1 ({len(result.normal_modes[0].displacements)} atoms)")

        if result.scf_iterations:
            print(f"SCF Iterations: {len(result.scf_iterations)} iterations")
            if result.scf_iterations:
                final = result.scf_iterations[-1]
                print(f"  Final: ΔE={final.delta_e:.2e}, RMSDP={final.rmsdp:.2e}")

        if result.timing_data:
            print(f"Total Time: {result.timing_data.total_time:.1f} sec")
            if result.timing_data.total_run_time > 0:
                mins = int(result.timing_data.total_run_time // 60)
                secs = result.timing_data.total_run_time % 60
                print(f"Total Run Time: {mins} min {secs:.1f} sec ({result.timing_data.total_run_time:.1f} sec)")
            if result.timing_data.fock_matrix_time > 0:
                print(f"  Fock Matrix: {result.timing_data.fock_matrix_time:.1f} sec ({result.timing_data.fock_matrix_time/result.timing_data.total_time*100:.1f}%)")

        if result.dft_grid_info:
            print(f"DFT Grid: {result.dft_grid_info.total_grid_points:,} points")
            print(f"  {result.dft_grid_info.angular_grid}, {result.dft_grid_info.avg_points_per_atom} pts/atom")

        if result.basis_set_info:
            print(f"Basis Set: {result.basis_set_info.name}")
            print(f"  {result.basis_set_info.num_basis_functions} functions, {result.basis_set_info.num_primitive_gaussians} primitives")
            if result.basis_set_info.element_basis:
                print(f"  Elements: {', '.join(result.basis_set_info.element_basis.keys())}")

        if result.energy_components:
            ec = result.energy_components
            print(f"Energy Components:")
            print(f"  Nuclear Repulsion: {ec.nuclear_repulsion:.6f} Eh")
            print(f"  Electronic Energy: {ec.electronic_energy:.6f} Eh")
            if ec.kinetic_energy != 0.0:
                print(f"  Kinetic Energy: {ec.kinetic_energy:.6f} Eh")
                print(f"  Potential Energy: {ec.potential_energy:.6f} Eh")
                print(f"  Virial Ratio: {ec.virial_ratio:.8f}")
            if ec.xc_energy != 0.0:
                print(f"  XC Energy: {ec.xc_energy:.6f} Eh (X: {ec.exchange_energy:.6f}, C: {ec.correlation_energy:.6f})")

        if result.cpcm_solvation:
            cpcm = result.cpcm_solvation
            print(f"CPCM Solvation:")
            print(f"  Surface Charge: {cpcm.surface_charge:.8f}")
            print(f"  Corrected Charge: {cpcm.corrected_charge:.8f}")
            print(f"  Outlying Charge Correction: {cpcm.outlying_charge_correction:.8f} Eh ({cpcm.outlying_charge_correction_ev:.5f} eV)")
            if cpcm.dielectric_energy != 0.0:
                print(f"  Dielectric Energy: {cpcm.dielectric_energy:.8f} Eh")

        if result.scf_convergence:
            conv = result.scf_convergence
            print(f"SCF Convergence:")
            print(f"  Energy Change: {conv.energy_change:.2e}")
            print(f"  MAX-Density Change: {conv.max_density_change:.2e}")
            print(f"  RMS-Density Change: {conv.rms_density_change:.2e}")
            print(f"  DIIS Error: {conv.diis_error:.2e}")
            print(f"  Orbital Gradient: {conv.orbital_gradient:.2e}")
            print(f"  Orbital Rotation: {conv.orbital_rotation:.2e}")
