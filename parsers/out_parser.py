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
from dataclasses import dataclass, field
from typing import Optional


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
class NMRShift:
    """NMR chemical shift for a nucleus."""
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
class OrcaOutput:
    """Complete parsed ORCA output."""
    job_info: JobInfo
    final_energy: Optional[float] = None
    scf_energies: list[float] = field(default_factory=list)
    optimization_energies: list[float] = field(default_factory=list)
    dipole_moment: Optional[DipoleMoment] = None
    polarizability: Optional[Polarizability] = None
    orbital_energies: list[OrbitalEnergy] = field(default_factory=list)
    frequencies: list[float] = field(default_factory=list)
    ir_spectrum: list[IRMode] = field(default_factory=list)
    raman_spectrum: list[RamanMode] = field(default_factory=list)
    dispersion_correction: Optional[DispersionCorrection] = None
    mulliken_charges: dict[int, tuple[str, float]] = field(default_factory=dict)
    loewdin_charges: dict[int, tuple[str, float]] = field(default_factory=dict)
    mayer_bond_orders: list[tuple[int, int, float]] = field(default_factory=list)
    thermochemistry: Optional[Thermochemistry] = None
    nmr_data: Optional[NMRData] = None

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
            'dispersion_correction': {
                'method': self.dispersion_correction.method,
                'total_correction': self.dispersion_correction.total_correction,
                'total_correction_kcal': self.dispersion_correction.total_correction_kcal,
                'e6_kcal': self.dispersion_correction.e6_kcal,
                'e8_kcal': self.dispersion_correction.e8_kcal
            } if self.dispersion_correction else None,
            'mulliken_charges': {
                str(k): {'element': v[0], 'charge': v[1]}
                for k, v in self.mulliken_charges.items()
            },
            'loewdin_charges': {
                str(k): {'element': v[0], 'charge': v[1]}
                for k, v in self.loewdin_charges.items()
            },
            'mayer_bond_orders': [
                {'atom1': b[0], 'atom2': b[1], 'order': b[2]}
                for b in self.mayer_bond_orders
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
            } if self.nmr_data else None
        }


def parse_out_file(filepath: str) -> OrcaOutput:
    """Parse ORCA output file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_out_content(content)


def parse_out_content(content: str) -> OrcaOutput:
    """Parse ORCA output content from string."""
    job_info = parse_job_info(content)
    final_energy = parse_final_energy(content)
    scf_energies = parse_scf_energies(content)
    opt_energies = parse_optimization_energies(content)
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
    thermo = parse_thermochemistry(content)
    nmr = parse_nmr_data(content)

    return OrcaOutput(
        job_info=job_info,
        final_energy=final_energy,
        scf_energies=scf_energies,
        optimization_energies=opt_energies,
        dipole_moment=dipole,
        polarizability=polarizability,
        orbital_energies=orbitals,
        frequencies=frequencies,
        ir_spectrum=ir_spectrum,
        raman_spectrum=raman_spectrum,
        dispersion_correction=dispersion,
        mulliken_charges=mulliken,
        loewdin_charges=loewdin,
        mayer_bond_orders=bond_orders,
        thermochemistry=thermo,
        nmr_data=nmr
    )


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
    """Extract static polarizability tensor."""
    pol_section = re.search(
        r'THE POLARIZABILITY TENSOR.*?'
        r'The raw cartesian tensor \(atomic units\):\s*'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*'
        r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*'
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
        return Polarizability(
            tensor=tensor,
            isotropic=float(pol_section.group(10))
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
        r'IR SPECTRUM\s*-+\s*(.*?)(?=-{20,}|\Z)',
        content, re.DOTALL
    )

    if ir_section:
        # Match: Mode, freq (cm^-1), eps, Int (km/mol), T^2, TX, TY, TZ
        matches = re.findall(
            r'^\s*(\d+):\s+(-?\d+\.?\d*)\s+\d+\.?\d*\s+(\d+\.?\d*)\s+\d+\.?\d*\s*'
            r'\(\s*(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\)',
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


if __name__ == '__main__':
    import sys
    import json

    if len(sys.argv) > 1:
        result = parse_out_file(sys.argv[1])
        print(f"Method: {result.job_info.method}")
        print(f"Basis: {result.job_info.basis_set}")
        print(f"Charge: {result.job_info.charge}, Mult: {result.job_info.multiplicity}")
        print(f"Final Energy: {result.final_energy:.6f} Eh")

        if result.dipole_moment:
            print(f"Dipole: {result.dipole_moment.magnitude_debye:.3f} Debye")

        if result.polarizability:
            print(f"Isotropic Polarizability: {result.polarizability.isotropic:.2f} a.u.")

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
            print(f"Bond Orders: {len(result.mayer_bond_orders)} bonds")

        if result.thermochemistry:
            print(f"Gibbs Free Energy: {result.thermochemistry.gibbs_free_energy:.6f} Eh")

        if result.nmr_data:
            print(f"NMR Shifts: {len(result.nmr_data.chemical_shifts)} nuclei")
            print(f"J-Couplings: {len(result.nmr_data.j_couplings)} pairs")

        if result.raman_spectrum:
            strongest = max(result.raman_spectrum, key=lambda x: x.activity)
            print(f"Raman Spectrum: {len(result.raman_spectrum)} modes")
            print(f"  Strongest: {strongest.frequency:.1f} cm^-1 (activity: {strongest.activity:.1f})")

        if result.dispersion_correction:
            print(f"Dispersion: {result.dispersion_correction.method}")
            print(f"  Energy: {result.dispersion_correction.total_correction_kcal:.2f} kcal/mol")
