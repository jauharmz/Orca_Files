"""
ORCA Main Output File Parser (.out)

Parses the main ORCA output file to extract:
- Job information (method, basis set, charge, multiplicity)
- Final energy
- SCF energies and convergence
- Dipole moment
- Vibrational frequencies
- Thermochemistry
- Mulliken charges
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
class Thermochemistry:
    """Thermochemistry data."""
    temperature: float = 298.15
    pressure: float = 1.0
    total_mass: float = 0.0
    zero_point_energy: Optional[float] = None
    inner_energy: Optional[float] = None
    enthalpy: Optional[float] = None
    entropy: Optional[float] = None
    gibbs_free_energy: Optional[float] = None


@dataclass
class OrcaOutput:
    """Complete parsed ORCA output."""
    job_info: JobInfo
    final_energy: Optional[float] = None
    scf_energies: list[float] = field(default_factory=list)
    optimization_energies: list[float] = field(default_factory=list)
    dipole_moment: Optional[DipoleMoment] = None
    frequencies: list[float] = field(default_factory=list)
    mulliken_charges: dict[int, tuple[str, float]] = field(default_factory=dict)
    thermochemistry: Optional[Thermochemistry] = None

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
            'frequencies': self.frequencies,
            'mulliken_charges': {
                str(k): {'element': v[0], 'charge': v[1]}
                for k, v in self.mulliken_charges.items()
            },
            'thermochemistry': {
                'temperature': self.thermochemistry.temperature,
                'pressure': self.thermochemistry.pressure,
                'total_mass': self.thermochemistry.total_mass,
                'zero_point_energy': self.thermochemistry.zero_point_energy,
                'inner_energy': self.thermochemistry.inner_energy,
                'enthalpy': self.thermochemistry.enthalpy,
                'gibbs_free_energy': self.thermochemistry.gibbs_free_energy
            } if self.thermochemistry else None
        }


def parse_out_file(filepath: str) -> OrcaOutput:
    """
    Parse ORCA output file.

    Args:
        filepath: Path to .out file

    Returns:
        OrcaOutput object with parsed data
    """
    with open(filepath, 'r') as f:
        content = f.read()

    return parse_out_content(content)


def parse_out_content(content: str) -> OrcaOutput:
    """
    Parse ORCA output content from string.

    Args:
        content: Output file content as string

    Returns:
        OrcaOutput object with parsed data
    """
    job_info = parse_job_info(content)
    final_energy = parse_final_energy(content)
    scf_energies = parse_scf_energies(content)
    opt_energies = parse_optimization_energies(content)
    dipole = parse_dipole_moment(content)
    frequencies = parse_frequencies(content)
    mulliken = parse_mulliken_charges(content)
    thermo = parse_thermochemistry(content)

    return OrcaOutput(
        job_info=job_info,
        final_energy=final_energy,
        scf_energies=scf_energies,
        optimization_energies=opt_energies,
        dipole_moment=dipole,
        frequencies=frequencies,
        mulliken_charges=mulliken,
        thermochemistry=thermo
    )


def parse_job_info(content: str) -> JobInfo:
    """Extract job information from output."""
    info = JobInfo()

    # Extract method and job types from input line
    # Look for lines like "! B3LYP D3BJ 6-311++G(d,p) OPT FREQ"
    input_match = re.search(r'\|\s*>\s*!\s*(.+)', content)
    if input_match:
        keywords = input_match.group(1).split()
        # Common methods
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
    # Get the last final energy (after optimization completes)
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
    # Find the last dipole moment section
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


def parse_frequencies(content: str) -> list[float]:
    """Extract vibrational frequencies in cm^-1."""
    frequencies = []

    # Look for frequency section
    freq_section = re.search(
        r'VIBRATIONAL FREQUENCIES\s*-+\s*(.*?)(?=\s*-{5,})',
        content, re.DOTALL
    )

    if freq_section:
        # Extract frequency values (skip imaginary frequencies marked with *)
        freq_matches = re.findall(
            r'^\s*\d+:\s+(-?\d+\.?\d*)\s+cm\*\*-1',
            freq_section.group(1), re.MULTILINE
        )
        for match in freq_matches:
            freq = float(match)
            if freq > 0:  # Only positive frequencies
                frequencies.append(freq)

    return frequencies


def parse_mulliken_charges(content: str) -> dict[int, tuple[str, float]]:
    """Extract Mulliken atomic charges."""
    charges = {}

    # Find the last Mulliken charges section
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

    # Total enthalpy
    enthalpy_match = re.search(r'Total enthalpy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if enthalpy_match:
        thermo.enthalpy = float(enthalpy_match.group(1))

    # Final Gibbs free energy
    gibbs_match = re.search(r'Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.?\d*)\s+Eh', content)
    if gibbs_match:
        thermo.gibbs_free_energy = float(gibbs_match.group(1))

    # Only return if we found thermochemistry data
    if thermo.total_mass > 0:
        return thermo
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

        if result.frequencies:
            print(f"Frequencies: {len(result.frequencies)} modes")
            print(f"  Lowest: {min(result.frequencies):.1f} cm^-1")
            print(f"  Highest: {max(result.frequencies):.1f} cm^-1")

        if result.thermochemistry:
            print(f"Gibbs Free Energy: {result.thermochemistry.gibbs_free_energy:.6f} Eh")
