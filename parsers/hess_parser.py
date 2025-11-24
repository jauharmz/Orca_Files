"""
ORCA Hessian File Parser (.hess)

Parses Hessian files containing:
- Vibrational frequencies
- Normal modes
- IR spectrum data
- Atomic structure
- Dipole derivatives
- Eigenvalues
"""

import re
from dataclasses import dataclass, field


@dataclass
class HessAtom:
    """Atom data from Hessian file."""
    element: str
    mass: float
    x: float  # Bohr
    y: float
    z: float


@dataclass
class IRMode:
    """Single IR vibrational mode."""
    frequency: float  # cm^-1
    intensity: float  # km/mol
    dipole_x: float = 0.0
    dipole_y: float = 0.0
    dipole_z: float = 0.0


@dataclass
class HessianData:
    """Parsed Hessian file data."""
    multiplicity: int = 1
    num_atoms: int = 0
    num_modes: int = 0
    temperature: float = 0.0
    frequency_scale_factor: float = 1.0
    atoms: list[HessAtom] = field(default_factory=list)
    frequencies: list[float] = field(default_factory=list)
    eigenvalues: list[float] = field(default_factory=list)
    ir_spectrum: list[IRMode] = field(default_factory=list)
    dipole_derivatives: list[list[float]] = field(default_factory=list)  # [mode][x,y,z]

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'multiplicity': self.multiplicity,
            'num_atoms': self.num_atoms,
            'num_modes': self.num_modes,
            'temperature': self.temperature,
            'frequency_scale_factor': self.frequency_scale_factor,
            'atoms': [
                {'element': a.element, 'mass': a.mass, 'x': a.x, 'y': a.y, 'z': a.z}
                for a in self.atoms
            ],
            'frequencies': self.frequencies,
            'eigenvalues': self.eigenvalues,
            'ir_spectrum': [
                {
                    'frequency': m.frequency,
                    'intensity': m.intensity,
                    'dipole_x': m.dipole_x,
                    'dipole_y': m.dipole_y,
                    'dipole_z': m.dipole_z
                }
                for m in self.ir_spectrum
            ],
            'dipole_derivatives': self.dipole_derivatives,
            'num_real_modes': len([f for f in self.frequencies if f > 0]),
            'frequency_range': {
                'min': min([f for f in self.frequencies if f > 0], default=0),
                'max': max(self.frequencies, default=0)
            }
        }


def parse_hess_file(filepath: str) -> HessianData:
    """Parse ORCA Hessian file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_hess_content(content)


def parse_hess_content(content: str) -> HessianData:
    """Parse Hessian content from string."""
    data = HessianData()

    # Parse multiplicity
    mult_match = re.search(r'\$multiplicity\s+(\d+)', content)
    if mult_match:
        data.multiplicity = int(mult_match.group(1))

    # Parse temperature
    temp_match = re.search(r'\$actual_temperature\s+(-?\d+\.?\d*)', content)
    if temp_match:
        data.temperature = float(temp_match.group(1))

    # Parse frequency scale factor
    scale_match = re.search(r'\$frequency_scale_factor\s+(-?\d+\.?\d*)', content)
    if scale_match:
        data.frequency_scale_factor = float(scale_match.group(1))

    # Parse atoms
    atoms_section = re.search(
        r'\$atoms\s+(\d+)\s+(.*?)(?=\$)',
        content, re.DOTALL
    )

    if atoms_section:
        data.num_atoms = int(atoms_section.group(1))
        atom_lines = atoms_section.group(2).strip().split('\n')

        for line in atom_lines:
            parts = line.split()
            if len(parts) >= 5:
                atom = HessAtom(
                    element=parts[0],
                    mass=float(parts[1]),
                    x=float(parts[2]),
                    y=float(parts[3]),
                    z=float(parts[4])
                )
                data.atoms.append(atom)

    # Parse vibrational frequencies
    freq_section = re.search(
        r'\$vibrational_frequencies\s+(\d+)\s+(.*?)(?=\$)',
        content, re.DOTALL
    )

    if freq_section:
        data.num_modes = int(freq_section.group(1))
        freq_lines = freq_section.group(2).strip().split('\n')

        for line in freq_lines:
            parts = line.split()
            if len(parts) >= 2:
                freq = float(parts[1])
                data.frequencies.append(freq)

    # Parse eigenvalues
    eigen_section = re.search(
        r'\$eigenvalues_hessian\s+\d+\s+(.*?)(?=\$)',
        content, re.DOTALL
    )

    if eigen_section:
        eigen_lines = eigen_section.group(1).strip().split('\n')
        for line in eigen_lines:
            parts = line.split()
            if len(parts) >= 2:
                data.eigenvalues.append(float(parts[1]))

    # Parse dipole derivatives
    dipole_section = re.search(
        r'\$dipole_derivatives\s+\d+\s+(.*?)(?=\$)',
        content, re.DOTALL
    )

    if dipole_section:
        dipole_lines = dipole_section.group(1).strip().split('\n')
        for line in dipole_lines:
            parts = line.split()
            if len(parts) >= 3:
                data.dipole_derivatives.append([
                    float(parts[0]),
                    float(parts[1]),
                    float(parts[2])
                ])

    # Parse IR spectrum
    ir_section = re.search(
        r'\$ir_spectrum\s+\d+\s+(.*?)(?=\$)',
        content, re.DOTALL
    )

    if ir_section:
        ir_lines = ir_section.group(1).strip().split('\n')

        for line in ir_lines:
            parts = line.split()
            if len(parts) >= 6:
                mode = IRMode(
                    frequency=float(parts[0]),
                    intensity=float(parts[2]),  # km/mol
                    dipole_x=float(parts[3]),
                    dipole_y=float(parts[4]),
                    dipole_z=float(parts[5])
                )
                data.ir_spectrum.append(mode)

    return data


def get_ir_spectrum_data(hess_data: HessianData) -> tuple[list[float], list[float]]:
    """Get IR spectrum as (frequencies, intensities) for plotting."""
    frequencies = []
    intensities = []

    for mode in hess_data.ir_spectrum:
        if mode.frequency > 0:
            frequencies.append(mode.frequency)
            intensities.append(mode.intensity)

    return frequencies, intensities


def find_strongest_peaks(hess_data: HessianData, n: int = 10) -> list[IRMode]:
    """Find the N strongest IR peaks."""
    real_modes = [m for m in hess_data.ir_spectrum if m.frequency > 0]
    sorted_modes = sorted(real_modes, key=lambda m: m.intensity, reverse=True)
    return sorted_modes[:n]


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        hess = parse_hess_file(sys.argv[1])
        print(f"Number of atoms: {hess.num_atoms}")
        print(f"Number of modes: {hess.num_modes}")
        print(f"Temperature: {hess.temperature} K")
        print(f"Scale factor: {hess.frequency_scale_factor}")

        if hess.atoms:
            print(f"\nAtoms: {len(hess.atoms)}")
            for i, a in enumerate(hess.atoms[:3]):
                print(f"  {i}: {a.element} ({a.mass:.2f} amu)")

        if hess.frequencies:
            real_freqs = [f for f in hess.frequencies if f > 0]
            print(f"\nReal frequencies: {len(real_freqs)}")
            print(f"Range: {min(real_freqs):.1f} - {max(real_freqs):.1f} cm^-1")

        if hess.eigenvalues:
            print(f"\nEigenvalues: {len(hess.eigenvalues)}")

        if hess.dipole_derivatives:
            print(f"Dipole derivatives: {len(hess.dipole_derivatives)} modes")

        if hess.ir_spectrum:
            print(f"\nTop 5 IR peaks:")
            peaks = find_strongest_peaks(hess, 5)
            for i, peak in enumerate(peaks):
                print(f"  {i+1}. {peak.frequency:.1f} cm^-1, {peak.intensity:.1f} km/mol")
