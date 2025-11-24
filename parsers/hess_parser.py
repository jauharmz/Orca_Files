"""
ORCA Hessian File Parser (.hess)

Parses Hessian files containing:
- Vibrational frequencies
- Normal modes
- IR spectrum data
"""

import re
from dataclasses import dataclass, field


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
    num_modes: int = 0
    frequencies: list[float] = field(default_factory=list)
    ir_spectrum: list[IRMode] = field(default_factory=list)
    # Normal modes could be added but they're large matrices

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'multiplicity': self.multiplicity,
            'num_modes': self.num_modes,
            'frequencies': self.frequencies,
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
            'num_real_modes': len([f for f in self.frequencies if f > 0]),
            'frequency_range': {
                'min': min([f for f in self.frequencies if f > 0], default=0),
                'max': max(self.frequencies, default=0)
            }
        }


def parse_hess_file(filepath: str) -> HessianData:
    """
    Parse ORCA Hessian file.

    Args:
        filepath: Path to .hess file

    Returns:
        HessianData object with parsed data
    """
    with open(filepath, 'r') as f:
        content = f.read()

    return parse_hess_content(content)


def parse_hess_content(content: str) -> HessianData:
    """
    Parse Hessian content from string.

    Args:
        content: Hessian file content as string

    Returns:
        HessianData object with parsed data
    """
    data = HessianData()

    # Parse multiplicity
    mult_match = re.search(r'\$multiplicity\s+(\d+)', content)
    if mult_match:
        data.multiplicity = int(mult_match.group(1))

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
    """
    Get IR spectrum as (frequencies, intensities) for plotting.

    Args:
        hess_data: HessianData object

    Returns:
        Tuple of (frequencies, intensities) lists, excluding zero modes
    """
    frequencies = []
    intensities = []

    for mode in hess_data.ir_spectrum:
        if mode.frequency > 0:
            frequencies.append(mode.frequency)
            intensities.append(mode.intensity)

    return frequencies, intensities


def find_strongest_peaks(hess_data: HessianData, n: int = 10) -> list[IRMode]:
    """
    Find the N strongest IR peaks.

    Args:
        hess_data: HessianData object
        n: Number of peaks to return

    Returns:
        List of IRMode objects sorted by intensity
    """
    real_modes = [m for m in hess_data.ir_spectrum if m.frequency > 0]
    sorted_modes = sorted(real_modes, key=lambda m: m.intensity, reverse=True)
    return sorted_modes[:n]


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        hess = parse_hess_file(sys.argv[1])
        print(f"Number of modes: {hess.num_modes}")
        print(f"Real frequencies: {len([f for f in hess.frequencies if f > 0])}")

        if hess.frequencies:
            real_freqs = [f for f in hess.frequencies if f > 0]
            print(f"Frequency range: {min(real_freqs):.1f} - {max(real_freqs):.1f} cm^-1")

        if hess.ir_spectrum:
            print(f"\nTop 5 IR peaks:")
            peaks = find_strongest_peaks(hess, 5)
            for i, peak in enumerate(peaks):
                print(f"  {i+1}. {peak.frequency:.1f} cm^-1, {peak.intensity:.1f} km/mol")
