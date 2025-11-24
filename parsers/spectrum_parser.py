"""
ORCA Spectrum File Parser (.spectrum)

Parses spectrum files containing vibronic spectrum data.
Format: Tab-separated with columns for Energy and various intensities.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class SpectrumData:
    """Parsed spectrum data."""
    energy: list[float]  # Energy in cm^-1
    total_spectrum: list[float]
    intensity_fc: Optional[list[float]] = None  # Franck-Condon
    intensity_ht: Optional[list[float]] = None  # Herzberg-Teller
    columns: list[str] = None

    @property
    def energy_ev(self) -> list[float]:
        """Convert energy from cm^-1 to eV."""
        # 1 eV = 8065.54 cm^-1
        return [e / 8065.54 for e in self.energy]

    @property
    def energy_nm(self) -> list[float]:
        """Convert energy from cm^-1 to wavelength in nm."""
        # wavelength (nm) = 1e7 / wavenumber (cm^-1)
        return [1e7 / e if e > 0 else 0 for e in self.energy]

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'energy_cm': self.energy,
            'energy_ev': self.energy_ev,
            'energy_nm': self.energy_nm,
            'total_spectrum': self.total_spectrum,
            'intensity_fc': self.intensity_fc,
            'intensity_ht': self.intensity_ht,
            'columns': self.columns,
            'num_points': len(self.energy),
            'energy_range': {
                'min_cm': min(self.energy),
                'max_cm': max(self.energy),
                'min_nm': min(self.energy_nm) if self.energy_nm else 0,
                'max_nm': max(self.energy_nm) if self.energy_nm else 0
            }
        }


def parse_spectrum_file(filepath: str) -> SpectrumData:
    """
    Parse ORCA spectrum file.

    Args:
        filepath: Path to .spectrum file

    Returns:
        SpectrumData object with parsed data
    """
    with open(filepath, 'r') as f:
        content = f.read()

    return parse_spectrum_content(content)


def parse_spectrum_content(content: str) -> SpectrumData:
    """
    Parse spectrum content from string.

    Args:
        content: Spectrum file content as string

    Returns:
        SpectrumData object with parsed data
    """
    lines = content.strip().split('\n')

    # First line is header
    header = lines[0].split()
    columns = header

    # Parse data
    energy = []
    total_spectrum = []
    intensity_fc = []
    intensity_ht = []

    for line in lines[1:]:
        if not line.strip():
            continue

        parts = line.split()
        if len(parts) >= 2:
            energy.append(float(parts[0]))
            total_spectrum.append(float(parts[1]))

            if len(parts) >= 3:
                intensity_fc.append(float(parts[2]))
            if len(parts) >= 4:
                intensity_ht.append(float(parts[3]))

    return SpectrumData(
        energy=energy,
        total_spectrum=total_spectrum,
        intensity_fc=intensity_fc if intensity_fc else None,
        intensity_ht=intensity_ht if intensity_ht else None,
        columns=columns
    )


def find_peaks(spectrum: SpectrumData, threshold: float = 0.1) -> list[dict]:
    """
    Find peaks in spectrum.

    Args:
        spectrum: SpectrumData object
        threshold: Minimum relative intensity for peak detection

    Returns:
        List of peak dictionaries with energy and intensity
    """
    intensities = spectrum.total_spectrum
    energies = spectrum.energy

    if not intensities:
        return []

    max_intensity = max(intensities)
    min_threshold = threshold * max_intensity

    peaks = []
    for i in range(1, len(intensities) - 1):
        if (intensities[i] > intensities[i-1] and
            intensities[i] > intensities[i+1] and
            intensities[i] > min_threshold):
            peaks.append({
                'energy_cm': energies[i],
                'energy_nm': 1e7 / energies[i] if energies[i] > 0 else 0,
                'intensity': intensities[i]
            })

    # Sort by intensity
    peaks.sort(key=lambda x: x['intensity'], reverse=True)

    return peaks


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        spectrum = parse_spectrum_file(sys.argv[1])
        print(f"Data points: {len(spectrum.energy)}")
        print(f"Energy range: {min(spectrum.energy):.1f} - {max(spectrum.energy):.1f} cm^-1")
        print(f"Wavelength range: {min(spectrum.energy_nm):.1f} - {max(spectrum.energy_nm):.1f} nm")

        peaks = find_peaks(spectrum)
        if peaks:
            print(f"\nTop 5 peaks:")
            for i, peak in enumerate(peaks[:5]):
                print(f"  {i+1}. {peak['energy_nm']:.1f} nm ({peak['energy_cm']:.1f} cm^-1)")
