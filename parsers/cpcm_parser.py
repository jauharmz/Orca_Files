"""
ORCA CPCM Solvation File Parser (.cpcm)

Parses CPCM solvation model data.
"""

import re
from dataclasses import dataclass, field


@dataclass
class CPCMData:
    """Parsed CPCM file data."""
    num_atoms: int = 0
    num_surface_points: int = 0
    volume: float = 0.0
    area: float = 0.0
    dielectric_energy: float = 0.0
    one_electron_energy: float = 0.0

    def to_dict(self) -> dict:
        return {
            'num_atoms': self.num_atoms,
            'num_surface_points': self.num_surface_points,
            'volume': self.volume,
            'area': self.area,
            'dielectric_energy': self.dielectric_energy,
            'one_electron_energy': self.one_electron_energy,
            'dielectric_energy_kcal': self.dielectric_energy * 627.509  # Eh to kcal/mol
        }


def parse_cpcm_file(filepath: str) -> CPCMData:
    """Parse ORCA CPCM file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_cpcm_content(content)


def parse_cpcm_content(content: str) -> CPCMData:
    """Parse CPCM content from string."""
    data = CPCMData()
    lines = content.split('\n')

    for line in lines:
        if '# Number of atoms' in line:
            data.num_atoms = int(line.split()[0])
        elif '# Number of surface points' in line:
            data.num_surface_points = int(line.split()[0])
        elif '# Volume' in line:
            data.volume = float(line.split()[0])
        elif '# Area' in line:
            data.area = float(line.split()[0])
        elif '# CPCM dielectric energy' in line:
            data.dielectric_energy = float(line.split()[0])
        elif '# One-electron operator energy' in line:
            data.one_electron_energy = float(line.split()[0])

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_cpcm_file(sys.argv[1])
        print(f"Atoms: {data.num_atoms}")
        print(f"Surface points: {data.num_surface_points}")
        print(f"Volume: {data.volume:.2f} bohr³")
        print(f"Area: {data.area:.2f} bohr²")
        print(f"Dielectric energy: {data.dielectric_energy:.6f} Eh")
        print(f"                 : {data.dielectric_energy * 627.509:.2f} kcal/mol")
