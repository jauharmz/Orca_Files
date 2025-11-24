"""
ORCA CPCM Solvation File Parser (.cpcm)

Parses CPCM solvation model data including:
- Summary energies and volumes
- Atomic coordinates with radii
- Surface points with charges and potentials
"""

import re
from dataclasses import dataclass, field


@dataclass
class CPCMAtom:
    """Atom in CPCM cavity."""
    x: float
    y: float
    z: float
    radius: float
    atomic_number: int


@dataclass
class SurfacePoint:
    """Single CPCM surface point."""
    x: float
    y: float
    z: float
    area: float
    potential: float
    charge: float
    atom_index: int


@dataclass
class CPCMData:
    """Parsed CPCM file data."""
    num_atoms: int = 0
    num_surface_points: int = 0
    volume: float = 0.0
    area: float = 0.0
    dielectric_energy: float = 0.0
    one_electron_energy: float = 0.0
    epsilon: float = 0.0
    atoms: list[CPCMAtom] = field(default_factory=list)
    surface_points: list[SurfacePoint] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'num_atoms': self.num_atoms,
            'num_surface_points': self.num_surface_points,
            'volume': self.volume,
            'area': self.area,
            'dielectric_energy': self.dielectric_energy,
            'one_electron_energy': self.one_electron_energy,
            'dielectric_energy_kcal': self.dielectric_energy * 627.509,
            'atoms': [
                {'x': a.x, 'y': a.y, 'z': a.z, 'radius': a.radius, 'Z': a.atomic_number}
                for a in self.atoms
            ],
            'surface_points_count': len(self.surface_points),
            'charge_distribution': {
                'min': min([sp.charge for sp in self.surface_points], default=0),
                'max': max([sp.charge for sp in self.surface_points], default=0),
                'total': sum([sp.charge for sp in self.surface_points])
            } if self.surface_points else None
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

    # Parse header
    in_atoms = False
    in_surface = False
    atom_count = 0
    surface_count = 0

    for i, line in enumerate(lines):
        # Header parsing
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

        # Detect sections
        if line.strip() and not line.startswith('#'):
            parts = line.split()

            # After header, first data is atomic coordinates
            if len(parts) == 5 and atom_count < data.num_atoms:
                try:
                    atom = CPCMAtom(
                        x=float(parts[0]),
                        y=float(parts[1]),
                        z=float(parts[2]),
                        radius=float(parts[3]),
                        atomic_number=int(parts[4])
                    )
                    data.atoms.append(atom)
                    atom_count += 1
                except ValueError:
                    pass

            # Surface points have 10 columns
            elif len(parts) >= 10 and atom_count == data.num_atoms:
                try:
                    sp = SurfacePoint(
                        x=float(parts[0]),
                        y=float(parts[1]),
                        z=float(parts[2]),
                        area=float(parts[3]),
                        potential=float(parts[4]),
                        charge=float(parts[5]),
                        atom_index=int(parts[9])
                    )
                    data.surface_points.append(sp)
                except (ValueError, IndexError):
                    pass

    return data


def get_charge_by_atom(data: CPCMData) -> dict[int, float]:
    """Get total surface charge per atom."""
    charges = {}
    for sp in data.surface_points:
        if sp.atom_index not in charges:
            charges[sp.atom_index] = 0.0
        charges[sp.atom_index] += sp.charge
    return charges


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

        if data.atoms:
            print(f"\nAtoms loaded: {len(data.atoms)}")
            for i, a in enumerate(data.atoms[:3]):
                print(f"  {i}: Z={a.atomic_number}, r={a.radius:.3f}")

        if data.surface_points:
            print(f"\nSurface points: {len(data.surface_points)}")
            charges = [sp.charge for sp in data.surface_points]
            print(f"  Charge range: {min(charges):.6f} to {max(charges):.6f}")
            print(f"  Total charge: {sum(charges):.6f}")

            # Charge per atom
            atom_charges = get_charge_by_atom(data)
            print(f"\nCharge by atom (top 5):")
            sorted_charges = sorted(atom_charges.items(), key=lambda x: abs(x[1]), reverse=True)
            for idx, charge in sorted_charges[:5]:
                print(f"  Atom {idx}: {charge:.6f}")
