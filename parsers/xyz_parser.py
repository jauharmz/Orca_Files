"""
XYZ File Parser for ORCA output

Parses standard XYZ molecular geometry files.
Format:
  Line 1: Number of atoms
  Line 2: Comment (often contains energy)
  Line 3+: Element X Y Z
"""

import re
from dataclasses import dataclass
from typing import Optional


@dataclass
class Atom:
    """Represents a single atom with its coordinates."""
    symbol: str
    x: float
    y: float
    z: float


@dataclass
class MoleculeGeometry:
    """Represents a molecular geometry from XYZ file."""
    num_atoms: int
    comment: str
    atoms: list[Atom]
    energy: Optional[float] = None

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'num_atoms': self.num_atoms,
            'comment': self.comment,
            'energy': self.energy,
            'atoms': [
                {'symbol': a.symbol, 'x': a.x, 'y': a.y, 'z': a.z}
                for a in self.atoms
            ]
        }


def parse_xyz_file(filepath: str) -> MoleculeGeometry:
    """
    Parse a single XYZ file.

    Args:
        filepath: Path to the XYZ file

    Returns:
        MoleculeGeometry object with parsed data
    """
    with open(filepath, 'r') as f:
        content = f.read()

    return parse_xyz_content(content)


def parse_xyz_content(content: str) -> MoleculeGeometry:
    """
    Parse XYZ content from string.

    Args:
        content: XYZ file content as string

    Returns:
        MoleculeGeometry object with parsed data
    """
    lines = content.strip().split('\n')

    # Line 1: Number of atoms
    num_atoms = int(lines[0].strip())

    # Line 2: Comment (may contain energy)
    comment = lines[1].strip()
    energy = extract_energy(comment)

    # Lines 3+: Atom coordinates
    atoms = []
    for i in range(2, 2 + num_atoms):
        parts = lines[i].split()
        atom = Atom(
            symbol=parts[0],
            x=float(parts[1]),
            y=float(parts[2]),
            z=float(parts[3])
        )
        atoms.append(atom)

    return MoleculeGeometry(
        num_atoms=num_atoms,
        comment=comment,
        atoms=atoms,
        energy=energy
    )


def parse_trajectory_xyz(filepath: str) -> list[MoleculeGeometry]:
    """
    Parse a trajectory XYZ file containing multiple frames.

    Args:
        filepath: Path to trajectory XYZ file

    Returns:
        List of MoleculeGeometry objects, one per frame
    """
    with open(filepath, 'r') as f:
        content = f.read()

    frames = []
    lines = content.strip().split('\n')

    i = 0
    while i < len(lines):
        # Read number of atoms
        num_atoms = int(lines[i].strip())

        # Extract this frame's content
        frame_lines = lines[i:i + num_atoms + 2]
        frame_content = '\n'.join(frame_lines)

        frame = parse_xyz_content(frame_content)
        frames.append(frame)

        i += num_atoms + 2

    return frames


def extract_energy(comment: str) -> Optional[float]:
    """
    Extract energy value from comment line.

    ORCA format: "Coordinates from ORCA-job p1xs0 E -662.998375154159"

    Args:
        comment: Comment line from XYZ file

    Returns:
        Energy in Hartree or None if not found
    """
    # Pattern for ORCA energy format
    match = re.search(r'E\s+(-?\d+\.?\d*)', comment)
    if match:
        return float(match.group(1))
    return None


if __name__ == '__main__':
    # Quick test
    import sys
    if len(sys.argv) > 1:
        geom = parse_xyz_file(sys.argv[1])
        print(f"Atoms: {geom.num_atoms}")
        print(f"Energy: {geom.energy} Hartree")
        for atom in geom.atoms[:5]:
            print(f"  {atom.symbol}: ({atom.x:.4f}, {atom.y:.4f}, {atom.z:.4f})")
