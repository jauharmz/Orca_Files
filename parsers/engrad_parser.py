"""
ORCA Energy and Gradient File Parser (.engrad)

Parses engrad files containing energy and gradient information.
"""

import re
from dataclasses import dataclass, field


@dataclass
class AtomGradient:
    """Atom with gradient information."""
    atomic_number: int
    x: float  # Bohr
    y: float
    z: float
    grad_x: float  # Eh/bohr
    grad_y: float
    grad_z: float


@dataclass
class EngradData:
    """Parsed engrad file data."""
    num_atoms: int = 0
    energy: float = 0.0
    gradients: list[float] = field(default_factory=list)
    atoms: list[AtomGradient] = field(default_factory=list)

    @property
    def max_gradient(self) -> float:
        """Maximum absolute gradient component."""
        if not self.gradients:
            return 0.0
        return max(abs(g) for g in self.gradients)

    @property
    def rms_gradient(self) -> float:
        """RMS gradient."""
        if not self.gradients:
            return 0.0
        return (sum(g**2 for g in self.gradients) / len(self.gradients)) ** 0.5

    def to_dict(self) -> dict:
        return {
            'num_atoms': self.num_atoms,
            'energy': self.energy,
            'max_gradient': self.max_gradient,
            'rms_gradient': self.rms_gradient,
            'gradients': self.gradients,
            'atoms': [
                {
                    'atomic_number': a.atomic_number,
                    'x': a.x, 'y': a.y, 'z': a.z,
                    'grad_x': a.grad_x, 'grad_y': a.grad_y, 'grad_z': a.grad_z
                }
                for a in self.atoms
            ]
        }


def parse_engrad_file(filepath: str) -> EngradData:
    """Parse ORCA engrad file."""
    with open(filepath, 'r') as f:
        content = f.read()
    return parse_engrad_content(content)


def parse_engrad_content(content: str) -> EngradData:
    """Parse engrad content from string."""
    data = EngradData()
    lines = content.strip().split('\n')

    # Parse sections
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if 'Number of atoms' in line:
            i += 1
            while i < len(lines) and lines[i].strip().startswith('#'):
                i += 1
            data.num_atoms = int(lines[i].strip())

        elif 'current total energy' in line:
            i += 1
            while i < len(lines) and lines[i].strip().startswith('#'):
                i += 1
            data.energy = float(lines[i].strip())

        elif 'current gradient' in line:
            i += 1
            while i < len(lines) and lines[i].strip().startswith('#'):
                i += 1
            for _ in range(data.num_atoms * 3):
                data.gradients.append(float(lines[i].strip()))
                i += 1
            i -= 1  # Adjust for loop increment

        elif 'atomic numbers and current coordinates' in line:
            i += 1
            while i < len(lines) and lines[i].strip().startswith('#'):
                i += 1
            grad_idx = 0
            for _ in range(data.num_atoms):
                parts = lines[i].split()
                atom = AtomGradient(
                    atomic_number=int(parts[0]),
                    x=float(parts[1]),
                    y=float(parts[2]),
                    z=float(parts[3]),
                    grad_x=data.gradients[grad_idx] if grad_idx < len(data.gradients) else 0,
                    grad_y=data.gradients[grad_idx+1] if grad_idx+1 < len(data.gradients) else 0,
                    grad_z=data.gradients[grad_idx+2] if grad_idx+2 < len(data.gradients) else 0
                )
                data.atoms.append(atom)
                grad_idx += 3
                i += 1
            i -= 1

        i += 1

    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_engrad_file(sys.argv[1])
        print(f"Atoms: {data.num_atoms}")
        print(f"Energy: {data.energy:.6f} Eh")
        print(f"Max gradient: {data.max_gradient:.6e} Eh/bohr")
        print(f"RMS gradient: {data.rms_gradient:.6e} Eh/bohr")
