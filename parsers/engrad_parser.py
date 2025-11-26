"""
ORCA Energy and Gradient File Parser (.engrad)

Parses engrad files containing energy and gradient information.
"""

import re
import logging
from dataclasses import dataclass, field

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)


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
    logger.info(f"Parsing engrad file: {filepath}")
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        logger.debug(f"Successfully read file {filepath}, size: {len(content)} bytes")
        result = parse_engrad_content(content)
        logger.info(f"Successfully parsed engrad file with {result.num_atoms} atoms, "
                   f"energy: {result.energy:.6f} Eh")
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


def parse_engrad_content(content: str) -> EngradData:
    """Parse engrad content from string."""
    data = EngradData()
    lines = content.strip().split('\n')
    logger.debug(f"Parsing engrad content with {len(lines)} lines")

    # Parse sections
    i = 0
    try:
        while i < len(lines):
            line = lines[i].strip()

            if 'Number of atoms' in line:
                i += 1
                while i < len(lines) and lines[i].strip().startswith('#'):
                    i += 1
                if i >= len(lines):
                    logger.error("Unexpected end of file while parsing number of atoms")
                    break
                try:
                    data.num_atoms = int(lines[i].strip())
                    logger.debug(f"Found {data.num_atoms} atoms")
                except ValueError as e:
                    logger.error(f"Failed to parse number of atoms from line {i}: '{lines[i]}': {e}")

            elif 'current total energy' in line:
                i += 1
                while i < len(lines) and lines[i].strip().startswith('#'):
                    i += 1
                if i >= len(lines):
                    logger.error("Unexpected end of file while parsing energy")
                    break
                try:
                    data.energy = float(lines[i].strip())
                    logger.debug(f"Parsed energy: {data.energy:.6f} Eh")
                except ValueError as e:
                    logger.error(f"Failed to parse energy from line {i}: '{lines[i]}': {e}")

            elif 'current gradient' in line:
                i += 1
                while i < len(lines) and lines[i].strip().startswith('#'):
                    i += 1
                if i >= len(lines):
                    logger.error("Unexpected end of file while parsing gradients")
                    break
                expected_gradients = data.num_atoms * 3
                logger.debug(f"Parsing {expected_gradients} gradient components")
                for grad_num in range(expected_gradients):
                    if i >= len(lines):
                        logger.warning(f"Ran out of lines while parsing gradients "
                                     f"(expected {expected_gradients}, got {grad_num})")
                        break
                    try:
                        data.gradients.append(float(lines[i].strip()))
                    except ValueError as e:
                        logger.warning(f"Failed to parse gradient at line {i}: '{lines[i]}': {e}")
                        data.gradients.append(0.0)
                    i += 1
                i -= 1  # Adjust for loop increment

            elif 'atomic numbers and current coordinates' in line:
                i += 1
                while i < len(lines) and lines[i].strip().startswith('#'):
                    i += 1
                if i >= len(lines):
                    logger.error("Unexpected end of file while parsing coordinates")
                    break
                grad_idx = 0
                logger.debug(f"Parsing coordinates for {data.num_atoms} atoms")
                for atom_num in range(data.num_atoms):
                    if i >= len(lines):
                        logger.warning(f"Ran out of lines while parsing atoms "
                                     f"(expected {data.num_atoms}, got {atom_num})")
                        break
                    try:
                        parts = lines[i].split()
                        if len(parts) < 4:
                            logger.warning(f"Insufficient data at line {i}: '{lines[i]}' "
                                         f"(expected: atomic_number x y z)")
                            i += 1
                            continue

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
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Failed to parse atom at line {i}: '{lines[i]}': {e}")
                    i += 1
                i -= 1

            i += 1

    except Exception as e:
        logger.error(f"Error parsing engrad content at line {i}: {e}", exc_info=True)
        raise

    logger.debug(f"Parsed {len(data.atoms)} atoms with {len(data.gradients)} gradient components")
    return data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        data = parse_engrad_file(sys.argv[1])
        print(f"Atoms: {data.num_atoms}")
        print(f"Energy: {data.energy:.6f} Eh")
        print(f"Max gradient: {data.max_gradient:.6e} Eh/bohr")
        print(f"RMS gradient: {data.rms_gradient:.6e} Eh/bohr")
