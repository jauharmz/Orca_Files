"""
XYZ File Parser for ORCA output

Parses standard XYZ molecular geometry files.
Format:
  Line 1: Number of atoms
  Line 2: Comment (often contains energy)
  Line 3+: Element X Y Z
"""

import re
import logging
from dataclasses import dataclass
from typing import Optional

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)


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
    logger.info(f"Parsing XYZ file: {filepath}")
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        logger.debug(f"Successfully read file {filepath}, size: {len(content)} bytes")
        result = parse_xyz_content(content)
        logger.info(f"Successfully parsed {result.num_atoms} atoms from {filepath}")
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


def parse_xyz_content(content: str) -> MoleculeGeometry:
    """
    Parse XYZ content from string.

    Args:
        content: XYZ file content as string

    Returns:
        MoleculeGeometry object with parsed data
    """
    try:
        lines = content.strip().split('\n')

        if len(lines) < 3:
            raise ValueError(f"Invalid XYZ format: expected at least 3 lines, got {len(lines)}")

        # Line 1: Number of atoms
        try:
            num_atoms = int(lines[0].strip())
            logger.debug(f"Parsing XYZ with {num_atoms} atoms")
        except ValueError as e:
            logger.error(f"Failed to parse number of atoms from line 1: '{lines[0]}': {e}")
            raise

        # Line 2: Comment (may contain energy)
        comment = lines[1].strip()
        energy = extract_energy(comment)

        # Lines 3+: Atom coordinates
        atoms = []
        expected_lines = 2 + num_atoms
        if len(lines) < expected_lines:
            logger.warning(f"Expected {expected_lines} lines for {num_atoms} atoms, got {len(lines)}")

        for i in range(2, min(2 + num_atoms, len(lines))):
            try:
                parts = lines[i].split()
                if len(parts) < 4:
                    logger.warning(f"Line {i+1} has insufficient data: '{lines[i]}' (expected: element x y z)")
                    continue
                atom = Atom(
                    symbol=parts[0],
                    x=float(parts[1]),
                    y=float(parts[2]),
                    z=float(parts[3])
                )
                atoms.append(atom)
            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse atom at line {i+1}: '{lines[i]}': {e}")

        if len(atoms) != num_atoms:
            logger.warning(f"Expected {num_atoms} atoms, but parsed {len(atoms)}")

        return MoleculeGeometry(
            num_atoms=num_atoms,
            comment=comment,
            atoms=atoms,
            energy=energy
        )
    except Exception as e:
        logger.error(f"Error parsing XYZ content: {e}", exc_info=True)
        raise


def parse_trajectory_xyz(filepath: str) -> list[MoleculeGeometry]:
    """
    Parse a trajectory XYZ file containing multiple frames.

    Args:
        filepath: Path to trajectory XYZ file

    Returns:
        List of MoleculeGeometry objects, one per frame
    """
    logger.info(f"Parsing trajectory XYZ file: {filepath}")
    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        logger.error(f"File not found: {filepath}")
        raise
    except IOError as e:
        logger.error(f"Error reading file {filepath}: {e}")
        raise

    frames = []
    lines = content.strip().split('\n')

    i = 0
    frame_num = 0
    while i < len(lines):
        try:
            # Read number of atoms
            if i >= len(lines):
                logger.warning(f"Reached end of file while parsing frame {frame_num}")
                break

            num_atoms = int(lines[i].strip())
            frame_num += 1
            logger.debug(f"Parsing frame {frame_num} with {num_atoms} atoms")

            # Extract this frame's content
            frame_end = i + num_atoms + 2
            if frame_end > len(lines):
                logger.warning(f"Frame {frame_num} incomplete: expected {num_atoms + 2} lines, "
                             f"only {len(lines) - i} remaining")
                break

            frame_lines = lines[i:frame_end]
            frame_content = '\n'.join(frame_lines)

            frame = parse_xyz_content(frame_content)
            frames.append(frame)

            i += num_atoms + 2
        except (ValueError, IndexError) as e:
            logger.error(f"Error parsing frame {frame_num} at line {i}: {e}")
            break

    logger.info(f"Successfully parsed {len(frames)} frames from {filepath}")
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
