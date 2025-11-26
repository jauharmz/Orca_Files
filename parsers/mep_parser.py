"""
Molecular Electrostatic Potential (MEP) Parser

Parses MEP data from ORCA outputs including:
- .cube files (Gaussian cube format)
- .property.txt files
- Integration with orca_2json for .gbw parsing

MEP is critical for:
- Understanding charge distribution
- Identifying interaction sites
- Predicting reactivity
- Degradation pathway analysis
"""

import re
import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, List
import struct

logger = logging.getLogger(__name__)


@dataclass
class MEPGridPoint:
    """Single point on MEP grid."""
    x: float  # Å
    y: float  # Å
    z: float  # Å
    potential: float  # a.u. (Hartree/electron)


@dataclass
class MEPData:
    """Molecular Electrostatic Potential data."""
    # Grid properties
    grid_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0)  # Å
    grid_dimensions: Tuple[int, int, int] = (0, 0, 0)  # nx, ny, nz
    grid_spacing: Tuple[float, float, float] = (0.0, 0.0, 0.0)  # dx, dy, dz in Å

    # Atomic positions (for reference)
    atoms: List[Tuple[str, float, float, float]] = field(default_factory=list)  # (element, x, y, z)

    # MEP values on 3D grid
    potential_grid: Optional[np.ndarray] = None  # 3D array (nx, ny, nz)

    # Statistical properties
    min_potential: Optional[float] = None  # a.u.
    max_potential: Optional[float] = None  # a.u.
    mean_potential: Optional[float] = None  # a.u.

    # Critical points (minima = nucleophilic sites, maxima = electrophilic sites)
    local_minima: List[Tuple[float, float, float, float]] = field(default_factory=list)  # (x, y, z, V)
    local_maxima: List[Tuple[float, float, float, float]] = field(default_factory=list)  # (x, y, z, V)

    # Surface MEP (van der Waals surface)
    vdw_surface_mep: Optional[List[Tuple[float, float, float, float]]] = None  # (x, y, z, V)

    def calculate_statistics(self):
        """Calculate statistical properties of the MEP grid."""
        if self.potential_grid is not None:
            self.min_potential = float(np.min(self.potential_grid))
            self.max_potential = float(np.max(self.potential_grid))
            self.mean_potential = float(np.mean(self.potential_grid))

    def find_critical_points(self, threshold: float = 0.01):
        """
        Find local minima and maxima in MEP grid.

        Args:
            threshold: Minimum potential difference to consider as critical point (a.u.)
        """
        if self.potential_grid is None:
            logger.warning("No potential grid available for critical point analysis")
            return

        nx, ny, nz = self.grid_dimensions
        ox, oy, oz = self.grid_origin
        dx, dy, dz = self.grid_spacing

        # Simple gradient-based approach for finding extrema
        # For production use, consider more sophisticated methods
        grid = self.potential_grid

        self.local_minima = []
        self.local_maxima = []

        # Check interior points (avoiding edges)
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                for k in range(1, nz - 1):
                    val = grid[i, j, k]

                    # Check if local minimum
                    is_min = all([
                        val < grid[i-1, j, k],
                        val < grid[i+1, j, k],
                        val < grid[i, j-1, k],
                        val < grid[i, j+1, k],
                        val < grid[i, j, k-1],
                        val < grid[i, j, k+1]
                    ])

                    # Check if local maximum
                    is_max = all([
                        val > grid[i-1, j, k],
                        val > grid[i+1, j, k],
                        val > grid[i, j-1, k],
                        val > grid[i, j+1, k],
                        val > grid[i, j, k-1],
                        val > grid[i, j, k+1]
                    ])

                    # Convert grid indices to Cartesian coordinates
                    x = ox + i * dx
                    y = oy + j * dy
                    z = oz + k * dz

                    if is_min and abs(val - self.min_potential) < threshold:
                        self.local_minima.append((x, y, z, val))
                    elif is_max and abs(val - self.max_potential) < threshold:
                        self.local_maxima.append((x, y, z, val))


class MEPParser:
    """Parser for Molecular Electrostatic Potential data."""

    @staticmethod
    def parse_cube_file(cube_file: str) -> MEPData:
        """
        Parse ORCA .cube file containing MEP data.

        Cube file format (Gaussian):
        - Line 1-2: Comments
        - Line 3: Number of atoms, origin (x, y, z) in Bohr
        - Lines 4-6: Grid points in each direction and spacing
        - Following lines: Atomic data (atomic number, charge, x, y, z)
        - Remaining lines: MEP values on grid

        Args:
            cube_file: Path to .cube file

        Returns:
            MEPData object
        """
        logger.info(f"Parsing cube file: {cube_file}")

        mep_data = MEPData()

        try:
            with open(cube_file, 'r') as f:
                lines = f.readlines()

            # Skip comment lines (first 2)
            line_idx = 2

            # Parse atom count and origin
            parts = lines[line_idx].split()
            n_atoms = int(parts[0])
            origin_bohr = [float(parts[i]) for i in range(1, 4)]
            # Convert Bohr to Angstrom
            BOHR_TO_ANGSTROM = 0.529177
            mep_data.grid_origin = tuple(o * BOHR_TO_ANGSTROM for o in origin_bohr)
            line_idx += 1

            # Parse grid dimensions and spacing
            grid_dims = []
            grid_spacing = []
            for _ in range(3):
                parts = lines[line_idx].split()
                n_points = int(parts[0])
                spacing_bohr = float(parts[1])  # Spacing in that direction
                grid_dims.append(abs(n_points))
                grid_spacing.append(abs(spacing_bohr) * BOHR_TO_ANGSTROM)
                line_idx += 1

            mep_data.grid_dimensions = tuple(grid_dims)
            mep_data.grid_spacing = tuple(grid_spacing)

            # Parse atomic positions
            atoms = []
            for _ in range(n_atoms):
                parts = lines[line_idx].split()
                atomic_num = int(parts[0])
                # Skip charge (parts[1])
                x, y, z = [float(parts[i]) * BOHR_TO_ANGSTROM for i in range(2, 5)]

                # Convert atomic number to element symbol (simplified)
                element_map = {
                    1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                    15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'
                }
                element = element_map.get(atomic_num, f"X{atomic_num}")
                atoms.append((element, x, y, z))
                line_idx += 1

            mep_data.atoms = atoms

            # Parse MEP grid data
            nx, ny, nz = mep_data.grid_dimensions
            potential_grid = np.zeros((nx, ny, nz))

            # Read all remaining data values
            data_values = []
            for line in lines[line_idx:]:
                data_values.extend([float(val) for val in line.split()])

            # Reshape into 3D grid (order: z fastest, then y, then x)
            if len(data_values) != nx * ny * nz:
                logger.warning(
                    f"Data size mismatch: expected {nx*ny*nz}, got {len(data_values)}"
                )

            idx = 0
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if idx < len(data_values):
                            potential_grid[i, j, k] = data_values[idx]
                            idx += 1

            mep_data.potential_grid = potential_grid

            # Calculate statistics
            mep_data.calculate_statistics()

            logger.info(f"Parsed MEP grid: {nx}×{ny}×{nz} points")
            logger.info(f"MEP range: {mep_data.min_potential:.6f} to {mep_data.max_potential:.6f} a.u.")

            return mep_data

        except Exception as e:
            logger.error(f"Error parsing cube file: {e}")
            raise

    @staticmethod
    def parse_property_file(property_file: str) -> Optional[dict]:
        """
        Parse MEP data from ORCA .property.txt file.

        The property.txt file may contain:
        - ESP charges (electrostatic potential derived charges)
        - Statistical data about MEP

        Args:
            property_file: Path to .property.txt file

        Returns:
            Dictionary with available MEP-related data
        """
        logger.info(f"Parsing property file for MEP data: {property_file}")

        mep_info = {
            'esp_charges': {},
            'has_mep_data': False
        }

        try:
            with open(property_file, 'r') as f:
                content = f.read()

            # Look for ESP charges section
            esp_pattern = r'\$esp_charges\s*\n(.*?)\$end'
            esp_match = re.search(esp_pattern, content, re.DOTALL)

            if esp_match:
                esp_section = esp_match.group(1)
                for line in esp_section.strip().split('\n'):
                    parts = line.split()
                    if len(parts) >= 3:
                        atom_idx = int(parts[0])
                        element = parts[1]
                        charge = float(parts[2])
                        mep_info['esp_charges'][atom_idx] = (element, charge)

            # Check for MEP-related keywords
            if 'ESP' in content or 'electrostatic potential' in content.lower():
                mep_info['has_mep_data'] = True

            return mep_info

        except Exception as e:
            logger.error(f"Error parsing property file: {e}")
            return None

    @staticmethod
    def calculate_mep_at_point(
        mep_data: MEPData,
        x: float,
        y: float,
        z: float
    ) -> Optional[float]:
        """
        Calculate MEP at a specific point using trilinear interpolation.

        Args:
            mep_data: MEPData object
            x, y, z: Coordinates in Angstrom

        Returns:
            Interpolated potential value in a.u., or None if out of bounds
        """
        if mep_data.potential_grid is None:
            return None

        ox, oy, oz = mep_data.grid_origin
        dx, dy, dz = mep_data.grid_spacing
        nx, ny, nz = mep_data.grid_dimensions

        # Convert to grid coordinates
        i = (x - ox) / dx
        j = (y - oy) / dy
        k = (z - oz) / dz

        # Check bounds
        if i < 0 or i >= nx - 1 or j < 0 or j >= ny - 1 or k < 0 or k >= nz - 1:
            return None

        # Trilinear interpolation
        i0, j0, k0 = int(i), int(j), int(k)
        i1, j1, k1 = i0 + 1, j0 + 1, k0 + 1

        # Fractional parts
        fi = i - i0
        fj = j - j0
        fk = k - k0

        grid = mep_data.potential_grid

        # Interpolate
        c000 = grid[i0, j0, k0]
        c100 = grid[i1, j0, k0]
        c010 = grid[i0, j1, k0]
        c110 = grid[i1, j1, k0]
        c001 = grid[i0, j0, k1]
        c101 = grid[i1, j0, k1]
        c011 = grid[i0, j1, k1]
        c111 = grid[i1, j1, k1]

        c00 = c000 * (1 - fi) + c100 * fi
        c01 = c001 * (1 - fi) + c101 * fi
        c10 = c010 * (1 - fi) + c110 * fi
        c11 = c011 * (1 - fi) + c111 * fi

        c0 = c00 * (1 - fj) + c10 * fj
        c1 = c01 * (1 - fj) + c11 * fj

        value = c0 * (1 - fk) + c1 * fk

        return value

    @staticmethod
    def extract_vdw_surface_mep(
        mep_data: MEPData,
        vdw_radii: Optional[dict] = None
    ) -> List[Tuple[float, float, float, float]]:
        """
        Extract MEP values on van der Waals surface.

        This is useful for visualizing reactive sites.

        Args:
            mep_data: MEPData object
            vdw_radii: Custom van der Waals radii {element: radius_in_angstrom}

        Returns:
            List of (x, y, z, potential) on vdW surface
        """
        # Default van der Waals radii (Angstrom)
        default_vdw = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
            'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98
        }

        vdw = vdw_radii if vdw_radii else default_vdw

        surface_points = []

        # For each atom, sample points on its vdW surface
        for element, ax, ay, az in mep_data.atoms:
            radius = vdw.get(element, 1.70)  # Default to carbon radius

            # Sample sphere surface using Fibonacci lattice
            n_points = 50  # Points per atom
            phi = (1 + np.sqrt(5)) / 2  # Golden ratio

            for i in range(n_points):
                theta = 2 * np.pi * i / phi
                cos_phi = 1 - (2 * i + 1) / n_points
                sin_phi = np.sqrt(1 - cos_phi**2)

                # Point on unit sphere
                x = radius * sin_phi * np.cos(theta) + ax
                y = radius * sin_phi * np.sin(theta) + ay
                z = radius * cos_phi + az

                # Get MEP at this point
                potential = MEPParser.calculate_mep_at_point(mep_data, x, y, z)
                if potential is not None:
                    surface_points.append((x, y, z, potential))

        return surface_points

    @staticmethod
    def print_summary(mep_data: MEPData) -> str:
        """Generate text summary of MEP analysis."""
        lines = []
        lines.append("=" * 80)
        lines.append("MOLECULAR ELECTROSTATIC POTENTIAL (MEP) ANALYSIS")
        lines.append("=" * 80)
        lines.append("")

        lines.append("GRID PROPERTIES:")
        lines.append("-" * 80)
        lines.append(f"  Grid dimensions: {mep_data.grid_dimensions[0]} × "
                    f"{mep_data.grid_dimensions[1]} × {mep_data.grid_dimensions[2]}")
        lines.append(f"  Grid origin (Å): ({mep_data.grid_origin[0]:.3f}, "
                    f"{mep_data.grid_origin[1]:.3f}, {mep_data.grid_origin[2]:.3f})")
        lines.append(f"  Grid spacing (Å): ({mep_data.grid_spacing[0]:.3f}, "
                    f"{mep_data.grid_spacing[1]:.3f}, {mep_data.grid_spacing[2]:.3f})")
        lines.append(f"  Total grid points: {np.prod(mep_data.grid_dimensions):,}")
        lines.append("")

        lines.append("MEP STATISTICS:")
        lines.append("-" * 80)
        if mep_data.min_potential is not None:
            lines.append(f"  Minimum potential: {mep_data.min_potential:12.6f} a.u. "
                        f"({mep_data.min_potential * 27.2114:10.3f} eV)")
        if mep_data.max_potential is not None:
            lines.append(f"  Maximum potential: {mep_data.max_potential:12.6f} a.u. "
                        f"({mep_data.max_potential * 27.2114:10.3f} eV)")
        if mep_data.mean_potential is not None:
            lines.append(f"  Mean potential:    {mep_data.mean_potential:12.6f} a.u. "
                        f"({mep_data.mean_potential * 27.2114:10.3f} eV)")
        lines.append("")

        if mep_data.local_minima:
            lines.append(f"NUCLEOPHILIC SITES (MEP minima, {len(mep_data.local_minima)} found):")
            lines.append("-" * 80)
            lines.append(f"{'X (Å)':<12} {'Y (Å)':<12} {'Z (Å)':<12} {'V (a.u.)':<12} {'V (eV)':<12}")
            lines.append("-" * 80)
            for x, y, z, v in sorted(mep_data.local_minima, key=lambda p: p[3])[:10]:
                lines.append(f"{x:<12.4f} {y:<12.4f} {z:<12.4f} {v:<12.6f} {v*27.2114:<12.3f}")
            lines.append("")

        if mep_data.local_maxima:
            lines.append(f"ELECTROPHILIC SITES (MEP maxima, {len(mep_data.local_maxima)} found):")
            lines.append("-" * 80)
            lines.append(f"{'X (Å)':<12} {'Y (Å)':<12} {'Z (Å)':<12} {'V (a.u.)':<12} {'V (eV)':<12}")
            lines.append("-" * 80)
            for x, y, z, v in sorted(mep_data.local_maxima, key=lambda p: p[3], reverse=True)[:10]:
                lines.append(f"{x:<12.4f} {y:<12.4f} {z:<12.4f} {v:<12.6f} {v*27.2114:<12.3f}")
            lines.append("")

        lines.append("=" * 80)
        lines.append("INTERPRETATION:")
        lines.append("  Negative MEP (blue):    Nucleophilic sites (electron-rich, attracts +)")
        lines.append("  Positive MEP (red):     Electrophilic sites (electron-poor, attracts -)")
        lines.append("  MEP minima:             Preferred sites for electrophile attack")
        lines.append("  MEP maxima:             Preferred sites for nucleophile attack")
        lines.append("=" * 80)

        return "\n".join(lines)
