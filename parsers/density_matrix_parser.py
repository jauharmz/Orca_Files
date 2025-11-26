"""
Density Matrix and MO Coefficient Parser

Parses density matrices and molecular orbital coefficients from ORCA outputs.

Data sources:
- .gbw files (binary, requires conversion)
- .molden files (MO coefficients)
- .out files (partial density data)
- orca_2mkl or orca_2json for .gbw extraction

Density matrices are needed for:
- Fukui function calculation
- Population analysis
- Bond order analysis
- Advanced property calculations
"""

import re
import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Tuple
import struct

logger = logging.getLogger(__name__)


@dataclass
class MOCoefficients:
    """Molecular orbital coefficients."""
    mo_index: int  # MO number (0-indexed)
    energy: float  # eV
    occupation: float  # 0, 1, or 2
    symmetry: Optional[str] = None  # e.g., "A1", "E", etc.
    coefficients: List[float] = field(default_factory=list)  # Coefficients for each AO
    ao_labels: List[str] = field(default_factory=list)  # Labels like "C1_2s", "O2_2pz"


@dataclass
class DensityMatrix:
    """Electron density matrix."""
    matrix_type: str = "total"  # "total", "alpha", "beta", "spin"
    dimension: int = 0  # Number of basis functions
    matrix: Optional[np.ndarray] = None  # 2D array (symmetric)

    # Metadata
    num_electrons: Optional[float] = None
    trace: Optional[float] = None  # Should equal num_electrons


@dataclass
class FockMatrix:
    """Fock matrix (effective one-electron Hamiltonian)."""
    matrix_type: str = "total"  # "total", "alpha", "beta"
    dimension: int = 0
    matrix: Optional[np.ndarray] = None  # 2D array (symmetric)


@dataclass
class OverlapMatrix:
    """Overlap matrix between basis functions."""
    dimension: int = 0
    matrix: Optional[np.ndarray] = None  # 2D array (symmetric)

    def calculate_orthogonality(self) -> float:
        """Calculate deviation from orthogonality."""
        if self.matrix is None:
            return 0.0

        # Perfect orthogonality: S = I (identity)
        identity = np.eye(self.dimension)
        deviation = np.linalg.norm(self.matrix - identity)
        return deviation


@dataclass
class ElectronicStructureData:
    """Complete electronic structure data."""
    # MO coefficients
    mo_coefficients: List[MOCoefficients] = field(default_factory=list)

    # Matrices
    density_matrix: Optional[DensityMatrix] = None
    density_matrix_alpha: Optional[DensityMatrix] = None
    density_matrix_beta: Optional[DensityMatrix] = None
    spin_density_matrix: Optional[DensityMatrix] = None
    fock_matrix: Optional[FockMatrix] = None
    overlap_matrix: Optional[OverlapMatrix] = None

    # Basis set information
    num_basis_functions: int = 0
    basis_function_labels: List[str] = field(default_factory=list)


class DensityMatrixParser:
    """Parser for density matrices and MO coefficients."""

    @staticmethod
    def parse_molden_file(molden_file: str) -> ElectronicStructureData:
        """
        Parse MO coefficients from .molden file.

        ORCA generates .molden files with:
            %output Print[P_MOs] 1 end

        Args:
            molden_file: Path to .molden file

        Returns:
            ElectronicStructureData object
        """
        logger.info(f"Parsing MOLDEN file: {molden_file}")

        data = ElectronicStructureData()

        try:
            with open(molden_file, 'r') as f:
                content = f.read()

            # Parse MO section
            mo_section_pattern = r'\[MO\](.*?)(?=\[|$)'
            mo_match = re.search(mo_section_pattern, content, re.DOTALL)

            if not mo_match:
                logger.warning("No [MO] section found in molden file")
                return data

            mo_section = mo_match.group(1)

            # Split into individual MOs
            # Each MO starts with "Sym=", "Ene=", "Spin=", "Occup="
            mo_blocks = re.split(r'(?=Sym\s*=)', mo_section.strip())

            for mo_idx, block in enumerate(mo_blocks):
                if not block.strip():
                    continue

                # Parse MO properties
                sym_match = re.search(r'Sym\s*=\s*(\S+)', block)
                ene_match = re.search(r'Ene\s*=\s*([-\d.Ee+]+)', block)
                occup_match = re.search(r'Occup\s*=\s*([-\d.Ee+]+)', block)

                if not (ene_match and occup_match):
                    continue

                symmetry = sym_match.group(1) if sym_match else None
                energy = float(ene_match.group(1))  # Already in eV or a.u. (check format)
                occupation = float(occup_match.group(1))

                # Parse coefficients
                # Format: "   1   0.12345"
                coef_pattern = r'\s*(\d+)\s+([-\d.Ee+]+)'
                coefficients = []
                for match in re.finditer(coef_pattern, block):
                    ao_index = int(match.group(1))
                    coef = float(match.group(2))
                    coefficients.append(coef)

                mo = MOCoefficients(
                    mo_index=mo_idx,
                    energy=energy,
                    occupation=occupation,
                    symmetry=symmetry,
                    coefficients=coefficients
                )
                data.mo_coefficients.append(mo)

            data.num_basis_functions = len(data.mo_coefficients[0].coefficients) if data.mo_coefficients else 0

            logger.info(f"Parsed {len(data.mo_coefficients)} molecular orbitals")
            logger.info(f"Basis functions: {data.num_basis_functions}")

            return data

        except Exception as e:
            logger.error(f"Error parsing MOLDEN file: {e}")
            raise

    @staticmethod
    def parse_density_from_output(output_file: str) -> Optional[ElectronicStructureData]:
        """
        Parse density matrix information from ORCA .out file.

        Note: The .out file typically doesn't contain full density matrices,
        but may contain partial data or references.

        Args:
            output_file: Path to ORCA .out file

        Returns:
            ElectronicStructureData object if data found, else None
        """
        logger.info(f"Searching for density matrix data in: {output_file}")

        try:
            with open(output_file, 'r') as f:
                content = f.read()

            data = ElectronicStructureData()

            # Look for density matrix printing (if enabled)
            # This is usually not printed by default; user must request it

            # Check for basis set size
            nbasis_match = re.search(r'Basis Dimension\s+Dim\s*=\s*(\d+)', content)
            if nbasis_match:
                data.num_basis_functions = int(nbasis_match.group(1))
                logger.info(f"Found basis dimension: {data.num_basis_functions}")

            # Look for "DENSITY MATRIX" section (rarely printed)
            density_pattern = r'DENSITY MATRIX\s*\n(.*?)(?=\n\s*\n|\Z)'
            density_match = re.search(density_pattern, content, re.DOTALL)

            if density_match:
                logger.info("Found density matrix section")
                # Parse matrix elements (implementation depends on ORCA output format)
                # This is a placeholder - actual parsing would need the exact format
                pass

            # If no density data found
            if data.num_basis_functions == 0:
                logger.warning("No density matrix data found in output file")
                return None

            return data

        except Exception as e:
            logger.error(f"Error parsing output file: {e}")
            return None

    @staticmethod
    def calculate_density_from_mos(
        mo_coefficients: List[MOCoefficients]
    ) -> DensityMatrix:
        """
        Calculate density matrix from MO coefficients.

        P_μν = Σ_i n_i C_μi C_νi

        where:
        - P is the density matrix
        - n_i is the occupation of MO i
        - C_μi is the coefficient of AO μ in MO i

        Args:
            mo_coefficients: List of MO coefficient objects

        Returns:
            DensityMatrix object
        """
        if not mo_coefficients:
            raise ValueError("No MO coefficients provided")

        nbasis = len(mo_coefficients[0].coefficients)

        # Initialize density matrix
        density = np.zeros((nbasis, nbasis))

        # Calculate density: P_μν = Σ_i n_i C_μi C_νi
        for mo in mo_coefficients:
            if len(mo.coefficients) != nbasis:
                logger.warning(f"MO {mo.mo_index} has inconsistent size")
                continue

            occupation = mo.occupation
            coef = np.array(mo.coefficients)

            # Outer product: C_μi C_νi
            density += occupation * np.outer(coef, coef)

        # Calculate trace (should equal total number of electrons)
        trace = np.trace(density)

        density_matrix = DensityMatrix(
            matrix_type="total",
            dimension=nbasis,
            matrix=density,
            num_electrons=trace,
            trace=trace
        )

        logger.info(f"Calculated density matrix ({nbasis}×{nbasis})")
        logger.info(f"Trace (electrons): {trace:.6f}")

        return density_matrix

    @staticmethod
    def calculate_spin_density(
        mo_alpha: List[MOCoefficients],
        mo_beta: List[MOCoefficients]
    ) -> DensityMatrix:
        """
        Calculate spin density matrix from alpha and beta MO coefficients.

        Spin density = P_alpha - P_beta

        Args:
            mo_alpha: Alpha MO coefficients
            mo_beta: Beta MO coefficients

        Returns:
            Spin density matrix
        """
        density_alpha = DensityMatrixParser.calculate_density_from_mos(mo_alpha)
        density_beta = DensityMatrixParser.calculate_density_from_mos(mo_beta)

        if density_alpha.matrix is None or density_beta.matrix is None:
            raise ValueError("Failed to calculate alpha/beta densities")

        spin_density = density_alpha.matrix - density_beta.matrix
        trace = np.trace(spin_density)

        return DensityMatrix(
            matrix_type="spin",
            dimension=density_alpha.dimension,
            matrix=spin_density,
            num_electrons=trace,  # Unpaired electrons
            trace=trace
        )

    @staticmethod
    def calculate_bond_order_from_density(
        density_matrix: DensityMatrix,
        overlap_matrix: OverlapMatrix,
        atom_i_aos: List[int],
        atom_j_aos: List[int]
    ) -> float:
        """
        Calculate bond order between two atoms from density and overlap matrices.

        Mayer bond order: B_IJ = Σ_μ∈I Σ_ν∈J (PS)_μν (PS)_νμ

        Args:
            density_matrix: Density matrix
            overlap_matrix: Overlap matrix
            atom_i_aos: List of AO indices for atom I
            atom_j_aos: List of AO indices for atom J

        Returns:
            Bond order value
        """
        if density_matrix.matrix is None or overlap_matrix.matrix is None:
            raise ValueError("Density or overlap matrix not available")

        P = density_matrix.matrix
        S = overlap_matrix.matrix

        # Calculate PS matrix
        PS = np.dot(P, S)

        # Calculate bond order
        bond_order = 0.0
        for mu in atom_i_aos:
            for nu in atom_j_aos:
                bond_order += PS[mu, nu] * PS[nu, mu]

        return bond_order

    @staticmethod
    def print_mo_summary(mo_coefficients: List[MOCoefficients], n_orbitals: int = 10) -> str:
        """Generate summary of MO coefficients."""
        lines = []
        lines.append("=" * 80)
        lines.append("MOLECULAR ORBITAL COEFFICIENTS SUMMARY")
        lines.append("=" * 80)
        lines.append("")

        # Find HOMO/LUMO
        occupied = [mo for mo in mo_coefficients if mo.occupation > 0.1]
        virtual = [mo for mo in mo_coefficients if mo.occupation < 0.1]

        if occupied:
            homo = max(occupied, key=lambda mo: mo.energy)
            lines.append(f"HOMO: MO {homo.mo_index}, Energy = {homo.energy:.4f} eV, "
                        f"Occupation = {homo.occupation:.2f}")

        if virtual:
            lumo = min(virtual, key=lambda mo: mo.energy)
            lines.append(f"LUMO: MO {lumo.mo_index}, Energy = {lumo.energy:.4f} eV, "
                        f"Occupation = {lumo.occupation:.2f}")

            if occupied:
                gap = lumo.energy - homo.energy
                lines.append(f"HOMO-LUMO Gap: {gap:.4f} eV")

        lines.append("")
        lines.append(f"OCCUPIED ORBITALS (showing up to {n_orbitals}):")
        lines.append("-" * 80)
        lines.append(f"{'MO':<6} {'Energy (eV)':<14} {'Occupation':<12} {'Symmetry':<12}")
        lines.append("-" * 80)

        for mo in sorted(occupied, key=lambda m: m.energy, reverse=True)[:n_orbitals]:
            sym = mo.symmetry if mo.symmetry else "N/A"
            lines.append(f"{mo.mo_index:<6} {mo.energy:<14.4f} {mo.occupation:<12.2f} {sym:<12}")

        lines.append("")
        lines.append(f"VIRTUAL ORBITALS (showing up to {n_orbitals}):")
        lines.append("-" * 80)
        lines.append(f"{'MO':<6} {'Energy (eV)':<14} {'Occupation':<12} {'Symmetry':<12}")
        lines.append("-" * 80)

        for mo in sorted(virtual, key=lambda m: m.energy)[:n_orbitals]:
            sym = mo.symmetry if mo.symmetry else "N/A"
            lines.append(f"{mo.mo_index:<6} {mo.energy:<14.4f} {mo.occupation:<12.2f} {sym:<12}")

        lines.append("")
        lines.append("=" * 80)

        return "\n".join(lines)
