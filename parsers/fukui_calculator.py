"""
Fukui Function Calculator

Calculates Fukui functions for reactivity analysis from ORCA outputs.
Fukui functions describe the sensitivity of molecular electron density to changes
in the number of electrons, indicating reactive sites.

For atom k:
- f⁺(k) = q(N) - q(N-1) : nucleophilic Fukui (electrophilic attack sites)
- f⁻(k) = q(N+1) - q(N) : electrophilic Fukui (nucleophilic attack sites)
- f⁰(k) = [q(N+1) - q(N-1)] / 2 : radical Fukui (radical attack sites)

Where q is the atomic charge and N is the number of electrons.
"""

import logging
from dataclasses import dataclass, field
from typing import Optional
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class AtomicFukui:
    """Fukui indices for a single atom."""
    atom_index: int
    element: str
    f_plus: float  # Nucleophilic Fukui (electrophilic attack)
    f_minus: float  # Electrophilic Fukui (nucleophilic attack)
    f_zero: float  # Radical Fukui
    dual_descriptor: float  # Δf = f⁻ - f⁺ (local softness)

    # Global reactivity indices (if available)
    local_softness_plus: Optional[float] = None
    local_softness_minus: Optional[float] = None
    local_electrophilicity: Optional[float] = None


@dataclass
class FukuiIndices:
    """Complete Fukui function analysis."""
    method: str = "mulliken"  # "mulliken" or "loewdin"
    atomic_fukui: list[AtomicFukui] = field(default_factory=list)

    # Global reactivity descriptors
    ionization_potential: Optional[float] = None  # IP = E(N-1) - E(N)
    electron_affinity: Optional[float] = None  # EA = E(N) - E(N+1)
    chemical_hardness: Optional[float] = None  # η = (IP - EA) / 2
    chemical_potential: Optional[float] = None  # μ = -(IP + EA) / 2
    electronegativity: Optional[float] = None  # χ = (IP + EA) / 2
    global_softness: Optional[float] = None  # S = 1 / (2η)
    electrophilicity_index: Optional[float] = None  # ω = μ² / (2η)

    # Energies used for calculation
    energy_neutral: Optional[float] = None
    energy_cation: Optional[float] = None
    energy_anion: Optional[float] = None

    def get_max_f_plus(self) -> Optional[AtomicFukui]:
        """Get atom with maximum f⁺ (most susceptible to electrophilic attack)."""
        if not self.atomic_fukui:
            return None
        return max(self.atomic_fukui, key=lambda x: x.f_plus)

    def get_max_f_minus(self) -> Optional[AtomicFukui]:
        """Get atom with maximum f⁻ (most susceptible to nucleophilic attack)."""
        if not self.atomic_fukui:
            return None
        return max(self.atomic_fukui, key=lambda x: x.f_minus)

    def get_max_dual_descriptor(self) -> Optional[AtomicFukui]:
        """Get atom with maximum dual descriptor (most ambiphilic)."""
        if not self.atomic_fukui:
            return None
        return max(self.atomic_fukui, key=lambda x: abs(x.dual_descriptor))


class FukuiCalculator:
    """Calculate Fukui functions from ORCA outputs."""

    @staticmethod
    def calculate_from_charges(
        neutral_charges: dict[int, tuple[str, float]],
        cation_charges: dict[int, tuple[str, float]],
        anion_charges: dict[int, tuple[str, float]],
        method: str = "mulliken",
        energy_neutral: Optional[float] = None,
        energy_cation: Optional[float] = None,
        energy_anion: Optional[float] = None
    ) -> FukuiIndices:
        """
        Calculate Fukui functions from atomic charges of N, N-1, N+1 systems.

        Args:
            neutral_charges: Charges from neutral molecule {atom_idx: (element, charge)}
            cation_charges: Charges from cation (N-1 electrons)
            anion_charges: Charges from anion (N+1 electrons)
            method: "mulliken" or "loewdin"
            energy_neutral: Total energy of neutral molecule (Hartree)
            energy_cation: Total energy of cation (Hartree)
            energy_anion: Total energy of anion (Hartree)

        Returns:
            FukuiIndices object with atomic and global descriptors
        """
        if not neutral_charges:
            raise ValueError("Neutral charges are required")

        # Check if all systems have the same atoms
        if set(neutral_charges.keys()) != set(cation_charges.keys()) or \
           set(neutral_charges.keys()) != set(anion_charges.keys()):
            logger.warning("Atom indices differ between charge sets")

        atomic_fukui_list = []

        # Calculate global reactivity descriptors if energies provided
        ionization_potential = None
        electron_affinity = None
        chemical_hardness = None
        chemical_potential = None
        electronegativity = None
        global_softness = None
        electrophilicity_index = None

        if energy_neutral is not None and energy_cation is not None:
            # IP = E(N-1) - E(N) in eV
            ionization_potential = (energy_cation - energy_neutral) * 27.2114  # Eh to eV

        if energy_neutral is not None and energy_anion is not None:
            # EA = E(N) - E(N+1) in eV
            electron_affinity = (energy_neutral - energy_anion) * 27.2114  # Eh to eV

        if ionization_potential is not None and electron_affinity is not None:
            # Chemical hardness: η = (IP - EA) / 2
            chemical_hardness = (ionization_potential - electron_affinity) / 2.0

            # Electronegativity: χ = (IP + EA) / 2
            electronegativity = (ionization_potential + electron_affinity) / 2.0

            # Chemical potential: μ = -χ
            chemical_potential = -electronegativity

            # Global softness: S = 1 / (2η)
            if chemical_hardness > 0:
                global_softness = 1.0 / (2.0 * chemical_hardness)

                # Electrophilicity index: ω = μ² / (2η)
                electrophilicity_index = (chemical_potential ** 2) / (2.0 * chemical_hardness)

        # Calculate atomic Fukui functions
        for atom_idx in sorted(neutral_charges.keys()):
            element, q_neutral = neutral_charges[atom_idx]

            # Get charges, defaulting to neutral charge if missing
            _, q_cation = cation_charges.get(atom_idx, (element, q_neutral))
            _, q_anion = anion_charges.get(atom_idx, (element, q_neutral))

            # Calculate Fukui functions
            # f⁺ = q(N) - q(N-1) for nucleophilic Fukui
            f_plus = q_neutral - q_cation

            # f⁻ = q(N+1) - q(N) for electrophilic Fukui
            f_minus = q_anion - q_neutral

            # f⁰ = [q(N+1) - q(N-1)] / 2 for radical Fukui
            f_zero = (q_anion - q_cation) / 2.0

            # Dual descriptor: Δf = f⁻ - f⁺
            # Positive: nucleophilic (electron accepting)
            # Negative: electrophilic (electron donating)
            dual_descriptor = f_minus - f_plus

            # Local softness if global softness available
            local_softness_plus = None
            local_softness_minus = None
            if global_softness is not None:
                local_softness_plus = f_plus * global_softness
                local_softness_minus = f_minus * global_softness

            # Local electrophilicity
            local_electrophilicity = None
            if electrophilicity_index is not None:
                local_electrophilicity = f_minus * electrophilicity_index

            atomic_fukui = AtomicFukui(
                atom_index=atom_idx,
                element=element,
                f_plus=f_plus,
                f_minus=f_minus,
                f_zero=f_zero,
                dual_descriptor=dual_descriptor,
                local_softness_plus=local_softness_plus,
                local_softness_minus=local_softness_minus,
                local_electrophilicity=local_electrophilicity
            )
            atomic_fukui_list.append(atomic_fukui)

        return FukuiIndices(
            method=method,
            atomic_fukui=atomic_fukui_list,
            ionization_potential=ionization_potential,
            electron_affinity=electron_affinity,
            chemical_hardness=chemical_hardness,
            chemical_potential=chemical_potential,
            electronegativity=electronegativity,
            global_softness=global_softness,
            electrophilicity_index=electrophilicity_index,
            energy_neutral=energy_neutral,
            energy_cation=energy_cation,
            energy_anion=energy_anion
        )

    @staticmethod
    def calculate_from_orca_outputs(
        neutral_output,
        cation_output,
        anion_output,
        use_loewdin: bool = False
    ) -> FukuiIndices:
        """
        Calculate Fukui functions from three OrcaOutput objects.

        Args:
            neutral_output: OrcaOutput from neutral molecule calculation
            cation_output: OrcaOutput from cation calculation (charge +1)
            anion_output: OrcaOutput from anion calculation (charge -1)
            use_loewdin: Use Löwdin charges instead of Mulliken (default: False)

        Returns:
            FukuiIndices object
        """
        method = "loewdin" if use_loewdin else "mulliken"

        # Extract charges based on method
        if use_loewdin:
            neutral_charges = neutral_output.loewdin_charges
            cation_charges = cation_output.loewdin_charges
            anion_charges = anion_output.loewdin_charges
        else:
            neutral_charges = neutral_output.mulliken_charges
            cation_charges = cation_output.mulliken_charges
            anion_charges = anion_output.mulliken_charges

        # Check if charges are available
        if not neutral_charges:
            raise ValueError(f"No {method.capitalize()} charges found in neutral output")
        if not cation_charges:
            raise ValueError(f"No {method.capitalize()} charges found in cation output")
        if not anion_charges:
            raise ValueError(f"No {method.capitalize()} charges found in anion output")

        # Extract energies
        energy_neutral = neutral_output.final_energy
        energy_cation = cation_output.final_energy
        energy_anion = anion_output.final_energy

        return FukuiCalculator.calculate_from_charges(
            neutral_charges=neutral_charges,
            cation_charges=cation_charges,
            anion_charges=anion_charges,
            method=method,
            energy_neutral=energy_neutral,
            energy_cation=energy_cation,
            energy_anion=energy_anion
        )

    @staticmethod
    def print_summary(fukui: FukuiIndices, top_n: int = 5) -> str:
        """
        Generate a text summary of Fukui analysis.

        Args:
            fukui: FukuiIndices object
            top_n: Number of top atoms to show for each index

        Returns:
            Formatted summary string
        """
        lines = []
        lines.append("=" * 80)
        lines.append(f"FUKUI FUNCTION ANALYSIS ({fukui.method.upper()} CHARGES)")
        lines.append("=" * 80)
        lines.append("")

        # Global descriptors
        lines.append("GLOBAL REACTIVITY DESCRIPTORS:")
        lines.append("-" * 80)
        if fukui.ionization_potential:
            lines.append(f"  Ionization Potential (IP):        {fukui.ionization_potential:10.4f} eV")
        if fukui.electron_affinity:
            lines.append(f"  Electron Affinity (EA):           {fukui.electron_affinity:10.4f} eV")
        if fukui.electronegativity:
            lines.append(f"  Electronegativity (χ):            {fukui.electronegativity:10.4f} eV")
        if fukui.chemical_hardness:
            lines.append(f"  Chemical Hardness (η):            {fukui.chemical_hardness:10.4f} eV")
        if fukui.chemical_potential:
            lines.append(f"  Chemical Potential (μ):           {fukui.chemical_potential:10.4f} eV")
        if fukui.global_softness:
            lines.append(f"  Global Softness (S):              {fukui.global_softness:10.4f} eV⁻¹")
        if fukui.electrophilicity_index:
            lines.append(f"  Electrophilicity Index (ω):       {fukui.electrophilicity_index:10.4f} eV")
        lines.append("")

        # Top electrophilic attack sites (high f⁺)
        lines.append(f"TOP {top_n} ELECTROPHILIC ATTACK SITES (high f⁺):")
        lines.append("-" * 80)
        lines.append(f"{'Atom':<8} {'Element':<8} {'f⁺':<12} {'f⁻':<12} {'f⁰':<12} {'Δf':<12}")
        lines.append("-" * 80)
        sorted_by_f_plus = sorted(fukui.atomic_fukui, key=lambda x: x.f_plus, reverse=True)
        for atom in sorted_by_f_plus[:top_n]:
            lines.append(
                f"{atom.atom_index:<8} {atom.element:<8} "
                f"{atom.f_plus:<12.6f} {atom.f_minus:<12.6f} "
                f"{atom.f_zero:<12.6f} {atom.dual_descriptor:<12.6f}"
            )
        lines.append("")

        # Top nucleophilic attack sites (high f⁻)
        lines.append(f"TOP {top_n} NUCLEOPHILIC ATTACK SITES (high f⁻):")
        lines.append("-" * 80)
        lines.append(f"{'Atom':<8} {'Element':<8} {'f⁺':<12} {'f⁻':<12} {'f⁰':<12} {'Δf':<12}")
        lines.append("-" * 80)
        sorted_by_f_minus = sorted(fukui.atomic_fukui, key=lambda x: x.f_minus, reverse=True)
        for atom in sorted_by_f_minus[:top_n]:
            lines.append(
                f"{atom.atom_index:<8} {atom.element:<8} "
                f"{atom.f_plus:<12.6f} {atom.f_minus:<12.6f} "
                f"{atom.f_zero:<12.6f} {atom.dual_descriptor:<12.6f}"
            )
        lines.append("")

        # Top radical attack sites (high f⁰)
        lines.append(f"TOP {top_n} RADICAL ATTACK SITES (high f⁰):")
        lines.append("-" * 80)
        lines.append(f"{'Atom':<8} {'Element':<8} {'f⁺':<12} {'f⁻':<12} {'f⁰':<12} {'Δf':<12}")
        lines.append("-" * 80)
        sorted_by_f_zero = sorted(fukui.atomic_fukui, key=lambda x: x.f_zero, reverse=True)
        for atom in sorted_by_f_zero[:top_n]:
            lines.append(
                f"{atom.atom_index:<8} {atom.element:<8} "
                f"{atom.f_plus:<12.6f} {atom.f_minus:<12.6f} "
                f"{atom.f_zero:<12.6f} {atom.dual_descriptor:<12.6f}"
            )
        lines.append("")

        lines.append("=" * 80)
        lines.append("INTERPRETATION:")
        lines.append("  f⁺ (nucleophilic Fukui):  High values → electrophilic attack sites")
        lines.append("  f⁻ (electrophilic Fukui): High values → nucleophilic attack sites")
        lines.append("  f⁰ (radical Fukui):       High values → radical attack sites")
        lines.append("  Δf (dual descriptor):     Positive → nucleophilic, Negative → electrophilic")
        lines.append("=" * 80)

        return "\n".join(lines)
