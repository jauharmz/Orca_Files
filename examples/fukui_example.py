"""
Example: Fukui Function Analysis

Demonstrates how to calculate and visualize Fukui functions for reactivity analysis.

Fukui functions identify reactive sites in molecules:
- f⁺: Sites susceptible to electrophilic attack
- f⁻: Sites susceptible to nucleophilic attack
- f⁰: Sites susceptible to radical attack

REQUIREMENTS:
To calculate Fukui functions, you need THREE ORCA calculations:
1. Neutral molecule (charge 0)
2. Cation (charge +1, N-1 electrons)
3. Anion (charge -1, N+1 electrons)

All three should use the SAME geometry and be single-point energy calculations.
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from parsers.out_parser import parse_orca_output
from parsers.fukui_calculator import FukuiCalculator, FukuiIndices
from visualization.fukui_visualization import (
    create_fukui_bar_plot,
    create_fukui_comparison_plot,
    create_global_descriptors_plot,
    create_fukui_heatmap
)
import matplotlib.pyplot as plt


def example_fukui_from_three_files():
    """
    Example: Calculate Fukui functions from three ORCA output files.

    This requires you to have performed three calculations:
    1. neutral.out - Neutral molecule
    2. cation.out - Cation (charge +1)
    3. anion.out - Anion (charge -1)
    """
    print("=" * 80)
    print("FUKUI FUNCTION ANALYSIS FROM THREE ORCA OUTPUTS")
    print("=" * 80)
    print()

    # File paths (UPDATE THESE to your actual file paths)
    neutral_file = "path/to/neutral.out"
    cation_file = "path/to/cation.out"
    anion_file = "path/to/anion.out"

    # Check if files exist
    for filepath in [neutral_file, cation_file, anion_file]:
        if not os.path.exists(filepath):
            print(f"ERROR: File not found: {filepath}")
            print("\nPlease update the file paths in this script to point to your ORCA outputs.")
            print("You need three calculations: neutral, cation (charge +1), and anion (charge -1)")
            return

    print("Parsing ORCA output files...")
    neutral_output = parse_orca_output(neutral_file)
    cation_output = parse_orca_output(cation_file)
    anion_output = parse_orca_output(anion_file)

    # Verify charges are correct
    print(f"  Neutral molecule: charge = {neutral_output.job_info.charge}")
    print(f"  Cation:          charge = {cation_output.job_info.charge}")
    print(f"  Anion:           charge = {anion_output.job_info.charge}")
    print()

    if neutral_output.job_info.charge != 0:
        print("WARNING: Neutral molecule should have charge = 0")
    if cation_output.job_info.charge != 1:
        print("WARNING: Cation should have charge = +1")
    if anion_output.job_info.charge != -1:
        print("WARNING: Anion should have charge = -1")

    # Calculate Fukui functions
    print("Calculating Fukui functions using Mulliken charges...")
    fukui_mulliken = FukuiCalculator.calculate_from_orca_outputs(
        neutral_output=neutral_output,
        cation_output=cation_output,
        anion_output=anion_output,
        use_loewdin=False
    )

    # Print summary
    summary = FukuiCalculator.print_summary(fukui_mulliken, top_n=5)
    print(summary)
    print()

    # Create visualizations
    print("Creating visualizations...")
    output_dir = "fukui_analysis"
    os.makedirs(output_dir, exist_ok=True)

    # 1. Bar plot of all Fukui indices
    fig1 = create_fukui_bar_plot(
        fukui_mulliken,
        output_file=f"{output_dir}/fukui_bar_plot.png",
        show_labels=True,
        top_n=20  # Show top 20 most reactive atoms
    )
    print(f"  Saved: {output_dir}/fukui_bar_plot.png")

    # 2. Comparison plot
    fig2 = create_fukui_comparison_plot(
        fukui_mulliken,
        output_file=f"{output_dir}/fukui_comparison.png",
        top_n=15
    )
    print(f"  Saved: {output_dir}/fukui_comparison.png")

    # 3. Global descriptors
    fig3 = create_global_descriptors_plot(
        fukui_mulliken,
        output_file=f"{output_dir}/global_descriptors.png"
    )
    print(f"  Saved: {output_dir}/global_descriptors.png")

    # 4. Heatmap
    fig4 = create_fukui_heatmap(
        fukui_mulliken,
        output_file=f"{output_dir}/fukui_heatmap.png",
        max_atoms=25
    )
    print(f"  Saved: {output_dir}/fukui_heatmap.png")

    plt.close('all')

    # Also calculate with Löwdin charges if available
    if neutral_output.loewdin_charges:
        print("\nCalculating Fukui functions using Löwdin charges...")
        fukui_loewdin = FukuiCalculator.calculate_from_orca_outputs(
            neutral_output=neutral_output,
            cation_output=cation_output,
            anion_output=anion_output,
            use_loewdin=True
        )

        summary_loewdin = FukuiCalculator.print_summary(fukui_loewdin, top_n=5)
        print(summary_loewdin)

        # Create comparison visualization
        fig5 = create_fukui_bar_plot(
            fukui_loewdin,
            output_file=f"{output_dir}/fukui_bar_plot_loewdin.png",
            show_labels=True,
            top_n=20
        )
        print(f"  Saved: {output_dir}/fukui_bar_plot_loewdin.png")
        plt.close('all')

    print()
    print("=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"Results saved to: {output_dir}/")
    print()
    print("INTERPRETATION GUIDE:")
    print("  • High f⁺ (red bars):  Electrophilic attack sites (electron-rich sites)")
    print("  • High f⁻ (blue bars): Nucleophilic attack sites (electron-poor sites)")
    print("  • High f⁰ (green bars): Radical attack sites")
    print("  • Positive Δf: Nucleophilic character (accepts electrons)")
    print("  • Negative Δf: Electrophilic character (donates electrons)")
    print()
    print("FOR DEGRADATION STUDIES:")
    print("  • Atoms with high f⁺ or f⁻ are most reactive")
    print("  • These sites are where degradation is likely to initiate")
    print("  • Compare with experimental degradation patterns")
    print("=" * 80)


def example_fukui_from_charges():
    """
    Example: Calculate Fukui functions directly from charge dictionaries.

    This is useful if you already have parsed charge data.
    """
    print("\n" + "=" * 80)
    print("FUKUI FUNCTION CALCULATION FROM CHARGE DATA")
    print("=" * 80)
    print()

    # Example charge data (replace with your actual data)
    # Format: {atom_index: (element, charge)}
    neutral_charges = {
        0: ('C', -0.123),
        1: ('C', -0.089),
        2: ('O', -0.456),
        3: ('H', 0.234),
        4: ('H', 0.221),
        5: ('H', 0.213)
    }

    cation_charges = {
        0: ('C', 0.012),
        1: ('C', 0.045),
        2: ('O', -0.234),
        3: ('H', 0.267),
        4: ('H', 0.254),
        5: ('H', 0.246)
    }

    anion_charges = {
        0: ('C', -0.245),
        1: ('C', -0.212),
        2: ('O', -0.612),
        3: ('H', 0.201),
        4: ('H', 0.189),
        5: ('H', 0.182)
    }

    # Energies in Hartree
    energy_neutral = -115.456789
    energy_cation = -115.123456
    energy_anion = -115.567890

    # Calculate Fukui functions
    fukui = FukuiCalculator.calculate_from_charges(
        neutral_charges=neutral_charges,
        cation_charges=cation_charges,
        anion_charges=anion_charges,
        method="mulliken",
        energy_neutral=energy_neutral,
        energy_cation=energy_cation,
        energy_anion=energy_anion
    )

    # Print summary
    summary = FukuiCalculator.print_summary(fukui, top_n=10)
    print(summary)

    # Print individual atom data
    print("\nDETAILED ATOMIC FUKUI INDICES:")
    print("-" * 80)
    print(f"{'Atom':<6} {'Element':<8} {'f⁺':<10} {'f⁻':<10} {'f⁰':<10} {'Δf':<10}")
    print("-" * 80)
    for atom in fukui.atomic_fukui:
        print(f"{atom.atom_index:<6} {atom.element:<8} "
              f"{atom.f_plus:<10.6f} {atom.f_minus:<10.6f} "
              f"{atom.f_zero:<10.6f} {atom.dual_descriptor:<10.6f}")
    print("-" * 80)


def example_identify_reactive_sites():
    """
    Example: Identify most reactive sites for different attack types.
    """
    print("\n" + "=" * 80)
    print("IDENTIFYING REACTIVE SITES")
    print("=" * 80)
    print()

    # Example: Using the charge data from previous example
    neutral_charges = {i: (f'C{i}', -0.1 + i*0.02) for i in range(10)}
    cation_charges = {i: (f'C{i}', 0.1 + i*0.01) for i in range(10)}
    anion_charges = {i: (f'C{i}', -0.3 + i*0.03) for i in range(10)}

    fukui = FukuiCalculator.calculate_from_charges(
        neutral_charges=neutral_charges,
        cation_charges=cation_charges,
        anion_charges=anion_charges,
        method="example"
    )

    # Get most reactive sites
    max_f_plus = fukui.get_max_f_plus()
    max_f_minus = fukui.get_max_f_minus()
    max_dual = fukui.get_max_dual_descriptor()

    print("MOST REACTIVE SITES:")
    print("-" * 80)
    if max_f_plus:
        print(f"Electrophilic attack site: Atom {max_f_plus.atom_index} ({max_f_plus.element})")
        print(f"  f⁺ = {max_f_plus.f_plus:.6f}")

    if max_f_minus:
        print(f"\nNucleophilic attack site: Atom {max_f_minus.atom_index} ({max_f_minus.element})")
        print(f"  f⁻ = {max_f_minus.f_minus:.6f}")

    if max_dual:
        print(f"\nMost ambiphilic site: Atom {max_dual.atom_index} ({max_dual.element})")
        print(f"  Δf = {max_dual.dual_descriptor:.6f}")
        if max_dual.dual_descriptor > 0:
            print("  Character: Nucleophilic (electron accepting)")
        else:
            print("  Character: Electrophilic (electron donating)")
    print("-" * 80)


if __name__ == "__main__":
    print("FUKUI FUNCTION ANALYSIS EXAMPLES")
    print("=" * 80)
    print()
    print("This script demonstrates how to calculate and visualize Fukui functions.")
    print()
    print("Available examples:")
    print("  1. Calculate from three ORCA output files (neutral, cation, anion)")
    print("  2. Calculate from charge dictionaries")
    print("  3. Identify most reactive sites")
    print()

    # Run example 2 (doesn't require files)
    example_fukui_from_charges()

    # Run example 3
    example_identify_reactive_sites()

    print("\n" + "=" * 80)
    print("To run the full analysis with real ORCA files:")
    print("  1. Update file paths in example_fukui_from_three_files()")
    print("  2. Uncomment the function call below")
    print("=" * 80)

    # Uncomment to run with real files (after updating paths)
    # example_fukui_from_three_files()
