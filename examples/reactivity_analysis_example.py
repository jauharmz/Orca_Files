"""
Example: Comprehensive Reactivity Analysis for Degradation Studies

Combines Fukui functions and MEP analysis to predict degradation pathways.

This example demonstrates:
1. Fukui function calculation (frontier orbital reactivity)
2. MEP analysis (electrostatic reactivity)
3. Combined analysis for degradation site prediction
4. Correlation between Fukui and MEP

REQUIREMENTS:
For Fukui functions:
  - Three ORCA calculations: neutral, cation (+1), anion (-1)
  - Same geometry for all three
  - Single-point energy calculations

For MEP:
  - ORCA calculation with MEP cube file generation:
    %plots ElPot("mep.cube"); end

WORKFLOW FOR DEGRADATION STUDIES:
1. Perform geometry optimization
2. Run single-point calculations: neutral, cation, anion
3. Generate MEP cube file
4. Run this analysis script
5. Identify most reactive sites
6. Correlate with experimental degradation patterns
"""

import sys
import os
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from parsers.out_parser import parse_orca_output
from parsers.fukui_calculator import FukuiCalculator, FukuiIndices
from parsers.mep_parser import MEPParser, MEPData
from visualization.fukui_visualization import (
    create_fukui_bar_plot,
    create_global_descriptors_plot
)
from visualization.mep_visualization import (
    create_mep_slice_plot,
    create_mep_critical_points_plot
)
import matplotlib.pyplot as plt


def comprehensive_reactivity_analysis(
    neutral_file: str,
    cation_file: str,
    anion_file: str,
    mep_cube_file: str,
    output_dir: str = "reactivity_analysis"
):
    """
    Perform comprehensive reactivity analysis combining Fukui and MEP.

    Args:
        neutral_file: Path to neutral molecule .out file
        cation_file: Path to cation .out file
        anion_file: Path to anion .out file
        mep_cube_file: Path to MEP .cube file
        output_dir: Directory to save results
    """
    print("=" * 80)
    print("COMPREHENSIVE REACTIVITY ANALYSIS FOR DEGRADATION STUDIES")
    print("=" * 80)
    print()

    os.makedirs(output_dir, exist_ok=True)

    # ============================================================================
    # STEP 1: FUKUI FUNCTION ANALYSIS
    # ============================================================================
    print("STEP 1: FUKUI FUNCTION ANALYSIS")
    print("-" * 80)

    print("Parsing ORCA output files...")
    neutral_output = parse_orca_output(neutral_file)
    cation_output = parse_orca_output(cation_file)
    anion_output = parse_orca_output(anion_file)

    print("Calculating Fukui functions...")
    fukui = FukuiCalculator.calculate_from_orca_outputs(
        neutral_output=neutral_output,
        cation_output=cation_output,
        anion_output=anion_output,
        use_loewdin=False
    )

    # Print Fukui summary
    fukui_summary = FukuiCalculator.print_summary(fukui, top_n=10)
    print(fukui_summary)

    # Save Fukui summary to file
    with open(f"{output_dir}/fukui_summary.txt", 'w') as f:
        f.write(fukui_summary)

    # Create Fukui visualizations
    print("\nCreating Fukui visualizations...")
    fig1 = create_fukui_bar_plot(fukui, output_file=f"{output_dir}/fukui_indices.png")
    fig2 = create_global_descriptors_plot(fukui, output_file=f"{output_dir}/global_descriptors.png")
    plt.close('all')
    print(f"  Saved: {output_dir}/fukui_indices.png")
    print(f"  Saved: {output_dir}/global_descriptors.png")

    # ============================================================================
    # STEP 2: MEP ANALYSIS
    # ============================================================================
    print("\nSTEP 2: MEP ANALYSIS")
    print("-" * 80)

    print("Parsing MEP cube file...")
    mep_data = MEPParser.parse_cube_file(mep_cube_file)

    print("Finding MEP critical points...")
    mep_data.find_critical_points(threshold=0.01)

    # Print MEP summary
    mep_summary = MEPParser.print_summary(mep_data)
    print(mep_summary)

    # Save MEP summary to file
    with open(f"{output_dir}/mep_summary.txt", 'w') as f:
        f.write(mep_summary)

    # Create MEP visualizations
    print("\nCreating MEP visualizations...")
    fig3 = create_mep_slice_plot(mep_data, slice_axis='z',
                                 output_file=f"{output_dir}/mep_slice.png")
    fig4 = create_mep_critical_points_plot(mep_data,
                                          output_file=f"{output_dir}/mep_critical_points.png")
    plt.close('all')
    print(f"  Saved: {output_dir}/mep_slice.png")
    print(f"  Saved: {output_dir}/mep_critical_points.png")

    # ============================================================================
    # STEP 3: COMBINED ANALYSIS
    # ============================================================================
    print("\nSTEP 3: COMBINED REACTIVITY ANALYSIS")
    print("-" * 80)

    # Extract surface MEP for correlation with Fukui
    print("Extracting MEP on van der Waals surface...")
    surface_mep = MEPParser.extract_vdw_surface_mep(mep_data)

    # Identify consensus reactive sites
    print("\nIdentifying consensus reactive sites...")

    # Get top Fukui sites
    top_f_minus = sorted(fukui.atomic_fukui, key=lambda x: x.f_minus, reverse=True)[:5]
    top_f_plus = sorted(fukui.atomic_fukui, key=lambda x: x.f_plus, reverse=True)[:5]

    # Get coordinates for these atoms
    coords = {i: (elem, x, y, z) for i, (elem, x, y, z) in enumerate(neutral_output.coordinates)}

    # Calculate MEP near each Fukui site
    print("\nCORRELATION BETWEEN FUKUI AND MEP:")
    print("=" * 80)

    print("\nTOP NUCLEOPHILIC SITES (high f⁻, for nucleophilic attack):")
    print("-" * 80)
    print(f"{'Atom':<6} {'Element':<8} {'f⁻':<12} {'MEP (eV)':<12} {'Reactivity':<15}")
    print("-" * 80)

    for atom in top_f_minus:
        if atom.atom_index in coords:
            elem, x, y, z = coords[atom.atom_index]
            mep_val = MEPParser.calculate_mep_at_point(mep_data, x, y, z)

            if mep_val is not None:
                mep_ev = mep_val * 27.2114

                # Assess reactivity
                if atom.f_minus > 0.1 and mep_val < -0.01:
                    reactivity = "VERY HIGH"
                elif atom.f_minus > 0.05:
                    reactivity = "HIGH"
                else:
                    reactivity = "MODERATE"

                print(f"{atom.atom_index:<6} {elem:<8} {atom.f_minus:<12.6f} "
                      f"{mep_ev:<12.3f} {reactivity:<15}")

    print("\nTOP ELECTROPHILIC SITES (high f⁺, for electrophilic attack):")
    print("-" * 80)
    print(f"{'Atom':<6} {'Element':<8} {'f⁺':<12} {'MEP (eV)':<12} {'Reactivity':<15}")
    print("-" * 80)

    for atom in top_f_plus:
        if atom.atom_index in coords:
            elem, x, y, z = coords[atom.atom_index]
            mep_val = MEPParser.calculate_mep_at_point(mep_data, x, y, z)

            if mep_val is not None:
                mep_ev = mep_val * 27.2114

                # Assess reactivity
                if atom.f_plus > 0.1 and mep_val > 0.01:
                    reactivity = "VERY HIGH"
                elif atom.f_plus > 0.05:
                    reactivity = "HIGH"
                else:
                    reactivity = "MODERATE"

                print(f"{atom.atom_index:<6} {elem:<8} {atom.f_plus:<12.6f} "
                      f"{mep_ev:<12.3f} {reactivity:<15}")

    # ============================================================================
    # STEP 4: DEGRADATION PATHWAY PREDICTIONS
    # ============================================================================
    print("\n\nSTEP 4: DEGRADATION PATHWAY PREDICTIONS")
    print("=" * 80)

    degradation_report = []
    degradation_report.append("PREDICTED DEGRADATION SITES")
    degradation_report.append("=" * 80)
    degradation_report.append("")
    degradation_report.append("OXIDATIVE DEGRADATION (electron loss):")
    degradation_report.append("-" * 80)
    degradation_report.append("Atoms with high f⁺ are most susceptible to oxidation.")
    degradation_report.append("These are electron-rich sites that can easily lose electrons.")
    degradation_report.append("")

    for i, atom in enumerate(top_f_plus[:3], 1):
        elem = coords[atom.atom_index][0] if atom.atom_index in coords else "?"
        degradation_report.append(
            f"{i}. Atom {atom.atom_index} ({elem}): f⁺ = {atom.f_plus:.4f}"
        )

    degradation_report.append("")
    degradation_report.append("REDUCTIVE DEGRADATION (electron gain):")
    degradation_report.append("-" * 80)
    degradation_report.append("Atoms with high f⁻ are most susceptible to reduction.")
    degradation_report.append("These are electron-poor sites that can easily accept electrons.")
    degradation_report.append("")

    for i, atom in enumerate(top_f_minus[:3], 1):
        elem = coords[atom.atom_index][0] if atom.atom_index in coords else "?"
        degradation_report.append(
            f"{i}. Atom {atom.atom_index} ({elem}): f⁻ = {atom.f_minus:.4f}"
        )

    degradation_report.append("")
    degradation_report.append("ELECTROPHILIC ATTACK SITES (from MEP):")
    degradation_report.append("-" * 80)
    degradation_report.append("Regions with negative MEP attract electrophiles (H⁺, metal ions).")
    degradation_report.append("")

    if mep_data.local_minima:
        for i, (x, y, z, v) in enumerate(sorted(mep_data.local_minima, key=lambda p: p[3])[:3], 1):
            degradation_report.append(
                f"{i}. Position ({x:.2f}, {y:.2f}, {z:.2f} Å): "
                f"MEP = {v*27.2114:.2f} eV"
            )

    degradation_report.append("")
    degradation_report.append("NUCLEOPHILIC ATTACK SITES (from MEP):")
    degradation_report.append("-" * 80)
    degradation_report.append("Regions with positive MEP attract nucleophiles (OH⁻, H₂O).")
    degradation_report.append("")

    if mep_data.local_maxima:
        for i, (x, y, z, v) in enumerate(sorted(mep_data.local_maxima, key=lambda p: p[3], reverse=True)[:3], 1):
            degradation_report.append(
                f"{i}. Position ({x:.2f}, {y:.2f}, {z:.2f} Å): "
                f"MEP = {v*27.2114:.2f} eV"
            )

    degradation_report.append("")
    degradation_report.append("=" * 80)
    degradation_report.append("RECOMMENDATIONS:")
    degradation_report.append("1. Compare predicted sites with experimental degradation data")
    degradation_report.append("2. Sites with BOTH high Fukui AND extreme MEP are most reactive")
    degradation_report.append("3. Consider pH and environmental conditions for actual degradation")
    degradation_report.append("4. Use MD simulations to study dynamic effects")
    degradation_report.append("5. Perform transition state calculations for key pathways")
    degradation_report.append("=" * 80)

    degradation_text = "\n".join(degradation_report)
    print(degradation_text)

    # Save degradation report
    with open(f"{output_dir}/degradation_predictions.txt", 'w') as f:
        f.write(degradation_text)

    print(f"\n  Saved: {output_dir}/degradation_predictions.txt")

    # ============================================================================
    # SUMMARY
    # ============================================================================
    print("\n\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nAll results saved to: {output_dir}/")
    print("\nGenerated files:")
    print(f"  - fukui_summary.txt: Fukui function analysis")
    print(f"  - fukui_indices.png: Fukui indices visualization")
    print(f"  - global_descriptors.png: Global reactivity descriptors")
    print(f"  - mep_summary.txt: MEP analysis")
    print(f"  - mep_slice.png: MEP spatial distribution")
    print(f"  - mep_critical_points.png: MEP critical points")
    print(f"  - degradation_predictions.txt: Degradation pathway predictions")
    print("=" * 80)


def example_usage():
    """Example of how to use the comprehensive analysis."""

    # File paths (UPDATE THESE)
    neutral_file = "path/to/neutral.out"
    cation_file = "path/to/cation.out"
    anion_file = "path/to/anion.out"
    mep_cube_file = "path/to/mep.cube"

    # Check if files exist
    files = [neutral_file, cation_file, anion_file, mep_cube_file]
    missing = [f for f in files if not os.path.exists(f)]

    if missing:
        print("ERROR: The following files were not found:")
        for f in missing:
            print(f"  - {f}")
        print("\nPlease update the file paths in this script.")
        print("\nREQUIRED ORCA CALCULATIONS:")
        print("1. Neutral molecule (charge 0)")
        print("2. Cation (charge +1)")
        print("3. Anion (charge -1)")
        print("4. MEP cube file (add to input: %plots ElPot(\"mep.cube\"); end)")
        return

    # Run analysis
    comprehensive_reactivity_analysis(
        neutral_file=neutral_file,
        cation_file=cation_file,
        anion_file=anion_file,
        mep_cube_file=mep_cube_file,
        output_dir="degradation_analysis"
    )


if __name__ == "__main__":
    print("COMPREHENSIVE REACTIVITY ANALYSIS FOR DEGRADATION STUDIES")
    print("=" * 80)
    print()
    print("This script combines Fukui functions and MEP to predict degradation sites.")
    print()
    print("REQUIREMENTS:")
    print("  • Three ORCA outputs: neutral, cation (+1), anion (-1)")
    print("  • MEP cube file from ORCA")
    print()
    print("Update file paths in example_usage() function and run.")
    print("=" * 80)
    print()

    # Uncomment to run with real files
    # example_usage()
