#!/usr/bin/env python3
"""
Example: Using ORCA Visualization Modules

Demonstrates how to create publication-quality visualizations
from parsed ORCA output data.

Usage:
    python examples/visualization_example.py
"""

import sys
import os
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from parsers.out_parser import parse_out_file
from raman_visualization import create_stacked_raman_plot, create_single_raman_with_regions
from ir_visualization import create_stacked_ir_plot, create_dual_ir_plot
from orbital_visualization import create_orbital_comparison, create_simple_orbital_diagram

def example_single_raman():
    """Example: Single Raman spectrum with regions."""
    print("\n" + "="*60)
    print("Example 1: Single Raman Spectrum with Regional Annotations")
    print("="*60)

    # Parse ORCA output
    test_file = Path(__file__).parent.parent / "p1xs0p.out"
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return

    print(f"📄 Parsing: {test_file.name}")
    result = parse_out_file(str(test_file))

    # Extract Raman data
    if not result.raman_spectrum:
        print("❌ No Raman data in file")
        return

    import numpy as np
    frequencies = np.array([mode.frequency for mode in result.raman_spectrum])
    activities = np.array([mode.activity for mode in result.raman_spectrum])

    print(f"✓ Found {len(frequencies)} Raman modes")

    # Create visualization
    fig = create_single_raman_with_regions(
        frequencies=frequencies,
        activities=activities,
        shift_range=(0, 3500),
        fwhm=17.0,
        show_regions=True,
        show_labels=True,
        use_temperature_correction=True,
        title=f"Raman Spectrum - {result.job_info.basis_set}",
        save_path="raman_single_example.png"
    )

    print("✅ Saved: raman_single_example.png")
    print("   Features: Temperature-corrected intensity, regional boundaries, functional group labels")


def example_stacked_raman():
    """Example: Multi-dataset Raman stacking."""
    print("\n" + "="*60)
    print("Example 2: Multi-Dataset Stacked Raman Spectra")
    print("="*60)

    # Parse ORCA output
    test_file = Path(__file__).parent.parent / "p1xs0p.out"
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return

    result = parse_out_file(str(test_file))

    if not result.raman_spectrum:
        print("❌ No Raman data")
        return

    # Simulate multiple datasets (in real use, parse multiple files)
    # Here we'll create variations for demonstration
    import numpy as np
    base_freqs = np.array([mode.frequency for mode in result.raman_spectrum])
    base_acts = np.array([mode.activity for mode in result.raman_spectrum])

    datasets = [
        {'frequencies': base_freqs, 'activities': base_acts},
        {'frequencies': base_freqs * 0.96, 'activities': base_acts * 1.1},  # Scaled
        {'frequencies': base_freqs * 0.98, 'activities': base_acts * 0.9},  # Different method
    ]

    labels = [
        "Original",
        "Scaled (0.96x)",
        "Method B"
    ]

    print(f"✓ Preparing {len(datasets)} datasets for stacking")

    # Create stacked plot
    fig = create_stacked_raman_plot(
        datasets=datasets,
        labels=labels,
        shift_range=(200, 3400),
        fwhm=17.0,
        y_space=0.15,
        show_regions=True,
        show_labels=True,
        use_temperature_correction=True,
        title="Stacked Raman Spectra - Multi-Dataset Comparison",
        save_path="raman_stacked_example.png"
    )

    print("✅ Saved: raman_stacked_example.png")
    print("   Features: Automatic vertical spacing, stick spectra, dataset labels")


def example_dual_ir():
    """Example: Dual IR plot (absorbance + transmittance)."""
    print("\n" + "="*60)
    print("Example 3: Dual IR Spectrum (Absorbance + Transmittance)")
    print("="*60)

    test_file = Path(__file__).parent.parent / "p1xs0p.out"
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return

    result = parse_out_file(str(test_file))

    if not result.ir_spectrum:
        print("❌ No IR data")
        return

    import numpy as np
    frequencies = np.array([mode.frequency for mode in result.ir_spectrum])
    intensities = np.array([mode.intensity for mode in result.ir_spectrum])

    print(f"✓ Found {len(frequencies)} IR modes")

    # Create dual plot
    fig = create_dual_ir_plot(
        frequencies=frequencies,
        intensities=intensities,
        wavenumber_range=(500, 4000),
        fwhm=20.0,
        scale_factor=0.96,  # DFT frequency scaling
        shift_cm=0.0,
        show_regions=True,
        show_labels=True,
        title=f"IR Spectrum - {result.job_info.basis_set}",
        save_path="ir_dual_example.png"
    )

    print("✅ Saved: ir_dual_example.png")
    print("   Features: Beer-Lambert absorbance/transmittance, frequency scaling (0.96x), regional markers")


def example_stacked_ir():
    """Example: Multi-dataset IR stacking."""
    print("\n" + "="*60)
    print("Example 4: Multi-Dataset Stacked IR Spectra")
    print("="*60)

    test_file = Path(__file__).parent.parent / "p1xs0p.out"
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return

    result = parse_out_file(str(test_file))

    if not result.ir_spectrum:
        print("❌ No IR data")
        return

    import numpy as np
    base_freqs = np.array([mode.frequency for mode in result.ir_spectrum])
    base_ints = np.array([mode.intensity for mode in result.ir_spectrum])

    # Simulate multiple datasets
    datasets = [
        {'frequencies': base_freqs, 'intensities': base_ints},
        {'frequencies': base_freqs * 0.96, 'intensities': base_ints * 1.2},
        {'frequencies': base_freqs * 0.98, 'intensities': base_ints * 0.8},
    ]

    labels = ["Original", "Scaled (0.96x)", "Method B"]

    print(f"✓ Preparing {len(datasets)} datasets")

    # Create stacked plot (transmittance)
    fig = create_stacked_ir_plot(
        datasets=datasets,
        labels=labels,
        wavenumber_range=(600, 4000),
        fwhm=20.0,
        y_space=15.0,
        plot_type='transmittance',
        show_regions=True,
        show_labels=True,
        title="Stacked IR Spectra - Transmittance",
        save_path="ir_stacked_example.png"
    )

    print("✅ Saved: ir_stacked_example.png")
    print("   Features: Transmittance mode, automatic spacing, stick spectra above baseline")


def example_orbital_diagram():
    """Example: Simple orbital energy diagram."""
    print("\n" + "="*60)
    print("Example 5: Simple Molecular Orbital Energy Diagram")
    print("="*60)

    test_file = Path(__file__).parent.parent / "p1xs0p.out"
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return

    result = parse_out_file(str(test_file))

    if not result.orbital_energies:
        print("❌ No orbital data")
        return

    # Convert to dict format
    orbitals = [
        {
            'index': orb.index,
            'occupation': orb.occupation,
            'energy_ev': orb.energy_ev
        }
        for orb in result.orbital_energies
    ]

    print(f"✓ Found {len(orbitals)} molecular orbitals")

    # Create diagram
    fig = create_simple_orbital_diagram(
        orbital_energies=orbitals,
        n_orbitals=12,
        title=f"MO Diagram - {result.job_info.basis_set}",
        save_path="orbital_simple_example.png"
    )

    print("✅ Saved: orbital_simple_example.png")
    print("   Features: HOMO/LUMO marked, gap value, color-coded occupied/virtual")


def example_orbital_comparison():
    """Example: Multi-dataset orbital comparison."""
    print("\n" + "="*60)
    print("Example 6: Multi-Dataset Orbital Energy Comparison")
    print("="*60)

    test_file = Path(__file__).parent.parent / "p1xs0p.out"
    if not test_file.exists():
        print(f"❌ Test file not found: {test_file}")
        return

    result = parse_out_file(str(test_file))

    if not result.orbital_energies:
        print("❌ No orbital data")
        return

    # Convert to dict format
    base_orbitals = [
        {
            'index': orb.index,
            'occupation': orb.occupation,
            'energy_ev': orb.energy_ev
        }
        for orb in result.orbital_energies
    ]

    # Simulate multiple datasets (in real use, parse multiple files)
    import numpy as np
    datasets = [
        {'orbital_energies': base_orbitals},
        {'orbital_energies': [
            {**orb, 'energy_ev': orb['energy_ev'] + np.random.normal(0, 0.1)}
            for orb in base_orbitals
        ]},
        {'orbital_energies': [
            {**orb, 'energy_ev': orb['energy_ev'] + np.random.normal(0, 0.15)}
            for orb in base_orbitals
        ]},
    ]

    labels = ["Method A", "Method B", "Method C"]

    print(f"✓ Preparing {len(datasets)} datasets for comparison")

    # Create comparison
    fig = create_orbital_comparison(
        datasets=datasets,
        labels=labels,
        n_orbitals=8,
        show_gap_values=True,
        title="Orbital Energy Comparison - Multiple Methods",
        save_path="orbital_comparison_example.png"
    )

    print("✅ Saved: orbital_comparison_example.png")
    print("   Features: Side-by-side comparison, connection lines, gap annotations")


def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("ORCA Visualization Examples")
    print("="*70)
    print("\nThis script demonstrates all visualization capabilities.")
    print("Output files will be saved to the current directory.\n")

    try:
        # Run all examples
        example_single_raman()
        example_stacked_raman()
        example_dual_ir()
        example_stacked_ir()
        example_orbital_diagram()
        example_orbital_comparison()

        print("\n" + "="*70)
        print("✅ All examples completed successfully!")
        print("="*70)
        print("\nGenerated files:")
        print("  - raman_single_example.png")
        print("  - raman_stacked_example.png")
        print("  - ir_dual_example.png")
        print("  - ir_stacked_example.png")
        print("  - orbital_simple_example.png")
        print("  - orbital_comparison_example.png")
        print("\n")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
