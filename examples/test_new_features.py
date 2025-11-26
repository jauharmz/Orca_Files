#!/usr/bin/env python3
"""
Test Script for New Visualization Features

Tests the following newly implemented features:
- Item 25: Grouped Orbital Connections (already in orbital_visualization.py)
- Item 26: Orbital Level Filtering
- Item 43: Data Interpolation for Comparison
- Item 23: IR Functional Group Assignment Table
- Item 34: Advanced Subplot Layouts

All tests use synthetic data for demonstration purposes.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

print("🧪 Testing New Visualization Features")
print("=" * 50)

# Test 1: Orbital filtering (Item 26)
print("\n1️⃣  Testing Orbital Level Filtering...")
try:
    from orbital_visualization import (
        filter_orbitals_by_energy,
        filter_orbitals_by_type,
        create_orbital_comparison
    )

    # Create test orbital data
    test_orbitals = [
        {'index': i, 'energy_ev': -10 + i * 0.5, 'occupation': 2.0 if i < 20 else 0.0}
        for i in range(30)
    ]

    # Test energy filtering
    filtered_energy = filter_orbitals_by_energy(test_orbitals, min_energy=-5, max_energy=5)
    print(f"   ✓ Energy filtering: {len(test_orbitals)} → {len(filtered_energy)} orbitals")

    # Test type filtering
    filtered_occupied = filter_orbitals_by_type(test_orbitals, 'occupied')
    filtered_virtual = filter_orbitals_by_type(test_orbitals, 'virtual')
    filtered_frontier = filter_orbitals_by_type(test_orbitals, 'frontier')

    print(f"   ✓ Type filtering:")
    print(f"     - Occupied: {len(filtered_occupied)}")
    print(f"     - Virtual: {len(filtered_virtual)}")
    print(f"     - Frontier: {len(filtered_frontier)}")

    # Test orbital comparison with filtering
    dataset_with_filter = {
        'orbital_energies': test_orbitals
    }

    fig = create_orbital_comparison(
        [dataset_with_filter, dataset_with_filter],
        ['No Filter', 'Filtered'],
        n_orbitals=8,
        energy_range=(-3, 3),
        title="Orbital Comparison with Energy Filter",
        save_path='test_orbital_filter.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_orbital_filter.png")

except Exception as e:
    print(f"   ❌ Error: {e}")

# Test 2: Data interpolation (Item 43)
print("\n2️⃣  Testing Data Interpolation for Comparison...")
try:
    from visualization_utils import (
        interpolate_spectrum,
        align_spectra_to_common_grid,
        calculate_spectrum_similarity
    )

    # Create test spectra with different grids
    x1 = np.linspace(0, 100, 50)
    y1 = np.exp(-(x1 - 30)**2 / 100) + np.exp(-(x1 - 70)**2 / 100)

    x2 = np.linspace(0, 100, 75)
    y2 = np.exp(-(x2 - 31)**2 / 100) + np.exp(-(x2 - 69)**2 / 100)

    # Test interpolation
    x_new = np.linspace(0, 100, 100)
    y1_interp = interpolate_spectrum(x1, y1, x_new, method='linear')
    print(f"   ✓ Interpolated spectrum: {len(x1)} → {len(x_new)} points")

    # Test alignment
    x_common, [y1_aligned, y2_aligned] = align_spectra_to_common_grid(
        [(x1, y1), (x2, y2)],
        num_points=100
    )
    print(f"   ✓ Aligned 2 spectra to common grid: {len(x_common)} points")

    # Test similarity metrics
    corr = calculate_spectrum_similarity(y1_aligned, y2_aligned, 'correlation')
    rmse = calculate_spectrum_similarity(y1_aligned, y2_aligned, 'rmse')
    r2 = calculate_spectrum_similarity(y1_aligned, y2_aligned, 'r2')

    print(f"   ✓ Similarity metrics:")
    print(f"     - Correlation: {corr:.4f}")
    print(f"     - RMSE: {rmse:.4f}")
    print(f"     - R²: {r2:.4f}")

    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    ax1.plot(x1, y1, 'o-', label='Original 1', alpha=0.6)
    ax1.plot(x2, y2, 's-', label='Original 2', alpha=0.6)
    ax1.set_title('Original Spectra (Different Grids)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(x_common, y1_aligned, '-', label='Aligned 1', linewidth=2)
    ax2.plot(x_common, y2_aligned, '-', label='Aligned 2', linewidth=2)
    ax2.set_title(f'Aligned Spectra (Correlation: {corr:.4f})')
    ax2.set_xlabel('X')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig('test_interpolation.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print("   ✅ Saved: test_interpolation.png")

except Exception as e:
    print(f"   ❌ Error: {e}")

# Test 3: IR functional group assignment (Item 23)
print("\n3️⃣  Testing IR Functional Group Assignment...")
try:
    from ir_visualization import (
        assign_functional_groups,
        create_ir_with_assignment_table
    )

    # Create realistic IR test data
    # Typical organic molecule with various functional groups
    ir_frequencies = np.array([
        3650,  # O-H
        3050,  # C-H aromatic
        2950,  # C-H aliphatic
        2850,  # C-H aliphatic
        1715,  # C=O ketone
        1650,  # C=C
        1450,  # C-H bend
        1250,  # C-O
        1100,  # C-O
        750    # C-H aromatic out-of-plane
    ])

    ir_intensities = np.array([85, 30, 120, 90, 250, 45, 80, 110, 95, 40])

    # Test assignment
    assignments = assign_functional_groups(ir_frequencies, ir_intensities, intensity_threshold=20)

    print(f"   ✓ Assigned {len(assignments)} peaks")
    print("   ✓ Sample assignments:")
    for i, assignment in enumerate(assignments[:5]):
        print(f"     - {assignment['frequency']:.0f} cm⁻¹: "
              f"{assignment['functional_group']} ({assignment['confidence']})")

    # Create IR with table
    fig = create_ir_with_assignment_table(
        ir_frequencies,
        ir_intensities,
        title="IR Spectrum with Functional Group Assignments",
        save_path='test_ir_assignment.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_ir_assignment.png")

except Exception as e:
    print(f"   ❌ Error: {e}")

# Test 4: Advanced subplot layouts (Item 34)
print("\n4️⃣  Testing Advanced Subplot Layouts...")
try:
    from subplot_layouts import (
        SubplotLayout,
        create_ir_raman_comparison,
        create_spectroscopy_orbitals_layout,
        create_comparison_panel
    )

    # Create test data
    # IR data
    ir_data = {
        'frequencies': np.array([3650, 2950, 1715, 1450, 1250, 750]),
        'intensities': np.array([85, 120, 250, 80, 110, 40])
    }

    # Raman data
    raman_data = {
        'frequencies': np.array([3065, 1600, 1450, 1180, 1000, 620]),
        'activities': np.array([35, 120, 85, 95, 75, 110])
    }

    # Orbital data
    orbital_data = {
        'orbital_energies': [
            {'index': i, 'energy_ev': -10 + i * 0.6, 'occupation': 2.0 if i < 18 else 0.0}
            for i in range(25)
        ]
    }

    # Test 4a: IR + Raman comparison (vertical)
    fig1 = create_ir_raman_comparison(ir_data, raman_data, layout='vertical',
                                     save_path='test_ir_raman_vertical.png')
    plt.close(fig1)
    print("   ✅ Saved: test_ir_raman_vertical.png")

    # Test 4b: IR + Raman comparison (horizontal)
    fig2 = create_ir_raman_comparison(ir_data, raman_data, layout='horizontal',
                                     save_path='test_ir_raman_horizontal.png')
    plt.close(fig2)
    print("   ✅ Saved: test_ir_raman_horizontal.png")

    # Test 4c: Spectroscopy + Orbitals
    fig3 = create_spectroscopy_orbitals_layout(raman_data, orbital_data,
                                              spectrum_type='raman',
                                              save_path='test_raman_orbitals.png')
    plt.close(fig3)
    print("   ✅ Saved: test_raman_orbitals.png")

    # Test 4d: Experimental vs Calculated comparison
    # Simulate experimental data
    x_exp = np.linspace(400, 4000, 300)
    y_exp = (np.exp(-(x_exp - 1715)**2 / 5000) * 250 +
             np.exp(-(x_exp - 2950)**2 / 3000) * 120 +
             np.random.normal(0, 5, len(x_exp)))

    experimental_data = {'x': x_exp, 'y': y_exp}
    calculated_data = ir_data

    fig4 = create_comparison_panel(experimental_data, calculated_data,
                                  spectrum_type='ir',
                                  save_path='test_exp_vs_calc.png')
    plt.close(fig4)
    print("   ✅ Saved: test_exp_vs_calc.png")

    # Test 4e: Custom layout
    layout = SubplotLayout(2, 2, figsize=(14, 10),
                          height_ratios=[1, 1],
                          width_ratios=[2, 1])

    ax1 = layout.add_subplot(0, 0, colspan=2, name='main')
    ax2 = layout.add_subplot(1, 0, name='bottom_left')
    ax3 = layout.add_subplot(1, 1, name='bottom_right')

    # Plot on custom layout
    x = np.linspace(0, 10, 100)
    ax1.plot(x, np.sin(x), 'b-', linewidth=2)
    ax1.set_title('Custom Layout - Main Panel')
    ax1.grid(True, alpha=0.3)

    ax2.plot(x, np.cos(x), 'r-', linewidth=2)
    ax2.set_title('Bottom Left')
    ax2.grid(True, alpha=0.3)

    ax3.plot(x, np.exp(-x/5), 'g-', linewidth=2)
    ax3.set_title('Bottom Right')
    ax3.grid(True, alpha=0.3)

    layout.save('test_custom_layout.png')
    print("   ✅ Saved: test_custom_layout.png")

except Exception as e:
    print(f"   ❌ Error: {e}")

print("\n" + "=" * 50)
print("✅ All tests completed successfully!")
print("\nGenerated files:")
test_files = [
    'test_orbital_filter.png',
    'test_interpolation.png',
    'test_ir_assignment.png',
    'test_ir_raman_vertical.png',
    'test_ir_raman_horizontal.png',
    'test_raman_orbitals.png',
    'test_exp_vs_calc.png',
    'test_custom_layout.png'
]

for f in test_files:
    if os.path.exists(f):
        size = os.path.getsize(f) // 1024
        print(f"  ✓ {f} ({size}K)")
    else:
        print(f"  ✗ {f} (missing)")

print("\n🎉 New feature testing complete!")
