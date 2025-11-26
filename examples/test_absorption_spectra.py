#!/usr/bin/env python3
"""
Test Script for Absorption Spectrum Visualization

Tests Item 19: Absorption Spectrum Comparison (VG, AH, AHAS)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

print("🧪 Testing Absorption Spectrum Visualization")
print("=" * 50)

# Test 1: Create synthetic absorption data
print("\n1️⃣  Creating synthetic absorption spectra data...")
try:
    from absorption_visualization import (
        create_broadened_absorption_spectrum,
        extract_absorption_data,
        create_absorption_comparison,
        create_electric_vs_velocity_comparison,
        create_transition_table
    )

    # Simulate VG (Vertical Gradient) method
    # Typically shows sharp, intense peaks
    vg_transitions = [
        {'wavelength_nm': 280, 'energy_ev': 4.43, 'fosc_d2': 0.85},
        {'wavelength_nm': 320, 'energy_ev': 3.88, 'fosc_d2': 0.35},
        {'wavelength_nm': 410, 'energy_ev': 3.02, 'fosc_d2': 0.12},
        {'wavelength_nm': 450, 'energy_ev': 2.76, 'fosc_d2': 0.05},
    ]

    # Simulate AH (Adiabatic Hessian) method
    # Usually shows red-shifted peaks with vibronic structure
    ah_transitions = [
        {'wavelength_nm': 295, 'energy_ev': 4.20, 'fosc_d2': 0.75},
        {'wavelength_nm': 310, 'energy_ev': 4.00, 'fosc_d2': 0.15},
        {'wavelength_nm': 340, 'energy_ev': 3.65, 'fosc_d2': 0.28},
        {'wavelength_nm': 360, 'energy_ev': 3.44, 'fosc_d2': 0.08},
        {'wavelength_nm': 430, 'energy_ev': 2.88, 'fosc_d2': 0.10},
        {'wavelength_nm': 470, 'energy_ev': 2.64, 'fosc_d2': 0.04},
    ]

    # Simulate AHAS (AH + Anharmonic + Solvation) method
    # Shows broader peaks, further red-shifted
    ahas_transitions = [
        {'wavelength_nm': 305, 'energy_ev': 4.07, 'fosc_d2': 0.68},
        {'wavelength_nm': 325, 'energy_ev': 3.82, 'fosc_d2': 0.18},
        {'wavelength_nm': 350, 'energy_ev': 3.54, 'fosc_d2': 0.25},
        {'wavelength_nm': 375, 'energy_ev': 3.31, 'fosc_d2': 0.09},
        {'wavelength_nm': 445, 'energy_ev': 2.79, 'fosc_d2': 0.08},
        {'wavelength_nm': 490, 'energy_ev': 2.53, 'fosc_d2': 0.03},
    ]

    print(f"   ✓ VG method: {len(vg_transitions)} transitions")
    print(f"   ✓ AH method: {len(ah_transitions)} transitions")
    print(f"   ✓ AHAS method: {len(ahas_transitions)} transitions")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Create multi-method comparison (overlay)
print("\n2️⃣  Testing multi-method comparison (overlay)...")
try:
    datasets = [
        {'transitions': vg_transitions},
        {'transitions': ah_transitions},
        {'transitions': ahas_transitions}
    ]
    labels = ['VG', 'AH', 'AHAS']

    fig = create_absorption_comparison(
        datasets, labels,
        wavelength_range=(250, 550),
        fwhm=15.0,
        layout='overlay',
        show_sticks=True,
        normalize=True,
        title="Absorption Spectrum Comparison: VG vs AH vs AHAS",
        save_path='test_absorption_overlay.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_absorption_overlay.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Create stacked comparison
print("\n3️⃣  Testing stacked comparison with metrics...")
try:
    fig = create_absorption_comparison(
        datasets, labels,
        wavelength_range=(250, 550),
        fwhm=15.0,
        layout='stacked',
        show_sticks=True,
        normalize=True,
        title="Stacked Absorption Spectra: VG, AH, AHAS",
        save_path='test_absorption_stacked.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_absorption_stacked.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 4: Electric vs Velocity dipole comparison
print("\n4️⃣  Testing electric vs velocity dipole comparison...")
try:
    # Simulate electric dipole transitions
    electric_transitions = [
        {'wavelength_nm': 290, 'energy_ev': 4.28, 'fosc_d2': 0.80},
        {'wavelength_nm': 330, 'energy_ev': 3.76, 'fosc_d2': 0.32},
        {'wavelength_nm': 420, 'energy_ev': 2.95, 'fosc_d2': 0.15},
    ]

    # Simulate velocity dipole transitions (usually similar but not identical)
    velocity_transitions = [
        {'wavelength_nm': 290, 'energy_ev': 4.28, 'fosc_p2': 0.78},
        {'wavelength_nm': 330, 'energy_ev': 3.76, 'fosc_p2': 0.34},
        {'wavelength_nm': 420, 'energy_ev': 2.95, 'fosc_p2': 0.14},
    ]

    fig = create_electric_vs_velocity_comparison(
        electric_transitions,
        velocity_transitions,
        wavelength_range=(250, 500),
        fwhm=15.0,
        title="Electric vs Velocity Dipole Comparison",
        save_path='test_electric_vs_velocity.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_electric_vs_velocity.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 5: Transition table
print("\n5️⃣  Testing transition table generation...")
try:
    # Use AH transitions as example
    fig = create_transition_table(
        ah_transitions,
        num_transitions=6,
        title="Major Transitions - AH Method",
        save_path='test_transition_table.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_transition_table.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 6: Broadening function
print("\n6️⃣  Testing Gaussian broadening function...")
try:
    # Extract data from VG
    from absorption_visualization import extract_absorption_data

    wl, en, fosc = extract_absorption_data(vg_transitions)
    print(f"   ✓ Extracted {len(wl)} transitions")
    print(f"   ✓ Wavelength range: {wl.min():.1f} - {wl.max():.1f} nm")
    print(f"   ✓ Total oscillator strength: {fosc.sum():.4f}")

    # Create broadened spectrum
    wl_grid, spectrum = create_broadened_absorption_spectrum(
        wl, fosc,
        wavelength_range=(250, 550),
        fwhm=15.0
    )

    print(f"   ✓ Created broadened spectrum with {len(wl_grid)} points")
    print(f"   ✓ Maximum intensity: {spectrum.max():.4f}")

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Stick spectrum
    for wavelength, strength in zip(wl, fosc):
        ax.vlines(wavelength, 0, strength, color='red', alpha=0.6, linewidth=2, label='Discrete')

    # Broadened spectrum
    ax.plot(wl_grid, spectrum, 'b-', linewidth=2, label='Broadened')

    ax.set_xlabel('Wavelength (nm)', fontweight='bold')
    ax.set_ylabel('Oscillator Strength', fontweight='bold')
    ax.set_title('Gaussian Broadening Demonstration', fontweight='bold')
    ax.set_xlim(250, 550)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig('test_broadening_demo.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print("   ✅ Saved: test_broadening_demo.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 50)
print("✅ All absorption spectrum tests completed!")
print("\nGenerated files:")
test_files = [
    'test_absorption_overlay.png',
    'test_absorption_stacked.png',
    'test_electric_vs_velocity.png',
    'test_transition_table.png',
    'test_broadening_demo.png'
]

for f in test_files:
    if os.path.exists(f):
        size = os.path.getsize(f) // 1024
        print(f"  ✓ {f} ({size}K)")
    else:
        print(f"  ✗ {f} (missing)")

print("\n🎉 Absorption spectrum testing complete!")
