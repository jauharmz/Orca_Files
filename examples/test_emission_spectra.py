#!/usr/bin/env python3
"""
Test Script for Emission Spectrum Visualization

Tests Item 20: Fluorescence vs Phosphorescence
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

print("🧪 Testing Emission Spectrum Visualization")
print("=" * 50)

# Test 1: Create synthetic emission data
print("\n1️⃣  Creating synthetic fluorescence and phosphorescence data...")
try:
    from emission_visualization import (
        create_emission_spectrum,
        calculate_stokes_shift,
        create_fluorescence_phosphorescence_comparison,
        create_absorption_emission_mirror
    )

    # Simulate fluorescence (S1 → S0)
    # Typically: narrow peaks, higher energy, short lifetime
    fluorescence_data = {
        'wavelengths': np.array([420, 450, 480]),  # nm
        'intensities': np.array([0.85, 1.0, 0.45])  # Relative intensities
    }

    # Simulate phosphorescence (T1 → S0)
    # Typically: broader peaks, lower energy (red-shifted), long lifetime
    phosphorescence_data = {
        'wavelengths': np.array([520, 560, 600]),  # nm (red-shifted)
        'intensities': np.array([0.65, 1.0, 0.35])  # Relative intensities
    }

    print(f"   ✓ Fluorescence: {len(fluorescence_data['wavelengths'])} transitions")
    print(f"   ✓ Phosphorescence: {len(phosphorescence_data['wavelengths'])} transitions")
    print(f"   ✓ Expected red shift: ~100 nm (T1 lower energy than S1)")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Calculate Stokes shift
print("\n2️⃣  Testing Stokes shift calculation...")
try:
    abs_max = 380.0  # Absorption maximum
    em_max = 450.0   # Emission maximum

    shift_nm, shift_cm1 = calculate_stokes_shift(abs_max, em_max)

    print(f"   ✓ Absorption max: {abs_max:.1f} nm")
    print(f"   ✓ Emission max: {em_max:.1f} nm")
    print(f"   ✓ Stokes shift: {shift_nm:.1f} nm ({shift_cm1:.0f} cm⁻¹)")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Create fluorescence vs phosphorescence comparison
print("\n3️⃣  Testing fluorescence vs phosphorescence comparison...")
try:
    fig = create_fluorescence_phosphorescence_comparison(
        fluorescence_data,
        phosphorescence_data,
        wavelength_range=(350, 700),
        fwhm_fluor=18.0,  # Narrower for fluorescence
        fwhm_phos=25.0,   # Broader for phosphorescence
        show_sticks=True,
        normalize=True,
        title="Fluorescence vs Phosphorescence Comparison",
        save_path='test_fluorescence_phosphorescence.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_fluorescence_phosphorescence.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 4: Create absorption-emission mirror image
print("\n4️⃣  Testing absorption-emission mirror image...")
try:
    # Simulate absorption spectrum
    absorption_data = {
        'wavelengths': np.array([350, 380, 410]),
        'intensities': np.array([0.45, 1.0, 0.35])
    }

    # Simulate emission spectrum (red-shifted)
    emission_data = {
        'wavelengths': np.array([420, 450, 480]),
        'intensities': np.array([0.40, 1.0, 0.50])
    }

    fig = create_absorption_emission_mirror(
        absorption_data,
        emission_data,
        wavelength_range=(300, 550),
        fwhm=20.0,
        title="Absorption-Emission Mirror Image (Stokes Shift)",
        save_path='test_absorption_emission_mirror.png'
    )
    plt.close(fig)

    print("   ✅ Saved: test_absorption_emission_mirror.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 5: Realistic organic molecule example
print("\n5️⃣  Testing realistic organic dye example...")
try:
    # Simulate a typical organic fluorophore (e.g., coumarin-like)
    # Absorption around 350-400 nm, fluorescence around 420-470 nm

    realistic_absorption = {
        'wavelengths': np.array([360, 375, 395]),
        'intensities': np.array([0.60, 1.0, 0.40])
    }

    realistic_fluorescence = {
        'wavelengths': np.array([430, 450, 475]),
        'intensities': np.array([0.70, 1.0, 0.50])
    }

    # Phosphorescence would be even more red-shifted (550-650 nm range)
    realistic_phosphorescence = {
        'wavelengths': np.array([560, 590, 620]),
        'intensities': np.array([0.50, 1.0, 0.45])
    }

    # Create comprehensive comparison
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Panel 1: Absorption
    wl_grid, abs_spec = create_emission_spectrum(
        np.array(realistic_absorption['wavelengths']),
        np.array(realistic_absorption['intensities']),
        (300, 700), 20.0
    )
    abs_spec = abs_spec / np.max(abs_spec)

    axes[0, 0].plot(wl_grid, abs_spec, 'b-', linewidth=2)
    axes[0, 0].fill_between(wl_grid, 0, abs_spec, color='blue', alpha=0.2)
    axes[0, 0].set_title('Absorption (S₀ → S₁)', fontweight='bold')
    axes[0, 0].set_xlabel('Wavelength (nm)')
    axes[0, 0].set_ylabel('Normalized Intensity')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_xlim(300, 700)

    # Panel 2: Fluorescence
    wl_grid, fluor_spec = create_emission_spectrum(
        np.array(realistic_fluorescence['wavelengths']),
        np.array(realistic_fluorescence['intensities']),
        (300, 700), 18.0
    )
    fluor_spec = fluor_spec / np.max(fluor_spec)

    axes[0, 1].plot(wl_grid, fluor_spec, color='#00FF00', linewidth=2)
    axes[0, 1].fill_between(wl_grid, 0, fluor_spec, color='green', alpha=0.2)
    axes[0, 1].set_title('Fluorescence (S₁ → S₀)', fontweight='bold')
    axes[0, 1].set_xlabel('Wavelength (nm)')
    axes[0, 1].set_ylabel('Normalized Intensity')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].set_xlim(300, 700)

    # Panel 3: Phosphorescence
    wl_grid, phos_spec = create_emission_spectrum(
        np.array(realistic_phosphorescence['wavelengths']),
        np.array(realistic_phosphorescence['intensities']),
        (300, 700), 25.0
    )
    phos_spec = phos_spec / np.max(phos_spec)

    axes[1, 0].plot(wl_grid, phos_spec, color='#FF4500', linewidth=2)
    axes[1, 0].fill_between(wl_grid, 0, phos_spec, color='red', alpha=0.2)
    axes[1, 0].set_title('Phosphorescence (T₁ → S₀)', fontweight='bold')
    axes[1, 0].set_xlabel('Wavelength (nm)')
    axes[1, 0].set_ylabel('Normalized Intensity')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_xlim(300, 700)

    # Panel 4: All overlaid
    axes[1, 1].plot(wl_grid, abs_spec, 'b-', linewidth=2, label='Absorption', alpha=0.7)
    axes[1, 1].plot(wl_grid, fluor_spec, color='#00FF00', linewidth=2,
                   label='Fluorescence', alpha=0.7)
    axes[1, 1].plot(wl_grid, phos_spec, color='#FF4500', linewidth=2,
                   label='Phosphorescence', alpha=0.7)

    # Calculate and show shifts
    abs_max = wl_grid[np.argmax(abs_spec)]
    fluor_max = wl_grid[np.argmax(fluor_spec)]
    phos_max = wl_grid[np.argmax(phos_spec)]

    axes[1, 1].axvline(abs_max, color='blue', linestyle='--', alpha=0.3)
    axes[1, 1].axvline(fluor_max, color='green', linestyle='--', alpha=0.3)
    axes[1, 1].axvline(phos_max, color='red', linestyle='--', alpha=0.3)

    axes[1, 1].set_title(f'Combined (Shifts: F={fluor_max-abs_max:.0f} nm, P={phos_max-abs_max:.0f} nm)',
                        fontweight='bold')
    axes[1, 1].set_xlabel('Wavelength (nm)')
    axes[1, 1].set_ylabel('Normalized Intensity')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].set_xlim(300, 700)

    plt.suptitle('Realistic Organic Fluorophore: Complete Photophysical Profile',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    fig.savefig('test_realistic_photophysics.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print("   ✅ Saved: test_realistic_photophysics.png")
    print(f"   ✓ Absorption max: {abs_max:.1f} nm")
    print(f"   ✓ Fluorescence max: {fluor_max:.1f} nm (Stokes: {fluor_max-abs_max:.0f} nm)")
    print(f"   ✓ Phosphorescence max: {phos_max:.1f} nm (Shift: {phos_max-abs_max:.0f} nm)")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 50)
print("✅ All emission spectrum tests completed!")
print("\nGenerated files:")
test_files = [
    'test_fluorescence_phosphorescence.png',
    'test_absorption_emission_mirror.png',
    'test_realistic_photophysics.png'
]

for f in test_files:
    if os.path.exists(f):
        size = os.path.getsize(f) // 1024
        print(f"  ✓ {f} ({size}K)")
    else:
        print(f"  ✗ {f} (missing)")

print("\n🎉 Emission spectrum testing complete!")
print("\n📊 Key Concepts Demonstrated:")
print("  • Fluorescence: Fast emission (ns), S₁ → S₀, narrow peaks")
print("  • Phosphorescence: Slow emission (ms-s), T₁ → S₀, broader peaks")
print("  • Stokes Shift: Energy loss between absorption and emission")
print("  • Red Shift: Phosphorescence occurs at longer wavelengths than fluorescence")
