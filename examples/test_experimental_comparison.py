#!/usr/bin/env python3
"""
Test Script for Experimental vs DFT Comparison Features

Tests:
- Item 15: Experimental vs DFT IR Comparison
- Item 30: Smart Dataset Label Positioning
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

print("🧪 Testing Experimental Comparison Features")
print("=" * 50)

# Test 1: Experimental vs DFT Comparison (Item 15)
print("\n1️⃣  Testing Experimental vs DFT IR Comparison...")
try:
    from experimental_comparison import (
        optimize_frequency_scaling,
        match_peaks,
        create_exp_vs_dft_comparison,
        create_peak_assignment_table
    )

    # Simulate realistic DFT IR data
    dft_frequencies = np.array([
        3750, 3100, 2980, 2920, 1740, 1645, 1480, 1370, 1260, 1090, 820
    ])
    dft_intensities = np.array([
        90, 35, 125, 110, 280, 60, 85, 70, 115, 100, 45
    ])

    # Simulate experimental data (with typical scaling offset and noise)
    x_exp = np.linspace(400, 4000, 500)
    y_exp = np.zeros_like(x_exp)

    # Add peaks at scaled DFT frequencies (typical scale factor ~0.97)
    true_scale = 0.97
    for freq, inten in zip(dft_frequencies, dft_intensities):
        scaled_freq = freq * true_scale
        # Add Gaussian peak with some noise
        y_exp += inten * np.exp(-(x_exp - scaled_freq)**2 / 400)

    # Add baseline and noise
    y_exp += 5 + np.random.normal(0, 2, len(x_exp))
    y_exp = np.maximum(y_exp, 0)  # No negative values

    experimental_data = {'x': x_exp, 'y': y_exp}

    # Test 1a: Optimize scaling factor
    print("   Testing scaling factor optimization...")
    best_scale, best_corr = optimize_frequency_scaling(
        x_exp, y_exp, dft_frequencies, dft_intensities,
        scale_range=(0.94, 1.00),
        num_trials=30
    )
    print(f"   ✓ Optimized scale: {best_scale:.4f} (true: {true_scale:.4f})")
    print(f"   ✓ Correlation: {best_corr:.4f}")

    # Test 1b: Peak matching
    print("   Testing peak matching...")
    # Find experimental peaks (simplified - just use known positions)
    exp_peaks = [(freq * true_scale, inten) for freq, inten in zip(dft_frequencies, dft_intensities)]
    calc_peaks = list(zip(dft_frequencies, dft_intensities))

    matches = match_peaks(exp_peaks, calc_peaks, scale_factor=best_scale, tolerance=20)
    print(f"   ✓ Matched {len(matches)} peaks")
    print(f"   ✓ Sample match: Exp {matches[0]['exp_freq']:.1f} cm⁻¹ → "
          f"Calc {matches[0]['calc_freq_scaled']:.1f} cm⁻¹ "
          f"(shift: {matches[0]['shift']:.1f} cm⁻¹)")

    # Test 1c: Create comparison plot
    print("   Creating comparison visualization...")
    fig = create_exp_vs_dft_comparison(
        experimental_data,
        dft_frequencies,
        dft_intensities,
        auto_scale=True,
        title="Test: Experimental vs DFT IR",
        save_path='test_exp_vs_dft_ir.png'
    )
    plt.close(fig)
    print("   ✅ Saved: test_exp_vs_dft_ir.png")

    # Test 1d: Create assignment table
    print("   Creating peak assignment table...")
    fig = create_peak_assignment_table(
        exp_peaks,
        calc_peaks,
        scale_factor=best_scale,
        tolerance=20,
        title="Test: Peak Assignment Table",
        save_path='test_peak_assignment.png'
    )
    plt.close(fig)
    print("   ✅ Saved: test_peak_assignment.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Smart Label Positioning (Item 30)
print("\n2️⃣  Testing Smart Dataset Label Positioning...")
try:
    from visualization_utils import (
        gaussian_broadening,
        normalize_spectrum,
        calculate_stack_offsets,
        find_optimal_label_positions
    )

    # Create multiple synthetic spectra with different peak patterns
    x = np.linspace(50, 3400, 2000)

    # Dataset 1: Peaks at low frequencies
    peaks1 = [(400, 1.0), (600, 0.8), (1000, 1.2)]
    spectrum1 = gaussian_broadening(x, peaks1, 17.0)
    spectrum1 = normalize_spectrum(spectrum1)

    # Dataset 2: Peaks at mid frequencies
    peaks2 = [(1200, 1.0), (1600, 1.5), (2000, 0.7)]
    spectrum2 = gaussian_broadening(x, peaks2, 17.0)
    spectrum2 = normalize_spectrum(spectrum2)

    # Dataset 3: Peaks at high frequencies
    peaks3 = [(2500, 0.9), (2900, 1.3), (3200, 0.6)]
    spectrum3 = gaussian_broadening(x, peaks3, 17.0)
    spectrum3 = normalize_spectrum(spectrum3)

    spectra = [spectrum1, spectrum2, spectrum3]
    labels = ['Low Freq', 'Mid Freq', 'High Freq']

    # Calculate offsets
    y_offsets = calculate_stack_offsets(spectra, y_space=0.1)

    # Test different label positioning methods
    methods = ['peak', 'flat', 'right']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, method in zip(axes, methods):
        # Find label positions
        positions = find_optimal_label_positions(
            spectra, y_offsets, (50, 3400), x, method=method
        )

        # Plot spectra
        for i, (spectrum, offset, label, color) in enumerate(
            zip(spectra, y_offsets, labels, ['red', 'blue', 'green'])):
            offset_spectrum = spectrum + offset
            ax.plot(x, offset_spectrum, color=color, linewidth=1.5)

            # Add label at optimal position
            x_pos, y_pos = positions[i]
            ax.annotate(label, xy=(x_pos, y_pos),
                       xytext=(10, 0), textcoords='offset points',
                       fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.4',
                                facecolor='white', edgecolor=color, linewidth=2),
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

        ax.invert_xaxis()
        ax.set_xlabel('Raman Shift (cm⁻¹)', fontweight='bold')
        ax.set_ylabel('Intensity (offset)', fontweight='bold')
        ax.set_title(f'Label Positioning: {method.upper()}', fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig('test_label_positioning.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"   ✓ Tested 3 label positioning methods: {', '.join(methods)}")
    print("   ✅ Saved: test_label_positioning.png")

except Exception as e:
    print(f"   ❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 50)
print("✅ All experimental comparison tests completed!")
print("\nGenerated files:")
test_files = [
    'test_exp_vs_dft_ir.png',
    'test_peak_assignment.png',
    'test_label_positioning.png'
]

for f in test_files:
    if os.path.exists(f):
        size = os.path.getsize(f) // 1024
        print(f"  ✓ {f} ({size}K)")
    else:
        print(f"  ✗ {f} (missing)")

print("\n🎉 Experimental comparison testing complete!")
