"""
Experimental vs DFT Spectroscopy Comparison Tools

Specialized tools for comparing experimental spectroscopy data with
DFT-calculated spectra, including:
- Frequency scaling optimization
- Peak alignment and matching
- Statistical comparison metrics
- Side-by-side and overlay visualizations

Implements Item 15: Experimental vs DFT IR Comparison
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Callable
import logging

# Import local utilities
from visualization_utils import (
    gaussian_broadening,
    normalize_spectrum,
    align_spectra_to_common_grid,
    calculate_spectrum_similarity,
    interpolate_spectrum
)

logger = logging.getLogger(__name__)


def optimize_frequency_scaling(
    experimental_x: np.ndarray,
    experimental_y: np.ndarray,
    calculated_frequencies: np.ndarray,
    calculated_intensities: np.ndarray,
    scale_range: Tuple[float, float] = (0.94, 1.00),
    num_trials: int = 50,
    fwhm: float = 20.0
) -> Tuple[float, float]:
    """
    Find optimal frequency scaling factor to match experimental spectrum.

    Args:
        experimental_x: Experimental wavenumbers (cm⁻¹)
        experimental_y: Experimental intensities
        calculated_frequencies: DFT frequencies (cm⁻¹)
        calculated_intensities: DFT intensities
        scale_range: Range of scaling factors to test (min, max)
        num_trials: Number of scaling factors to test
        fwhm: FWHM for Gaussian broadening

    Returns:
        (best_scale_factor, best_correlation)
    """
    logger.info(f"Optimizing frequency scaling in range {scale_range}")

    # Normalize experimental
    exp_y_norm = normalize_spectrum(experimental_y)

    # Test different scaling factors
    scale_factors = np.linspace(scale_range[0], scale_range[1], num_trials)
    correlations = []

    # Common x-axis for comparison
    x_min = max(experimental_x.min(), calculated_frequencies.min() * scale_range[0])
    x_max = min(experimental_x.max(), calculated_frequencies.max() * scale_range[1])
    x_common = np.linspace(x_min, x_max, 1000)

    for scale in scale_factors:
        # Scale calculated frequencies
        scaled_freqs = calculated_frequencies * scale

        # Generate spectrum
        peaks = list(zip(scaled_freqs, calculated_intensities))
        calc_spectrum = gaussian_broadening(x_common, peaks, fwhm)
        calc_spectrum_norm = normalize_spectrum(calc_spectrum)

        # Interpolate experimental to common grid
        exp_interp = interpolate_spectrum(experimental_x, exp_y_norm, x_common)

        # Calculate correlation
        corr = calculate_spectrum_similarity(exp_interp, calc_spectrum_norm, 'correlation')
        correlations.append(corr)

    # Find best scaling factor
    best_idx = np.argmax(correlations)
    best_scale = scale_factors[best_idx]
    best_corr = correlations[best_idx]

    logger.info(f"Optimal scaling factor: {best_scale:.4f} (correlation: {best_corr:.4f})")

    return best_scale, best_corr


def match_peaks(
    experimental_peaks: List[Tuple[float, float]],
    calculated_peaks: List[Tuple[float, float]],
    scale_factor: float = 1.0,
    tolerance: float = 20.0
) -> List[Dict]:
    """
    Match experimental and calculated peaks within tolerance.

    Args:
        experimental_peaks: List of (frequency, intensity) for experimental
        calculated_peaks: List of (frequency, intensity) for DFT
        scale_factor: Scaling factor to apply to calculated frequencies
        tolerance: Maximum frequency difference (cm⁻¹) for matching

    Returns:
        List of matched peak dicts with 'exp_freq', 'calc_freq', 'shift', etc.
    """
    logger.info(f"Matching peaks with tolerance {tolerance} cm⁻¹")

    matches = []

    for exp_freq, exp_inten in experimental_peaks:
        best_match = None
        best_diff = float('inf')

        for calc_freq, calc_inten in calculated_peaks:
            scaled_freq = calc_freq * scale_factor
            diff = abs(exp_freq - scaled_freq)

            if diff < tolerance and diff < best_diff:
                best_diff = diff
                best_match = {
                    'exp_freq': exp_freq,
                    'exp_intensity': exp_inten,
                    'calc_freq': calc_freq,
                    'calc_freq_scaled': scaled_freq,
                    'calc_intensity': calc_inten,
                    'shift': exp_freq - scaled_freq,
                    'shift_abs': abs(exp_freq - scaled_freq)
                }

        if best_match is not None:
            matches.append(best_match)

    logger.info(f"Matched {len(matches)} peaks out of {len(experimental_peaks)} experimental peaks")

    return matches


def create_exp_vs_dft_comparison(
    experimental_data: Dict,
    calculated_frequencies: np.ndarray,
    calculated_intensities: np.ndarray,
    auto_scale: bool = True,
    manual_scale: float = 0.97,
    wavenumber_range: Tuple[float, float] = (400, 4000),
    fwhm: float = 20.0,
    figsize: Tuple[float, float] = (14, 12),
    title: str = "Experimental vs DFT IR Comparison",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create comprehensive experimental vs DFT IR comparison.

    Args:
        experimental_data: Dict with 'x' (wavenumbers) and 'y' (intensities)
        calculated_frequencies: DFT frequencies (cm⁻¹)
        calculated_intensities: DFT intensities (km/mol)
        auto_scale: If True, optimize scaling factor; if False, use manual_scale
        manual_scale: Manual scaling factor (typically 0.96-0.97)
        wavenumber_range: Range to plot
        fwhm: FWHM for Gaussian broadening
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info("Creating experimental vs DFT IR comparison")

    exp_x = np.array(experimental_data['x'])
    exp_y = np.array(experimental_data['y'])

    # Optimize scaling if requested
    if auto_scale:
        scale_factor, correlation = optimize_frequency_scaling(
            exp_x, exp_y, calculated_frequencies, calculated_intensities,
            fwhm=fwhm
        )
        scale_label = f"Optimized: {scale_factor:.4f}"
    else:
        scale_factor = manual_scale
        correlation = None
        scale_label = f"Manual: {scale_factor:.4f}"

    # Create figure with 4 panels
    fig = plt.figure(figsize=figsize)

    # Panel 1: Experimental spectrum (top)
    ax1 = plt.subplot(4, 1, 1)
    exp_y_norm = normalize_spectrum(exp_y)
    ax1.plot(exp_x, exp_y_norm, 'k-', linewidth=1.5, label='Experimental')
    ax1.invert_xaxis()
    ax1.set_ylabel('Normalized Intensity', fontweight='bold')
    ax1.set_title(f'{title} - Experimental', fontweight='bold')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(wavenumber_range)

    # Panel 2: Calculated spectrum (second)
    ax2 = plt.subplot(4, 1, 2)

    # Generate broadened calculated spectrum
    x_calc = np.linspace(wavenumber_range[0], wavenumber_range[1], 2000)
    scaled_freqs = calculated_frequencies * scale_factor
    peaks = list(zip(scaled_freqs, calculated_intensities))
    calc_spectrum = gaussian_broadening(x_calc, peaks, fwhm)
    calc_spectrum_norm = normalize_spectrum(calc_spectrum)

    ax2.plot(x_calc, calc_spectrum_norm, 'b-', linewidth=1.5,
            label=f'DFT (scaled: {scale_factor:.4f})')

    # Add stick spectrum
    for freq, inten in zip(scaled_freqs, calculated_intensities):
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            ax2.vlines(freq, 0, inten / calculated_intensities.max(),
                      color='red', alpha=0.4, linewidth=1.0)

    ax2.invert_xaxis()
    ax2.set_ylabel('Normalized Intensity', fontweight='bold')
    ax2.set_title(f'DFT Calculated - Scale Factor: {scale_label}', fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(wavenumber_range)

    # Panel 3: Overlay (third)
    ax3 = plt.subplot(4, 1, 3)

    # Align spectra for direct comparison
    x_common, [exp_aligned, calc_aligned] = align_spectra_to_common_grid(
        [(exp_x, exp_y_norm), (x_calc, calc_spectrum_norm)],
        x_min=wavenumber_range[0],
        x_max=wavenumber_range[1],
        num_points=1000
    )

    # Calculate metrics if not already done
    if correlation is None:
        correlation = calculate_spectrum_similarity(exp_aligned, calc_aligned, 'correlation')

    rmse = calculate_spectrum_similarity(exp_aligned, calc_aligned, 'rmse')
    r2 = calculate_spectrum_similarity(exp_aligned, calc_aligned, 'r2')

    ax3.plot(x_common, exp_aligned, 'k-', linewidth=1.5, alpha=0.7, label='Experimental')
    ax3.plot(x_common, calc_aligned, 'b-', linewidth=1.5, alpha=0.7, label='DFT')
    ax3.fill_between(x_common, exp_aligned, calc_aligned, alpha=0.2, color='gray')

    ax3.invert_xaxis()
    ax3.set_ylabel('Normalized Intensity', fontweight='bold')
    ax3.set_title(f'Overlay - Correlation: {correlation:.4f}, R²: {r2:.4f}, RMSE: {rmse:.4f}',
                 fontweight='bold')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(wavenumber_range)

    # Panel 4: Residuals (bottom)
    ax4 = plt.subplot(4, 1, 4)

    residuals = exp_aligned - calc_aligned
    ax4.plot(x_common, residuals, 'r-', linewidth=1.0)
    ax4.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax4.fill_between(x_common, 0, residuals, alpha=0.3, color='red')

    ax4.invert_xaxis()
    ax4.set_xlabel('Wavenumber (cm⁻¹)', fontweight='bold')
    ax4.set_ylabel('Residuals', fontweight='bold')
    ax4.set_title('Difference (Exp - DFT)', fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(wavenumber_range)

    # Style all axes
    for ax in [ax1, ax2, ax3, ax4]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


def create_peak_assignment_table(
    experimental_peaks: List[Tuple[float, float]],
    calculated_peaks: List[Tuple[float, float]],
    scale_factor: float = 0.97,
    tolerance: float = 20.0,
    figsize: Tuple[float, float] = (12, 8),
    title: str = "Peak Assignment Table",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create table showing experimental vs calculated peak assignments.

    Args:
        experimental_peaks: List of (frequency, intensity) for experimental
        calculated_peaks: List of (frequency, intensity) for DFT
        scale_factor: Scaling factor for calculated frequencies
        tolerance: Maximum frequency difference for matching
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info("Creating peak assignment table")

    # Match peaks
    matches = match_peaks(experimental_peaks, calculated_peaks, scale_factor, tolerance)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')

    # Prepare table data
    table_data = [['Exp. Freq.\n(cm⁻¹)', 'Exp. Inten.', 'Calc. Freq.\n(cm⁻¹)',
                   'Calc. Inten.\n(km/mol)', 'Scaled Freq.\n(cm⁻¹)', 'Shift\n(cm⁻¹)']]

    for match in matches[:20]:  # Show top 20 matches
        table_data.append([
            f"{match['exp_freq']:.1f}",
            f"{match['exp_intensity']:.2f}",
            f"{match['calc_freq']:.1f}",
            f"{match['calc_intensity']:.1f}",
            f"{match['calc_freq_scaled']:.1f}",
            f"{match['shift']:+.1f}"
        ])

    # Create table
    table = ax.table(cellText=table_data, cellLoc='center',
                    loc='center',
                    colWidths=[0.15, 0.15, 0.15, 0.18, 0.17, 0.15])

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2.5)

    # Style header row
    for i in range(6):
        cell = table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white')

    # Color-code shifts
    for i, match in enumerate(matches[:20], start=1):
        shift_abs = abs(match['shift'])
        if shift_abs < 5:
            color = '#C6EFCE'  # Green - excellent match
        elif shift_abs < 10:
            color = '#FFEB9C'  # Yellow - good match
        elif shift_abs < 20:
            color = '#FFC7CE'  # Red - acceptable match
        else:
            color = '#E0E0E0'  # Gray - poor match

        table[(i, 5)].set_facecolor(color)

    ax.set_title(f'{title}\nScale Factor: {scale_factor:.4f}, Tolerance: ±{tolerance} cm⁻¹',
                fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


# Export main functions
__all__ = [
    'optimize_frequency_scaling',
    'match_peaks',
    'create_exp_vs_dft_comparison',
    'create_peak_assignment_table'
]
