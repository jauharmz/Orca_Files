"""
Multi-Dataset IR Spectrum Stacking and Advanced Visualization

Creates publication-quality stacked IR spectra with:
- Absorbance and transmittance plots
- Frequency scaling for DFT corrections
- Multi-dataset stacking with vertical offsets
- Functional group region markers

Based on implementation from 0cbz.ipynb.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple
import logging

# Import local utilities
from visualization_utils import (
    gaussian_broadening,
    calculate_stack_offsets,
    calculate_absorbance_transmittance,
    apply_frequency_scaling,
    get_color_palette,
    get_ir_region_boundaries,
    get_ir_region_labels,
    normalize_spectrum
)

logger = logging.getLogger(__name__)


def create_stacked_ir_plot(
    datasets: List[Dict],
    labels: List[str],
    wavenumber_range: Tuple[float, float] = (400, 4000),
    fwhm: float = 20.0,
    y_space: float = 15.0,
    scale_factor: float = 1.0,
    shift_cm: float = 0.0,
    figsize: Tuple[float, float] = (12, 8),
    colors: Optional[List[str]] = None,
    show_regions: bool = True,
    show_labels: bool = True,
    plot_type: str = 'absorbance',  # 'absorbance' or 'transmittance'
    concentration_mol_l: float = 0.01,
    path_length_cm: float = 1.0,
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create stacked IR spectra plot for multiple datasets.

    Args:
        datasets: List of dicts with 'frequencies' and 'intensities' keys
        labels: List of dataset labels
        wavenumber_range: (min, max) wavenumber range in cm⁻¹
        fwhm: Full width at half maximum for Gaussian broadening
        y_space: Vertical spacing between spectra
        scale_factor: Frequency scaling factor (e.g., 0.96 for B3LYP)
        shift_cm: Frequency shift in cm⁻¹
        figsize: Figure size (width, height)
        colors: Custom color list
        show_regions: Show regional boundary lines
        show_labels: Show dataset labels on right side
        plot_type: 'absorbance' or 'transmittance'
        concentration_mol_l: Sample concentration for Beer-Lambert
        path_length_cm: Path length for Beer-Lambert
        title: Plot title (auto-generated if None)
        save_path: Path to save figure

    Returns:
        matplotlib Figure object
    """
    logger.info(f"Creating stacked IR plot for {len(datasets)} datasets")

    # Validate inputs
    if len(datasets) != len(labels):
        raise ValueError(f"Number of datasets ({len(datasets)}) must match labels ({len(labels)})")

    # Generate wavenumber grid
    x = np.linspace(wavenumber_range[0], wavenumber_range[1], 5000)

    # Get colors
    if colors is None:
        colors = get_color_palette(len(datasets))

    # Process all datasets
    spectra = []
    for i, (dataset, label) in enumerate(zip(datasets, labels)):
        try:
            frequencies = np.array(dataset['frequencies'])
            intensities = np.array(dataset['intensities'])

            # Apply frequency scaling/shifting
            if scale_factor != 1.0 or shift_cm != 0.0:
                frequencies = apply_frequency_scaling(frequencies, scale_factor, shift_cm)
                logger.debug(f"Applied scaling: {scale_factor}x + {shift_cm} cm⁻¹")

            # Create peaks list
            peaks = list(zip(frequencies, intensities))

            # Apply Gaussian broadening
            spectrum = gaussian_broadening(x, peaks, fwhm)

            # Calculate absorbance/transmittance
            if plot_type == 'transmittance':
                _, spectrum = calculate_absorbance_transmittance(
                    spectrum, concentration_mol_l, path_length_cm
                )
            else:  # absorbance
                spectrum, _ = calculate_absorbance_transmittance(
                    spectrum, concentration_mol_l, path_length_cm
                )

            spectra.append(spectrum)
            logger.debug(f"Processed dataset {i}: {label}, {len(frequencies)} peaks")

        except Exception as e:
            logger.error(f"Error processing dataset {i} ({label}): {e}")
            raise

    # Calculate vertical offsets
    y_offsets = calculate_stack_offsets(spectra, y_space=y_space)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot each spectrum with offset
    for i, (spectrum, offset, label, color) in enumerate(zip(spectra, y_offsets, labels, colors)):
        offset_spectrum = spectrum + offset
        ax.plot(x, offset_spectrum, color=color, linewidth=1.5, label=label)

        # Add stick spectrum (above 100% for transmittance, above baseline for absorbance)
        if 'frequencies' in datasets[i]:
            freqs = datasets[i]['frequencies']
            if scale_factor != 1.0 or shift_cm != 0.0:
                freqs = apply_frequency_scaling(np.array(freqs), scale_factor, shift_cm)

            for freq in freqs:
                if wavenumber_range[0] <= freq <= wavenumber_range[1]:
                    if plot_type == 'transmittance':
                        stick_height = 5.0
                        stick_base = offset + np.max(spectrum[np.abs(x - freq) < 50])
                    else:
                        stick_height = np.max(spectrum) * 0.05
                        stick_base = offset

                    ax.vlines(freq, stick_base, stick_base + stick_height,
                            color=color, alpha=0.6, linewidth=0.8)

    # Add regional boundaries
    if show_regions:
        boundaries = get_ir_region_boundaries()
        y_min, y_max = ax.get_ylim()
        for boundary in boundaries:
            if wavenumber_range[0] <= boundary <= wavenumber_range[1]:
                ax.axvline(boundary, color='gray', linestyle='--',
                          alpha=0.5, linewidth=0.8, zorder=0)

    # Add dataset labels on right side
    if show_labels:
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(y_offsets)
        ax2.set_yticklabels(labels, fontsize=10)
        ax2.tick_params(axis='y', length=0)

    # Formatting
    ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12, fontweight='bold')

    if plot_type == 'transmittance':
        ylabel = 'Transmittance (%, offset)'
    else:
        ylabel = 'Absorbance (offset)'
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')

    # Auto-generate title if not provided
    if title is None:
        correction_str = f" (scaled: {scale_factor}×)" if scale_factor != 1.0 else ""
        title = f"Stacked IR Spectra - {plot_type.capitalize()}{correction_str}"
    ax.set_title(title, fontsize=14, fontweight='bold')

    ax.set_xlim(wavenumber_range)

    # Invert x-axis (high to low wavenumber)
    ax.invert_xaxis()

    # Hide left y-axis (arbitrary offsets)
    ax.get_yaxis().set_visible(False)
    ax.spines['left'].set_visible(False)

    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

    plt.tight_layout()

    # Save if requested
    if save_path:
        logger.info(f"Saving figure to {save_path}")
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    logger.info("Stacked IR plot created successfully")
    return fig


def create_dual_ir_plot(
    frequencies: np.ndarray,
    intensities: np.ndarray,
    wavenumber_range: Tuple[float, float] = (400, 4000),
    fwhm: float = 20.0,
    scale_factor: float = 1.0,
    shift_cm: float = 0.0,
    concentration_mol_l: float = 0.01,
    path_length_cm: float = 1.0,
    figsize: Tuple[float, float] = (12, 10),
    show_regions: bool = True,
    show_labels: bool = True,
    title: str = "IR Spectrum",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create dual-panel IR spectrum (absorbance + transmittance).

    Args:
        frequencies: IR frequencies (cm⁻¹)
        intensities: IR intensities (km/mol)
        wavenumber_range: Plot range
        fwhm: Gaussian broadening width
        scale_factor: Frequency scaling
        shift_cm: Frequency shift
        concentration_mol_l: Concentration for Beer-Lambert
        path_length_cm: Path length for Beer-Lambert
        figsize: Figure size
        show_regions: Show boundary lines
        show_labels: Show functional group labels
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    # Apply frequency scaling
    if scale_factor != 1.0 or shift_cm != 0.0:
        frequencies = apply_frequency_scaling(frequencies, scale_factor, shift_cm)

    # Generate grid
    x = np.linspace(wavenumber_range[0], wavenumber_range[1], 5000)

    # Build continuous spectrum
    peaks = list(zip(frequencies, intensities))
    spectrum = gaussian_broadening(x, peaks, fwhm)

    # Calculate absorbance and transmittance
    absorbance, transmittance = calculate_absorbance_transmittance(
        spectrum, concentration_mol_l, path_length_cm
    )

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    # Plot absorbance (top)
    ax1.plot(x, absorbance, color='blue', linewidth=1.5)
    ax1.fill_between(x, absorbance, alpha=0.2, color='blue')
    ax1.set_ylabel('Absorbance', fontsize=11, fontweight='bold')
    ax1.set_title(f"{title} - Absorbance", fontsize=12, fontweight='bold')

    # Stick spectrum for absorbance
    max_abs = np.max(absorbance)
    for freq, inten in zip(frequencies, intensities):
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            ax1.vlines(freq, max_abs * 1.05, max_abs * 1.15,
                      color='red', alpha=0.6, linewidth=1.0)

    # Plot transmittance (bottom)
    ax2.plot(x, transmittance, color='green', linewidth=1.5)
    ax2.fill_between(x, transmittance, alpha=0.2, color='green')
    ax2.set_ylabel('Transmittance (%)', fontsize=11, fontweight='bold')
    ax2.set_title(f"{title} - Transmittance", fontsize=12, fontweight='bold')
    ax2.set_xlabel('Wavenumber (cm⁻¹)', fontsize=11, fontweight='bold')

    # Stick spectrum for transmittance
    for freq, inten in zip(frequencies, intensities):
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            ax2.vlines(freq, 100, 105,
                      color='red', alpha=0.6, linewidth=1.0)

    # Add regional boundaries to both
    if show_regions:
        boundaries = get_ir_region_boundaries()
        for ax in [ax1, ax2]:
            for boundary in boundaries:
                if wavenumber_range[0] <= boundary <= wavenumber_range[1]:
                    ax.axvline(boundary, color='gray', linestyle='--',
                             alpha=0.5, linewidth=0.8)

    # Add functional group labels
    if show_labels:
        labels = get_ir_region_labels()
        y_pos_abs = max_abs * 0.95
        y_pos_trans = 95

        for freq, label in labels:
            if wavenumber_range[0] <= freq <= wavenumber_range[1]:
                ax1.text(freq, y_pos_abs, label, fontsize=8, ha='center',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                edgecolor='gray', alpha=0.7))
                ax2.text(freq, y_pos_trans, label, fontsize=8, ha='center',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                edgecolor='gray', alpha=0.7))

    # Formatting
    for ax in [ax1, ax2]:
        ax.invert_xaxis()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, alpha=0.3)

    ax1.set_ylim(bottom=0, top=max_abs * 1.2)
    ax2.set_ylim(0, 110)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig


# Export main functions
__all__ = ['create_stacked_ir_plot', 'create_dual_ir_plot']
