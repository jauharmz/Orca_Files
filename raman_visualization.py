"""
Multi-Dataset Raman Spectrum Stacking

Creates publication-quality stacked Raman spectra with vertical offset
for comparing multiple datasets. Based on implementation from 0cbz.ipynb.

Features:
- Automatic vertical spacing with minimal overlap
- Gaussian broadening of discrete peaks
- Temperature-dependent intensity conversion
- Regional boundary markers
- Right-side dataset labels
- Customizable colors and styling
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple
import logging

# Import local utilities
from visualization_utils import (
    gaussian_broadening,
    calculate_stack_offsets,
    raman_activity_to_intensity,
    get_color_palette,
    get_raman_region_boundaries,
    get_raman_region_labels,
    normalize_spectrum
)

logger = logging.getLogger(__name__)


def create_stacked_raman_plot(
    datasets: List[Dict],
    labels: List[str],
    shift_range: Tuple[float, float] = (50, 3400),
    fwhm: float = 17.0,
    y_space: float = 0.1,
    figsize: Tuple[float, float] = (12, 8),
    colors: Optional[List[str]] = None,
    show_regions: bool = True,
    show_labels: bool = True,
    use_temperature_correction: bool = True,
    laser_wavelength_nm: float = 532,
    low_freq_cutoff: float = 0.0,
    title: str = "Stacked Raman Spectra",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create stacked Raman spectra plot for multiple datasets.

    Args:
        datasets: List of dicts with 'frequencies' and 'activities' keys
        labels: List of dataset labels
        shift_range: (min, max) Raman shift range in cm⁻¹
        fwhm: Full width at half maximum for Gaussian broadening
        y_space: Vertical spacing between spectra
        figsize: Figure size (width, height)
        colors: Custom color list (default: auto-generated)
        show_regions: Show regional boundary lines
        show_labels: Show dataset labels on right side
        use_temperature_correction: Apply temperature-dependent intensity
        laser_wavelength_nm: Laser wavelength for intensity correction
        low_freq_cutoff: Suppress modes below this frequency
        title: Plot title
        save_path: Path to save figure (optional)

    Returns:
        matplotlib Figure object
    """
    logger.info(f"Creating stacked Raman plot for {len(datasets)} datasets")

    # Validate inputs
    if len(datasets) != len(labels):
        raise ValueError(f"Number of datasets ({len(datasets)}) must match labels ({len(labels)})")

    # Generate Raman shift grid
    shifts = np.linspace(shift_range[0], shift_range[1], 5000)

    # Get colors
    if colors is None:
        colors = get_color_palette(len(datasets))

    # Process all datasets
    spectra = []
    for i, (dataset, label) in enumerate(zip(datasets, labels)):
        try:
            frequencies = np.array(dataset['frequencies'])
            activities = np.array(dataset['activities'])

            # Apply temperature correction if requested
            if use_temperature_correction:
                intensities = raman_activity_to_intensity(
                    frequencies, activities,
                    laser_wavelength_nm=laser_wavelength_nm,
                    low_freq_cutoff=low_freq_cutoff
                )
            else:
                intensities = activities

            # Create peaks list
            peaks = list(zip(frequencies, intensities))

            # Apply Gaussian broadening
            spectrum = gaussian_broadening(shifts, peaks, fwhm)

            # Normalize
            spectrum = normalize_spectrum(spectrum, method='max')

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
        ax.plot(shifts, offset_spectrum, color=color, linewidth=1.5, label=label)

        # Add stick spectrum below baseline
        if 'frequencies' in datasets[i]:
            freqs = datasets[i]['frequencies']
            for freq in freqs:
                if shift_range[0] <= freq <= shift_range[1]:
                    stick_height = 0.05  # Relative to normalized spectrum
                    ax.vlines(freq, offset - stick_height, offset,
                            color=color, alpha=0.6, linewidth=0.8)

    # Add regional boundaries
    if show_regions:
        boundaries = get_raman_region_boundaries()
        y_min, y_max = ax.get_ylim()
        for boundary in boundaries:
            if shift_range[0] <= boundary <= shift_range[1]:
                ax.axvline(boundary, color='gray', linestyle='--',
                          alpha=0.5, linewidth=0.8, zorder=0)

    # Add dataset labels on right side
    if show_labels:
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(y_offsets)
        ax2.set_yticklabels(labels, fontsize=10)
        ax2.tick_params(axis='y', length=0)  # Hide tick marks

    # Formatting
    ax.set_xlabel('Raman Shift (cm⁻¹)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Normalized Intensity (offset)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlim(shift_range)

    # Invert x-axis (high to low wavenumber - spectroscopy convention)
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

    logger.info("Stacked Raman plot created successfully")
    return fig


def create_single_raman_with_regions(
    frequencies: np.ndarray,
    activities: np.ndarray,
    shift_range: Tuple[float, float] = (0, 4000),
    fwhm: float = 17.0,
    figsize: Tuple[float, float] = (12, 6),
    show_regions: bool = True,
    show_labels: bool = True,
    use_temperature_correction: bool = True,
    title: str = "Raman Spectrum",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create publication-quality Raman spectrum with regional annotations.

    Args:
        frequencies: Raman shift values (cm⁻¹)
        activities: Raman activities
        shift_range: (min, max) range for plot
        fwhm: Gaussian broadening width
        figsize: Figure size
        show_regions: Show boundary lines
        show_labels: Show region labels
        use_temperature_correction: Apply physical intensity formula
        title: Plot title
        save_path: Path to save figure

    Returns:
        matplotlib Figure object
    """
    # Generate grid
    shifts = np.linspace(shift_range[0], shift_range[1], 8000)

    # Apply temperature correction if requested
    if use_temperature_correction:
        intensities = raman_activity_to_intensity(frequencies, activities)
    else:
        intensities = activities

    # Build continuous spectrum
    peaks = list(zip(frequencies, intensities))
    spectrum = gaussian_broadening(shifts, peaks, fwhm)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot continuous spectrum
    ax.plot(shifts, spectrum, color='black', linewidth=1.5)
    ax.fill_between(shifts, spectrum, alpha=0.2, color='blue')

    # Plot stick spectrum below
    for freq, inten in zip(frequencies, intensities):
        if shift_range[0] <= freq <= shift_range[1]:
            max_spectrum = np.max(spectrum)
            stick_height = max_spectrum * 0.1
            ax.vlines(freq, -stick_height, 0, color='red', alpha=0.7, linewidth=1.2)

    # Add regional boundaries
    if show_regions:
        boundaries = get_raman_region_boundaries()
        y_max = np.max(spectrum) * 1.1
        for boundary in boundaries:
            if shift_range[0] <= boundary <= shift_range[1]:
                ax.axvline(boundary, color='gray', linestyle='--',
                          alpha=0.5, linewidth=0.8)

    # Add region labels
    if show_labels:
        region_labels = get_raman_region_labels()
        y_pos = np.max(spectrum) * 0.95
        for shift, label in region_labels:
            if shift_range[0] <= shift <= shift_range[1]:
                ax.text(shift, y_pos, label, fontsize=9, ha='center',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                               edgecolor='gray', alpha=0.8))

    # Formatting
    ax.set_xlabel('Raman Shift (cm⁻¹)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Intensity (arb. units)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlim(shift_range)
    ax.set_ylim(bottom=-np.max(spectrum) * 0.15)

    # Invert x-axis
    ax.invert_xaxis()

    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig


# Export main functions
__all__ = ['create_stacked_raman_plot', 'create_single_raman_with_regions']
