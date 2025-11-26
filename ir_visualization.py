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


def assign_functional_groups(
    frequencies: np.ndarray,
    intensities: np.ndarray,
    intensity_threshold: float = 20.0
) -> List[Dict[str, any]]:
    """
    Assign functional groups based on IR peak positions.

    Args:
        frequencies: Peak frequencies (cm⁻¹)
        intensities: Peak intensities
        intensity_threshold: Minimum intensity to consider

    Returns:
        List of dicts with 'frequency', 'intensity', 'functional_group', 'confidence'
    """
    # Comprehensive functional group assignment table
    # Format: (min_freq, max_freq, functional_group, confidence_level)
    assignment_table = [
        # O-H and N-H stretches
        (3200, 3700, 'O-H stretch (alcohol/phenol)', 'high'),
        (3100, 3500, 'N-H stretch (amine/amide)', 'medium'),
        (2500, 3300, 'O-H stretch (carboxylic acid, broad)', 'high'),

        # C-H stretches
        (3000, 3100, 'C-H stretch (aromatic)', 'high'),
        (2850, 3000, 'C-H stretch (aliphatic)', 'high'),
        (2720, 2820, 'C-H stretch (aldehyde)', 'medium'),

        # Triple bonds
        (2100, 2260, 'C≡C stretch (alkyne)', 'high'),
        (2210, 2260, 'C≡N stretch (nitrile)', 'high'),

        # C=O stretches
        (1800, 1850, 'C=O stretch (acid chloride)', 'high'),
        (1735, 1750, 'C=O stretch (ester)', 'high'),
        (1700, 1725, 'C=O stretch (carboxylic acid)', 'high'),
        (1680, 1720, 'C=O stretch (ketone)', 'high'),
        (1670, 1690, 'C=O stretch (amide I)', 'high'),
        (1650, 1680, 'C=C stretch (alkene)', 'medium'),

        # Aromatic rings
        (1450, 1600, 'C=C stretch (aromatic ring)', 'high'),

        # C-H bending
        (1340, 1470, 'C-H bending (alkanes)', 'medium'),
        (1000, 1300, 'C-O stretch (alcohols/ethers/esters)', 'medium'),

        # Fingerprint region
        (690, 900, 'C-H out-of-plane bending (aromatic)', 'low'),
        (650, 1000, 'Fingerprint region (varies)', 'low'),

        # Halogen stretches
        (500, 800, 'C-Cl stretch', 'low'),
        (400, 700, 'C-Br stretch', 'low'),

        # Other important groups
        (1500, 1560, 'N-O stretch (nitro)', 'high'),
        (1020, 1220, 'S=O stretch (sulfone)', 'medium'),
    ]

    assignments = []

    for freq, intensity in zip(frequencies, intensities):
        if intensity < intensity_threshold:
            continue

        # Find matching functional groups
        matches = []
        for min_freq, max_freq, group, confidence in assignment_table:
            if min_freq <= freq <= max_freq:
                matches.append({
                    'functional_group': group,
                    'confidence': confidence
                })

        if matches:
            # Use highest confidence match
            confidence_order = {'high': 3, 'medium': 2, 'low': 1}
            best_match = max(matches, key=lambda x: confidence_order[x['confidence']])

            assignments.append({
                'frequency': freq,
                'intensity': intensity,
                'functional_group': best_match['functional_group'],
                'confidence': best_match['confidence']
            })
        else:
            # No clear assignment
            assignments.append({
                'frequency': freq,
                'intensity': intensity,
                'functional_group': 'Unassigned',
                'confidence': 'none'
            })

    logger.info(f"Assigned {len(assignments)} peaks to functional groups")
    return assignments


def create_ir_with_assignment_table(
    frequencies: np.ndarray,
    intensities: np.ndarray,
    wavenumber_range: Tuple[float, float] = (400, 4000),
    fwhm: float = 20.0,
    scale_factor: float = 1.0,
    shift_cm: float = 0.0,
    concentration_mol_l: float = 0.01,
    path_length_cm: float = 1.0,
    intensity_threshold: float = 20.0,
    figsize: Tuple[float, float] = (14, 10),
    title: str = "IR Spectrum with Functional Group Assignment",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create IR spectrum with functional group assignment table.

    Args:
        frequencies: Peak frequencies (cm⁻¹)
        intensities: Peak intensities (km/mol)
        wavenumber_range: Range of wavenumbers to plot
        fwhm: Full width at half maximum for Gaussian broadening
        scale_factor: DFT frequency scaling factor (typically 0.96-0.97)
        shift_cm: Additional frequency shift (cm⁻¹)
        concentration_mol_l: Concentration for Beer-Lambert law (mol/L)
        path_length_cm: Path length for Beer-Lambert law (cm)
        intensity_threshold: Minimum intensity for assignment
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info("Creating IR spectrum with functional group assignment table")

    # Apply frequency scaling
    frequencies = apply_frequency_scaling(frequencies, scale_factor, shift_cm)

    # Generate broadened spectrum
    x = np.linspace(wavenumber_range[0], wavenumber_range[1], 2000)
    peaks = list(zip(frequencies, intensities))
    spectrum = gaussian_broadening(x, peaks, fwhm)

    # Calculate absorbance
    absorbance, transmittance = calculate_absorbance_transmittance(
        spectrum, concentration_mol_l, path_length_cm
    )

    # Get functional group assignments
    assignments = assign_functional_groups(frequencies, intensities, intensity_threshold)

    # Create figure with spectrum and table
    fig = plt.figure(figsize=figsize)

    # Spectrum plot (top 60%)
    ax_spec = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    ax_spec.plot(x, absorbance, 'b-', linewidth=1.5, label='Absorbance')
    ax_spec.invert_xaxis()
    ax_spec.set_ylabel('Absorbance', fontsize=11, fontweight='bold')
    ax_spec.set_title(title, fontsize=12, fontweight='bold')
    ax_spec.spines['top'].set_visible(False)
    ax_spec.spines['right'].set_visible(False)
    ax_spec.grid(True, alpha=0.3)

    # Mark peaks in spectrum
    for assignment in assignments:
        freq = assignment['frequency']
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            ax_spec.axvline(freq, color='red', linestyle=':', alpha=0.5, linewidth=1)

    # Assignment table (bottom 40%)
    ax_table = plt.subplot2grid((5, 1), (3, 0), rowspan=2)
    ax_table.axis('off')

    # Prepare table data
    table_data = [['Frequency\n(cm⁻¹)', 'Intensity\n(km/mol)', 'Functional Group', 'Confidence']]

    # Sort by frequency (descending)
    assignments_sorted = sorted(assignments, key=lambda x: x['frequency'], reverse=True)

    for assignment in assignments_sorted[:15]:  # Show top 15
        table_data.append([
            f"{assignment['frequency']:.1f}",
            f"{assignment['intensity']:.1f}",
            assignment['functional_group'],
            assignment['confidence'].capitalize()
        ])

    # Create table
    table = ax_table.table(cellText=table_data, cellLoc='left',
                          loc='center', colWidths=[0.15, 0.15, 0.50, 0.20])

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Style header row
    for i in range(4):
        cell = table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white')

    # Color-code confidence levels
    confidence_colors = {
        'High': '#C6EFCE',
        'Medium': '#FFEB9C',
        'Low': '#FFC7CE',
        'None': '#E0E0E0'
    }

    for i, assignment in enumerate(assignments_sorted[:15], start=1):
        conf_level = assignment['confidence'].capitalize()
        color = confidence_colors.get(conf_level, 'white')
        table[(i, 3)].set_facecolor(color)

    ax_spec.set_xlabel('Wavenumber (cm⁻¹)', fontsize=11, fontweight='bold')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


# Export main functions
__all__ = [
    'create_stacked_ir_plot',
    'create_dual_ir_plot',
    'assign_functional_groups',
    'create_ir_with_assignment_table'
]
