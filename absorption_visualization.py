"""
Absorption Spectrum Comparison and Visualization

Creates publication-quality absorption spectra comparisons including:
- Multi-method comparison (VG, AH, AHAS, etc.)
- Electric vs Velocity dipole comparison
- Gaussian broadened spectra from discrete transitions
- Oscillator strength vs wavelength plots
- Side-by-side and overlay comparisons

Implements Item 19: Absorption Spectrum Comparison
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Union
import logging

# Import local utilities
from visualization_utils import (
    gaussian_broadening,
    normalize_spectrum,
    get_color_palette,
    calculate_stack_offsets
)

logger = logging.getLogger(__name__)


def create_broadened_absorption_spectrum(
    wavelengths: np.ndarray,
    oscillator_strengths: np.ndarray,
    wavelength_range: Tuple[float, float] = (200, 800),
    fwhm: float = 20.0,
    num_points: int = 2000
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create Gaussian-broadened absorption spectrum from discrete transitions.

    Args:
        wavelengths: Transition wavelengths (nm)
        oscillator_strengths: Oscillator strengths (dimensionless)
        wavelength_range: Range of wavelengths to plot (nm)
        fwhm: Full width at half maximum for Gaussian broadening (nm)
        num_points: Number of points in output spectrum

    Returns:
        (wavelength_grid, broadened_spectrum)
    """
    # Create wavelength grid
    wl_grid = np.linspace(wavelength_range[0], wavelength_range[1], num_points)

    # Create peaks list
    peaks = list(zip(wavelengths, oscillator_strengths))

    # Apply Gaussian broadening
    spectrum = gaussian_broadening(wl_grid, peaks, fwhm)

    logger.debug(f"Created broadened spectrum from {len(wavelengths)} transitions")

    return wl_grid, spectrum


def extract_absorption_data(transitions: List[Dict]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract wavelength, energy, and oscillator strength from transition list.

    Args:
        transitions: List of transition dicts with 'wavelength_nm', 'energy_ev',
                    and 'fosc_d2' or 'fosc_p2' keys

    Returns:
        (wavelengths_nm, energies_ev, oscillator_strengths)
    """
    wavelengths = []
    energies = []
    oscillator_strengths = []

    for trans in transitions:
        wl = trans.get('wavelength_nm')
        en = trans.get('energy_ev')
        # Try both electric and velocity dipole keys
        fosc = trans.get('fosc_d2', trans.get('fosc_p2', 0.0))

        if wl is not None and en is not None:
            wavelengths.append(wl)
            energies.append(en)
            oscillator_strengths.append(fosc)

    return np.array(wavelengths), np.array(energies), np.array(oscillator_strengths)


def create_absorption_comparison(
    datasets: List[Dict],
    labels: List[str],
    wavelength_range: Tuple[float, float] = (200, 800),
    fwhm: float = 20.0,
    figsize: Tuple[float, float] = (14, 10),
    colors: Optional[List[str]] = None,
    layout: str = 'overlay',  # 'overlay' or 'stacked'
    show_sticks: bool = True,
    normalize: bool = True,
    title: str = "Absorption Spectrum Comparison",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create multi-dataset absorption spectrum comparison.

    Args:
        datasets: List of dicts with 'transitions' key containing transition data
        labels: Dataset labels (e.g., ['VG', 'AH', 'AHAS'] or ['Electric', 'Velocity'])
        wavelength_range: Range of wavelengths to plot (nm)
        fwhm: FWHM for Gaussian broadening (nm)
        figsize: Figure size
        colors: Custom color list
        layout: 'overlay' (all on same axes) or 'stacked' (vertical offset)
        show_sticks: Show stick spectrum for transitions
        normalize: Normalize spectra to same maximum
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info(f"Creating absorption comparison for {len(datasets)} datasets")

    if len(datasets) != len(labels):
        raise ValueError(f"Number of datasets ({len(datasets)}) must match labels ({len(labels)})")

    # Get colors
    if colors is None:
        colors = get_color_palette(len(datasets))

    # Process all datasets
    spectra_data = []
    for i, (dataset, label) in enumerate(zip(datasets, labels)):
        try:
            transitions = dataset.get('transitions', [])
            if not transitions:
                logger.warning(f"No transitions found in dataset {i}: {label}")
                continue

            # Extract data
            wavelengths, energies, fosc = extract_absorption_data(transitions)

            if len(wavelengths) == 0:
                logger.warning(f"No valid transitions in dataset {i}: {label}")
                continue

            # Create broadened spectrum
            wl_grid, spectrum = create_broadened_absorption_spectrum(
                wavelengths, fosc, wavelength_range, fwhm
            )

            if normalize:
                spectrum = spectrum / np.max(spectrum) if np.max(spectrum) > 0 else spectrum

            spectra_data.append({
                'label': label,
                'wavelengths': wavelengths,
                'energies': energies,
                'fosc': fosc,
                'wl_grid': wl_grid,
                'spectrum': spectrum,
                'color': colors[i]
            })

            logger.debug(f"Processed dataset {i}: {label}, {len(wavelengths)} transitions")

        except Exception as e:
            logger.error(f"Error processing dataset {i} ({label}): {e}")
            raise

    if not spectra_data:
        raise ValueError("No valid spectra data to plot")

    # Create figure
    if layout == 'stacked':
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1])
    else:
        fig, ax1 = plt.subplots(figsize=figsize)
        ax2 = None

    # Plot broadened spectra
    if layout == 'overlay':
        # All spectra on same axes
        for data in spectra_data:
            ax1.plot(data['wl_grid'], data['spectrum'],
                    color=data['color'], linewidth=2.0, label=data['label'])

            # Add stick spectrum
            if show_sticks:
                for wl, fosc in zip(data['wavelengths'], data['fosc']):
                    if wavelength_range[0] <= wl <= wavelength_range[1]:
                        ax1.vlines(wl, 0, fosc * 0.8,
                                 color=data['color'], alpha=0.3, linewidth=1.0)

    elif layout == 'stacked':
        # Calculate vertical offsets
        spectra_list = [data['spectrum'] for data in spectra_data]
        y_offsets = calculate_stack_offsets(spectra_list, y_space=0.1)

        for data, offset in zip(spectra_data, y_offsets):
            offset_spectrum = data['spectrum'] + offset
            ax1.plot(data['wl_grid'], offset_spectrum,
                    color=data['color'], linewidth=2.0, label=data['label'])

            # Add stick spectrum
            if show_sticks:
                for wl, fosc in zip(data['wavelengths'], data['fosc']):
                    if wavelength_range[0] <= wl <= wavelength_range[1]:
                        ax1.vlines(wl, offset, offset + fosc * 0.3,
                                 color=data['color'], alpha=0.3, linewidth=1.0)

    # Formatting for main plot
    ax1.set_xlabel('Wavelength (nm)', fontsize=12, fontweight='bold')
    if normalize:
        ax1.set_ylabel('Normalized Intensity', fontsize=12, fontweight='bold')
    else:
        ax1.set_ylabel('Oscillator Strength', fontsize=12, fontweight='bold')

    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_xlim(wavelength_range)
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add secondary x-axis for energy
    ax1_energy = ax1.twiny()
    # Convert wavelength to energy: E(eV) = 1240/λ(nm)
    energy_min = 1240 / wavelength_range[1]
    energy_max = 1240 / wavelength_range[0]
    ax1_energy.set_xlim(energy_max, energy_min)  # Inverted for increasing wavelength
    ax1_energy.set_xlabel('Energy (eV)', fontsize=11, fontweight='bold')

    # Comparison metrics panel (for overlay layout)
    if layout == 'stacked' and ax2 is not None:
        ax2.axis('off')

        # Create comparison table
        table_data = [['Method', 'Max λ (nm)', 'Max E (eV)', 'Total f_osc']]

        for data in spectra_data:
            # Find maximum
            max_idx = np.argmax(data['spectrum'])
            max_wl = data['wl_grid'][max_idx]
            max_e = 1240 / max_wl if max_wl > 0 else 0

            # Total oscillator strength
            total_fosc = np.sum(data['fosc'])

            table_data.append([
                data['label'],
                f"{max_wl:.1f}",
                f"{max_e:.2f}",
                f"{total_fosc:.4f}"
            ])

        table = ax2.table(cellText=table_data, cellLoc='center',
                         loc='center', colWidths=[0.3, 0.25, 0.25, 0.20])

        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)

        # Style header
        for i in range(4):
            cell = table[(0, i)]
            cell.set_facecolor('#4472C4')
            cell.set_text_props(weight='bold', color='white')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


def create_electric_vs_velocity_comparison(
    electric_transitions: List[Dict],
    velocity_transitions: List[Dict],
    wavelength_range: Tuple[float, float] = (200, 800),
    fwhm: float = 20.0,
    figsize: Tuple[float, float] = (14, 10),
    title: str = "Electric vs Velocity Dipole Comparison",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create specialized comparison between electric and velocity dipole spectra.

    Args:
        electric_transitions: Electric dipole transition list
        velocity_transitions: Velocity dipole transition list
        wavelength_range: Range to plot (nm)
        fwhm: FWHM for broadening (nm)
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info("Creating electric vs velocity dipole comparison")

    datasets = [
        {'transitions': electric_transitions},
        {'transitions': velocity_transitions}
    ]

    labels = ['Electric Dipole', 'Velocity Dipole']

    fig = create_absorption_comparison(
        datasets, labels,
        wavelength_range=wavelength_range,
        fwhm=fwhm,
        figsize=figsize,
        colors=['blue', 'red'],
        layout='overlay',
        show_sticks=True,
        normalize=True,
        title=title,
        save_path=save_path
    )

    return fig


def create_transition_table(
    transitions: List[Dict],
    num_transitions: int = 10,
    figsize: Tuple[float, float] = (10, 6),
    title: str = "Major Transitions",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create table showing major transitions.

    Args:
        transitions: List of transition dicts
        num_transitions: Number of transitions to show
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info(f"Creating transition table for {len(transitions)} transitions")

    # Extract data
    wavelengths, energies, fosc = extract_absorption_data(transitions)

    # Sort by oscillator strength (descending)
    sorted_indices = np.argsort(fosc)[::-1][:num_transitions]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')

    # Prepare table data
    table_data = [['#', 'λ (nm)', 'E (eV)', 'f_osc', 'Intensity']]

    for i, idx in enumerate(sorted_indices, start=1):
        # Calculate relative intensity
        intensity = fosc[idx] / fosc[sorted_indices[0]] * 100

        table_data.append([
            str(i),
            f"{wavelengths[idx]:.1f}",
            f"{energies[idx]:.3f}",
            f"{fosc[idx]:.4f}",
            f"{intensity:.1f}%"
        ])

    # Create table
    table = ax.table(cellText=table_data, cellLoc='center',
                    loc='center', colWidths=[0.10, 0.20, 0.20, 0.25, 0.25])

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.5)

    # Style header row
    for i in range(5):
        cell = table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white')

    # Color-code by intensity
    for i in range(1, len(table_data)):
        intensity_val = float(table_data[i][4].rstrip('%'))
        if intensity_val > 80:
            color = '#C6EFCE'  # Green - strong
        elif intensity_val > 50:
            color = '#FFEB9C'  # Yellow - medium
        elif intensity_val > 20:
            color = '#FFC7CE'  # Orange - weak
        else:
            color = '#E0E0E0'  # Gray - very weak

        table[(i, 4)].set_facecolor(color)

    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


# Export main functions
__all__ = [
    'create_broadened_absorption_spectrum',
    'extract_absorption_data',
    'create_absorption_comparison',
    'create_electric_vs_velocity_comparison',
    'create_transition_table'
]
