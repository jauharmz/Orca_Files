"""
Fluorescence and Phosphorescence Emission Spectrum Visualization

Creates publication-quality emission spectra including:
- Fluorescence (singlet → singlet, fast emission)
- Phosphorescence (triplet → singlet, slow emission)
- Side-by-side and overlay comparisons
- Stokes shift analysis
- Mirror image absorption vs emission

Implements Item 20: Fluorescence vs Phosphorescence
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


def create_emission_spectrum(
    wavelengths: np.ndarray,
    intensities: np.ndarray,
    wavelength_range: Tuple[float, float] = (300, 800),
    fwhm: float = 20.0,
    num_points: int = 2000
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create Gaussian-broadened emission spectrum from discrete transitions.

    Args:
        wavelengths: Emission wavelengths (nm)
        intensities: Emission intensities (arbitrary units or Einstein A coefficients)
        wavelength_range: Range of wavelengths to plot (nm)
        fwhm: Full width at half maximum for Gaussian broadening (nm)
        num_points: Number of points in output spectrum

    Returns:
        (wavelength_grid, broadened_spectrum)
    """
    # Create wavelength grid
    wl_grid = np.linspace(wavelength_range[0], wavelength_range[1], num_points)

    # Create peaks list
    peaks = list(zip(wavelengths, intensities))

    # Apply Gaussian broadening
    spectrum = gaussian_broadening(wl_grid, peaks, fwhm)

    logger.debug(f"Created emission spectrum from {len(wavelengths)} transitions")

    return wl_grid, spectrum


def calculate_stokes_shift(
    absorption_max: float,
    emission_max: float
) -> Tuple[float, float]:
    """
    Calculate Stokes shift between absorption and emission.

    Args:
        absorption_max: Maximum of absorption spectrum (nm)
        emission_max: Maximum of emission spectrum (nm)

    Returns:
        (stokes_shift_nm, stokes_shift_cm1)
    """
    # Stokes shift in nm (always positive by convention)
    shift_nm = emission_max - absorption_max

    # Convert to wavenumber (cm⁻¹)
    # Δν̃ = 1/λ_abs - 1/λ_em (multiply by 10^7 for nm to cm⁻¹)
    shift_cm1 = (1.0 / absorption_max - 1.0 / emission_max) * 1e7

    logger.debug(f"Stokes shift: {shift_nm:.1f} nm ({shift_cm1:.0f} cm⁻¹)")

    return shift_nm, shift_cm1


def create_fluorescence_phosphorescence_comparison(
    fluorescence_data: Dict,
    phosphorescence_data: Dict,
    wavelength_range: Tuple[float, float] = (300, 800),
    fwhm_fluor: float = 20.0,
    fwhm_phos: float = 30.0,
    figsize: Tuple[float, float] = (14, 10),
    show_sticks: bool = True,
    normalize: bool = True,
    title: str = "Fluorescence vs Phosphorescence",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create fluorescence vs phosphorescence comparison.

    Args:
        fluorescence_data: Dict with 'wavelengths' and 'intensities' for fluorescence
        phosphorescence_data: Dict with 'wavelengths' and 'intensities' for phosphorescence
        wavelength_range: Range of wavelengths to plot (nm)
        fwhm_fluor: FWHM for fluorescence broadening (nm, typically narrower)
        fwhm_phos: FWHM for phosphorescence broadening (nm, typically broader)
        figsize: Figure size
        show_sticks: Show stick spectrum for transitions
        normalize: Normalize spectra to same maximum
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info("Creating fluorescence vs phosphorescence comparison")

    # Extract data
    fluor_wl = np.array(fluorescence_data['wavelengths'])
    fluor_int = np.array(fluorescence_data['intensities'])

    phos_wl = np.array(phosphorescence_data['wavelengths'])
    phos_int = np.array(phosphorescence_data['intensities'])

    # Create broadened spectra
    wl_grid_f, fluor_spectrum = create_emission_spectrum(
        fluor_wl, fluor_int, wavelength_range, fwhm_fluor
    )

    wl_grid_p, phos_spectrum = create_emission_spectrum(
        phos_wl, phos_int, wavelength_range, fwhm_phos
    )

    if normalize:
        fluor_spectrum = fluor_spectrum / np.max(fluor_spectrum) if np.max(fluor_spectrum) > 0 else fluor_spectrum
        phos_spectrum = phos_spectrum / np.max(phos_spectrum) if np.max(phos_spectrum) > 0 else phos_spectrum

    # Find maxima
    fluor_max_idx = np.argmax(fluor_spectrum)
    fluor_max_wl = wl_grid_f[fluor_max_idx]

    phos_max_idx = np.argmax(phos_spectrum)
    phos_max_wl = wl_grid_p[phos_max_idx]

    # Create figure with 2 subplots
    fig = plt.figure(figsize=figsize)

    # Main comparison (top 70%)
    ax1 = plt.subplot2grid((10, 1), (0, 0), rowspan=7)

    # Plot fluorescence
    ax1.plot(wl_grid_f, fluor_spectrum, color='#00FF00', linewidth=2.5,
            label=f'Fluorescence (max: {fluor_max_wl:.1f} nm)', alpha=0.7)
    ax1.fill_between(wl_grid_f, 0, fluor_spectrum, color='#00FF00', alpha=0.2)

    # Plot phosphorescence
    ax1.plot(wl_grid_p, phos_spectrum, color='#FF4500', linewidth=2.5,
            label=f'Phosphorescence (max: {phos_max_wl:.1f} nm)', alpha=0.7)
    ax1.fill_between(wl_grid_p, 0, phos_spectrum, color='#FF4500', alpha=0.2)

    # Add stick spectra
    if show_sticks:
        for wl, intensity in zip(fluor_wl, fluor_int):
            if wavelength_range[0] <= wl <= wavelength_range[1]:
                ax1.vlines(wl, 0, intensity * 0.15 if normalize else intensity,
                          color='green', alpha=0.5, linewidth=1.5)

        for wl, intensity in zip(phos_wl, phos_int):
            if wavelength_range[0] <= wl <= wavelength_range[1]:
                ax1.vlines(wl, 0, intensity * 0.15 if normalize else intensity,
                          color='red', alpha=0.5, linewidth=1.5)

    # Mark maxima
    ax1.axvline(fluor_max_wl, color='green', linestyle='--', alpha=0.5, linewidth=1)
    ax1.axvline(phos_max_wl, color='red', linestyle='--', alpha=0.5, linewidth=1)

    # Add wavelength difference annotation
    red_shift = phos_max_wl - fluor_max_wl
    if red_shift > 0:
        ax1.annotate('', xy=(phos_max_wl, 0.9), xytext=(fluor_max_wl, 0.9),
                    arrowprops=dict(arrowstyle='<->', color='black', lw=2))
        ax1.text((fluor_max_wl + phos_max_wl) / 2, 0.92,
                f'Red shift: {red_shift:.1f} nm',
                fontsize=10, ha='center', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7))

    ax1.set_ylabel('Normalized Intensity' if normalize else 'Intensity',
                  fontsize=12, fontweight='bold')
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(wavelength_range)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add secondary x-axis for energy
    ax1_energy = ax1.twiny()
    energy_min = 1240 / wavelength_range[1]
    energy_max = 1240 / wavelength_range[0]
    ax1_energy.set_xlim(energy_max, energy_min)
    ax1_energy.set_xlabel('Energy (eV)', fontsize=11, fontweight='bold')

    # Information panel (bottom 30%)
    ax2 = plt.subplot2grid((10, 1), (7, 0), rowspan=3)
    ax2.axis('off')

    # Create information table
    table_data = [
        ['Property', 'Fluorescence', 'Phosphorescence'],
        ['Type', 'S₁ → S₀', 'T₁ → S₀'],
        ['Spin', 'Singlet → Singlet', 'Triplet → Singlet'],
        ['Lifetime', '~10⁻⁹ s (ns)', '~10⁻³ to 10² s (ms-s)'],
        ['Max λ (nm)', f'{fluor_max_wl:.1f}', f'{phos_max_wl:.1f}'],
        ['Max E (eV)', f'{1240/fluor_max_wl:.2f}', f'{1240/phos_max_wl:.2f}'],
        ['Rel. Position', 'Higher energy', f'Red-shifted by {red_shift:.1f} nm']
    ]

    table = ax2.table(cellText=table_data, cellLoc='left',
                     loc='center', colWidths=[0.30, 0.35, 0.35])

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Style header row
    for i in range(3):
        cell = table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white')

    # Color code fluorescence and phosphorescence columns
    for i in range(1, len(table_data)):
        table[(i, 1)].set_facecolor('#E8F5E9')  # Light green
        table[(i, 2)].set_facecolor('#FFE0B2')  # Light orange

    ax1.set_xlabel('Wavelength (nm)', fontsize=12, fontweight='bold')

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


def create_absorption_emission_mirror(
    absorption_data: Dict,
    emission_data: Dict,
    wavelength_range: Tuple[float, float] = (250, 700),
    fwhm: float = 20.0,
    figsize: Tuple[float, float] = (12, 8),
    title: str = "Absorption-Emission Mirror Image",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create mirror image plot showing absorption and emission spectra.

    Demonstrates Stokes shift and mirror symmetry relationship.

    Args:
        absorption_data: Dict with 'wavelengths' and 'intensities' for absorption
        emission_data: Dict with 'wavelengths' and 'intensities' for emission
        wavelength_range: Range to plot (nm)
        fwhm: FWHM for broadening (nm)
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info("Creating absorption-emission mirror image")

    # Extract data
    abs_wl = np.array(absorption_data['wavelengths'])
    abs_int = np.array(absorption_data['intensities'])

    em_wl = np.array(emission_data['wavelengths'])
    em_int = np.array(emission_data['intensities'])

    # Create spectra
    wl_grid, abs_spectrum = create_emission_spectrum(
        abs_wl, abs_int, wavelength_range, fwhm
    )

    _, em_spectrum = create_emission_spectrum(
        em_wl, em_int, wavelength_range, fwhm
    )

    # Normalize
    abs_spectrum = abs_spectrum / np.max(abs_spectrum) if np.max(abs_spectrum) > 0 else abs_spectrum
    em_spectrum = em_spectrum / np.max(em_spectrum) if np.max(em_spectrum) > 0 else em_spectrum

    # Find maxima
    abs_max_wl = wl_grid[np.argmax(abs_spectrum)]
    em_max_wl = wl_grid[np.argmax(em_spectrum)]

    # Calculate Stokes shift
    stokes_nm, stokes_cm1 = calculate_stokes_shift(abs_max_wl, em_max_wl)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot absorption (upward)
    ax.plot(wl_grid, abs_spectrum, color='blue', linewidth=2.5,
           label=f'Absorption (max: {abs_max_wl:.1f} nm)')
    ax.fill_between(wl_grid, 0, abs_spectrum, color='blue', alpha=0.2)

    # Plot emission (downward - mirror image)
    ax.plot(wl_grid, -em_spectrum, color='red', linewidth=2.5,
           label=f'Emission (max: {em_max_wl:.1f} nm)')
    ax.fill_between(wl_grid, 0, -em_spectrum, color='red', alpha=0.2)

    # Mark maxima
    ax.axvline(abs_max_wl, color='blue', linestyle='--', alpha=0.5)
    ax.axvline(em_max_wl, color='red', linestyle='--', alpha=0.5)

    # Add Stokes shift annotation
    ax.annotate('', xy=(em_max_wl, 0), xytext=(abs_max_wl, 0),
               arrowprops=dict(arrowstyle='<->', color='black', lw=2.5))
    ax.text((abs_max_wl + em_max_wl) / 2, 0.1,
           f'Stokes Shift\n{stokes_nm:.1f} nm\n({stokes_cm1:.0f} cm⁻¹)',
           fontsize=11, ha='center', fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8))

    ax.axhline(0, color='black', linewidth=1, linestyle='-', alpha=0.5)

    ax.set_xlabel('Wavelength (nm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Normalized Intensity', fontsize=12, fontweight='bold')
    ax.set_title(f'{title}\n(Absorption ↑, Emission ↓)', fontsize=14, fontweight='bold')
    ax.set_xlim(wavelength_range)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='x')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved figure to {save_path}")

    return fig


# Export main functions
__all__ = [
    'create_emission_spectrum',
    'calculate_stokes_shift',
    'create_fluorescence_phosphorescence_comparison',
    'create_absorption_emission_mirror'
]
