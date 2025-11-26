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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
) -> go.Figure:
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
        plotly Figure object
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

    # Create figure with subplots
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.7, 0.3],
        vertical_spacing=0.15,
        specs=[[{"secondary_y": False}], [{"type": "table"}]],
        subplot_titles=(title, "")
    )

    # Plot fluorescence
    fig.add_trace(go.Scatter(
        x=wl_grid_f, y=fluor_spectrum,
        mode='lines', name=f'Fluorescence (max: {fluor_max_wl:.1f} nm)',
        line=dict(color='#00FF00', width=2.5),
        fill='tozeroy', fillcolor='rgba(0, 255, 0, 0.2)'
    ), row=1, col=1)

    # Plot phosphorescence
    fig.add_trace(go.Scatter(
        x=wl_grid_p, y=phos_spectrum,
        mode='lines', name=f'Phosphorescence (max: {phos_max_wl:.1f} nm)',
        line=dict(color='#FF4500', width=2.5),
        fill='tozeroy', fillcolor='rgba(255, 69, 0, 0.2)'
    ), row=1, col=1)

    # Add stick spectra
    if show_sticks:
        for wl, intensity in zip(fluor_wl, fluor_int):
            if wavelength_range[0] <= wl <= wavelength_range[1]:
                normalized_intensity = intensity * 0.15 if normalize else intensity
                fig.add_trace(go.Scatter(
                    x=[wl, wl], y=[0, normalized_intensity],
                    mode='lines', line=dict(color='green', width=1.5),
                    opacity=0.5, showlegend=False, hoverinfo='skip'
                ), row=1, col=1)

        for wl, intensity in zip(phos_wl, phos_int):
            if wavelength_range[0] <= wl <= wavelength_range[1]:
                normalized_intensity = intensity * 0.15 if normalize else intensity
                fig.add_trace(go.Scatter(
                    x=[wl, wl], y=[0, normalized_intensity],
                    mode='lines', line=dict(color='red', width=1.5),
                    opacity=0.5, showlegend=False, hoverinfo='skip'
                ), row=1, col=1)

    # Mark maxima
    fig.add_vline(x=fluor_max_wl, line=dict(color='green', dash='dash', width=1),
                  opacity=0.5, row=1, col=1)
    fig.add_vline(x=phos_max_wl, line=dict(color='red', dash='dash', width=1),
                  opacity=0.5, row=1, col=1)

    # Add wavelength difference annotation
    red_shift = phos_max_wl - fluor_max_wl
    if red_shift > 0:
        fig.add_annotation(
            x=(fluor_max_wl + phos_max_wl) / 2, y=0.92,
            text=f'Red shift: {red_shift:.1f} nm',
            showarrow=False,
            bgcolor='yellow', opacity=0.7,
            row=1, col=1
        )

    # Add secondary x-axis for energy
    energy_min = 1240 / wavelength_range[1]
    energy_max = 1240 / wavelength_range[0]

    # Information panel
    table_data = [
        ['Property', 'Fluorescence', 'Phosphorescence'],
        ['Type', 'S₁ → S₀', 'T₁ → S₀'],
        ['Spin', 'Singlet → Singlet', 'Triplet → Singlet'],
        ['Lifetime', '~10⁻⁹ s (ns)', '~10⁻³ to 10² s (ms-s)'],
        ['Max λ (nm)', f'{fluor_max_wl:.1f}', f'{phos_max_wl:.1f}'],
        ['Max E (eV)', f'{1240/fluor_max_wl:.2f}', f'{1240/phos_max_wl:.2f}'],
        ['Rel. Position', 'Higher energy', f'Red-shifted by {red_shift:.1f} nm']
    ]

    # Add table
    fig.add_trace(go.Table(
        header=dict(
            values=['Property', 'Fluorescence', 'Phosphorescence'],
            fill_color='#4472C4',
            font=dict(color='white', size=11),
            align='left'
        ),
        cells=dict(
            values=[[row[0] for row in table_data[1:]],
                   [row[1] for row in table_data[1:]],
                   [row[2] for row in table_data[1:]]],
            fill_color=[['white'] * 6, ['#E8F5E9'] * 6, ['#FFE0B2'] * 6],
            align='left',
            font=dict(size=10)
        )
    ), row=2, col=1)

    # Update layout
    ylabel = 'Normalized Intensity' if normalize else 'Intensity'

    fig.update_xaxes(title_text='Wavelength (nm)', row=1, col=1, range=wavelength_range, showgrid=True)
    fig.update_yaxes(title_text=ylabel, row=1, col=1, showgrid=True)

    fig.update_layout(
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=True,
        xaxis2=dict(
            overlaying='x',
            side='top',
            range=[energy_max, energy_min],
            title='Energy (eV)'
        )
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
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
) -> go.Figure:
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
        plotly Figure object
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
    fig = go.Figure()

    # Plot absorption (upward)
    fig.add_trace(go.Scatter(
        x=wl_grid, y=abs_spectrum,
        mode='lines', name=f'Absorption (max: {abs_max_wl:.1f} nm)',
        line=dict(color='blue', width=2.5),
        fill='tozeroy', fillcolor='rgba(0, 0, 255, 0.2)'
    ))

    # Plot emission (downward - mirror image)
    fig.add_trace(go.Scatter(
        x=wl_grid, y=-em_spectrum,
        mode='lines', name=f'Emission (max: {em_max_wl:.1f} nm)',
        line=dict(color='red', width=2.5),
        fill='tozeroy', fillcolor='rgba(255, 0, 0, 0.2)'
    ))

    # Mark maxima
    fig.add_vline(x=abs_max_wl, line=dict(color='blue', dash='dash'), opacity=0.5)
    fig.add_vline(x=em_max_wl, line=dict(color='red', dash='dash'), opacity=0.5)

    # Add Stokes shift annotation
    fig.add_annotation(
        x=(abs_max_wl + em_max_wl) / 2, y=0.1,
        text=f'Stokes Shift<br>{stokes_nm:.1f} nm<br>({stokes_cm1:.0f} cm⁻¹)',
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=2.5,
        arrowcolor='black',
        ax=(em_max_wl - abs_max_wl) * 40,
        ay=0,
        bgcolor='yellow',
        opacity=0.8
    )

    fig.add_hline(y=0, line=dict(color='black', width=1), opacity=0.5)

    fig.update_layout(
        xaxis_title='Wavelength (nm)',
        yaxis_title='Normalized Intensity',
        title=f'{title}<br>(Absorption ↑, Emission ↓)',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        xaxis=dict(range=wavelength_range, showgrid=True),
        yaxis=dict(showgrid=True, gridcolor='lightgray', zeroline=True)
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved figure to {save_path}")

    return fig


# Export main functions
__all__ = [
    'create_emission_spectrum',
    'calculate_stokes_shift',
    'create_fluorescence_phosphorescence_comparison',
    'create_absorption_emission_mirror'
]
