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
) -> go.Figure:
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
        plotly Figure object
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
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.75, 0.25],
            vertical_spacing=0.1,
            specs=[[{"secondary_y": False}], [{"type": "table"}]]
        )
    else:
        fig = go.Figure()

    # Plot broadened spectra
    if layout == 'overlay':
        # All spectra on same axes
        for data in spectra_data:
            fig.add_trace(go.Scatter(
                x=data['wl_grid'], y=data['spectrum'],
                mode='lines', name=data['label'],
                line=dict(color=data['color'], width=2.0)
            ))

            # Add stick spectrum
            if show_sticks:
                for wl, fosc in zip(data['wavelengths'], data['fosc']):
                    if wavelength_range[0] <= wl <= wavelength_range[1]:
                        fig.add_trace(go.Scatter(
                            x=[wl, wl], y=[0, fosc * 0.8],
                            mode='lines', line=dict(color=data['color'], width=1.0),
                            opacity=0.3, showlegend=False, hoverinfo='skip'
                        ))

    elif layout == 'stacked':
        # Calculate vertical offsets
        spectra_list = [data['spectrum'] for data in spectra_data]
        y_offsets = calculate_stack_offsets(spectra_list, y_space=0.1)

        for data, offset in zip(spectra_data, y_offsets):
            offset_spectrum = data['spectrum'] + offset
            fig.add_trace(go.Scatter(
                x=data['wl_grid'], y=offset_spectrum,
                mode='lines', name=data['label'],
                line=dict(color=data['color'], width=2.0)
            ), row=1, col=1)

            # Add stick spectrum
            if show_sticks:
                for wl, fosc in zip(data['wavelengths'], data['fosc']):
                    if wavelength_range[0] <= wl <= wavelength_range[1]:
                        fig.add_trace(go.Scatter(
                            x=[wl, wl], y=[offset, offset + fosc * 0.3],
                            mode='lines', line=dict(color=data['color'], width=1.0),
                            opacity=0.3, showlegend=False, hoverinfo='skip'
                        ), row=1, col=1)

    # Formatting for main plot
    ylabel = 'Normalized Intensity' if normalize else 'Oscillator Strength'

    # Add secondary x-axis for energy
    energy_min = 1240 / wavelength_range[1]
    energy_max = 1240 / wavelength_range[0]

    fig.update_layout(
        xaxis_title='Wavelength (nm)',
        yaxis_title=ylabel,
        title=title,
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=True,
        xaxis=dict(range=wavelength_range, showgrid=True),
        xaxis2=dict(
            overlaying='x',
            side='top',
            range=[energy_max, energy_min],
            title='Energy (eV)'
        )
    )

    # Comparison metrics panel (for stacked layout)
    if layout == 'stacked':
        # Create comparison table
        table_data = [['Method', 'Max λ (nm)', 'Max E (eV)', 'Total f_osc']]

        for data in spectra_data:
            # Find maximum
            max_idx = np.argmax(data['spectrum'])
            max_wl = data['wl_grid'][max_idx]
            max_e = 1240 / max_wl if max_wl > 0 else 0

            # Total oscillator strength
            total_fosc = np.sum(data['fosc'])

            table_data[0].append(data['label'])
            table_data.append([
                'Max λ (nm)', f"{max_wl:.1f}", '', ''
            ])

        # Add table to subplot
        header = table_data[0]
        cells = [row[1:] for row in table_data[1:]]

        fig.add_trace(go.Table(
            header=dict(values=header,
                       fill_color='#4472C4',
                       font=dict(color='white', size=11),
                       align='center'),
            cells=dict(values=[[data['label'],
                               f"{data['wl_grid'][np.argmax(data['spectrum'])]:.1f}",
                               f"{1240/data['wl_grid'][np.argmax(data['spectrum'])]:.2f}",
                               f"{np.sum(data['fosc']):.4f}"] for data in spectra_data],
                      align='center')
        ), row=2, col=1)

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
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
) -> go.Figure:
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
        plotly Figure object
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
) -> go.Figure:
    """
    Create table showing major transitions.

    Args:
        transitions: List of transition dicts
        num_transitions: Number of transitions to show
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        plotly Figure object
    """
    logger.info(f"Creating transition table for {len(transitions)} transitions")

    # Extract data
    wavelengths, energies, fosc = extract_absorption_data(transitions)

    # Sort by oscillator strength (descending)
    sorted_indices = np.argsort(fosc)[::-1][:num_transitions]

    # Prepare table data
    header = ['#', 'λ (nm)', 'E (eV)', 'f_osc', 'Intensity']
    cells_data = [[], [], [], [], []]

    for i, idx in enumerate(sorted_indices, start=1):
        # Calculate relative intensity
        intensity = fosc[idx] / fosc[sorted_indices[0]] * 100

        cells_data[0].append(str(i))
        cells_data[1].append(f"{wavelengths[idx]:.1f}")
        cells_data[2].append(f"{energies[idx]:.3f}")
        cells_data[3].append(f"{fosc[idx]:.4f}")
        cells_data[4].append(f"{intensity:.1f}%")

    # Create color coding for intensity
    fill_colors = []
    for i, idx in enumerate(sorted_indices):
        intensity_val = fosc[idx] / fosc[sorted_indices[0]] * 100
        if intensity_val > 80:
            color = '#C6EFCE'  # Green - strong
        elif intensity_val > 50:
            color = '#FFEB9C'  # Yellow - medium
        elif intensity_val > 20:
            color = '#FFC7CE'  # Orange - weak
        else:
            color = '#E0E0E0'  # Gray - very weak
        fill_colors.append(['white', 'white', 'white', 'white', color])

    # Create figure
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=header,
            fill_color='#4472C4',
            font=dict(color='white', size=12, family='Arial Bold'),
            align='center',
            height=40
        ),
        cells=dict(
            values=cells_data,
            fill_color=[['white'] * num_transitions] * 4 + [[fill_colors[i][4] for i in range(num_transitions)]],
            align='center',
            height=35,
            font=dict(size=11)
        )
    )])

    fig.update_layout(
        title=title,
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        template='plotly_white'
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
    'create_broadened_absorption_spectrum',
    'extract_absorption_data',
    'create_absorption_comparison',
    'create_electric_vs_velocity_comparison',
    'create_transition_table'
]
