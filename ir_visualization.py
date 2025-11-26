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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
    normalize_spectrum,
    find_optimal_label_positions
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
) -> go.Figure:
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
        Plotly Figure object
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
    fig = go.Figure()

    # Plot each spectrum with offset
    for i, (spectrum, offset, label, color) in enumerate(zip(spectra, y_offsets, labels, colors)):
        offset_spectrum = spectrum + offset

        # Add main spectrum line
        fig.add_trace(go.Scatter(
            x=x,
            y=offset_spectrum,
            mode='lines',
            name=label,
            line=dict(color=color, width=1.5),
            hovertemplate='Wavenumber: %{x:.1f} cm⁻¹<br>Intensity: %{y:.3f}<extra></extra>'
        ))

        # Add stick spectrum
        if 'frequencies' in datasets[i]:
            freqs = datasets[i]['frequencies']
            if scale_factor != 1.0 or shift_cm != 0.0:
                freqs = apply_frequency_scaling(np.array(freqs), scale_factor, shift_cm)

            stick_x = []
            stick_y = []
            for freq in freqs:
                if wavenumber_range[0] <= freq <= wavenumber_range[1]:
                    if plot_type == 'transmittance':
                        stick_height = 5.0
                        stick_base = offset + np.max(spectrum[np.abs(x - freq) < 50])
                    else:
                        stick_height = np.max(spectrum) * 0.05
                        stick_base = offset

                    stick_x.extend([freq, freq, None])
                    stick_y.extend([stick_base, stick_base + stick_height, None])

            if stick_x:
                fig.add_trace(go.Scatter(
                    x=stick_x,
                    y=stick_y,
                    mode='lines',
                    line=dict(color=color, width=0.8),
                    opacity=0.6,
                    showlegend=False,
                    hoverinfo='skip'
                ))

    # Add regional boundaries
    if show_regions:
        boundaries = get_ir_region_boundaries()
        for boundary in boundaries:
            if wavenumber_range[0] <= boundary <= wavenumber_range[1]:
                fig.add_vline(
                    x=boundary,
                    line_dash="dash",
                    line_color="gray",
                    opacity=0.5,
                    line_width=0.8
                )

    # Auto-generate title if not provided
    if title is None:
        correction_str = f" (scaled: {scale_factor}×)" if scale_factor != 1.0 else ""
        title = f"Stacked IR Spectra - {plot_type.capitalize()}{correction_str}"

    # Formatting
    ylabel = 'Transmittance (%, offset)' if plot_type == 'transmittance' else 'Absorbance (offset)'

    fig.update_layout(
        title=dict(text=title, font=dict(size=14, family='Arial Black')),
        xaxis_title=dict(text='Wavenumber (cm⁻¹)', font=dict(size=12, family='Arial Black')),
        yaxis_title=dict(text=ylabel, font=dict(size=12, family='Arial Black')),
        template='plotly_white',
        hovermode='closest',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=True,
        legend=dict(
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02
        )
    )

    # Invert x-axis (high to low wavenumber - spectroscopy convention)
    fig.update_xaxes(autorange='reversed', range=[wavenumber_range[1], wavenumber_range[0]])

    # Hide y-axis ticks (arbitrary offsets)
    fig.update_yaxes(showticklabels=False)

    # Save if requested
    if save_path:
        logger.info(f"Saving figure to {save_path}")
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=1400, height=800, scale=3)

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
) -> go.Figure:
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
        Plotly Figure object
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
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=(f"{title} - Absorbance", f"{title} - Transmittance"),
        vertical_spacing=0.12,
        shared_xaxes=True
    )

    # Plot absorbance (top)
    max_abs = np.max(absorbance)

    fig.add_trace(go.Scatter(
        x=x, y=absorbance,
        mode='lines',
        line=dict(color='blue', width=1.5),
        fill='tozeroy',
        fillcolor='rgba(0, 0, 255, 0.2)',
        name='Absorbance',
        hovertemplate='Wavenumber: %{x:.1f} cm⁻¹<br>Absorbance: %{y:.3f}<extra></extra>'
    ), row=1, col=1)

    # Stick spectrum for absorbance
    stick_x_abs = []
    stick_y_abs = []
    for freq, inten in zip(frequencies, intensities):
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            stick_x_abs.extend([freq, freq, None])
            stick_y_abs.extend([max_abs * 1.05, max_abs * 1.15, None])

    if stick_x_abs:
        fig.add_trace(go.Scatter(
            x=stick_x_abs, y=stick_y_abs,
            mode='lines',
            line=dict(color='red', width=1.0),
            opacity=0.6,
            showlegend=False,
            hoverinfo='skip'
        ), row=1, col=1)

    # Plot transmittance (bottom)
    fig.add_trace(go.Scatter(
        x=x, y=transmittance,
        mode='lines',
        line=dict(color='green', width=1.5),
        fill='tozeroy',
        fillcolor='rgba(0, 255, 0, 0.2)',
        name='Transmittance',
        hovertemplate='Wavenumber: %{x:.1f} cm⁻¹<br>Transmittance: %{y:.2f}%<extra></extra>'
    ), row=2, col=1)

    # Stick spectrum for transmittance
    stick_x_trans = []
    stick_y_trans = []
    for freq, inten in zip(frequencies, intensities):
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            stick_x_trans.extend([freq, freq, None])
            stick_y_trans.extend([100, 105, None])

    if stick_x_trans:
        fig.add_trace(go.Scatter(
            x=stick_x_trans, y=stick_y_trans,
            mode='lines',
            line=dict(color='red', width=1.0),
            opacity=0.6,
            showlegend=False,
            hoverinfo='skip'
        ), row=2, col=1)

    # Add regional boundaries to both
    if show_regions:
        boundaries = get_ir_region_boundaries()
        for boundary in boundaries:
            if wavenumber_range[0] <= boundary <= wavenumber_range[1]:
                fig.add_vline(
                    x=boundary,
                    line_dash="dash",
                    line_color="gray",
                    opacity=0.5,
                    line_width=0.8,
                    row='all'
                )

    # Add functional group labels
    if show_labels:
        label_data = get_ir_region_labels()
        for freq, label_text in label_data:
            if wavenumber_range[0] <= freq <= wavenumber_range[1]:
                # Add to absorbance plot
                fig.add_annotation(
                    x=freq,
                    y=max_abs * 0.95,
                    text=label_text,
                    showarrow=False,
                    bgcolor='white',
                    bordercolor='gray',
                    borderwidth=1,
                    borderpad=2,
                    font=dict(size=8),
                    row=1, col=1
                )
                # Add to transmittance plot
                fig.add_annotation(
                    x=freq,
                    y=95,
                    text=label_text,
                    showarrow=False,
                    bgcolor='white',
                    bordercolor='gray',
                    borderwidth=1,
                    borderpad=2,
                    font=dict(size=8),
                    row=2, col=1
                )

    # Update layout
    fig.update_layout(
        template='plotly_white',
        hovermode='closest',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=False
    )

    # Update axes
    fig.update_xaxes(
        title_text='Wavenumber (cm⁻¹)',
        title_font=dict(size=11, family='Arial Black'),
        autorange='reversed',
        row=2, col=1
    )
    fig.update_yaxes(
        title_text='Absorbance',
        title_font=dict(size=11, family='Arial Black'),
        range=[0, max_abs * 1.2],
        row=1, col=1
    )
    fig.update_yaxes(
        title_text='Transmittance (%)',
        title_font=dict(size=11, family='Arial Black'),
        range=[0, 110],
        row=2, col=1
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=1200, height=1000, scale=3)

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
) -> go.Figure:
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
        Plotly Figure object
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
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.6, 0.4],
        subplot_titles=(title, "Functional Group Assignments"),
        vertical_spacing=0.15,
        specs=[[{"type": "scatter"}], [{"type": "table"}]]
    )

    # Spectrum plot (top)
    fig.add_trace(go.Scatter(
        x=x,
        y=absorbance,
        mode='lines',
        line=dict(color='blue', width=1.5),
        name='Absorbance',
        hovertemplate='Wavenumber: %{x:.1f} cm⁻¹<br>Absorbance: %{y:.3f}<extra></extra>'
    ), row=1, col=1)

    # Mark peaks in spectrum
    for assignment in assignments:
        freq = assignment['frequency']
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            fig.add_vline(
                x=freq,
                line_dash="dot",
                line_color="red",
                opacity=0.5,
                line_width=1,
                row=1
            )

    # Assignment table (bottom)
    # Sort by frequency (descending)
    assignments_sorted = sorted(assignments, key=lambda x: x['frequency'], reverse=True)[:15]

    # Prepare table data
    table_header = ['Frequency<br>(cm⁻¹)', 'Intensity<br>(km/mol)', 'Functional Group', 'Confidence']
    table_cells = [
        [f"{a['frequency']:.1f}" for a in assignments_sorted],
        [f"{a['intensity']:.1f}" for a in assignments_sorted],
        [a['functional_group'] for a in assignments_sorted],
        [a['confidence'].capitalize() for a in assignments_sorted]
    ]

    # Color-code confidence levels
    confidence_colors = {
        'High': '#C6EFCE',
        'Medium': '#FFEB9C',
        'Low': '#FFC7CE',
        'None': '#E0E0E0'
    }
    cell_colors = [
        ['white'] * len(assignments_sorted),
        ['white'] * len(assignments_sorted),
        ['white'] * len(assignments_sorted),
        [confidence_colors.get(a['confidence'].capitalize(), 'white') for a in assignments_sorted]
    ]

    fig.add_trace(go.Table(
        header=dict(
            values=table_header,
            fill_color='#4472C4',
            align='left',
            font=dict(color='white', size=10, family='Arial Black')
        ),
        cells=dict(
            values=table_cells,
            fill_color=cell_colors,
            align='left',
            font=dict(size=9),
            height=25
        )
    ), row=2, col=1)

    # Update layout
    fig.update_layout(
        template='plotly_white',
        hovermode='closest',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=False
    )

    fig.update_xaxes(
        title_text='Wavenumber (cm⁻¹)',
        title_font=dict(size=11, family='Arial Black'),
        autorange='reversed',
        row=1, col=1
    )
    fig.update_yaxes(
        title_text='Absorbance',
        title_font=dict(size=11, family='Arial Black'),
        row=1, col=1
    )

    if save_path:
        logger.info(f"Saved figure to {save_path}")
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=1400, height=1000, scale=3)

    return fig


# Export main functions
__all__ = [
    'create_stacked_ir_plot',
    'create_dual_ir_plot',
    'assign_functional_groups',
    'create_ir_with_assignment_table'
]
