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
import plotly.graph_objects as go
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
    normalize_spectrum,
    find_optimal_label_positions
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
    label_position_method: str = 'right',
    use_temperature_correction: bool = True,
    laser_wavelength_nm: float = 532,
    low_freq_cutoff: float = 0.0,
    title: str = "Stacked Raman Spectra",
    save_path: Optional[str] = None
) -> go.Figure:
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
        show_labels: Show dataset labels
        label_position_method: 'right' (edge), 'peak' (at max), or 'flat' (at flattest region)
        use_temperature_correction: Apply temperature-dependent intensity
        laser_wavelength_nm: Laser wavelength for intensity correction
        low_freq_cutoff: Suppress modes below this frequency
        title: Plot title
        save_path: Path to save figure (optional)

    Returns:
        Plotly Figure object
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
    fig = go.Figure()

    # Plot each spectrum with offset
    for i, (spectrum, offset, label, color) in enumerate(zip(spectra, y_offsets, labels, colors)):
        offset_spectrum = spectrum + offset

        # Add main spectrum line
        fig.add_trace(go.Scatter(
            x=shifts,
            y=offset_spectrum,
            mode='lines',
            name=label,
            line=dict(color=color, width=1.5),
            hovertemplate='Raman Shift: %{x:.1f} cm⁻¹<br>Intensity: %{y:.3f}<extra></extra>'
        ))

        # Add stick spectrum below baseline
        if 'frequencies' in datasets[i]:
            freqs = datasets[i]['frequencies']
            stick_x = []
            stick_y = []
            for freq in freqs:
                if shift_range[0] <= freq <= shift_range[1]:
                    stick_height = 0.05  # Relative to normalized spectrum
                    stick_x.extend([freq, freq, None])
                    stick_y.extend([offset - stick_height, offset, None])

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
        boundaries = get_raman_region_boundaries()
        for boundary in boundaries:
            if shift_range[0] <= boundary <= shift_range[1]:
                fig.add_vline(
                    x=boundary,
                    line_dash="dash",
                    line_color="gray",
                    opacity=0.5,
                    line_width=0.8
                )

    # Formatting
    fig.update_layout(
        title=dict(text=title, font=dict(size=14, family='Arial Black')),
        xaxis_title=dict(text='Raman Shift (cm⁻¹)', font=dict(size=12, family='Arial Black')),
        yaxis_title=dict(text='Normalized Intensity (offset)', font=dict(size=12, family='Arial Black')),
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
    fig.update_xaxes(autorange='reversed', range=[shift_range[1], shift_range[0]])

    # Hide y-axis ticks (arbitrary offsets)
    fig.update_yaxes(showticklabels=False)

    # Save if requested
    if save_path:
        logger.info(f"Saving figure to {save_path}")
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=1200, height=800, scale=3)

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
) -> go.Figure:
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
        Plotly Figure object
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

    max_spectrum = np.max(spectrum)

    # Create figure
    fig = go.Figure()

    # Plot continuous spectrum
    fig.add_trace(go.Scatter(
        x=shifts,
        y=spectrum,
        mode='lines',
        line=dict(color='black', width=1.5),
        fill='tozeroy',
        fillcolor='rgba(0, 0, 255, 0.2)',
        name='Raman Spectrum',
        hovertemplate='Raman Shift: %{x:.1f} cm⁻¹<br>Intensity: %{y:.3f}<extra></extra>'
    ))

    # Plot stick spectrum below
    stick_x = []
    stick_y = []
    for freq, inten in zip(frequencies, intensities):
        if shift_range[0] <= freq <= shift_range[1]:
            stick_height = max_spectrum * 0.1
            stick_x.extend([freq, freq, None])
            stick_y.extend([-stick_height, 0, None])

    if stick_x:
        fig.add_trace(go.Scatter(
            x=stick_x,
            y=stick_y,
            mode='lines',
            line=dict(color='red', width=1.2),
            opacity=0.7,
            showlegend=False,
            hoverinfo='skip'
        ))

    # Add regional boundaries
    if show_regions:
        boundaries = get_raman_region_boundaries()
        for boundary in boundaries:
            if shift_range[0] <= boundary <= shift_range[1]:
                fig.add_vline(
                    x=boundary,
                    line_dash="dash",
                    line_color="gray",
                    opacity=0.5,
                    line_width=0.8
                )

    # Add region labels
    if show_labels:
        region_labels = get_raman_region_labels()
        for shift, label_text in region_labels:
            if shift_range[0] <= shift <= shift_range[1]:
                fig.add_annotation(
                    x=shift,
                    y=max_spectrum * 0.95,
                    text=label_text,
                    showarrow=False,
                    bgcolor='white',
                    bordercolor='gray',
                    borderwidth=1,
                    borderpad=3,
                    font=dict(size=9)
                )

    # Formatting
    fig.update_layout(
        title=dict(text=title, font=dict(size=14, family='Arial Black')),
        xaxis_title=dict(text='Raman Shift (cm⁻¹)', font=dict(size=12, family='Arial Black')),
        yaxis_title=dict(text='Intensity (arb. units)', font=dict(size=12, family='Arial Black')),
        template='plotly_white',
        hovermode='closest',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=False
    )

    # Invert x-axis
    fig.update_xaxes(autorange='reversed', range=[shift_range[1], shift_range[0]])
    fig.update_yaxes(range=[-max_spectrum * 0.15, max_spectrum * 1.1])

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=1200, height=600, scale=3)

    return fig


# Export main functions
__all__ = ['create_stacked_raman_plot', 'create_single_raman_with_regions']
