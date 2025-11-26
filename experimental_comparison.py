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
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
) -> go.Figure:
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
        plotly Figure object
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
    fig = make_subplots(
        rows=4, cols=1,
        row_heights=[1, 1, 1, 1],
        vertical_spacing=0.08,
        subplot_titles=(
            f'{title} - Experimental',
            f'DFT Calculated - Scale Factor: {scale_label}',
            'Overlay',
            'Difference (Exp - DFT)'
        )
    )

    # Panel 1: Experimental spectrum (top)
    exp_y_norm = normalize_spectrum(exp_y)
    fig.add_trace(go.Scatter(
        x=exp_x, y=exp_y_norm,
        mode='lines', name='Experimental',
        line=dict(color='black', width=1.5)
    ), row=1, col=1)

    # Panel 2: Calculated spectrum (second)
    # Generate broadened calculated spectrum
    x_calc = np.linspace(wavenumber_range[0], wavenumber_range[1], 2000)
    scaled_freqs = calculated_frequencies * scale_factor
    peaks = list(zip(scaled_freqs, calculated_intensities))
    calc_spectrum = gaussian_broadening(x_calc, peaks, fwhm)
    calc_spectrum_norm = normalize_spectrum(calc_spectrum)

    fig.add_trace(go.Scatter(
        x=x_calc, y=calc_spectrum_norm,
        mode='lines', name=f'DFT (scaled: {scale_factor:.4f})',
        line=dict(color='blue', width=1.5)
    ), row=2, col=1)

    # Add stick spectrum
    for freq, inten in zip(scaled_freqs, calculated_intensities):
        if wavenumber_range[0] <= freq <= wavenumber_range[1]:
            fig.add_trace(go.Scatter(
                x=[freq, freq], y=[0, inten / calculated_intensities.max()],
                mode='lines',
                line=dict(color='red', width=1.0),
                opacity=0.4,
                showlegend=False,
                hoverinfo='skip'
            ), row=2, col=1)

    # Panel 3: Overlay (third)
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

    fig.add_trace(go.Scatter(
        x=x_common, y=exp_aligned,
        mode='lines', name='Experimental',
        line=dict(color='black', width=1.5),
        opacity=0.7
    ), row=3, col=1)

    fig.add_trace(go.Scatter(
        x=x_common, y=calc_aligned,
        mode='lines', name='DFT',
        line=dict(color='blue', width=1.5),
        opacity=0.7
    ), row=3, col=1)

    # Update title with metrics
    fig.layout.annotations[2].update(
        text=f'Overlay - Correlation: {correlation:.4f}, R²: {r2:.4f}, RMSE: {rmse:.4f}'
    )

    # Panel 4: Residuals (bottom)
    residuals = exp_aligned - calc_aligned
    fig.add_trace(go.Scatter(
        x=x_common, y=residuals,
        mode='lines', name='Residuals',
        line=dict(color='red', width=1.0),
        fill='tozeroy',
        fillcolor='rgba(255, 0, 0, 0.3)'
    ), row=4, col=1)

    fig.add_hline(y=0, line=dict(color='black', dash='dash', width=1),
                  opacity=0.5, row=4, col=1)

    # Update all axes
    for i in range(1, 5):
        fig.update_xaxes(autorange='reversed', showgrid=True, row=i, col=1, range=wavenumber_range)
        fig.update_yaxes(showgrid=True, row=i, col=1)

    fig.update_yaxes(title_text='Normalized Intensity', row=1, col=1)
    fig.update_yaxes(title_text='Normalized Intensity', row=2, col=1)
    fig.update_yaxes(title_text='Normalized Intensity', row=3, col=1)
    fig.update_yaxes(title_text='Residuals', row=4, col=1)
    fig.update_xaxes(title_text='Wavenumber (cm⁻¹)', row=4, col=1)

    fig.update_layout(
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=False
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
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
) -> go.Figure:
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
        plotly Figure object
    """
    logger.info("Creating peak assignment table")

    # Match peaks
    matches = match_peaks(experimental_peaks, calculated_peaks, scale_factor, tolerance)

    # Prepare table data (limit to top 20 matches)
    matches = matches[:20]

    header = ['Exp. Freq.<br>(cm⁻¹)', 'Exp. Inten.', 'Calc. Freq.<br>(cm⁻¹)',
              'Calc. Inten.<br>(km/mol)', 'Scaled Freq.<br>(cm⁻¹)', 'Shift<br>(cm⁻¹)']

    # Create cell values
    cell_values = [
        [f"{m['exp_freq']:.1f}" for m in matches],
        [f"{m['exp_intensity']:.2f}" for m in matches],
        [f"{m['calc_freq']:.1f}" for m in matches],
        [f"{m['calc_intensity']:.1f}" for m in matches],
        [f"{m['calc_freq_scaled']:.1f}" for m in matches],
        [f"{m['shift']:+.1f}" for m in matches]
    ]

    # Color-code shifts
    shift_colors = []
    for match in matches:
        shift_abs = abs(match['shift'])
        if shift_abs < 5:
            color = '#C6EFCE'  # Green - excellent match
        elif shift_abs < 10:
            color = '#FFEB9C'  # Yellow - good match
        elif shift_abs < 20:
            color = '#FFC7CE'  # Red - acceptable match
        else:
            color = '#E0E0E0'  # Gray - poor match
        shift_colors.append(color)

    # Create figure
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=header,
            fill_color='#4472C4',
            font=dict(color='white', size=11),
            align='center',
            height=40
        ),
        cells=dict(
            values=cell_values,
            fill_color=[['white'] * len(matches)] * 5 + [shift_colors],
            align='center',
            height=35,
            font=dict(size=10)
        )
    )])

    fig.update_layout(
        title=f'{title}<br>Scale Factor: {scale_factor:.4f}, Tolerance: ±{tolerance} cm⁻¹',
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
    'optimize_frequency_scaling',
    'match_peaks',
    'create_exp_vs_dft_comparison',
    'create_peak_assignment_table'
]
