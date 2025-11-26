"""
Advanced Subplot Layout Utilities for Multi-Panel Visualizations

Provides flexible grid layouts for combining different plot types:
- IR + Raman side-by-side or stacked
- Spectroscopy + Orbital diagrams
- Comparison panels (Experimental vs Calculated)
- Custom grid arrangements with shared axes

Based on advanced visualization needs from ORCA output analysis.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import List, Dict, Optional, Tuple, Callable, Any
import logging

logger = logging.getLogger(__name__)


class SubplotLayout:
    """
    Manager for complex subplot layouts with multiple plot types.

    Features:
    - Flexible grid layouts (rows x columns)
    - Shared axes between subplots
    - Different subplot sizes
    - Easy integration with existing plot functions
    """

    def __init__(self, rows: int = 2, cols: int = 2,
                 figsize: Tuple[float, float] = (14, 10),
                 height_ratios: Optional[List[float]] = None,
                 width_ratios: Optional[List[float]] = None,
                 hspace: float = 0.3,
                 wspace: float = 0.3):
        """
        Initialize subplot layout manager.

        Args:
            rows: Number of rows
            cols: Number of columns
            figsize: Figure size (width, height)
            height_ratios: Relative heights of rows (e.g., [2, 1] for 2:1 ratio)
            width_ratios: Relative widths of columns
            hspace: Vertical spacing between subplots
            wspace: Horizontal spacing between subplots
        """
        self.rows = rows
        self.cols = cols
        self.figsize = figsize

        # Create subplots with specs
        specs = [[{"secondary_y": False} for _ in range(cols)] for _ in range(rows)]

        self.fig = make_subplots(
            rows=rows, cols=cols,
            row_heights=height_ratios,
            column_widths=width_ratios,
            vertical_spacing=hspace,
            horizontal_spacing=wspace,
            specs=specs
        )

        self.axes = {}
        self.current_row = 1
        self.current_col = 1

        logger.info(f"Created {rows}x{cols} subplot layout")

    def add_subplot(self, row: int, col: int,
                   rowspan: int = 1, colspan: int = 1,
                   name: Optional[str] = None) -> Tuple[int, int]:
        """
        Add a subplot to the layout.

        Args:
            row: Starting row (0-indexed, will be converted to 1-indexed for Plotly)
            col: Starting column (0-indexed, will be converted to 1-indexed for Plotly)
            rowspan: Number of rows to span
            colspan: Number of columns to span
            name: Optional name for referencing this subplot

        Returns:
            (row, col) tuple for plotly (1-indexed)
        """
        # Convert to 1-indexed for Plotly
        plotly_row = row + 1
        plotly_col = col + 1

        if name is None:
            name = f"ax_{row}_{col}"

        self.axes[name] = (plotly_row, plotly_col)
        logger.debug(f"Added subplot '{name}' at ({plotly_row},{plotly_col})")

        return plotly_row, plotly_col

    def get_axis(self, name: str) -> Optional[Tuple[int, int]]:
        """Get axis location by name."""
        return self.axes.get(name)

    def save(self, path: str, dpi: int = 300):
        """Save figure to file."""
        if path.endswith('.html'):
            self.fig.write_html(path)
        else:
            self.fig.write_image(path, width=self.figsize[0] * 80, height=self.figsize[1] * 80, scale=3)
        logger.info(f"Saved figure to {path}")

    def show(self):
        """Display figure."""
        self.fig.show()


def create_ir_raman_comparison(
    ir_data: Dict,
    raman_data: Dict,
    layout: str = 'vertical',
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> go.Figure:
    """
    Create side-by-side or stacked IR and Raman comparison.

    Args:
        ir_data: Dict with 'frequencies', 'intensities' for IR
        raman_data: Dict with 'frequencies', 'activities' for Raman
        layout: 'vertical' (stacked) or 'horizontal' (side-by-side)
        figsize: Figure size (auto-calculated if None)
        save_path: Save path

    Returns:
        plotly Figure object
    """
    from visualization_utils import gaussian_broadening, normalize_spectrum

    logger.info(f"Creating IR/Raman comparison with {layout} layout")

    # Determine layout
    if layout == 'vertical':
        rows, cols = 2, 1
        if figsize is None:
            figsize = (10, 12)
    else:  # horizontal
        rows, cols = 1, 2
        if figsize is None:
            figsize = (16, 6)

    fig = make_subplots(
        rows=rows, cols=cols,
        vertical_spacing=0.15,
        horizontal_spacing=0.1,
        subplot_titles=("IR Spectrum", "Raman Spectrum") if layout == 'horizontal' else None
    )

    # Process IR data
    ir_freq = np.array(ir_data['frequencies'])
    ir_inten = np.array(ir_data['intensities'])

    x_ir = np.linspace(400, 4000, 2000)
    ir_spectrum = gaussian_broadening(x_ir, list(zip(ir_freq, ir_inten)), fwhm=20.0)
    ir_spectrum = normalize_spectrum(ir_spectrum)

    # IR subplot
    ir_row, ir_col = (1, 1) if layout == 'vertical' else (1, 1)
    fig.add_trace(go.Scatter(
        x=x_ir, y=ir_spectrum,
        mode='lines', name='IR',
        line=dict(color='blue', width=1.5)
    ), row=ir_row, col=ir_col)

    # Process Raman data
    raman_freq = np.array(raman_data['frequencies'])
    raman_act = np.array(raman_data['activities'])

    x_raman = np.linspace(50, 3400, 2000)
    raman_spectrum = gaussian_broadening(x_raman, list(zip(raman_freq, raman_act)), fwhm=17.0)
    raman_spectrum = normalize_spectrum(raman_spectrum)

    # Raman subplot
    raman_row, raman_col = (2, 1) if layout == 'vertical' else (1, 2)
    fig.add_trace(go.Scatter(
        x=x_raman, y=raman_spectrum,
        mode='lines', name='Raman',
        line=dict(color='red', width=1.5)
    ), row=raman_row, col=raman_col)

    # Update axes
    fig.update_xaxes(title_text='Wavenumber (cm⁻¹)', row=ir_row, col=ir_col, autorange='reversed', showgrid=True)
    fig.update_yaxes(title_text='Normalized Intensity', row=ir_row, col=ir_col, showgrid=True)

    fig.update_xaxes(title_text='Raman Shift (cm⁻¹)', row=raman_row, col=raman_col, autorange='reversed', showgrid=True)
    fig.update_yaxes(title_text='Normalized Intensity', row=raman_row, col=raman_col, showgrid=True)

    fig.update_layout(
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=False,
        title_text="IR and Raman Spectroscopy Comparison" if layout == 'vertical' else None
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)

    return fig


def create_spectroscopy_orbitals_layout(
    spectrum_data: Dict,
    orbital_data: Dict,
    spectrum_type: str = 'raman',
    figsize: Tuple[float, float] = (16, 10),
    save_path: Optional[str] = None
) -> go.Figure:
    """
    Create layout with spectrum on left and orbital diagram on right.

    Args:
        spectrum_data: Dict with spectroscopy data
        orbital_data: Dict with orbital energies
        spectrum_type: 'raman' or 'ir'
        figsize: Figure size
        save_path: Save path

    Returns:
        plotly Figure object
    """
    logger.info(f"Creating {spectrum_type} + orbitals layout")

    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.67, 0.33],
        horizontal_spacing=0.15,
        subplot_titles=(f'{spectrum_type.upper()} Spectrum', 'Molecular Orbitals')
    )

    # Spectrum subplot (left, larger)
    from visualization_utils import gaussian_broadening, normalize_spectrum

    frequencies = np.array(spectrum_data['frequencies'])

    if spectrum_type == 'raman':
        intensities = np.array(spectrum_data['activities'])
        x = np.linspace(50, 3400, 2000)
        fwhm = 17.0
        xlabel = 'Raman Shift (cm⁻¹)'
        color = 'red'
    else:  # ir
        intensities = np.array(spectrum_data['intensities'])
        x = np.linspace(400, 4000, 2000)
        fwhm = 20.0
        xlabel = 'Wavenumber (cm⁻¹)'
        color = 'blue'

    spectrum = gaussian_broadening(x, list(zip(frequencies, intensities)), fwhm)
    spectrum = normalize_spectrum(spectrum)

    fig.add_trace(go.Scatter(
        x=x, y=spectrum,
        mode='lines', name=spectrum_type.upper(),
        line=dict(color=color, width=1.5)
    ), row=1, col=1)

    # Orbital diagram (right, smaller)
    from orbital_visualization import find_homo_lumo

    orbitals = orbital_data['orbital_energies']
    homo_idx, lumo_idx = find_homo_lumo(orbitals)

    # Select orbitals around HOMO-LUMO
    n_orbitals = 10
    n_below = n_orbitals // 2
    n_above = n_orbitals - n_below

    start_idx = max(0, homo_idx - n_below + 1)
    end_idx = min(len(orbitals), lumo_idx + n_above)

    selected_orbitals = orbitals[start_idx:end_idx]

    energies = [orb.get('energy_ev', 0.0) for orb in selected_orbitals]
    occupations = [orb.get('occupation', 0.0) for orb in selected_orbitals]

    # Draw orbital levels
    for energy, occupation in zip(energies, occupations):
        color = 'blue' if occupation > 0.5 else 'red'
        fig.add_trace(go.Scatter(
            x=[0.3, 0.7], y=[energy, energy],
            mode='lines',
            line=dict(color=color, width=2.5),
            showlegend=False,
            hoverinfo='skip'
        ), row=1, col=2)

    # Mark HOMO and LUMO
    homo_energy = orbitals[homo_idx].get('energy_ev', 0.0)
    lumo_energy = orbitals[lumo_idx].get('energy_ev', 0.0)
    gap = lumo_energy - homo_energy

    # Add HOMO/LUMO annotations
    fig.add_annotation(x=0.1, y=homo_energy, text='HOMO', showarrow=False,
                      font=dict(size=9), xref='x2', yref='y2')
    fig.add_annotation(x=0.1, y=lumo_energy, text='LUMO', showarrow=False,
                      font=dict(size=9), xref='x2', yref='y2')
    fig.add_annotation(x=0.8, y=(homo_energy + lumo_energy) / 2,
                      text=f'{gap:.2f} eV', showarrow=False,
                      font=dict(size=9), xref='x2', yref='y2')

    # Update axes
    fig.update_xaxes(title_text=xlabel, row=1, col=1, autorange='reversed', showgrid=True)
    fig.update_yaxes(title_text='Normalized Intensity', row=1, col=1, showgrid=True)

    fig.update_xaxes(title_text='', row=1, col=2, range=[0, 1], showticklabels=False, showgrid=False)
    fig.update_yaxes(title_text='Orbital Energy (eV)', row=1, col=2, showgrid=True)

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

    return fig


def create_comparison_panel(
    experimental_data: Dict,
    calculated_data: Dict,
    spectrum_type: str = 'ir',
    figsize: Tuple[float, float] = (14, 10),
    save_path: Optional[str] = None
) -> go.Figure:
    """
    Create 3-panel comparison: Experimental, Calculated, and Overlay.

    Args:
        experimental_data: Dict with 'x' and 'y' for experimental spectrum
        calculated_data: Dict with 'frequencies' and 'intensities' for DFT
        spectrum_type: 'ir' or 'raman'
        figsize: Figure size
        save_path: Save path

    Returns:
        plotly Figure object
    """
    logger.info(f"Creating experimental vs calculated {spectrum_type} comparison")

    fig = make_subplots(
        rows=3, cols=1,
        row_heights=[1, 1, 1],
        vertical_spacing=0.1,
        subplot_titles=('Experimental Spectrum', 'Calculated Spectrum (DFT)', 'Overlay')
    )

    from visualization_utils import (
        gaussian_broadening, normalize_spectrum,
        align_spectra_to_common_grid, calculate_spectrum_similarity
    )

    # Experimental (top)
    x_exp = np.array(experimental_data['x'])
    y_exp = np.array(experimental_data['y'])
    y_exp_norm = normalize_spectrum(y_exp)

    fig.add_trace(go.Scatter(
        x=x_exp, y=y_exp_norm,
        mode='lines', name='Experimental',
        line=dict(color='black', width=1.5)
    ), row=1, col=1)

    # Calculated (middle)
    frequencies = np.array(calculated_data['frequencies'])

    if spectrum_type == 'raman':
        intensities = np.array(calculated_data['activities'])
        x_range = (50, 3400)
        fwhm = 17.0
        xlabel = 'Raman Shift (cm⁻¹)'
    else:
        intensities = np.array(calculated_data['intensities'])
        x_range = (400, 4000)
        fwhm = 20.0
        xlabel = 'Wavenumber (cm⁻¹)'

    x_calc = np.linspace(x_range[0], x_range[1], 2000)
    y_calc = gaussian_broadening(x_calc, list(zip(frequencies, intensities)), fwhm)
    y_calc_norm = normalize_spectrum(y_calc)

    fig.add_trace(go.Scatter(
        x=x_calc, y=y_calc_norm,
        mode='lines', name='Calculated',
        line=dict(color='blue', width=1.5)
    ), row=2, col=1)

    # Overlay (bottom)
    # Align spectra for comparison
    x_common, [y_exp_interp, y_calc_interp] = align_spectra_to_common_grid(
        [(x_exp, y_exp_norm), (x_calc, y_calc_norm)],
        x_min=max(x_exp.min(), x_calc.min()),
        x_max=min(x_exp.max(), x_calc.max()),
        num_points=1000
    )

    # Calculate similarity
    correlation = calculate_spectrum_similarity(y_exp_interp, y_calc_interp, 'correlation')
    r2 = calculate_spectrum_similarity(y_exp_interp, y_calc_interp, 'r2')

    fig.add_trace(go.Scatter(
        x=x_common, y=y_exp_interp,
        mode='lines', name='Experimental',
        line=dict(color='black', width=1.5),
        opacity=0.7
    ), row=3, col=1)

    fig.add_trace(go.Scatter(
        x=x_common, y=y_calc_interp,
        mode='lines', name='Calculated',
        line=dict(color='blue', width=1.5),
        opacity=0.7
    ), row=3, col=1)

    # Update all axes
    for i in range(1, 4):
        fig.update_xaxes(autorange='reversed', showgrid=True, row=i, col=1)
        fig.update_yaxes(title_text='Normalized Intensity', showgrid=True, row=i, col=1)

    fig.update_xaxes(title_text=xlabel, row=3, col=1)

    # Update subplot title for overlay with metrics
    fig.layout.annotations[2].update(text=f'Overlay (Correlation: {correlation:.3f}, R²: {r2:.3f})')

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

    return fig


# Export main functions
__all__ = [
    'SubplotLayout',
    'create_ir_raman_comparison',
    'create_spectroscopy_orbitals_layout',
    'create_comparison_panel'
]
