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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

        self.fig = plt.figure(figsize=figsize)
        self.gs = gridspec.GridSpec(rows, cols, figure=self.fig,
                                   height_ratios=height_ratios,
                                   width_ratios=width_ratios,
                                   hspace=hspace, wspace=wspace)

        self.axes = {}
        logger.info(f"Created {rows}x{cols} subplot layout")

    def add_subplot(self, row: int, col: int,
                   rowspan: int = 1, colspan: int = 1,
                   name: Optional[str] = None) -> plt.Axes:
        """
        Add a subplot to the layout.

        Args:
            row: Starting row (0-indexed)
            col: Starting column (0-indexed)
            rowspan: Number of rows to span
            colspan: Number of columns to span
            name: Optional name for referencing this subplot

        Returns:
            matplotlib Axes object
        """
        ax = self.fig.add_subplot(self.gs[row:row+rowspan, col:col+colspan])

        if name is None:
            name = f"ax_{row}_{col}"

        self.axes[name] = ax
        logger.debug(f"Added subplot '{name}' at ({row},{col}) with span ({rowspan},{colspan})")

        return ax

    def get_axis(self, name: str) -> Optional[plt.Axes]:
        """Get axis by name."""
        return self.axes.get(name)

    def save(self, path: str, dpi: int = 300):
        """Save figure to file."""
        self.fig.savefig(path, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved figure to {path}")

    def show(self):
        """Display figure."""
        plt.tight_layout()
        plt.show()


def create_ir_raman_comparison(
    ir_data: Dict,
    raman_data: Dict,
    layout: str = 'vertical',
    figsize: Optional[Tuple[float, float]] = None,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create side-by-side or stacked IR and Raman comparison.

    Args:
        ir_data: Dict with 'frequencies', 'intensities' for IR
        raman_data: Dict with 'frequencies', 'activities' for Raman
        layout: 'vertical' (stacked) or 'horizontal' (side-by-side)
        figsize: Figure size (auto-calculated if None)
        save_path: Save path

    Returns:
        matplotlib Figure object
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

    subplot_layout = SubplotLayout(rows, cols, figsize=figsize,
                                   hspace=0.3, wspace=0.3)

    # IR subplot
    if layout == 'vertical':
        ax_ir = subplot_layout.add_subplot(0, 0, name='ir')
    else:
        ax_ir = subplot_layout.add_subplot(0, 0, name='ir')

    ir_freq = np.array(ir_data['frequencies'])
    ir_inten = np.array(ir_data['intensities'])

    x_ir = np.linspace(400, 4000, 2000)
    ir_spectrum = gaussian_broadening(x_ir, list(zip(ir_freq, ir_inten)), fwhm=20.0)
    ir_spectrum = normalize_spectrum(ir_spectrum)

    ax_ir.plot(x_ir, ir_spectrum, 'b-', linewidth=1.5)
    ax_ir.invert_xaxis()
    ax_ir.set_xlabel('Wavenumber (cm⁻¹)', fontweight='bold')
    ax_ir.set_ylabel('Normalized Intensity', fontweight='bold')
    ax_ir.set_title('IR Spectrum', fontweight='bold')
    ax_ir.grid(True, alpha=0.3)
    ax_ir.spines['top'].set_visible(False)
    ax_ir.spines['right'].set_visible(False)

    # Raman subplot
    if layout == 'vertical':
        ax_raman = subplot_layout.add_subplot(1, 0, name='raman')
    else:
        ax_raman = subplot_layout.add_subplot(0, 1, name='raman')

    raman_freq = np.array(raman_data['frequencies'])
    raman_act = np.array(raman_data['activities'])

    x_raman = np.linspace(50, 3400, 2000)
    raman_spectrum = gaussian_broadening(x_raman, list(zip(raman_freq, raman_act)), fwhm=17.0)
    raman_spectrum = normalize_spectrum(raman_spectrum)

    ax_raman.plot(x_raman, raman_spectrum, 'r-', linewidth=1.5)
    ax_raman.invert_xaxis()
    ax_raman.set_xlabel('Raman Shift (cm⁻¹)', fontweight='bold')
    ax_raman.set_ylabel('Normalized Intensity', fontweight='bold')
    ax_raman.set_title('Raman Spectrum', fontweight='bold')
    ax_raman.grid(True, alpha=0.3)
    ax_raman.spines['top'].set_visible(False)
    ax_raman.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        subplot_layout.save(save_path)

    return subplot_layout.fig


def create_spectroscopy_orbitals_layout(
    spectrum_data: Dict,
    orbital_data: Dict,
    spectrum_type: str = 'raman',
    figsize: Tuple[float, float] = (16, 10),
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create layout with spectrum on left and orbital diagram on right.

    Args:
        spectrum_data: Dict with spectroscopy data
        orbital_data: Dict with orbital energies
        spectrum_type: 'raman' or 'ir'
        figsize: Figure size
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info(f"Creating {spectrum_type} + orbitals layout")

    subplot_layout = SubplotLayout(1, 2, figsize=figsize,
                                   width_ratios=[2, 1], wspace=0.3)

    # Spectrum subplot (left, larger)
    ax_spec = subplot_layout.add_subplot(0, 0, name='spectrum')

    from visualization_utils import gaussian_broadening, normalize_spectrum

    frequencies = np.array(spectrum_data['frequencies'])

    if spectrum_type == 'raman':
        intensities = np.array(spectrum_data['activities'])
        x = np.linspace(50, 3400, 2000)
        fwhm = 17.0
        xlabel = 'Raman Shift (cm⁻¹)'
        title = 'Raman Spectrum'
        color = 'red'
    else:  # ir
        intensities = np.array(spectrum_data['intensities'])
        x = np.linspace(400, 4000, 2000)
        fwhm = 20.0
        xlabel = 'Wavenumber (cm⁻¹)'
        title = 'IR Spectrum'
        color = 'blue'

    spectrum = gaussian_broadening(x, list(zip(frequencies, intensities)), fwhm)
    spectrum = normalize_spectrum(spectrum)

    ax_spec.plot(x, spectrum, color=color, linewidth=1.5)
    ax_spec.invert_xaxis()
    ax_spec.set_xlabel(xlabel, fontweight='bold')
    ax_spec.set_ylabel('Normalized Intensity', fontweight='bold')
    ax_spec.set_title(title, fontweight='bold')
    ax_spec.grid(True, alpha=0.3)
    ax_spec.spines['top'].set_visible(False)
    ax_spec.spines['right'].set_visible(False)

    # Orbital diagram (right, smaller)
    ax_orb = subplot_layout.add_subplot(0, 1, name='orbitals')

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

    for energy, occupation in zip(energies, occupations):
        color = 'blue' if occupation > 0.5 else 'red'
        ax_orb.hlines(energy, 0.3, 0.7, colors=color, linewidth=2.5)

    # Mark HOMO and LUMO
    homo_energy = orbitals[homo_idx].get('energy_ev', 0.0)
    lumo_energy = orbitals[lumo_idx].get('energy_ev', 0.0)
    gap = lumo_energy - homo_energy

    ax_orb.text(0.1, homo_energy, 'HOMO', fontsize=9, fontweight='bold', va='center')
    ax_orb.text(0.1, lumo_energy, 'LUMO', fontsize=9, fontweight='bold', va='center')
    ax_orb.text(0.8, (homo_energy + lumo_energy) / 2, f'{gap:.2f} eV',
               fontsize=9, fontweight='bold', va='center')

    ax_orb.set_ylabel('Orbital Energy (eV)', fontweight='bold')
    ax_orb.set_title('Molecular Orbitals', fontweight='bold')
    ax_orb.set_xlim(0, 1)
    ax_orb.set_xticks([])
    ax_orb.grid(True, axis='y', alpha=0.3)
    ax_orb.spines['top'].set_visible(False)
    ax_orb.spines['right'].set_visible(False)
    ax_orb.spines['bottom'].set_visible(False)

    plt.tight_layout()

    if save_path:
        subplot_layout.save(save_path)

    return subplot_layout.fig


def create_comparison_panel(
    experimental_data: Dict,
    calculated_data: Dict,
    spectrum_type: str = 'ir',
    figsize: Tuple[float, float] = (14, 10),
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create 3-panel comparison: Experimental, Calculated, and Overlay.

    Args:
        experimental_data: Dict with 'x' and 'y' for experimental spectrum
        calculated_data: Dict with 'frequencies' and 'intensities' for DFT
        spectrum_type: 'ir' or 'raman'
        figsize: Figure size
        save_path: Save path

    Returns:
        matplotlib Figure object
    """
    logger.info(f"Creating experimental vs calculated {spectrum_type} comparison")

    subplot_layout = SubplotLayout(3, 1, figsize=figsize,
                                   height_ratios=[1, 1, 1], hspace=0.3)

    from visualization_utils import (
        gaussian_broadening, normalize_spectrum,
        align_spectra_to_common_grid, calculate_spectrum_similarity
    )

    # Experimental (top)
    ax_exp = subplot_layout.add_subplot(0, 0, name='experimental')

    x_exp = np.array(experimental_data['x'])
    y_exp = np.array(experimental_data['y'])
    y_exp_norm = normalize_spectrum(y_exp)

    ax_exp.plot(x_exp, y_exp_norm, 'k-', linewidth=1.5, label='Experimental')
    ax_exp.set_ylabel('Normalized Intensity', fontweight='bold')
    ax_exp.set_title('Experimental Spectrum', fontweight='bold')
    ax_exp.legend()
    ax_exp.grid(True, alpha=0.3)
    ax_exp.invert_xaxis()

    # Calculated (middle)
    ax_calc = subplot_layout.add_subplot(1, 0, name='calculated')

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

    ax_calc.plot(x_calc, y_calc_norm, 'b-', linewidth=1.5, label='Calculated')
    ax_calc.set_ylabel('Normalized Intensity', fontweight='bold')
    ax_calc.set_title('Calculated Spectrum (DFT)', fontweight='bold')
    ax_calc.legend()
    ax_calc.grid(True, alpha=0.3)
    ax_calc.invert_xaxis()

    # Overlay (bottom)
    ax_overlay = subplot_layout.add_subplot(2, 0, name='overlay')

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

    ax_overlay.plot(x_common, y_exp_interp, 'k-', linewidth=1.5, alpha=0.7, label='Experimental')
    ax_overlay.plot(x_common, y_calc_interp, 'b-', linewidth=1.5, alpha=0.7, label='Calculated')
    ax_overlay.set_xlabel(xlabel, fontweight='bold')
    ax_overlay.set_ylabel('Normalized Intensity', fontweight='bold')
    ax_overlay.set_title(f'Overlay (Correlation: {correlation:.3f}, R²: {r2:.3f})',
                        fontweight='bold')
    ax_overlay.legend()
    ax_overlay.grid(True, alpha=0.3)
    ax_overlay.invert_xaxis()

    for ax in [ax_exp, ax_calc, ax_overlay]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        subplot_layout.save(save_path)

    return subplot_layout.fig


# Export main functions
__all__ = [
    'SubplotLayout',
    'create_ir_raman_comparison',
    'create_spectroscopy_orbitals_layout',
    'create_comparison_panel'
]
