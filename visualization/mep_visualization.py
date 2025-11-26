"""
Molecular Electrostatic Potential (MEP) Visualization

Create visualizations for MEP analysis including:
- 2D slices through the MEP grid
- Surface MEP on van der Waals surface
- Isosurface contours
- Critical points (nucleophilic/electrophilic sites)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, TwoSlopeNorm
from mpl_toolkits.mplot3d import Axes3D
from typing import Optional, Tuple, List
import logging

logger = logging.getLogger(__name__)


def create_mep_slice_plot(
    mep_data,
    slice_axis: str = 'z',
    slice_position: Optional[float] = None,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 10),
    dpi: int = 300,
    cmap: str = 'RdBu_r',
    show_atoms: bool = True
) -> plt.Figure:
    """
    Create 2D slice visualization through MEP grid.

    Args:
        mep_data: MEPData object
        slice_axis: Axis perpendicular to slice ('x', 'y', or 'z')
        slice_position: Position along slice_axis (Å), or None for middle
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        cmap: Colormap ('RdBu_r' shows negative=blue, positive=red)
        show_atoms: Show atomic positions

    Returns:
        Matplotlib figure object
    """
    if mep_data.potential_grid is None:
        raise ValueError("No MEP grid data available")

    grid = mep_data.potential_grid
    nx, ny, nz = mep_data.grid_dimensions
    ox, oy, oz = mep_data.grid_origin
    dx, dy, dz = mep_data.grid_spacing

    # Determine slice index
    if slice_axis.lower() == 'x':
        axis_idx = 0
        axis_size = nx
        axis_origin = ox
        axis_spacing = dx
        extent_axes = ('y', 'z')
        extent = [oy, oy + ny*dy, oz, oz + nz*dz]
    elif slice_axis.lower() == 'y':
        axis_idx = 1
        axis_size = ny
        axis_origin = oy
        axis_spacing = dy
        extent_axes = ('x', 'z')
        extent = [ox, ox + nx*dx, oz, oz + nz*dz]
    elif slice_axis.lower() == 'z':
        axis_idx = 2
        axis_size = nz
        axis_origin = oz
        axis_spacing = dz
        extent_axes = ('x', 'y')
        extent = [ox, ox + nx*dx, oy, oy + ny*dy]
    else:
        raise ValueError(f"Invalid slice_axis: {slice_axis}")

    if slice_position is None:
        slice_idx = axis_size // 2
    else:
        slice_idx = int((slice_position - axis_origin) / axis_spacing)
        slice_idx = max(0, min(slice_idx, axis_size - 1))

    # Extract slice
    if axis_idx == 0:
        slice_data = grid[slice_idx, :, :]
    elif axis_idx == 1:
        slice_data = grid[:, slice_idx, :]
    else:
        slice_data = grid[:, :, slice_idx]

    # Transpose for correct orientation
    slice_data = slice_data.T

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Create symmetric colormap normalization centered at zero
    vmax = max(abs(mep_data.min_potential), abs(mep_data.max_potential))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Plot 1: MEP slice
    im1 = ax1.imshow(slice_data, cmap=cmap, norm=norm, origin='lower',
                     extent=extent, aspect='auto', interpolation='bilinear')
    ax1.set_xlabel(f'{extent_axes[0].upper()} (Å)', fontsize=11, fontweight='bold')
    ax1.set_ylabel(f'{extent_axes[1].upper()} (Å)', fontsize=11, fontweight='bold')
    ax1.set_title(f'MEP Slice ({slice_axis.upper()}-plane at '
                  f'{axis_origin + slice_idx * axis_spacing:.2f} Å)',
                  fontsize=12, fontweight='bold')

    # Add colorbar
    cbar1 = plt.colorbar(im1, ax=ax1, pad=0.02)
    cbar1.set_label('MEP (a.u.)', rotation=270, labelpad=20, fontsize=10, fontweight='bold')

    # Show atomic positions if requested
    if show_atoms and mep_data.atoms:
        for element, x, y, z in mep_data.atoms:
            # Project onto slice plane
            if slice_axis.lower() == 'x':
                if abs(x - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    ax1.plot(y, z, 'o', color='white', markersize=8,
                            markeredgecolor='black', markeredgewidth=1.5)
                    ax1.text(y, z, element, ha='center', va='center',
                            fontsize=7, fontweight='bold')
            elif slice_axis.lower() == 'y':
                if abs(y - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    ax1.plot(x, z, 'o', color='white', markersize=8,
                            markeredgecolor='black', markeredgewidth=1.5)
                    ax1.text(x, z, element, ha='center', va='center',
                            fontsize=7, fontweight='bold')
            elif slice_axis.lower() == 'z':
                if abs(z - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    ax1.plot(x, y, 'o', color='white', markersize=8,
                            markeredgecolor='black', markeredgewidth=1.5)
                    ax1.text(x, y, element, ha='center', va='center',
                            fontsize=7, fontweight='bold')

    # Plot 2: Contour plot
    im2 = ax2.contourf(slice_data, levels=20, cmap=cmap, norm=norm,
                       extent=extent, origin='lower')
    ax2.contour(slice_data, levels=10, colors='black', alpha=0.3,
                linewidths=0.5, extent=extent, origin='lower')
    ax2.set_xlabel(f'{extent_axes[0].upper()} (Å)', fontsize=11, fontweight='bold')
    ax2.set_ylabel(f'{extent_axes[1].upper()} (Å)', fontsize=11, fontweight='bold')
    ax2.set_title(f'MEP Contours ({slice_axis.upper()}-plane)',
                  fontsize=12, fontweight='bold')

    # Add colorbar
    cbar2 = plt.colorbar(im2, ax=ax2, pad=0.02)
    cbar2.set_label('MEP (a.u.)', rotation=270, labelpad=20, fontsize=10, fontweight='bold')

    # Show atoms on contour plot too
    if show_atoms and mep_data.atoms:
        for element, x, y, z in mep_data.atoms:
            if slice_axis.lower() == 'x':
                if abs(x - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    ax2.plot(y, z, 'o', color='white', markersize=8,
                            markeredgecolor='black', markeredgewidth=1.5)
            elif slice_axis.lower() == 'y':
                if abs(y - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    ax2.plot(x, z, 'o', color='white', markersize=8,
                            markeredgecolor='black', markeredgewidth=1.5)
            elif slice_axis.lower() == 'z':
                if abs(z - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    ax2.plot(x, y, 'o', color='white', markersize=8,
                            markeredgecolor='black', markeredgewidth=1.5)

    plt.suptitle('Molecular Electrostatic Potential (MEP) Analysis',
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved MEP slice plot to {output_file}")

    return fig


def create_mep_multi_slice_plot(
    mep_data,
    axis: str = 'z',
    n_slices: int = 4,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (16, 10),
    dpi: int = 300,
    cmap: str = 'RdBu_r'
) -> plt.Figure:
    """
    Create multiple slices through MEP grid.

    Args:
        mep_data: MEPData object
        axis: Axis perpendicular to slices ('x', 'y', or 'z')
        n_slices: Number of slices to show
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        cmap: Colormap

    Returns:
        Matplotlib figure object
    """
    if mep_data.potential_grid is None:
        raise ValueError("No MEP grid data available")

    grid = mep_data.potential_grid
    nx, ny, nz = mep_data.grid_dimensions
    ox, oy, oz = mep_data.grid_origin
    dx, dy, dz = mep_data.grid_spacing

    # Determine slices
    if axis.lower() == 'x':
        axis_size = nx
        axis_origin = ox
        axis_spacing = dx
        extent_axes = ('y', 'z')
        extent = [oy, oy + ny*dy, oz, oz + nz*dz]
    elif axis.lower() == 'y':
        axis_size = ny
        axis_origin = oy
        axis_spacing = dy
        extent_axes = ('x', 'z')
        extent = [ox, ox + nx*dx, oz, oz + nz*dz]
    elif axis.lower() == 'z':
        axis_size = nz
        axis_origin = oz
        axis_spacing = dz
        extent_axes = ('x', 'y')
        extent = [ox, ox + nx*dx, oy, oy + ny*dy]
    else:
        raise ValueError(f"Invalid axis: {axis}")

    # Calculate slice indices
    slice_indices = np.linspace(0, axis_size - 1, n_slices, dtype=int)

    # Create symmetric colormap normalization
    vmax = max(abs(mep_data.min_potential), abs(mep_data.max_potential))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Create subplots
    ncols = min(n_slices, 4)
    nrows = (n_slices + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n_slices == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, idx in enumerate(slice_indices):
        if axis.lower() == 'x':
            slice_data = grid[idx, :, :].T
        elif axis.lower() == 'y':
            slice_data = grid[:, idx, :].T
        else:
            slice_data = grid[:, :, idx].T

        im = axes[i].imshow(slice_data, cmap=cmap, norm=norm, origin='lower',
                           extent=extent, aspect='auto', interpolation='bilinear')
        axes[i].set_xlabel(f'{extent_axes[0].upper()} (Å)', fontsize=9)
        axes[i].set_ylabel(f'{extent_axes[1].upper()} (Å)', fontsize=9)
        position = axis_origin + idx * axis_spacing
        axes[i].set_title(f'{axis.upper()} = {position:.2f} Å', fontsize=10, fontweight='bold')

    # Hide extra subplots
    for i in range(n_slices, len(axes)):
        axes[i].axis('off')

    # Add colorbar
    fig.colorbar(im, ax=axes, label='MEP (a.u.)', pad=0.02, aspect=30)

    plt.suptitle(f'MEP Multiple Slices ({axis.upper()}-axis)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved MEP multi-slice plot to {output_file}")

    return fig


def create_mep_3d_isosurface_plot(
    mep_data,
    isovalues: Optional[List[float]] = None,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 10),
    dpi: int = 300,
    show_atoms: bool = True,
    max_points: int = 5000
) -> plt.Figure:
    """
    Create 3D visualization of MEP isosurfaces.

    Args:
        mep_data: MEPData object
        isovalues: MEP isovalues to plot (a.u.), or None for automatic
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        show_atoms: Show atomic positions
        max_points: Maximum points to plot (for performance)

    Returns:
        Matplotlib figure object
    """
    if mep_data.potential_grid is None:
        raise ValueError("No MEP grid data available")

    grid = mep_data.potential_grid
    nx, ny, nz = mep_data.grid_dimensions
    ox, oy, oz = mep_data.grid_origin
    dx, dy, dz = mep_data.grid_spacing

    # Create coordinate grids
    x = np.linspace(ox, ox + (nx-1)*dx, nx)
    y = np.linspace(oy, oy + (ny-1)*dy, ny)
    z = np.linspace(oz, oz + (nz-1)*dz, nz)

    # Default isovalues if not provided
    if isovalues is None:
        isovalues = [
            mep_data.min_potential * 0.5,  # Negative (nucleophilic)
            0.0,
            mep_data.max_potential * 0.5   # Positive (electrophilic)
        ]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Sample grid points for visualization (for performance)
    step = max(1, int((nx * ny * nz / max_points) ** (1/3)))
    x_sample = x[::step]
    y_sample = y[::step]
    z_sample = z[::step]
    grid_sample = grid[::step, ::step, ::step]

    # Create meshgrid for sampled points
    X, Y, Z = np.meshgrid(x_sample, y_sample, z_sample, indexing='ij')

    # Flatten for scatter plot
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    Z_flat = Z.flatten()
    V_flat = grid_sample.flatten()

    # Color based on MEP value
    norm = TwoSlopeNorm(vmin=mep_data.min_potential, vcenter=0, vmax=mep_data.max_potential)
    colors = cm.RdBu_r(norm(V_flat))

    # Plot points
    scatter = ax.scatter(X_flat, Y_flat, Z_flat, c=V_flat, cmap='RdBu_r',
                        norm=norm, s=10, alpha=0.3, marker='.')

    # Plot atoms
    if show_atoms and mep_data.atoms:
        atom_x = [atom[1] for atom in mep_data.atoms]
        atom_y = [atom[2] for atom in mep_data.atoms]
        atom_z = [atom[3] for atom in mep_data.atoms]
        ax.scatter(atom_x, atom_y, atom_z, c='black', s=200, marker='o',
                  edgecolors='white', linewidths=2, alpha=0.8)

        # Label atoms
        for element, x, y, z in mep_data.atoms:
            ax.text(x, y, z, element, color='white', fontsize=8,
                   fontweight='bold', ha='center', va='center')

    ax.set_xlabel('X (Å)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Y (Å)', fontsize=10, fontweight='bold')
    ax.set_zlabel('Z (Å)', fontsize=10, fontweight='bold')
    ax.set_title('MEP 3D Isosurface Visualization', fontsize=12, fontweight='bold', pad=20)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('MEP (a.u.)', rotation=270, labelpad=20, fontsize=10, fontweight='bold')

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved MEP 3D isosurface plot to {output_file}")

    return fig


def create_mep_critical_points_plot(
    mep_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 6),
    dpi: int = 300
) -> plt.Figure:
    """
    Create visualization of MEP critical points (reactive sites).

    Args:
        mep_data: MEPData object with critical points calculated
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure

    Returns:
        Matplotlib figure object
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Plot 1: Nucleophilic sites (minima)
    if mep_data.local_minima:
        minima_sorted = sorted(mep_data.local_minima, key=lambda p: p[3])[:10]
        sites = [f"Site {i+1}" for i in range(len(minima_sorted))]
        potentials = [p[3] * 27.2114 for p in minima_sorted]  # Convert to eV

        bars1 = ax1.barh(sites, potentials, color='#3498db', alpha=0.7, edgecolor='black')
        ax1.set_xlabel('MEP (eV)', fontsize=11, fontweight='bold')
        ax1.set_title('Top Nucleophilic Sites (MEP Minima)', fontsize=12, fontweight='bold')
        ax1.axvline(x=0, color='gray', linestyle='--', linewidth=1)
        ax1.grid(True, alpha=0.3, axis='x')

        # Add value labels
        for bar, val in zip(bars1, potentials):
            width = bar.get_width()
            ax1.text(width, bar.get_y() + bar.get_height()/2,
                    f'{val:.2f}',
                    ha='left' if width < 0 else 'left', va='center',
                    fontsize=8, fontweight='bold', color='black')
    else:
        ax1.text(0.5, 0.5, 'No critical points calculated\nUse mep_data.find_critical_points()',
                ha='center', va='center', transform=ax1.transAxes, fontsize=11)

    # Plot 2: Electrophilic sites (maxima)
    if mep_data.local_maxima:
        maxima_sorted = sorted(mep_data.local_maxima, key=lambda p: p[3], reverse=True)[:10]
        sites = [f"Site {i+1}" for i in range(len(maxima_sorted))]
        potentials = [p[3] * 27.2114 for p in maxima_sorted]  # Convert to eV

        bars2 = ax2.barh(sites, potentials, color='#e74c3c', alpha=0.7, edgecolor='black')
        ax2.set_xlabel('MEP (eV)', fontsize=11, fontweight='bold')
        ax2.set_title('Top Electrophilic Sites (MEP Maxima)', fontsize=12, fontweight='bold')
        ax2.axvline(x=0, color='gray', linestyle='--', linewidth=1)
        ax2.grid(True, alpha=0.3, axis='x')

        # Add value labels
        for bar, val in zip(bars2, potentials):
            width = bar.get_width()
            ax2.text(width, bar.get_y() + bar.get_height()/2,
                    f'{val:.2f}',
                    ha='left' if width > 0 else 'right', va='center',
                    fontsize=8, fontweight='bold', color='black')
    else:
        ax2.text(0.5, 0.5, 'No critical points calculated\nUse mep_data.find_critical_points()',
                ha='center', va='center', transform=ax2.transAxes, fontsize=11)

    plt.suptitle('MEP Critical Points Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved MEP critical points plot to {output_file}")

    return fig
