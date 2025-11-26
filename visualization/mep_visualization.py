"""
Molecular Electrostatic Potential (MEP) Visualization

Create visualizations for MEP analysis including:
- 2D slices through the MEP grid
- Surface MEP on van der Waals surface
- Isosurface contours
- Critical points (nucleophilic/electrophilic sites)
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
) -> go.Figure:
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
        Plotly figure object
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
    elif slice_axis.lower() == 'y':
        axis_idx = 1
        axis_size = ny
        axis_origin = oy
        axis_spacing = dy
        extent_axes = ('x', 'z')
    elif slice_axis.lower() == 'z':
        axis_idx = 2
        axis_size = nz
        axis_origin = oz
        axis_spacing = dz
        extent_axes = ('x', 'y')
    else:
        raise ValueError(f"Invalid slice_axis: {slice_axis}")

    if slice_position is None:
        slice_idx = axis_size // 2
    else:
        slice_idx = int((slice_position - axis_origin) / axis_spacing)
        slice_idx = max(0, min(slice_idx, axis_size - 1))

    # Extract slice
    if axis_idx == 0:
        slice_data = grid[slice_idx, :, :].T
    elif axis_idx == 1:
        slice_data = grid[:, slice_idx, :].T
    else:
        slice_data = grid[:, :, slice_idx].T

    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            f'MEP Slice ({slice_axis.upper()}-plane at {axis_origin + slice_idx * axis_spacing:.2f} Å)',
            f'MEP Contours ({slice_axis.upper()}-plane)'
        ),
        horizontal_spacing=0.1
    )

    # Create symmetric colormap centered at zero
    vmax = max(abs(mep_data.min_potential), abs(mep_data.max_potential))

    # Plot 1: MEP slice (heatmap)
    fig.add_trace(go.Heatmap(
        z=slice_data,
        colorscale=cmap,
        zmid=0,
        zmin=-vmax,
        zmax=vmax,
        colorbar=dict(title='MEP (a.u.)', x=0.45)
    ), row=1, col=1)

    # Plot 2: Contour plot
    fig.add_trace(go.Contour(
        z=slice_data,
        colorscale=cmap,
        contours=dict(
            coloring='fill',
            showlabels=True
        ),
        line=dict(width=0.5),
        zmid=0,
        zmin=-vmax,
        zmax=vmax,
        colorbar=dict(title='MEP (a.u.)', x=1.0)
    ), row=1, col=2)

    # Show atoms on both plots if requested
    if show_atoms and mep_data.atoms:
        for element, x, y, z in mep_data.atoms:
            # Project onto slice plane
            if slice_axis.lower() == 'x':
                if abs(x - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    # Add marker to both plots
                    for col in [1, 2]:
                        fig.add_trace(go.Scatter(
                            x=[y], y=[z],
                            mode='markers+text',
                            marker=dict(size=8, color='white', line=dict(color='black', width=1.5)),
                            text=element,
                            textposition='middle center',
                            textfont=dict(size=7, color='black'),
                            showlegend=False,
                            hoverinfo='skip'
                        ), row=1, col=col)
            elif slice_axis.lower() == 'y':
                if abs(y - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    for col in [1, 2]:
                        fig.add_trace(go.Scatter(
                            x=[x], y=[z],
                            mode='markers+text',
                            marker=dict(size=8, color='white', line=dict(color='black', width=1.5)),
                            text=element,
                            textposition='middle center',
                            textfont=dict(size=7, color='black'),
                            showlegend=False,
                            hoverinfo='skip'
                        ), row=1, col=col)
            elif slice_axis.lower() == 'z':
                if abs(z - (axis_origin + slice_idx * axis_spacing)) < 2.0:
                    for col in [1, 2]:
                        fig.add_trace(go.Scatter(
                            x=[x], y=[y],
                            mode='markers+text',
                            marker=dict(size=8, color='white', line=dict(color='black', width=1.5)),
                            text=element,
                            textposition='middle center',
                            textfont=dict(size=7, color='black'),
                            showlegend=False,
                            hoverinfo='skip'
                        ), row=1, col=col)

    # Update axes
    fig.update_xaxes(title_text=f'{extent_axes[0].upper()} (Å)', row=1, col=1)
    fig.update_yaxes(title_text=f'{extent_axes[1].upper()} (Å)', row=1, col=1)
    fig.update_xaxes(title_text=f'{extent_axes[0].upper()} (Å)', row=1, col=2)
    fig.update_yaxes(title_text=f'{extent_axes[1].upper()} (Å)', row=1, col=2)

    fig.update_layout(
        title_text='Molecular Electrostatic Potential (MEP) Analysis',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
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
) -> go.Figure:
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
        Plotly figure object
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
    elif axis.lower() == 'y':
        axis_size = ny
        axis_origin = oy
        axis_spacing = dy
    elif axis.lower() == 'z':
        axis_size = nz
        axis_origin = oz
        axis_spacing = dz
    else:
        raise ValueError(f"Invalid axis: {axis}")

    # Calculate slice indices
    slice_indices = np.linspace(0, axis_size - 1, n_slices, dtype=int)

    # Create symmetric colormap normalization
    vmax = max(abs(mep_data.min_potential), abs(mep_data.max_potential))

    # Create subplots
    ncols = min(n_slices, 4)
    nrows = (n_slices + ncols - 1) // ncols

    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=[f'{axis.upper()} = {axis_origin + idx * axis_spacing:.2f} Å'
                       for idx in slice_indices]
    )

    for i, idx in enumerate(slice_indices):
        row = i // ncols + 1
        col = i % ncols + 1

        if axis.lower() == 'x':
            slice_data = grid[idx, :, :].T
        elif axis.lower() == 'y':
            slice_data = grid[:, idx, :].T
        else:
            slice_data = grid[:, :, idx].T

        fig.add_trace(go.Heatmap(
            z=slice_data,
            colorscale=cmap,
            zmid=0,
            zmin=-vmax,
            zmax=vmax,
            showscale=(i == 0),
            colorbar=dict(title='MEP (a.u.)') if i == 0 else None
        ), row=row, col=col)

    fig.update_layout(
        title_text=f'MEP Multiple Slices ({axis.upper()}-axis)',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
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
) -> go.Figure:
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
        Plotly figure object
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

    fig = go.Figure()

    # Plot points colored by MEP value
    fig.add_trace(go.Scatter3d(
        x=X_flat, y=Y_flat, z=Z_flat,
        mode='markers',
        marker=dict(
            size=2,
            color=V_flat,
            colorscale='RdBu_r',
            cmid=0,
            cmin=mep_data.min_potential,
            cmax=mep_data.max_potential,
            colorbar=dict(title='MEP (a.u.)', x=1.0),
            opacity=0.3
        ),
        showlegend=False,
        hoverinfo='skip'
    ))

    # Plot atoms
    if show_atoms and mep_data.atoms:
        atom_x = [atom[1] for atom in mep_data.atoms]
        atom_y = [atom[2] for atom in mep_data.atoms]
        atom_z = [atom[3] for atom in mep_data.atoms]
        atom_elements = [atom[0] for atom in mep_data.atoms]

        fig.add_trace(go.Scatter3d(
            x=atom_x, y=atom_y, z=atom_z,
            mode='markers+text',
            marker=dict(size=8, color='black', line=dict(color='white', width=2)),
            text=atom_elements,
            textposition='middle center',
            textfont=dict(size=8, color='white'),
            showlegend=False,
            hoverinfo='text',
            hovertext=[f'{elem}' for elem in atom_elements]
        ))

    fig.update_layout(
        title='MEP 3D Isosurface Visualization',
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            aspectmode='data'
        ),
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved MEP 3D isosurface plot to {output_file}")

    return fig


def create_mep_critical_points_plot(
    mep_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 6),
    dpi: int = 300
) -> go.Figure:
    """
    Create visualization of MEP critical points (reactive sites).

    Args:
        mep_data: MEPData object with critical points calculated
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure

    Returns:
        Plotly figure object
    """
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('Top Nucleophilic Sites (MEP Minima)',
                       'Top Electrophilic Sites (MEP Maxima)')
    )

    # Plot 1: Nucleophilic sites (minima)
    if mep_data.local_minima:
        minima_sorted = sorted(mep_data.local_minima, key=lambda p: p[3])[:10]
        sites = [f"Site {i+1}" for i in range(len(minima_sorted))]
        potentials = [p[3] * 27.2114 for p in minima_sorted]  # Convert to eV

        fig.add_trace(go.Bar(
            y=sites, x=potentials,
            orientation='h',
            marker=dict(color='#3498db', line=dict(color='black', width=1)),
            opacity=0.7,
            text=[f'{v:.2f}' for v in potentials],
            textposition='outside',
            showlegend=False
        ), row=1, col=1)

        fig.add_vline(x=0, line=dict(color='gray', dash='dash', width=1), row=1, col=1)
    else:
        fig.add_annotation(
            x=0.25, y=0.5,
            text='No critical points calculated<br>Use mep_data.find_critical_points()',
            showarrow=False,
            xref='paper', yref='paper',
            font=dict(size=11)
        )

    # Plot 2: Electrophilic sites (maxima)
    if mep_data.local_maxima:
        maxima_sorted = sorted(mep_data.local_maxima, key=lambda p: p[3], reverse=True)[:10]
        sites = [f"Site {i+1}" for i in range(len(maxima_sorted))]
        potentials = [p[3] * 27.2114 for p in maxima_sorted]  # Convert to eV

        fig.add_trace(go.Bar(
            y=sites, x=potentials,
            orientation='h',
            marker=dict(color='#e74c3c', line=dict(color='black', width=1)),
            opacity=0.7,
            text=[f'{v:.2f}' for v in potentials],
            textposition='outside',
            showlegend=False
        ), row=1, col=2)

        fig.add_vline(x=0, line=dict(color='gray', dash='dash', width=1), row=1, col=2)
    else:
        fig.add_annotation(
            x=0.75, y=0.5,
            text='No critical points calculated<br>Use mep_data.find_critical_points()',
            showarrow=False,
            xref='paper', yref='paper',
            font=dict(size=11)
        )

    # Update axes
    fig.update_xaxes(title_text='MEP (eV)', showgrid=True, gridcolor='lightgray', row=1, col=1)
    fig.update_xaxes(title_text='MEP (eV)', showgrid=True, gridcolor='lightgray', row=1, col=2)

    fig.update_layout(
        title_text='MEP Critical Points Analysis',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved MEP critical points plot to {output_file}")

    return fig
