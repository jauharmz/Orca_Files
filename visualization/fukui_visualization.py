"""
Fukui Function Visualization

Create visualizations for Fukui function analysis including:
- Atomic Fukui indices (f⁺, f⁻, f⁰) as bar plots
- Dual descriptor plots
- Reactivity site identification
- Global reactivity descriptor summaries
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def create_fukui_bar_plot(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 10),
    dpi: int = 300,
    show_labels: bool = True,
    top_n: Optional[int] = None
) -> go.Figure:
    """
    Create bar plot visualization of Fukui indices.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        show_labels: Show atom labels on x-axis
        top_n: Only show top N most reactive atoms (None = show all)

    Returns:
        Plotly figure object
    """
    if not fukui_data.atomic_fukui:
        raise ValueError("No Fukui data available")

    # Sort by f⁰ (overall reactivity) and optionally limit
    sorted_atoms = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_zero, reverse=True)
    if top_n:
        sorted_atoms = sorted_atoms[:top_n]

    n_atoms = len(sorted_atoms)
    atom_indices = [atom.atom_index for atom in sorted_atoms]
    elements = [atom.element for atom in sorted_atoms]
    f_plus = np.array([atom.f_plus for atom in sorted_atoms])
    f_minus = np.array([atom.f_minus for atom in sorted_atoms])
    f_zero = np.array([atom.f_zero for atom in sorted_atoms])
    dual = np.array([atom.dual_descriptor for atom in sorted_atoms])

    x = np.arange(n_atoms)

    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Electrophilic Attack Sites (High f⁺)',
            'Nucleophilic Attack Sites (High f⁻)',
            'Radical Attack Sites (High f⁰)',
            'Dual Descriptor (Δf = f⁻ - f⁺)'
        ),
        vertical_spacing=0.15,
        horizontal_spacing=0.1
    )

    # Plot 1: f⁺ (nucleophilic Fukui - electrophilic attack sites)
    colors_fplus = ['#c0392b' if i in np.argsort(f_plus)[-3:] else '#e74c3c' for i in range(n_atoms)]
    fig.add_trace(go.Bar(
        x=x, y=f_plus,
        marker_color=colors_fplus,
        marker_line=dict(color='black', width=1),
        opacity=0.7,
        showlegend=False,
        hovertemplate='%{y:.4f}<extra></extra>'
    ), row=1, col=1)

    # Plot 2: f⁻ (electrophilic Fukui - nucleophilic attack sites)
    colors_fminus = ['#2980b9' if i in np.argsort(f_minus)[-3:] else '#3498db' for i in range(n_atoms)]
    fig.add_trace(go.Bar(
        x=x, y=f_minus,
        marker_color=colors_fminus,
        marker_line=dict(color='black', width=1),
        opacity=0.7,
        showlegend=False,
        hovertemplate='%{y:.4f}<extra></extra>'
    ), row=1, col=2)

    # Plot 3: f⁰ (radical Fukui - radical attack sites)
    colors_fzero = ['#27ae60' if i in np.argsort(f_zero)[-3:] else '#2ecc71' for i in range(n_atoms)]
    fig.add_trace(go.Bar(
        x=x, y=f_zero,
        marker_color=colors_fzero,
        marker_line=dict(color='black', width=1),
        opacity=0.7,
        showlegend=False,
        hovertemplate='%{y:.4f}<extra></extra>'
    ), row=2, col=1)

    # Plot 4: Dual descriptor (Δf = f⁻ - f⁺)
    colors_dual = ['#3498db' if d > 0 else '#e74c3c' for d in dual]
    fig.add_trace(go.Bar(
        x=x, y=dual,
        marker_color=colors_dual,
        marker_line=dict(color='black', width=1),
        opacity=0.7,
        showlegend=False,
        hovertemplate='%{y:.4f}<extra></extra>'
    ), row=2, col=2)

    # Add legend to dual descriptor
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color='#3498db', opacity=0.7),
        name='Nucleophilic (Δf > 0)',
        showlegend=True
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color='#e74c3c', opacity=0.7),
        name='Electrophilic (Δf < 0)',
        showlegend=True
    ))

    # Update axes
    if show_labels:
        ticktext = [f"{idx}<br>{elem}" for idx, elem in zip(atom_indices, elements)]
    else:
        ticktext = [str(idx) for idx in atom_indices]

    for row in range(1, 3):
        for col in range(1, 3):
            fig.update_xaxes(
                title_text='Atom Index', row=row, col=col,
                tickmode='array', tickvals=x, ticktext=ticktext,
                tickfont=dict(size=8)
            )
            fig.update_yaxes(showgrid=True, gridcolor='lightgray', row=row, col=col)

    fig.update_yaxes(title_text='f⁺ (Nucleophilic Fukui)', row=1, col=1)
    fig.update_yaxes(title_text='f⁻ (Electrophilic Fukui)', row=1, col=2)
    fig.update_yaxes(title_text='f⁰ (Radical Fukui)', row=2, col=1)
    fig.update_yaxes(title_text='Δf (Dual Descriptor)', row=2, col=2)

    # Add horizontal lines at y=0
    for row, col in [(1, 1), (1, 2), (2, 1)]:
        fig.add_hline(y=0, line=dict(color='gray', dash='dash', width=0.8), row=row, col=col)
    fig.add_hline(y=0, line=dict(color='black', width=1.5), row=2, col=2)

    fig.update_layout(
        title_text=f'Fukui Function Analysis ({fukui_data.method.upper()} Charges)',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=True,
        legend=dict(x=0.75, y=0.25)
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved Fukui bar plot to {output_file}")

    return fig


def create_fukui_comparison_plot(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 6),
    dpi: int = 300,
    top_n: int = 15
) -> go.Figure:
    """
    Create side-by-side comparison of all three Fukui indices.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        top_n: Show top N atoms by f⁰

    Returns:
        Plotly figure object
    """
    if not fukui_data.atomic_fukui:
        raise ValueError("No Fukui data available")

    # Sort by f⁰ and limit
    sorted_atoms = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_zero, reverse=True)[:top_n]

    n_atoms = len(sorted_atoms)
    atom_labels = [f"{atom.atom_index}-{atom.element}" for atom in sorted_atoms]
    f_plus = np.array([atom.f_plus for atom in sorted_atoms])
    f_minus = np.array([atom.f_minus for atom in sorted_atoms])
    f_zero = np.array([atom.f_zero for atom in sorted_atoms])

    x = np.arange(n_atoms)
    width = 0.25

    fig = go.Figure()

    # Plot all three Fukui indices
    fig.add_trace(go.Bar(
        x=x - width, y=f_plus, width=width,
        name='f⁺ (Electrophilic attack)',
        marker=dict(color='#e74c3c', line=dict(color='black', width=1)),
        opacity=0.7
    ))

    fig.add_trace(go.Bar(
        x=x, y=f_minus, width=width,
        name='f⁻ (Nucleophilic attack)',
        marker=dict(color='#3498db', line=dict(color='black', width=1)),
        opacity=0.7
    ))

    fig.add_trace(go.Bar(
        x=x + width, y=f_zero, width=width,
        name='f⁰ (Radical attack)',
        marker=dict(color='#2ecc71', line=dict(color='black', width=1)),
        opacity=0.7
    ))

    fig.add_hline(y=0, line=dict(color='gray', dash='dash', width=0.8))

    fig.update_layout(
        xaxis=dict(
            title='Atom',
            tickmode='array',
            tickvals=x,
            ticktext=atom_labels,
            tickangle=45,
            tickfont=dict(size=9)
        ),
        yaxis=dict(
            title='Fukui Index Value',
            showgrid=True,
            gridcolor='lightgray'
        ),
        title=f'Fukui Indices Comparison (Top {top_n} Reactive Atoms)',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        barmode='group',
        legend=dict(x=1.0, y=1.0, xanchor='right')
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved Fukui comparison plot to {output_file}")

    return fig


def create_global_descriptors_plot(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 8),
    dpi: int = 300
) -> go.Figure:
    """
    Create visualization of global reactivity descriptors.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure

    Returns:
        Plotly figure object
    """
    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.4, 0.6],
        specs=[[{"type": "bar"}, {"type": "table"}]],
        subplot_titles=('Global Reactivity Descriptors', 'Top Reactive Sites')
    )

    # Collect global descriptors
    descriptors = []
    values = []
    colors = []

    if fukui_data.ionization_potential is not None:
        descriptors.append('IP<br>(Ionization<br>Potential)')
        values.append(fukui_data.ionization_potential)
        colors.append('#e74c3c')

    if fukui_data.electron_affinity is not None:
        descriptors.append('EA<br>(Electron<br>Affinity)')
        values.append(fukui_data.electron_affinity)
        colors.append('#3498db')

    if fukui_data.electronegativity is not None:
        descriptors.append('χ<br>(Electro-<br>negativity)')
        values.append(fukui_data.electronegativity)
        colors.append('#9b59b6')

    if fukui_data.chemical_hardness is not None:
        descriptors.append('η<br>(Chemical<br>Hardness)')
        values.append(fukui_data.chemical_hardness)
        colors.append('#f39c12')

    if fukui_data.electrophilicity_index is not None:
        descriptors.append('ω<br>(Electro-<br>philicity)')
        values.append(fukui_data.electrophilicity_index)
        colors.append('#e67e22')

    # Plot 1: Global descriptors bar chart
    if descriptors:
        fig.add_trace(go.Bar(
            x=descriptors, y=values,
            marker=dict(color=colors, line=dict(color='black', width=1.5)),
            opacity=0.7,
            text=[f'{v:.3f}' for v in values],
            textposition='outside',
            textfont=dict(size=9),
            showlegend=False
        ), row=1, col=1)

        fig.update_yaxes(title_text='Value (eV)', showgrid=True, gridcolor='lightgray', row=1, col=1)

    # Plot 2: Top reactive sites table
    if fukui_data.atomic_fukui:
        sorted_f_plus = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_plus, reverse=True)[:5]
        sorted_f_minus = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_minus, reverse=True)[:5]

        # Prepare table data
        header = ['Type', 'Atom', 'Index Value']
        cell_values = [
            ['f⁺'] * 5 + ['f⁻'] * 5,
            [f"{atom.atom_index} ({atom.element})" for atom in sorted_f_plus] +
            [f"{atom.atom_index} ({atom.element})" for atom in sorted_f_minus],
            [f"{atom.f_plus:.4f}" for atom in sorted_f_plus] +
            [f"{atom.f_minus:.4f}" for atom in sorted_f_minus]
        ]

        colors_col = ['#ffe6e6'] * 5 + ['#e6f2ff'] * 5

        fig.add_trace(go.Table(
            header=dict(
                values=header,
                fill_color='#4472C4',
                font=dict(color='white', size=11),
                align='center'
            ),
            cells=dict(
                values=cell_values,
                fill_color=[colors_col, 'white', 'white'],
                align='center',
                font=dict(size=9, family='monospace')
            )
        ), row=1, col=2)

    fig.update_layout(
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        showlegend=False
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved global descriptors plot to {output_file}")

    return fig


def create_fukui_heatmap(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 8),
    dpi: int = 300,
    max_atoms: int = 30
) -> go.Figure:
    """
    Create heatmap of Fukui indices across all atoms.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        max_atoms: Maximum number of atoms to display

    Returns:
        Plotly figure object
    """
    if not fukui_data.atomic_fukui:
        raise ValueError("No Fukui data available")

    # Sort by f⁰ and limit
    sorted_atoms = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_zero, reverse=True)[:max_atoms]
    n_atoms = len(sorted_atoms)

    # Prepare data matrix
    data = np.zeros((4, n_atoms))
    for i, atom in enumerate(sorted_atoms):
        data[0, i] = atom.f_plus
        data[1, i] = atom.f_minus
        data[2, i] = atom.f_zero
        data[3, i] = atom.dual_descriptor

    atom_labels = [f"{atom.atom_index}-{atom.element}" for atom in sorted_atoms]
    fukui_labels = ['f⁺<br>(Electrophilic<br>attack)', 'f⁻<br>(Nucleophilic<br>attack)',
                    'f⁰<br>(Radical<br>attack)', 'Δf<br>(Dual<br>descriptor)']

    # Create annotations for values
    annotations = []
    for i in range(4):
        for j in range(n_atoms):
            annotations.append(
                dict(
                    x=j, y=i,
                    text=f'{data[i, j]:.3f}',
                    showarrow=False,
                    font=dict(size=6, color='black')
                )
            )

    fig = go.Figure(data=go.Heatmap(
        z=data,
        x=atom_labels,
        y=fukui_labels,
        colorscale='RdYlBu_r',
        colorbar=dict(title='Fukui Index Value')
    ))

    fig.update_layout(
        title=f'Fukui Indices Heatmap ({fukui_data.method.upper()} Charges)',
        xaxis=dict(title='Atom Index', tickangle=45, tickfont=dict(size=8)),
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        annotations=annotations
    )

    if output_file:
        if output_file.endswith('.html'):
            fig.write_html(output_file)
        else:
            fig.write_image(output_file, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saved Fukui heatmap to {output_file}")

    return fig
