"""
Multi-Dataset Orbital Energy Comparison

Creates side-by-side orbital energy diagrams with connection lines
showing how orbital energies change across different calculations.

Features:
- Side-by-side orbital diagrams
- Dashed connection lines between corresponding orbitals
- HOMO-LUMO gap annotations
- Support for flat and nested dataset lists
- Color-coded occupied/virtual orbitals
- Grouped comparisons (only connect within groups)
- Orbital filtering by energy range
- Orbital filtering by type (occupied/virtual/frontier)

Based on implementation from 0cbz.ipynb.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import List, Dict, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)


def find_homo_lumo(orbital_energies: List[Dict]) -> Tuple[int, int]:
    """
    Find HOMO and LUMO indices from orbital energies.

    Args:
        orbital_energies: List of orbital dicts with 'occupation' and 'energy_ev' keys

    Returns:
        (homo_index, lumo_index) tuple
    """
    homo_idx = None
    lumo_idx = None

    for i, orb in enumerate(orbital_energies):
        occupation = orb.get('occupation', 0.0)
        if occupation > 0.5:  # Occupied
            homo_idx = i
        elif occupation < 0.5 and lumo_idx is None:  # First virtual
            lumo_idx = i
            break

    if homo_idx is None:
        logger.warning("HOMO not found, using highest energy orbital")
        homo_idx = len(orbital_energies) // 2

    if lumo_idx is None:
        logger.warning("LUMO not found, using orbital after HOMO")
        lumo_idx = homo_idx + 1

    return homo_idx, lumo_idx


def filter_orbitals_by_energy(
    orbital_energies: List[Dict],
    min_energy: Optional[float] = None,
    max_energy: Optional[float] = None
) -> List[Dict]:
    """
    Filter orbitals by energy range.

    Args:
        orbital_energies: List of orbital dicts
        min_energy: Minimum energy in eV (None = no lower limit)
        max_energy: Maximum energy in eV (None = no upper limit)

    Returns:
        Filtered list of orbitals
    """
    filtered = []
    for orb in orbital_energies:
        energy = orb.get('energy_ev', 0.0)
        if min_energy is not None and energy < min_energy:
            continue
        if max_energy is not None and energy > max_energy:
            continue
        filtered.append(orb)

    logger.debug(f"Filtered {len(orbital_energies)} orbitals to {len(filtered)} by energy range")
    return filtered


def filter_orbitals_by_type(
    orbital_energies: List[Dict],
    orbital_type: str = 'all'
) -> List[Dict]:
    """
    Filter orbitals by occupation type.

    Args:
        orbital_energies: List of orbital dicts
        orbital_type: 'all', 'occupied', 'virtual', or 'frontier'
                     'frontier' returns HOMO-2 to LUMO+2

    Returns:
        Filtered list of orbitals
    """
    if orbital_type == 'all':
        return orbital_energies

    if orbital_type == 'occupied':
        filtered = [orb for orb in orbital_energies if orb.get('occupation', 0.0) > 0.5]
    elif orbital_type == 'virtual':
        filtered = [orb for orb in orbital_energies if orb.get('occupation', 0.0) < 0.5]
    elif orbital_type == 'frontier':
        homo_idx, lumo_idx = find_homo_lumo(orbital_energies)
        start_idx = max(0, homo_idx - 2)
        end_idx = min(len(orbital_energies), lumo_idx + 3)
        filtered = orbital_energies[start_idx:end_idx]
    else:
        logger.warning(f"Unknown orbital_type '{orbital_type}', returning all")
        return orbital_energies

    logger.debug(f"Filtered to {len(filtered)} {orbital_type} orbitals")
    return filtered


def create_orbital_comparison(
    datasets: Union[List[List[Dict]], List[Dict]],
    labels: List[str],
    n_orbitals: int = 6,
    figsize: Tuple[float, float] = (14, 8),
    bar_width: float = 0.6,
    connection_color: str = 'gray',
    connection_style: str = 'dash',
    connection_alpha: float = 0.4,
    show_gap_values: bool = True,
    energy_range: Optional[Tuple[float, float]] = None,
    orbital_filter: str = 'all',
    title: str = "Orbital Energy Comparison",
    save_path: Optional[str] = None
) -> go.Figure:
    """
    Create multi-dataset orbital energy comparison diagram.

    Args:
        datasets: List of datasets (flat) or list of groups (nested).
                 Each dataset should have 'orbital_energies' key.
                 Flat: [dataset1, dataset2, dataset3] - connects all adjacent
                 Nested: [[dataset1, dataset2], [dataset3, dataset4]] - connects within groups
        labels: Dataset labels matching structure
        n_orbitals: Number of orbitals to show around HOMO-LUMO gap (n/2 above, n/2 below)
        figsize: Figure size
        bar_width: Width of orbital energy bars (0-1)
        connection_color: Color for connection lines
        connection_style: Line style for connections
        connection_alpha: Transparency of connection lines
        show_gap_values: Show HOMO-LUMO gap values
        energy_range: Optional (min_eV, max_eV) tuple to filter by energy
        orbital_filter: 'all', 'occupied', 'virtual', or 'frontier' to filter orbital type
        title: Plot title
        save_path: Save path

    Returns:
        plotly Figure object
    """
    logger.info(f"Creating orbital comparison for {len(labels)} datasets")

    # Determine if datasets is flat or nested
    is_nested = isinstance(datasets[0], (list, tuple))

    if is_nested:
        # Flatten for processing but remember groups
        flat_datasets = []
        flat_labels = []
        groups = []
        for group in datasets:
            group_start = len(flat_datasets)
            for dataset in group:
                flat_datasets.append(dataset)
            groups.append(list(range(group_start, len(flat_datasets))))
        flat_labels = labels
    else:
        flat_datasets = datasets
        flat_labels = labels
        # Single group with all datasets
        groups = [list(range(len(datasets)))]

    n_datasets = len(flat_datasets)

    # Extract orbital energies and find HOMO/LUMO for each dataset
    orbital_data = []
    for i, dataset in enumerate(flat_datasets):
        try:
            orbitals = dataset.get('orbital_energies', [])
            if not orbitals:
                raise ValueError(f"No orbital_energies in dataset {i}")

            # Apply filters if specified
            if energy_range is not None:
                orbitals = filter_orbitals_by_energy(orbitals, energy_range[0], energy_range[1])
            if orbital_filter != 'all':
                orbitals = filter_orbitals_by_type(orbitals, orbital_filter)

            if not orbitals:
                raise ValueError(f"No orbitals remaining after filtering in dataset {i}")

            homo_idx, lumo_idx = find_homo_lumo(orbitals)

            # Extract n_orbitals/2 around gap
            n_below = n_orbitals // 2
            n_above = n_orbitals - n_below

            start_idx = max(0, homo_idx - n_below + 1)
            end_idx = min(len(orbitals), lumo_idx + n_above)

            selected_orbitals = orbitals[start_idx:end_idx]

            # Get energies in eV
            energies = [orb.get('energy_ev', 0.0) for orb in selected_orbitals]
            occupations = [orb.get('occupation', 0.0) for orb in selected_orbitals]

            homo_energy = orbitals[homo_idx].get('energy_ev', 0.0)
            lumo_energy = orbitals[lumo_idx].get('energy_ev', 0.0)
            gap = lumo_energy - homo_energy

            orbital_data.append({
                'energies': energies,
                'occupations': occupations,
                'homo_energy': homo_energy,
                'lumo_energy': lumo_energy,
                'gap': gap,
                'homo_idx_in_selection': homo_idx - start_idx,
                'lumo_idx_in_selection': lumo_idx - start_idx
            })

            logger.debug(f"Dataset {i}: HOMO={homo_energy:.3f} eV, LUMO={lumo_energy:.3f} eV, Gap={gap:.3f} eV")

        except Exception as e:
            logger.error(f"Error processing dataset {i}: {e}")
            raise

    # Create figure
    fig = go.Figure()

    # X positions for datasets
    x_positions = np.arange(n_datasets) * 2  # Space them out

    # Plot orbital levels for each dataset
    for i, (data, x_pos, label) in enumerate(zip(orbital_data, x_positions, flat_labels)):
        energies = data['energies']
        occupations = data['occupations']

        for j, (energy, occupation) in enumerate(zip(energies, occupations)):
            # Color based on occupation
            color = 'blue' if occupation > 0.5 else 'red'

            # Draw horizontal bar (simulated with scatter)
            fig.add_trace(go.Scatter(
                x=[x_pos - bar_width/2, x_pos + bar_width/2],
                y=[energy, energy],
                mode='lines',
                line=dict(color=color, width=2.5),
                showlegend=False,
                hoverinfo='text',
                hovertext=f'{label}<br>Energy: {energy:.3f} eV'
            ))

        # Add HOMO-LUMO gap arrow and label
        if show_gap_values:
            gap = data['gap']
            homo_e = data['homo_energy']
            lumo_e = data['lumo_energy']

            # Add annotation for gap
            fig.add_annotation(
                x=x_pos + bar_width/2 + 0.5,
                y=(homo_e + lumo_e) / 2,
                text=f'{gap:.2f} eV',
                showarrow=True,
                arrowhead=3,
                arrowsize=1,
                arrowwidth=1.5,
                arrowcolor='black',
                ax=0,
                ay=(lumo_e - homo_e) * 30,
                bgcolor='white',
                bordercolor='black',
                borderwidth=1
            )

    # Draw connection lines between corresponding orbitals
    for group in groups:
        for i in range(len(group) - 1):
            idx1 = group[i]
            idx2 = group[i + 1]

            data1 = orbital_data[idx1]
            data2 = orbital_data[idx2]

            x1 = x_positions[idx1]
            x2 = x_positions[idx2]

            # Connect corresponding orbital levels
            n_to_connect = min(len(data1['energies']), len(data2['energies']))

            for j in range(n_to_connect):
                e1 = data1['energies'][j]
                e2 = data2['energies'][j]

                fig.add_trace(go.Scatter(
                    x=[x1 + bar_width/2, x2 - bar_width/2],
                    y=[e1, e2],
                    mode='lines',
                    line=dict(color=connection_color, dash=connection_style, width=1.0),
                    opacity=connection_alpha,
                    showlegend=False,
                    hoverinfo='skip'
                ))

    # Add legend manually
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='blue', width=2.5),
        name='Occupied'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='red', width=2.5),
        name='Virtual'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color=connection_color, dash=connection_style, width=1.0),
        opacity=connection_alpha,
        name='Connection'
    ))

    # Update layout
    fig.update_layout(
        title=title,
        xaxis_title='Dataset',
        yaxis_title='Orbital Energy (eV)',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        xaxis=dict(
            tickmode='array',
            tickvals=x_positions,
            ticktext=flat_labels,
            tickangle=45
        ),
        yaxis=dict(showgrid=True, gridcolor='lightgray'),
        showlegend=True
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)
        logger.info(f"Saving figure to {save_path}")

    logger.info("Orbital comparison created successfully")
    return fig


def create_simple_orbital_diagram(
    orbital_energies: List[Dict],
    n_orbitals: int = 10,
    figsize: Tuple[float, float] = (6, 8),
    title: str = "Molecular Orbital Diagram",
    save_path: Optional[str] = None
) -> go.Figure:
    """
    Create a simple molecular orbital energy diagram for a single dataset.

    Args:
        orbital_energies: List of orbital dicts
        n_orbitals: Number of orbitals to show around HOMO-LUMO
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        plotly Figure object
    """
    # Find HOMO/LUMO
    homo_idx, lumo_idx = find_homo_lumo(orbital_energies)

    # Select orbitals around gap
    n_below = n_orbitals // 2
    n_above = n_orbitals - n_below

    start_idx = max(0, homo_idx - n_below + 1)
    end_idx = min(len(orbital_energies), lumo_idx + n_above)

    selected_orbitals = orbital_energies[start_idx:end_idx]

    # Extract energies
    energies = [orb.get('energy_ev', 0.0) for orb in selected_orbitals]
    occupations = [orb.get('occupation', 0.0) for orb in selected_orbitals]
    indices = [orb.get('index', start_idx + i) for i, orb in enumerate(selected_orbitals)]

    # Create figure
    fig = go.Figure()

    # Plot orbital levels
    for i, (energy, occupation, idx) in enumerate(zip(energies, occupations, indices)):
        color = 'blue' if occupation > 0.5 else 'red'

        # Draw horizontal line
        fig.add_trace(go.Scatter(
            x=[0.3, 0.7],
            y=[energy, energy],
            mode='lines',
            line=dict(color=color, width=3),
            showlegend=False,
            hoverinfo='text',
            hovertext=f'MO {idx}<br>Energy: {energy:.3f} eV'
        ))

        # Add orbital index
        fig.add_annotation(
            x=0.75, y=energy,
            text=f'MO {idx}',
            showarrow=False,
            xanchor='left',
            font=dict(size=9)
        )

    # Mark HOMO and LUMO
    homo_energy = orbital_energies[homo_idx].get('energy_ev', 0.0)
    lumo_energy = orbital_energies[lumo_idx].get('energy_ev', 0.0)
    gap = lumo_energy - homo_energy

    fig.add_hline(y=homo_energy, line=dict(color='blue', dash='dash', width=1), opacity=0.5)
    fig.add_hline(y=lumo_energy, line=dict(color='red', dash='dash', width=1), opacity=0.5)

    # Add HOMO/LUMO labels
    fig.add_annotation(x=0.15, y=homo_energy, text='HOMO', showarrow=False,
                      font=dict(size=10, color='blue'), yshift=-10)
    fig.add_annotation(x=0.15, y=lumo_energy, text='LUMO', showarrow=False,
                      font=dict(size=10, color='red'), yshift=10)

    # Add gap annotation
    fig.add_annotation(
        x=0.1, y=(homo_energy + lumo_energy) / 2,
        text=f'{gap:.2f} eV',
        showarrow=True,
        arrowhead=3,
        arrowsize=1,
        arrowwidth=1.5,
        arrowcolor='black',
        ax=0,
        ay=(lumo_energy - homo_energy) * 30,
        font=dict(size=10)
    )

    # Add legend
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='blue', width=3),
        name='Occupied'
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='red', width=3),
        name='Virtual'
    ))

    # Update layout
    fig.update_layout(
        title=title,
        yaxis_title='Orbital Energy (eV)',
        template='plotly_white',
        width=figsize[0] * 80,
        height=figsize[1] * 80,
        xaxis=dict(range=[0, 1], showticklabels=False, showgrid=False),
        yaxis=dict(showgrid=True, gridcolor='lightgray'),
        showlegend=True
    )

    if save_path:
        if save_path.endswith('.html'):
            fig.write_html(save_path)
        else:
            fig.write_image(save_path, width=figsize[0] * 80, height=figsize[1] * 80, scale=3)

    return fig


# Export main functions
__all__ = [
    'create_orbital_comparison',
    'create_simple_orbital_diagram',
    'find_homo_lumo',
    'filter_orbitals_by_energy',
    'filter_orbitals_by_type'
]
