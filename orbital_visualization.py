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
import matplotlib.pyplot as plt
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
    connection_style: str = '--',
    connection_alpha: float = 0.4,
    show_gap_values: bool = True,
    energy_range: Optional[Tuple[float, float]] = None,
    orbital_filter: str = 'all',
    title: str = "Orbital Energy Comparison",
    save_path: Optional[str] = None
) -> plt.Figure:
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
        matplotlib Figure object
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
    fig, ax = plt.subplots(figsize=figsize)

    # X positions for datasets
    x_positions = np.arange(n_datasets) * 2  # Space them out

    # Plot orbital levels for each dataset
    for i, (data, x_pos, label) in enumerate(zip(orbital_data, x_positions, flat_labels)):
        energies = data['energies']
        occupations = data['occupations']

        for j, (energy, occupation) in enumerate(zip(energies, occupations)):
            # Color based on occupation
            color = 'blue' if occupation > 0.5 else 'red'

            # Draw horizontal bar
            ax.hlines(energy, x_pos - bar_width/2, x_pos + bar_width/2,
                     colors=color, linewidth=2.5)

        # Add HOMO-LUMO gap arrow and label
        if show_gap_values:
            gap = data['gap']
            homo_e = data['homo_energy']
            lumo_e = data['lumo_energy']

            # Arrow
            ax.annotate('', xy=(x_pos + bar_width/2 + 0.2, lumo_e),
                       xytext=(x_pos + bar_width/2 + 0.2, homo_e),
                       arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))

            # Gap value
            gap_text = f'{gap:.2f} eV'
            ax.text(x_pos + bar_width/2 + 0.5, (homo_e + lumo_e) / 2, gap_text,
                   fontsize=9, va='center',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black'))

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

                ax.plot([x1 + bar_width/2, x2 - bar_width/2], [e1, e2],
                       linestyle=connection_style, color=connection_color,
                       alpha=connection_alpha, linewidth=1.0, zorder=0)

    # Formatting
    ax.set_xlabel('Dataset', fontsize=12, fontweight='bold')
    ax.set_ylabel('Orbital Energy (eV)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Set x-axis labels
    ax.set_xticks(x_positions)
    ax.set_xticklabels(flat_labels, fontsize=10, rotation=45, ha='right')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', label='Occupied'),
        Patch(facecolor='red', label='Virtual'),
        plt.Line2D([0], [0], linestyle=connection_style, color=connection_color,
                   alpha=connection_alpha, label='Connection')
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=10)

    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()

    if save_path:
        logger.info(f"Saving figure to {save_path}")
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    logger.info("Orbital comparison created successfully")
    return fig


def create_simple_orbital_diagram(
    orbital_energies: List[Dict],
    n_orbitals: int = 10,
    figsize: Tuple[float, float] = (6, 8),
    title: str = "Molecular Orbital Diagram",
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create a simple molecular orbital energy diagram for a single dataset.

    Args:
        orbital_energies: List of orbital dicts
        n_orbitals: Number of orbitals to show around HOMO-LUMO
        figsize: Figure size
        title: Plot title
        save_path: Save path

    Returns:
        matplotlib Figure object
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
    fig, ax = plt.subplots(figsize=figsize)

    # Plot orbital levels
    for i, (energy, occupation, idx) in enumerate(zip(energies, occupations, indices)):
        color = 'blue' if occupation > 0.5 else 'red'
        ax.hlines(energy, 0.3, 0.7, colors=color, linewidth=3)

        # Add orbital index
        ax.text(0.75, energy, f'MO {idx}', fontsize=9, va='center')

    # Mark HOMO and LUMO
    homo_energy = orbital_energies[homo_idx].get('energy_ev', 0.0)
    lumo_energy = orbital_energies[lumo_idx].get('energy_ev', 0.0)
    gap = lumo_energy - homo_energy

    ax.axhline(homo_energy, color='blue', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(lumo_energy, color='red', linestyle='--', alpha=0.5, linewidth=1)

    ax.text(0.15, homo_energy, 'HOMO', fontsize=10, fontweight='bold', va='bottom')
    ax.text(0.15, lumo_energy, 'LUMO', fontsize=10, fontweight='bold', va='top')

    # Add gap annotation
    ax.annotate('', xy=(0.2, lumo_energy), xytext=(0.2, homo_energy),
               arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(0.1, (homo_energy + lumo_energy) / 2, f'{gap:.2f} eV',
           fontsize=10, va='center', ha='right', fontweight='bold')

    # Formatting
    ax.set_ylabel('Orbital Energy (eV)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1)
    ax.set_xticks([])

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', label='Occupied'),
        Patch(facecolor='red', label='Virtual')
    ]
    ax.legend(handles=legend_elements, loc='best')

    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig


# Export main functions
__all__ = [
    'create_orbital_comparison',
    'create_simple_orbital_diagram',
    'find_homo_lumo',
    'filter_orbitals_by_energy',
    'filter_orbitals_by_type'
]
