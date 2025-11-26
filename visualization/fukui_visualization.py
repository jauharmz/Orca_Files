"""
Fukui Function Visualization

Create visualizations for Fukui function analysis including:
- Atomic Fukui indices (f⁺, f⁻, f⁰) as bar plots
- Dual descriptor plots
- Reactivity site identification
- Global reactivity descriptor summaries
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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
) -> plt.Figure:
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
        Matplotlib figure object
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
    width = 0.2

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)

    # Plot 1: f⁺ (nucleophilic Fukui - electrophilic attack sites)
    bars1 = ax1.bar(x, f_plus, width=0.6, color='#e74c3c', alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Atom Index', fontsize=11, fontweight='bold')
    ax1.set_ylabel('f⁺ (Nucleophilic Fukui)', fontsize=11, fontweight='bold')
    ax1.set_title('Electrophilic Attack Sites (High f⁺)', fontsize=12, fontweight='bold')
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax1.grid(True, alpha=0.3, axis='y')

    if show_labels:
        ax1.set_xticks(x)
        ax1.set_xticklabels([f"{idx}\n{elem}" for idx, elem in zip(atom_indices, elements)],
                           rotation=0, fontsize=8)
    else:
        ax1.set_xticks(x)
        ax1.set_xticklabels(atom_indices, rotation=45)

    # Highlight top 3
    top_3_indices = np.argsort(f_plus)[-3:]
    for idx in top_3_indices:
        bars1[idx].set_color('#c0392b')
        bars1[idx].set_linewidth(2)

    # Plot 2: f⁻ (electrophilic Fukui - nucleophilic attack sites)
    bars2 = ax2.bar(x, f_minus, width=0.6, color='#3498db', alpha=0.7, edgecolor='black')
    ax2.set_xlabel('Atom Index', fontsize=11, fontweight='bold')
    ax2.set_ylabel('f⁻ (Electrophilic Fukui)', fontsize=11, fontweight='bold')
    ax2.set_title('Nucleophilic Attack Sites (High f⁻)', fontsize=12, fontweight='bold')
    ax2.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax2.grid(True, alpha=0.3, axis='y')

    if show_labels:
        ax2.set_xticks(x)
        ax2.set_xticklabels([f"{idx}\n{elem}" for idx, elem in zip(atom_indices, elements)],
                           rotation=0, fontsize=8)
    else:
        ax2.set_xticks(x)
        ax2.set_xticklabels(atom_indices, rotation=45)

    # Highlight top 3
    top_3_indices = np.argsort(f_minus)[-3:]
    for idx in top_3_indices:
        bars2[idx].set_color('#2980b9')
        bars2[idx].set_linewidth(2)

    # Plot 3: f⁰ (radical Fukui - radical attack sites)
    bars3 = ax3.bar(x, f_zero, width=0.6, color='#2ecc71', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Atom Index', fontsize=11, fontweight='bold')
    ax3.set_ylabel('f⁰ (Radical Fukui)', fontsize=11, fontweight='bold')
    ax3.set_title('Radical Attack Sites (High f⁰)', fontsize=12, fontweight='bold')
    ax3.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax3.grid(True, alpha=0.3, axis='y')

    if show_labels:
        ax3.set_xticks(x)
        ax3.set_xticklabels([f"{idx}\n{elem}" for idx, elem in zip(atom_indices, elements)],
                           rotation=0, fontsize=8)
    else:
        ax3.set_xticks(x)
        ax3.set_xticklabels(atom_indices, rotation=45)

    # Highlight top 3
    top_3_indices = np.argsort(f_zero)[-3:]
    for idx in top_3_indices:
        bars3[idx].set_color('#27ae60')
        bars3[idx].set_linewidth(2)

    # Plot 4: Dual descriptor (Δf = f⁻ - f⁺)
    colors = ['#3498db' if d > 0 else '#e74c3c' for d in dual]
    bars4 = ax4.bar(x, dual, width=0.6, color=colors, alpha=0.7, edgecolor='black')
    ax4.set_xlabel('Atom Index', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Δf (Dual Descriptor)', fontsize=11, fontweight='bold')
    ax4.set_title('Dual Descriptor (Δf = f⁻ - f⁺)', fontsize=12, fontweight='bold')
    ax4.axhline(y=0, color='black', linestyle='-', linewidth=1.5)
    ax4.grid(True, alpha=0.3, axis='y')

    if show_labels:
        ax4.set_xticks(x)
        ax4.set_xticklabels([f"{idx}\n{elem}" for idx, elem in zip(atom_indices, elements)],
                           rotation=0, fontsize=8)
    else:
        ax4.set_xticks(x)
        ax4.set_xticklabels(atom_indices, rotation=45)

    # Add legend to dual descriptor
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#3498db', alpha=0.7, label='Nucleophilic (Δf > 0)'),
        Patch(facecolor='#e74c3c', alpha=0.7, label='Electrophilic (Δf < 0)')
    ]
    ax4.legend(handles=legend_elements, loc='upper right', fontsize=9)

    # Add global title with method
    fig.suptitle(f'Fukui Function Analysis ({fukui_data.method.upper()} Charges)',
                 fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout(rect=[0, 0, 1, 0.99])

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved Fukui bar plot to {output_file}")

    return fig


def create_fukui_comparison_plot(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 6),
    dpi: int = 300,
    top_n: int = 15
) -> plt.Figure:
    """
    Create side-by-side comparison of all three Fukui indices.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        top_n: Show top N atoms by f⁰

    Returns:
        Matplotlib figure object
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

    fig, ax = plt.subplots(figsize=figsize)

    # Plot all three Fukui indices
    bars1 = ax.bar(x - width, f_plus, width, label='f⁺ (Electrophilic attack)',
                   color='#e74c3c', alpha=0.7, edgecolor='black')
    bars2 = ax.bar(x, f_minus, width, label='f⁻ (Nucleophilic attack)',
                   color='#3498db', alpha=0.7, edgecolor='black')
    bars3 = ax.bar(x + width, f_zero, width, label='f⁰ (Radical attack)',
                   color='#2ecc71', alpha=0.7, edgecolor='black')

    ax.set_xlabel('Atom', fontsize=12, fontweight='bold')
    ax.set_ylabel('Fukui Index Value', fontsize=12, fontweight='bold')
    ax.set_title(f'Fukui Indices Comparison (Top {top_n} Reactive Atoms)',
                 fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(atom_labels, rotation=45, ha='right', fontsize=9)
    ax.legend(fontsize=10, loc='upper right')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved Fukui comparison plot to {output_file}")

    return fig


def create_global_descriptors_plot(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (12, 8),
    dpi: int = 300
) -> plt.Figure:
    """
    Create visualization of global reactivity descriptors.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure

    Returns:
        Matplotlib figure object
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Collect global descriptors
    descriptors = []
    values = []
    colors = []

    if fukui_data.ionization_potential is not None:
        descriptors.append('IP\n(Ionization\nPotential)')
        values.append(fukui_data.ionization_potential)
        colors.append('#e74c3c')

    if fukui_data.electron_affinity is not None:
        descriptors.append('EA\n(Electron\nAffinity)')
        values.append(fukui_data.electron_affinity)
        colors.append('#3498db')

    if fukui_data.electronegativity is not None:
        descriptors.append('χ\n(Electro-\nnegativity)')
        values.append(fukui_data.electronegativity)
        colors.append('#9b59b6')

    if fukui_data.chemical_hardness is not None:
        descriptors.append('η\n(Chemical\nHardness)')
        values.append(fukui_data.chemical_hardness)
        colors.append('#f39c12')

    if fukui_data.electrophilicity_index is not None:
        descriptors.append('ω\n(Electro-\nphilicity)')
        values.append(fukui_data.electrophilicity_index)
        colors.append('#e67e22')

    # Plot 1: Global descriptors bar chart
    if descriptors:
        x = np.arange(len(descriptors))
        bars = ax1.bar(x, values, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax1.set_ylabel('Value (eV)', fontsize=12, fontweight='bold')
        ax1.set_title('Global Reactivity Descriptors', fontsize=13, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(descriptors, fontsize=10)
        ax1.grid(True, alpha=0.3, axis='y')

        # Add value labels on bars
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2., height,
                    f'{val:.3f}',
                    ha='center', va='bottom', fontsize=9, fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'No global descriptors available',
                ha='center', va='center', transform=ax1.transAxes, fontsize=12)
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)

    # Plot 2: Top reactive sites summary
    if fukui_data.atomic_fukui:
        top_f_plus = fukui_data.get_max_f_plus()
        top_f_minus = fukui_data.get_max_f_minus()

        # Get top 5 for each
        sorted_f_plus = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_plus, reverse=True)[:5]
        sorted_f_minus = sorted(fukui_data.atomic_fukui, key=lambda x: x.f_minus, reverse=True)[:5]

        ax2.axis('off')

        # Title
        ax2.text(0.5, 0.95, 'Top Reactive Sites', ha='center', va='top',
                fontsize=13, fontweight='bold', transform=ax2.transAxes)

        # Electrophilic attack sites (f⁺)
        y_pos = 0.82
        ax2.text(0.05, y_pos, 'Electrophilic Attack Sites (f⁺):', ha='left', va='top',
                fontsize=11, fontweight='bold', color='#e74c3c', transform=ax2.transAxes)
        y_pos -= 0.06
        for i, atom in enumerate(sorted_f_plus, 1):
            text = f"{i}. Atom {atom.atom_index} ({atom.element}): f⁺ = {atom.f_plus:.4f}"
            ax2.text(0.08, y_pos, text, ha='left', va='top',
                    fontsize=9, transform=ax2.transAxes, family='monospace')
            y_pos -= 0.05

        # Nucleophilic attack sites (f⁻)
        y_pos -= 0.04
        ax2.text(0.05, y_pos, 'Nucleophilic Attack Sites (f⁻):', ha='left', va='top',
                fontsize=11, fontweight='bold', color='#3498db', transform=ax2.transAxes)
        y_pos -= 0.06
        for i, atom in enumerate(sorted_f_minus, 1):
            text = f"{i}. Atom {atom.atom_index} ({atom.element}): f⁻ = {atom.f_minus:.4f}"
            ax2.text(0.08, y_pos, text, ha='left', va='top',
                    fontsize=9, transform=ax2.transAxes, family='monospace')
            y_pos -= 0.05

        # Add interpretation guide
        y_pos -= 0.04
        ax2.text(0.05, y_pos, 'Interpretation:', ha='left', va='top',
                fontsize=10, fontweight='bold', style='italic', transform=ax2.transAxes)
        y_pos -= 0.05
        guide_text = [
            "• High f⁺: Susceptible to electrophilic attack",
            "• High f⁻: Susceptible to nucleophilic attack",
            "• For degradation: High f⁺ or f⁻ = reactive sites"
        ]
        for line in guide_text:
            ax2.text(0.08, y_pos, line, ha='left', va='top',
                    fontsize=8, transform=ax2.transAxes)
            y_pos -= 0.04

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved global descriptors plot to {output_file}")

    return fig


def create_fukui_heatmap(
    fukui_data,
    output_file: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 8),
    dpi: int = 300,
    max_atoms: int = 30
) -> plt.Figure:
    """
    Create heatmap of Fukui indices across all atoms.

    Args:
        fukui_data: FukuiIndices object
        output_file: Path to save figure (optional)
        figsize: Figure size (width, height)
        dpi: Resolution for saved figure
        max_atoms: Maximum number of atoms to display

    Returns:
        Matplotlib figure object
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
    fukui_labels = ['f⁺\n(Electrophilic\nattack)', 'f⁻\n(Nucleophilic\nattack)',
                    'f⁰\n(Radical\nattack)', 'Δf\n(Dual\ndescriptor)']

    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap
    im = ax.imshow(data, cmap='RdYlBu_r', aspect='auto', interpolation='nearest')

    # Set ticks
    ax.set_xticks(np.arange(n_atoms))
    ax.set_yticks(np.arange(4))
    ax.set_xticklabels(atom_labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(fukui_labels, fontsize=10)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Fukui Index Value', rotation=270, labelpad=20, fontsize=11, fontweight='bold')

    # Add value annotations
    for i in range(4):
        for j in range(n_atoms):
            text = ax.text(j, i, f'{data[i, j]:.3f}',
                          ha="center", va="center", color="black", fontsize=6)

    ax.set_title(f'Fukui Indices Heatmap ({fukui_data.method.upper()} Charges)',
                 fontsize=13, fontweight='bold', pad=10)
    ax.set_xlabel('Atom Index', fontsize=11, fontweight='bold')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
        logger.info(f"Saved Fukui heatmap to {output_file}")

    return fig
