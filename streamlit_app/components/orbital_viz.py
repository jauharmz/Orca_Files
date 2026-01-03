"""
Orbital Visualization Component - Energy levels and diagram

Features:
- Orbital energy level diagram
- HOMO/LUMO highlighting
- Spin-up/down separation
- N-orbital slider
- Multi-molecule comparison
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List


def render_orbital_viz(df: pd.DataFrame):
    """Render orbital visualization tab."""
    
    st.header("🔮 Orbital Analysis")
    
    # Molecule selector
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    selected_mol = st.selectbox("Select Molecule", mol_ids, key="orbital_mol_select")
    
    mol_row = df[df["molecule_id"] == selected_mol]
    if mol_row.empty:
        st.warning("Select a molecule")
        return
    
    mol_data = mol_row.iloc[0]
    orbitals = mol_data.get("orbitals")
    
    if orbitals is None or (hasattr(orbitals, 'empty') and orbitals.empty):
        st.warning("No orbital data available for this molecule")
        return
    
    # Tabs
    tabs = st.tabs(["📊 Energy Levels", "📈 Comparison"])
    
    with tabs[0]:
        render_orbital_levels(orbitals, selected_mol, mol_data)
    
    with tabs[1]:
        render_orbital_comparison(df, mol_ids)


def render_orbital_levels(orbitals: pd.DataFrame, mol_id: str, mol_data: pd.Series):
    """Render orbital energy level diagram."""
    
    st.subheader(f"📊 Orbital Energy Levels - {mol_id}")
    
    # Customization
    with st.expander("⚙️ Customization", expanded=True):
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            n_orbitals = st.slider("Orbitals Around HOMO", 3, 20, 8, key="n_orb")
        with col2:
            show_labels = st.checkbox("Show Labels", True, key="orb_labels")
        with col3:
            energy_range = st.slider("Energy Range (eV)", -20.0, 10.0, (-15.0, 5.0), key="orb_range")
        with col4:
            separate_spin = st.checkbox("Separate Spin", False, key="orb_spin")
    
    # Get HOMO/LUMO info
    homo_energy = mol_data.get("homo_energy")
    lumo_energy = mol_data.get("lumo_energy")
    
    # Prepare orbital data
    if "eV" not in orbitals.columns and "Eh" in orbitals.columns:
        orbitals = orbitals.copy()
        orbitals["eV"] = orbitals["Eh"] * 27.2114
    
    energy_col = "eV" if "eV" in orbitals.columns else orbitals.columns[1]
    
    # Filter by energy range
    orbitals_filtered = orbitals[
        (orbitals[energy_col] >= energy_range[0]) & 
        (orbitals[energy_col] <= energy_range[1])
    ].copy()
    
    # Find HOMO index
    if "OCC" in orbitals.columns:
        occupied = orbitals[orbitals["OCC"] > 0.5]
        if not occupied.empty:
            homo_idx = occupied[energy_col].idxmax()
        else:
            homo_idx = None
    else:
        homo_idx = None
    
    # Create plot
    fig = go.Figure()
    
    # Group by spin if available and requested
    if separate_spin and "spin" in orbitals_filtered.columns:
        spins = orbitals_filtered["spin"].unique()
        x_offset = {"up": -0.2, "down": 0.2, "alpha": -0.2, "beta": 0.2}
        
        for spin in spins:
            spin_data = orbitals_filtered[orbitals_filtered["spin"] == spin]
            offset = x_offset.get(str(spin).lower(), 0)
            
            for idx, row in spin_data.iterrows():
                energy = row[energy_col]
                is_homo = homo_energy and abs(energy - homo_energy) < 0.01
                is_lumo = lumo_energy and abs(energy - lumo_energy) < 0.01
                
                color = "red" if is_lumo else ("blue" if is_homo else "gray")
                
                fig.add_trace(go.Scatter(
                    x=[offset - 0.15, offset + 0.15],
                    y=[energy, energy],
                    mode="lines",
                    line=dict(color=color, width=3),
                    showlegend=False,
                    hovertemplate=f"Orbital {idx}<br>Energy: {energy:.4f} eV<br>Spin: {spin}<extra></extra>"
                ))
    else:
        for i, (idx, row) in enumerate(orbitals_filtered.iterrows()):
            energy = row[energy_col]
            is_homo = homo_energy and abs(energy - homo_energy) < 0.01
            is_lumo = lumo_energy and abs(energy - lumo_energy) < 0.01
            
            color = "red" if is_lumo else ("blue" if is_homo else "gray")
            width = 4 if (is_homo or is_lumo) else 2
            
            fig.add_trace(go.Scatter(
                x=[-0.3, 0.3],
                y=[energy, energy],
                mode="lines",
                line=dict(color=color, width=width),
                showlegend=False,
                hovertemplate=f"Orbital {idx}<br>Energy: {energy:.4f} eV<extra></extra>"
            ))
            
            if show_labels and (is_homo or is_lumo):
                label = "HOMO" if is_homo else "LUMO"
                fig.add_annotation(
                    x=0.4,
                    y=energy,
                    text=label,
                    showarrow=False,
                    font=dict(size=12, color=color)
                )
    
    # Add HOMO-LUMO gap annotation
    if homo_energy and lumo_energy:
        gap = lumo_energy - homo_energy
        fig.add_annotation(
            x=0.6,
            y=(homo_energy + lumo_energy) / 2,
            text=f"Gap: {gap:.2f} eV",
            showarrow=False,
            font=dict(size=11)
        )
    
    fig.update_layout(
        title="Orbital Energy Levels",
        xaxis=dict(showticklabels=False, range=[-1, 1]),
        yaxis_title="Energy (eV)",
        yaxis=dict(range=[energy_range[0], energy_range[1]]),
        height=600
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 View Orbital Data"):
        from components.data_table import render_simple_table
        render_simple_table(orbitals_filtered, "Orbital Data")


def render_orbital_comparison(df: pd.DataFrame, mol_ids: List[str]):
    """Render orbital comparison across molecules."""
    
    st.subheader("📈 Orbital Comparison")
    
    selected = st.multiselect(
        "Select Molecules to Compare",
        mol_ids,
        default=mol_ids[:min(5, len(mol_ids))],
        key="orbital_compare_select"
    )
    
    if not selected:
        st.warning("Select molecules to compare")
        return
    
    subset = df[df["molecule_id"].isin(selected)]
    subset = subset[subset["homo_energy"].notna()]
    
    if subset.empty:
        st.warning("No orbital data for selected molecules")
        return
    
    # Create comparison chart
    fig = go.Figure()
    
    for i, (_, row) in enumerate(subset.iterrows()):
        mol_id = row["molecule_id"]
        homo = row.get("homo_energy", 0)
        lumo = row.get("lumo_energy", 0)
        
        # HOMO level
        fig.add_trace(go.Scatter(
            x=[i - 0.2, i + 0.2],
            y=[homo, homo],
            mode="lines",
            line=dict(color="blue", width=4),
            name=f"{mol_id} HOMO" if i == 0 else None,
            showlegend=(i == 0),
            hovertemplate=f"<b>{mol_id}</b><br>HOMO: {homo:.3f} eV<extra></extra>"
        ))
        
        # LUMO level
        fig.add_trace(go.Scatter(
            x=[i - 0.2, i + 0.2],
            y=[lumo, lumo],
            mode="lines",
            line=dict(color="red", width=4),
            name=f"{mol_id} LUMO" if i == 0 else None,
            showlegend=(i == 0),
            hovertemplate=f"<b>{mol_id}</b><br>LUMO: {lumo:.3f} eV<extra></extra>"
        ))
        
        # Gap line
        fig.add_trace(go.Scatter(
            x=[i, i],
            y=[homo, lumo],
            mode="lines",
            line=dict(color="gray", width=1, dash="dash"),
            showlegend=False
        ))
    
    fig.update_layout(
        title="HOMO-LUMO Comparison",
        xaxis=dict(
            tickmode="array",
            tickvals=list(range(len(subset))),
            ticktext=subset["molecule_id"].tolist()
        ),
        yaxis_title="Energy (eV)",
        legend=dict(x=0.9, y=0.95)
    )
    
    st.plotly_chart(fig, use_container_width=True)
