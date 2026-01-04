"""
Orbital Visualization Component - Complete Implementation

Features:
- Multi-molecule orbital comparison with state labels
- Orbital energy level diagram
- HOMO/LUMO highlighting
- Occupied vs virtual distinction
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List


def render_orbital_viz(df: pd.DataFrame):
    """Render orbital visualization with molecule+state filter."""
    
    st.subheader("🔮 Orbital Analysis")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build molecule options
    mol_options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        mol_options.append({"label": label, "idx": idx})
    
    unique_labels = list(dict.fromkeys([m["label"] for m in mol_options]))
    
    # Selection
    col1, col2 = st.columns([4, 2])
    with col1:
        selected = st.multiselect(
            "Select Molecules",
            unique_labels,
            default=unique_labels[:min(4, len(unique_labels))],
            key="orbital_mol_select"
        )
    with col2:
        n_orbitals = st.slider("Orbitals Around HOMO", 3, 20, 10, key="n_orbitals")
    
    if not selected:
        st.warning("Select molecules to view")
        return
    
    # Get selected data
    selected_indices = [m["idx"] for m in mol_options if m["label"] in selected]
    selected_df = df.loc[selected_indices]
    
    render_orbital_levels(selected_df, n_orbitals)


def render_orbital_levels(df: pd.DataFrame, n_orbitals: int):
    """Render orbital energy level diagram."""
    
    fig = go.Figure()
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    x_offset = 0
    all_labels = []
    all_data = []
    
    for df_idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        
        orbitals = row.get("orbitals")
        homo_energy = row.get("homo_energy")
        
        if orbitals is None or (hasattr(orbitals, 'empty') and orbitals.empty):
            continue
        
        # Get orbital energies
        if "energy" in orbitals.columns:
            energies = orbitals["energy"].dropna().sort_values()
        elif "Energy" in orbitals.columns:
            energies = orbitals["Energy"].dropna().sort_values()
        else:
            continue
        
        if len(energies) == 0:
            continue
        
        # Filter to N orbitals around HOMO
        if homo_energy is not None and not pd.isna(homo_energy):
            diff = np.abs(energies.values - homo_energy)
            sorted_indices = np.argsort(diff)
            closest_idx = sorted_indices[:min(n_orbitals * 2, len(sorted_indices))]
            energies = energies.iloc[closest_idx].sort_values()
        else:
            # Just take last N*2 orbitals
            energies = energies.tail(n_orbitals * 2)
        
        color_idx = x_offset % len(colors)
        color = colors[color_idx]
        
        for i, e in enumerate(energies.values):
            is_occupied = homo_energy is not None and e <= homo_energy
            is_homo = homo_energy is not None and abs(e - homo_energy) < 0.01
            is_lumo = homo_energy is not None and e > homo_energy and i == len([x for x in energies.values if x <= homo_energy])
            
            line_width = 5 if (is_homo or is_lumo) else 2
            line_style = "solid" if is_occupied else "dash"
            
            fig.add_trace(go.Scatter(
                x=[x_offset - 0.3, x_offset + 0.3],
                y=[e, e],
                mode="lines",
                line=dict(color=color, width=line_width, dash=line_style),
                showlegend=False,
                hovertemplate=f'{label}<br>Energy: {e:.4f} eV<br>{"HOMO" if is_homo else "LUMO" if is_lumo else "Occupied" if is_occupied else "Virtual"}<extra></extra>'
            ))
            
            # Add HOMO/LUMO labels
            if is_homo:
                fig.add_annotation(
                    x=x_offset + 0.35, y=e,
                    text="HOMO",
                    showarrow=False,
                    font=dict(size=9, color=color)
                )
            if is_lumo:
                fig.add_annotation(
                    x=x_offset + 0.35, y=e,
                    text="LUMO",
                    showarrow=False,
                    font=dict(size=9, color=color)
                )
        
        all_labels.append(label)
        all_data.append({
            "label": label,
            "homo": homo_energy,
            "n_orbitals": len(orbitals)
        })
        x_offset += 1
    
    if x_offset == 0:
        st.warning("No orbital data available for selected molecules")
        return
    
    fig.update_layout(
        title="Orbital Energy Levels",
        xaxis=dict(
            tickvals=list(range(len(all_labels))),
            ticktext=all_labels,
            title="Molecule"
        ),
        yaxis_title="Energy (eV)",
        showlegend=False,
        hovermode="closest"
    )
    
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Solid = Occupied, Dashed = Virtual. Bold = HOMO/LUMO.")
    
    # Data table
    with st.expander("📋 Orbital Summary"):
        st.dataframe(pd.DataFrame(all_data), use_container_width=True, hide_index=True)
