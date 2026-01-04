"""
Orbital Visualization Component

Features:
- Molecule + State filter
- Orbital energy levels
- HOMO/LUMO highlighting
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from typing import List


def render_orbital_viz(df: pd.DataFrame):
    """Render orbital visualization with molecule+state filter."""
    
    st.subheader("🔮 Orbital Analysis")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build labels
    labels = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        labels.append({"label": label, "idx": idx})
    
    unique_labels = list(dict.fromkeys([l["label"] for l in labels]))
    
    # Filter
    c1, c2 = st.columns([4, 2])
    with c1:
        selected = st.multiselect("Compare Orbitals", unique_labels, 
                                 default=unique_labels[:min(4, len(unique_labels))], key="orbital_sel")
    with c2:
        n_orbitals = st.slider("Orbitals Around HOMO", 3, 20, 8, key="n_orb")
    
    if not selected:
        st.warning("Select molecules to compare")
        return
    
    # Get data
    sel_indices = [l["idx"] for l in labels if l["label"] in selected]
    sel_df = df.loc[sel_indices]
    
    render_orbital_levels(sel_df, selected, n_orbitals)


def render_orbital_levels(df: pd.DataFrame, labels: List[str], n_orbitals: int):
    """Render orbital energy level diagram."""
    
    fig = go.Figure()
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A']
    
    x_offset = 0
    all_data = []
    
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        
        orbitals = row.get("orbitals")
        if orbitals is None or (hasattr(orbitals, 'empty') and orbitals.empty):
            continue
        
        homo_energy = row.get("homo_energy")
        
        # Get orbital energies
        if "energy" in orbitals.columns:
            energies = orbitals["energy"].dropna().sort_values()
        else:
            continue
        
        # Filter to N orbitals around HOMO
        if homo_energy is not None:
            diff = (energies - homo_energy).abs()
            closest = diff.nsmallest(n_orbitals * 2).index
            energies = energies.loc[closest].sort_values()
        else:
            energies = energies.tail(n_orbitals * 2)
        
        color_idx = list(df.index).index(idx) % len(colors)
        
        for i, e in enumerate(energies):
            is_occupied = homo_energy is not None and e <= homo_energy
            is_homo = homo_energy is not None and abs(e - homo_energy) < 0.01
            is_lumo = homo_energy is not None and e > homo_energy and i == len([x for x in energies if x <= homo_energy])
            
            line_width = 4 if (is_homo or is_lumo) else 2
            line_style = "solid" if is_occupied else "dash"
            
            fig.add_trace(go.Scatter(
                x=[x_offset - 0.3, x_offset + 0.3],
                y=[e, e],
                mode="lines",
                line=dict(color=colors[color_idx], width=line_width, dash=line_style),
                showlegend=False,
                hovertemplate=f"{label}<br>Energy: {e:.3f} eV<extra></extra>"
            ))
            
            if is_homo:
                fig.add_annotation(x=x_offset + 0.4, y=e, text="HOMO", showarrow=False, font=dict(size=9))
            if is_lumo:
                fig.add_annotation(x=x_offset + 0.4, y=e, text="LUMO", showarrow=False, font=dict(size=9))
        
        all_data.append({"label": label, "homo": homo_energy, "n_orbitals": len(orbitals)})
        x_offset += 1
    
    if x_offset == 0:
        st.warning("No orbital data available")
        return
    
    # Add molecule labels
    fig.update_layout(
        title="Orbital Energy Levels",
        xaxis=dict(
            tickvals=list(range(x_offset)),
            ticktext=[d["label"] for d in all_data],
            title="Molecule"
        ),
        yaxis_title="Energy (eV)",
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Solid = occupied, Dashed = virtual. Bold = HOMO/LUMO.")
    
    with st.expander("📋 Data"):
        st.dataframe(pd.DataFrame(all_data), use_container_width=True, hide_index=True)
