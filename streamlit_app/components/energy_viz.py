"""
Energy Visualization Component - Complete Implementation

Features:
- Multi-molecule energy comparison with state labels
- Relative energy calculation (kcal/mol)
- HOMO-LUMO energy level diagram
- Bar charts and data tables
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List


def render_energy_viz(df: pd.DataFrame):
    """Render energy visualization with molecule+state filter."""
    
    st.subheader("⚡ Energy Analysis")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build molecule options with state labels
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
    
    if not unique_labels:
        st.warning("No molecules available")
        return
    
    # Selection - default top 10
    col1, col2 = st.columns([4, 1])
    with col1:
        selected = st.multiselect(
            "Select Molecules to Compare",
            unique_labels,
            default=unique_labels[:min(10, len(unique_labels))],
            key="energy_mol_select"
        )
    with col2:
        relative = st.checkbox("Relative", True, key="energy_relative")
    
    if not selected:
        st.warning("Select molecules to compare")
        return
    
    # Get selected data
    selected_indices = [m["idx"] for m in mol_options if m["label"] in selected]
    selected_df = df.loc[selected_indices]
    
    # Tabs
    energy_tabs = st.tabs(["📊 Energy Comparison", "🔋 HOMO-LUMO Diagram"])
    
    with energy_tabs[0]:
        render_energy_comparison(selected_df, relative)
    with energy_tabs[1]:
        render_homo_lumo_diagram(selected_df)


def render_energy_comparison(df: pd.DataFrame, relative: bool):
    """Render energy comparison bar chart."""
    
    # Collect energy data
    data = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        
        gibbs = row.get("gibbs_Eh")
        sp = row.get("single_point_Eh")
        energy = gibbs if gibbs is not None and not pd.isna(gibbs) else sp
        
        if energy is not None and not pd.isna(energy):
            data.append({
                "label": label,
                "energy_Eh": energy,
                "type": "Gibbs" if gibbs is not None and not pd.isna(gibbs) else "SP"
            })
    
    if not data:
        st.warning("No energy data available for selected molecules")
        return
    
    energy_df = pd.DataFrame(data)
    
    # Calculate relative energies
    if relative:
        min_e = energy_df["energy_Eh"].min()
        energy_df["rel_kcal"] = (energy_df["energy_Eh"] - min_e) * 627.509  # Hartree to kcal/mol
        y_col = "rel_kcal"
        y_title = "Relative Energy (kcal/mol)"
        
        # Find reference
        ref_mol = energy_df.loc[energy_df["energy_Eh"] == min_e, "label"].values[0]
        st.info(f"📌 Reference: {ref_mol} (lowest energy)")
    else:
        y_col = "energy_Eh"
        y_title = "Energy (Hartree)"
    
    # Create bar chart
    fig = go.Figure()
    
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    fig.add_trace(go.Bar(
        x=energy_df["label"],
        y=energy_df[y_col],
        marker_color=[colors[i % len(colors)] for i in range(len(energy_df))],
        text=[f"{v:.2f}" for v in energy_df[y_col]],
        textposition="outside",
        hovertemplate='%{x}<br>%{y:.4f}<extra></extra>'
    ))
    
    fig.update_layout(
        title="Energy Comparison",
        xaxis_title="Molecule",
        yaxis_title=y_title,
        showlegend=False,
        xaxis_tickangle=-45
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 Energy Data"):
        display_df = energy_df.copy()
        if relative:
            display_df["energy_kcal"] = display_df["energy_Eh"] * 627.509
        st.dataframe(display_df, use_container_width=True, hide_index=True)


def render_homo_lumo_diagram(df: pd.DataFrame):
    """Render HOMO-LUMO energy level diagram."""
    
    # Collect orbital data
    data = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        
        homo = row.get("homo_energy")
        lumo = row.get("lumo_energy")
        
        if homo is not None and not pd.isna(homo):
            data.append({
                "label": label,
                "homo": homo,
                "lumo": lumo if lumo is not None and not pd.isna(lumo) else None,
                "gap": lumo - homo if lumo is not None and not pd.isna(lumo) else None
            })
    
    if not data:
        st.warning("No HOMO/LUMO data available")
        return
    
    hl_df = pd.DataFrame(data)
    
    # Create diagram
    fig = go.Figure()
    
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    for i, row in hl_df.iterrows():
        x_pos = i
        color = colors[i % len(colors)]
        
        # HOMO level (solid line)
        if row["homo"] is not None:
            fig.add_trace(go.Scatter(
                x=[x_pos - 0.3, x_pos + 0.3],
                y=[row["homo"], row["homo"]],
                mode="lines",
                line=dict(color=color, width=5),
                name=f"{row['label']} HOMO",
                showlegend=False,
                hovertemplate=f'{row["label"]}<br>HOMO: {row["homo"]:.3f} eV<extra></extra>'
            ))
        
        # LUMO level (dashed line)
        if row["lumo"] is not None:
            fig.add_trace(go.Scatter(
                x=[x_pos - 0.3, x_pos + 0.3],
                y=[row["lumo"], row["lumo"]],
                mode="lines",
                line=dict(color=color, width=5, dash="dash"),
                name=f"{row['label']} LUMO",
                showlegend=False,
                hovertemplate=f'{row["label"]}<br>LUMO: {row["lumo"]:.3f} eV<extra></extra>'
            ))
        
        # Gap annotation
        if row["homo"] is not None and row["lumo"] is not None:
            mid_y = (row["homo"] + row["lumo"]) / 2
            fig.add_annotation(
                x=x_pos,
                y=mid_y,
                text=f"{row['gap']:.2f} eV",
                showarrow=False,
                font=dict(size=10, color="#333"),
                bgcolor="white",
                borderpad=2
            )
            
            # Connecting line
            fig.add_trace(go.Scatter(
                x=[x_pos, x_pos],
                y=[row["homo"], row["lumo"]],
                mode="lines",
                line=dict(color=color, width=1, dash="dot"),
                showlegend=False,
                hoverinfo="skip"
            ))
    
    fig.update_layout(
        title="HOMO-LUMO Energy Levels",
        xaxis=dict(
            tickvals=list(range(len(hl_df))),
            ticktext=hl_df["label"].tolist(),
            title="Molecule"
        ),
        yaxis_title="Energy (eV)",
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Solid = HOMO (occupied), Dashed = LUMO (virtual)")
    
    # Data table
    with st.expander("📋 HOMO-LUMO Data"):
        st.dataframe(hl_df, use_container_width=True, hide_index=True)
