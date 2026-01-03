"""
Energy Visualization Component - Bar charts, comparisons, pathways

Features:
- Energy comparison bar chart
- Multi-molecule grouped bars
- Energy pathway diagram
- Unit conversion (Eh, kcal/mol, kJ/mol, eV)
- Relative energy reference
- Data table with export
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List, Optional


# Unit conversion factors from Hartree
UNIT_FACTORS = {
    "Eh": 1.0,
    "kcal/mol": 627.509,
    "kJ/mol": 2625.5,
    "eV": 27.2114
}


def render_energy_viz(df: pd.DataFrame):
    """Render energy visualization tab."""
    
    st.header("⚡ Energy Analysis")
    
    # Molecule selector
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    selected_mols = st.multiselect(
        "Select Molecules",
        mol_ids,
        default=mol_ids[:min(10, len(mol_ids))],
        key="energy_mol_select"
    )
    
    if not selected_mols:
        st.warning("Select at least one molecule")
        return
    
    # Tabs
    tabs = st.tabs(["📊 Comparison", "🛤️ Pathway", "📈 HOMO-LUMO"])
    
    with tabs[0]:
        render_energy_comparison(df, selected_mols)
    
    with tabs[1]:
        render_energy_pathway(df, selected_mols)
    
    with tabs[2]:
        render_homo_lumo(df, selected_mols)


def render_energy_comparison(df: pd.DataFrame, selected_mols: List[str]):
    """Render energy comparison chart."""
    
    st.subheader("📊 Energy Comparison")
    
    # Customization
    with st.expander("⚙️ Customization", expanded=True):
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            energy_type = st.selectbox(
                "Energy Type",
                ["gibbs_Eh", "single_point_Eh"],
                key="energy_type"
            )
        with col2:
            unit = st.selectbox("Unit", list(UNIT_FACTORS.keys()), key="energy_unit")
        with col3:
            relative = st.checkbox("Relative to Minimum", True, key="energy_relative")
        with col4:
            sort_by = st.selectbox("Sort", ["Molecule ID", "Energy"], key="energy_sort")
    
    # Filter data
    subset = df[df["molecule_id"].isin(selected_mols)].copy()
    
    if subset.empty or energy_type not in subset.columns:
        st.warning("No energy data available")
        return
    
    subset = subset[subset[energy_type].notna()]
    
    if subset.empty:
        st.warning(f"No {energy_type} data for selected molecules")
        return
    
    # Convert units
    factor = UNIT_FACTORS[unit]
    energy_col = f"energy_{unit}"
    
    if relative:
        min_energy = subset[energy_type].min()
        subset[energy_col] = (subset[energy_type] - min_energy) * factor
        y_title = f"Relative Energy ({unit})"
    else:
        subset[energy_col] = subset[energy_type] * factor
        y_title = f"Energy ({unit})"
    
    # Sort
    if sort_by == "Energy":
        subset = subset.sort_values(energy_col)
    else:
        subset = subset.sort_values("molecule_id")
    
    # Create plot
    fig = go.Figure()
    
    # Color by state if available
    colors = []
    state_colors = {"S0": "#636EFA", "S1": "#EF553B", "T1": "#00CC96"}
    for _, row in subset.iterrows():
        state = row.get("optimized_state", "S0")
        colors.append(state_colors.get(state, "#636EFA"))
    
    fig.add_trace(go.Bar(
        x=subset["molecule_id"],
        y=subset[energy_col],
        marker_color=colors,
        text=[f"{e:.2f}" for e in subset[energy_col]],
        textposition="outside",
        hovertemplate="<b>%{x}</b><br>Energy: %{y:.4f} " + unit + "<extra></extra>"
    ))
    
    fig.update_layout(
        title="Energy Comparison",
        xaxis_title="Molecule",
        yaxis_title=y_title,
        hovermode="x unified",
        showlegend=False
    )
    
    if relative:
        fig.update_yaxes(range=[0, subset[energy_col].max() * 1.15])
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 View Data Table"):
        table_cols = ["molecule_id", energy_type, energy_col]
        if "optimized_state" in subset.columns:
            table_cols.append("optimized_state")
        if "method_id" in subset.columns:
            table_cols.append("method_id")
        
        from components.data_table import render_simple_table
        render_simple_table(subset[table_cols], "Energy Data")


def render_energy_pathway(df: pd.DataFrame, selected_mols: List[str]):
    """Render energy pathway diagram."""
    
    st.subheader("🛤️ Energy Pathway")
    
    # Customization
    col1, col2 = st.columns(2)
    with col1:
        unit = st.selectbox("Unit", list(UNIT_FACTORS.keys()), index=1, key="pathway_unit")
    with col2:
        show_barriers = st.checkbox("Show Barriers", True, key="pathway_barriers")
    
    # Filter and sort
    subset = df[df["molecule_id"].isin(selected_mols)].copy()
    energy_col = "gibbs_Eh" if "gibbs_Eh" in subset.columns else "single_point_Eh"
    subset = subset[subset[energy_col].notna()]
    
    if len(subset) < 2:
        st.warning("Need at least 2 molecules with energy for pathway")
        return
    
    # Sort by molecule_id to get pathway order
    subset = subset.sort_values("molecule_id")
    
    # Convert to relative energy
    factor = UNIT_FACTORS[unit]
    min_energy = subset[energy_col].min()
    subset["rel_energy"] = (subset[energy_col] - min_energy) * factor
    
    # Create pathway plot
    fig = go.Figure()
    
    x_positions = list(range(len(subset)))
    mol_ids = subset["molecule_id"].tolist()
    rel_energies = subset["rel_energy"].tolist()
    
    # Draw energy levels as horizontal lines
    for i, (mol, energy) in enumerate(zip(mol_ids, rel_energies)):
        fig.add_trace(go.Scatter(
            x=[i - 0.3, i + 0.3],
            y=[energy, energy],
            mode="lines",
            line=dict(color="#636EFA", width=4),
            showlegend=False,
            hovertemplate=f"<b>{mol}</b><br>Energy: {energy:.2f} {unit}<extra></extra>"
        ))
    
    # Draw connecting lines
    for i in range(len(mol_ids) - 1):
        fig.add_trace(go.Scatter(
            x=[i + 0.3, i + 0.7],
            y=[rel_energies[i], rel_energies[i + 1]],
            mode="lines",
            line=dict(color="gray", width=1, dash="dash"),
            showlegend=False
        ))
        
        # Barrier annotation
        if show_barriers:
            delta = rel_energies[i + 1] - rel_energies[i]
            mid_x = i + 0.5
            mid_y = (rel_energies[i] + rel_energies[i + 1]) / 2
            fig.add_annotation(
                x=mid_x,
                y=mid_y,
                text=f"Δ{delta:+.2f}",
                showarrow=False,
                font=dict(size=10, color="red" if delta > 0 else "green")
            )
    
    fig.update_layout(
        title="Energy Pathway",
        xaxis=dict(
            tickmode="array",
            tickvals=x_positions,
            ticktext=mol_ids,
            title="Reaction Coordinate"
        ),
        yaxis_title=f"Relative Energy ({unit})",
        hovermode="x unified"
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_homo_lumo(df: pd.DataFrame, selected_mols: List[str]):
    """Render HOMO-LUMO comparison."""
    
    st.subheader("📈 HOMO-LUMO Comparison")
    
    subset = df[df["molecule_id"].isin(selected_mols)].copy()
    subset = subset[subset["homo_energy"].notna() & subset["lumo_energy"].notna()]
    
    if subset.empty:
        st.warning("No HOMO-LUMO data available")
        return
    
    # Create grouped bar chart
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        name="HOMO",
        x=subset["molecule_id"],
        y=subset["homo_energy"],
        marker_color="blue",
        hovertemplate="<b>%{x}</b><br>HOMO: %{y:.3f} eV<extra></extra>"
    ))
    
    fig.add_trace(go.Bar(
        name="LUMO",
        x=subset["molecule_id"],
        y=subset["lumo_energy"],
        marker_color="red",
        hovertemplate="<b>%{x}</b><br>LUMO: %{y:.3f} eV<extra></extra>"
    ))
    
    # Add gap annotations
    for _, row in subset.iterrows():
        gap = row["lumo_energy"] - row["homo_energy"]
        fig.add_annotation(
            x=row["molecule_id"],
            y=(row["homo_energy"] + row["lumo_energy"]) / 2,
            text=f"Δ{gap:.2f}",
            showarrow=False,
            font=dict(size=10)
        )
    
    fig.update_layout(
        title="HOMO-LUMO Comparison",
        xaxis_title="Molecule",
        yaxis_title="Energy (eV)",
        barmode="group",
        legend=dict(x=0.9, y=0.95)
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 View Data Table"):
        table_cols = ["molecule_id", "homo_energy", "lumo_energy"]
        if "homo_lumo_gap" in subset.columns:
            table_cols.append("homo_lumo_gap")
        from components.data_table import render_simple_table
        render_simple_table(subset[table_cols], "HOMO-LUMO Data")
