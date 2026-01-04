"""
Energy Visualization Component

Features:
- Molecule + State comparison
- Energy bar charts
- HOMO-LUMO diagrams
- Relative energy calculation
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from typing import List


def render_energy_viz(df: pd.DataFrame):
    """Render energy visualization with molecule+state filter."""
    
    st.subheader("⚡ Energy Analysis")
    
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
    c1, c2, c3 = st.columns([4, 1, 1])
    with c1:
        selected = st.multiselect("Compare", unique_labels, 
                                 default=unique_labels[:min(6, len(unique_labels))], key="energy_sel")
    with c2:
        if st.button("✅ All", key="energy_all"):
            selected = unique_labels
    with c3:
        relative = st.checkbox("Relative", True, key="energy_rel")
    
    if not selected:
        st.warning("Select molecules to compare")
        return
    
    # Get data
    sel_indices = [l["idx"] for l in labels if l["label"] in selected]
    sel_df = df.loc[sel_indices]
    
    # Tabs
    tabs = st.tabs(["📊 Energy Comparison", "🔋 HOMO-LUMO"])
    
    with tabs[0]:
        render_energy_comparison(sel_df, selected, relative)
    with tabs[1]:
        render_homo_lumo(sel_df, selected)


def render_energy_comparison(df: pd.DataFrame, labels: List[str], relative: bool):
    """Render energy comparison bar chart."""
    
    # Collect energies
    data = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        
        gibbs = row.get("gibbs_Eh")
        sp = row.get("single_point_Eh")
        energy = gibbs if gibbs is not None else sp
        
        if energy is not None:
            data.append({"label": label, "energy": energy, "type": "Gibbs" if gibbs else "SP"})
    
    if not data:
        st.warning("No energy data available")
        return
    
    energy_df = pd.DataFrame(data)
    
    # Calculate relative energies
    if relative:
        min_e = energy_df["energy"].min()
        energy_df["rel_kcal"] = (energy_df["energy"] - min_e) * 627.509  # Eh to kcal/mol
        y_col = "rel_kcal"
        y_title = "Relative Energy (kcal/mol)"
    else:
        energy_df["energy_kcal"] = energy_df["energy"] * 627.509
        y_col = "energy"
        y_title = "Energy (Eh)"
    
    # Plot
    fig = go.Figure()
    
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3']
    
    fig.add_trace(go.Bar(
        x=energy_df["label"],
        y=energy_df[y_col],
        marker_color=[colors[i % len(colors)] for i in range(len(energy_df))],
        text=[f"{v:.2f}" for v in energy_df[y_col]],
        textposition="outside"
    ))
    
    fig.update_layout(
        title="Energy Comparison",
        xaxis_title="Molecule",
        yaxis_title=y_title,
        showlegend=False
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 Data"):
        st.dataframe(energy_df, use_container_width=True, hide_index=True)


def render_homo_lumo(df: pd.DataFrame, labels: List[str]):
    """Render HOMO-LUMO diagram."""
    
    data = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        
        homo = row.get("homo_energy")
        lumo = row.get("lumo_energy")
        
        if homo is not None:
            data.append({"label": label, "homo": homo, "lumo": lumo, 
                        "gap": lumo - homo if lumo else None})
    
    if not data:
        st.warning("No orbital energy data")
        return
    
    hl_df = pd.DataFrame(data)
    
    fig = go.Figure()
    
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A']
    
    for i, row in hl_df.iterrows():
        x_pos = i
        
        # HOMO
        if row["homo"] is not None:
            fig.add_trace(go.Scatter(
                x=[x_pos - 0.3, x_pos + 0.3],
                y=[row["homo"], row["homo"]],
                mode="lines",
                line=dict(color=colors[i % len(colors)], width=4),
                name=f"{row['label']} HOMO",
                showlegend=False
            ))
        
        # LUMO
        if row["lumo"] is not None:
            fig.add_trace(go.Scatter(
                x=[x_pos - 0.3, x_pos + 0.3],
                y=[row["lumo"], row["lumo"]],
                mode="lines",
                line=dict(color=colors[i % len(colors)], width=4, dash="dash"),
                name=f"{row['label']} LUMO",
                showlegend=False
            ))
        
        # Gap arrow
        if row["homo"] is not None and row["lumo"] is not None:
            fig.add_annotation(
                x=x_pos,
                y=(row["homo"] + row["lumo"]) / 2,
                text=f"{row['gap']:.2f} eV",
                showarrow=False,
                font=dict(size=10)
            )
    
    fig.update_layout(
        title="HOMO-LUMO Energy Levels",
        xaxis=dict(
            tickvals=list(range(len(hl_df))),
            ticktext=hl_df["label"].tolist()
        ),
        yaxis_title="Energy (eV)",
        hovermode="x unified"
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    with st.expander("📋 Data"):
        st.dataframe(hl_df, use_container_width=True, hide_index=True)
