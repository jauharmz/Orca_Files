"""
Energy Visualization Component - Enhanced Implementation

Features:
- Multi-molecule energy comparison with state labels
- Relative energy calculation (kcal/mol)
- HOMO-LUMO energy level diagram
- Energy Pathway with reaction corrections (add/remove species)
- Species energy definitions (OH, H2O, CO2, etc.)
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List, Dict, Optional, Tuple


# =============================================================================
# PRESET SPECIES ENERGIES (Hartrees at B3LYP/6-31G* level, approximate)
# =============================================================================

DEFAULT_SPECIES_ENERGIES = {
    "H2O": -76.4088,      # Water
    "OH": -75.7277,       # Hydroxide
    "H": -0.5000,         # Hydrogen atom
    "H2": -1.1754,        # Hydrogen gas
    "CO2": -188.5908,     # Carbon dioxide
    "HCO3": -264.4392,    # Bicarbonate
    "OMe": -115.0634,     # Methoxide
    "MeOH": -115.7234,    # Methanol
    "NH3": -56.5477,      # Ammonia
    "O2": -150.3176,      # Oxygen
}


# =============================================================================
# MAIN RENDERER
# =============================================================================

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
        mol_options.append({"label": label, "idx": idx, "mol_id": mol_id})
    
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
    energy_tabs = st.tabs(["📊 Energy Comparison", "🔋 HOMO-LUMO Diagram", "🛤️ Energy Pathway"])
    
    with energy_tabs[0]:
        render_energy_comparison(selected_df, relative)
    with energy_tabs[1]:
        render_homo_lumo_diagram(selected_df)
    with energy_tabs[2]:
        render_energy_pathway(df, mol_options, selected)


# =============================================================================
# ENERGY COMPARISON
# =============================================================================

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


# =============================================================================
# HOMO-LUMO DIAGRAM
# =============================================================================

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
            safe_label = str(row['label']).replace('<', '&lt;').replace('>', '&gt;')
            fig.add_trace(go.Scatter(
                x=[x_pos - 0.3, x_pos + 0.3],
                y=[row["homo"], row["homo"]],
                mode="lines",
                line=dict(color=color, width=5),
                name=f"{row['label']} HOMO",
                showlegend=False,
                hovertemplate=f'{safe_label}<br>HOMO: {row["homo"]:.3f} eV<extra></extra>'
            ))
        
        # LUMO level (dashed line)
        if row["lumo"] is not None:
            safe_label = str(row['label']).replace('<', '&lt;').replace('>', '&gt;')
            fig.add_trace(go.Scatter(
                x=[x_pos - 0.3, x_pos + 0.3],
                y=[row["lumo"], row["lumo"]],
                mode="lines",
                line=dict(color=color, width=5, dash="dash"),
                name=f"{row['label']} LUMO",
                showlegend=False,
                hovertemplate=f'{safe_label}<br>LUMO: {row["lumo"]:.3f} eV<extra></extra>'
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


# =============================================================================
# ENERGY PATHWAY
# =============================================================================

def render_energy_pathway(df: pd.DataFrame, mol_options: List[Dict], selected: List[str]):
    """Render energy pathway diagram with reaction corrections."""
    
    st.markdown("### 🛤️ Energy Pathway")
    st.caption("Define reaction pathway and apply species corrections")
    
    # Settings expander
    with st.expander("⚙️ Pathway Settings", expanded=True):
        # Pathway definition
        st.markdown("**Define Pathway Order**")
        pathway_labels = st.multiselect(
            "Select molecules in reaction order (first = reference)",
            selected,
            default=selected[:min(6, len(selected))],
            key="pathway_order"
        )
        
        if len(pathway_labels) < 2:
            st.warning("Select at least 2 molecules to define a pathway")
            return
        
        # Species energies
        st.markdown("**Species Energies (Hartree)**")
        st.caption("Values used for add/remove corrections in reactions")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            e_h2o = st.number_input("H₂O", value=DEFAULT_SPECIES_ENERGIES["H2O"], format="%.4f", key="e_h2o")
        with col2:
            e_oh = st.number_input("OH⁻", value=DEFAULT_SPECIES_ENERGIES["OH"], format="%.4f", key="e_oh")
        with col3:
            e_co2 = st.number_input("CO₂", value=DEFAULT_SPECIES_ENERGIES["CO2"], format="%.4f", key="e_co2")
        with col4:
            e_hco3 = st.number_input("HCO₃⁻", value=DEFAULT_SPECIES_ENERGIES["HCO3"], format="%.4f", key="e_hco3")
        
        species_energies = {
            "H2O": e_h2o,
            "OH": e_oh,
            "CO2": e_co2,
            "HCO3": e_hco3,
        }
        
        # Step corrections
        st.markdown("**Step Corrections**")
        st.caption("Define species added/removed at each step (e.g., '+2 OH, -1 H₂O')")
        
        step_corrections = {}
        for i in range(len(pathway_labels) - 1):
            step_from = pathway_labels[i]
            step_to = pathway_labels[i + 1]
            
            col1, col2 = st.columns([1, 3])
            with col1:
                st.markdown(f"**Step {i+1}:** {step_from.split('[')[0].strip()} → {step_to.split('[')[0].strip()}")
            with col2:
                correction_str = st.text_input(
                    f"Correction (e.g., '+2 OH, -1 H2O')",
                    key=f"step_corr_{i}",
                    label_visibility="collapsed",
                    placeholder="e.g., +2 OH, -1 H2O"
                )
                step_corrections[i] = parse_correction_string(correction_str)
    
    # Get pathway data
    pathway_data = []
    label_to_idx = {m["label"]: m["idx"] for m in mol_options}
    
    for label in pathway_labels:
        if label not in label_to_idx:
            continue
        idx = label_to_idx[label]
        row = df.loc[idx]
        
        gibbs = row.get("gibbs_Eh")
        sp = row.get("single_point_Eh")
        energy = gibbs if gibbs is not None and not pd.isna(gibbs) else sp
        
        if energy is not None and not pd.isna(energy):
            pathway_data.append({
                "label": label,
                "short_label": label.split("[")[0].strip(),
                "energy_Eh": energy,
                "state": row.get("optimized_state", "")
            })
    
    if len(pathway_data) < 2:
        st.warning("Not enough energy data for pathway")
        return
    
    # Calculate corrected energies
    ref_energy = pathway_data[0]["energy_Eh"]
    cumulative_correction = 0.0
    
    for i, pd_item in enumerate(pathway_data):
        # Apply step correction (from previous step)
        if i > 0 and (i - 1) in step_corrections:
            corr = step_corrections[i - 1]
            for species, count in corr.get("add", {}).items():
                if species in species_energies:
                    cumulative_correction += count * species_energies[species]
            for species, count in corr.get("remove", {}).items():
                if species in species_energies:
                    cumulative_correction -= count * species_energies[species]
        
        # Corrected relative energy
        raw_rel = (pd_item["energy_Eh"] - ref_energy) * 627.509
        corrected_rel = raw_rel + (cumulative_correction * 627.509)
        
        pd_item["raw_rel_kcal"] = raw_rel
        pd_item["corrected_rel_kcal"] = corrected_rel
        pd_item["step"] = i
    
    # Create pathway diagram
    fig = go.Figure()
    
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    # Draw pathway line
    x_vals = [p["step"] for p in pathway_data]
    y_vals = [p["corrected_rel_kcal"] for p in pathway_data]
    
    fig.add_trace(go.Scatter(
        x=x_vals,
        y=y_vals,
        mode="lines+markers",
        line=dict(color="#636EFA", width=2),
        marker=dict(size=12, color="#636EFA"),
        name="Corrected",
        hovertemplate='%{text}<br>ΔG: %{y:.2f} kcal/mol<extra></extra>',
        text=[p["label"] for p in pathway_data]
    ))
    
    # Draw horizontal bars at each point
    for i, p in enumerate(pathway_data):
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(
            x=[p["step"] - 0.3, p["step"] + 0.3],
            y=[p["corrected_rel_kcal"], p["corrected_rel_kcal"]],
            mode="lines",
            line=dict(color=color, width=6),
            showlegend=False,
            hoverinfo="skip"
        ))
        
        # Add label
        fig.add_annotation(
            x=p["step"],
            y=p["corrected_rel_kcal"],
            text=f"{p['corrected_rel_kcal']:.1f}",
            showarrow=True,
            arrowhead=0,
            arrowsize=1,
            arrowwidth=1,
            arrowcolor="#666",
            ax=0,
            ay=-30,
            font=dict(size=10),
            bgcolor="white",
            borderpad=2
        )
    
    fig.update_layout(
        title="Energy Pathway (with corrections)",
        xaxis=dict(
            tickvals=x_vals,
            ticktext=[p["short_label"] for p in pathway_data],
            title="Reaction Step"
        ),
        yaxis_title="Relative Energy (kcal/mol)",
        showlegend=True,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    # Add reference line at 0
    fig.add_hline(y=0, line=dict(color="gray", width=1, dash="dash"), opacity=0.5)
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 Pathway Data"):
        pathway_df = pd.DataFrame(pathway_data)
        display_cols = ["step", "label", "energy_Eh", "raw_rel_kcal", "corrected_rel_kcal"]
        available_cols = [c for c in display_cols if c in pathway_df.columns]
        st.dataframe(pathway_df[available_cols], use_container_width=True, hide_index=True)


def parse_correction_string(s: str) -> Dict:
    """Parse correction string like '+2 OH, -1 H2O' into dict.
    
    Returns:
        Dict with 'add' and 'remove' keys containing species counts
    """
    result = {"add": {}, "remove": {}}
    
    if not s or not s.strip():
        return result
    
    # Split by comma
    parts = s.split(",")
    
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        # Parse sign, count, and species
        # Format: +2 OH or -1 H2O or +OH or -H2O
        import re
        match = re.match(r'([+-]?)(\d*)\s*(\w+)', part)
        if match:
            sign = match.group(1) or "+"
            count_str = match.group(2) or "1"
            species = match.group(3)
            
            count = int(count_str)
            
            if sign == "+":
                result["add"][species] = count
            else:
                result["remove"][species] = count
    
    return result
