"""
Orbital Visualization Component - Enhanced Implementation

Features:
- Multi-molecule orbital comparison with state labels
- Orbital energy level diagram
- HOMO/LUMO highlighting
- Open-shell support (SOMO/SUMO for alpha/beta spins)
- Gap mode selection (HL, SL, SS)
- Connector lines between molecules
- Occupied vs virtual distinction
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List, Optional, Dict


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
    col1, col2, col3 = st.columns([3, 1, 1])
    with col1:
        selected = st.multiselect(
            "Select Molecules",
            unique_labels,
            default=unique_labels[:min(4, len(unique_labels))],
            key="orbital_mol_select"
        )
    with col2:
        n_orbitals = st.slider("Orbitals", 3, 20, 10, key="n_orbitals",
                              help="Number of orbitals around HOMO to display")
    with col3:
        gap_mode = st.selectbox("Gap Mode", ["HL", "SL", "SS"], key="gap_mode",
                               help="HL=HOMO-LUMO, SL=SOMO-LUMO, SS=SOMO-SUMO")
    
    if not selected:
        st.warning("Select molecules to view")
        return
    
    # Get selected data
    selected_indices = [m["idx"] for m in mol_options if m["label"] in selected]
    selected_df = df.loc[selected_indices]
    
    # Settings expander
    with st.expander("⚙️ Orbital Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            show_connectors = st.checkbox("Show Connector Lines", False, key="orbital_connectors",
                                         help="Draw lines connecting corresponding orbitals")
        with col2:
            show_labels = st.checkbox("Show Orbital Labels", True, key="orbital_labels")
    
    render_orbital_levels(selected_df, n_orbitals, gap_mode, show_connectors, show_labels)


def detect_open_shell(orbitals: pd.DataFrame) -> Dict:
    """Detect if system is open-shell and identify SOMO/SUMO.
    
    Returns:
        Dict with 'is_open_shell', 'n_alpha', 'n_beta', 'somo_idx', 'sumo_idx'
    """
    result = {
        "is_open_shell": False,
        "n_alpha": 0,
        "n_beta": 0,
        "somo_idx": None,
        "sumo_idx": None
    }
    
    if orbitals is None or orbitals.empty:
        return result
    
    # Check for spin column
    if "spin" in orbitals.columns:
        # Check for alpha/beta spin labels
        alpha = orbitals[orbitals["spin"].isin(["alpha", "up", "a"])]
        beta = orbitals[orbitals["spin"].isin(["beta", "down", "b"])]
        
        if len(alpha) > 0 and len(beta) > 0:
            # Check occupancy using OCC or occupancy column
            occ_col = "OCC" if "OCC" in orbitals.columns else ("occupancy" if "occupancy" in orbitals.columns else None)
            
            if occ_col:
                alpha_occ = alpha[alpha[occ_col] > 0.5]
                beta_occ = beta[beta[occ_col] > 0.5]
                
                if len(alpha_occ) != len(beta_occ):
                    result["is_open_shell"] = True
                    result["n_alpha"] = len(alpha_occ)
                    result["n_beta"] = len(beta_occ)
                    
                    # SOMO is highest occupied with unpaired electron
                    # SUMO is lowest unoccupied in the minority spin
                    if len(alpha_occ) > len(beta_occ):
                        result["somo_idx"] = len(alpha_occ) - 1
                        beta_virt = beta[beta[occ_col] <= 0.5]
                        if len(beta_virt) > 0:
                            result["sumo_idx"] = len(beta_occ)
    
    return result


def render_orbital_levels(
    df: pd.DataFrame, 
    n_orbitals: int,
    gap_mode: str = "HL",
    show_connectors: bool = False,
    show_labels: bool = True
):
    """Render orbital energy level diagram with enhanced features."""
    
    fig = go.Figure()
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    x_offset = 0
    all_labels = []
    all_data = []
    orbital_positions = []  # For connector lines
    
    for df_idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        
        orbitals = row.get("orbitals")
        homo_energy = row.get("homo_energy")
        lumo_energy = row.get("lumo_energy")
        
        if orbitals is None or (hasattr(orbitals, 'empty') and orbitals.empty):
            continue
        
        # Detect open-shell
        shell_info = detect_open_shell(orbitals)
        
        # Get orbital energies - check multiple possible column names
        if "eV" in orbitals.columns:
            energies = orbitals["eV"].dropna().sort_values()
        elif "energy" in orbitals.columns:
            energies = orbitals["energy"].dropna().sort_values()
        elif "Energy" in orbitals.columns:
            energies = orbitals["Energy"].dropna().sort_values()
        elif "Eh" in orbitals.columns:
            # Convert Eh to eV (1 Eh = 27.2114 eV)
            energies = (orbitals["Eh"] * 27.2114).dropna().sort_values()
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
        
        mol_orbital_positions = []
        
        for i, e in enumerate(energies.values):
            is_occupied = homo_energy is not None and e <= homo_energy
            is_homo = homo_energy is not None and abs(e - homo_energy) < 0.01
            is_lumo = lumo_energy is not None and abs(e - lumo_energy) < 0.01
            
            # Check for SOMO (open-shell)
            is_somo = shell_info["is_open_shell"] and shell_info["somo_idx"] is not None and i == shell_info["somo_idx"]
            is_sumo = shell_info["is_open_shell"] and shell_info["sumo_idx"] is not None and i == shell_info["sumo_idx"]
            
            line_width = 5 if (is_homo or is_lumo or is_somo or is_sumo) else 2
            line_style = "solid" if is_occupied else "dash"
            
            # Special styling for SOMO/SUMO
            if is_somo or is_sumo:
                line_style = "dashdot"
            
            # Build hover text safely
            orbital_type = "SOMO" if is_somo else "SUMO" if is_sumo else "HOMO" if is_homo else "LUMO" if is_lumo else "Occupied" if is_occupied else "Virtual"
            safe_label = str(label).replace('<', '&lt;').replace('>', '&gt;')
            
            fig.add_trace(go.Scatter(
                x=[x_offset - 0.3, x_offset + 0.3],
                y=[e, e],
                mode="lines",
                line=dict(color=color, width=line_width, dash=line_style),
                showlegend=False,
                hovertemplate=f'{safe_label}<br>Energy: {e:.4f} eV<br>{orbital_type}<extra></extra>'
            ))
            
            # Store position for connectors
            mol_orbital_positions.append({
                "x": x_offset,
                "y": e,
                "is_homo": is_homo,
                "is_lumo": is_lumo,
                "is_somo": is_somo,
                "is_sumo": is_sumo,
                "is_occupied": is_occupied
            })
            
            # Add orbital labels
            if show_labels:
                if is_homo:
                    fig.add_annotation(
                        x=x_offset + 0.35, y=e,
                        text="HOMO",
                        showarrow=False,
                        font=dict(size=9, color=color)
                    )
                elif is_lumo:
                    fig.add_annotation(
                        x=x_offset + 0.35, y=e,
                        text="LUMO",
                        showarrow=False,
                        font=dict(size=9, color=color)
                    )
                elif is_somo:
                    fig.add_annotation(
                        x=x_offset + 0.35, y=e,
                        text="SOMO",
                        showarrow=False,
                        font=dict(size=9, color=color)
                    )
                elif is_sumo:
                    fig.add_annotation(
                        x=x_offset + 0.35, y=e,
                        text="SUMO",
                        showarrow=False,
                        font=dict(size=9, color=color)
                    )
        
        orbital_positions.append(mol_orbital_positions)
        
        # Calculate and display gap based on mode
        gap_value = None
        gap_label = ""
        
        if gap_mode == "HL" and homo_energy and lumo_energy:
            gap_value = lumo_energy - homo_energy
            gap_label = "HL Gap"
        elif gap_mode == "SL" and shell_info["is_open_shell"]:
            # SOMO-LUMO gap
            if shell_info["somo_idx"] is not None and lumo_energy:
                somo_e = energies.values[min(shell_info["somo_idx"], len(energies)-1)]
                gap_value = lumo_energy - somo_e
                gap_label = "SL Gap"
        elif gap_mode == "SS" and shell_info["is_open_shell"]:
            # SOMO-SUMO gap
            if shell_info["somo_idx"] is not None and shell_info["sumo_idx"] is not None:
                somo_e = energies.values[min(shell_info["somo_idx"], len(energies)-1)]
                sumo_e = energies.values[min(shell_info["sumo_idx"], len(energies)-1)]
                gap_value = sumo_e - somo_e
                gap_label = "SS Gap"
        
        # Fallback to HL if other modes not available
        if gap_value is None and homo_energy and lumo_energy:
            gap_value = lumo_energy - homo_energy
            gap_label = "HL Gap"
        
        all_labels.append(label)
        all_data.append({
            "label": label,
            "homo": homo_energy,
            "lumo": lumo_energy,
            "gap": gap_value,
            "gap_type": gap_label,
            "n_orbitals": len(orbitals),
            "open_shell": shell_info["is_open_shell"]
        })
        x_offset += 1
    
    # Draw connector lines between molecules
    if show_connectors and len(orbital_positions) > 1:
        for i in range(len(orbital_positions) - 1):
            mol1_pos = orbital_positions[i]
            mol2_pos = orbital_positions[i + 1]
            
            # Connect HOMO to HOMO, LUMO to LUMO
            for p1 in mol1_pos:
                for p2 in mol2_pos:
                    if (p1["is_homo"] and p2["is_homo"]) or (p1["is_lumo"] and p2["is_lumo"]):
                        fig.add_trace(go.Scatter(
                            x=[p1["x"] + 0.3, p2["x"] - 0.3],
                            y=[p1["y"], p2["y"]],
                            mode="lines",
                            line=dict(color="#888", width=1, dash="dot"),
                            showlegend=False,
                            hoverinfo="skip"
                        ))
    
    if x_offset == 0:
        st.warning("No orbital data available for selected molecules")
        return
    
    fig.update_layout(
        title=f"Orbital Energy Levels ({gap_mode} gap mode)",
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
    
    # Legend
    legend_text = "**Legend:** Solid = Occupied, Dashed = Virtual, Dot-Dash = SOMO/SUMO. Bold = HOMO/LUMO."
    st.caption(legend_text)
    
    # Check for open-shell systems
    open_shell_mols = [d["label"] for d in all_data if d.get("open_shell")]
    if open_shell_mols:
        st.info(f"🔄 Open-shell system detected: {', '.join(open_shell_mols)}")
    
    # Data table
    with st.expander("📋 Orbital Summary"):
        summary_df = pd.DataFrame(all_data)
        st.dataframe(summary_df, use_container_width=True, hide_index=True)
