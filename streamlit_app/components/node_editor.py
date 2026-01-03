"""
Node Editor Component - Reaction pathway visualization

A drag-and-drop node editor for visualizing molecular reaction pathways.
Similar to ComfyUI/n8n style node-based interfaces.

Features:
- Molecules as clickable nodes
- Connect nodes with arrows
- Click node for details/charts
- Energy labels on connections
"""

import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple, Optional
import plotly.graph_objects as go
import json


def render_node_editor(df: pd.DataFrame):
    """Render node-based reaction editor."""
    
    st.header("🔗 Reaction Pathway Editor")
    
    st.info("""
    **How to use:**
    1. Select molecules to add as nodes
    2. Click "Add Connection" to link molecules
    3. Click on nodes to view molecule details
    4. Use the pathway visualization below
    """)
    
    # Initialize reaction state
    if "reaction_nodes" not in st.session_state:
        st.session_state.reaction_nodes = []
    if "reaction_edges" not in st.session_state:
        st.session_state.reaction_edges = []
    
    # Tabs
    tabs = st.tabs(["🎨 Editor", "📊 Visualization", "📋 Details"])
    
    with tabs[0]:
        render_editor_panel(df)
    
    with tabs[1]:
        render_pathway_viz(df)
    
    with tabs[2]:
        render_node_details(df)


def render_editor_panel(df: pd.DataFrame):
    """Render node editor panel."""
    
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📦 Add Nodes")
        
        new_nodes = st.multiselect(
            "Select molecules to add",
            [m for m in mol_ids if m not in st.session_state.reaction_nodes],
            key="new_nodes"
        )
        
        if st.button("➕ Add Selected Nodes", type="primary"):
            for node in new_nodes:
                if node not in st.session_state.reaction_nodes:
                    st.session_state.reaction_nodes.append(node)
            st.rerun()
        
        # Current nodes
        st.markdown("**Current Nodes:**")
        for node in st.session_state.reaction_nodes:
            col_a, col_b = st.columns([3, 1])
            with col_a:
                st.write(f"• {node}")
            with col_b:
                if st.button("❌", key=f"remove_{node}"):
                    st.session_state.reaction_nodes.remove(node)
                    # Remove related edges
                    st.session_state.reaction_edges = [
                        e for e in st.session_state.reaction_edges
                        if e[0] != node and e[1] != node
                    ]
                    st.rerun()
    
    with col2:
        st.subheader("🔗 Add Connections")
        
        nodes = st.session_state.reaction_nodes
        if len(nodes) >= 2:
            from_node = st.selectbox("From", nodes, key="edge_from")
            to_node = st.selectbox("To", [n for n in nodes if n != from_node], key="edge_to")
            edge_label = st.text_input("Label (optional)", "", key="edge_label")
            
            if st.button("➕ Add Connection", type="primary"):
                edge = (from_node, to_node, edge_label)
                if edge not in st.session_state.reaction_edges:
                    st.session_state.reaction_edges.append(edge)
                st.rerun()
        else:
            st.info("Add at least 2 nodes to create connections")
        
        # Current edges
        st.markdown("**Current Connections:**")
        for i, (src, dst, label) in enumerate(st.session_state.reaction_edges):
            col_a, col_b = st.columns([3, 1])
            with col_a:
                label_str = f" ({label})" if label else ""
                st.write(f"• {src} → {dst}{label_str}")
            with col_b:
                if st.button("❌", key=f"remove_edge_{i}"):
                    st.session_state.reaction_edges.pop(i)
                    st.rerun()
    
    # Quick actions
    st.divider()
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🔄 Clear All"):
            st.session_state.reaction_nodes = []
            st.session_state.reaction_edges = []
            st.rerun()
    
    with col2:
        if st.button("📥 Auto-detect Pathway"):
            # Auto-detect pathway from naming (p1x → p2x → p3x)
            auto_detect_pathway(mol_ids)
            st.rerun()
    
    with col3:
        # Export pathway
        if st.session_state.reaction_nodes:
            pathway_data = {
                "nodes": st.session_state.reaction_nodes,
                "edges": st.session_state.reaction_edges
            }
            st.download_button(
                "📤 Export",
                json.dumps(pathway_data, indent=2),
                "pathway.json",
                "application/json"
            )


def auto_detect_pathway(mol_ids: List[str]):
    """Auto-detect pathway from molecule naming patterns."""
    import re
    
    # Group by base name
    pattern = re.compile(r'^(p\d+)([a-z]*).*$', re.I)
    groups = {}
    
    for mol in mol_ids:
        match = pattern.match(mol)
        if match:
            base = match.group(1).lower()
            if base not in groups:
                groups[base] = []
            groups[base].append(mol)
    
    # Sort by base name
    sorted_bases = sorted(groups.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
    
    # Take first variant from each base
    nodes = []
    for base in sorted_bases[:10]:  # Limit to 10
        variants = sorted(groups[base])
        if variants:
            nodes.append(variants[0])
    
    # Create edges
    edges = []
    for i in range(len(nodes) - 1):
        edges.append((nodes[i], nodes[i + 1], ""))
    
    st.session_state.reaction_nodes = nodes
    st.session_state.reaction_edges = edges


def render_pathway_viz(df: pd.DataFrame):
    """Render pathway visualization."""
    
    nodes = st.session_state.reaction_nodes
    edges = st.session_state.reaction_edges
    
    if not nodes:
        st.info("Add molecules to visualize pathway")
        return
    
    st.subheader("📊 Pathway Visualization")
    
    # Options
    col1, col2 = st.columns(2)
    with col1:
        show_energy = st.checkbox("Show Energy", True, key="pathway_energy")
    with col2:
        layout = st.selectbox("Layout", ["Horizontal", "Vertical"], key="pathway_layout")
    
    # Get energy data
    energy_data = {}
    for node in nodes:
        mol_row = df[df["molecule_id"] == node]
        if not mol_row.empty:
            energy = mol_row.iloc[0].get("gibbs_Eh") or mol_row.iloc[0].get("single_point_Eh")
            energy_data[node] = energy
    
    # Calculate positions
    n = len(nodes)
    if layout == "Horizontal":
        positions = {node: (i * 2, 0) for i, node in enumerate(nodes)}
    else:
        positions = {node: (0, -i * 2) for i, node in enumerate(nodes)}
    
    # If we have energy data and showing energy, use energy for Y
    if show_energy and energy_data:
        min_energy = min(e for e in energy_data.values() if e is not None)
        for node in nodes:
            if energy_data.get(node) is not None:
                rel_energy = (energy_data[node] - min_energy) * 627.509  # kcal/mol
                x, y = positions[node]
                if layout == "Horizontal":
                    positions[node] = (x, rel_energy)
                else:
                    positions[node] = (rel_energy, y)
    
    # Create plot
    fig = go.Figure()
    
    # Draw edges first
    for src, dst, label in edges:
        if src in positions and dst in positions:
            x0, y0 = positions[src]
            x1, y1 = positions[dst]
            
            # Arrow line
            fig.add_trace(go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(color="gray", width=2),
                showlegend=False,
                hoverinfo="skip"
            ))
            
            # Arrow head
            fig.add_annotation(
                x=x1, y=y1,
                ax=x0, ay=y0,
                xref="x", yref="y",
                axref="x", ayref="y",
                showarrow=True,
                arrowhead=2,
                arrowsize=1.5,
                arrowwidth=2,
                arrowcolor="gray"
            )
            
            # Label
            if label:
                fig.add_annotation(
                    x=(x0 + x1) / 2,
                    y=(y0 + y1) / 2,
                    text=label,
                    showarrow=False,
                    font=dict(size=10)
                )
    
    # Draw nodes
    for node in nodes:
        if node in positions:
            x, y = positions[node]
            energy = energy_data.get(node)
            
            hover_text = f"<b>{node}</b>"
            if energy:
                rel_e = (energy - min_energy) * 627.509 if min_energy else 0
                hover_text += f"<br>Energy: {rel_e:.2f} kcal/mol"
            
            fig.add_trace(go.Scatter(
                x=[x],
                y=[y],
                mode="markers+text",
                marker=dict(size=40, color="#636EFA", symbol="square"),
                text=[node],
                textposition="middle center",
                textfont=dict(color="white", size=10),
                name=node,
                hovertemplate=hover_text + "<extra></extra>",
                showlegend=False
            ))
    
    fig.update_layout(
        title="Reaction Pathway",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False) if layout == "Horizontal" else dict(title="Relative Energy (kcal/mol)"),
        yaxis=dict(title="Relative Energy (kcal/mol)") if layout == "Horizontal" and show_energy else dict(showgrid=False, zeroline=False, showticklabels=False),
        height=500,
        hovermode="closest"
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_node_details(df: pd.DataFrame):
    """Render details for selected node."""
    
    nodes = st.session_state.reaction_nodes
    
    if not nodes:
        st.info("Add molecules to view details")
        return
    
    selected = st.selectbox("Select Node to View", nodes, key="node_detail_select")
    
    mol_row = df[df["molecule_id"] == selected]
    if mol_row.empty:
        st.warning("No data for selected molecule")
        return
    
    mol_data = mol_row.iloc[0]
    
    st.subheader(f"📋 Details: {selected}")
    
    # Properties
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Method", mol_data.get("method_id", "N/A"))
    with col2:
        gibbs = mol_data.get("gibbs_Eh")
        st.metric("Gibbs (Eh)", f"{gibbs:.6f}" if gibbs else "N/A")
    with col3:
        homo = mol_data.get("homo_energy")
        st.metric("HOMO (eV)", f"{homo:.3f}" if homo else "N/A")
    
    # Charts
    chart_type = st.selectbox(
        "View Chart",
        ["None", "IR Spectrum", "Raman Spectrum", "Orbital Levels"],
        key="node_chart_select"
    )
    
    if chart_type == "IR Spectrum":
        ir = mol_data.get("ir")
        if ir is not None and not ir.empty:
            from components.spectra_viz import render_ir_spectra
            # Simple IR plot
            import plotly.graph_objects as go
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=ir["freq_cm-1"], 
                y=ir["intensity_km/mol"] if "intensity_km/mol" in ir.columns else ir.iloc[:, 1],
                mode="lines",
                fill="tozeroy"
            ))
            fig.update_layout(
                title=f"IR Spectrum - {selected}",
                xaxis_title="Frequency (cm⁻¹)",
                xaxis=dict(autorange="reversed")
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No IR data available")
    
    elif chart_type == "Raman Spectrum":
        raman = mol_data.get("raman")
        if raman is not None and not raman.empty:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=raman["freq_cm-1"],
                y=raman["activity"] if "activity" in raman.columns else raman.iloc[:, 1],
                mode="lines",
                line=dict(color="green"),
                fill="tozeroy",
                fillcolor="rgba(0,255,0,0.1)"
            ))
            fig.update_layout(
                title=f"Raman Spectrum - {selected}",
                xaxis_title="Frequency (cm⁻¹)"
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No Raman data available")
    
    elif chart_type == "Orbital Levels":
        orbitals = mol_data.get("orbitals")
        if orbitals is not None and not orbitals.empty:
            homo = mol_data.get("homo_energy", 0)
            lumo = mol_data.get("lumo_energy", 0)
            st.metric("HOMO-LUMO Gap", f"{lumo - homo:.3f} eV" if homo and lumo else "N/A")
            
            from components.data_table import render_simple_table
            render_simple_table(orbitals.head(20), "Orbital Data")
        else:
            st.info("No orbital data available")
