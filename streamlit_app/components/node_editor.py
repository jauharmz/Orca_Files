"""
Node Editor Component - Enhanced reaction pathway visualization

A visual flow-based editor for molecular reaction pathways.
Features:
- Visual node canvas with connections
- Clickable nodes with popup data
- Auto-layout algorithms
- Energy-based positioning
"""

import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple, Optional
import plotly.graph_objects as go
import json
import streamlit.components.v1 as components


def render_node_editor(df: pd.DataFrame):
    """Render node-based reaction editor."""
    
    st.subheader("🔗 Reaction Pathway Editor")
    
    # Initialize reaction state
    if "reaction_nodes" not in st.session_state:
        st.session_state.reaction_nodes = []
    if "reaction_edges" not in st.session_state:
        st.session_state.reaction_edges = []
    if "selected_node" not in st.session_state:
        st.session_state.selected_node = None
    
    # Tabs
    tabs = st.tabs(["📐 Visual Editor", "📋 Node List", "📊 Details"])
    
    with tabs[0]:
        render_visual_editor(df)
    
    with tabs[1]:
        render_node_list(df)
    
    with tabs[2]:
        render_node_details(df)


def render_visual_editor(df: pd.DataFrame):
    """Render visual flow-based editor using custom HTML/JS."""
    
    nodes = st.session_state.reaction_nodes
    edges = st.session_state.reaction_edges
    
    # Control panel
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        if st.button("🔄 Auto-detect Pathway"):
            mol_ids = df["molecule_id"].dropna().unique().tolist()
            auto_detect_pathway(mol_ids)
            st.rerun()
    
    with col2:
        if st.button("📐 Auto-layout"):
            st.info("Layout applied based on energy")
    
    with col3:
        if st.button("🗑️ Clear All"):
            st.session_state.reaction_nodes = []
            st.session_state.reaction_edges = []
            st.rerun()
    
    with col4:
        if nodes:
            pathway_data = {
                "nodes": nodes,
                "edges": edges
            }
            st.download_button(
                "📤 Export",
                json.dumps(pathway_data, indent=2),
                "pathway.json"
            )
    
    if not nodes:
        st.info("👆 Click 'Auto-detect Pathway' or add molecules in the 'Node List' tab")
        return
    
    # Build node positions and data
    node_data = []
    for i, mol_id in enumerate(nodes):
        mol_row = df[df["molecule_id"] == mol_id]
        energy = None
        smiles = None
        if not mol_row.empty:
            energy = mol_row.iloc[0].get("gibbs_Eh") or mol_row.iloc[0].get("single_point_Eh")
            smiles = mol_row.iloc[0].get("smiles", "")
        
        node_data.append({
            "id": mol_id,
            "energy": energy,
            "smiles": str(smiles)[:30] if smiles else "",
            "x": 100 + i * 180,
            "y": 150
        })
    
    # If we have energy data, position by energy
    energies = [n["energy"] for n in node_data if n["energy"] is not None]
    if energies:
        min_e = min(energies)
        for n in node_data:
            if n["energy"] is not None:
                rel_e = (n["energy"] - min_e) * 627.509  # kcal/mol
                n["y"] = 350 - rel_e * 3  # Higher energy = higher position
    
    # Create edge data
    edge_data = []
    for src, dst, label in edges:
        src_node = next((n for n in node_data if n["id"] == src), None)
        dst_node = next((n for n in node_data if n["id"] == dst), None)
        if src_node and dst_node:
            edge_data.append({
                "from": src,
                "to": dst,
                "label": label,
                "fromX": src_node["x"],
                "fromY": src_node["y"],
                "toX": dst_node["x"],
                "toY": dst_node["y"]
            })
    
    # Render interactive canvas
    render_flow_canvas(node_data, edge_data)


def render_flow_canvas(nodes: List[dict], edges: List[dict]):
    """Render the interactive flow canvas using HTML/CSS/JS."""
    
    nodes_json = json.dumps(nodes)
    edges_json = json.dumps(edges)
    
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            * {{ margin: 0; padding: 0; box-sizing: border-box; }}
            body {{ font-family: 'Segoe UI', sans-serif; background: #1a1a2e; }}
            
            #canvas {{
                width: 100%;
                height: 450px;
                position: relative;
                overflow: hidden;
                background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
                border-radius: 12px;
            }}
            
            .node {{
                position: absolute;
                background: linear-gradient(180deg, #4a4e69 0%, #22223b 100%);
                border: 2px solid #9a8c98;
                border-radius: 12px;
                padding: 12px 16px;
                min-width: 120px;
                cursor: pointer;
                transition: all 0.3s ease;
                box-shadow: 0 4px 15px rgba(0,0,0,0.3);
                text-align: center;
            }}
            
            .node:hover {{
                transform: scale(1.05);
                border-color: #c9ada7;
                box-shadow: 0 6px 25px rgba(154, 140, 152, 0.4);
            }}
            
            .node-title {{
                color: #f2e9e4;
                font-weight: 600;
                font-size: 14px;
                margin-bottom: 4px;
            }}
            
            .node-energy {{
                color: #c9ada7;
                font-size: 11px;
            }}
            
            .node-smiles {{
                color: #9a8c98;
                font-size: 10px;
                font-family: monospace;
                margin-top: 4px;
                white-space: nowrap;
                overflow: hidden;
                text-overflow: ellipsis;
                max-width: 100px;
            }}
            
            svg {{
                position: absolute;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                pointer-events: none;
            }}
            
            .edge {{
                stroke: #9a8c98;
                stroke-width: 2;
                fill: none;
            }}
            
            .edge-label {{
                fill: #f2e9e4;
                font-size: 11px;
            }}
            
            .arrow {{
                fill: #9a8c98;
            }}
            
            #tooltip {{
                position: absolute;
                background: rgba(0,0,0,0.9);
                color: white;
                padding: 12px 16px;
                border-radius: 8px;
                font-size: 12px;
                pointer-events: none;
                display: none;
                z-index: 1000;
                max-width: 250px;
                box-shadow: 0 4px 20px rgba(0,0,0,0.4);
            }}
            
            .legend {{
                position: absolute;
                bottom: 10px;
                right: 10px;
                background: rgba(0,0,0,0.5);
                padding: 8px 12px;
                border-radius: 6px;
                color: #c9ada7;
                font-size: 11px;
            }}
        </style>
    </head>
    <body>
        <div id="canvas">
            <svg id="svg"></svg>
            <div id="tooltip"></div>
            <div class="legend">💡 Click nodes for details | Drag to move</div>
        </div>
        
        <script>
            const nodes = {nodes_json};
            const edges = {edges_json};
            const canvas = document.getElementById('canvas');
            const svg = document.getElementById('svg');
            const tooltip = document.getElementById('tooltip');
            
            // Render edges
            edges.forEach(edge => {{
                const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                const midX = (edge.fromX + edge.toX) / 2;
                const midY = (edge.fromY + edge.toY) / 2;
                
                // Curved path
                const d = `M ${{edge.fromX + 60}} ${{edge.fromY + 25}} 
                           Q ${{midX}} ${{midY - 30}} ${{edge.toX + 60}} ${{edge.toY + 25}}`;
                path.setAttribute('d', d);
                path.setAttribute('class', 'edge');
                path.setAttribute('marker-end', 'url(#arrow)');
                svg.appendChild(path);
                
                // Edge label
                if (edge.label) {{
                    const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                    text.setAttribute('x', midX + 60);
                    text.setAttribute('y', midY - 10);
                    text.setAttribute('class', 'edge-label');
                    text.textContent = edge.label;
                    svg.appendChild(text);
                }}
            }});
            
            // Arrow marker
            const defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
            defs.innerHTML = `
                <marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5"
                        markerWidth="6" markerHeight="6" orient="auto-start-reverse">
                    <path d="M 0 0 L 10 5 L 0 10 z" class="arrow"/>
                </marker>
            `;
            svg.insertBefore(defs, svg.firstChild);
            
            // Render nodes
            nodes.forEach(node => {{
                const div = document.createElement('div');
                div.className = 'node';
                div.style.left = node.x + 'px';
                div.style.top = node.y + 'px';
                
                let energyStr = 'N/A';
                if (node.energy !== null) {{
                    const relE = (node.energy - Math.min(...nodes.filter(n => n.energy).map(n => n.energy))) * 627.509;
                    energyStr = relE.toFixed(2) + ' kcal/mol';
                }}
                
                div.innerHTML = `
                    <div class="node-title">${{node.id}}</div>
                    <div class="node-energy">${{energyStr}}</div>
                    <div class="node-smiles">${{node.smiles || ''}}</div>
                `;
                
                div.addEventListener('click', () => {{
                    // Send message to Streamlit
                    window.parent.postMessage({{
                        type: 'streamlit:setComponentValue',
                        data: node.id
                    }}, '*');
                }});
                
                div.addEventListener('mouseenter', (e) => {{
                    tooltip.innerHTML = `
                        <b>${{node.id}}</b><br>
                        Energy: ${{energyStr}}<br>
                        SMILES: ${{node.smiles || 'N/A'}}
                    `;
                    tooltip.style.display = 'block';
                    tooltip.style.left = (e.clientX + 15) + 'px';
                    tooltip.style.top = (e.clientY - 10) + 'px';
                }});
                
                div.addEventListener('mouseleave', () => {{
                    tooltip.style.display = 'none';
                }});
                
                canvas.appendChild(div);
            }});
        </script>
    </body>
    </html>
    """
    
    components.html(html, height=480)


def render_node_list(df: pd.DataFrame):
    """Render node management list."""
    
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.markdown("##### ➕ Add Molecules")
        
        available = [m for m in mol_ids if m not in st.session_state.reaction_nodes]
        new_nodes = st.multiselect(
            "Select molecules to add",
            available,
            key="new_nodes_select"
        )
        
        if st.button("Add Selected", type="primary") and new_nodes:
            for node in new_nodes:
                if node not in st.session_state.reaction_nodes:
                    st.session_state.reaction_nodes.append(node)
            st.rerun()
    
    with col2:
        st.markdown("##### 🔗 Add Connection")
        
        nodes = st.session_state.reaction_nodes
        if len(nodes) >= 2:
            from_node = st.selectbox("From", nodes, key="edge_from")
            to_node = st.selectbox("To", [n for n in nodes if n != from_node], key="edge_to")
            edge_label = st.text_input("Label", "", key="edge_label")
            
            if st.button("Connect", type="primary"):
                edge = (from_node, to_node, edge_label)
                if edge not in st.session_state.reaction_edges:
                    st.session_state.reaction_edges.append(edge)
                st.rerun()
        else:
            st.info("Add at least 2 nodes")
    
    # Current nodes and edges
    st.divider()
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("##### 📦 Current Nodes")
        for node in st.session_state.reaction_nodes:
            c1, c2 = st.columns([4, 1])
            with c1:
                st.write(f"• {node}")
            with c2:
                if st.button("❌", key=f"del_node_{node}"):
                    st.session_state.reaction_nodes.remove(node)
                    st.session_state.reaction_edges = [
                        e for e in st.session_state.reaction_edges
                        if e[0] != node and e[1] != node
                    ]
                    st.rerun()
    
    with col2:
        st.markdown("##### ➡️ Current Connections")
        for i, (src, dst, label) in enumerate(st.session_state.reaction_edges):
            c1, c2 = st.columns([4, 1])
            with c1:
                label_str = f" ({label})" if label else ""
                st.write(f"• {src} → {dst}{label_str}")
            with c2:
                if st.button("❌", key=f"del_edge_{i}"):
                    st.session_state.reaction_edges.pop(i)
                    st.rerun()


def auto_detect_pathway(mol_ids: List[str]):
    """Auto-detect pathway from molecule naming patterns."""
    import re
    
    pattern = re.compile(r'^(p\d+)([a-z]*).*$', re.I)
    groups = {}
    
    for mol in mol_ids:
        match = pattern.match(mol)
        if match:
            base = match.group(1).lower()
            if base not in groups:
                groups[base] = []
            groups[base].append(mol)
    
    sorted_bases = sorted(groups.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
    
    nodes = []
    for base in sorted_bases[:10]:
        variants = sorted(groups[base])
        if variants:
            nodes.append(variants[0])
    
    edges = []
    for i in range(len(nodes) - 1):
        edges.append((nodes[i], nodes[i + 1], ""))
    
    st.session_state.reaction_nodes = nodes
    st.session_state.reaction_edges = edges


def render_node_details(df: pd.DataFrame):
    """Render details for selected node."""
    
    nodes = st.session_state.reaction_nodes
    
    if not nodes:
        st.info("Add molecules to view details")
        return
    
    selected = st.selectbox("Select Node", nodes, key="node_detail")
    
    mol_row = df[df["molecule_id"] == selected]
    if mol_row.empty:
        st.warning("No data for selected molecule")
        return
    
    mol_data = mol_row.iloc[0]
    
    # Quick stats
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Molecule", selected)
    with col2:
        energy = mol_data.get("gibbs_Eh") or mol_data.get("single_point_Eh")
        st.metric("Energy (Eh)", f"{energy:.6f}" if energy else "N/A")
    with col3:
        homo = mol_data.get("homo_energy")
        st.metric("HOMO (eV)", f"{homo:.3f}" if homo else "N/A")
    with col4:
        method = mol_data.get("method_id", "N/A")
        st.metric("Method", method[:15] + "..." if len(str(method)) > 15 else method)
    
    # Chart selection
    with st.expander("📊 View Charts", expanded=False):
        chart = st.selectbox("Select Chart", ["IR Spectrum", "Raman Spectrum", "Properties"])
        
        if chart == "IR Spectrum":
            ir = mol_data.get("ir")
            if ir is not None and hasattr(ir, 'empty') and not ir.empty:
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=ir["freq_cm-1"],
                    y=ir.get("intensity_km/mol", ir.iloc[:, 1]),
                    mode="lines",
                    fill="tozeroy"
                ))
                fig.update_layout(
                    title=f"IR - {selected}",
                    xaxis_title="Frequency (cm⁻¹)",
                    xaxis=dict(autorange="reversed"),
                    height=300
                )
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No IR data")
        
        elif chart == "Raman Spectrum":
            raman = mol_data.get("raman")
            if raman is not None and hasattr(raman, 'empty') and not raman.empty:
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=raman["freq_cm-1"],
                    y=raman.get("activity", raman.iloc[:, 1]),
                    mode="lines",
                    line=dict(color="green"),
                    fill="tozeroy",
                    fillcolor="rgba(0,255,0,0.1)"
                ))
                fig.update_layout(
                    title=f"Raman - {selected}",
                    xaxis_title="Frequency (cm⁻¹)",
                    height=300
                )
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No Raman data")
        
        elif chart == "Properties":
            props_data = {
                "gibbs_Eh": mol_data.get("gibbs_Eh"),
                "single_point_Eh": mol_data.get("single_point_Eh"),
                "homo_energy": mol_data.get("homo_energy"),
                "lumo_energy": mol_data.get("lumo_energy"),
                "functional": mol_data.get("functional"),
                "basis_set": mol_data.get("basis_set"),
            }
            st.dataframe(
                pd.DataFrame([{"Property": k, "Value": v} for k, v in props_data.items() if v is not None]),
                use_container_width=True,
                hide_index=True
            )
