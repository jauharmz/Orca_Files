"""
Node Editor Component - React Flow based reaction pathway visualization

Uses streamlit-flow-component for interactive node-based editing.
Install: pip install streamlit-flow-component

Features:
- Drag-and-drop nodes
- Connect with edges
- Click nodes for details popup
- Energy-based vertical positioning
- ComfyUI/n8n style interface
"""
from __future__ import annotations  # Enable string-based type hints

import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple, TYPE_CHECKING, Any
import json

# Type hints only for static analysis, not runtime
if TYPE_CHECKING:
    from streamlit_flow.state import StreamlitFlowState

# Try to import streamlit-flow-component at runtime
FLOW_AVAILABLE = False
streamlit_flow = None
StreamlitFlowNode = None
StreamlitFlowEdge = None
StreamlitFlowStateClass = None
TreeLayout = None
LayeredLayout = None

try:
    from streamlit_flow import streamlit_flow as sf_func
    from streamlit_flow.elements import StreamlitFlowNode as SFNode, StreamlitFlowEdge as SFEdge
    from streamlit_flow.state import StreamlitFlowState as SFState
    from streamlit_flow.layouts import TreeLayout as TLayout, LayeredLayout as LLayout
    
    streamlit_flow = sf_func
    StreamlitFlowNode = SFNode
    StreamlitFlowEdge = SFEdge
    StreamlitFlowStateClass = SFState
    TreeLayout = TLayout
    LayeredLayout = LLayout
    FLOW_AVAILABLE = True
except ImportError:
    pass


def render_node_editor(df: pd.DataFrame):
    """Render node-based reaction editor."""
    
    st.subheader("🔗 Reaction Pathway Editor")
    
    if not FLOW_AVAILABLE:
        st.error("""
        **streamlit-flow-component not installed!**
        
        Install it with:
        ```
        pip install streamlit-flow-component
        ```
        """)
        render_fallback_editor(df)
        return
    
    # Initialize flow state in session
    if "flow_state" not in st.session_state:
        st.session_state.flow_state = None
    
    # Tabs
    tabs = st.tabs(["📐 Visual Editor", "📋 Setup", "📊 Details"])
    
    with tabs[0]:
        render_react_flow_editor(df)
    
    with tabs[1]:
        render_flow_setup(df)
    
    with tabs[2]:
        render_node_details(df)


def render_react_flow_editor(df: pd.DataFrame):
    """Render React Flow based editor."""
    
    # Control buttons
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("🔄 Auto-detect Pathway"):
            mol_ids = df["molecule_id"].dropna().unique().tolist()
            state = create_auto_pathway(df, mol_ids)
            st.session_state.flow_state = state
            st.rerun()
    
    with col2:
        layout_type = st.selectbox("Layout", ["Layered", "Tree"], key="flow_layout")
    
    with col3:
        if st.button("🗑️ Clear All"):
            st.session_state.flow_state = None
            st.rerun()
    
    # Get or create flow state
    if st.session_state.flow_state is None:
        # Create empty state
        st.session_state.flow_state = StreamlitFlowStateClass(nodes=[], edges=[])
    
    # Layout configuration
    if layout_type == "Tree":
        layout = TreeLayout(direction="right")
    else:
        layout = LayeredLayout(direction="right")
    
    # Render the flow component
    state = streamlit_flow(
        key="reaction_flow",
        state=st.session_state.flow_state,
        layout=layout,
        fit_view=True,
        height=500,
        enable_node_menu=True,
        enable_edge_menu=True,
        enable_pane_menu=True,
        allow_new_edges=True,
        animate_new_edges=True,
        style={"backgroundColor": "#1a1a2e"}
    )
    
    # Update state
    st.session_state.flow_state = state
    
    # Show selected node info
    if state.selected_id:
        st.info(f"Selected: **{state.selected_id}**")
        show_node_popup(df, state.selected_id)


def create_auto_pathway(df: pd.DataFrame, mol_ids: List[str]):
    """Create automatic pathway from molecule naming."""
    import re
    
    # Group by base molecule name
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
    sorted_bases = sorted(groups.keys(), key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0)
    
    # Take first variant from each base (up to 10)
    selected_mols = []
    for base in sorted_bases[:10]:
        variants = sorted(groups[base])
        if variants:
            selected_mols.append(variants[0])
    
    # Create nodes with energy data
    nodes = []
    energies = []
    
    for i, mol_id in enumerate(selected_mols):
        mol_row = df[df["molecule_id"] == mol_id]
        energy = None
        state = ""
        if not mol_row.empty:
            energy = mol_row.iloc[0].get("gibbs_Eh") or mol_row.iloc[0].get("single_point_Eh")
            state = mol_row.iloc[0].get("optimized_state", "")
        
        energies.append(energy)
        
        # Calculate y position based on energy
        y_pos = 100
        if energy is not None and energies[0] is not None:
            rel_e = (energy - min(e for e in energies if e)) * 627.509  # kcal/mol
            y_pos = 100 + rel_e * 2
        
        energy_str = f"{energy:.4f} Eh" if energy else "N/A"
        
        nodes.append(StreamlitFlowNode(
            id=mol_id,
            pos=(i * 200, y_pos),
            data={
                "content": f"**{mol_id}**\n\n{state}\n\n{energy_str}"
            },
            node_type="default",
            source_position="right",
            target_position="left",
            draggable=True
        ))
    
    # Create edges
    edges = []
    for i in range(len(selected_mols) - 1):
        edges.append(StreamlitFlowEdge(
            id=f"e{i}",
            source=selected_mols[i],
            target=selected_mols[i + 1],
            animated=True,
            edge_type="smoothstep"
        ))
    
    return StreamlitFlowStateClass(nodes=nodes, edges=edges)


def render_flow_setup(df: pd.DataFrame):
    """Setup nodes and edges manually."""
    
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("##### ➕ Add Nodes")
        selected_mols = st.multiselect("Select molecules", mol_ids, key="flow_add_mols")
        
        if st.button("Add Selected", type="primary") and selected_mols:
            if st.session_state.flow_state is None:
                st.session_state.flow_state = StreamlitFlowStateClass(nodes=[], edges=[])
            
            existing_ids = [n.id for n in st.session_state.flow_state.nodes]
            
            for i, mol_id in enumerate(selected_mols):
                if mol_id not in existing_ids:
                    mol_row = df[df["molecule_id"] == mol_id]
                    energy = None
                    if not mol_row.empty:
                        energy = mol_row.iloc[0].get("gibbs_Eh")
                    
                    st.session_state.flow_state.nodes.append(
                        StreamlitFlowNode(
                            id=mol_id,
                            pos=(len(existing_ids) * 200, 100),
                            data={"content": f"**{mol_id}**\n\n{energy:.4f} Eh" if energy else f"**{mol_id}**"},
                            node_type="default",
                            draggable=True
                        )
                    )
            st.rerun()
    
    with col2:
        st.markdown("##### Current Nodes")
        if st.session_state.flow_state:
            for node in st.session_state.flow_state.nodes:
                st.write(f"• {node.id}")
        else:
            st.info("No nodes yet")


def show_node_popup(df: pd.DataFrame, node_id: str):
    """Show popup with node details."""
    
    mol_row = df[df["molecule_id"] == node_id]
    if mol_row.empty:
        return
    
    data = mol_row.iloc[0]
    
    with st.expander(f"📊 Details: {node_id}", expanded=True):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            energy = data.get("gibbs_Eh") or data.get("single_point_Eh")
            st.metric("Energy (Eh)", f"{energy:.6f}" if energy else "N/A")
        
        with col2:
            homo = data.get("homo_energy")
            st.metric("HOMO (eV)", f"{homo:.3f}" if homo else "N/A")
        
        with col3:
            state = data.get("optimized_state", "N/A")
            st.metric("State", state)


def render_node_details(df: pd.DataFrame):
    """Render detailed view for nodes."""
    
    if st.session_state.flow_state and st.session_state.flow_state.nodes:
        node_ids = [n.id for n in st.session_state.flow_state.nodes]
        selected = st.selectbox("Select Node", node_ids, key="detail_node")
        
        mol_row = df[df["molecule_id"] == selected]
        if not mol_row.empty:
            data = mol_row.iloc[0]
            
            # Show properties table
            props = {
                "molecule_id": data.get("molecule_id"),
                "optimized_state": data.get("optimized_state"),
                "gibbs_Eh": data.get("gibbs_Eh"),
                "single_point_Eh": data.get("single_point_Eh"),
                "homo_energy": data.get("homo_energy"),
                "lumo_energy": data.get("lumo_energy"),
            }
            st.dataframe(
                pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items() if v is not None]),
                use_container_width=True,
                hide_index=True
            )
    else:
        st.info("Add nodes to view details")


def render_fallback_editor(df: pd.DataFrame):
    """Fallback editor when streamlit-flow-component is not available."""
    
    st.warning("Using fallback editor. Install `streamlit-flow-component` for the full experience.")
    
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    # Simple node management
    if "simple_nodes" not in st.session_state:
        st.session_state.simple_nodes = []
    if "simple_edges" not in st.session_state:
        st.session_state.simple_edges = []
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("##### Add Nodes")
        new_nodes = st.multiselect("Select", [m for m in mol_ids if m not in st.session_state.simple_nodes])
        if st.button("Add") and new_nodes:
            st.session_state.simple_nodes.extend(new_nodes)
            st.rerun()
    
    with col2:
        st.markdown("##### Add Edge")
        if len(st.session_state.simple_nodes) >= 2:
            src = st.selectbox("From", st.session_state.simple_nodes, key="fb_src")
            dst = st.selectbox("To", [n for n in st.session_state.simple_nodes if n != src], key="fb_dst")
            if st.button("Connect"):
                st.session_state.simple_edges.append((src, dst))
                st.rerun()
    
    # Display
    st.markdown("##### Current Pathway")
    st.write("Nodes:", st.session_state.simple_nodes)
    st.write("Edges:", st.session_state.simple_edges)
