"""
Molecule Viewer Component - Complete 2D/3D visualization

Features:
- Molecule + State selector (or combined selector)
- 3D viewer with py3Dmol (ball-stick, stick, sphere styles)
- 2D viewer with RDKit SVG
- Property tables with coordinates and Mulliken charges
- Hover tooltips with atom info
"""

import streamlit as st
import pandas as pd
from typing import Optional
import streamlit.components.v1 as components


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer with molecule+state selection."""
    
    st.subheader("🧬 Molecule Viewer")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build options with molecule_id + state combined labels
    options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        options.append({"label": label, "idx": idx, "mol_id": mol_id, "state": state})
    
    # Remove duplicates by label
    seen = set()
    unique_options = []
    for opt in options:
        if opt["label"] not in seen:
            seen.add(opt["label"])
            unique_options.append(opt)
    
    if not unique_options:
        st.warning("No molecules available")
        return
    
    labels = [opt["label"] for opt in unique_options]
    
    # Selection
    selected_label = st.selectbox("Select Molecule", labels, key="mol_viewer_select")
    
    # Get selected row
    selected_opt = next((opt for opt in unique_options if opt["label"] == selected_label), None)
    if selected_opt is None:
        return
    
    mol_row = df.loc[selected_opt["idx"]]
    
    # Info metrics row
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        st.metric("Molecule", mol_row.get("molecule_id", "N/A"))
    with c2:
        state = mol_row.get("optimized_state", "N/A")
        st.metric("State", state if state and str(state) != "nan" else "N/A")
    with c3:
        energy = mol_row.get("gibbs_Eh") or mol_row.get("single_point_Eh")
        st.metric("Energy (Eh)", f"{energy:.4f}" if energy else "N/A")
    with c4:
        homo = mol_row.get("homo_energy")
        st.metric("HOMO (eV)", f"{homo:.2f}" if homo else "N/A")
    with c5:
        lumo = mol_row.get("lumo_energy")
        st.metric("LUMO (eV)", f"{lumo:.2f}" if lumo else "N/A")
    
    # View tabs
    view_tabs = st.tabs(["🔮 3D View", "📐 2D View", "📊 Properties"])
    
    with view_tabs[0]:
        render_3d_view(mol_row)
    with view_tabs[1]:
        render_2d_view(mol_row)
    with view_tabs[2]:
        render_properties(mol_row)


def render_3d_view(mol_row: pd.Series):
    """Render interactive 3D molecular structure."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available for 3D view")
        return
    
    # Settings in expander
    with st.expander("⚙️ Visualization Settings", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            style = st.selectbox("Style", ["Ball-Stick", "Stick", "Sphere", "Line"], key="view3d_style")
        with c2:
            bg_color = st.color_picker("Background", "#1e1e2e", key="view3d_bg")
        with c3:
            spin = st.checkbox("Spin Animation", False, key="view3d_spin")
    
    # Build XYZ string
    mol_name = mol_row.get("molecule_id", "molecule")
    n_atoms = len(coords)
    xyz_lines = [str(n_atoms), str(mol_name)]
    
    for _, atom in coords.iterrows():
        el = str(atom.get("atom", atom.get("element", "C")))
        x = float(atom.get("x", 0))
        y = float(atom.get("y", 0))
        z = float(atom.get("z", 0))
        xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
    
    xyz_str = "\\n".join(xyz_lines)
    
    # Get Mulliken charges for hover tooltips
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    charges_js = str(charges) if charges else "[]"
    
    # Style mapping
    style_map = {
        "Ball-Stick": '{"stick":{"radius":0.12},"sphere":{"scale":0.25}}',
        "Stick": '{"stick":{"radius":0.12}}',
        "Sphere": '{"sphere":{"scale":0.35}}',
        "Line": '{"line":{}}'
    }
    style_js = style_map.get(style, style_map["Ball-Stick"])
    
    spin_js = "setInterval(function(){viewer.rotate(1,{x:0,y:1,z:0});viewer.render();},40);" if spin else ""
    
    html = f'''<!DOCTYPE html>
<html><head>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body {{ margin: 0; padding: 0; font-family: sans-serif; }}
#container {{ width: 100%; height: 460px; position: relative; }}
#tooltip {{ 
    position: absolute; 
    background: rgba(0,0,0,0.85); 
    color: white; 
    padding: 8px 12px; 
    border-radius: 6px; 
    font-size: 12px; 
    display: none; 
    pointer-events: none;
    z-index: 100;
    max-width: 200px;
}}
</style>
</head>
<body>
<div id="container"></div>
<div id="tooltip"></div>
<script>
var viewer = $3Dmol.createViewer("container", {{backgroundColor: "{bg_color}"}});
var xyz = "{xyz_str}";
viewer.addModel(xyz, "xyz");
viewer.setStyle({{}}, {style_js});

var charges = {charges_js};

viewer.setHoverCallback(function(atom, v, event) {{
    var tip = document.getElementById('tooltip');
    if (atom) {{
        var html = '<b>' + atom.elem + ' #' + atom.serial + '</b>';
        html += '<br>X: ' + atom.x.toFixed(4);
        html += '<br>Y: ' + atom.y.toFixed(4);
        html += '<br>Z: ' + atom.z.toFixed(4);
        if (charges.length > 0 && charges[atom.serial - 1] !== undefined) {{
            html += '<br>Charge: ' + charges[atom.serial - 1].toFixed(4);
        }}
        tip.innerHTML = html;
        tip.style.display = 'block';
        tip.style.left = (event.offsetX + 15) + 'px';
        tip.style.top = (event.offsetY - 10) + 'px';
    }} else {{
        tip.style.display = 'none';
    }}
}});

viewer.zoomTo();
viewer.render();
{spin_js}
</script>
</body></html>'''
    
    components.html(html, height=480)
    st.caption("💡 Hover over atoms for details. Scroll to zoom, drag to rotate.")


def render_2d_view(mol_row: pd.Series):
    """Render 2D molecular structure using RDKit."""
    
    smiles = mol_row.get("smiles")
    
    if not smiles or str(smiles) == "nan":
        st.warning("No SMILES data available for 2D view")
        return
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            st.warning(f"Cannot parse SMILES: {smiles}")
            return
        
        # Compute 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Draw as SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(450, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        # Center the SVG
        st.markdown(f'<div style="display:flex;justify-content:center;background:#f8f9fa;padding:20px;border-radius:10px;">{svg}</div>', unsafe_allow_html=True)
        
        # Show SMILES string
        with st.expander("📋 SMILES String"):
            st.code(smiles, language=None)
        
    except ImportError:
        st.info("RDKit not installed. Install with: pip install rdkit")
        st.markdown("**SMILES:**")
        st.code(smiles)


def render_properties(mol_row: pd.Series):
    """Render property tables."""
    
    prop_tabs = st.tabs(["📍 Coordinates", "⚡ Properties", "💫 Mulliken Charges", "🔧 Method"])
    
    with prop_tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            st.dataframe(coords, use_container_width=True, height=350)
        else:
            st.info("No coordinate data")
    
    with prop_tabs[1]:
        # Scalar properties table
        props = {}
        for key in ["molecule_id", "optimized_state", "smiles", "charge", "multiplicity",
                    "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy", 
                    "homo_lumo_gap", "calc_class", "esd_type"]:
            val = mol_row.get(key)
            if val is not None and str(val) != "nan":
                props[key] = val
        
        if props:
            st.dataframe(
                pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items()]),
                hide_index=True,
                use_container_width=True
            )
        else:
            st.info("No property data")
    
    with prop_tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            st.dataframe(mulliken, use_container_width=True, height=350)
        else:
            st.info("No Mulliken data")
    
    with prop_tabs[3]:
        # Method details
        method_props = {}
        for key in ["functional", "basis_set", "dispersion", "solvent", "method_id"]:
            val = mol_row.get(key)
            if val is not None and str(val) != "nan":
                method_props[key] = val
        
        if method_props:
            st.dataframe(
                pd.DataFrame([{"Property": k, "Value": v} for k, v in method_props.items()]),
                hide_index=True,
                use_container_width=True
            )
        else:
            st.info("No method data")
