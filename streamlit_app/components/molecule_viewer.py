"""
Molecule Viewer Component - 2D and 3D visualization

Features:
- Per-molecule filter with state
- 3D viewer using py3Dmol (simplified)
- 2D viewer using RDKit
- Data tables
"""

import streamlit as st
import pandas as pd
from typing import Optional
import streamlit.components.v1 as components


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer with its own filters."""
    
    st.subheader("🧬 Molecule Viewer")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build options: molecule_id + state
    options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        options.append({"label": label, "idx": idx, "mol_id": mol_id, "state": state})
    
    # Unique labels
    seen = set()
    unique_options = []
    for opt in options:
        if opt["label"] not in seen:
            seen.add(opt["label"])
            unique_options.append(opt)
    
    if not unique_options:
        st.warning("No molecules available")
        return
    
    # Selector
    labels = [opt["label"] for opt in unique_options]
    selected_label = st.selectbox("Select Molecule", labels, key="mol_viewer_select")
    
    # Get selected row
    selected_opt = next((opt for opt in unique_options if opt["label"] == selected_label), None)
    if selected_opt is None:
        return
    
    mol_row = df.loc[selected_opt["idx"]]
    
    # Info
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Molecule", mol_row.get("molecule_id", "N/A"))
    with c2:
        st.metric("State", mol_row.get("optimized_state", "N/A"))
    with c3:
        e = mol_row.get("gibbs_Eh") or mol_row.get("single_point_Eh")
        st.metric("Energy (Eh)", f"{e:.4f}" if e else "N/A")
    with c4:
        h = mol_row.get("homo_energy")
        st.metric("HOMO (eV)", f"{h:.2f}" if h else "N/A")
    
    # View tabs
    tabs = st.tabs(["🔮 3D View", "📐 2D View", "📊 Data"])
    
    with tabs[0]:
        render_3d_simple(mol_row)
    with tabs[1]:
        render_2d(mol_row)
    with tabs[2]:
        render_data(mol_row)


def render_3d_simple(mol_row: pd.Series):
    """Simple 3D viewer - minimal settings to avoid rerun issues."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available")
        return
    
    # Build XYZ string
    mol_name = mol_row.get("molecule_id", "mol")
    n_atoms = len(coords)
    xyz_lines = [str(n_atoms), str(mol_name)]
    
    for _, atom in coords.iterrows():
        el = str(atom.get("atom", atom.get("element", "C")))
        x = float(atom.get("x", 0))
        y = float(atom.get("y", 0))
        z = float(atom.get("z", 0))
        xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
    
    xyz_str = "\\n".join(xyz_lines)
    
    # Get Mulliken charges for hover
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    charges_js = str(charges) if charges else "[]"
    
    # Simple HTML - no reactive settings
    html = f'''<!DOCTYPE html>
<html>
<head>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; }}
        #container {{ width: 100%; height: 450px; position: relative; }}
        #tip {{ 
            position: absolute; 
            background: rgba(0,0,0,0.85); 
            color: white; 
            padding: 8px 12px; 
            border-radius: 6px; 
            font-size: 12px; 
            display: none; 
            pointer-events: none;
            font-family: sans-serif;
            z-index: 100;
        }}
    </style>
</head>
<body>
    <div id="container"></div>
    <div id="tip"></div>
    <script>
        var viewer = $3Dmol.createViewer("container", {{backgroundColor: "#1e1e2e"}});
        var xyz = "{xyz_str}";
        viewer.addModel(xyz, "xyz");
        viewer.setStyle({{}}, {{stick: {{radius: 0.12}}, sphere: {{scale: 0.25}}}});
        
        var charges = {charges_js};
        
        viewer.setHoverCallback(function(atom, v, event) {{
            var tip = document.getElementById('tip');
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
    </script>
</body>
</html>'''
    
    components.html(html, height=470)
    st.caption("💡 Hover over atoms for details. Scroll to zoom, drag to rotate.")


def render_2d(mol_row: pd.Series):
    """Render 2D structure."""
    
    smiles = mol_row.get("smiles")
    
    if not smiles or str(smiles) == "nan":
        st.warning("No SMILES data available")
        return
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            st.warning(f"Cannot parse SMILES: {smiles}")
            return
        
        AllChem.Compute2DCoords(mol)
        
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        st.markdown(f'<div style="display:flex;justify-content:center;">{svg}</div>', unsafe_allow_html=True)
        st.code(smiles, language=None)
        
    except ImportError:
        st.info("RDKit not installed. SMILES string:")
        st.code(smiles)


def render_data(mol_row: pd.Series):
    """Render data tables."""
    
    tabs = st.tabs(["Coordinates", "Properties", "Mulliken"])
    
    with tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            st.dataframe(coords, use_container_width=True, height=300)
        else:
            st.info("No coordinate data")
    
    with tabs[1]:
        props = {}
        for key in ["molecule_id", "optimized_state", "method_id", "functional", "basis_set",
                    "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy", "homo_lumo_gap",
                    "charge", "multiplicity"]:
            val = mol_row.get(key)
            if val is not None and str(val) != "nan":
                props[key] = val
        
        if props:
            st.dataframe(
                pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items()]),
                use_container_width=True,
                hide_index=True
            )
        else:
            st.info("No property data")
    
    with tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            st.dataframe(mulliken, use_container_width=True, height=300)
        else:
            st.info("No Mulliken data")
