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
    """Render interactive 3D molecular structure with atom labels."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available for 3D view")
        return
    
    # Get molecule info for display
    mol_id = mol_row.get("molecule_id", "molecule")
    state = mol_row.get("optimized_state", "")
    energy = mol_row.get("gibbs_Eh") or mol_row.get("single_point_Eh")
    homo = mol_row.get("homo_energy")
    lumo = mol_row.get("lumo_energy")
    
    # Settings in expander - use unique stable keys
    with st.expander("⚙️ 3D Visualization Settings", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            style = st.selectbox("Style", ["Ball-Stick", "Stick", "Sphere", "Line"], 
                               key="mol3d_style", index=0)
            show_labels = st.checkbox("Show Atom Labels", False, key="mol3d_labels")
        with c2:
            bg_color = st.color_picker("Background", "#1e1e2e", key="mol3d_bg")
            spin = st.checkbox("Spin Animation", False, key="mol3d_spin")
        
        c3, c4 = st.columns(2)
        with c3:
            label_size = st.slider("Label Size", 8, 20, 12, key="mol3d_label_size") if show_labels else 12
        with c4:
            show_indices = st.checkbox("Show Atom Indices", False, key="mol3d_indices") if show_labels else False
    
    # Build XYZ string
    n_atoms = len(coords)
    xyz_lines = [str(n_atoms), str(mol_id)]
    atom_data = []
    
    for i, (_, atom) in enumerate(coords.iterrows()):
        el = str(atom.get("atom", atom.get("element", "C")))
        x = float(atom.get("x", 0))
        y = float(atom.get("y", 0))
        z = float(atom.get("z", 0))
        xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
        atom_data.append({"el": el, "x": x, "y": y, "z": z, "idx": i + 1})
    
    xyz_str = "\\n".join(xyz_lines)
    
    # Get Mulliken charges for hover tooltips
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    charges_js = str(charges) if charges else "[]"
    
    # Atom data for labels
    atom_data_js = str(atom_data).replace("'", '"')
    
    # Style mapping
    style_map = {
        "Ball-Stick": '{"stick":{"radius":0.12},"sphere":{"scale":0.25}}',
        "Stick": '{"stick":{"radius":0.12}}',
        "Sphere": '{"sphere":{"scale":0.35}}',
        "Line": '{"line":{}}'
    }
    style_js = style_map.get(style, style_map["Ball-Stick"])
    
    spin_js = "setInterval(function(){viewer.rotate(1,{x:0,y:1,z:0});viewer.render();},40);" if spin else ""
    
    # Label generation JavaScript
    if show_labels:
        label_js = f'''
        var atomData = {atom_data_js};
        for (var i = 0; i < atomData.length; i++) {{
            var a = atomData[i];
            var labelText = {"'a.el + a.idx'" if show_indices else "'a.el'"};
            viewer.addLabel(a.el{' + a.idx' if show_indices else ''}, {{
                position: {{x: a.x, y: a.y, z: a.z}},
                fontSize: {label_size},
                fontColor: 'white',
                backgroundColor: 'rgba(0,0,0,0.5)',
                backgroundOpacity: 0.5,
                showBackground: true
            }});
        }}
        '''
    else:
        label_js = ""
    
    # Build energy info for display in viewer
    energy_info = ""
    if energy:
        energy_info += f"Energy: {energy:.6f} Eh"
    if homo:
        energy_info += f" | HOMO: {homo:.2f} eV"
    if lumo:
        energy_info += f" | LUMO: {lumo:.2f} eV"
    
    html = f'''<!DOCTYPE html>
<html><head>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body {{ margin: 0; padding: 0; font-family: sans-serif; }}
#container {{ width: 100%; height: 440px; position: relative; }}
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
#info-bar {{
    position: absolute;
    bottom: 5px;
    left: 5px;
    background: rgba(0,0,0,0.7);
    color: #00ff88;
    padding: 4px 10px;
    border-radius: 4px;
    font-size: 11px;
    font-family: monospace;
}}
</style>
</head>
<body>
<div id="container"></div>
<div id="tooltip"></div>
<div id="info-bar">{energy_info if energy_info else f"{n_atoms} atoms"}</div>
<script>
var viewer = $3Dmol.createViewer("container", {{backgroundColor: "{bg_color}"}});
var xyz = "{xyz_str}";
viewer.addModel(xyz, "xyz");
viewer.setStyle({{}}, {style_js});

var charges = {charges_js};

{label_js}

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
    """Render 2D molecular structure - uses SMILES if available, otherwise coordinates."""
    
    smiles = mol_row.get("smiles")
    coords = mol_row.get("cart_coords")
    
    # Settings
    with st.expander("⚙️ 2D View Settings", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            img_size = st.slider("Image Size", 300, 600, 450, key="mol2d_size")
            show_atom_indices = st.checkbox("Show Atom Indices", False, key="mol2d_indices")
        with c2:
            bg_color = st.color_picker("Background", "#f8f9fa", key="mol2d_bg")
            bond_width = st.slider("Bond Width", 1, 5, 2, key="mol2d_bond")
    
    # Try SMILES first
    if smiles and str(smiles) != "nan":
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw, AllChem
            from rdkit.Chem.Draw import rdMolDraw2D
            
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                # Compute 2D coordinates
                AllChem.Compute2DCoords(mol)
                
                # Draw options
                drawer = rdMolDraw2D.MolDraw2DSVG(img_size, int(img_size * 0.88))
                opts = drawer.drawOptions()
                opts.addAtomIndices = show_atom_indices
                opts.bondLineWidth = bond_width
                
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                
                # Display
                st.markdown(f'<div style="display:flex;justify-content:center;background:{bg_color};padding:20px;border-radius:10px;">{svg}</div>', unsafe_allow_html=True)
                
                with st.expander("📋 SMILES String"):
                    st.code(smiles, language=None)
                return
        except ImportError:
            st.info("RDKit not installed. Using coordinate-based 2D view.")
        except Exception as e:
            st.warning(f"SMILES parsing failed: {e}. Using coordinate-based 2D view.")
    
    # Fallback: Use coordinates for 2D view
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No molecular data available for 2D view (no SMILES or coordinates)")
        return
    
    st.info("📐 2D projection from 3D coordinates (SMILES not available)")
    
    # Create SVG from coordinates (project XY plane)
    render_2d_from_coords(coords, img_size, bg_color, show_atom_indices)


def render_2d_from_coords(coords, img_size: int, bg_color: str, show_indices: bool):
    """Render 2D view by projecting 3D coordinates onto XY plane."""
    
    # Element colors
    element_colors = {
        'H': '#FFFFFF', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D',
        'S': '#FFFF30', 'P': '#FF8000', 'F': '#90E050', 'Cl': '#1FF01F',
        'Br': '#A62929', 'I': '#940094', 'B': '#FFB5B5', 'Si': '#F0C8A0'
    }
    
    # Extract coordinates
    atoms = []
    x_vals, y_vals = [], []
    
    for _, atom in coords.iterrows():
        el = str(atom.get("atom", atom.get("element", "C")))
        x = float(atom.get("x", 0))
        y = float(atom.get("y", 0))
        atoms.append({"el": el, "x": x, "y": y})
        x_vals.append(x)
        y_vals.append(y)
    
    if not atoms:
        st.warning("No atom data")
        return
    
    # Scale to fit
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    
    margin = 40
    width = img_size
    height = int(img_size * 0.88)
    
    x_range = max(x_max - x_min, 0.1)
    y_range = max(y_max - y_min, 0.1)
    scale = min((width - 2*margin) / x_range, (height - 2*margin) / y_range)
    
    def scale_x(x):
        return margin + (x - x_min) * scale + (width - 2*margin - x_range * scale) / 2
    
    def scale_y(y):
        return height - margin - (y - y_min) * scale - (height - 2*margin - y_range * scale) / 2
    
    # Build SVG
    svg_parts = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    svg_parts.append(f'<rect width="100%" height="100%" fill="{bg_color}"/>')
    
    # Draw bonds (simple distance-based)
    bond_threshold = 2.0  # Angstrom
    for i, a1 in enumerate(atoms):
        for j, a2 in enumerate(atoms):
            if i < j:
                dx = atoms[i]["x"] - atoms[j]["x"]
                dy = atoms[i]["y"] - atoms[j]["y"]
                dist = (dx*dx + dy*dy) ** 0.5
                if dist < bond_threshold:
                    x1, y1 = scale_x(a1["x"]), scale_y(a1["y"])
                    x2, y2 = scale_x(a2["x"]), scale_y(a2["y"])
                    svg_parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#666" stroke-width="2"/>')
    
    # Draw atoms
    for i, atom in enumerate(atoms):
        x, y = scale_x(atom["x"]), scale_y(atom["y"])
        color = element_colors.get(atom["el"], '#808080')
        radius = 12 if atom["el"] != 'H' else 8
        
        svg_parts.append(f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{color}" stroke="#333" stroke-width="1"/>')
        
        # Label
        label = f"{atom['el']}{i+1}" if show_indices else atom["el"]
        font_size = 10 if show_indices else 11
        text_color = '#000' if color in ['#FFFFFF', '#FFFF30', '#FFB5B5', '#F0C8A0'] else '#FFF'
        svg_parts.append(f'<text x="{x}" y="{y+4}" text-anchor="middle" font-size="{font_size}" font-family="sans-serif" fill="{text_color}">{label}</text>')
    
    svg_parts.append('</svg>')
    svg = '\n'.join(svg_parts)
    
    # Display
    st.markdown(f'<div style="display:flex;justify-content:center;padding:10px;">{svg}</div>', unsafe_allow_html=True)
    st.caption("📐 XY projection of 3D coordinates. Bond detection: distance < 2.0 Å")


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
