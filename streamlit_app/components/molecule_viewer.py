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
    """Render interactive 3D molecular structure with enhanced visualization options."""
    
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
    
    # Get Mulliken charges for electrostatic coloring
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    
    # Settings in expander - all features enabled by default
    with st.expander("⚙️ 3D Visualization Settings", expanded=True):
        # Row 1: Style and Color
        c1, c2, c3 = st.columns(3)
        with c1:
            style = st.selectbox("Style", 
                ["Ball & Stick", "Stick", "Spacefill (CPK)", "Line (Wire)", "Cross"], 
                key="mol3d_style", index=0)
        with c2:
            color_options = ["Standard (Jmol)", "Professional (Light Carbon)", "Vibrant (Teal Carbon)"]
            if charges:
                color_options.insert(0, "Electrostatic (Charge)")
            color_mode = st.selectbox("Color Mode", color_options, key="mol3d_color", index=0)
        with c3:
            material = st.selectbox("Material", ["Shiny (Plastic)", "Matte"], key="mol3d_material", index=0)
        
        # Row 2: Labels and Surface
        c4, c5, c6 = st.columns(3)
        with c4:
            show_labels = st.checkbox("Atom Labels (Heavy)", True, key="mol3d_labels")
        with c5:
            show_surface = st.checkbox("Show Surface", False, key="mol3d_surface")
        with c6:
            spin = st.checkbox("Spin Animation", False, key="mol3d_spin")
        
        # Row 3: Surface options
        c7, c8, c9 = st.columns(3)
        with c7:
            if show_surface:
                surface_type = st.selectbox("Surface Type", 
                    ["VDW (Space-filling)", "SAS (Solvent Accessible)", "SES (Molecular)"],
                    key="mol3d_surftype")
            else:
                surface_type = "VDW (Space-filling)"
        with c8:
            if show_surface:
                surface_opacity = st.slider("Surface Opacity", 0.1, 1.0, 0.5, key="mol3d_opacity")
            else:
                surface_opacity = 0.5
        with c9:
            if show_labels:
                label_size = st.slider("Label Size", 8, 16, 10, key="mol3d_label_size")
            else:
                label_size = 10
        
        # Row 4: Advanced settings
        c10, c11, c12 = st.columns(3)
        with c10:
            bg_color = st.color_picker("Background", "#1e1e2e", key="mol3d_bg")
        with c11:
            stick_radius = st.slider("Bond Radius", 0.05, 0.25, 0.12, key="mol3d_radius")
        with c12:
            sphere_scale = st.slider("Atom Scale", 0.1, 1.0, 0.25, key="mol3d_scale")
    
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
        charge = charges[i] if i < len(charges) else 0
        atom_data.append({"el": el, "x": x, "y": y, "z": z, "idx": i + 1, "charge": charge})
    
    xyz_str = "\\n".join(xyz_lines)
    
    # Build charges array for JS
    charges_js = str(charges) if charges else "[]"
    
    # Atom data for labels (only heavy atoms)
    heavy_atoms = [a for a in atom_data if a["el"] != "H"]
    atom_data_js = str(heavy_atoms).replace("'", '"')
    
    # Style and color JavaScript
    material_spec = ", specular: 1, shininess: 50" if material == "Shiny (Plastic)" else ""
    
    # Color scheme based on mode
    if color_mode == "Electrostatic (Charge)" and charges:
        # Use charge-based coloring
        color_js = f'''
        viewer.setStyle({{}}, {{{get_style_type(style, stick_radius, sphere_scale)}: {{
            colorscheme: {{prop: 'partialCharge', gradient: 'rwb', min: -0.25, max: 0.25}}{material_spec}
        }}}});
        '''
    elif color_mode == "Professional (Light Carbon)":
        color_scheme = "{'C': '#d3d3d3', 'H': 'white', 'N': '#3050F8', 'O': '#FF0D0D', 'S': '#FFFF30', 'F': '#90E050', 'Cl': '#1FF01F'}"
        color_js = f'''
        viewer.setStyle({{}}, {{{get_style_type(style, stick_radius, sphere_scale)}: {{
            colorscheme: {{colorfunc: function(atom) {{
                var colors = {color_scheme};
                return colors[atom.elem] || '#808080';
            }}}}{material_spec}
        }}}});
        '''
    elif color_mode == "Vibrant (Teal Carbon)":
        color_scheme = "{'C': '#008080', 'H': 'white', 'N': '#483D8B', 'O': '#FF4500', 'S': '#FFD700'}"
        color_js = f'''
        viewer.setStyle({{}}, {{{get_style_type(style, stick_radius, sphere_scale)}: {{
            colorscheme: {{colorfunc: function(atom) {{
                var colors = {color_scheme};
                return colors[atom.elem] || '#808080';
            }}}}{material_spec}
        }}}});
        '''
    else:  # Standard Jmol
        style_spec = get_style_spec(style, material_spec, stick_radius, sphere_scale)
        color_js = f'viewer.setStyle({{}}, {style_spec});'
    
    # Surface JavaScript with surface type
    if show_surface:
        surf_type_map = {
            "VDW (Space-filling)": "$3Dmol.SurfaceType.VDW",
            "SAS (Solvent Accessible)": "$3Dmol.SurfaceType.SAS",
            "SES (Molecular)": "$3Dmol.SurfaceType.SES"
        }
        surf_type_js = surf_type_map.get(surface_type, "$3Dmol.SurfaceType.VDW")
        surface_js = f'viewer.addSurface({surf_type_js}, {{opacity: {surface_opacity}, color: "white"}});'
    else:
        surface_js = ""
    
    # Spin JavaScript
    spin_js = "setInterval(function(){viewer.rotate(1,{x:0,y:1,z:0});viewer.render();},40);" if spin else ""
    
    # Label JavaScript (heavy atoms only)
    text_color = 'black' if bg_color in ['#ffffff', '#f8f9fa', '#e8e8e8'] else 'white'
    if show_labels:
        label_js = f'''
        var atomData = {atom_data_js};
        for (var i = 0; i < atomData.length; i++) {{
            var a = atomData[i];
            viewer.addLabel(a.el, {{
                position: {{x: a.x, y: a.y, z: a.z}},
                fontSize: {label_size},
                fontColor: '{text_color}',
                showBackground: false
            }});
        }}
        '''
    else:
        label_js = ""
    
    # Build energy info for display
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
#container {{ width: 100%; height: 450px; position: relative; }}
#tooltip {{ 
    position: absolute; 
    background: rgba(0,0,0,0.9); 
    color: white; 
    padding: 10px 14px; 
    border-radius: 8px; 
    font-size: 12px; 
    display: none; 
    pointer-events: none;
    z-index: 100;
    max-width: 220px;
    box-shadow: 0 2px 10px rgba(0,0,0,0.3);
}}
#info-bar {{
    position: absolute;
    bottom: 8px;
    left: 8px;
    background: rgba(0,0,0,0.75);
    color: #00ff88;
    padding: 6px 12px;
    border-radius: 6px;
    font-size: 11px;
    font-family: monospace;
}}
</style>
</head>
<body>
<div id="container"></div>
<div id="tooltip"></div>
<div id="info-bar">{energy_info if energy_info else f"{n_atoms} atoms | {len(heavy_atoms)} heavy"}</div>
<script>
var viewer = $3Dmol.createViewer("container", {{backgroundColor: "{bg_color}"}});
var xyz = "{xyz_str}";
viewer.addModel(xyz, "xyz");

{color_js}

{surface_js}

var charges = {charges_js};

{label_js}

viewer.setHoverCallback(function(atom, v, event) {{
    var tip = document.getElementById('tooltip');
    if (atom) {{
        var html = '<b style="font-size:14px">' + atom.elem + ' #' + atom.serial + '</b>';
        html += '<br><span style="color:#888">Position:</span>';
        html += '<br>X: ' + atom.x.toFixed(4);
        html += '<br>Y: ' + atom.y.toFixed(4);
        html += '<br>Z: ' + atom.z.toFixed(4);
        if (charges.length > 0 && charges[atom.serial - 1] !== undefined) {{
            var q = charges[atom.serial - 1];
            var qColor = q > 0 ? '#ff6b6b' : '#4ecdc4';
            html += '<br><span style="color:' + qColor + '">Charge: ' + q.toFixed(4) + '</span>';
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
    
    components.html(html, height=490)
    st.caption("💡 Hover over atoms for details. Scroll to zoom, drag to rotate.")


def get_style_type(style: str, stick_radius: float = 0.12, sphere_scale: float = 0.25) -> str:
    """Get the 3Dmol style type string with parameters."""
    if style == "Ball & Stick":
        return f'"stick": {{"radius": {stick_radius}}}, "sphere": {{"scale": {sphere_scale}}}'
    elif style == "Spacefill (CPK)":
        return f'"sphere": {{"scale": 0.85}}'
    elif style == "Line (Wire)":
        return f'"line": {{"linewidth": 2}}'
    elif style == "Cross":
        return f'"cross": {{"radius": 0.3}}'
    else:  # Stick
        return f'"stick": {{"radius": {stick_radius}}}'


def get_style_spec(style: str, material_spec: str, stick_radius: float = 0.12, sphere_scale: float = 0.25) -> str:
    """Get the full style specification for 3Dmol."""
    if style == "Ball & Stick":
        return f'{{"stick": {{"radius": {stick_radius}{material_spec}}}, "sphere": {{"scale": {sphere_scale}{material_spec}}}}}'
    elif style == "Spacefill (CPK)":
        return f'{{"sphere": {{"scale": 0.85{material_spec}}}}}'
    elif style == "Line (Wire)":
        return f'{{"line": {{"linewidth": 2}}}}'
    elif style == "Cross":
        return f'{{"cross": {{"radius": 0.3{material_spec}}}}}'
    else:  # Stick
        return f'{{"stick": {{"radius": {stick_radius}{material_spec}}}}}'


def render_2d_view(mol_row: pd.Series):
    """Render 2D molecular structure using RDKit with extensive customization."""
    
    smiles = mol_row.get("smiles")
    coords = mol_row.get("cart_coords")
    
    # Settings with more RDKit options
    with st.expander("⚙️ 2D View Settings", expanded=True):
        # Row 1: Basic settings
        c1, c2, c3 = st.columns(3)
        with c1:
            img_size = st.slider("Image Size", 300, 700, 450, key="mol2d_size")
        with c2:
            bg_color = st.color_picker("Background", "#ffffff", key="mol2d_bg")
        with c3:
            color_palette = st.selectbox("Color Palette", 
                ["Standard (RDKit)", "Black & White", "Avalon", "Dark Mode"],
                key="mol2d_palette")
        
        # Row 2: Display options
        c4, c5, c6 = st.columns(3)
        with c4:
            show_atom_indices = st.checkbox("Atom Indices", False, key="mol2d_indices")
        with c5:
            show_stereo = st.checkbox("Stereo Labels", True, key="mol2d_stereo")
        with c6:
            kekulize = st.checkbox("Kekulize Bonds", False, key="mol2d_kekulize")
        
        # Row 3: Line settings
        c7, c8, c9 = st.columns(3)
        with c7:
            bond_width = st.slider("Bond Width", 1.0, 4.0, 2.0, key="mol2d_bond")
        with c8:
            atom_font = st.slider("Atom Font Size", 10, 24, 14, key="mol2d_fontsize")
        with c9:
            padding = st.slider("Padding", 0.0, 0.2, 0.05, key="mol2d_padding")
    
    mol = None
    source = None
    
    # Try SMILES first
    if smiles and str(smiles) != "nan":
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                source = "SMILES"
        except ImportError:
            pass
        except Exception:
            pass
    
    # Fallback: Create molecule from XYZ coordinates using RDKit
    if mol is None and coords is not None and not (hasattr(coords, 'empty') and coords.empty):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # Build XYZ block
            n_atoms = len(coords)
            xyz_lines = [str(n_atoms), "molecule"]
            
            for _, atom in coords.iterrows():
                el = str(atom.get("atom", atom.get("element", "C")))
                x = float(atom.get("x", 0))
                y = float(atom.get("y", 0))
                z = float(atom.get("z", 0))
                xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
            
            xyz_block = "\n".join(xyz_lines)
            
            # Create molecule from XYZ using rdDetermineBonds
            raw_mol = Chem.MolFromXYZBlock(xyz_block)
            if raw_mol is not None:
                try:
                    # Try to determine bonds
                    from rdkit.Chem import rdDetermineBonds
                    rdDetermineBonds.DetermineBonds(raw_mol)
                    mol = raw_mol
                    source = "XYZ coordinates"
                except Exception:
                    # Fallback: use raw mol without bonds
                    mol = raw_mol
                    source = "XYZ coordinates (no bonds)"
        except ImportError:
            st.warning("RDKit not installed. Install with: pip install rdkit")
            return
        except Exception as e:
            st.warning(f"Failed to create molecule from coordinates: {e}")
    
    if mol is None:
        st.warning("No molecular data available for 2D view")
        return
    
    # Render with RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        
        # Kekulize if requested
        if kekulize:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception:
                pass
        
        # Compute 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Setup drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(img_size, int(img_size * 0.88))
        opts = drawer.drawOptions()
        
        # Apply display options
        opts.addAtomIndices = show_atom_indices
        opts.addStereoAnnotation = show_stereo
        opts.bondLineWidth = bond_width
        opts.padding = padding
        
        # Font size - set min and max to same value for consistent sizing
        opts.minFontSize = atom_font
        opts.maxFontSize = atom_font + 4
        
        # Apply color palette
        if color_palette == "Black & White":
            opts.useBWAtomPalette()
        elif color_palette == "Avalon":
            try:
                opts.useAvalonAtomPalette()
            except Exception:
                pass  # Fallback to default if not available
        elif color_palette == "Dark Mode":
            # Custom dark mode palette
            opts.setBackgroundColour((0.1, 0.1, 0.15, 1.0))
            # Set atom colors for dark background
            opts.updateAtomPalette({
                6: (0.9, 0.9, 0.9),   # Carbon - light gray
                7: (0.3, 0.5, 1.0),   # Nitrogen - blue
                8: (1.0, 0.3, 0.3),   # Oxygen - red
                16: (1.0, 1.0, 0.3),  # Sulfur - yellow
                9: (0.3, 1.0, 0.3),   # Fluorine - green
                17: (0.3, 0.9, 0.3),  # Chlorine - green
                35: (0.6, 0.3, 0.3),  # Bromine - brown
                15: (1.0, 0.5, 0.0),  # Phosphorus - orange
            })
        
        # Draw molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        # Apply background color (override for non-dark mode)
        if color_palette != "Dark Mode":
            svg = svg.replace('style="background-color:#FFFFFF;', f'style="background-color:{bg_color};')
        
        # Display
        display_bg = "#1a1a24" if color_palette == "Dark Mode" else bg_color
        st.markdown(f'<div style="display:flex;justify-content:center;background:{display_bg};padding:20px;border-radius:10px;">{svg}</div>', unsafe_allow_html=True)
        
        st.caption(f"📐 Source: {source}")
        
        if smiles and str(smiles) != "nan":
            with st.expander("📋 SMILES String"):
                st.code(smiles, language=None)
                
    except Exception as e:
        st.error(f"RDKit rendering failed: {e}")


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
