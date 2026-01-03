"""
Molecule Viewer Component - 2D and 3D visualization with data hover

Features:
- 2D structure from SMILES (RDKit or fallback)
- 3D structure from coordinates (py3Dmol)
- Hover/click for atom properties
- Animation and styling options
"""

import streamlit as st
import pandas as pd
from typing import Optional, Dict, List
import streamlit.components.v1 as components


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer tab."""
    
    st.subheader("🧬 Molecule Viewer")
    
    # Molecule selector
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    if not mol_ids:
        st.warning("No molecules available")
        return
    
    # Top row: molecule selector
    selected_mol = st.selectbox("Select Molecule", mol_ids, key="mol_viewer_select")
    
    # Get molecule data
    mol_row = df[df["molecule_id"] == selected_mol].iloc[0] if selected_mol else None
    
    if mol_row is None:
        st.warning("Select a molecule to view")
        return
    
    # Show molecule info
    info_col1, info_col2, info_col3, info_col4 = st.columns(4)
    with info_col1:
        st.metric("Molecule", selected_mol)
    with info_col2:
        smiles = mol_row.get("smiles", "N/A")
        st.metric("SMILES", smiles[:15] + "..." if smiles and len(str(smiles)) > 15 else smiles)
    with info_col3:
        energy = mol_row.get("gibbs_Eh") or mol_row.get("single_point_Eh")
        st.metric("Energy (Eh)", f"{energy:.4f}" if energy else "N/A")
    with info_col4:
        homo = mol_row.get("homo_energy")
        st.metric("HOMO (eV)", f"{homo:.2f}" if homo else "N/A")
    
    # View mode as tabs (instead of radio)
    view_tabs = st.tabs(["🔮 3D Structure", "📐 2D Structure"])
    
    with view_tabs[0]:
        render_3d_viewer(mol_row, selected_mol)
    
    with view_tabs[1]:
        render_2d_viewer(mol_row, selected_mol)
    
    # Data section - collapsed by default
    with st.expander("📊 Molecule Data", expanded=False):
        render_molecule_data(mol_row, selected_mol, df)


def render_3d_viewer(mol_row: pd.Series, mol_id: str):
    """Render 3D molecular structure with persistent settings."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available")
        return
    
    # Settings in expander for cleaner UI
    with st.expander("⚙️ 3D Settings", expanded=False):
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            style = st.selectbox("Style", ["ball-stick", "stick", "sphere", "line"], 
                               key=f"3d_style_{mol_id}", index=0)
        with col2:
            bg_color = st.color_picker("Background", "#FFFFFF", key=f"3d_bg_{mol_id}")
        with col3:
            show_labels = st.checkbox("Atom Labels", False, key=f"3d_labels_{mol_id}")
        with col4:
            spin = st.checkbox("Spin Animation", False, key=f"3d_spin_{mol_id}")
    
    # Create XYZ string
    xyz_lines = [str(len(coords)), mol_id]
    atom_data = []
    for i, (_, atom) in enumerate(coords.iterrows()):
        element = atom.get("atom", atom.get("element", "X"))
        x = atom.get("x", 0)
        y = atom.get("y", 0)
        z = atom.get("z", 0)
        xyz_lines.append(f"{element}  {x:.6f}  {y:.6f}  {z:.6f}")
        atom_data.append({"element": element, "x": x, "y": y, "z": z})
    
    xyz_str = "\\n".join(xyz_lines)
    
    # Mulliken charges for hover (if available)
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    
    # Style configuration
    style_config = {
        "ball-stick": '{"stick": {"radius": 0.15, "colorscheme": "Jmol"}, "sphere": {"scale": 0.3, "colorscheme": "Jmol"}}',
        "stick": '{"stick": {"radius": 0.15, "colorscheme": "Jmol"}}',
        "sphere": '{"sphere": {"scale": 0.4, "colorscheme": "Jmol"}}',
        "line": '{"line": {"colorscheme": "Jmol"}}'
    }[style]
    
    # Spin animation JS
    spin_js = """
        function spin() {
            viewer.rotate(1, {x:0, y:1, z:0});
            viewer.render();
            requestAnimationFrame(spin);
        }
        spin();
    """ if spin else ""
    
    # Label JS
    label_js = """
        var atoms = viewer.selectedAtoms({});
        for (var i = 0; i < atoms.length; i++) {
            viewer.addLabel(atoms[i].elem, {
                position: atoms[i], 
                backgroundColor: 'rgba(0,0,0,0.8)', 
                fontColor: 'white', 
                fontSize: 14,
                borderRadius: 4
            });
        }
    """ if show_labels else ""
    
    # Build charges array for JS
    charges_js = str(charges) if charges else "[]"
    
    # py3Dmol HTML with enhanced features
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
            body {{ margin: 0; padding: 0; }}
            #viewer {{ width: 100%; height: 500px; position: relative; }}
            #tooltip {{
                position: fixed;
                background: rgba(0,0,0,0.9);
                color: white;
                padding: 8px 12px;
                border-radius: 6px;
                font-size: 13px;
                pointer-events: none;
                display: none;
                z-index: 1000;
                font-family: 'Segoe UI', sans-serif;
                line-height: 1.4;
                box-shadow: 0 2px 8px rgba(0,0,0,0.3);
            }}
            #info {{
                position: absolute;
                top: 10px;
                left: 10px;
                background: rgba(255,255,255,0.9);
                padding: 8px 12px;
                border-radius: 6px;
                font-size: 12px;
                font-family: 'Segoe UI', sans-serif;
                box-shadow: 0 1px 4px rgba(0,0,0,0.2);
            }}
        </style>
    </head>
    <body>
        <div id="viewer"></div>
        <div id="tooltip"></div>
        <script>
            var viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "{bg_color}"}});
            var xyz = "{xyz_str}";
            viewer.addModel(xyz, "xyz");
            viewer.setStyle({{}}, {style_config});
            {label_js}
            
            // Charges data
            var charges = {charges_js};
            
            // Enhanced hover callback with popup data
            viewer.setHoverCallback(function(atom, viewer, event, container) {{
                var tooltip = document.getElementById('tooltip');
                if (atom) {{
                    var html = '<b>' + atom.elem + ' #' + (atom.serial) + '</b>';
                    html += '<br>X: ' + atom.x.toFixed(4) + ' Å';
                    html += '<br>Y: ' + atom.y.toFixed(4) + ' Å';
                    html += '<br>Z: ' + atom.z.toFixed(4) + ' Å';
                    if (charges.length > 0 && charges[atom.serial - 1] !== undefined) {{
                        var charge = charges[atom.serial - 1];
                        var chargeColor = charge > 0 ? '#ff6b6b' : '#4ecdc4';
                        html += '<br>Charge: <span style="color:' + chargeColor + '">' + charge.toFixed(4) + '</span>';
                    }}
                    tooltip.innerHTML = html;
                    tooltip.style.display = 'block';
                    tooltip.style.left = (event.clientX + 15) + 'px';
                    tooltip.style.top = (event.clientY - 10) + 'px';
                }} else {{
                    tooltip.style.display = 'none';
                }}
            }});
            
            // Click callback for atom selection
            viewer.setClickCallback(function(atom, viewer, event) {{
                if (atom) {{
                    viewer.addLabel(atom.elem + ' (' + atom.serial + ')', {{
                        position: atom,
                        backgroundColor: 'rgba(50,50,50,0.9)',
                        fontColor: 'white',
                        fontSize: 14,
                        borderRadius: 4
                    }});
                    viewer.render();
                }}
            }});
            
            viewer.zoomTo();
            viewer.render();
            
            {spin_js}
        </script>
    </body>
    </html>
    """
    
    components.html(html, height=520)
    
    st.caption("💡 Hover over atoms to see details. Click to add labels.")


def render_2d_viewer(mol_row: pd.Series, mol_id: str):
    """Render 2D molecular structure."""
    
    smiles = mol_row.get("smiles")
    
    if not smiles or str(smiles) == "nan":
        st.warning("No SMILES data available for 2D rendering")
        return
    
    # Settings
    with st.expander("⚙️ 2D Settings", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1:
            show_atom_nums = st.checkbox("Atom Numbers", False, key=f"2d_nums_{mol_id}")
        with col2:
            show_stereo = st.checkbox("Stereo", True, key=f"2d_stereo_{mol_id}")
        with col3:
            img_size = st.slider("Size", 300, 600, 400, key=f"2d_size_{mol_id}")
    
    # Try RDKit first
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        import io
        import base64
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            st.warning(f"Could not parse SMILES: {smiles}")
            return
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Atom numbering
        if show_atom_nums:
            for atom in mol.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
        
        # Draw with enhanced styling
        drawer = rdMolDraw2D.MolDraw2DSVG(img_size, img_size)
        drawer.drawOptions().addAtomIndices = show_atom_nums
        drawer.drawOptions().addStereoAnnotation = show_stereo
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        # Display SVG
        st.markdown(f"""
        <div style="display: flex; justify-content: center; padding: 20px;">
            {svg}
        </div>
        """, unsafe_allow_html=True)
        
        st.caption(f"**SMILES:** `{smiles}`")
        
    except ImportError:
        st.info("RDKit not installed. Showing SMILES text instead.")
        st.code(smiles, language=None)
        st.markdown(f"[🔗 View on PubChem](https://pubchem.ncbi.nlm.nih.gov/#query={smiles})")


def render_molecule_data(mol_row: pd.Series, mol_id: str, df: pd.DataFrame):
    """Render molecule data tables."""
    
    data_tabs = st.tabs(["📍 Coordinates", "📊 Properties", "⚡ Mulliken"])
    
    with data_tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            st.dataframe(coords, use_container_width=True, height=300)
        else:
            st.info("No coordinate data")
    
    with data_tabs[1]:
        props = {
            "molecule_id": mol_row.get("molecule_id"),
            "method_id": mol_row.get("method_id"),
            "functional": mol_row.get("functional"),
            "basis_set": mol_row.get("basis_set"),
            "gibbs_Eh": mol_row.get("gibbs_Eh"),
            "single_point_Eh": mol_row.get("single_point_Eh"),
            "homo_energy": mol_row.get("homo_energy"),
            "lumo_energy": mol_row.get("lumo_energy"),
            "homo_lumo_gap": mol_row.get("homo_lumo_gap"),
            "charge": mol_row.get("charge"),
            "multiplicity": mol_row.get("multiplicity"),
        }
        props_df = pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items() if v is not None])
        st.dataframe(props_df, use_container_width=True, hide_index=True)
    
    with data_tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            st.dataframe(mulliken, use_container_width=True, height=300)
        else:
            st.info("No Mulliken charge data")
