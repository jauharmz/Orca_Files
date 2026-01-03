"""
Molecule Viewer Component - 2D and 3D visualization with data hover

Features:
- 2D structure from SMILES (RDKit or fallback)
- 3D structure from coordinates (py3Dmol)
- Hover/click for atom properties
- Multi-molecule comparison
"""

import streamlit as st
import pandas as pd
from typing import Optional, Dict, List
import streamlit.components.v1 as components


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer tab."""
    
    st.header("🧬 Molecule Viewer")
    
    # Molecule selector
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    col1, col2 = st.columns([1, 3])
    
    with col1:
        selected_mol = st.selectbox("Select Molecule", mol_ids, key="mol_viewer_select")
        view_mode = st.radio("View Mode", ["3D Structure", "2D Structure", "Both"], key="mol_view_mode")
    
    # Get molecule data
    mol_row = df[df["molecule_id"] == selected_mol].iloc[0] if selected_mol else None
    
    if mol_row is None:
        st.warning("Select a molecule to view")
        return
    
    with col2:
        # Show molecule info
        info_col1, info_col2, info_col3 = st.columns(3)
        with info_col1:
            st.metric("Molecule", selected_mol)
        with info_col2:
            smiles = mol_row.get("smiles", "N/A")
            st.metric("SMILES", smiles[:20] + "..." if smiles and len(str(smiles)) > 20 else smiles)
        with info_col3:
            charge = mol_row.get("charge", "N/A")
            mult = mol_row.get("multiplicity", "N/A")
            st.metric("Charge/Mult", f"{charge}/{mult}")
    
    st.divider()
    
    # Render based on view mode
    if view_mode in ["3D Structure", "Both"]:
        render_3d_viewer(mol_row, selected_mol)
    
    if view_mode in ["2D Structure", "Both"]:
        render_2d_viewer(mol_row, selected_mol)
    
    # Data section
    st.divider()
    render_molecule_data(mol_row, selected_mol, df)


def render_3d_viewer(mol_row: pd.Series, mol_id: str):
    """Render 3D molecular structure."""
    
    st.subheader("🔮 3D Structure")
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available")
        return
    
    # Style options
    col1, col2, col3 = st.columns(3)
    with col1:
        style = st.selectbox("Style", ["ball-stick", "stick", "sphere", "line"], key="3d_style")
    with col2:
        bg_color = st.color_picker("Background", "#FFFFFF", key="3d_bg")
    with col3:
        show_labels = st.checkbox("Atom Labels", False, key="3d_labels")
    
    # Create XYZ string
    xyz_lines = [str(len(coords)), mol_id]
    for _, atom in coords.iterrows():
        element = atom.get("atom", atom.get("element", "X"))
        x = atom.get("x", 0)
        y = atom.get("y", 0)
        z = atom.get("z", 0)
        xyz_lines.append(f"{element}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_str = "\\n".join(xyz_lines)
    
    # Mulliken charges for hover (if available)
    mulliken = mol_row.get("mulliken")
    charges_js = "null"
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
            charges_js = str(charges)
    
    # Style configuration
    style_config = {
        "ball-stick": '{"stick": {"radius": 0.15}, "sphere": {"scale": 0.3}}',
        "stick": '{"stick": {"radius": 0.15}}',
        "sphere": '{"sphere": {"scale": 0.4}}',
        "line": '{"line": {}}'
    }[style]
    
    label_js = ""
    if show_labels:
        label_js = """
        var atoms = viewer.selectedAtoms({});
        for (var i = 0; i < atoms.length; i++) {
            viewer.addLabel(atoms[i].elem, {position: atoms[i], backgroundColor: 'white', fontColor: 'black', fontSize: 12});
        }
        """
    
    # py3Dmol HTML
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
            #viewer {{ width: 100%; height: 500px; position: relative; }}
            .tooltip {{
                position: absolute;
                background: rgba(0,0,0,0.8);
                color: white;
                padding: 5px 10px;
                border-radius: 4px;
                font-size: 12px;
                pointer-events: none;
                display: none;
            }}
        </style>
    </head>
    <body>
        <div id="viewer"></div>
        <div id="tooltip" class="tooltip"></div>
        <script>
            var viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "{bg_color}"}});
            var xyz = "{xyz_str}";
            viewer.addModel(xyz, "xyz");
            viewer.setStyle({{}}, {style_config});
            {label_js}
            
            // Hover callback
            var charges = {charges_js};
            viewer.setHoverCallback(function(atom, viewer, event, container) {{
                var tooltip = document.getElementById('tooltip');
                if (atom) {{
                    var text = atom.elem + " (" + atom.serial + ")";
                    text += "\\nX: " + atom.x.toFixed(3);
                    text += "\\nY: " + atom.y.toFixed(3);
                    text += "\\nZ: " + atom.z.toFixed(3);
                    if (charges && charges[atom.serial - 1] !== undefined) {{
                        text += "\\nCharge: " + charges[atom.serial - 1].toFixed(4);
                    }}
                    tooltip.innerHTML = text.replace(/\\n/g, '<br>');
                    tooltip.style.display = 'block';
                    tooltip.style.left = event.clientX + 'px';
                    tooltip.style.top = event.clientY + 'px';
                }} else {{
                    tooltip.style.display = 'none';
                }}
            }});
            
            viewer.zoomTo();
            viewer.render();
        </script>
    </body>
    </html>
    """
    
    components.html(html, height=550)


def render_2d_viewer(mol_row: pd.Series, mol_id: str):
    """Render 2D molecular structure."""
    
    st.subheader("📐 2D Structure")
    
    smiles = mol_row.get("smiles")
    
    if not smiles:
        st.warning("No SMILES data available for 2D rendering")
        return
    
    # Try to use RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        import io
        import base64
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            st.warning(f"Could not parse SMILES: {smiles}")
            return
        
        # Options
        col1, col2 = st.columns(2)
        with col1:
            show_atom_nums = st.checkbox("Show Atom Numbers", False, key="2d_atom_nums")
        with col2:
            img_size = st.slider("Image Size", 200, 600, 400, key="2d_size")
        
        # Draw
        if show_atom_nums:
            for atom in mol.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
        
        img = Draw.MolToImage(mol, size=(img_size, img_size))
        
        # Convert to base64
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        img_str = base64.b64encode(buffer.getvalue()).decode()
        
        st.image(f"data:image/png;base64,{img_str}", caption=f"{mol_id}: {smiles}")
        
    except ImportError:
        st.info("RDKit not installed. Showing SMILES text instead.")
        st.code(smiles, language=None)
        st.markdown(f"[View on PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{smiles})")


def render_molecule_data(mol_row: pd.Series, mol_id: str, df: pd.DataFrame):
    """Render molecule data tables."""
    
    st.subheader(f"📊 Data for {mol_id}")
    
    tabs = st.tabs(["Coordinates", "Properties", "Mulliken"])
    
    with tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            from components.data_table import render_simple_table
            render_simple_table(coords, "Cartesian Coordinates")
        else:
            st.info("No coordinate data")
    
    with tabs[1]:
        # Scalar properties
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
    
    with tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            from components.data_table import render_simple_table
            render_simple_table(mulliken, "Mulliken Charges")
        else:
            st.info("No Mulliken charge data")
