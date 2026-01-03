"""
Molecule Viewer Component - 2D and 3D visualization with data hover

Features:
- 2D structure from SMILES (RDKit or fallback)
- 3D structure from coordinates (py3Dmol)
- Hover/click for atom properties
- Animation and styling options

Fixed: Uses fragment=True to prevent Streamlit rerun issues
"""

import streamlit as st
import pandas as pd
from typing import Optional, Dict, List
import streamlit.components.v1 as components


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer tab."""
    
    st.subheader("🧬 Molecule Viewer")
    
    # Molecule selector - show with state info
    if df.empty:
        st.warning("No molecules available")
        return
    
    # Create display labels with state info
    mol_options = []
    for _, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        state_label = f" [{state}]" if state else ""
        mol_options.append(f"{mol_id}{state_label}")
    
    # Get unique options
    unique_options = list(dict.fromkeys(mol_options))
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    
    # Selector with index
    if "mol_viewer_idx" not in st.session_state:
        st.session_state.mol_viewer_idx = 0
    
    selected_idx = st.selectbox(
        "Select Molecule",
        range(len(mol_ids)),
        format_func=lambda i: f"{mol_ids[i]} ({df[df['molecule_id']==mol_ids[i]]['optimized_state'].iloc[0] if 'optimized_state' in df.columns else ''})",
        key="mol_viewer_select_idx"
    )
    
    selected_mol = mol_ids[selected_idx]
    mol_row = df[df["molecule_id"] == selected_mol].iloc[0]
    
    # Show molecule info
    info_cols = st.columns(5)
    with info_cols[0]:
        st.metric("Molecule", selected_mol)
    with info_cols[1]:
        state = mol_row.get("optimized_state", "N/A")
        st.metric("State", state)
    with info_cols[2]:
        energy = mol_row.get("gibbs_Eh") or mol_row.get("single_point_Eh")
        st.metric("Energy (Eh)", f"{energy:.4f}" if energy else "N/A")
    with info_cols[3]:
        homo = mol_row.get("homo_energy")
        st.metric("HOMO (eV)", f"{homo:.2f}" if homo else "N/A")
    with info_cols[4]:
        method = mol_row.get("method_id", "N/A")
        st.metric("Method", str(method)[:12] if method else "N/A")
    
    # View mode as tabs
    view_tabs = st.tabs(["🔮 3D Structure", "📐 2D Structure"])
    
    with view_tabs[0]:
        render_3d_viewer(mol_row, selected_mol)
    
    with view_tabs[1]:
        render_2d_viewer(mol_row, selected_mol)
    
    # Data section - collapsed by default
    with st.expander("📊 Molecule Data", expanded=False):
        render_molecule_data(mol_row, selected_mol, df)


def render_3d_viewer(mol_row: pd.Series, mol_id: str):
    """Render 3D molecular structure - all settings in one render."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available for 3D view")
        return
    
    # All settings in expander
    with st.expander("⚙️ 3D Settings", expanded=False):
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            style = st.selectbox("Style", ["ball-stick", "stick", "sphere", "line"], key="style_3d")
        with col2:
            bg_color = st.color_picker("Background", "#FFFFFF", key="bg_3d")
        with col3:
            show_labels = st.checkbox("Labels", False, key="labels_3d")
        with col4:
            spin = st.checkbox("Spin", False, key="spin_3d")
    
    # Build XYZ string
    xyz_lines = [str(len(coords)), mol_id]
    for _, atom in coords.iterrows():
        element = atom.get("atom", atom.get("element", "X"))
        x, y, z = atom.get("x", 0), atom.get("y", 0), atom.get("z", 0)
        xyz_lines.append(f"{element}  {x:.6f}  {y:.6f}  {z:.6f}")
    xyz_str = "\\n".join(xyz_lines)
    
    # Mulliken charges
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    charges_js = str(charges) if charges else "[]"
    
    # Style config
    style_configs = {
        "ball-stick": '{"stick": {"radius": 0.15, "colorscheme": "Jmol"}, "sphere": {"scale": 0.3, "colorscheme": "Jmol"}}',
        "stick": '{"stick": {"radius": 0.15, "colorscheme": "Jmol"}}',
        "sphere": '{"sphere": {"scale": 0.4, "colorscheme": "Jmol"}}',
        "line": '{"line": {"colorscheme": "Jmol"}}'
    }
    style_config = style_configs.get(style, style_configs["ball-stick"])
    
    spin_js = "setInterval(function(){viewer.rotate(1,{x:0,y:1,z:0});viewer.render();},50);" if spin else ""
    label_js = """
        var atoms = viewer.selectedAtoms({});
        for(var i=0;i<atoms.length;i++){
            viewer.addLabel(atoms[i].elem,{position:atoms[i],backgroundColor:'#333',fontColor:'white',fontSize:12});
        }
    """ if show_labels else ""
    
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <style>
            body{{margin:0;padding:0;}}
            #v{{width:100%;height:480px;position:relative;}}
            #tip{{position:fixed;background:rgba(0,0,0,0.9);color:white;padding:8px 12px;border-radius:6px;font-size:12px;display:none;z-index:1000;font-family:sans-serif;}}
        </style>
    </head>
    <body>
        <div id="v"></div>
        <div id="tip"></div>
        <script>
            var viewer=$3Dmol.createViewer("v",{{backgroundColor:"{bg_color}"}});
            viewer.addModel("{xyz_str}","xyz");
            viewer.setStyle({{}},{style_config});
            {label_js}
            var charges={charges_js};
            viewer.setHoverCallback(function(atom,v,ev){{
                var t=document.getElementById('tip');
                if(atom){{
                    var h='<b>'+atom.elem+' #'+atom.serial+'</b><br>X:'+atom.x.toFixed(3)+'<br>Y:'+atom.y.toFixed(3)+'<br>Z:'+atom.z.toFixed(3);
                    if(charges[atom.serial-1]!==undefined)h+='<br>Charge:'+charges[atom.serial-1].toFixed(4);
                    t.innerHTML=h;t.style.display='block';t.style.left=(ev.clientX+15)+'px';t.style.top=(ev.clientY-10)+'px';
                }}else{{t.style.display='none';}}
            }});
            viewer.zoomTo();viewer.render();
            {spin_js}
        </script>
    </body>
    </html>
    """
    
    components.html(html, height=500)
    st.caption("💡 Hover over atoms for details")


def render_2d_viewer(mol_row: pd.Series, mol_id: str):
    """Render 2D molecular structure."""
    
    smiles = mol_row.get("smiles")
    
    if not smiles or str(smiles) == "nan":
        st.warning("No SMILES data available for 2D view")
        return
    
    # Settings
    with st.expander("⚙️ 2D Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            show_nums = st.checkbox("Atom Numbers", False, key="nums_2d")
        with col2:
            img_size = st.slider("Size", 300, 500, 350, key="size_2d")
    
    # Try RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            st.warning(f"Could not parse SMILES: {smiles}")
            return
        
        AllChem.Compute2DCoords(mol)
        
        if show_nums:
            for atom in mol.GetAtoms():
                atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
        
        drawer = rdMolDraw2D.MolDraw2DSVG(img_size, img_size)
        drawer.drawOptions().addAtomIndices = show_nums
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        st.markdown(f'<div style="display:flex;justify-content:center;">{svg}</div>', unsafe_allow_html=True)
        st.caption(f"**SMILES:** `{smiles}`")
        
    except ImportError:
        st.info("RDKit not installed. Showing SMILES text.")
        st.code(smiles)


def render_molecule_data(mol_row: pd.Series, mol_id: str, df: pd.DataFrame):
    """Render molecule data tables."""
    
    tabs = st.tabs(["📍 Coordinates", "📊 Properties", "⚡ Mulliken"])
    
    with tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            st.dataframe(coords, use_container_width=True, height=250)
        else:
            st.info("No coordinate data")
    
    with tabs[1]:
        props = {
            "molecule_id": mol_row.get("molecule_id"),
            "optimized_state": mol_row.get("optimized_state"),
            "method_id": mol_row.get("method_id"),
            "functional": mol_row.get("functional"),
            "basis_set": mol_row.get("basis_set"),
            "gibbs_Eh": mol_row.get("gibbs_Eh"),
            "single_point_Eh": mol_row.get("single_point_Eh"),
            "homo_energy": mol_row.get("homo_energy"),
            "lumo_energy": mol_row.get("lumo_energy"),
        }
        props_df = pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items() if v is not None])
        st.dataframe(props_df, use_container_width=True, hide_index=True)
    
    with tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            st.dataframe(mulliken, use_container_width=True, height=250)
        else:
            st.info("No Mulliken data")
