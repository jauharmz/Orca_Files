"""
Molecule Viewer Component - 2D and 3D visualization

Uses static HTML rendering to avoid Streamlit rerun issues.
Settings are applied directly without causing widget conflicts.
"""

import streamlit as st
import pandas as pd
from typing import Optional
import streamlit.components.v1 as components
import hashlib


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer tab."""
    
    st.subheader("🧬 Molecule Viewer")
    
    if df.empty:
        st.warning("No molecules available")
        return
    
    # Molecule selector
    mol_ids = sorted(df["molecule_id"].dropna().unique().tolist())
    
    if not mol_ids:
        st.warning("No molecules available")
        return
    
    # Build options with state
    options = []
    for _, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        options.append((label, row.name))  # (label, index)
    
    # Remove duplicates
    seen = set()
    unique_options = []
    for label, idx in options:
        if label not in seen:
            seen.add(label)
            unique_options.append((label, idx))
    
    labels = [opt[0] for opt in unique_options]
    indices = [opt[1] for opt in unique_options]
    
    selected_label = st.selectbox("Select Molecule", labels, key="mol_view_select")
    selected_idx = indices[labels.index(selected_label)]
    mol_row = df.loc[selected_idx]
    
    # Info row
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Molecule", mol_row.get("molecule_id", "N/A"))
    with col2:
        st.metric("State", mol_row.get("optimized_state", "N/A"))
    with col3:
        energy = mol_row.get("gibbs_Eh") or mol_row.get("single_point_Eh")
        st.metric("Energy (Eh)", f"{energy:.4f}" if energy else "N/A")
    with col4:
        homo = mol_row.get("homo_energy")
        st.metric("HOMO (eV)", f"{homo:.2f}" if homo else "N/A")
    
    # View tabs
    view_tabs = st.tabs(["🔮 3D Structure", "📐 2D Structure", "📊 Data"])
    
    with view_tabs[0]:
        render_3d_static(mol_row)
    
    with view_tabs[1]:
        render_2d_static(mol_row)
    
    with view_tabs[2]:
        render_data_tables(mol_row)


def render_3d_static(mol_row: pd.Series):
    """Render 3D viewer with inline settings (no widget-based updates)."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data available")
        return
    
    # Settings in columns (values directly embedded in HTML)
    st.markdown("**Settings:**")
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        style = st.selectbox("Style", ["Ball-Stick", "Stick", "Sphere", "Line"], key="s3d")
    with c2:
        bg = st.color_picker("Background", "#1a1a2e", key="bg3d")
    with c3:
        spin = st.checkbox("Spin", False, key="spin3d")
    with c4:
        labels = st.checkbox("Labels", False, key="lbl3d")
    
    # Build XYZ
    mol_id = mol_row.get("molecule_id", "molecule")
    xyz_lines = [str(len(coords)), str(mol_id)]
    for _, atom in coords.iterrows():
        el = atom.get("atom", atom.get("element", "C"))
        x, y, z = atom.get("x", 0), atom.get("y", 0), atom.get("z", 0)
        xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
    xyz_str = "\\n".join(xyz_lines)
    
    # Mulliken charges
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    
    # Style map
    style_map = {
        "Ball-Stick": '{"stick":{"radius":0.12},"sphere":{"scale":0.25}}',
        "Stick": '{"stick":{"radius":0.12}}',
        "Sphere": '{"sphere":{"scale":0.35}}',
        "Line": '{"line":{}}'
    }
    style_js = style_map.get(style, style_map["Ball-Stick"])
    
    spin_js = "setInterval(()=>{v.rotate(1,{x:0,y:1,z:0});v.render();},40);" if spin else ""
    label_js = """
        v.selectedAtoms({}).forEach(a=>{
            v.addLabel(a.elem,{position:a,backgroundColor:'#333',fontColor:'#fff',fontSize:11});
        });
    """ if labels else ""
    
    charges_json = str(charges) if charges else "[]"
    
    # Create unique key for this molecule+settings to force iframe refresh
    settings_hash = hashlib.md5(f"{mol_id}{style}{bg}{spin}{labels}".encode()).hexdigest()[:8]
    
    html = f"""<!DOCTYPE html>
<html><head>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body{{margin:0;font-family:sans-serif;}}
#m{{width:100%;height:460px;}}
#t{{position:fixed;background:#222;color:#fff;padding:6px 10px;border-radius:4px;font-size:11px;display:none;z-index:99;}}
</style>
</head><body>
<div id="m"></div>
<div id="t"></div>
<script>
var v=$3Dmol.createViewer("m",{{backgroundColor:"{bg}"}});
v.addModel("{xyz_str}","xyz");
v.setStyle({{}},{style_js});
{label_js}
var q={charges_json};
v.setHoverCallback(function(a,vw,e){{
    var t=document.getElementById('t');
    if(a){{
        var h='<b>'+a.elem+' #'+a.serial+'</b><br>'+a.x.toFixed(3)+', '+a.y.toFixed(3)+', '+a.z.toFixed(3);
        if(q[a.serial-1]!==undefined)h+='<br>Charge: '+q[a.serial-1].toFixed(4);
        t.innerHTML=h;t.style.display='block';
        t.style.left=(e.clientX+10)+'px';t.style.top=(e.clientY-10)+'px';
    }}else{{t.style.display='none';}}
}});
v.zoomTo();v.render();
{spin_js}
</script>
</body></html>"""
    
    components.html(html, height=480)
    st.caption("Hover atoms for details")


def render_2d_static(mol_row: pd.Series):
    """Render 2D structure."""
    
    smiles = mol_row.get("smiles")
    
    if not smiles or str(smiles) == "nan":
        st.warning("No SMILES data")
        return
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            st.warning(f"Invalid SMILES: {smiles}")
            return
        
        AllChem.Compute2DCoords(mol)
        
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        st.markdown(f'<div style="display:flex;justify-content:center;">{svg}</div>', unsafe_allow_html=True)
        st.code(smiles, language=None)
        
    except ImportError:
        st.info("RDKit not available. SMILES:")
        st.code(smiles)


def render_data_tables(mol_row: pd.Series):
    """Render data tables for molecule."""
    
    tabs = st.tabs(["Coordinates", "Properties", "Mulliken"])
    
    with tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            st.dataframe(coords, use_container_width=True, height=250)
        else:
            st.info("No coordinates")
    
    with tabs[1]:
        props = {k: mol_row.get(k) for k in [
            "molecule_id", "optimized_state", "method_id", "functional", 
            "basis_set", "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy"
        ] if mol_row.get(k) is not None}
        if props:
            st.dataframe(pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items()]), 
                        hide_index=True, use_container_width=True)
        else:
            st.info("No properties")
    
    with tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            st.dataframe(mulliken, use_container_width=True, height=250)
        else:
            st.info("No Mulliken data")
