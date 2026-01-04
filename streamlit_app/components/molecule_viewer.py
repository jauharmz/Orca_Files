"""
Molecule Viewer Component - 2D/3D molecular visualization

Features:
- Molecule + State selector
- 3D viewer with py3Dmol
- 2D viewer with RDKit
- Property tables
"""

import streamlit as st
import pandas as pd
import streamlit.components.v1 as components


def render_molecule_viewer(df: pd.DataFrame):
    """Render molecule viewer with molecule+state filter."""
    
    st.subheader("🧬 Molecule Viewer")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build unique (molecule_id, state) pairs for selector
    pairs = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        pairs.append({"label": label, "idx": idx, "mol_id": mol_id, "state": state})
    
    # Unique labels
    seen = set()
    unique = [p for p in pairs if p["label"] not in seen and not seen.add(p["label"])]
    
    if not unique:
        st.warning("No molecules available")
        return
    
    # Dual filter
    c1, c2 = st.columns(2)
    with c1:
        mol_ids = sorted(set(p["mol_id"] for p in unique))
        sel_mol = st.selectbox("Molecule", mol_ids, key="mol_view_mol")
    with c2:
        mol_states = sorted([p["state"] for p in unique if p["mol_id"] == sel_mol])
        sel_state = st.selectbox("State", mol_states, key="mol_view_state") if mol_states else None
    
    # Get selected row
    match = [p for p in unique if p["mol_id"] == sel_mol and p["state"] == sel_state]
    if not match:
        st.warning("No data for selection")
        return
    
    mol_row = df.loc[match[0]["idx"]]
    
    # Info metrics
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
        render_3d(mol_row)
    with tabs[1]:
        render_2d(mol_row)
    with tabs[2]:
        render_data(mol_row)


def render_3d(mol_row: pd.Series):
    """Render 3D molecular structure."""
    
    coords = mol_row.get("cart_coords")
    
    if coords is None or (hasattr(coords, 'empty') and coords.empty):
        st.warning("No coordinate data")
        return
    
    # Build XYZ
    mol_name = mol_row.get("molecule_id", "mol")
    xyz_lines = [str(len(coords)), str(mol_name)]
    
    for _, atom in coords.iterrows():
        el = str(atom.get("atom", atom.get("element", "C")))
        x, y, z = float(atom.get("x", 0)), float(atom.get("y", 0)), float(atom.get("z", 0))
        xyz_lines.append(f"{el} {x:.6f} {y:.6f} {z:.6f}")
    
    xyz_str = "\\n".join(xyz_lines)
    
    # Mulliken charges for hover
    mulliken = mol_row.get("mulliken")
    charges = []
    if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
        if "Charge" in mulliken.columns:
            charges = mulliken["Charge"].tolist()
    charges_js = str(charges) if charges else "[]"
    
    html = f'''<!DOCTYPE html><html><head>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>body{{margin:0;}}#v{{width:100%;height:450px;}}
#t{{position:absolute;background:#222;color:#fff;padding:6px 10px;border-radius:4px;font-size:11px;display:none;z-index:99;}}</style>
</head><body><div id="v"></div><div id="t"></div><script>
var v=$3Dmol.createViewer("v",{{backgroundColor:"#1e1e2e"}});
v.addModel("{xyz_str}","xyz");
v.setStyle({{}},{{stick:{{radius:0.12}},sphere:{{scale:0.25}}}});
var q={charges_js};
v.setHoverCallback(function(a,vw,e){{var t=document.getElementById('t');
if(a){{var h='<b>'+a.elem+' #'+a.serial+'</b><br>'+a.x.toFixed(4)+', '+a.y.toFixed(4)+', '+a.z.toFixed(4);
if(q[a.serial-1]!==undefined)h+='<br>Charge: '+q[a.serial-1].toFixed(4);
t.innerHTML=h;t.style.display='block';t.style.left=(e.offsetX+15)+'px';t.style.top=(e.offsetY-10)+'px';
}}else{{t.style.display='none';}}}});
v.zoomTo();v.render();
</script></body></html>'''
    
    components.html(html, height=470)
    st.caption("Hover atoms for details. Scroll=zoom, drag=rotate.")


def render_2d(mol_row: pd.Series):
    """Render 2D structure."""
    
    smiles = mol_row.get("smiles")
    
    if not smiles or str(smiles) == "nan":
        st.warning("No SMILES data")
        return
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
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
        st.info("RDKit not installed")
        st.code(smiles)


def render_data(mol_row: pd.Series):
    """Render data tables."""
    
    tabs = st.tabs(["Coordinates", "Properties", "Mulliken"])
    
    with tabs[0]:
        coords = mol_row.get("cart_coords")
        if coords is not None and hasattr(coords, 'empty') and not coords.empty:
            st.dataframe(coords, use_container_width=True, height=300)
        else:
            st.info("No coordinates")
    
    with tabs[1]:
        props = {k: mol_row.get(k) for k in [
            "molecule_id", "optimized_state", "functional", "basis_set",
            "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy"
        ] if mol_row.get(k) is not None and str(mol_row.get(k)) != "nan"}
        
        if props:
            st.dataframe(pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items()]),
                        hide_index=True, use_container_width=True)
        else:
            st.info("No properties")
    
    with tabs[2]:
        mulliken = mol_row.get("mulliken")
        if mulliken is not None and hasattr(mulliken, 'empty') and not mulliken.empty:
            st.dataframe(mulliken, use_container_width=True, height=300)
        else:
            st.info("No Mulliken data")
