"""CPCM solvation file preview."""

import streamlit as st
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.cpcm_parser import parse_cpcm_file, parse_cpcm_content


def cpcm_preview_page():
    st.title("CPCM Solvation Viewer")

    uploaded_file = st.file_uploader("Upload CPCM file", type=['cpcm'])
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.cpcm')
        data = parse_cpcm_file(filepath)
    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        data = parse_cpcm_content(content)
    else:
        st.info("Please upload a CPCM file or select a sample file.")
        return

    st.header("Solvation Summary")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Surface Points", data.num_surface_points)
        st.metric("Volume (bohr³)", f"{data.volume:.2f}")
    with col2:
        st.metric("Area (bohr²)", f"{data.area:.2f}")
        st.metric("Atoms", data.num_atoms)

    st.header("Solvation Energy")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Dielectric Energy (Eh)", f"{data.dielectric_energy:.6f}")
        st.metric("Dielectric Energy (kcal/mol)", f"{data.dielectric_energy * 627.509:.2f}")
    with col2:
        st.metric("One-electron Energy (Eh)", f"{data.one_electron_energy:.6f}")

    if st.button("Export as JSON"):
        import json
        st.download_button("Download JSON", json.dumps(data.to_dict(), indent=2), "cpcm.json", "application/json")


if __name__ == '__main__':
    cpcm_preview_page()
