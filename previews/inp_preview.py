"""
Input file preview using Streamlit.
"""

import streamlit as st
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.inp_parser import parse_inp_file, parse_inp_content


def inp_preview_page():
    """Main Streamlit page for input file preview."""
    st.title("ORCA Input File Viewer")

    uploaded_file = st.file_uploader("Upload input file", type=['inp'])
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.inp'
        )
        data = parse_inp_file(filepath)
    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        data = parse_inp_content(content)
    else:
        st.info("Please upload an ORCA input file or select a sample file.")
        return

    # Job info
    st.header("Calculation Settings")
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Method", data.method or "N/A")
        st.metric("Charge", data.charge)
        st.metric("MaxCore (MB)", data.maxcore or "Default")
    with col2:
        st.metric("Basis Set", data.basis_set or "N/A")
        st.metric("Multiplicity", data.multiplicity)
        st.metric("Processors", data.nprocs)

    # Job types
    if data.job_types:
        st.write(f"**Job Types:** {', '.join(data.job_types)}")

    # Keywords
    st.header("Keywords")
    st.code(' '.join(data.keywords), language=None)

    # Coordinates
    if data.coordinates:
        st.header("Molecular Geometry")
        st.metric("Number of Atoms", len(data.coordinates))

        coord_data = []
        for i, (elem, x, y, z) in enumerate(data.coordinates):
            coord_data.append({
                'Index': i,
                'Element': elem,
                'X (Å)': f"{x:.6f}",
                'Y (Å)': f"{y:.6f}",
                'Z (Å)': f"{z:.6f}"
            })

        df = pd.DataFrame(coord_data)
        st.dataframe(df, use_container_width=True)

    # Export
    if st.button("Export as JSON"):
        import json
        st.download_button(
            "Download JSON",
            json.dumps(data.to_dict(), indent=2),
            "input.json",
            "application/json"
        )


if __name__ == '__main__':
    inp_preview_page()
