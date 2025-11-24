"""
ORCA Output Viewer - Main Streamlit Application

Run with: streamlit run app.py
"""

import streamlit as st
import os
import sys

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

st.set_page_config(
    page_title="ORCA Output Viewer",
    page_icon="🔬",
    layout="wide"
)

st.title("ORCA Output Viewer")
st.markdown("Parse and visualize ORCA quantum chemistry output files")

# Sidebar navigation
st.sidebar.title("File Types")

page = st.sidebar.radio(
    "Select viewer:",
    [
        "XYZ - Molecular Geometry",
        "OUT - Main Output",
        "SPECTRUM - Spectra",
        "HESS - IR Spectrum",
        "INP - Input File",
        "ENGRAD - Energy & Gradient",
        "PROPERTY - Charges",
    ]
)

st.sidebar.markdown("---")
st.sidebar.markdown("### Status")
st.sidebar.markdown("✅ All 7 parsers complete")

# Main content
if page == "XYZ - Molecular Geometry":
    from previews.xyz_preview import xyz_preview_page
    xyz_preview_page()
elif page == "OUT - Main Output":
    from previews.out_preview import out_preview_page
    out_preview_page()
elif page == "SPECTRUM - Spectra":
    from previews.spectrum_preview import spectrum_preview_page
    spectrum_preview_page()
elif page == "HESS - IR Spectrum":
    from previews.hess_preview import hess_preview_page
    hess_preview_page()
elif page == "INP - Input File":
    from previews.inp_preview import inp_preview_page
    inp_preview_page()
elif page == "ENGRAD - Energy & Gradient":
    from previews.engrad_preview import engrad_preview_page
    engrad_preview_page()
elif page == "PROPERTY - Charges":
    from previews.property_preview import property_preview_page
    property_preview_page()

# Footer
st.markdown("---")
st.markdown("*ORCA Output Parser & Viewer - Work in Progress*")
