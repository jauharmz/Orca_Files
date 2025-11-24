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
        # Future parsers will be added here
        # "OUT - Main Output",
        # "HESS - Frequencies",
        # "SPECTRUM - Spectra",
    ]
)

st.sidebar.markdown("---")
st.sidebar.markdown("### Status")
st.sidebar.markdown("✅ XYZ Parser")
st.sidebar.markdown("⏳ More parsers coming...")

# Main content
if page == "XYZ - Molecular Geometry":
    from previews.xyz_preview import xyz_preview_page
    xyz_preview_page()

# Footer
st.markdown("---")
st.markdown("*ORCA Output Parser & Viewer - Work in Progress*")
