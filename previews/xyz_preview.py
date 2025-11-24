"""
XYZ file preview using Streamlit and py3Dmol.

This module provides visualization for molecular geometries.
"""

import streamlit as st
import py3Dmol
from streamlit.components.v1 import html
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.xyz_parser import parse_xyz_file, parse_xyz_content, parse_trajectory_xyz


def show_molecule_3d(xyz_content: str, width: int = 600, height: int = 400):
    """
    Display 3D molecular structure using py3Dmol.

    Args:
        xyz_content: XYZ format string
        width: Viewer width in pixels
        height: Viewer height in pixels
    """
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(xyz_content, 'xyz')
    viewer.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
    viewer.zoomTo()

    # Generate HTML
    viewer_html = viewer._make_html()

    html(viewer_html, width=width, height=height)


def show_atom_table(geom):
    """Display atom information in a table."""
    import pandas as pd

    data = []
    for i, atom in enumerate(geom.atoms):
        data.append({
            'Index': i + 1,
            'Element': atom.symbol,
            'X (Å)': f"{atom.x:.6f}",
            'Y (Å)': f"{atom.y:.6f}",
            'Z (Å)': f"{atom.z:.6f}"
        })

    df = pd.DataFrame(data)
    st.dataframe(df, use_container_width=True)


def xyz_preview_page():
    """Main Streamlit page for XYZ file preview."""
    st.title("XYZ File Viewer")

    # File upload
    uploaded_file = st.file_uploader("Upload XYZ file", type=['xyz'])

    # Or use sample file
    sample_files = [
        'p1xs0.xyz',
        'p1xs0_trj.xyz'
    ]

    use_sample = st.checkbox("Use sample file")

    if use_sample:
        selected_sample = st.selectbox("Select sample", sample_files)
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            selected_sample
        )

        with open(filepath, 'r') as f:
            content = f.read()

        # Check if trajectory
        if '_trj' in selected_sample:
            frames = parse_trajectory_xyz(filepath)
            st.info(f"Trajectory file with {len(frames)} frames")

            frame_idx = st.slider("Frame", 0, len(frames) - 1, 0)
            geom = frames[frame_idx]

            # Reconstruct XYZ content for this frame
            lines = [str(geom.num_atoms), geom.comment]
            for atom in geom.atoms:
                lines.append(f"{atom.symbol} {atom.x} {atom.y} {atom.z}")
            content = '\n'.join(lines)
        else:
            geom = parse_xyz_file(filepath)

    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        geom = parse_xyz_content(content)
    else:
        st.info("Please upload an XYZ file or select a sample file.")
        return

    # Display info
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Number of Atoms", geom.num_atoms)
    with col2:
        if geom.energy:
            st.metric("Energy (Hartree)", f"{geom.energy:.6f}")

    st.write(f"**Comment:** {geom.comment}")

    # 3D viewer
    st.subheader("3D Structure")
    show_molecule_3d(content)

    # Atom table
    st.subheader("Atom Coordinates")
    show_atom_table(geom)

    # Export option
    if st.button("Export as JSON"):
        import json
        json_data = json.dumps(geom.to_dict(), indent=2)
        st.download_button(
            "Download JSON",
            json_data,
            "molecule.json",
            "application/json"
        )


if __name__ == '__main__':
    xyz_preview_page()
