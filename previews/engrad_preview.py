"""
Engrad file preview using Streamlit.
"""

import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.engrad_parser import parse_engrad_file, parse_engrad_content


def show_gradient_plot(data):
    """Plot gradient components per atom."""
    if not data.atoms:
        return

    atoms = []
    grad_x = []
    grad_y = []
    grad_z = []

    for i, atom in enumerate(data.atoms):
        atoms.append(f"{i}")
        grad_x.append(atom.grad_x)
        grad_y.append(atom.grad_y)
        grad_z.append(atom.grad_z)

    fig = go.Figure()
    fig.add_trace(go.Bar(name='X', x=atoms, y=grad_x))
    fig.add_trace(go.Bar(name='Y', x=atoms, y=grad_y))
    fig.add_trace(go.Bar(name='Z', x=atoms, y=grad_z))

    fig.update_layout(
        title="Gradient Components per Atom",
        xaxis_title="Atom Index",
        yaxis_title="Gradient (Eh/bohr)",
        barmode='group',
        template="plotly_white"
    )

    st.plotly_chart(fig, use_container_width=True)


def engrad_preview_page():
    """Main Streamlit page for engrad file preview."""
    st.title("Energy & Gradient Viewer")

    uploaded_file = st.file_uploader("Upload engrad file", type=['engrad'])
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.engrad'
        )
        data = parse_engrad_file(filepath)
    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        data = parse_engrad_content(content)
    else:
        st.info("Please upload an engrad file or select a sample file.")
        return

    # Display metrics
    st.header("Energy & Gradient Summary")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Energy (Eh)", f"{data.energy:.6f}")
    with col2:
        st.metric("Max Gradient", f"{data.max_gradient:.2e}")
    with col3:
        st.metric("RMS Gradient", f"{data.rms_gradient:.2e}")

    st.metric("Number of Atoms", data.num_atoms)

    # Gradient plot
    st.header("Gradient Analysis")
    show_gradient_plot(data)

    # Atom table
    st.header("Atomic Data")
    atom_data = []
    for i, atom in enumerate(data.atoms):
        atom_data.append({
            'Index': i,
            'Z': atom.atomic_number,
            'X (bohr)': f"{atom.x:.4f}",
            'Y (bohr)': f"{atom.y:.4f}",
            'Z (bohr)': f"{atom.z:.4f}",
            'Grad X': f"{atom.grad_x:.2e}",
            'Grad Y': f"{atom.grad_y:.2e}",
            'Grad Z': f"{atom.grad_z:.2e}"
        })

    df = pd.DataFrame(atom_data)
    st.dataframe(df, use_container_width=True)

    # Export
    if st.button("Export as JSON"):
        import json
        st.download_button(
            "Download JSON",
            json.dumps(data.to_dict(), indent=2),
            "engrad.json",
            "application/json"
        )


if __name__ == '__main__':
    engrad_preview_page()
