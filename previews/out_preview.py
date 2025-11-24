"""
ORCA Output file preview using Streamlit and Plotly.

This module provides visualization for ORCA calculation results.
"""

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.out_parser import parse_out_file, parse_out_content


def show_energy_plot(energies: list[float], title: str = "Energy Convergence"):
    """Plot energy convergence during optimization."""
    if not energies:
        st.warning("No energy data to plot")
        return

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=list(range(1, len(energies) + 1)),
        y=energies,
        mode='lines+markers',
        name='Energy'
    ))

    fig.update_layout(
        title=title,
        xaxis_title="Cycle",
        yaxis_title="Energy (Hartree)",
        template="plotly_white"
    )

    st.plotly_chart(fig, use_container_width=True)


def show_frequency_histogram(frequencies: list[float]):
    """Plot histogram of vibrational frequencies."""
    if not frequencies:
        st.warning("No frequency data to plot")
        return

    fig = px.histogram(
        x=frequencies,
        nbins=50,
        title="Vibrational Frequency Distribution"
    )
    fig.update_layout(
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Count",
        template="plotly_white"
    )

    st.plotly_chart(fig, use_container_width=True)


def show_mulliken_charges(charges: dict):
    """Display Mulliken charges as a bar chart."""
    if not charges:
        st.warning("No Mulliken charge data")
        return

    # Prepare data
    data = []
    for idx, (element, charge) in charges.items():
        data.append({
            'Atom': f"{element}{idx}",
            'Element': element,
            'Charge': charge
        })

    df = pd.DataFrame(data)

    fig = px.bar(
        df, x='Atom', y='Charge',
        color='Element',
        title="Mulliken Atomic Charges"
    )
    fig.update_layout(template="plotly_white")

    st.plotly_chart(fig, use_container_width=True)


def out_preview_page():
    """Main Streamlit page for ORCA output file preview."""
    st.title("ORCA Output Viewer")

    # File upload
    uploaded_file = st.file_uploader("Upload ORCA output file", type=['out'])

    # Or use sample file
    sample_files = ['p1xs0.out']
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        selected_sample = st.selectbox("Select sample", sample_files)
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            selected_sample
        )

        result = parse_out_file(filepath)

    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        result = parse_out_content(content)
    else:
        st.info("Please upload an ORCA output file or select a sample file.")
        return

    # Job Information
    st.header("Job Information")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Basis Set", result.job_info.basis_set or "N/A")
    with col2:
        st.metric("Charge", result.job_info.charge)
    with col3:
        st.metric("Multiplicity", result.job_info.multiplicity)

    # Energy Summary
    st.header("Energy")
    col1, col2 = st.columns(2)
    with col1:
        if result.final_energy:
            st.metric("Final Energy (Eh)", f"{result.final_energy:.6f}")
    with col2:
        if result.thermochemistry and result.thermochemistry.gibbs_free_energy:
            st.metric("Gibbs Free Energy (Eh)",
                     f"{result.thermochemistry.gibbs_free_energy:.6f}")

    # Dipole Moment
    if result.dipole_moment:
        st.header("Dipole Moment")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("X (a.u.)", f"{result.dipole_moment.x:.4f}")
        with col2:
            st.metric("Y (a.u.)", f"{result.dipole_moment.y:.4f}")
        with col3:
            st.metric("Z (a.u.)", f"{result.dipole_moment.z:.4f}")
        with col4:
            st.metric("Total (Debye)", f"{result.dipole_moment.magnitude_debye:.3f}")

    # Energy Convergence Plot
    if result.optimization_energies:
        st.header("Optimization Convergence")
        show_energy_plot(result.optimization_energies, "Energy during Optimization")

    # Vibrational Frequencies
    if result.frequencies:
        st.header("Vibrational Analysis")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Number of Modes", len(result.frequencies))
        with col2:
            st.metric("Lowest (cm⁻¹)", f"{min(result.frequencies):.1f}")
        with col3:
            st.metric("Highest (cm⁻¹)", f"{max(result.frequencies):.1f}")

        show_frequency_histogram(result.frequencies)

        # Show frequencies table
        with st.expander("View all frequencies"):
            freq_df = pd.DataFrame({
                'Mode': list(range(1, len(result.frequencies) + 1)),
                'Frequency (cm⁻¹)': result.frequencies
            })
            st.dataframe(freq_df, use_container_width=True)

    # Mulliken Charges
    if result.mulliken_charges:
        st.header("Mulliken Charges")
        show_mulliken_charges(result.mulliken_charges)

    # Thermochemistry
    if result.thermochemistry:
        st.header("Thermochemistry")
        thermo = result.thermochemistry
        col1, col2 = st.columns(2)
        with col1:
            st.write(f"**Temperature:** {thermo.temperature} K")
            st.write(f"**Pressure:** {thermo.pressure} atm")
            st.write(f"**Total Mass:** {thermo.total_mass} AMU")
        with col2:
            if thermo.zero_point_energy:
                st.write(f"**ZPE:** {thermo.zero_point_energy:.6f} Eh")
            if thermo.enthalpy:
                st.write(f"**Enthalpy:** {thermo.enthalpy:.6f} Eh")
            if thermo.gibbs_free_energy:
                st.write(f"**Gibbs Free Energy:** {thermo.gibbs_free_energy:.6f} Eh")

    # Export option
    if st.button("Export as JSON"):
        import json
        json_data = json.dumps(result.to_dict(), indent=2)
        st.download_button(
            "Download JSON",
            json_data,
            "orca_output.json",
            "application/json"
        )


if __name__ == '__main__':
    out_preview_page()
