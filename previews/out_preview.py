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


def show_ir_spectrum(ir_modes):
    """Plot IR spectrum as stick plot."""
    if not ir_modes:
        return

    freqs = [m.frequency for m in ir_modes]
    intensities = [m.intensity for m in ir_modes]

    fig = go.Figure()
    for f, i in zip(freqs, intensities):
        fig.add_trace(go.Scatter(
            x=[f, f],
            y=[0, i],
            mode='lines',
            line=dict(color='red', width=1),
            showlegend=False
        ))

    fig.update_layout(
        title="IR Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Intensity (km/mol)",
        template="plotly_white"
    )

    st.plotly_chart(fig, use_container_width=True)


def show_orbital_energies(orbitals):
    """Plot orbital energy diagram."""
    if not orbitals:
        return

    occupied = [o for o in orbitals if o.occupation > 0]
    virtual = [o for o in orbitals if o.occupation == 0]

    # Only show last few occupied and first few virtual
    homo_idx = len(occupied) - 1 if occupied else -1
    lumo_idx = 0 if virtual else -1

    # Get energies around HOMO-LUMO gap
    show_occ = occupied[-10:] if len(occupied) > 10 else occupied
    show_virt = virtual[:10] if len(virtual) > 10 else virtual

    fig = go.Figure()

    # Occupied orbitals
    for o in show_occ:
        color = 'blue' if o.index != homo_idx else 'green'
        fig.add_trace(go.Scatter(
            x=[0.3, 0.7],
            y=[o.energy_ev, o.energy_ev],
            mode='lines',
            line=dict(color=color, width=3),
            name=f"{'HOMO' if o.index == homo_idx else 'Occ'}: {o.energy_ev:.2f} eV",
            showlegend=(o.index == homo_idx)
        ))

    # Virtual orbitals
    for o in show_virt:
        color = 'red' if o == virtual[0] else 'orange'
        fig.add_trace(go.Scatter(
            x=[0.3, 0.7],
            y=[o.energy_ev, o.energy_ev],
            mode='lines',
            line=dict(color=color, width=3),
            name=f"{'LUMO' if o == virtual[0] else 'Virt'}: {o.energy_ev:.2f} eV",
            showlegend=(o == virtual[0])
        ))

    fig.update_layout(
        title="Frontier Orbitals",
        yaxis_title="Energy (eV)",
        xaxis=dict(showticklabels=False, showgrid=False),
        template="plotly_white",
        showlegend=True
    )

    st.plotly_chart(fig, use_container_width=True)


def show_charges_comparison(mulliken: dict, loewdin: dict):
    """Compare Mulliken and Loewdin charges."""
    if not mulliken:
        return

    data = []
    for idx, (element, m_charge) in mulliken.items():
        l_charge = loewdin.get(idx, (element, 0))[1] if loewdin else 0
        data.append({
            'Atom': f"{element}{idx}",
            'Mulliken': m_charge,
            'Loewdin': l_charge
        })

    df = pd.DataFrame(data)

    fig = go.Figure()
    fig.add_trace(go.Bar(name='Mulliken', x=df['Atom'], y=df['Mulliken']))
    if loewdin:
        fig.add_trace(go.Bar(name='Loewdin', x=df['Atom'], y=df['Loewdin']))

    fig.update_layout(
        title="Atomic Charges",
        xaxis_title="Atom",
        yaxis_title="Charge",
        barmode='group',
        template="plotly_white"
    )

    st.plotly_chart(fig, use_container_width=True)


def show_nmr_shifts(nmr_data):
    """Display NMR chemical shifts."""
    if not nmr_data or not nmr_data.chemical_shifts:
        return

    # Group by element
    elements = {}
    for shift in nmr_data.chemical_shifts:
        if shift.element not in elements:
            elements[shift.element] = []
        elements[shift.element].append(shift)

    for elem, shifts in elements.items():
        st.write(f"**{elem} Chemical Shifts:**")
        data = [{
            'Atom': s.atom_index,
            'Isotropic (ppm)': f"{s.isotropic:.2f}",
            'Anisotropy (ppm)': f"{s.anisotropy:.2f}"
        } for s in shifts]
        st.dataframe(pd.DataFrame(data), use_container_width=True)


def out_preview_page():
    """Main Streamlit page for ORCA output file preview."""
    st.title("ORCA Output Viewer")

    # File upload
    uploaded_file = st.file_uploader("Upload ORCA output file", type=['out'])

    # Or use sample file
    sample_files = ['p1xs0.out', 'p1xs0p.out', 'p1xs0vg.out']
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
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Method", result.job_info.method or "N/A")
    with col2:
        st.metric("Basis Set", result.job_info.basis_set or "N/A")
    with col3:
        st.metric("Charge", result.job_info.charge)
    with col4:
        st.metric("Multiplicity", result.job_info.multiplicity)

    # Energy Summary
    st.header("Energy")
    col1, col2, col3 = st.columns(3)
    with col1:
        if result.final_energy:
            st.metric("Final Energy (Eh)", f"{result.final_energy:.6f}")
    with col2:
        if result.thermochemistry and result.thermochemistry.gibbs_free_energy:
            st.metric("Gibbs Free Energy (Eh)",
                     f"{result.thermochemistry.gibbs_free_energy:.6f}")
    with col3:
        if result.thermochemistry and result.thermochemistry.zero_point_energy:
            st.metric("ZPE (Eh)", f"{result.thermochemistry.zero_point_energy:.6f}")

    # Dipole and Polarizability
    col1, col2 = st.columns(2)
    with col1:
        if result.dipole_moment:
            st.subheader("Dipole Moment")
            st.write(f"**X:** {result.dipole_moment.x:.4f} a.u.")
            st.write(f"**Y:** {result.dipole_moment.y:.4f} a.u.")
            st.write(f"**Z:** {result.dipole_moment.z:.4f} a.u.")
            st.write(f"**Total:** {result.dipole_moment.magnitude_debye:.3f} Debye")

    with col2:
        if result.polarizability:
            st.subheader("Polarizability")
            st.write(f"**Isotropic:** {result.polarizability.isotropic:.2f} a.u.")
            if result.polarizability.tensor:
                st.write("**Tensor (a.u.):**")
                df = pd.DataFrame(result.polarizability.tensor,
                                  columns=['X', 'Y', 'Z'],
                                  index=['X', 'Y', 'Z'])
                st.dataframe(df.style.format("{:.2f}"))

    # Orbital Energies
    if result.orbital_energies:
        st.header("Orbital Energies")
        occupied = [o for o in result.orbital_energies if o.occupation > 0]
        virtual = [o for o in result.orbital_energies if o.occupation == 0]

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Occupied", len(occupied))
        with col2:
            st.metric("Virtual", len(virtual))
        with col3:
            if occupied and virtual:
                gap = virtual[0].energy_ev - occupied[-1].energy_ev
                st.metric("HOMO-LUMO Gap (eV)", f"{gap:.2f}")

        if occupied and virtual:
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**HOMO:** {occupied[-1].energy_ev:.3f} eV")
            with col2:
                st.write(f"**LUMO:** {virtual[0].energy_ev:.3f} eV")

        show_orbital_energies(result.orbital_energies)

    # Energy Convergence Plot
    if result.optimization_energies and len(result.optimization_energies) > 1:
        st.header("Optimization Convergence")
        show_energy_plot(result.optimization_energies, "Energy during Optimization")

    # Vibrational Frequencies and IR Spectrum
    if result.frequencies or result.ir_spectrum:
        st.header("Vibrational Analysis")

        if result.frequencies:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Number of Modes", len(result.frequencies))
            with col2:
                st.metric("Lowest (cm⁻¹)", f"{min(result.frequencies):.1f}")
            with col3:
                st.metric("Highest (cm⁻¹)", f"{max(result.frequencies):.1f}")

        if result.ir_spectrum:
            show_ir_spectrum(result.ir_spectrum)

            # Show top IR peaks
            sorted_modes = sorted(result.ir_spectrum, key=lambda x: x.intensity, reverse=True)
            st.write("**Top IR Peaks:**")
            for i, mode in enumerate(sorted_modes[:5]):
                st.write(f"{i+1}. {mode.frequency:.1f} cm⁻¹ ({mode.intensity:.1f} km/mol)")

        # Frequencies table
        if result.frequencies:
            with st.expander("View all frequencies"):
                freq_df = pd.DataFrame({
                    'Mode': list(range(1, len(result.frequencies) + 1)),
                    'Frequency (cm⁻¹)': result.frequencies
                })
                st.dataframe(freq_df, use_container_width=True)

    # Atomic Charges
    if result.mulliken_charges:
        st.header("Atomic Charges")
        show_charges_comparison(result.mulliken_charges, result.loewdin_charges)

    # Mayer Bond Orders
    if result.mayer_bond_orders:
        st.header("Mayer Bond Orders")
        st.write(f"**Total bonds:** {len(result.mayer_bond_orders)}")

        # Show top bond orders
        sorted_bonds = sorted(result.mayer_bond_orders, key=lambda x: x[2], reverse=True)
        with st.expander("View bond orders"):
            data = [{
                'Atom 1': b[0],
                'Atom 2': b[1],
                'Order': f"{b[2]:.3f}"
            } for b in sorted_bonds[:20]]
            st.dataframe(pd.DataFrame(data), use_container_width=True)

    # NMR Data
    if result.nmr_data:
        st.header("NMR Data")
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Chemical Shifts", len(result.nmr_data.chemical_shifts))
        with col2:
            st.metric("J-Couplings", len(result.nmr_data.j_couplings))

        show_nmr_shifts(result.nmr_data)

        if result.nmr_data.j_couplings:
            st.write("**J-Couplings (Hz):**")
            data = [{
                'Atom 1': j.atom1_index,
                'Atom 2': j.atom2_index,
                'J (Hz)': f"{j.value:.2f}"
            } for j in result.nmr_data.j_couplings]
            st.dataframe(pd.DataFrame(data), use_container_width=True)

    # Thermochemistry
    if result.thermochemistry:
        st.header("Thermochemistry")
        thermo = result.thermochemistry

        col1, col2 = st.columns(2)
        with col1:
            st.write(f"**Temperature:** {thermo.temperature} K")
            st.write(f"**Pressure:** {thermo.pressure} atm")
            st.write(f"**Total Mass:** {thermo.total_mass:.2f} AMU")

        with col2:
            if thermo.zero_point_energy:
                st.write(f"**ZPE:** {thermo.zero_point_energy:.6f} Eh")
            if thermo.enthalpy:
                st.write(f"**Enthalpy:** {thermo.enthalpy:.6f} Eh")
            if thermo.gibbs_free_energy:
                st.write(f"**Gibbs Free Energy:** {thermo.gibbs_free_energy:.6f} Eh")

        # Entropy breakdown
        if any([thermo.entropy_electronic, thermo.entropy_vibrational,
                thermo.entropy_rotational, thermo.entropy_translational]):
            with st.expander("Entropy breakdown"):
                if thermo.entropy_electronic is not None:
                    st.write(f"Electronic: {thermo.entropy_electronic:.6f} Eh")
                if thermo.entropy_vibrational is not None:
                    st.write(f"Vibrational: {thermo.entropy_vibrational:.6f} Eh")
                if thermo.entropy_rotational is not None:
                    st.write(f"Rotational: {thermo.entropy_rotational:.6f} Eh")
                if thermo.entropy_translational is not None:
                    st.write(f"Translational: {thermo.entropy_translational:.6f} Eh")
                if thermo.entropy is not None:
                    st.write(f"**Total:** {thermo.entropy:.6f} Eh")

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
