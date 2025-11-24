"""
Unified file preview with automatic type detection.
"""

import streamlit as st
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def detect_file_type(filename: str) -> str:
    """Detect ORCA file type from filename."""
    name = filename.lower()

    if name.endswith('.xyz'):
        if '_trj' in name:
            return 'trajectory'
        return 'xyz'
    elif name.endswith('.out'):
        return 'out'
    elif name.endswith('.hess'):
        return 'hess'
    elif name.endswith('.spectrum'):
        return 'spectrum'
    elif name.endswith('.engrad'):
        return 'engrad'
    elif name.endswith('.property.txt'):
        return 'property'
    elif name.endswith('.inp'):
        return 'inp'
    elif name.endswith('.cpcm_corr'):
        return 'cpcm_corr'
    elif name.endswith('.cpcm'):
        return 'cpcm'
    elif name.endswith('.bibtex') or name.endswith('.bib'):
        return 'bibtex'
    else:
        return 'unknown'


def unified_preview_page():
    """Unified file upload with automatic type detection."""
    st.title("ORCA File Analyzer")
    st.markdown("Upload any ORCA file - type will be automatically detected")

    uploaded_file = st.file_uploader(
        "Upload ORCA file",
        type=['xyz', 'out', 'hess', 'spectrum', 'engrad', 'txt', 'inp', 'cpcm', 'cpcm_corr', 'bibtex', 'bib']
    )

    if not uploaded_file:
        st.info("Supported formats: .xyz, .out, .hess, .spectrum, .engrad, .property.txt, .inp, .cpcm, .cpcm_corr, .bibtex")
        return

    file_type = detect_file_type(uploaded_file.name)
    content = uploaded_file.read().decode('utf-8')

    st.success(f"Detected file type: **{file_type.upper()}**")

    if file_type == 'xyz' or file_type == 'trajectory':
        from parsers.xyz_parser import parse_xyz_content, parse_trajectory_xyz
        if file_type == 'trajectory':
            # For trajectory, we need to handle differently
            st.info("Trajectory file detected")
        geom = parse_xyz_content(content)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Atoms", geom.num_atoms)
        with col2:
            if geom.energy:
                st.metric("Energy (Eh)", f"{geom.energy:.6f}")

        # Show coordinates table
        import pandas as pd
        data = [{'Element': a.symbol, 'X': f"{a.x:.4f}", 'Y': f"{a.y:.4f}", 'Z': f"{a.z:.4f}"}
                for a in geom.atoms]
        st.dataframe(pd.DataFrame(data), use_container_width=True)

    elif file_type == 'out':
        from parsers.out_parser import parse_out_content
        result = parse_out_content(content)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Final Energy (Eh)", f"{result.final_energy:.6f}" if result.final_energy else "N/A")
        with col2:
            if result.dipole_moment:
                st.metric("Dipole (Debye)", f"{result.dipole_moment.magnitude_debye:.3f}")

        if result.frequencies:
            st.metric("Vibrational Modes", len(result.frequencies))

    elif file_type == 'hess':
        from parsers.hess_parser import parse_hess_content, find_strongest_peaks
        hess = parse_hess_content(content)

        real_freqs = [f for f in hess.frequencies if f > 0]
        st.metric("Vibrational Modes", len(real_freqs))

        peaks = find_strongest_peaks(hess, 5)
        if peaks:
            st.write("**Top IR Peaks:**")
            for i, p in enumerate(peaks):
                st.write(f"{i+1}. {p.frequency:.1f} cm⁻¹ ({p.intensity:.1f} km/mol)")

    elif file_type == 'spectrum':
        from parsers.spectrum_parser import parse_spectrum_content, find_peaks
        spectrum = parse_spectrum_content(content)

        st.metric("Data Points", len(spectrum.energy))
        peaks = find_peaks(spectrum)
        if peaks:
            st.write(f"**Strongest peak:** {peaks[0]['energy_nm']:.1f} nm")

    elif file_type == 'engrad':
        from parsers.engrad_parser import parse_engrad_content
        data = parse_engrad_content(content)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Energy (Eh)", f"{data.energy:.6f}")
        with col2:
            st.metric("Max Gradient", f"{data.max_gradient:.2e}")

    elif file_type == 'property':
        from parsers.property_parser import parse_property_content
        data = parse_property_content(content)

        st.metric("Status", data.status)
        st.metric("Atoms", data.num_atoms)

    elif file_type == 'inp':
        from parsers.inp_parser import parse_inp_content
        data = parse_inp_content(content)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Method", data.method or "N/A")
        with col2:
            st.metric("Basis Set", data.basis_set or "N/A")
        st.write(f"**Charge:** {data.charge}, **Multiplicity:** {data.multiplicity}")

    elif file_type == 'cpcm_corr':
        from parsers.cpcm_corr_parser import parse_cpcm_corr_content
        data = parse_cpcm_corr_content(content)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Corrected Energy (Eh)", f"{data.corrected_energy:.6f}")
        with col2:
            st.metric("Surface Charges", len(data.corrected_charges))

    elif file_type == 'cpcm':
        from parsers.cpcm_parser import parse_cpcm_content
        data = parse_cpcm_content(content)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Dielectric Energy (Eh)", f"{data.dielectric_energy:.6f}")
        with col2:
            st.metric("Volume (bohr³)", f"{data.volume:.1f}")

    elif file_type == 'bibtex':
        from parsers.bibtex_parser import parse_bibtex_content
        data = parse_bibtex_content(content)

        st.metric("Citations", len(data.citations))
        for c in data.citations[:5]:
            st.write(f"- {c.author.split(',')[0]} et al. ({c.year})")
    else:
        st.error("Unknown file type")


if __name__ == '__main__':
    unified_preview_page()
