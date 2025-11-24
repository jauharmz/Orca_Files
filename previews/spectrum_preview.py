"""
Spectrum file preview using Streamlit and Plotly.

This module provides visualization for ORCA spectrum data.
"""

import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.spectrum_parser import parse_spectrum_file, parse_spectrum_content, find_peaks


def show_spectrum_plot(spectrum, x_unit: str = 'nm'):
    """
    Plot spectrum with selected x-axis unit.

    Args:
        spectrum: SpectrumData object
        x_unit: 'nm' for wavelength or 'cm' for wavenumber
    """
    fig = go.Figure()

    if x_unit == 'nm':
        x_data = spectrum.energy_nm
        x_label = "Wavelength (nm)"
    else:
        x_data = spectrum.energy
        x_label = "Energy (cm⁻¹)"

    # Total spectrum
    fig.add_trace(go.Scatter(
        x=x_data,
        y=spectrum.total_spectrum,
        mode='lines',
        name='Total Spectrum',
        line=dict(color='blue')
    ))

    # FC component
    if spectrum.intensity_fc:
        fig.add_trace(go.Scatter(
            x=x_data,
            y=spectrum.intensity_fc,
            mode='lines',
            name='Franck-Condon',
            line=dict(color='green', dash='dash')
        ))

    # HT component
    if spectrum.intensity_ht:
        fig.add_trace(go.Scatter(
            x=x_data,
            y=spectrum.intensity_ht,
            mode='lines',
            name='Herzberg-Teller',
            line=dict(color='red', dash='dot')
        ))

    fig.update_layout(
        title="Vibronic Spectrum",
        xaxis_title=x_label,
        yaxis_title="Intensity",
        template="plotly_white",
        hovermode='x unified'
    )

    st.plotly_chart(fig, use_container_width=True)


def spectrum_preview_page():
    """Main Streamlit page for spectrum file preview."""
    st.title("Spectrum Viewer")

    # File upload
    uploaded_file = st.file_uploader("Upload spectrum file", type=['spectrum'])

    # Or use sample file
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0vg.spectrum'
        )
        spectrum = parse_spectrum_file(filepath)

    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        spectrum = parse_spectrum_content(content)
    else:
        st.info("Please upload a spectrum file or select a sample file.")
        return

    # Display info
    st.header("Spectrum Information")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Data Points", len(spectrum.energy))
    with col2:
        st.metric("Energy Range (cm⁻¹)",
                 f"{min(spectrum.energy):.0f} - {max(spectrum.energy):.0f}")
    with col3:
        st.metric("Wavelength Range (nm)",
                 f"{min(spectrum.energy_nm):.0f} - {max(spectrum.energy_nm):.0f}")

    # Plot options
    st.header("Spectrum Plot")
    x_unit = st.radio("X-axis unit:", ['nm', 'cm⁻¹'], horizontal=True)
    x_unit_key = 'nm' if x_unit == 'nm' else 'cm'

    show_spectrum_plot(spectrum, x_unit_key)

    # Peak analysis
    st.header("Peak Analysis")
    threshold = st.slider("Peak threshold (%)", 1, 50, 10) / 100
    peaks = find_peaks(spectrum, threshold)

    if peaks:
        st.write(f"Found {len(peaks)} peaks above {threshold*100:.0f}% threshold")

        # Show top peaks
        peak_data = []
        for i, peak in enumerate(peaks[:20]):
            peak_data.append({
                'Rank': i + 1,
                'Wavelength (nm)': f"{peak['energy_nm']:.1f}",
                'Energy (cm⁻¹)': f"{peak['energy_cm']:.1f}",
                'Intensity': f"{peak['intensity']:.2e}"
            })

        df = pd.DataFrame(peak_data)
        st.dataframe(df, use_container_width=True)
    else:
        st.info("No peaks found above threshold")

    # Export option
    if st.button("Export as JSON"):
        import json
        json_data = json.dumps(spectrum.to_dict(), indent=2)
        st.download_button(
            "Download JSON",
            json_data,
            "spectrum.json",
            "application/json"
        )


if __name__ == '__main__':
    spectrum_preview_page()
