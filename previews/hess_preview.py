"""
Hessian file preview using Streamlit and Plotly.

This module provides visualization for IR spectra from Hessian files.
"""

import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.hess_parser import (
    parse_hess_file,
    parse_hess_content,
    get_ir_spectrum_data,
    find_strongest_peaks
)


def show_ir_spectrum_plot(frequencies: list[float], intensities: list[float]):
    """Plot IR spectrum as stick plot."""
    fig = go.Figure()

    # Stick spectrum
    for freq, intensity in zip(frequencies, intensities):
        fig.add_trace(go.Scatter(
            x=[freq, freq],
            y=[0, intensity],
            mode='lines',
            line=dict(color='blue', width=1),
            showlegend=False,
            hovertemplate=f'{freq:.1f} cm⁻¹<br>{intensity:.1f} km/mol<extra></extra>'
        ))

    fig.update_layout(
        title="IR Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Intensity (km/mol)",
        template="plotly_white",
        xaxis=dict(range=[0, max(frequencies) * 1.1] if frequencies else [0, 4000])
    )

    st.plotly_chart(fig, use_container_width=True)


def hess_preview_page():
    """Main Streamlit page for Hessian file preview."""
    st.title("Hessian / IR Spectrum Viewer")

    # File upload
    uploaded_file = st.file_uploader("Upload Hessian file", type=['hess'])

    # Or use sample file
    use_sample = st.checkbox("Use sample file")

    if use_sample:
        sample_files = ['p1xs0.hess', 'p1xs0p_job2.hess']
        selected = st.selectbox("Select sample", sample_files)
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            selected
        )
        hess = parse_hess_file(filepath)

    elif uploaded_file:
        content = uploaded_file.read().decode('utf-8')
        hess = parse_hess_content(content)
    else:
        st.info("Please upload a Hessian file or select a sample file.")
        return

    # Display info
    st.header("Vibrational Analysis")
    real_freqs = [f for f in hess.frequencies if f > 0]

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Modes", hess.num_modes)
    with col2:
        st.metric("Real Modes", len(real_freqs))
    with col3:
        if real_freqs:
            st.metric("Frequency Range",
                     f"{min(real_freqs):.0f} - {max(real_freqs):.0f}")

    # IR Spectrum Plot
    if hess.ir_spectrum:
        st.header("IR Spectrum")
        frequencies, intensities = get_ir_spectrum_data(hess)
        show_ir_spectrum_plot(frequencies, intensities)

        # Strongest peaks table
        st.header("Strongest IR Peaks")
        num_peaks = st.slider("Number of peaks", 5, 20, 10)
        peaks = find_strongest_peaks(hess, num_peaks)

        peak_data = []
        for i, peak in enumerate(peaks):
            peak_data.append({
                'Rank': i + 1,
                'Frequency (cm⁻¹)': f"{peak.frequency:.1f}",
                'Intensity (km/mol)': f"{peak.intensity:.1f}",
            })

        df = pd.DataFrame(peak_data)
        st.dataframe(df, use_container_width=True)

    # All frequencies table
    with st.expander("View all frequencies"):
        freq_data = []
        for i, freq in enumerate(hess.frequencies):
            if freq > 0:
                # Find matching IR data
                intensity = 0
                for mode in hess.ir_spectrum:
                    if abs(mode.frequency - freq) < 0.1:
                        intensity = mode.intensity
                        break
                freq_data.append({
                    'Mode': i + 1,
                    'Frequency (cm⁻¹)': f"{freq:.1f}",
                    'Intensity (km/mol)': f"{intensity:.1f}"
                })

        df = pd.DataFrame(freq_data)
        st.dataframe(df, use_container_width=True)

    # Export option
    if st.button("Export as JSON"):
        import json
        json_data = json.dumps(hess.to_dict(), indent=2)
        st.download_button(
            "Download JSON",
            json_data,
            "hessian.json",
            "application/json"
        )


if __name__ == '__main__':
    hess_preview_page()
