"""
Spectra Visualization Component - Complete IR/Raman/UV-Vis

Features:
- Multi-molecule comparison with state labels
- IR spectrum with frequency range, smoothing, normalization
- Raman spectrum with activity
- UV-Vis from TDDFT with broadening
- Select All / Clear All buttons
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter1d
from typing import List, Optional


def render_spectra_viz(df: pd.DataFrame):
    """Render spectra visualization with molecule+state selector."""
    
    st.subheader("📈 Spectroscopy")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build molecule options with state labels
    mol_options = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        if state and str(state) != "nan":
            label = f"{mol_id} [{state}]"
        else:
            label = mol_id
        mol_options.append({"label": label, "idx": idx})
    
    # Unique labels
    unique_labels = list(dict.fromkeys([m["label"] for m in mol_options]))
    
    if not unique_labels:
        st.warning("No molecules available")
        return
    
    # Molecule selector - default top 10
    if "spectra_selected_mols" not in st.session_state:
        st.session_state.spectra_selected_mols = unique_labels[:min(10, len(unique_labels))]
    
    # Validate against current options
    valid_selection = [m for m in st.session_state.spectra_selected_mols if m in unique_labels]
    
    selected_mols = st.multiselect(
        "Select Molecules to Compare",
        unique_labels,
        default=valid_selection if valid_selection else unique_labels[:min(10, len(unique_labels))],
        key="spectra_mol_select"
    )
    st.session_state.spectra_selected_mols = selected_mols
    
    if not selected_mols:
        st.warning("Select at least one molecule")
        return
    
    # Get selected data
    selected_indices = [m["idx"] for m in mol_options if m["label"] in selected_mols]
    selected_df = df.loc[selected_indices]
    
    st.caption(f"Comparing {len(selected_mols)} spectra")
    
    # Spectra tabs
    spec_tabs = st.tabs(["🔴 IR Spectrum", "🟢 Raman Spectrum", "🟣 UV-Vis Absorption"])
    
    with spec_tabs[0]:
        render_ir_spectrum(selected_df, selected_mols)
    with spec_tabs[1]:
        render_raman_spectrum(selected_df, selected_mols)
    with spec_tabs[2]:
        render_uvvis_spectrum(selected_df, selected_mols)


def smooth_spectrum(x: np.ndarray, y: np.ndarray, sigma: float) -> tuple:
    """Apply Gaussian smoothing to spectrum."""
    if len(x) == 0 or len(y) == 0:
        return x, y
    if sigma <= 0:
        return x, y
    
    # Sort by x
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    y_sorted = y[sort_idx]
    
    # Interpolate to regular grid
    x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), 1000)
    y_interp = np.interp(x_smooth, x_sorted, y_sorted)
    
    # Apply Gaussian filter
    y_smooth = gaussian_filter1d(y_interp, sigma=sigma)
    
    return x_smooth, y_smooth


def render_ir_spectrum(df: pd.DataFrame, labels: List[str]):
    """Render IR spectrum comparison."""
    
    # Settings
    with st.expander("⚙️ IR Settings", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1:
            freq_range = st.slider("Frequency Range (cm⁻¹)", 0, 4000, (400, 4000), key="ir_freq_range")
        with col2:
            smoothing = st.slider("Gaussian Smoothing", 0, 50, 10, key="ir_smoothing")
        with col3:
            normalize = st.checkbox("Normalize Intensity", False, key="ir_normalize")
    
    # Collect IR data
    ir_data = {}
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        ir = row.get("ir")
        if ir is not None and hasattr(ir, 'empty') and not ir.empty:
            ir_data[label] = ir
    
    if not ir_data:
        st.warning("No IR data available for selected molecules")
        return
    
    # Create plot
    fig = go.Figure()
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692']
    
    for i, (label, ir) in enumerate(ir_data.items()):
        # Get frequency and intensity columns
        if "freq_cm-1" in ir.columns:
            x = ir["freq_cm-1"].values
        elif "frequency" in ir.columns:
            x = ir["frequency"].values
        else:
            x = ir.iloc[:, 0].values
        
        if "intensity_km/mol" in ir.columns:
            y = ir["intensity_km/mol"].values
        elif "intensity" in ir.columns:
            y = ir["intensity"].values
        else:
            y = ir.iloc[:, 1].values
        
        # Filter by frequency range
        mask = (x >= freq_range[0]) & (x <= freq_range[1])
        x, y = x[mask], y[mask]
        
        if len(x) == 0:
            continue
        
        # Apply smoothing
        if smoothing > 0:
            x, y = smooth_spectrum(x, y, smoothing)
        
        # Normalize
        if normalize and len(y) > 0 and np.max(np.abs(y)) > 0:
            y = y / np.max(np.abs(y))
        
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(
            x=x, y=y,
            mode='lines',
            name=label,
            line=dict(color=color, width=1.5),
            hovertemplate='%{x:.0f} cm⁻¹<br>%{y:.2f}<extra>%{fullData.name}</extra>'
        ))
    
    fig.update_layout(
        title="IR Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Intensity (km/mol)" if not normalize else "Normalized Intensity",
        xaxis=dict(autorange="reversed"),  # IR convention
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_raman_spectrum(df: pd.DataFrame, labels: List[str]):
    """Render Raman spectrum comparison."""
    
    with st.expander("⚙️ Raman Settings", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1:
            freq_range = st.slider("Frequency Range (cm⁻¹)", 0, 4000, (400, 4000), key="raman_freq_range")
        with col2:
            smoothing = st.slider("Gaussian Smoothing", 0, 50, 10, key="raman_smoothing")
        with col3:
            normalize = st.checkbox("Normalize Activity", False, key="raman_normalize")
    
    # Collect Raman data
    raman_data = {}
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        raman = row.get("raman")
        if raman is not None and hasattr(raman, 'empty') and not raman.empty:
            raman_data[label] = raman
    
    if not raman_data:
        st.warning("No Raman data available for selected molecules")
        return
    
    fig = go.Figure()
    colors = ['#00CC96', '#636EFA', '#EF553B', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692']
    
    for i, (label, raman) in enumerate(raman_data.items()):
        if "freq_cm-1" in raman.columns:
            x = raman["freq_cm-1"].values
        else:
            x = raman.iloc[:, 0].values
        
        if "activity" in raman.columns:
            y = raman["activity"].values
        else:
            y = raman.iloc[:, 1].values
        
        mask = (x >= freq_range[0]) & (x <= freq_range[1])
        x, y = x[mask], y[mask]
        
        if len(x) == 0:
            continue
        
        if smoothing > 0:
            x, y = smooth_spectrum(x, y, smoothing)
        
        if normalize and len(y) > 0 and np.max(np.abs(y)) > 0:
            y = y / np.max(np.abs(y))
        
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name=label,
                                line=dict(color=color, width=1.5)))
    
    fig.update_layout(
        title="Raman Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Activity" if not normalize else "Normalized Activity",
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_uvvis_spectrum(df: pd.DataFrame, labels: List[str]):
    """Render UV-Vis spectrum from TDDFT data."""
    
    with st.expander("⚙️ UV-Vis Settings", expanded=False):
        col1, col2, col3 = st.columns(3)
        with col1:
            wl_range = st.slider("Wavelength Range (nm)", 100, 900, (200, 700), key="uvvis_wl_range")
        with col2:
            sigma = st.slider("Gaussian Broadening (nm)", 5, 50, 20, key="uvvis_sigma")
        with col3:
            show_sticks = st.checkbox("Show Stick Spectrum", True, key="uvvis_sticks")
    
    # Collect TDDFT data
    tddft_data = {}
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state and str(state) != "nan" else mol_id
        
        tddft = row.get("tddft_states")
        if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
            if "energy_ev" in tddft.columns:
                tddft_data[label] = tddft
    
    if not tddft_data:
        st.warning("No TDDFT data available for selected molecules")
        return
    
    fig = go.Figure()
    colors = ['#AB63FA', '#636EFA', '#EF553B', '#00CC96', '#FFA15A', '#19D3F3', '#FF6692']
    
    for i, (label, tddft) in enumerate(tddft_data.items()):
        # Convert energy (eV) to wavelength (nm)
        energy_ev = tddft["energy_ev"].values
        wavelength = 1239.84 / energy_ev  # nm from eV
        
        # Get oscillator strength
        if "weight" in tddft.columns:
            fosc = np.abs(tddft["weight"].values)
        elif "fosc" in tddft.columns:
            fosc = tddft["fosc"].values
        else:
            fosc = np.ones_like(wavelength) * 0.5
        
        # Filter by wavelength range
        mask = (wavelength >= wl_range[0]) & (wavelength <= wl_range[1])
        wl_filtered = wavelength[mask]
        fosc_filtered = fosc[mask]
        
        color = colors[i % len(colors)]
        
        # Add stick spectrum
        if show_sticks:
            for wl, f in zip(wl_filtered, fosc_filtered):
                fig.add_trace(go.Scatter(
                    x=[wl, wl], y=[0, f],
                    mode='lines',
                    line=dict(color=color, width=2),
                    showlegend=False,
                    hovertemplate=f'{label}<br>λ: {wl:.1f} nm<br>f: {f:.3f}<extra></extra>'
                ))
        
        # Create broadened spectrum
        wl_grid = np.linspace(wl_range[0], wl_range[1], 500)
        spectrum = np.zeros_like(wl_grid)
        
        for wl, f in zip(wl_filtered, fosc_filtered):
            spectrum += f * np.exp(-0.5 * ((wl_grid - wl) / sigma) ** 2)
        
        fig.add_trace(go.Scatter(
            x=wl_grid, y=spectrum,
            mode='lines',
            name=label,
            line=dict(color=color, width=1.5)
        ))
    
    fig.update_layout(
        title="UV-Vis Absorption Spectrum",
        xaxis_title="Wavelength (nm)",
        yaxis_title="Intensity / Oscillator Strength",
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    st.plotly_chart(fig, use_container_width=True)
