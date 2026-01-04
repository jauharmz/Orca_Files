"""
Spectra Visualization Component - IR/Raman/UV-Vis

Features:
- Molecule + State filter
- Multi-molecule comparison
- Smoothing and normalization
- Interactive Plotly charts
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter1d
from typing import List


def render_spectra_viz(df: pd.DataFrame):
    """Render spectra visualization with molecule+state filter."""
    
    st.subheader("📈 Spectroscopy")
    
    if df.empty:
        st.warning("No data available")
        return
    
    # Build unique labels
    labels = []
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        labels.append({"label": label, "idx": idx})
    
    unique_labels = list(dict.fromkeys([l["label"] for l in labels]))
    
    # Dual filter
    c1, c2, c3 = st.columns([4, 1, 1])
    with c1:
        selected = st.multiselect("Select Molecules", unique_labels, 
                                 default=unique_labels[:min(3, len(unique_labels))], key="spectra_sel")
    with c2:
        if st.button("✅ All", key="spectra_all"):
            st.session_state.spectra_sel = unique_labels
            st.rerun()
    with c3:
        if st.button("❌ None", key="spectra_none"):
            st.session_state.spectra_sel = []
            st.rerun()
    
    if not selected:
        st.warning("Select at least one molecule")
        return
    
    # Get selected indices
    sel_indices = [l["idx"] for l in labels if l["label"] in selected]
    sel_df = df.loc[sel_indices]
    
    # Spectra tabs
    tabs = st.tabs(["🔴 IR", "🟢 Raman", "🟣 UV-Vis"])
    
    with tabs[0]:
        render_ir(sel_df, selected)
    with tabs[1]:
        render_raman(sel_df, selected)
    with tabs[2]:
        render_uvvis(sel_df, selected)


def smooth_spectrum(x, y, sigma):
    """Apply Gaussian smoothing."""
    if sigma <= 0:
        return x, y
    x_smooth = np.linspace(x.min(), x.max(), 1000)
    y_interp = np.interp(x_smooth, np.sort(x), y[np.argsort(x)])
    return x_smooth, gaussian_filter1d(y_interp, sigma=sigma)


def render_ir(df: pd.DataFrame, labels: List[str]):
    """Render IR spectra."""
    
    # Settings
    with st.expander("⚙️ Settings", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            freq_range = st.slider("Frequency (cm⁻¹)", 0, 4000, (400, 4000), key="ir_freq")
        with c2:
            smoothing = st.slider("Smoothing", 0, 50, 0, key="ir_smooth")
        with c3:
            normalize = st.checkbox("Normalize", False, key="ir_norm")
    
    # Collect IR data
    ir_data = {}
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        ir = row.get("ir")
        if ir is not None and hasattr(ir, 'empty') and not ir.empty:
            ir_data[label] = ir
    
    if not ir_data:
        st.warning("No IR data available")
        return
    
    # Plot
    fig = go.Figure()
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A']
    
    for i, (label, ir) in enumerate(ir_data.items()):
        x = ir["freq_cm-1"].values
        y = ir["intensity_km/mol"].values if "intensity_km/mol" in ir.columns else ir.iloc[:, 1].values
        
        mask = (x >= freq_range[0]) & (x <= freq_range[1])
        x, y = x[mask], y[mask]
        
        if len(x) == 0:
            continue
        
        if smoothing > 0:
            x, y = smooth_spectrum(x, y, smoothing)
        
        if normalize and len(y) > 0 and y.max() > 0:
            y = y / y.max()
        
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=label,
                                line=dict(color=colors[i % len(colors)], width=1.5)))
    
    fig.update_layout(
        title="IR Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Intensity" if not normalize else "Normalized",
        xaxis=dict(autorange="reversed"),
        hovermode="x unified"
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_raman(df: pd.DataFrame, labels: List[str]):
    """Render Raman spectra."""
    
    with st.expander("⚙️ Settings", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            freq_range = st.slider("Frequency", 0, 4000, (400, 4000), key="raman_freq")
        with c2:
            smoothing = st.slider("Smoothing", 0, 50, 0, key="raman_smooth")
        with c3:
            normalize = st.checkbox("Normalize", False, key="raman_norm")
    
    raman_data = {}
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        raman = row.get("raman")
        if raman is not None and hasattr(raman, 'empty') and not raman.empty:
            raman_data[label] = raman
    
    if not raman_data:
        st.warning("No Raman data available")
        return
    
    fig = go.Figure()
    colors = ['#00CC96', '#636EFA', '#EF553B', '#AB63FA', '#FFA15A']
    
    for i, (label, raman) in enumerate(raman_data.items()):
        x = raman["freq_cm-1"].values
        y = raman["activity"].values if "activity" in raman.columns else raman.iloc[:, 1].values
        
        mask = (x >= freq_range[0]) & (x <= freq_range[1])
        x, y = x[mask], y[mask]
        
        if len(x) == 0:
            continue
        
        if smoothing > 0:
            x, y = smooth_spectrum(x, y, smoothing)
        
        if normalize and len(y) > 0 and y.max() > 0:
            y = y / y.max()
        
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=label,
                                line=dict(color=colors[i % len(colors)], width=1.5)))
    
    fig.update_layout(
        title="Raman Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Activity" if not normalize else "Normalized",
        hovermode="x unified"
    )
    
    st.plotly_chart(fig, use_container_width=True)


def render_uvvis(df: pd.DataFrame, labels: List[str]):
    """Render UV-Vis spectra from TDDFT."""
    
    with st.expander("⚙️ Settings", expanded=False):
        c1, c2, c3 = st.columns(3)
        with c1:
            wl_range = st.slider("Wavelength (nm)", 100, 900, (200, 700), key="uv_wl")
        with c2:
            sigma = st.slider("Broadening σ", 5, 50, 20, key="uv_sigma")
        with c3:
            show_sticks = st.checkbox("Show Sticks", True, key="uv_sticks")
    
    tddft_data = {}
    for idx, row in df.iterrows():
        mol_id = row.get("molecule_id", "unknown")
        state = row.get("optimized_state", "")
        label = f"{mol_id} [{state}]" if state else mol_id
        tddft = row.get("tddft_states")
        if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
            if "energy_ev" in tddft.columns:
                tddft_data[label] = tddft
    
    if not tddft_data:
        st.warning("No TDDFT data available")
        return
    
    fig = go.Figure()
    colors = ['#AB63FA', '#636EFA', '#EF553B', '#00CC96', '#FFA15A']
    
    for i, (label, tddft) in enumerate(tddft_data.items()):
        wavelength = 1239.84 / tddft["energy_ev"]
        fosc = tddft.get("weight", pd.Series([0.5] * len(tddft))).abs()
        
        mask = (wavelength >= wl_range[0]) & (wavelength <= wl_range[1])
        wl_f = wavelength[mask].values
        fosc_f = fosc[mask].values
        
        if show_sticks:
            for wl, f in zip(wl_f, fosc_f):
                fig.add_trace(go.Scatter(x=[wl, wl], y=[0, f], mode="lines",
                                        line=dict(color=colors[i % len(colors)], width=2),
                                        showlegend=False))
        
        # Broadened spectrum
        wl_grid = np.linspace(wl_range[0], wl_range[1], 500)
        spectrum = np.zeros_like(wl_grid)
        for wl, f in zip(wl_f, fosc_f):
            spectrum += f * np.exp(-0.5 * ((wl_grid - wl) / sigma) ** 2)
        
        fig.add_trace(go.Scatter(x=wl_grid, y=spectrum, mode="lines", name=label,
                                line=dict(color=colors[i % len(colors)], width=1.5)))
    
    fig.update_layout(
        title="UV-Vis Absorption",
        xaxis_title="Wavelength (nm)",
        yaxis_title="Intensity",
        hovermode="x unified"
    )
    
    st.plotly_chart(fig, use_container_width=True)
