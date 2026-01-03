"""
Spectra Visualization Component - IR/Raman/UV-Vis with customization

Features:
- Interactive Plotly charts
- Frequency range sliders
- Peak width/broadening controls
- Multi-molecule overlay
- Normalization options
- Data table with export
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from typing import Dict, List, Optional


def render_spectra_viz(df: pd.DataFrame):
    """Render spectra visualization tab."""
    
    st.header("📈 Spectroscopy")
    
    # Molecule selector (multi-select for comparison)
    mol_ids = df["molecule_id"].dropna().unique().tolist()
    selected_mols = st.multiselect(
        "Select Molecules to Compare",
        mol_ids,
        default=mol_ids[:min(3, len(mol_ids))],
        key="spectra_mol_select"
    )
    
    if not selected_mols:
        st.warning("Select at least one molecule")
        return
    
    # Spectra type tabs
    tabs = st.tabs(["🔴 IR", "🟢 Raman", "🟣 UV-Vis", "📊 IR-Raman Comparison"])
    
    with tabs[0]:
        render_ir_spectra(df, selected_mols)
    
    with tabs[1]:
        render_raman_spectra(df, selected_mols)
    
    with tabs[2]:
        render_uvvis_spectra(df, selected_mols)
    
    with tabs[3]:
        render_ir_raman_comparison(df, selected_mols)


def render_ir_spectra(df: pd.DataFrame, selected_mols: List[str]):
    """Render IR spectra with customization."""
    
    st.subheader("🔴 IR Spectrum")
    
    # Customization options
    with st.expander("⚙️ Customization", expanded=True):
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            freq_range = st.slider(
                "Frequency Range (cm⁻¹)",
                0, 4000, (400, 4000),
                key="ir_freq_range"
            )
        with col2:
            normalize = st.checkbox("Normalize", False, key="ir_normalize")
        with col3:
            invert_x = st.checkbox("Invert X-axis", True, key="ir_invert")
        with col4:
            show_peaks = st.checkbox("Show Peak Labels", False, key="ir_peaks")
    
    # Collect IR data
    ir_data = {}
    for mol_id in selected_mols:
        mol_row = df[df["molecule_id"] == mol_id]
        if not mol_row.empty:
            ir = mol_row.iloc[0].get("ir")
            if ir is not None and hasattr(ir, 'empty') and not ir.empty:
                ir_data[mol_id] = ir
    
    if not ir_data:
        st.warning("No IR data available for selected molecules")
        return
    
    # Create plot
    fig = go.Figure()
    
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', 
              '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    
    for i, (mol_id, ir) in enumerate(ir_data.items()):
        freq_col = "freq_cm-1"
        int_col = "intensity_km/mol" if "intensity_km/mol" in ir.columns else ir.columns[1]
        
        x = ir[freq_col]
        y = ir[int_col]
        
        # Apply frequency filter
        mask = (x >= freq_range[0]) & (x <= freq_range[1])
        x, y = x[mask], y[mask]
        
        # Normalize if requested
        if normalize and len(y) > 0:
            y = y / y.max()
        
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=mol_id,
            line=dict(color=colors[i % len(colors)], width=1.5),
            fill="tozeroy" if len(ir_data) == 1 else None,
            fillcolor=f"rgba{tuple(list(int(colors[i % len(colors)][j:j+2], 16) for j in (1, 3, 5)) + [0.2])}" if len(ir_data) == 1 else None,
            hovertemplate=f"<b>{mol_id}</b><br>Freq: %{{x:.1f}} cm⁻¹<br>Intensity: %{{y:.2f}}<extra></extra>"
        ))
        
        # Peak labels
        if show_peaks and len(y) > 0:
            peak_idx = y.nlargest(5).index
            for idx in peak_idx:
                fig.add_annotation(
                    x=x.loc[idx],
                    y=y.loc[idx],
                    text=f"{x.loc[idx]:.0f}",
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=0.5,
                    font=dict(size=10)
                )
    
    fig.update_layout(
        title="IR Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Intensity (km/mol)" if not normalize else "Normalized Intensity",
        hovermode="x unified",
        legend=dict(x=0.02, y=0.98),
        xaxis=dict(autorange="reversed" if invert_x else True)
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 View Data Table"):
        from components.data_table import render_data_table
        render_data_table(
            ir_data,
            coordinate_col="freq_cm-1",
            data_cols=["intensity_km/mol", "eps"] if "eps" in next(iter(ir_data.values())).columns else ["intensity_km/mol"],
            title="IR Spectrum Data"
        )


def render_raman_spectra(df: pd.DataFrame, selected_mols: List[str]):
    """Render Raman spectra with customization."""
    
    st.subheader("🟢 Raman Spectrum")
    
    # Customization
    with st.expander("⚙️ Customization", expanded=True):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            freq_range = st.slider("Frequency Range", 0, 4000, (400, 4000), key="raman_freq")
        with col2:
            normalize = st.checkbox("Normalize", False, key="raman_norm")
        with col3:
            show_depol = st.checkbox("Show Depolarization", False, key="raman_depol")
    
    # Collect Raman data
    raman_data = {}
    for mol_id in selected_mols:
        mol_row = df[df["molecule_id"] == mol_id]
        if not mol_row.empty:
            raman = mol_row.iloc[0].get("raman")
            if raman is not None and hasattr(raman, 'empty') and not raman.empty:
                raman_data[mol_id] = raman
    
    if not raman_data:
        st.warning("No Raman data available for selected molecules")
        return
    
    # Create plot
    fig = go.Figure()
    
    colors = ['#00CC96', '#636EFA', '#EF553B', '#AB63FA', '#FFA15A']
    
    for i, (mol_id, raman) in enumerate(raman_data.items()):
        x = raman["freq_cm-1"]
        y = raman["activity"] if "activity" in raman.columns else raman.iloc[:, 1]
        
        mask = (x >= freq_range[0]) & (x <= freq_range[1])
        x, y = x[mask], y[mask]
        
        if normalize and len(y) > 0:
            y = y / y.max()
        
        fig.add_trace(go.Scatter(
            x=x, y=y,
            mode="lines",
            name=mol_id,
            line=dict(color=colors[i % len(colors)], width=1.5),
            hovertemplate=f"<b>{mol_id}</b><br>Freq: %{{x:.1f}} cm⁻¹<br>Activity: %{{y:.4f}}<extra></extra>"
        ))
    
    fig.update_layout(
        title="Raman Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Raman Activity" if not normalize else "Normalized Activity",
        hovermode="x unified",
        legend=dict(x=0.02, y=0.98)
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 View Data Table"):
        from components.data_table import render_data_table
        render_data_table(
            raman_data,
            coordinate_col="freq_cm-1",
            data_cols=["activity", "depolarization"] if show_depol else ["activity"],
            title="Raman Spectrum Data"
        )


def render_uvvis_spectra(df: pd.DataFrame, selected_mols: List[str]):
    """Render UV-Vis spectra from TDDFT data."""
    
    st.subheader("🟣 UV-Vis Spectrum")
    
    # Customization
    with st.expander("⚙️ Customization", expanded=True):
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            wl_range = st.slider("Wavelength Range (nm)", 100, 900, (200, 700), key="uv_wl")
        with col2:
            sigma = st.slider("Broadening σ (nm)", 5, 50, 20, key="uv_sigma")
        with col3:
            show_sticks = st.checkbox("Show Sticks", True, key="uv_sticks")
        with col4:
            show_broadened = st.checkbox("Show Broadened", True, key="uv_broad")
    
    # Collect TDDFT data
    tddft_data = {}
    for mol_id in selected_mols:
        mol_row = df[df["molecule_id"] == mol_id]
        if not mol_row.empty:
            tddft = mol_row.iloc[0].get("tddft_states")
            if tddft is not None and hasattr(tddft, 'empty') and not tddft.empty:
                if "energy_ev" in tddft.columns:
                    tddft_data[mol_id] = tddft
    
    if not tddft_data:
        st.warning("No TDDFT data available for selected molecules")
        return
    
    # Create plot
    fig = go.Figure()
    
    colors = ['#AB63FA', '#636EFA', '#EF553B', '#00CC96', '#FFA15A']
    
    for i, (mol_id, tddft) in enumerate(tddft_data.items()):
        # Convert eV to nm
        wavelength = 1239.84 / tddft["energy_ev"]
        fosc = tddft.get("weight", pd.Series([0.5] * len(tddft))).abs()
        
        # Filter by wavelength
        mask = (wavelength >= wl_range[0]) & (wavelength <= wl_range[1])
        wavelength_f = wavelength[mask]
        fosc_f = fosc[mask]
        
        # Stick spectrum
        if show_sticks:
            for wl, f in zip(wavelength_f, fosc_f):
                fig.add_trace(go.Scatter(
                    x=[wl, wl],
                    y=[0, f],
                    mode="lines",
                    line=dict(color=colors[i % len(colors)], width=2),
                    showlegend=False,
                    hovertemplate=f"<b>{mol_id}</b><br>λ: {wl:.1f} nm<br>f: {f:.4f}<extra></extra>"
                ))
        
        # Broadened spectrum
        if show_broadened:
            wl_grid = np.linspace(wl_range[0], wl_range[1], 500)
            spectrum = np.zeros_like(wl_grid)
            for wl, f in zip(wavelength_f, fosc_f):
                spectrum += f * np.exp(-0.5 * ((wl_grid - wl) / sigma) ** 2)
            
            fig.add_trace(go.Scatter(
                x=wl_grid,
                y=spectrum,
                mode="lines",
                name=mol_id,
                line=dict(color=colors[i % len(colors)], width=1.5),
                fill="tozeroy" if len(tddft_data) == 1 else None,
                hovertemplate=f"<b>{mol_id}</b><br>λ: %{{x:.1f}} nm<br>Intensity: %{{y:.4f}}<extra></extra>"
            ))
    
    fig.update_layout(
        title="UV-Vis Absorption Spectrum",
        xaxis_title="Wavelength (nm)",
        yaxis_title="Oscillator Strength / Intensity",
        hovermode="x unified",
        legend=dict(x=0.02, y=0.98),
        xaxis=dict(range=[wl_range[0], wl_range[1]])
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Data table
    with st.expander("📋 View TDDFT States"):
        from components.data_table import render_data_table
        render_data_table(
            tddft_data,
            coordinate_col="energy_ev",
            data_cols=["weight", "from_orb", "to_orb"] if "from_orb" in next(iter(tddft_data.values())).columns else ["weight"],
            title="TDDFT States"
        )


def render_ir_raman_comparison(df: pd.DataFrame, selected_mols: List[str]):
    """Render IR-Raman dual-axis comparison."""
    
    st.subheader("📊 IR-Raman Comparison")
    
    if len(selected_mols) < 1:
        st.warning("Select at least one molecule")
        return
    
    mol_id = st.selectbox("Select Molecule for Comparison", selected_mols, key="ir_raman_mol")
    
    mol_row = df[df["molecule_id"] == mol_id].iloc[0]
    
    ir = mol_row.get("ir")
    raman = mol_row.get("raman")
    
    if ir is None or raman is None:
        st.warning("Both IR and Raman data required for comparison")
        return
    
    # Create dual-axis plot
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    # IR trace
    fig.add_trace(
        go.Scatter(
            x=ir["freq_cm-1"],
            y=ir["intensity_km/mol"] if "intensity_km/mol" in ir.columns else ir.iloc[:, 1],
            mode="lines",
            name="IR",
            line=dict(color="blue", width=1.5),
            fill="tozeroy",
            fillcolor="rgba(0,0,255,0.1)",
            hovertemplate="<b>IR</b><br>Freq: %{x:.1f} cm⁻¹<br>Intensity: %{y:.2f}<extra></extra>"
        ),
        secondary_y=False
    )
    
    # Raman trace
    fig.add_trace(
        go.Scatter(
            x=raman["freq_cm-1"],
            y=raman["activity"] if "activity" in raman.columns else raman.iloc[:, 1],
            mode="lines",
            name="Raman",
            line=dict(color="green", width=1.5),
            fill="tozeroy",
            fillcolor="rgba(0,255,0,0.1)",
            hovertemplate="<b>Raman</b><br>Freq: %{x:.1f} cm⁻¹<br>Activity: %{y:.4f}<extra></extra>"
        ),
        secondary_y=True
    )
    
    fig.update_layout(
        title=f"IR-Raman Comparison - {mol_id}",
        xaxis_title="Frequency (cm⁻¹)",
        xaxis=dict(autorange="reversed"),
        hovermode="x unified",
        legend=dict(x=0.8, y=0.95)
    )
    fig.update_yaxes(title_text="IR Intensity (km/mol)", secondary_y=False, color="blue")
    fig.update_yaxes(title_text="Raman Activity", secondary_y=True, color="green")
    
    st.plotly_chart(fig, use_container_width=True)
