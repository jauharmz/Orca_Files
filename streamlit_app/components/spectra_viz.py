"""
Spectra Visualization Component - Enhanced IR/Raman/UV-Vis

Features:
- Multi-molecule comparison with state labels
- IR spectrum with FWHM broadening, peak labeling, stacked/overlay mode
- Raman spectrum with activity and peak labeling
- UV-Vis from TDDFT with broadening and stacked mode
- Region boundaries with vertical lines
- Invert x-axis toggle
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.signal import find_peaks
from typing import List, Optional, Dict, Tuple


# =============================================================================
# PRESET ASSIGNMENTS FOR PEAK LABELING
# =============================================================================

IR_ASSIGNMENTS = {
    (3200, 3600): "O-H",
    (3000, 3100): "=C-H",
    (2850, 2960): "C-H",
    (2100, 2260): "C≡C/C≡N",
    (1650, 1750): "C=O",
    (1600, 1680): "C=C",
    (1000, 1300): "C-O",
    (650, 900): "C-H bend"
}

RAMAN_ASSIGNMENTS = {
    (3000, 3100): "C-H stretch",
    (2800, 3000): "C-H alkyl",
    (1600, 1700): "C=C",
    (1400, 1500): "C-H bend",
    (1000, 1200): "C-O",
    (600, 800): "C-S",
}

IR_REGION_BOUNDARIES = [4000, 2700, 2000, 1600, 600]
RAMAN_REGION_BOUNDARIES = [4000, 2700, 2000, 1600, 600]


# =============================================================================
# MAIN RENDERER
# =============================================================================

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


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def apply_gaussian_broadening(
    freqs: np.ndarray, 
    intensities: np.ndarray, 
    fwhm: float,
    x_min: float,
    x_max: float,
    n_points: int = 2000
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply Gaussian broadening to discrete peaks.
    
    Args:
        freqs: Peak frequencies
        intensities: Peak intensities
        fwhm: Full Width at Half Maximum
        x_min: Minimum x value for output grid
        x_max: Maximum x value for output grid
        n_points: Number of points in output grid
        
    Returns:
        Tuple of (x_grid, broadened_spectrum)
    """
    if len(freqs) == 0:
        return np.array([]), np.array([])
    
    # Convert FWHM to sigma: σ = FWHM / (2√(2ln2))
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    
    x_grid = np.linspace(x_min, x_max, n_points)
    spectrum = np.zeros_like(x_grid)
    
    for freq, intensity in zip(freqs, intensities):
        if x_min <= freq <= x_max:
            spectrum += intensity * np.exp(-0.5 * ((x_grid - freq) / sigma) ** 2)
    
    return x_grid, spectrum


def find_peak_labels(
    x: np.ndarray, 
    y: np.ndarray, 
    mode: str = "top_n",
    top_n: int = 5,
    threshold: float = 0.3,
    min_distance: float = 30.0,
    assignments: Optional[Dict] = None
) -> List[Dict]:
    """Find peaks and generate labels.
    
    Args:
        x: X-axis values (frequencies)
        y: Y-axis values (intensities)
        mode: Labeling mode ("top_n", "threshold", "all", "none")
        top_n: Number of top peaks to label (for top_n mode)
        threshold: Minimum intensity threshold (for threshold mode)
        min_distance: Minimum distance between peaks
        assignments: Dictionary of (min, max) -> label for functional groups
        
    Returns:
        List of dicts with 'x', 'y', 'label', 'assignment' keys
    """
    if len(x) == 0 or len(y) == 0 or mode == "none":
        return []
    
    # Calculate minimum distance in samples
    dx = np.abs(x[1] - x[0]) if len(x) > 1 else 1.0
    min_dist_samples = max(1, int(min_distance / dx))
    
    # Find peaks
    peaks, properties = find_peaks(y, height=0.01, prominence=0.01, distance=min_dist_samples)
    
    if len(peaks) == 0:
        return []
    
    peak_x = x[peaks]
    peak_y = y[peaks]
    
    # Select peaks based on mode
    if mode == "top_n":
        # Sort by height and take top N
        sorted_idx = np.argsort(peak_y)[::-1][:top_n]
        selected_peaks = sorted(sorted_idx.tolist())
    elif mode == "threshold":
        # Take peaks above threshold
        max_y = np.max(peak_y)
        selected_peaks = [i for i, h in enumerate(peak_y) if h >= threshold * max_y]
    elif mode == "all":
        selected_peaks = list(range(len(peaks)))
    else:
        return []
    
    # Build labels
    labels = []
    for i in selected_peaks:
        freq = peak_x[i]
        intensity = peak_y[i]
        
        # Find assignment
        assignment = None
        if assignments:
            for (min_f, max_f), name in assignments.items():
                if min_f <= freq <= max_f:
                    assignment = name
                    break
        
        labels.append({
            "x": freq,
            "y": intensity,
            "label": f"{freq:.0f}",
            "assignment": assignment
        })
    
    return labels


def add_region_boundaries(fig: go.Figure, boundaries: List[float], y_max: float):
    """Add vertical region boundary lines to a figure."""
    for boundary in boundaries:
        fig.add_vline(
            x=boundary, 
            line=dict(color="gray", width=1, dash="dash"),
            opacity=0.5
        )


# =============================================================================
# IR SPECTRUM
# =============================================================================

def render_ir_spectrum(df: pd.DataFrame, labels: List[str]):
    """Render IR spectrum comparison with enhanced features."""
    
    # Settings
    with st.expander("⚙️ IR Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            freq_range = st.slider("Frequency Range (cm⁻¹)", 0, 4000, (400, 4000), key="ir_freq_range")
            fwhm = st.slider("FWHM Broadening (cm⁻¹)", 5, 100, 20, key="ir_fwhm",
                           help="Full Width at Half Maximum for Gaussian broadening")
        with col2:
            display_mode = st.radio("Display Mode", ["Overlay", "Stacked"], key="ir_display_mode", horizontal=True)
            invert_x = st.checkbox("Invert X-Axis (IR convention)", True, key="ir_invert_x")
        
        col3, col4 = st.columns(2)
        with col3:
            normalize = st.checkbox("Normalize Intensity", True, key="ir_normalize")
            show_regions = st.checkbox("Show Region Boundaries", False, key="ir_regions")
        with col4:
            peak_mode = st.selectbox("Peak Labels", ["none", "top_n", "threshold", "all"], key="ir_peak_mode")
            if peak_mode == "top_n":
                peak_n = st.number_input("Top N Peaks", 1, 20, 5, key="ir_peak_n")
            elif peak_mode == "threshold":
                peak_threshold = st.slider("Threshold (%)", 0, 100, 30, key="ir_peak_threshold") / 100
    
    # Set defaults
    if peak_mode == "top_n":
        peak_params = {"mode": "top_n", "top_n": peak_n}
    elif peak_mode == "threshold":
        peak_params = {"mode": "threshold", "threshold": peak_threshold}
    elif peak_mode == "all":
        peak_params = {"mode": "all"}
    else:
        peak_params = {"mode": "none"}
    
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
    colors = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    y_offset = 0
    y_spacing = 1.2 if display_mode == "Stacked" else 0
    all_annotations = []
    
    for i, (label, ir) in enumerate(ir_data.items()):
        # Get frequency and intensity columns
        if "freq_cm-1" in ir.columns:
            freqs = ir["freq_cm-1"].values
        elif "frequency" in ir.columns:
            freqs = ir["frequency"].values
        else:
            freqs = ir.iloc[:, 0].values
        
        if "intensity_km/mol" in ir.columns:
            intensities = ir["intensity_km/mol"].values
        elif "intensity" in ir.columns:
            intensities = ir["intensity"].values
        elif "eps" in ir.columns:
            intensities = ir["eps"].values
        else:
            intensities = ir.iloc[:, 1].values
        
        # Filter by frequency range
        mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
        freqs = freqs[mask]
        intensities = intensities[mask]
        
        if len(freqs) == 0:
            continue
        
        # Apply Gaussian broadening
        x, y = apply_gaussian_broadening(freqs, intensities, fwhm, freq_range[0], freq_range[1])
        
        # Normalize
        if normalize and len(y) > 0 and np.max(np.abs(y)) > 0:
            y = y / np.max(np.abs(y))
        
        # Apply stacking offset
        y_plot = y + y_offset
        
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(
            x=x, y=y_plot,
            mode='lines',
            name=label,
            line=dict(color=color, width=1.5),
            hovertemplate='%{x:.0f} cm⁻¹<br>%{y:.3f}<extra>%{fullData.name}</extra>'
        ))
        
        # Find and add peak labels
        if peak_params["mode"] != "none":
            peak_labels = find_peak_labels(
                x, y,
                mode=peak_params["mode"],
                top_n=peak_params.get("top_n", 5),
                threshold=peak_params.get("threshold", 0.3),
                assignments=IR_ASSIGNMENTS
            )
            for pl in peak_labels:
                all_annotations.append(dict(
                    x=pl["x"],
                    y=pl["y"] + y_offset + 0.05,
                    text=pl["label"],
                    showarrow=False,
                    font=dict(size=9, color=color),
                    textangle=-90
                ))
        
        # Update offset for stacking
        if display_mode == "Stacked":
            y_offset += (np.max(y) if len(y) > 0 else 1.0) + y_spacing
    
    # Add region boundaries
    if show_regions:
        for boundary in IR_REGION_BOUNDARIES:
            if freq_range[0] <= boundary <= freq_range[1]:
                fig.add_vline(x=boundary, line=dict(color="gray", width=1, dash="dash"), opacity=0.4)
    
    # Add annotations
    for ann in all_annotations:
        fig.add_annotation(**ann)
    
    fig.update_layout(
        title="IR Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Intensity" if not normalize else "Normalized Intensity",
        xaxis=dict(autorange="reversed" if invert_x else True),
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        showlegend=True
    )
    
    if display_mode == "Stacked":
        fig.update_yaxes(showticklabels=False)
    
    st.plotly_chart(fig, use_container_width=True)


# =============================================================================
# RAMAN SPECTRUM
# =============================================================================

def render_raman_spectrum(df: pd.DataFrame, labels: List[str]):
    """Render Raman spectrum comparison with enhanced features."""
    
    with st.expander("⚙️ Raman Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            freq_range = st.slider("Frequency Range (cm⁻¹)", 0, 4000, (400, 4000), key="raman_freq_range")
            fwhm = st.slider("FWHM Broadening (cm⁻¹)", 5, 100, 15, key="raman_fwhm")
        with col2:
            display_mode = st.radio("Display Mode", ["Overlay", "Stacked"], key="raman_display_mode", horizontal=True)
            invert_x = st.checkbox("Invert X-Axis", True, key="raman_invert_x")
        
        col3, col4 = st.columns(2)
        with col3:
            normalize = st.checkbox("Normalize Activity", True, key="raman_normalize")
            show_regions = st.checkbox("Show Region Boundaries", False, key="raman_regions")
        with col4:
            peak_mode = st.selectbox("Peak Labels", ["none", "top_n", "threshold", "all"], key="raman_peak_mode")
            if peak_mode == "top_n":
                peak_n = st.number_input("Top N Peaks", 1, 20, 5, key="raman_peak_n")
            elif peak_mode == "threshold":
                peak_threshold = st.slider("Threshold (%)", 0, 100, 30, key="raman_peak_threshold") / 100
    
    # Set defaults
    if peak_mode == "top_n":
        peak_params = {"mode": "top_n", "top_n": peak_n}
    elif peak_mode == "threshold":
        peak_params = {"mode": "threshold", "threshold": peak_threshold}
    elif peak_mode == "all":
        peak_params = {"mode": "all"}
    else:
        peak_params = {"mode": "none"}
    
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
    colors = ['#00CC96', '#636EFA', '#EF553B', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    y_offset = 0
    y_spacing = 1.2 if display_mode == "Stacked" else 0
    all_annotations = []
    
    for i, (label, raman) in enumerate(raman_data.items()):
        if "freq_cm-1" in raman.columns:
            freqs = raman["freq_cm-1"].values
        else:
            freqs = raman.iloc[:, 0].values
        
        if "activity" in raman.columns:
            activities = raman["activity"].values
        else:
            activities = raman.iloc[:, 1].values
        
        mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
        freqs = freqs[mask]
        activities = activities[mask]
        
        if len(freqs) == 0:
            continue
        
        # Apply Gaussian broadening
        x, y = apply_gaussian_broadening(freqs, activities, fwhm, freq_range[0], freq_range[1])
        
        if normalize and len(y) > 0 and np.max(np.abs(y)) > 0:
            y = y / np.max(np.abs(y))
        
        y_plot = y + y_offset
        
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(x=x, y=y_plot, mode='lines', name=label,
                                line=dict(color=color, width=1.5)))
        
        # Peak labels
        if peak_params["mode"] != "none":
            peak_labels = find_peak_labels(
                x, y,
                mode=peak_params["mode"],
                top_n=peak_params.get("top_n", 5),
                threshold=peak_params.get("threshold", 0.3),
                assignments=RAMAN_ASSIGNMENTS
            )
            for pl in peak_labels:
                all_annotations.append(dict(
                    x=pl["x"],
                    y=pl["y"] + y_offset + 0.05,
                    text=pl["label"],
                    showarrow=False,
                    font=dict(size=9, color=color),
                    textangle=-90
                ))
        
        if display_mode == "Stacked":
            y_offset += (np.max(y) if len(y) > 0 else 1.0) + y_spacing
    
    # Region boundaries
    if show_regions:
        for boundary in RAMAN_REGION_BOUNDARIES:
            if freq_range[0] <= boundary <= freq_range[1]:
                fig.add_vline(x=boundary, line=dict(color="gray", width=1, dash="dash"), opacity=0.4)
    
    for ann in all_annotations:
        fig.add_annotation(**ann)
    
    fig.update_layout(
        title="Raman Spectrum",
        xaxis_title="Frequency (cm⁻¹)",
        yaxis_title="Activity" if not normalize else "Normalized Activity",
        xaxis=dict(autorange="reversed" if invert_x else True),
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    if display_mode == "Stacked":
        fig.update_yaxes(showticklabels=False)
    
    st.plotly_chart(fig, use_container_width=True)


# =============================================================================
# UV-VIS SPECTRUM
# =============================================================================

def render_uvvis_spectrum(df: pd.DataFrame, labels: List[str]):
    """Render UV-Vis spectrum from TDDFT data with enhanced features."""
    
    with st.expander("⚙️ UV-Vis Settings", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            wl_range = st.slider("Wavelength Range (nm)", 100, 900, (200, 700), key="uvvis_wl_range")
            fwhm = st.slider("FWHM Broadening (nm)", 5, 80, 20, key="uvvis_fwhm",
                           help="Full Width at Half Maximum for Gaussian broadening")
        with col2:
            display_mode = st.radio("Display Mode", ["Overlay", "Stacked"], key="uvvis_display_mode", horizontal=True)
            show_sticks = st.checkbox("Show Stick Spectrum", True, key="uvvis_sticks")
        
        col3, col4 = st.columns(2)
        with col3:
            normalize = st.checkbox("Normalize Intensity", True, key="uvvis_normalize")
            show_regions = st.checkbox("Show UV/Vis Regions", False, key="uvvis_regions")
        with col4:
            peak_mode = st.selectbox("Peak Labels", ["none", "top_n", "threshold"], key="uvvis_peak_mode")
            if peak_mode == "top_n":
                peak_n = st.number_input("Top N Peaks", 1, 10, 3, key="uvvis_peak_n")
            elif peak_mode == "threshold":
                peak_threshold = st.slider("Threshold (%)", 0, 100, 30, key="uvvis_peak_threshold") / 100
    
    # Set defaults
    if peak_mode == "top_n":
        peak_params = {"mode": "top_n", "top_n": peak_n}
    elif peak_mode == "threshold":
        peak_params = {"mode": "threshold", "threshold": peak_threshold}
    else:
        peak_params = {"mode": "none"}
    
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
    colors = ['#AB63FA', '#636EFA', '#EF553B', '#00CC96', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880']
    
    y_offset = 0
    y_spacing = 1.2 if display_mode == "Stacked" else 0
    all_annotations = []
    
    # UV-Vis region boundaries
    uvvis_regions = [280, 315, 400, 450, 495, 570, 590, 620]
    
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
        
        # Create broadened spectrum
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        wl_grid = np.linspace(wl_range[0], wl_range[1], 500)
        spectrum = np.zeros_like(wl_grid)
        
        for wl, f in zip(wl_filtered, fosc_filtered):
            spectrum += f * np.exp(-0.5 * ((wl_grid - wl) / sigma) ** 2)
        
        if normalize and np.max(spectrum) > 0:
            spectrum = spectrum / np.max(spectrum)
        
        y_plot = spectrum + y_offset
        
        # Add stick spectrum
        if show_sticks:
            stick_y = fosc_filtered.copy()
            if normalize and np.max(stick_y) > 0:
                stick_y = stick_y / np.max(fosc_filtered) * np.max(spectrum)
            for wl, f in zip(wl_filtered, stick_y):
                fig.add_trace(go.Scatter(
                    x=[wl, wl], y=[y_offset, f + y_offset],
                    mode='lines',
                    line=dict(color=color, width=2),
                    showlegend=False,
                    hovertemplate=f'{label}<br>λ: {wl:.1f} nm<br>f: {f:.3f}<extra></extra>'
                ))
        
        fig.add_trace(go.Scatter(
            x=wl_grid, y=y_plot,
            mode='lines',
            name=label,
            line=dict(color=color, width=1.5)
        ))
        
        # Peak labels
        if peak_params["mode"] != "none":
            peak_labels = find_peak_labels(
                wl_grid, spectrum,
                mode=peak_params["mode"],
                top_n=peak_params.get("top_n", 3),
                threshold=peak_params.get("threshold", 0.3),
                min_distance=15.0
            )
            for pl in peak_labels:
                all_annotations.append(dict(
                    x=pl["x"],
                    y=pl["y"] + y_offset + 0.05,
                    text=f"{pl['x']:.0f}",
                    showarrow=False,
                    font=dict(size=9, color=color),
                    textangle=0
                ))
        
        if display_mode == "Stacked":
            y_offset += (np.max(spectrum) if len(spectrum) > 0 else 1.0) + y_spacing
    
    # Region boundaries
    if show_regions:
        for boundary in uvvis_regions:
            if wl_range[0] <= boundary <= wl_range[1]:
                fig.add_vline(x=boundary, line=dict(color="gray", width=1, dash="dash"), opacity=0.4)
    
    for ann in all_annotations:
        fig.add_annotation(**ann)
    
    fig.update_layout(
        title="UV-Vis Absorption Spectrum",
        xaxis_title="Wavelength (nm)",
        yaxis_title="Intensity / Oscillator Strength" if not normalize else "Normalized Intensity",
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01)
    )
    
    if display_mode == "Stacked":
        fig.update_yaxes(showticklabels=False)
    
    st.plotly_chart(fig, use_container_width=True)
