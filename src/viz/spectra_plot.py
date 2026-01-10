"""Spectra plot visualizer with logging."""

from typing import Optional
import pandas as pd
import plotly.graph_objects as go

from ..core.base_visualizer import BaseVisualizer
from ..logger import get_logger


class SpectraVisualizer(BaseVisualizer):
    """Create spectral plots (IR, Raman, UV-Vis)."""
    
    def __init__(self, df: pd.DataFrame):
        super().__init__(df)
        self.logger = get_logger("SpectraVisualizer")
    
    def create_ir_figure(
        self,
        freq_col: str = "freq_cm-1",
        intensity_col: str = "intensity_km/mol"
    ) -> go.Figure:
        """
        Create IR spectrum plot.
        
        Args:
            freq_col: Frequency column
            intensity_col: Intensity column
            
        Returns:
            Plotly figure
        """
        self.logger.info("Creating IR spectrum plot")
        
        spectrum = self.data
        
        if spectrum is None or spectrum.empty:
            self.logger.warning("No IR spectrum data")
            return go.Figure()
        
        self.logger.info(f"  Plotting {len(spectrum)} peaks")
        
        if freq_col in spectrum.columns:
            freq_range = f"{spectrum[freq_col].min():.0f} - {spectrum[freq_col].max():.0f} cm⁻¹"
            self.logger.debug(f"  Frequency range: {freq_range}")
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=spectrum[freq_col],
            y=spectrum[intensity_col],
            mode="lines",
            name="IR",
            line=dict(width=1.5, color="blue"),
            fill="tozeroy",
            fillcolor="rgba(0,0,255,0.1)"
        ))
        
        fig.update_layout(
            title="IR Spectrum",
            xaxis_title="Frequency (cm⁻¹)",
            yaxis_title="Intensity (km/mol)",
            xaxis=dict(autorange="reversed")  # IR convention
        )
        
        self._apply_theme(fig)
        self.logger.info("  IR spectrum plot created successfully")
        
        return fig
    
    def create_raman_figure(self) -> go.Figure:
        """Create Raman spectrum plot."""
        self.logger.info("Creating Raman spectrum plot")
        
        spectrum = self.data
        
        if spectrum is None or spectrum.empty:
            self.logger.warning("No Raman spectrum data")
            return go.Figure()
        
        self.logger.info(f"  Plotting {len(spectrum)} peaks")
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=spectrum["freq_cm-1"],
            y=spectrum["activity"],
            mode="lines",
            name="Raman",
            line=dict(width=1.5, color="green")
        ))
        
        fig.update_layout(
            title="Raman Spectrum",
            xaxis_title="Frequency (cm⁻¹)",
            yaxis_title="Activity"
        )
        
        self._apply_theme(fig)
        self.logger.info("  Raman spectrum plot created successfully")
        
        return fig
    
    def create_uv_vis_figure(
        self,
        wavelength_col: str = "wavelength_nm",
        intensity_col: str = "fosc"
    ) -> go.Figure:
        """Create UV-Vis absorption spectrum."""
        self.logger.info("Creating UV-Vis absorption spectrum")
        
        spectrum = self.data
        
        if spectrum is None or spectrum.empty:
            self.logger.warning("No UV-Vis spectrum data")
            return go.Figure()
        
        self.logger.info(f"  Plotting {len(spectrum)} transitions")
        
        if wavelength_col in spectrum.columns:
            wl_range = f"{spectrum[wavelength_col].min():.0f} - {spectrum[wavelength_col].max():.0f} nm"
            self.logger.debug(f"  Wavelength range: {wl_range}")
        
        fig = go.Figure()
        
        # Stick spectrum
        for _, row in spectrum.iterrows():
            fig.add_trace(go.Scatter(
                x=[row[wavelength_col], row[wavelength_col]],
                y=[0, row[intensity_col]],
                mode="lines",
                line=dict(width=2, color="purple"),
                showlegend=False,
                hovertemplate=f"{row[wavelength_col]:.1f} nm<br>f={row[intensity_col]:.4f}<extra></extra>"
            ))
        
        fig.update_layout(
            title="UV-Vis Absorption Spectrum",
            xaxis_title="Wavelength (nm)",
            yaxis_title="Oscillator Strength"
        )
        
        self._apply_theme(fig)
        self.logger.info("  UV-Vis spectrum plot created successfully")
        
        return fig
