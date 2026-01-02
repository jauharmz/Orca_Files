"""Orbital plot visualizer."""

from typing import Optional
import pandas as pd
import plotly.graph_objects as go

from ..core.base_visualizer import BaseVisualizer


class OrbitalPlotVisualizer(BaseVisualizer):
    """Create orbital energy level diagrams."""
    
    def create_figure(self, n_orbitals: int = 5) -> go.Figure:
        """
        Create orbital energy level diagram.
        
        Args:
            n_orbitals: Number of HOMO/LUMO orbitals to show
            
        Returns:
            Plotly figure
        """
        orbitals = self.data
        
        if orbitals is None or orbitals.empty:
            self.logger.warning("No orbital data")
            return go.Figure()
        
        # Filter to selected orbitals
        occupied = orbitals[orbitals["OCC"] > 0].nlargest(n_orbitals, "eV")
        virtual = orbitals[orbitals["OCC"] == 0].nsmallest(n_orbitals, "eV")
        
        fig = go.Figure()
        
        # Add occupied orbitals (HOMO)
        for _, orb in occupied.iterrows():
            label = orb.get("label", "HOMO")
            fig.add_trace(go.Bar(
                x=[label],
                y=[orb["eV"]],
                name=label,
                marker_color="blue",
                hovertemplate=f"{label}: {orb['eV']:.3f} eV<extra></extra>"
            ))
        
        # Add virtual orbitals (LUMO)
        for _, orb in virtual.iterrows():
            label = orb.get("label", "LUMO")
            fig.add_trace(go.Bar(
                x=[label],
                y=[orb["eV"]],
                name=label,
                marker_color="red",
                hovertemplate=f"{label}: {orb['eV']:.3f} eV<extra></extra>"
            ))
        
        fig.update_layout(
            title="Orbital Energy Levels",
            yaxis_title="Energy (eV)",
            xaxis_title="Orbital",
            showlegend=False,
            barmode="group"
        )
        
        self._apply_theme(fig)
        self._log_created("orbital plot")
        
        return fig
