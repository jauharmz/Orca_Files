"""Orbital plot visualizer with logging."""

from typing import Optional
import pandas as pd
import plotly.graph_objects as go

from ..core.base_visualizer import BaseVisualizer
from ..logger import get_logger


class OrbitalPlotVisualizer(BaseVisualizer):
    """Create orbital energy level diagrams."""
    
    def __init__(self, df: pd.DataFrame):
        super().__init__(df)
        self.logger = get_logger("OrbitalPlotVisualizer")
    
    def create_figure(self, n_orbitals: int = 5) -> go.Figure:
        """
        Create orbital energy level diagram.
        
        Args:
            n_orbitals: Number of HOMO/LUMO orbitals to show
            
        Returns:
            Plotly figure
        """
        self.logger.info(f"Creating orbital plot (n_orbitals={n_orbitals})")
        
        orbitals = self.data
        
        if orbitals is None or orbitals.empty:
            self.logger.warning("No orbital data")
            return go.Figure()
        
        # Filter to selected orbitals
        occupied = orbitals[orbitals["OCC"] > 0].nlargest(n_orbitals, "eV")
        virtual = orbitals[orbitals["OCC"] == 0].nsmallest(n_orbitals, "eV")
        
        self.logger.info(f"  Showing {len(occupied)} occupied, {len(virtual)} virtual orbitals")
        
        if not occupied.empty:
            homo = occupied["eV"].max()
            self.logger.debug(f"  HOMO: {homo:.3f} eV")
        if not virtual.empty:
            lumo = virtual["eV"].min()
            self.logger.debug(f"  LUMO: {lumo:.3f} eV")
        
        fig = go.Figure()
        
        # Add occupied orbitals (HOMO)
        for i, (_, orb) in enumerate(occupied.iterrows()):
            lvl = orb.get("lvl", i)
            label = "HOMO" if lvl == 0 else f"HOMO-{lvl}"
            fig.add_trace(go.Bar(
                x=[label],
                y=[orb["eV"]],
                name=label,
                marker_color="blue",
                hovertemplate=f"{label}: {orb['eV']:.3f} eV<extra></extra>"
            ))
        
        # Add virtual orbitals (LUMO)
        for i, (_, orb) in enumerate(virtual.iterrows()):
            lvl = orb.get("lvl", i)
            label = "LUMO" if lvl == 0 else f"LUMO+{lvl}"
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
        self.logger.info("  Orbital plot created successfully")
        
        return fig
