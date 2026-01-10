"""Energy diagram visualizer with logging."""

from typing import List, Optional
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

from ..core.base_visualizer import BaseVisualizer
from ..logger import get_logger


class EnergyDiagramVisualizer(BaseVisualizer):
    """Create energy comparison diagrams."""
    
    def __init__(self, df: pd.DataFrame, id_column: str = "molecule_id"):
        """
        Initialize visualizer.
        
        Args:
            df: DataFrame with energy data
            id_column: Column with molecule IDs
        """
        super().__init__(df)
        self.id_column = id_column
        self.logger = get_logger("EnergyDiagramVisualizer")
    
    def create_figure(self, energy_col: str = "gibbs_Eh") -> go.Figure:
        """
        Create energy bar chart.
        
        Args:
            energy_col: Column with energy values
            
        Returns:
            Plotly figure
        """
        self.logger.info(f"Creating energy diagram using '{energy_col}'")
        
        if energy_col not in self.data.columns:
            self.logger.warning(f"Column {energy_col} not found")
            return go.Figure()
        
        df = self.data.dropna(subset=[energy_col]).copy()
        
        if df.empty:
            self.logger.warning("No energy data available")
            return go.Figure()
        
        self.logger.info(f"  Plotting {len(df)} molecules")
        
        # Calculate relative energy in kcal/mol
        ref = df[energy_col].min()
        df["relative_kcal"] = (df[energy_col] - ref) * 627.509
        
        self.logger.debug(f"  Reference energy: {ref:.6f} Eh")
        self.logger.debug(f"  Energy range: 0 to {df['relative_kcal'].max():.2f} kcal/mol")
        
        # Create figure
        fig = px.bar(
            df,
            x=self.id_column,
            y="relative_kcal",
            title="Relative Energy",
            labels={
                self.id_column: "Molecule",
                "relative_kcal": "Relative Energy (kcal/mol)"
            },
            color="relative_kcal",
            color_continuous_scale="RdYlGn_r"
        )
        
        fig.update_layout(
            xaxis_tickangle=-45,
            showlegend=False
        )
        
        self._apply_theme(fig)
        self.logger.info("  Energy diagram created successfully")
        
        return fig
    
    def create_pathway_figure(
        self,
        pathway: List[str],
        energy_col: str = "gibbs_Eh"
    ) -> go.Figure:
        """
        Create energy pathway diagram.
        
        Args:
            pathway: List of molecule IDs in pathway order
            energy_col: Column with energy values
            
        Returns:
            Plotly figure
        """
        self.logger.info(f"Creating pathway diagram with {len(pathway)} steps")
        
        # Filter and order by pathway
        df = self.data[self.data[self.id_column].isin(pathway)].copy()
        
        if df.empty:
            self.logger.warning("No molecules found in pathway")
            return go.Figure()
        
        df["order"] = df[self.id_column].apply(lambda x: pathway.index(x))
        df = df.sort_values("order")
        
        # Calculate relative energy
        ref = df[energy_col].min()
        df["relative_kcal"] = (df[energy_col] - ref) * 627.509
        
        self.logger.info(f"  Pathway: {' -> '.join(df[self.id_column].tolist())}")
        
        fig = go.Figure()
        
        # Add line
        fig.add_trace(go.Scatter(
            x=df[self.id_column],
            y=df["relative_kcal"],
            mode="lines+markers",
            name="Energy",
            line=dict(width=3),
            marker=dict(size=12)
        ))
        
        fig.update_layout(
            title="Energy Pathway",
            xaxis_title="Reaction Step",
            yaxis_title="Relative Energy (kcal/mol)",
            xaxis_tickangle=-45
        )
        
        self._apply_theme(fig)
        self.logger.info("  Pathway diagram created successfully")
        
        return fig
