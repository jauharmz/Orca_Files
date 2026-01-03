"""Base visualizer class."""

from typing import Any, Dict, Optional
import plotly.graph_objects as go

from ..logger import get_logger


class BaseVisualizer:
    """
    Base class for all visualizers.
    
    Each visualizer creates specific Plotly figures from parsed ORCA data.
    """
    
    def __init__(self, data: Any, config: Optional[Dict] = None):
        """
        Initialize visualizer.
        
        Args:
            data: Parsed data for visualization
            config: Optional configuration dict
        """
        self.data = data
        self.config = config or {}
        self.logger = get_logger(self.__class__.__name__)
    
    def _apply_theme(self, fig: go.Figure) -> go.Figure:
        """Apply default theme to figure."""
        fig.update_layout(
            template=self.config.get("template", "plotly_white"),
            font=dict(family="Arial", size=12),
            margin=dict(l=60, r=40, t=60, b=40),
        )
        return fig
    
    def _log_created(self, plot_type: str):
        """Log that figure was created."""
        self.logger.info(f"Created {plot_type} figure")
