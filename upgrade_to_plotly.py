"""
Script to upgrade all ORCA visualization files from matplotlib to Plotly.

This script provides helper functions and templates for converting visualization code.
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from typing import List, Tuple, Optional, Union, Dict, Any
import numpy as np

# Template functions for common Plotly operations

def create_plotly_template(template_name: str = 'plotly_white'):
    """
    Get standard Plotly template for publication-quality figures.

    Args:
        template_name: Template name ('plotly', 'plotly_white', 'plotly_dark', 'simple_white')

    Returns:
        Template name
    """
    return template_name


def save_plotly_figure(fig: go.Figure, save_path: Optional[str] = None,
                       format: str = 'html', width: int = 1400, height: int = 800,
                       scale: float = 3.0):
    """
    Save Plotly figure to file.

    Args:
        fig: Plotly figure object
        save_path: Path to save file
        format: 'html', 'png', 'pdf', 'svg', 'jpeg'
        width: Width in pixels
        height: Height in pixels
        scale: Scale factor for static images (higher = better quality)
    """
    if not save_path:
        return

    if format == 'html':
        fig.write_html(save_path)
    else:
        # Requires kaleido package
        fig.write_image(save_path, format=format, width=width, height=height, scale=scale)


def add_vertical_line(fig: go.Figure, x: float, color: str = 'gray',
                      dash: str = 'dash', width: float = 1.0, opacity: float = 0.5,
                      row: Optional[int] = None, col: Optional[int] = None):
    """
    Add vertical line to Plotly figure.

    Args:
        fig: Plotly figure
        x: X position for vertical line
        color: Line color
        dash: Line style ('solid', 'dash', 'dot', 'dashdot')
        width: Line width
        opacity: Line opacity (0-1)
        row: Subplot row (if using subplots)
        col: Subplot column (if using subplots)
    """
    fig.add_vline(
        x=x,
        line_dash=dash,
        line_color=color,
        line_width=width,
        opacity=opacity,
        row=row,
        col=col
    )


def add_horizontal_line(fig: go.Figure, y: float, color: str = 'gray',
                        dash: str = 'dash', width: float = 1.0, opacity: float = 0.5,
                        row: Optional[int] = None, col: Optional[int] = None):
    """
    Add horizontal line to Plotly figure.

    Args:
        fig: Plotly figure
        y: Y position for horizontal line
        color: Line color
        dash: Line style
        width: Line width
        opacity: Line opacity (0-1)
        row: Subplot row (if using subplots)
        col: Subplot column (if using subplots)
    """
    fig.add_hline(
        y=y,
        line_dash=dash,
        line_color=color,
        line_width=width,
        opacity=opacity,
        row=row,
        col=col
    )


def update_layout_for_chemistry(fig: go.Figure,
                                title: str,
                                x_title: str,
                                y_title: str,
                                template: str = 'plotly_white',
                                font_size: int = 12,
                                title_font_size: int = 16,
                                show_grid: bool = True,
                                invert_xaxis: bool = False):
    """
    Apply standard layout settings for chemistry visualizations.

    Args:
        fig: Plotly figure
        title: Plot title
        x_title: X-axis title
        y_title: Y-axis title
        template: Plotly template
        font_size: Base font size
        title_font_size: Title font size
        show_grid: Show grid lines
        invert_xaxis: Invert x-axis (common for spectroscopy)
    """
    fig.update_layout(
        title={
            'text': title,
            'font': {'size': title_font_size, 'weight': 700},
            'x': 0.5,
            'xanchor': 'center'
        },
        xaxis_title={'text': x_title, 'font': {'size': font_size, 'weight': 600}},
        yaxis_title={'text': y_title, 'font': {'size': font_size, 'weight': 600}},
        font={'size': font_size, 'family': 'Arial, sans-serif'},
        template=template,
        hovermode='x unified',
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="gray",
            borderwidth=1
        )
    )

    fig.update_xaxes(
        showgrid=show_grid,
        gridwidth=0.5,
        gridcolor='rgba(128,128,128,0.2)',
        autorange='reversed' if invert_xaxis else True
    )

    fig.update_yaxes(
        showgrid=show_grid,
        gridwidth=0.5,
        gridcolor='rgba(128,128,128,0.2)'
    )


# Matplotlib to Plotly conversion helpers
MATPLOTLIB_TO_PLOTLY_COLORS = {
    'black': '#000000',
    'red': '#FF0000',
    'blue': '#0000FF',
    'green': '#00FF00',
    'orange': '#FFA500',
    'purple': '#800080',
    'brown': '#A52A2A',
    'pink': '#FFC0CB',
    'gray': '#808080',
    'grey': '#808080'
}

MATPLOTLIB_TO_PLOTLY_LINESTYLES = {
    '-': 'solid',
    '--': 'dash',
    '-.': 'dashdot',
    ':': 'dot'
}


def convert_color(mpl_color: str) -> str:
    """Convert matplotlib color to Plotly compatible color."""
    return MATPLOTLIB_TO_PLOTLY_COLORS.get(mpl_color, mpl_color)


def convert_linestyle(mpl_linestyle: str) -> str:
    """Convert matplotlib linestyle to Plotly linestyle."""
    return MATPLOTLIB_TO_PLOTLY_LINESTYLES.get(mpl_linestyle, 'solid')


if __name__ == '__main__':
    print("Plotly upgrade helper module loaded successfully!")
    print("This module provides utilities for converting matplotlib to Plotly visualizations.")
