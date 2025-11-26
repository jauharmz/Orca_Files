# 🚀 Plotly Visualization Upgrade

## Overview

All ORCA visualization modules have been upgraded from **matplotlib** to **Plotly** for modern, interactive, web-ready visualizations.

## ✨ What's New

### Interactive Features
- **Zoom and Pan**: Click and drag to zoom, double-click to reset
- **Hover Data**: Hover over peaks to see exact values
- **Toggle Traces**: Click legend items to show/hide datasets
- **Export Options**: Download plots as PNG, SVG, or interactive HTML
- **Responsive**: Automatically adapts to screen size

### Web Integration
- **Flask API**: New visualization endpoints serve interactive Plotly charts
- **HTML Export**: Generate standalone HTML visualizations
- **CDN Support**: Lightweight delivery via Plotly.js CDN

## 📦 Updated Dependencies

```bash
pip install plotly==5.18.0 kaleido==0.2.1 numpy scipy
```

See `requirements.txt` for complete list.

## 🔄 Migration Guide

### For Existing Code

All function signatures remain **identical** - just update your imports:

#### Before (matplotlib):
```python
from ir_visualization import create_stacked_ir_plot
import matplotlib.pyplot as plt

fig = create_stacked_ir_plot(datasets, labels)
fig.savefig('output.png', dpi=300)
plt.show()
```

#### After (Plotly):
```python
from ir_visualization import create_stacked_ir_plot

fig = create_stacked_ir_plot(datasets, labels)
# Save as interactive HTML (recommended)
fig.write_html('output.html')
# Or save as static image
fig.write_image('output.png', width=1400, height=800, scale=3)
# Display in notebook/browser
fig.show()
```

### Return Type Change

- **Before**: Functions returned `matplotlib.figure.Figure`
- **After**: Functions return `plotly.graph_objects.Figure`

All other parameters and functionality remain the same.

## 📁 Converted Modules

### Core Visualization
| Module | Functions | Status |
|--------|-----------|--------|
| `visualization_utils.py` | Core utilities | ✅ Updated |
| `ir_visualization.py` | IR spectra | ✅ Converted |
| `raman_visualization.py` | Raman spectra | ✅ Converted |
| `absorption_visualization.py` | UV-Vis absorption | ✅ Converted |
| `emission_visualization.py` | Fluorescence/phosphorescence | ✅ Converted |
| `orbital_visualization.py` | Orbital energy diagrams | ✅ Converted |
| `experimental_comparison.py` | Exp vs DFT comparison | ✅ Converted |
| `subplot_layouts.py` | Multi-panel layouts | ✅ Converted |

### Advanced Visualization
| Module | Functions | Status |
|--------|-----------|--------|
| `visualization/fukui_visualization.py` | Fukui indices | ✅ Converted |
| `visualization/mep_visualization.py` | MEP analysis | ✅ Converted |

### Web Application
| File | Status |
|------|--------|
| `app.py` | ✅ Updated with Plotly endpoints |

## 🌐 Flask API Endpoints

New visualization endpoints:

```
GET /api/visualize/ir          - Interactive IR spectrum
GET /api/visualize/raman       - Interactive Raman spectrum
GET /api/visualize/orbitals    - Orbital energy diagram
GET /api/visualize/absorption  - UV-Vis absorption spectrum
```

Each endpoint returns interactive HTML that can be embedded in web pages.

## 📖 Usage Examples

### 1. IR Spectra

```python
from ir_visualization import create_stacked_ir_plot, create_dual_ir_plot
import numpy as np

# Single IR spectrum with absorbance and transmittance
frequencies = np.array([1000, 1500, 2000, 2900, 3300])
intensities = np.array([50, 120, 80, 200, 150])

fig = create_dual_ir_plot(
    frequencies,
    intensities,
    title="My IR Spectrum",
    show_regions=True
)

# Save as interactive HTML
fig.write_html("ir_spectrum.html")

# Or display in Jupyter notebook
fig.show()
```

### 2. Raman Spectra

```python
from raman_visualization import create_stacked_raman_plot

datasets = [
    {'frequencies': [500, 1000, 1500], 'activities': [10, 50, 30]},
    {'frequencies': [520, 1020, 1480], 'activities': [12, 48, 35]}
]
labels = ['Sample A', 'Sample B']

fig = create_stacked_raman_plot(
    datasets,
    labels,
    use_temperature_correction=True,
    show_regions=True
)

fig.write_html("raman_comparison.html")
```

### 3. Orbital Energy Diagrams

```python
from orbital_visualization import create_orbital_comparison

# Data structure from ORCA parser
dataset = {
    'orbital_energies': [
        {'energy_ev': -10.5, 'occupation': 2.0},
        {'energy_ev': -8.2, 'occupation': 2.0},
        {'energy_ev': -5.5, 'occupation': 2.0},  # HOMO
        {'energy_ev': -2.1, 'occupation': 0.0},  # LUMO
        {'energy_ev': 1.5, 'occupation': 0.0},
    ]
}

fig = create_orbital_comparison(
    [dataset],
    ['DFT Calculation'],
    n_orbitals=10,
    show_gap_values=True
)

fig.write_html("orbitals.html")
```

### 4. Fukui Function Analysis

```python
from visualization.fukui_visualization import create_fukui_bar_plot

# fukui_data from parser
fig = create_fukui_bar_plot(
    fukui_data,
    top_n=15,
    show_labels=True
)

fig.write_html("fukui_indices.html")
```

### 5. MEP Visualization

```python
from visualization.mep_visualization import create_mep_slice_plot

fig = create_mep_slice_plot(
    mep_data,
    slice_axis='z',
    slice_position=0.0,
    show_atoms=True
)

fig.write_html("mep_slice.html")
```

## 🎨 Customization

All Plotly figures can be customized after creation:

```python
fig = create_stacked_ir_plot(datasets, labels)

# Customize layout
fig.update_layout(
    template='plotly_dark',  # Change theme
    font=dict(size=14, family='Arial'),
    height=900,
    width=1600
)

# Customize axes
fig.update_xaxes(
    range=[400, 4000],
    showgrid=True,
    gridcolor='rgba(200,200,200,0.3)'
)

# Add annotations
fig.add_annotation(
    x=1700,
    y=0.5,
    text="C=O stretch",
    showarrow=True,
    arrowhead=2
)

fig.write_html("customized.html")
```

## 🖼️ Saving Formats

### Interactive HTML (Recommended)
```python
fig.write_html("plot.html")
```
- Full interactivity preserved
- Self-contained file
- Can be opened in any browser
- Perfect for sharing and reports

### Static Images
```python
# High-resolution PNG
fig.write_image("plot.png", width=1400, height=800, scale=3)

# Vector formats (requires kaleido)
fig.write_image("plot.svg")
fig.write_image("plot.pdf")
fig.write_image("plot.eps")
```

## 📊 Features by Module

### IR Visualization
- ✅ Stacked multi-dataset comparison
- ✅ Dual absorbance/transmittance plots
- ✅ Functional group assignment
- ✅ Frequency scaling corrections
- ✅ Regional boundary markers
- ✅ Stick spectrum overlay
- ✅ Interactive hover with peak values

### Raman Visualization
- ✅ Temperature-dependent intensity
- ✅ Multi-dataset stacking
- ✅ Regional annotations
- ✅ Stick spectrum overlay
- ✅ Normalization options
- ✅ Custom color palettes

### Absorption/Emission
- ✅ Gaussian-broadened spectra
- ✅ Oscillator strength display
- ✅ Mirror image plots (Stokes shift)
- ✅ Multi-method comparison
- ✅ Electric vs velocity dipole
- ✅ Transition table overlays

### Orbital Diagrams
- ✅ Side-by-side comparisons
- ✅ HOMO-LUMO gap annotations
- ✅ Connection lines between orbitals
- ✅ Orbital filtering
- ✅ Energy range selection

### Fukui/Reactivity
- ✅ Bar plots for f⁺, f⁻, f⁰
- ✅ Dual descriptor visualization
- ✅ Top-N reactive sites
- ✅ Color-coded reactivity
- ✅ Heatmap visualizations

### MEP Analysis
- ✅ 2D slice plots
- ✅ 3D surface rendering
- ✅ Contour plots
- ✅ Critical point identification
- ✅ Atom position overlays

## 🔧 Troubleshooting

### Image Export Not Working

If `fig.write_image()` fails:
```bash
pip install kaleido==0.2.1
```

### Plots Not Showing in Jupyter

```python
import plotly.io as pio
pio.renderers.default = "notebook"  # or "browser"
```

### Memory Issues with Large Data

Use sampling or reduce resolution:
```python
fig = create_broadened_spectrum(
    wavelengths,
    intensities,
    num_points=1000  # Reduce from default 2000
)
```

## 📝 Breaking Changes

### Save Method
- **Old**: `fig.savefig(path, dpi=300)`
- **New**: `fig.write_html(path)` or `fig.write_image(path, scale=3)`

### Show Method
- **Old**: `plt.show()`
- **New**: `fig.show()`

### Return Type
- **Old**: `matplotlib.figure.Figure`
- **New**: `plotly.graph_objects.Figure`

## 🎯 Best Practices

1. **Use HTML for sharing**: Interactive HTML files are self-contained and work everywhere
2. **High DPI for publications**: Use `scale=3` or higher for publication-quality images
3. **Leverage interactivity**: Add hover templates with detailed information
4. **Responsive layouts**: Don't hardcode sizes, let Plotly adapt
5. **Use templates**: `plotly_white` for light backgrounds, `plotly_dark` for presentations

## 📚 Resources

- [Plotly Python Documentation](https://plotly.com/python/)
- [Plotly Graph Objects Reference](https://plotly.com/python/graph-objects/)
- [Plotly Styling Guide](https://plotly.com/python/templates/)
- [Kaleido for Static Export](https://github.com/plotly/Kaleido)

## 🤝 Contributing

When adding new visualizations:
1. Use `plotly.graph_objects` for maximum control
2. Return `go.Figure` objects
3. Support both `.write_html()` and `.write_image()` saving
4. Add hover templates with meaningful data
5. Follow the existing API patterns
6. Document with clear examples

## ⚡ Performance Notes

- Plotly handles up to ~100K points efficiently
- For larger datasets, consider binning or sampling
- HTML files with many traces can be large (10-50 MB)
- Static images render faster than HTML for very large datasets

## 🙋 Support

For issues or questions:
- Check the examples in `/examples/` directory
- Review function docstrings
- See Plotly documentation for advanced customization

---

**Upgrade completed**: January 2025
**Plotly version**: 5.18.0
**Maintains compatibility**: Yes - all function signatures preserved
