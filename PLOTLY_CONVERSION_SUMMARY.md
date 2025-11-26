# ORCA Visualization Modules - Plotly Conversion Summary

## Conversion Status

This document tracks the conversion of ORCA visualization modules from matplotlib to Plotly.

### Completed Conversions (2/9):
1. ✅ **ir_visualization.py** - Fully converted to Plotly
   - create_stacked_ir_plot() - Interactive stacked IR spectra with hover data
   - create_dual_ir_plot() - Dual-panel absorbance/transmittance with subplots
   - create_ir_with_assignment_table() - Spectrum with functional group table
   - All x-axes inverted for spectroscopy convention
   - Saves as .html or image formats

2. ✅ **raman_visualization.py** - Fully converted to Plotly
   - create_stacked_raman_plot() - Interactive stacked Raman spectra
   - create_single_raman_with_regions() - Single spectrum with annotations
   - Temperature correction support maintained
   - Interactive hover with spectral data

### Remaining Conversions (7/9):
3. ⏳ absorption_visualization.py
4. ⏳ emission_visualization.py
5. ⏳ orbital_visualization.py
6. ⏳ subplot_layouts.py
7. ⏳ experimental_comparison.py
8. ⏳ visualization/fukui_visualization.py
9. ⏳ visualization/mep_visualization.py

## Key Plotly Patterns Used

### Line Plots
```python
fig.add_trace(go.Scatter(
    x=x, y=y,
    mode='lines',
    line=dict(color='blue', width=2),
    hovertemplate='X: %{x:.2f}<br>Y: %{y:.3f}<extra></extra>'
))
```

### Stick Spectra
```python
stick_x = []
stick_y = []
for x_val in x_values:
    stick_x.extend([x_val, x_val, None])
    stick_y.extend([0, height, None])

fig.add_trace(go.Scatter(x=stick_x, y=stick_y, mode='lines'))
```

### Subplots
```python
from plotly.subplots import make_subplots
fig = make_subplots(rows=2, cols=1, subplot_titles=('Top', 'Bottom'))
fig.add_trace(trace1, row=1, col=1)
fig.add_trace(trace2, row=2, col=1)
```

### Inverted X-Axis (Spectroscopy)
```python
fig.update_xaxes(autorange='reversed', range=[max_val, min_val])
```

### Tables
```python
fig.add_trace(go.Table(
    header=dict(values=['Col1', 'Col2'], fill_color='#4472C4'),
    cells=dict(values=[[val1, val2], [val3, val4]])
), row=2, col=1)
```

### Saving
```python
if save_path.endswith('.html'):
    fig.write_html(save_path)
else:
    fig.write_image(save_path, width=1400, height=800, scale=3)
```

## Migration Notes

All converted functions:
- Maintain original function signatures
- Keep all parameters unchanged
- Preserve all docstrings
- Return `go.Figure` instead of `plt.Figure`
- Support both .html and image output formats
- Include interactive hover data with spectral information
- Use inverted x-axis for all spectroscopy plots
- Maintain backward compatibility

