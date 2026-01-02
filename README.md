# ORCA Quantum Chemistry Visualization Platform

**An interactive Python-based parser and visualization system for ORCA quantum chemistry output files.**

> **Stack**: Streamlit + Plotly + py3Dmol | **Deployment**: Local | **User**: Single-user

---

## 📋 Key Features

| Feature | Description |
|---------|-------------|
| **Modular Parser** | Data-type based parsers with logging |
| **Interactive Viz** | Plotly + py3Dmol visualizations |
| **Multi-Comparison** | Compare multiple molecules side-by-side |
| **Pathway Detection** | Auto-detect degradation pathways |
| **Spectral Scaling** | Linear & relative frequency scaling |
| **HTML Export** | Single-file interactive reports |

---

## 🏗️ System Architecture

### High-Level Overview

```mermaid
graph TB
    subgraph "📁 Input"
        F1[.out Files]
        F2[HuggingFace]
    end
    
    subgraph "⚙️ Parser Layer"
        PF[ParserFactory]
        GP[GeometryParser]
        EP[EnergyParser]
        OP[OrbitalParser]
        SP[SpectroscopyParser]
        TP[TDDFTParser]
    end
    
    subgraph "🔬 Analysis Layer"
        PD[PathwayDetector]
        CE[ComparisonEngine]
        SS[SpectralScaler]
    end
    
    subgraph "📊 Visualization Layer"
        VF[VisualizerFactory]
        M3D[Molecule3D]
        EDV[EnergyDiagram]
        SPV[SpectraPlot]
        PWV[PathwayVisualizer]
    end
    
    subgraph "📤 Export Layer"
        HE[HTMLExporter]
        PE[PlotExporter]
    end
    
    subgraph "🖥️ Streamlit App"
        UI[StreamlitUI]
    end
    
    F1 & F2 --> PF
    PF --> GP & EP & OP & SP & TP
    GP & EP & OP & SP & TP --> PD & CE & SS
    PD & CE & SS --> VF
    VF --> M3D & EDV & SPV & PWV
    M3D & EDV & SPV & PWV --> HE & PE
    HE & PE --> UI
```

---

### Pathway Detection & Comparison Architecture

```mermaid
graph TB
    subgraph "📊 Data Input"
        DF[(Parsed DataFrame)]
    end
    
    subgraph "🔍 Pathway Detection"
        PD[PathwayDetector]
        NM[NamingMatcher]
        GD[GraphBuilder]
        RR[ReactionRules]
    end
    
    subgraph "📐 Pathway Model"
        PW[Pathway]
        ED[Edge]
        RX[Reaction]
        SC[StepCorrection]
    end
    
    subgraph "🎨 Pathway Visualization"
        PWV[PathwayVisualizer]
        CS[ColorScheme]
        OD[OverlapDetector]
    end
    
    DF --> PD
    PD --> NM --> GD
    GD --> RR
    
    NM --> PW
    GD --> ED
    RR --> RX & SC
    
    PW & ED & RX --> PWV
    PWV --> CS & OD
```

---

### Spectral Scaling Architecture

```mermaid
graph LR
    subgraph "📊 Input"
        RAW[Raw Spectrum]
    end
    
    subgraph "🔧 Scaler"
        SS[SpectralScaler]
        LS[LinearScaler]
        RS[RelativeScaler]
    end
    
    subgraph "📐 Formulas"
        LF["ν_s = s × ν"]
        RF["ν_s = ν_min + s × (ν - ν_min)"]
    end
    
    subgraph "📊 Output"
        SCALED[Scaled Spectrum]
    end
    
    RAW --> SS
    SS --> LS --> LF
    SS --> RS --> RF
    LF & RF --> SCALED
```

---

### HTML Export Architecture

```mermaid
graph TB
    subgraph "📊 Visualizations"
        V1[3D Molecule]
        V2[Energy Diagram]
        V3[Orbital Plot]
        V4[Spectra Plot]
        V5[Pathway Diagram]
    end
    
    subgraph "📤 HTML Exporter"
        HE[HTMLExporter]
        TB[TemplateBuilder]
        JS[PlotlyJS Embed]
        CSS[Styling]
    end
    
    subgraph "📄 Output"
        HTML[Single HTML Report]
    end
    
    V1 & V2 & V3 & V4 & V5 --> HE
    HE --> TB --> JS & CSS
    JS & CSS --> HTML
```

---

## 📁 Project Structure

```
Orca_Files/
├── app.py
├── requirements.txt
│
├── src/
│   ├── core/
│   │   ├── base_parser.py
│   │   ├── base_visualizer.py
│   │   └── data_models.py
│   │
│   ├── parser/
│   │   ├── factory.py
│   │   ├── geometry.py
│   │   ├── energy.py
│   │   ├── orbitals.py
│   │   ├── spectroscopy.py
│   │   ├── tddft.py
│   │   └── batch.py
│   │
│   ├── analysis/                    # NEW: Analysis modules
│   │   ├── __init__.py
│   │   ├── pathway_detector.py      # Auto-detect pathways
│   │   ├── comparison_engine.py     # Multi-molecule comparison
│   │   ├── spectral_scaler.py       # Linear/relative scaling
│   │   └── reaction_rules.py        # Reaction definitions
│   │
│   ├── viz/
│   │   ├── factory.py
│   │   ├── molecule_3d.py
│   │   ├── energy_diagram.py
│   │   ├── orbital_plot.py
│   │   ├── spectra_plot.py
│   │   └── pathway_viz.py           # NEW: Pathway visualization
│   │
│   ├── export/                      # NEW: Export modules
│   │   ├── __init__.py
│   │   ├── html_exporter.py         # Single HTML report
│   │   ├── plot_exporter.py         # PNG/SVG export
│   │   └── templates/
│   │       └── report.html
│   │
│   └── ui/
│       ├── pages/
│       │   ├── home.py
│       │   ├── upload.py
│       │   ├── molecule.py
│       │   ├── compare.py           # NEW: Comparison page
│       │   └── pathway.py           # NEW: Pathway page
│       └── components/
│
└── tests/
```

---

## 🔬 Advanced Features

### 1. Pathway Detection

```python
# Auto-detect degradation pathways from molecule IDs
detector = PathwayDetector(df)
pathways = detector.detect()

# Output:
# [
#   Pathway(nodes=["p1x", "p2x", "p3x", "p4x", "p5x"]),
#   Pathway(nodes=["p1x", "p6x"]),
#   Pathway(nodes=["p1x", "p1a", "p2a", ...]),
# ]

# Define reaction rules
detector.set_reaction_rules({
    ("p1", "p2"): {"add": {"OH": 4}, "remove": {"H2O": 3}},
    ("p2", "p3"): {"add": {"OH": 2}, "remove": {"H2O": 1}},
    ("p3", "p4"): {"remove": {"CO2": 1}},
    ("p4", "p5"): {"add": {"OH": 1}, "remove": {"HCO3": 1}},
    ("p1", "p6"): {"add": {"OH": 1}, "remove": {"H2O": 1, "OMe": 1}},
})

# Color schemes
detector.set_color_scheme("by_variant")  # or "by_destination"
```

### 2. Multi-Comparison

```python
# Compare multiple molecules
engine = ComparisonEngine(df)
comparison = engine.compare(["p1x", "p2x", "p3x"], 
                            properties=["energy", "orbitals", "spectra"])

# Side-by-side visualization
fig = engine.create_comparison_figure()
```

### 3. Spectral Scaling

```python
scaler = SpectralScaler(spectrum_df)

# Linear scaling: ν_s = s × ν
scaled_linear = scaler.linear_scale(factor=0.97)

# Relative scaling: ν_s = ν_min + s × (ν - ν_min)
scaled_relative = scaler.relative_scale(factor=1.5)
```

| Scale Type | Formula | Use Case |
|------------|---------|----------|
| Linear | `ν_s = s × ν` | DFT frequency correction |
| Relative | `ν_s = ν_min + s(ν - ν_min)` | Visual comparison |

### 4. HTML Report Export

```python
exporter = HTMLExporter(df)
exporter.add_molecule_3d(mol_id="p1x")
exporter.add_energy_diagram(mol_ids=["p1x", "p2x"])
exporter.add_pathway_diagram(pathways)
exporter.add_spectra_plot(mol_id="p1x", spectrum_type="ir")

# Export single interactive HTML
exporter.export("report.html")
```

---

## 🧪 Tests

```bash
# Run all tests with output
pytest tests/ -v -s

# Run specific module
pytest tests/test_pathway_detector.py -v
```

---

## 🚀 Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jauharmz/Orca_Files/blob/main/ORCA_Test_v2.ipynb)

---

*Last Updated: 2026-01-02*
