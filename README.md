# ORCA Quantum Chemistry Visualization Platform

**An interactive Python-based parser and visualization system for ORCA quantum chemistry output files.**

> **Stack**: Streamlit + Plotly + py3Dmol | **Deployment**: Local

---

## 📋 Key Features

| Feature | Description |
|---------|-------------|
| **Modular Parser** | Data-type based parsers with logging |
| **Interactive Viz** | Plotly + py3Dmol visualizations |
| **Multi-Comparison** | Compare multiple molecules side-by-side |
| **Hierarchy Detection** | Auto-detect data hierarchy from naming |
| **Partition Detection** | Auto-detect data partitions (S0/S1/T1, OPT/SP) |
| **Pathway Detection** | Auto-detect degradation pathways |
| **Spectral Scaling** | Linear & relative frequency scaling |
| **Data Export** | Export parsed data (JSON, CSV, Parquet) |
| **HTML Export** | Single-file interactive reports |

---

## 🏗️ System Architecture

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
    end
    
    subgraph "🔬 Analysis Layer"
        HD[HierarchyDetector]
        PRT[PartitionDetector]
        PD[PathwayDetector]
        CE[ComparisonEngine]
        SS[SpectralScaler]
    end
    
    subgraph "📊 Visualization Layer"
        VF[VisualizerFactory]
        M3D[Molecule3D]
        EDV[EnergyDiagram]
        SPV[SpectraPlot]
        PWV[PathwayViz]
    end
    
    subgraph "📤 Export Layer"
        DE[DataExporter]
        HE[HTMLExporter]
    end
    
    F1 & F2 --> PF
    PF --> GP & EP & OP & SP
    GP & EP & OP & SP --> HD & PRT & PD
    HD & PRT --> CE & SS
    PD & CE & SS --> VF
    VF --> M3D & EDV & SPV & PWV
    M3D & EDV & SPV & PWV --> DE & HE
```

---

### Data Detection Architecture

```mermaid
graph TB
    subgraph "📊 Input"
        DF[(Parsed DataFrame)]
    end
    
    subgraph "🔍 Detection Modules"
        HD[HierarchyDetector]
        PRT[PartitionDetector]
        PD[PathwayDetector]
    end
    
    subgraph "📐 Hierarchy"
        H1[Root Nodes]
        H2[Variant Groups]
    end
    
    subgraph "📐 Partitions"
        P1[State: S0/S1/T1]
        P2[Calc Type: OPT/SP]
        P3[ESD Type: VG/AH]
    end
    
    subgraph "📐 Pathways"
        PW1[Edges]
        PW2[Reactions]
    end
    
    DF --> HD --> H1 & H2
    DF --> PRT --> P1 & P2 & P3
    DF & H1 --> PD --> PW1 & PW2
```

**Partition Detection Example:**
```python
detector = PartitionDetector(df)
partitions = detector.detect()

# Output:
# {
#   "by_state": {"S0": [...], "S1": [...], "T1": [...]},
#   "by_calc_type": {"OPT": [...], "SP": [...]},
#   "by_esd_type": {"VG": [...], "AH": [...], "AHAS": [...]}
# }
```

---

### Pathway Detection Architecture

```mermaid
graph LR
    subgraph "📊 Input"
        DF[(DataFrame)]
        HR[Hierarchy]
        PT[Partitions]
    end
    
    subgraph "🔍 Detector"
        PD[PathwayDetector]
        RR[ReactionRules]
        SC[StepCorrections]
    end
    
    subgraph "📐 Output"
        PW[Pathways]
        ED[Edges]
        CS[ColorSchemes]
    end
    
    DF & HR & PT --> PD
    PD --> RR & SC --> PW & ED
    PD --> CS
```

---

## 📁 Project Structure

```
Orca_Files/
├── app.py
├── src/
│   ├── parser/
│   │   ├── factory.py
│   │   ├── geometry.py
│   │   ├── energy.py
│   │   ├── orbitals.py
│   │   └── spectroscopy.py
│   │
│   ├── analysis/
│   │   ├── hierarchy_detector.py
│   │   ├── partition_detector.py    # NEW
│   │   ├── pathway_detector.py
│   │   ├── comparison_engine.py
│   │   └── spectral_scaler.py
│   │
│   ├── viz/
│   │   ├── molecule_3d.py
│   │   ├── energy_diagram.py
│   │   ├── spectra_plot.py
│   │   └── pathway_viz.py
│   │
│   └── export/
│       ├── data_exporter.py
│       └── html_exporter.py
│
└── tests/
```

---

## 🔬 Detection APIs

### Hierarchy
```python
hd = HierarchyDetector(df)
hierarchy = hd.detect()
# p1x, p1a → p1 (root) with variants
```

### Partition
```python
pd = PartitionDetector(df)
partitions = pd.detect()
# by_state, by_calc_type, by_esd_type
```

### Pathway
```python
pwd = PathwayDetector(df, hierarchy, partitions)
pwd.set_reaction_rules({...})
pathways = pwd.detect()
```

---

## 📤 Export APIs

### Data
```python
exporter = DataExporter(df)
exporter.to_json("data.json")
exporter.to_csv("data.csv")
exporter.to_parquet("data.parquet")
```

### HTML
```python
exporter = HTMLExporter(df)
exporter.add_all_visualizations()
exporter.export("report.html")
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
