# ORCA Quantum Chemistry Visualization Platform

**An interactive Python-based parser and visualization system for ORCA quantum chemistry output files.**

> **Stack**: Streamlit + Plotly + py3Dmol | **Deployment**: Local | **User**: Single-user project

---

## 📋 Project Overview

| Component | Technology | Status |
|-----------|------------|--------|
| **Parser** | Python, regex, pandas | 🔄 Refactoring |
| **Visualization** | Plotly, py3Dmol | ⏳ Migration |
| **Web UI** | Streamlit | ⏳ Planned |
| **Logging** | Python logging | ⏳ Planned |

---

## 🏗️ System Architecture

### High-Level Overview

```mermaid
graph TB
    subgraph "📁 Input Sources"
        F1[Local .out Files]
        F2[Folder Batch Upload]
        F3[HuggingFace Dataset]
    end
    
    subgraph "🔄 File Handler"
        FH[FileHandler]
        FV[FileValidator]
        FL[FileLoader]
    end
    
    subgraph "⚙️ Parser Registry"
        PR[ParserRegistry]
        PF[ParserFactory]
    end
    
    subgraph "📊 Data Store"
        DS[(MoleculeDataStore)]
        DF[(DataFrame Cache)]
    end
    
    subgraph "📈 Visualization Registry"
        VR[VisualizerRegistry]
        VF[VisualizerFactory]
    end
    
    subgraph "🖥️ Streamlit App"
        UI[StreamlitUI]
        SS[SessionState]
    end
    
    F1 & F2 & F3 --> FH
    FH --> FV --> FL
    FL --> PR --> PF
    PF --> DS --> DF
    DF --> VR --> VF
    VF --> UI
    UI <--> SS
```

---

### Parser Module Architecture

```mermaid
graph LR
    subgraph "📄 Input"
        TXT[Raw ORCA Text]
    end
    
    subgraph "🏭 Parser Factory"
        PF[ParserFactory]
    end
    
    subgraph "🔧 Base Parser"
        BP[BaseParser]
        BP --> LOG[Logger]
        BP --> UTIL[RegexUtils]
    end
    
    subgraph "📦 Data Parsers"
        GP[GeometryParser]
        EP[EnergyParser]
        OP[OrbitalParser]
        SP[SpectroscopyParser]
        TP[TDDFTParser]
        DP[DipoleParser]
        MP[MullikenParser]
    end
    
    subgraph "📊 Output Models"
        GM[GeometryData]
        EM[EnergyData]
        OM[OrbitalData]
        SM[SpectraData]
        TM[TDDFTData]
        DM[DipoleData]
        MM[ChargeData]
    end
    
    TXT --> PF
    PF --> GP & EP & OP & SP & TP & DP & MP
    
    GP -.-> BP
    EP -.-> BP
    OP -.-> BP
    SP -.-> BP
    TP -.-> BP
    DP -.-> BP
    MP -.-> BP
    
    GP --> GM
    EP --> EM
    OP --> OM
    SP --> SM
    TP --> TM
    DP --> DM
    MP --> MM
```

---

### Visualization Module Architecture

```mermaid
graph LR
    subgraph "📊 Data Input"
        DF[(DataFrame)]
    end
    
    subgraph "🏭 Visualizer Factory"
        VF[VisualizerFactory]
    end
    
    subgraph "🎨 Base Visualizer"
        BV[BaseVisualizer]
        BV --> CFG[PlotConfig]
        BV --> THM[ThemeManager]
    end
    
    subgraph "📈 Visualizer Modules"
        M3D[Molecule3DVisualizer]
        EV[EnergyDiagramVisualizer]
        OV[OrbitalPlotVisualizer]
        SV[SpectraVisualizer]
        TV[TDDFTVisualizer]
    end
    
    subgraph "🖼️ Output"
        P3D[py3Dmol View]
        PL3D[Plotly 3D Scatter]
        PLB[Plotly Bar Chart]
        PLL[Plotly Line Chart]
        PLS[Plotly Sankey]
    end
    
    DF --> VF
    VF --> M3D & EV & OV & SV & TV
    
    M3D -.-> BV
    EV -.-> BV
    OV -.-> BV
    SV -.-> BV
    TV -.-> BV
    
    M3D --> P3D & PL3D
    EV --> PLB
    OV --> PLB
    SV --> PLL
    TV --> PLS
```

---

### Streamlit App Architecture

```mermaid
graph TB
    subgraph "🖥️ Streamlit Pages"
        HP[HomePage]
        UP[UploadPage]
        MP[MoleculePage]
        CP[ComparePage]
        SP[SettingsPage]
    end
    
    subgraph "🧩 UI Components"
        FU[FileUploader]
        MS[MoleculeSelector]
        VT[VisualizationTabs]
        EX[ExportButtons]
        LG[LogViewer]
    end
    
    subgraph "📦 State Management"
        SS[SessionState]
        DC[DataCache]
        VC[ViewCache]
    end
    
    subgraph "🔧 Services"
        PS[ParserService]
        VS[VizService]
        LS[LogService]
    end
    
    HP --> FU
    UP --> FU --> MS
    MP --> VT --> EX
    SP --> LG
    
    FU --> PS
    VT --> VS
    LG --> LS
    
    PS & VS & LS <--> SS
    SS <--> DC & VC
```

---

### Data Flow Architecture

```mermaid
sequenceDiagram
    participant U as User
    participant ST as Streamlit
    participant FH as FileHandler
    participant PF as ParserFactory
    participant DS as DataStore
    participant VF as VizFactory
    participant LOG as Logger
    
    U->>ST: Upload .out file(s)
    ST->>FH: validate_files()
    FH->>LOG: log("Validating files...")
    FH-->>ST: valid_files[]
    
    loop For each file
        ST->>PF: parse(file)
        PF->>LOG: log("Parsing geometry...")
        PF->>LOG: log("Parsing energy...")
        PF->>LOG: log("Parsing orbitals...")
        PF-->>DS: MoleculeData
    end
    
    DS-->>ST: DataFrame
    ST->>LOG: log(f"Parsed {n} molecules")
    
    U->>ST: Select molecule
    ST->>VF: create_visualizations(mol_id)
    VF->>LOG: log("Creating 3D viewer...")
    VF->>LOG: log("Creating energy plot...")
    VF-->>ST: [Figure, Figure, ...]
    
    ST-->>U: Display visualizations
```

---

## 📁 Project Structure

```
Orca_Files/
├── README.md
├── requirements.txt
├── app.py                          # Streamlit entry point
│
├── src/
│   ├── __init__.py
│   ├── config.py                   # Configuration constants
│   ├── logger.py                   # Centralized logging
│   │
│   ├── core/                       # Core abstractions
│   │   ├── __init__.py
│   │   ├── base_parser.py          # Abstract parser
│   │   ├── base_visualizer.py      # Abstract visualizer
│   │   ├── data_models.py          # Pydantic/dataclass models
│   │   └── exceptions.py           # Custom exceptions
│   │
│   ├── parser/                     # Modular parsers
│   │   ├── __init__.py
│   │   ├── registry.py             # Parser registration
│   │   ├── factory.py              # ParserFactory
│   │   ├── geometry.py             # Coords, SMILES
│   │   ├── energy.py               # Gibbs, single-point
│   │   ├── orbitals.py             # HOMO/LUMO, spin
│   │   ├── spectroscopy.py         # IR, Raman, NMR
│   │   ├── tddft.py                # TD-DFT states
│   │   ├── dipole.py               # Electric/velocity
│   │   ├── mulliken.py             # Charges
│   │   └── batch.py                # Multi-file
│   │
│   ├── viz/                        # Modular visualizers
│   │   ├── __init__.py
│   │   ├── registry.py             # Visualizer registration
│   │   ├── factory.py              # VisualizerFactory
│   │   ├── config.py               # Plot themes/configs
│   │   ├── molecule_3d.py          # py3Dmol + Plotly 3D
│   │   ├── energy_diagram.py       # Energy pathways
│   │   ├── orbital_plot.py         # Orbital bars
│   │   ├── spectra_plot.py         # IR/Raman/UV-Vis
│   │   └── tddft_plot.py           # TD-DFT diagrams
│   │
│   ├── services/                   # Business logic
│   │   ├── __init__.py
│   │   ├── parser_service.py       # Parsing orchestration
│   │   ├── viz_service.py          # Viz orchestration
│   │   └── file_service.py         # File handling
│   │
│   ├── ui/                         # Streamlit components
│   │   ├── __init__.py
│   │   ├── components/
│   │   │   ├── file_uploader.py
│   │   │   ├── molecule_selector.py
│   │   │   ├── viz_tabs.py
│   │   │   └── export_panel.py
│   │   └── pages/
│   │       ├── home.py
│   │       ├── upload.py
│   │       ├── molecule.py
│   │       └── settings.py
│   │
│   └── utils/
│       ├── __init__.py
│       ├── converters.py           # Unit conversions
│       ├── validators.py           # Input validation
│       └── regex_patterns.py       # Shared regex
│
├── tests/
│   ├── __init__.py
│   ├── conftest.py                 # Pytest fixtures
│   ├── test_parsers/
│   │   ├── test_geometry.py
│   │   ├── test_energy.py
│   │   └── ...
│   ├── test_visualizers/
│   └── test_data/
│       └── *.out
│
├── notebooks/
│   ├── 01_parser_demo.ipynb
│   └── 02_viz_demo.ipynb
│
└── legacy/
    ├── orca_praser.py
    └── *.ipynb
```

---

## 🔧 Parser Module Details

| Parser | Input Regex Pattern | Output Fields |
|--------|---------------------|---------------|
| `GeometryParser` | `CARTESIAN COORDINATES` | `cart_coords`, `internal_coords`, `smiles` |
| `EnergyParser` | `Final Gibbs`, `FINAL SINGLE POINT` | `gibbs_Eh`, `single_point_Eh` |
| `OrbitalParser` | `ORBITAL ENERGIES`, `SPIN UP/DOWN` | `orbitals` DataFrame |
| `SpectroscopyParser` | `IR SPECTRUM`, `RAMAN`, `NMR` | `ir`, `raman`, `nmr` |
| `TDDFTParser` | `TD-DFT EXCITED STATES` | `tddft_states` DataFrame |
| `DipoleParser` | `ELECTRIC DIPOLE`, `VELOCITY DIPOLE` | `electric_dipole`, `velocity_dipole` |
| `MullikenParser` | `MULLIKEN POPULATION` | `mulliken_charges` |

---

## 📊 Visualizer Module Details

| Visualizer | Input Data | Output Type | Interactive Features |
|------------|------------|-------------|---------------------|
| `Molecule3DVisualizer` | `cart_coords` | py3Dmol, Plotly 3D | Rotate, zoom, style toggle |
| `EnergyDiagramVisualizer` | `energies` | Plotly bar | Hover values, pathway click |
| `OrbitalPlotVisualizer` | `orbitals` | Plotly bar | Hover labels, gap highlight |
| `SpectraVisualizer` | `ir`, `raman` | Plotly line | Peak labels, zoom |
| `TDDFTVisualizer` | `tddft_states` | Plotly sankey | Transition hover |

---

## 🪵 Logging Configuration

```python
# Log format
%(asctime)s | %(name)s | %(levelname)s | %(message)s

# Example output
2026-01-02 18:55:00 | GeometryParser | INFO | Found 42 atoms
2026-01-02 18:55:00 | EnergyParser | INFO | Gibbs energy: -854.123456 Eh
2026-01-02 18:55:01 | Molecule3DVisualizer | DEBUG | Creating py3Dmol view
```

---

## 🚀 Quick Start

```bash
# Install
pip install -r requirements.txt

# Run tests with output
pytest tests/ -v -s

# Run Streamlit app
streamlit run app.py
```

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jauharmz/Orca_Files/blob/main/ORCA_Test_v2.ipynb)

---

## 🤗 HuggingFace Sample Data

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="JauharMz/Orca",
    repo_type="dataset",
    local_dir="./data"
)
```

---

*Last Updated: 2026-01-02*
