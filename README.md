# ORCA Quantum Chemistry Visualization Platform

**An interactive Python-based parser and visualization system for ORCA quantum chemistry output files.**

> **Stack**: Streamlit + Plotly + py3Dmol | **Deployment**: Local | **Logging**: Full debug support

---

## 📋 Key Features

| Feature | Description |
|---------|-------------|
| **Modular Parser** | Data-type based parsers (geometry, energy, orbitals, spectroscopy, tddft, dipole, mulliken) |
| **Hierarchy Detection** | Auto-detect molecule hierarchy from naming patterns (p1x, p1a, p1b → p1 root) |
| **Partition Detection** | Auto-detect partitions by state (S0/S1/T1), calc type (OPT/SP), ESD type (VG/AH/AHAS) |
| **Pathway Detection** | Auto-detect degradation pathways with reaction rules and step corrections |
| **Multi-Comparison** | Compare multiple molecules side-by-side (Energies, Spectra, Orbitals) |
| **Spectral Scaling** | Linear (`ν_s = s × ν`) and relative (`ν_s = ν_min + s(ν - ν_min)`) scaling |
| **Data Export** | Export parsed data to JSON, CSV, Parquet, Pickle |
| **Interactive HTML** | **NEW:** Generate self-contained research papers with 3D views, spectra overlay, and embedded raw data. |
| **Interactive Viz** | Plotly charts + py3Dmol 3D molecular viewer + RDKit 2D Structures |

---

## 🏗️ System Architecture

### Complete System Overview

```mermaid
graph TB
    subgraph "📁 Input Sources"
        FS[Local Files]
        FB[Folder Batch]
        HF[HuggingFace Dataset]
    end
    
    subgraph "🔄 File Handler"
        FH[FileHandler]
        FV[FileValidator]
        FL[FileLoader]
    end
    
    subgraph "⚙️ Parser Layer"
        direction TB
        PR[ParserRegistry]
        PF[ParserFactory]
        
        subgraph "Parsers"
            GP[GeometryParser]
            EP[EnergyParser]
            OP[OrbitalParser]
            SCP[SpectroscopyParser]
            TP[TDDFTParser]
            DP[DipoleParser]
            MP[MullikenParser]
        end
    end
    
    subgraph "🔬 Analysis Layer"
        direction TB
        
        subgraph "Detection"
            HD[HierarchyDetector]
            PRT[PartitionDetector]
            PWD[PathwayDetector]
        end
        
        subgraph "Processing"
            CE[ComparisonEngine]
            SS[SpectralScaler]
            RR[ReactionRules]
        end
    end
    
    subgraph "📊 Data Layer"
        DS[(DataStore)]
        DC[DataCache]
        DM[DataModels]
    end
    
    subgraph "📈 Visualization Layer"
        direction TB
        VR[VisualizerRegistry]
        VF[VisualizerFactory]
        
        subgraph "Visualizers"
            M3D[Molecule3DVisualizer]
            EDV[EnergyDiagramVisualizer]
            OPV[OrbitalPlotVisualizer]
            SPV[SpectraVisualizer]
            PWV[PathwayVisualizer]
            HTV[HierarchyTreeVisualizer]
        end
    end
    
    subgraph "📤 Export Layer"
        DE[DataExporter]
        HE[HTMLExporter]
        PE[PlotExporter]
    end
    
    subgraph "🖥️ Streamlit Application"
        direction TB
        APP[streamlit_app/app.py]
        
        subgraph "Components"
            C_VIZ[Viz Components]
            C_EXP[ExportPanel]
            C_SEL[MoleculeSelector]
            C_3D[3D Viewer]
        end
    end
    
    subgraph "🪵 Logging"
        LOG[Logger]
        LF[LogFormatter]
        LH[LogHandlers]
    end
    
    FS & FB & HF --> FH
    FH --> FV --> FL
    FL --> PR --> PF
    PF --> GP & EP & OP & SCP & TP & DP & MP
    GP & EP & OP & SCP & TP & DP & MP --> DS
    DS --> HD & PRT
    HD & PRT --> PWD
    PWD --> RR
    DS --> CE & SS
    DS & HD & PRT & PWD --> DC
    DC --> VR --> VF
    VF --> M3D & EDV & OPV & SPV & PWV & HTV
    M3D & EDV & OPV & SPV & PWV & HTV --> DE & HE & PE
    DE & HE & PE --> APP
    APP --> C_SEL & C_VIZ & C_EXP & C_3D
    
    GP & EP & HD & PWD & M3D --> LOG
    LOG --> LF --> LH
```

---

### Parser Module Architecture

The modular parser is refactored from `orca_praser.py` into independent, focused modules:

```mermaid
graph TB
    subgraph "📄 Input"
        TXT[Raw ORCA Text]
        FN[Filename]
    end
    
    subgraph "🏭 Factory"
        PF[ParserFactory]
        BP[BatchParser]
    end
    
    subgraph "🔧 Core"
        BASE[BaseParser]
        LOG[Logger]
        RX[RegexPatterns]
        DM[DataModels]
    end
    
    subgraph "📦 Modular Parsers"
        GP[GeometryParser]
        EP[EnergyParser]
        OP[OrbitalParser]
        SCP[SpectroscopyParser]
        TP[TDDFTParser]
        SFP[SpectrumFileParser]
    end
    
    subgraph "📊 Data Models"
        GMD[GeometryData]
        EMD[EnergyData]
        OMD[OrbitalData]
        SMD[SpectraData]
        TMD[TDDFTData]
        MMD[MullikenData]
        IMD[InternalCoordsData]
    end
    
    subgraph "📋 Output"
        RES[ParseResult]
        DF[(DataFrame)]
        CSV[(CSV Files)]
        JSON[(JSON)]
    end
    
    TXT & FN --> PF
    PF --> GP & EP & OP & SCP & TP
    BP --> PF
    
    GP & EP & OP & SCP & TP -.-> BASE
    BASE --> LOG & RX
    
    GP --> GMD & IMD
    EP --> EMD
    OP --> OMD
    SCP --> SMD & MMD
    TP --> TMD
    
    GMD & EMD & OMD & SMD & TMD & MMD & IMD --> RES
    RES --> DF --> CSV & JSON
```

#### Data Components (28 fields)

| Category | Field | Source Module | Description |
|----------|-------|---------------|-------------|
| **Identity** | molecule_id | batch.py | Extracted from filename |
| | smiles | geometry.py | SMILES from coordinates |
| | charge | geometry.py | Molecular charge |
| | multiplicity | geometry.py | Spin multiplicity |
| **Energy** | gibbs_Eh | energy.py | Gibbs free energy (Eh) |
| | single_point_Eh | energy.py | Single-point energy (Eh) |
| **Orbitals** | homo_energy | orbitals.py | HOMO energy (eV) |
| | lumo_energy | orbitals.py | LUMO energy (eV) |
| | homo_lumo_gap | orbitals.py | HOMO-LUMO gap (eV) |
| | orbitals | orbitals.py | Full orbital DataFrame |
| **Geometry** | cart_coords | geometry.py | Cartesian coordinates |
| | bonds | geometry.py | Internal bond coords |
| | angles | geometry.py | Internal angle coords |
| | dihedrals | geometry.py | Internal dihedral coords |
| **Spectroscopy** | ir | spectroscopy.py | IR spectrum |
| | vibrations | spectroscopy.py | Vibrational frequencies |
| | raman | spectroscopy.py | Raman spectrum |
| | mulliken | spectroscopy.py | Mulliken charges |
| | nmr_shielding | spectroscopy.py | NMR chemical shielding |
| | nmr_coupling | spectroscopy.py | NMR J-coupling |
| **TD-DFT** | tddft_states | tddft.py | Excited states |
| | electric_dipole_abs | tddft.py | Electric dipole absorption |
| | electric_dipole_soc | tddft.py | Electric dipole SOC |
| | velocity_dipole_abs | tddft.py | Velocity dipole absorption |
| | velocity_dipole_soc | tddft.py | Velocity dipole SOC |
| **Method** | method_id | method.py | Composite method identifier |
| | functional | method.py | XC functional (B3LYP, PBE0...) |
| | basis_set | method.py | Basis set (def2-TZVP...) |
| | dispersion | method.py | Dispersion (D3BJ, D4...) |
| | solvent | method.py | Solvent (water, ethanol...) |
| **Metadata** | is_optimization | energy.py | Optimization calc? |
| | optimized_state | energy.py | S0, S1, T1 |
| | calc_class | energy.py | single_point/optimization/tddft |
| | esd_type | energy.py | VG/AH/AHAS spectrum type |

---

### 🧬 ORCA Data Architecture

> **Core Principle**: ORCA does not produce files — it produces solutions to Hamiltonians under specific approximations.

The data architecture ensures scientifically correct storage while providing ergonomic access for analysis.

#### The Five Architectural Layers

```mermaid
graph TB
    subgraph "Layer 1: Molecule"
        MOL[Molecule]
        MOL --> ID[molecule_id]
        MOL --> SMILES[smiles]
        MOL --> CHG["charge / multiplicity"]
    end
    
    subgraph "Layer 2: Method"
        MTH[Method Descriptor]
        MTH --> FUNC["functional<br/>B3LYP / PBE0"]
        MTH --> BASIS["basis_set<br/>def2-TZVP"]
        MTH --> DISP["dispersion<br/>D3BJ / D4"]
        MTH --> SOL["solvent<br/>water / gas"]
    end
    
    subgraph "Layer 3: State"
        STATE[Electronic State]
        STATE --> S0[S0]
        STATE --> S1[S1]
        STATE --> T1[T1]
    end
    
    subgraph "Layer 4: Task"
        TASK[Task]
        TASK --> OPT[OPT]
        TASK --> SP[SP]
        TASK --> TDDFT[TDDFT]
    end
    
    subgraph "Layer 5: Properties"
        PROP[Properties]
        PROP --> GEO[geometry]
        PROP --> ORB[orbitals]
        PROP --> SPEC[spectra]
    end
    
    MOL --> MTH --> STATE --> TASK --> PROP
```

#### Method Descriptor (Layer 2)

A method is defined by a **composite descriptor**, not a single keyword:

| Dimension | Examples |
|-----------|----------|
| Formalism | DFT, HF, MP2, CCSD, CASSCF |
| Functional | B3LYP, ωB97X, PBE0 |
| Basis set | def2-SVP, def2-TZVP, def2-QZVP |
| Dispersion | D3BJ, D4, none |
| Relativistic | none, ZORA, DKH, X2C |
| Environment | gas, CPCM, SMD |
| Solvent | water, ethanol, acetonitrile |

> **Changing any of these creates a new method.**

#### MoleculeStore: Hierarchical Storage

```mermaid
graph TD
    subgraph "MoleculeStore"
        STORE[MoleculeStore]
        
        subgraph "p1x"
            M1[p1x]
            subgraph "B3LYP/def2-TZVP/D3BJ"
            MTH1["Method 1"]
            S0_1[S0] --> OPT1[OPT] & SP1[SP]
            S1_1[S1] --> TDDFT1[TDDFT]
            end
        end
    end
    
    STORE --> M1
    M1 --> MTH1 --> S0_1 & S1_1
```

#### Storage vs Access

| Layer | Purpose | Example |
|-------|---------|---------|
| **Storage** | Full data, all methods, reproducibility | `store._data[mol][method][state][task]` |
| **Access** | Simple queries, canonical projection | `store.get("p1x")` → best result |

#### Canonical Projection

For simple analysis, the system auto-selects the "canonical" (best) result:

- **State priority**: S0 > S1 > T1
- **Task priority**: OPT > SP > TDDFT
- **Basis priority**: def2-QZVP > def2-TZVP > def2-SVP

```python
# Simple access (uses projection)
store = MoleculeStore()
result = store.get("p1x")  # Returns canonical result

# Explicit access (for comparison)
result = store.get("p1x", method_id="DFT/B3LYP/def2-TZVP/D3BJ", state="S0")
```

#### Architectural Principles

1. **Method identity is composite** - not a single keyword
2. **Filenames are never identity** - molecule_id is extracted
3. **States are not identity** - S0 from method A ≠ S0 from method B
4. **Storage reflects physics** - all methods preserved
5. **Access reflects thinking** - simple queries return projected view
6. **Projection is mandatory** - hides complexity by default
7. **Method awareness is opt-in** - explicit only when comparing

---

### Detection & Analysis Architecture

```mermaid
graph TB
    subgraph "📊 Input"
        DF[(Parsed DataFrame)]
        CFG[Config]
    end
    
    subgraph "🔍 Hierarchy Detection"
        HD[HierarchyDetector]
        NP[NamingParser]
        PT[PatternMatcher]
        GB[GroupBuilder]
        
        subgraph "Hierarchy Output"
        RT[RootNodes]
        VR[VariantGroups]
        TR[TreeStructure]
        end
    end
    
    subgraph "📊 Partition Detection"
        PRT[PartitionDetector]
        
        subgraph "Partition Types"
        PS[StatePartition]
        PC[CalcTypePartition]
        PE[ESDPartition]
        end
        
        subgraph "Partition Output"
        S0[S0 Group]
        S1[S1 Group]
        T1[T1 Group]
        OPT[OPT Group]
        SP[SP Group]
        end
    end
    
    subgraph "🛤️ Pathway Detection"
        PWD[PathwayDetector]
        RR[ReactionRules]
        SC[StepCorrections]
        CS[ColorSchemes]
        
        subgraph "Pathway Output"
        PW[Pathways]
        ED[Edges]
        RX[Reactions]
        end
    end
    
    subgraph "⚖️ Comparison Engine"
        CE[ComparisonEngine]
        
        subgraph "Compare Types"
        CEN[EnergyCompare]
        COR[OrbitalCompare]
        CSP[SpectraCompare]
        CGE[GeometryCompare]
        end
    end
    
    subgraph "📐 Spectral Scaler"
        SS[SpectralScaler]
        LS[LinearScaler]
        RS[RelativeScaler]
        end
    
    DF & CFG --> HD
    HD --> NP --> PT --> GB
    GB --> RT & VR --> TR
    
    DF --> PRT
    PRT --> PS & PC & PE
    PS --> S0 & S1 & T1
    PC --> OPT & SP
    
    DF & TR --> PWD
    PWD --> RR & SC
    RR & SC --> PW & ED & RX
    PWD --> CS
    
    DF --> CE
    CE --> CEN & COR & CSP & CGE
    
    DF --> SS
    SS --> LS & RS
```

---

### Visualization Architecture

```mermaid
graph TB
    subgraph "📊 Data Input"
        DF[(DataFrame)]
        HR[Hierarchy]
        PT[Partitions]
        PW[Pathways]
    end
    
    subgraph "🏭 Factory"
        VF[VisualizerFactory]
        VR[VisualizerRegistry]
    end
    
    subgraph "🎨 Base"
        BV[BaseVisualizer]
        CFG[PlotConfig]
        THM[ThemeManager]
        LOG[Logger]
    end
    
    subgraph "📈 Visualizers"
        M3D[Molecule3DVisualizer]
        EDV[EnergyDiagramVisualizer]
        OPV[OrbitalPlotVisualizer]
        SPV[SpectraVisualizer]
        PWV[PathwayVisualizer]
        HTV[HierarchyTreeVisualizer]
        CMP[ComparisonVisualizer]
    end
    
    subgraph "🖼️ Renderers"
        P3M[py3Dmol Renderer]
        PL3[Plotly 3D Scatter]
        PLB[Plotly Bar]
        PLL[Plotly Line]
        PLS[Plotly Sankey]
        PLT[Plotly Treemap]
    end
    
    subgraph "📤 Output"
        FIG[Plotly Figure]
        HTM[HTML Widget]
        IMG[Image]
    end
    
    DF & HR & PT & PW --> VF
    VF --> VR
    VR --> M3D & EDV & OPV & SPV & PWV & HTV & CMP
    
    M3D & EDV & OPV & SPV & PWV & HTV & CMP -.-> BV
    BV --> CFG & THM & LOG
    
    M3D --> P3M & PL3
    EDV --> PLB
    OPV --> PLB
    SPV --> PLL
    PWV --> PLS
    HTV --> PLT
    CMP --> PLB & PLL
    
    P3M & PL3 & PLB & PLL & PLS & PLT --> FIG & HTM & IMG
```

---

### Export Architecture

```mermaid
graph TB
    subgraph "📊 Input"
        DF[(DataFrame)]
        FIG[Figures]
        MD[Metadata]
        CFG[Config]
    end
    
    subgraph "📤 Data Exporter"
        DE[DataExporter]
        JE[JSONExporter]
        CE[CSVExporter]
        PE[ParquetExporter]
        PK[PickleExporter]
    end
    
    subgraph "🌐 HTML Exporter"
        HE[HTMLExporter]
        TB[TemplateBuilder]
        PJ[PlotlyJS Embedder]
        CSS[StyleInjector]
    end
    
    subgraph "🖼️ Plot Exporter"
        PLE[PlotExporter]
        PNG[PNGExporter]
        SVG[SVGExporter]
        PDF[PDFExporter]
    end
    
    subgraph "📁 Output"
        OJ[data.json]
        OC[data.csv]
        OP[data.parquet]
        OK[data.pkl]
        OH[report.html]
        OI[plots/]
    end
    
    DF & MD --> DE
    DE --> JE --> OJ
    DE --> CE --> OC
    DE --> PE --> OP
    DE --> PK --> OK
    
    DF & FIG & MD & CFG --> HE
    HE --> TB --> PJ & CSS
    PJ & CSS --> OH
    
    FIG --> PLE
    PLE --> PNG & SVG & PDF --> OI
```

---

### Streamlit Application Architecture (Updated)

```mermaid
graph TB
    subgraph "🖥️ Entry"
        APP[streamlit_app/app.py]
    end
    
    subgraph "🧩 Components"
        EXP[ExportPanel]
        UPL[FileUploader]
        MOL[MoleculeInfo]
        TAB[VizTabs]
    end
    
    subgraph "🧠 Logic"
        LOG[LogParser]
        VIZ[VisualizerFactory]
    end
    
    APP --> UPL
    UPL --> LOG
    LOG --> VIZ
    VIZ --> TAB
    VIZ --> EXP
```

---

## 🎨 Interactive HTML Report (New)

The system now generates a **self-contained Interactive Research Paper**. 

**Key Capabilities:**
1.  **Multi-Select Comparison**: Overlay spectra (IR, Raman, UV-Vis) and compare energies for multiple selected molecules simultaneously.
2.  **3D/2D Viewer Parity**: Full interactive 3D viewer (3Dmol.js) and RDKit-generated 2D structures.
3.  **Embedded Data**: The HTML file contains **all** parsed data (XYZ coordinates, orbital energies, spectral peaks) as embedded JSON. This means the report is **fully offline-capable**—no kernel or server needed to view.
4.  **UI Parity**: Mirroring the Streamlit dashboard experience, including interactive sliders for broadening, style switching, and dark mode.

---

## 📁 Project Structure

```
Orca_Files/
├── README.md                           # This file
├── requirements.txt                    # Dependencies
├── streamlit_app/                      # Streamlit Application
│   ├── app.py                          # Main Entry Point
│   ├── utils/                          # UI Utilities
│   └── components/                     # UI Components
│       ├── export_panel.py             # HTML Export logic
│       ├── file_uploader.py            # Upload widget
│       └── ...
│
├── src/                                # Core Library
│   ├── parser/                         # Modular Parsers (Geom, Energy, etc.)
│   ├── analysis/                       # Analysis Logic (Comparison, Pathways)
│   ├── viz/                            # Visualization Logic (Plotly, 3Dmol)
│   ├── export/                         # Data Exporter Logic
│   └── core/                           # Base Classes
│
├── tests/                              # Test Suite
├── notebooks/                          # Usage Demos
└── orca_praser.py                      # Original Single-file Parser
```

---

## 🔬 API Reference

### Parser

```python
from src.parser.factory import ParserFactory

factory = ParserFactory()
result = factory.parse("molecule.out")

# Access parsed data
coords = result.geometry.cart_coords
energy = result.energy.gibbs_Eh
orbitals = result.orbitals.homo_lumo
```

### Hierarchy Detection

```python
from src.analysis.hierarchy_detector import HierarchyDetector

detector = HierarchyDetector(df)
hierarchy = detector.detect()

# p1x, p1a, p1b → p1 (root) with variants x, a, b
print(hierarchy.to_tree())
```

### Partition Detection

```python
from src.analysis.partition_detector import PartitionDetector

detector = PartitionDetector(df)
partitions = detector.detect()

# {"by_state": {"S0": [...], "S1": [...]}, 
#  "by_calc_type": {"OPT": [...], "SP": [...]}}
```

### Pathway Detection

```python
from src.analysis.pathway_detector import PathwayDetector

detector = PathwayDetector(df, hierarchy)
detector.set_reaction_rules({
    ("p1", "p2"): {"add": {"OH": 4}, "remove": {"H2O": 3}},
    ("p2", "p3"): {"add": {"OH": 2}, "remove": {"H2O": 1}},
})
detector.set_step_corrections({
    ("p1x", "p1a"): {"add": {"OH": 2}, "remove": {"H2O": 1}},
})
detector.set_color_scheme("by_variant")
pathways = detector.detect()
```

### Spectral Scaling

```python
from src.analysis.spectral_scaler import SpectralScaler

scaler = SpectralScaler(spectrum_df)

# Linear: ν_s = s × ν
linear = scaler.linear_scale(factor=0.97)

# Relative: ν_s = ν_min + s × (ν - ν_min)
relative = scaler.relative_scale(factor=1.5)
```

### Data Export

```python
from src.export.data_exporter import DataExporter

exporter = DataExporter(df, metadata=metadata)
exporter.to_json("data.json")
exporter.to_csv("data.csv")
exporter.to_parquet("data.parquet")
exporter.export_bundle("results/")  # All formats + metadata
```

### HTML Export

```python
from src.export.html_exporter import HTMLExporter

exporter = HTMLExporter(df)
exporter.add_molecule_3d("p1x")
exporter.add_energy_diagram(["p1x", "p2x", "p3x"])
exporter.add_pathway_diagram(pathways)
exporter.add_spectra("p1x", spectrum_type="ir")
exporter.export("report.html")
```

---

## 🪵 Logging

```python
# src/logger.py configuration
LOG_FORMAT = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'

# Levels: DEBUG, INFO, WARNING, ERROR
# Output: Console + orca_viz.log file
```

---

## 🚀 Quick Start

```bash
# Install
pip install -r requirements.txt

# Run Streamlit app (New Entry Point)
streamlit run streamlit_app/app.py

# Parse with modular parser (exports 16 CSVs)
python tests/test_comprehensive.py

# Parse with original parser (for comparison)
python tests/test_original_parser.py

# Compare both parsers
python tests/test_comparison.py

# Generate visualization HTML
python tests/test_visualizations.py

# Run unit tests
pytest tests/ -v -s
```

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jauharmz/Orca_Files/blob/main/ORCA_Test_v2.ipynb)

---

## 🤗 HuggingFace

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="JauharMz/Orca",
    repo_type="dataset",
    local_dir="./data"
)
```

---

*Last Updated: 2026-01-04*
