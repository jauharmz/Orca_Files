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
| **Multi-Comparison** | Compare multiple molecules side-by-side |
| **Spectral Scaling** | Linear (`ν_s = s × ν`) and relative (`ν_s = ν_min + s(ν - ν_min)`) scaling |
| **Data Export** | Export parsed data to JSON, CSV, Parquet, Pickle |
| **HTML Export** | Single-file interactive reports with embedded Plotly.js |
| **Interactive Viz** | Plotly charts + py3Dmol 3D molecular viewer |

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
        APP[App Entry]
        SM[SessionManager]
        
        subgraph "Pages"
            PH[HomePage]
            PU[UploadPage]
            PMol[MoleculePage]
            PCmp[ComparePage]
            PPW[PathwayPage]
            PExp[ExportPage]
            PSet[SettingsPage]
        end
        
        subgraph "Components"
            CFU[FileUploader]
            CMS[MoleculeSelector]
            CVT[VisualizationTabs]
            CEP[ExportPanel]
            CLV[LogViewer]
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
    APP --> SM
    SM --> PH & PU & PMol & PCmp & PPW & PExp & PSet
    PH & PU --> CFU & CMS
    PMol & PCmp --> CVT
    PExp --> CEP
    PSet --> CLV
    
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

### Streamlit Application Architecture

```mermaid
graph TB
    subgraph "🖥️ Entry"
        APP[app.py]
        CFG[config.py]
    end
    
    subgraph "📦 State Management"
        SM[SessionManager]
        SS[SessionState]
        DC[DataCache]
        VC[ViewCache]
    end
    
    subgraph "🔧 Services"
        PS[ParserService]
        AS[AnalysisService]
        VS[VizService]
        ES[ExportService]
        LS[LogService]
    end
    
    subgraph "📄 Pages"
        PH[HomePage]
        PU[UploadPage]
        PMol[MoleculePage]
        PCmp[ComparePage]
        PPW[PathwayPage]
        PSpc[SpectraPage]
        PExp[ExportPage]
        PSet[SettingsPage]
    end
    
    subgraph "🧩 Components"
        CFU[FileUploader]
        CMS[MoleculeSelector]
        CPF[PartitionFilter]
        CVT[VisualizationTabs]
        CPW[PathwayEditor]
        CSS[ScalerSliders]
        CEP[ExportPanel]
        CLV[LogViewer]
    end
    
    APP --> CFG --> SM
    SM --> SS & DC & VC
    
    SM --> PS & AS & VS & ES & LS
    
    APP --> PH & PU & PMol & PCmp & PPW & PSpc & PExp & PSet
    
    PU --> CFU
    PMol --> CMS & CVT
    PCmp --> CMS & CPF & CVT
    PPW --> CPW & CVT
    PSpc --> CSS & CVT
    PExp --> CEP
    PSet --> CLV
```

---

### Data Flow Sequence

```mermaid
sequenceDiagram
    participant U as User
    participant ST as Streamlit
    participant FH as FileHandler
    participant PF as ParserFactory
    participant HD as HierarchyDetector
    participant PRT as PartitionDetector
    participant PWD as PathwayDetector
    participant VF as VizFactory
    participant EX as Exporter
    participant LOG as Logger
    
    U->>ST: Upload .out files
    ST->>FH: validate_files(files)
    FH->>LOG: log("Validating...")
    FH-->>ST: valid_files[]
    
    loop Each file
        ST->>PF: parse(file)
        PF->>LOG: log("Parsing geometry...")
        PF->>LOG: log("Parsing energy...")
        PF->>LOG: log("Parsing orbitals...")
        PF->>LOG: log("Parsing spectra...")
        PF-->>ST: ParseResult
    end
    
    ST->>HD: detect(df)
    HD->>LOG: log("Detecting hierarchy...")
    HD-->>ST: Hierarchy
    
    ST->>PRT: detect(df)
    PRT->>LOG: log("Detecting partitions...")
    PRT-->>ST: Partitions
    
    ST->>PWD: detect(df, hierarchy)
    PWD->>LOG: log("Detecting pathways...")
    PWD-->>ST: Pathways
    
    U->>ST: Select molecule
    ST->>VF: create_visualizations(mol_id)
    VF->>LOG: log("Creating 3D view...")
    VF->>LOG: log("Creating energy plot...")
    VF-->>ST: [Figure, ...]
    ST-->>U: Display visualizations
    
    U->>ST: Export
    ST->>EX: export(format)
    EX->>LOG: log("Exporting...")
    EX-->>ST: file_path
    ST-->>U: Download link
```

---

## 📁 Project Structure

```
Orca_Files/
├── README.md                           # This file
├── requirements.txt                    # Dependencies
├── app.py                              # Streamlit entry
│
├── src/
│   ├── __init__.py
│   ├── config.py                       # Global config
│   ├── logger.py                       # Logging setup
│   │
│   ├── core/                           # Core abstractions
│   │   ├── __init__.py
│   │   ├── base_parser.py              # BaseParser ABC
│   │   ├── base_visualizer.py          # BaseVisualizer ABC
│   │   ├── data_models.py              # Pydantic models
│   │   ├── exceptions.py               # Custom exceptions
│   │   └── registry.py                 # Registry pattern
│   │
│   ├── parser/                         # Parser modules (refactored from orca_praser.py)
│   │   ├── __init__.py                 # Package exports
│   │   ├── factory.py                  # ParserFactory - orchestrates all parsers
│   │   ├── batch.py                    # BatchParser - multi-file with molecule grouping
│   │   ├── geometry.py                 # Coords, SMILES, internal coords (B/A/D)
│   │   ├── energy.py                   # Gibbs, SP energy, ESD flag, calc type
│   │   ├── orbitals.py                 # HOMO/LUMO, spin handling, orbital levels
│   │   ├── spectroscopy.py             # IR, Raman, Mulliken, NMR
│   │   ├── tddft.py                    # TD-DFT states, dipole spectra, H→L labels
│   │   ├── spectrum_file.py            # External spectrum files (VG/AH/AHAS/FLUOR)
│   │   └── regex_patterns.py           # Shared regex patterns
│   │
│   ├── analysis/                       # Analysis modules
│   │   ├── __init__.py
│   │   ├── hierarchy_detector.py       # Root/variant groups
│   │   ├── partition_detector.py       # S0/S1, OPT/SP
│   │   ├── pathway_detector.py         # Degradation paths
│   │   ├── reaction_rules.py           # Add/remove rules
│   │   ├── comparison_engine.py        # Multi-compare
│   │   └── spectral_scaler.py          # Linear/relative
│   │
│   ├── viz/                            # Visualization modules
│   │   ├── __init__.py
│   │   ├── factory.py                  # VisualizerFactory
│   │   ├── config.py                   # Plot themes
│   │   ├── molecule_3d.py              # py3Dmol + Plotly 3D
│   │   ├── energy_diagram.py           # Energy bars
│   │   ├── orbital_plot.py             # Orbital levels
│   │   ├── spectra_plot.py             # IR/Raman/UV-Vis
│   │   ├── pathway_viz.py              # Pathway network
│   │   ├── hierarchy_tree.py           # Tree view
│   │   └── comparison_viz.py           # Side-by-side
│   │
│   ├── export/                         # Export modules
│   │   ├── __init__.py
│   │   ├── data_exporter.py            # JSON/CSV/Parquet
│   │   ├── html_exporter.py            # Interactive HTML
│   │   ├── plot_exporter.py            # PNG/SVG
│   │   └── templates/
│   │       └── report.html             # HTML template
│   │
│   ├── services/                       # Business logic
│   │   ├── __init__.py
│   │   ├── parser_service.py           # Parsing orchestration
│   │   ├── analysis_service.py         # Analysis orchestration
│   │   ├── viz_service.py              # Viz orchestration
│   │   └── export_service.py           # Export orchestration
│   │
│   ├── ui/                             # Streamlit UI
│   │   ├── __init__.py
│   │   ├── session_manager.py          # State management
│   │   ├── pages/
│   │   │   ├── __init__.py
│   │   │   ├── home.py
│   │   │   ├── upload.py
│   │   │   ├── molecule.py
│   │   │   ├── compare.py
│   │   │   ├── pathway.py
│   │   │   ├── spectra.py
│   │   │   ├── export.py
│   │   │   └── settings.py
│   │   └── components/
│   │       ├── __init__.py
│   │       ├── file_uploader.py
│   │       ├── molecule_selector.py
│   │       ├── partition_filter.py
│   │       ├── viz_tabs.py
│   │       ├── pathway_editor.py
│   │       ├── scaler_sliders.py
│   │       ├── export_panel.py
│   │       └── log_viewer.py
│   │
│   └── utils/                          # Utilities
│       ├── __init__.py
│       ├── converters.py               # Unit conversions
│       ├── validators.py               # Input validation
│       └── helpers.py                  # General helpers
│
├── tests/                              # Test suite
│   ├── __init__.py
│   ├── README.md                       # Test documentation
│   ├── test_comprehensive.py           # Full parser test → 16 CSV exports
│   ├── test_original_parser.py         # Original orca_praser.py test
│   ├── test_comparison.py              # Side-by-side parser comparison
│   ├── test_visualizations.py          # Visualization tests → HTML plots
│   ├── test_parsers.py                 # Unit tests for parsers
│   ├── test_analysis.py                # Analysis module tests
│   └── test_real.py                    # Integration tests
│
├── notebooks/                          # Demo notebooks
│   ├── 01_parser_demo.ipynb
│   ├── 02_analysis_demo.ipynb
│   └── 03_viz_demo.ipynb
│
├── orca_praser.py                      # Original reference parser
└── requirements.txt                    # Dependencies
```

### CSV Export Structure

When running `python tests/test_comprehensive.py`, up to 16 CSV files are generated:

| File | Contents | Rows |
|------|----------|------|
| `parsed_molecules.csv` | Scalar data + counts | 1 per molecule |
| `parsed_cart_coords.csv` | x, y, z coordinates | All atoms |
| `parsed_orbitals.csv` | OCC, Eh, eV, spin, lvl | All orbitals |
| `parsed_tddft_states.csv` | Excited state transitions | All states |
| `parsed_ir_spectra.csv` | freq, eps, intensity | All IR peaks |
| `parsed_vibrations.csv` | freq, imaginary flag | All modes |
| `parsed_raman.csv` | freq, activity, depolarization | All peaks |
| `parsed_mulliken.csv` | Nucleus, Element, Pop, Charge | All atoms |
| `parsed_nmr_shielding.csv` | Isotropic, Anisotropy | All nuclei |
| `parsed_nmr_coupling.csv` | Nucleus1, Nucleus2, J_Hz | All couplings |
| `parsed_internal_bonds.csv` | Bond definitions, values | All bonds |
| `parsed_internal_angles.csv` | Angle definitions, values | All angles |
| `parsed_internal_dihedrals.csv` | Dihedral definitions, values | All dihedrals |
| `parsed_electric_dipole.csv` | Absorption spectrum | All transitions |
| `parsed_velocity_dipole.csv` | Velocity dipole spectrum | All transitions |
| `parsed_data.json` | Complete nested data | All data |

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

# Run Streamlit app
streamlit run app.py
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

*Last Updated: 2026-01-03*
