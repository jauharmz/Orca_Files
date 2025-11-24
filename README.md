# ORCA Output Parser & Viewer

A Python-based tool to parse and visualize ORCA quantum chemistry output files.

## Project Overview

This project provides parsers for various ORCA output file formats and a UI to preview the parsed data.

---

## File Format Documentation

### Analyzed Files

| Extension | Format | Description | Status |
|-----------|--------|-------------|--------|
| `.xyz` | XYZ | Molecular geometry coordinates | ✅ Done |
| `.inp` | Text | ORCA input file (calculation settings) | ⏳ Pending |
| `.out` | Text | Main ORCA output (energies, SCF, etc.) | ✅ Done |
| `.hess` | Text | Hessian matrix (vibrational analysis) | ✅ Done |
| `.property.txt` | Text | Computed properties | ⏳ Pending |
| `.engrad` | Text | Energy and gradient | ⏳ Pending |
| `.opt` | Text | Geometry optimization trajectory | ⏳ Pending |
| `.cpcm` | Text | CPCM solvation model data | ⏳ Pending |
| `.cpcm_corr` | Text | CPCM corrections | ⏳ Pending |
| `.densitiesinfo` | Text | Electron density information | ⏳ Pending |
| `.bibtex` | BibTeX | Citation references | ⏳ Pending |
| `_trj.xyz` | XYZ | Optimization trajectory | ✅ Done |
| `.spectrum` | Text | Spectrum data (UV-Vis, IR, etc.) | ✅ Done |
| `.cis` | Binary/Text | CI Singles data | ⏳ Pending |
| `.ges` | Binary/Text | Ground/Excited state data | ⏳ Pending |

---

## Detailed File Descriptions

### 1. XYZ Files (`.xyz`, `_trj.xyz`)

**Format:** Standard XYZ molecular geometry format
- Line 1: Number of atoms
- Line 2: Comment/title (energy info)
- Lines 3+: `Element X Y Z` coordinates in Angstroms

**Data Contains:**
- Atomic symbols
- Cartesian coordinates (X, Y, Z)
- Total energy (in comment line)

**Preview:** 3D molecular structure viewer

---

### 2. Input Files (`.inp`)

**Format:** ORCA input format
- `!` line: Keywords (method, basis set, job type)
- `%` blocks: Additional settings
- `* xyz` block: Molecular geometry

**Data Contains:**
- Calculation method (e.g., B3LYP)
- Basis set (e.g., 6-311++G(d,p))
- Job types (OPT, FREQ, etc.)
- Solvent model
- Molecular geometry

**Preview:** Formatted text with syntax highlighting

---

## Roadmap

### Phase 1: Core Parsers (Current)
- [x] **1.1** Parse `.xyz` files → Extract atoms & coordinates
- [ ] **1.2** Parse `.inp` files → Extract calculation settings
- [x] **1.3** Parse `.out` files → Extract energies, SCF convergence
- [x] **1.4** Parse `.hess` files → Extract frequencies, normal modes
- [ ] **1.5** Parse `.property.txt` → Extract computed properties
- [ ] **1.6** Parse `.engrad` files → Extract gradients
- [ ] **1.7** Parse `.opt` files → Extract optimization steps
- [x] **1.8** Parse `.spectrum` files → Extract spectral data
- [ ] **1.9** Parse `.cpcm` / `.cpcm_corr` → Extract solvation data
- [ ] **1.10** Parse `.densitiesinfo` → Extract density info
- [x] **1.11** Parse `_trj.xyz` → Extract trajectory frames

### Phase 2: Data Preview Components
- [ ] **2.1** 3D Molecular viewer (using py3Dmol or similar)
- [ ] **2.2** Energy plot viewer
- [ ] **2.3** Spectrum plot viewer (UV-Vis, IR)
- [ ] **2.4** Table viewer for properties
- [ ] **2.5** Optimization trajectory animation

### Phase 3: Web UI Integration
- [ ] **3.1** File upload interface
- [ ] **3.2** Automatic file type detection
- [ ] **3.3** Unified dashboard for all previews
- [ ] **3.4** Export parsed data (JSON, CSV)

---

## Tech Stack

- **Backend:** Python 3.10+
- **Web Framework:** Flask or Streamlit (for quick prototyping)
- **Parsing:** Custom parsers with regex
- **Visualization:**
  - Plotly (for spectra and energy plots)
  - py3Dmol (for 3D molecular structures)
- **Testing:** pytest
- **Frontend:** Streamlit (single-page app) or Flask + HTML/JS

---

## Project Structure

```
Orca_Files/
├── README.md
├── requirements.txt
├── parsers/
│   ├── __init__.py
│   ├── xyz_parser.py
│   ├── inp_parser.py
│   ├── out_parser.py
│   ├── hess_parser.py
│   └── ...
├── previews/
│   ├── __init__.py
│   ├── molecule_viewer.py
│   ├── spectrum_viewer.py
│   └── ...
├── tests/
│   ├── test_xyz_parser.py
│   ├── test_inp_parser.py
│   └── ...
├── app.py              # Main Streamlit/Flask app
└── sample_data/        # Example ORCA files (current files)
```

---

## Getting Started

```bash
# Install dependencies
pip install -r requirements.txt

# Run tests
pytest tests/

# Launch preview app
streamlit run app.py
```

---

## Current Progress

**Working on:** File format analysis and documentation

**Next step:** Implement `.xyz` parser with tests and preview

