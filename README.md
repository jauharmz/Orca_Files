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
| `.inp` | Text | ORCA input file (calculation settings) | ✅ Done |
| `.out` | Text | Main ORCA output (energies, SCF, etc.) | ✅ Done |
| `.hess` | Text | Hessian matrix (vibrational analysis) | ✅ Done |
| `.property.txt` | Text | Computed properties | ✅ Done |
| `.engrad` | Text | Energy and gradient | ✅ Done |
| `.opt` | Binary | Geometry optimization trajectory | ⚠️ Binary |
| `.cpcm` | Text | CPCM solvation model data | ✅ Done |
| `.cpcm_corr` | Text | CPCM corrections | ✅ Done |
| `.densitiesinfo` | Binary | Electron density information | ⚠️ Binary |
| `.bibtex` | BibTeX | Citation references | ✅ Done |
| `_trj.xyz` | XYZ | Optimization trajectory | ✅ Done |
| `.spectrum` | Text | Spectrum data (UV-Vis, IR, etc.) | ✅ Done |
| `.cis` | Binary | CI Singles data | ⚠️ Binary |
| `.ges` | Binary | Ground/Excited state data | ⚠️ Binary |

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
- [x] **1.2** Parse `.inp` files → Extract calculation settings
- [x] **1.3** Parse `.out` files → Extract energies, SCF convergence
- [x] **1.4** Parse `.hess` files → Extract frequencies, normal modes
- [x] **1.5** Parse `.property.txt` → Extract computed properties
- [x] **1.6** Parse `.engrad` files → Extract gradients
- [x] **1.8** Parse `.spectrum` files → Extract spectral data
- [x] **1.11** Parse `_trj.xyz` → Extract trajectory frames

### Binary Files (Not Implemented)
- `.opt`, `.densitiesinfo`, `.cis`, `.ges` - Binary format, needs special handling

### Phase 2: Data Preview Components
- [x] **2.1** 3D Molecular viewer (using py3Dmol)
- [x] **2.2** Energy plot viewer
- [x] **2.3** Spectrum plot viewer (UV-Vis, IR)
- [x] **2.4** Table viewer for properties

### Phase 3: Web UI Integration
- [x] **3.1** File upload interface
- [x] **3.2** Automatic file type detection
- [x] **3.3** Unified dashboard for all previews
- [x] **3.4** Export parsed data (JSON)

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

