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

### 3. Output Files (`.out`)

**Format:** ORCA text output with multiple sections

**Currently Parsed (27 sections):**
1. **Job Info:** Method, basis set, charge, multiplicity
2. **Final Energy:** Total SCF energy
3. **SCF Energies:** Energy at each SCF iteration
4. **Optimization Energies:** Energy at each geometry step
5. **Coordinates:** Cartesian coordinates (Angstrom)
6. **Dipole Moment:** X, Y, Z components (a.u. and Debye)
7. **Polarizability:** Static polarizability tensor + eigenvalues
8. **Orbital Energies:** HOMO, LUMO, all occupied/virtual orbitals
9. **Vibrational Frequencies:** All normal modes (cm⁻¹)
10. **IR Spectrum:** Frequencies with intensities (km/mol)
11. **Raman Spectrum:** Frequencies with activities and depolarization
12. **Dispersion Correction:** DFTD3 (E6, E8 components)
13. **Mulliken Charges:** Atomic partial charges
14. **Loewdin Charges:** Alternative charge analysis
15. **Mayer Bond Orders:** Covalent bond strengths
16. **Thermochemistry:**
    - Zero-point energy (ZPE)
    - Thermal corrections
    - Enthalpy, Entropy (electronic/vibrational/rotational/translational)
    - Gibbs free energy
17. **NMR Chemical Shifts:** Isotropic shielding (ppm)
18. **NMR J-Couplings:** Isotropic coupling constants (Hz)
19. **Normal Modes:** Vibrational displacement vectors (partial)
20. **SCF Iterations:** Detailed convergence data per iteration
21. **Timing Data:** Computational timing breakdown
22. **DFT Grid Info:** Integration grid parameters
23. **Basis Set Info:** Basis set name, functions, primitives
24. **Energy Components:** Nuclear/electronic/kinetic/potential/virial/XC
25. **CPCM Solvation:** Surface charge, dielectric energy
26. **SCF Convergence:** Final convergence metrics
27. **Mulliken Orbital Populations:** s,p,d,f,g breakdown per atom

**Unparsed Sections Available (30 more):**

**HIGH PRIORITY (11 sections):**
- Mulliken/Loewdin Orbital Populations Per MO (~11k MOs - HUGE dataset)
- Mulliken Overlap Charges
- Loewdin Bond Orders
- Mulliken/Loewdin Orbital Charges
- J-Coupling Tensor Components (DSO/PSO/FC/SD/SD-FC)
- Chemical Shielding Tensors (full anisotropic)
- Chemical Shielding Summary

**MEDIUM PRIORITY (12 sections):**
- Internal Coordinates (bond lengths/angles/dihedrals)
- Cartesian Coordinates (atomic units)
- Geometric Perturbations
- Basis Set Details (contractions, exponents)
- SHARK Integral Package
- COSX Grid Generation
- Initial Guess Orbitals
- DIIS/SOSCF Details
- CPCM Detailed Parameters

**LOW PRIORITY (7 sections):**
- Normal Modes (complete all 63 modes)
- SCF Hessian Matrix
- MO Coefficients (MASSIVE dataset)
- Pople Linear Equation Solver
- SCF Settings (detailed)
- Citations
- Total Runtime

See `ORCA_SECTIONS_ROADMAP.md` for complete analysis.

**Preview:**
- Energy convergence plots
- Orbital energy diagrams
- IR/Raman spectrum plots
- Atomic charge bar charts
- NMR shift tables
- Thermochemistry summary
- Orbital population analysis

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

### Binary Files (Not Supported)
- `.opt`, `.densitiesinfo`, `.cis`, `.ges` - Proprietary ORCA binary format without public documentation

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

