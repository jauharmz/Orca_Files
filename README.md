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

**Currently Parsed (28 sections):**
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
16. **Loewdin Bond Orders:** Alternative bond order analysis
17. **Thermochemistry:**
    - Zero-point energy (ZPE)
    - Thermal corrections
    - Enthalpy, Entropy (electronic/vibrational/rotational/translational)
    - Gibbs free energy
18. **NMR Chemical Shifts:** Isotropic shielding (ppm)
19. **NMR J-Couplings:** Isotropic coupling constants (Hz)
20. **Normal Modes:** Vibrational displacement vectors (partial)
21. **SCF Iterations:** Detailed convergence data per iteration
22. **Timing Data:** Computational timing breakdown
23. **DFT Grid Info:** Integration grid parameters
24. **Basis Set Info:** Basis set name, functions, primitives
25. **Energy Components:** Nuclear/electronic/kinetic/potential/virial/XC
26. **CPCM Solvation:** Surface charge, dielectric energy
27. **SCF Convergence:** Final convergence metrics
28. **Mulliken Orbital Populations:** s,p,d,f,g breakdown per atom

**Unparsed Sections Available (29 more):**

**HIGH PRIORITY (10 sections):**
- Mulliken/Loewdin Orbital Populations Per MO (~11k MOs - HUGE dataset)
- Mulliken Overlap Charges
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

**Preview:**
- Energy convergence plots
- Orbital energy diagrams
- IR/Raman spectrum plots
- Atomic charge bar charts
- NMR shift tables
- Thermochemistry summary
- Orbital population analysis

---

## Detailed ORCA Sections Analysis

**Last Updated:** 2025-11-25
**Test File:** p1xs0p.out (113,234 lines)
**Current Coverage:** 28/57 sections parsed (49%)

### Status Summary

| Category | Parsed | Remaining | Total |
|----------|--------|-----------|-------|
| **HIGH Priority** | 6 | 10 | 16 |
| **MEDIUM Priority** | 15 | 12 | 27 |
| **LOW Priority** | 7 | 7 | 14 |
| **TOTAL** | **28** | **29** | **57** |

### High Priority Unparsed Sections (10)

#### Electronic Structure - Charge Analysis (4 sections)

**28. MULLIKEN OVERLAP CHARGES** ⭐
- **Location:** Line 75366
- **Data:** Charge overlap between atom pairs
- **Size:** 23×23 = 529 values
- **Use Case:** Alternative bonding analysis
- **Estimated Effort:** 2-3 hours

**29. MULLIKEN ORBITAL CHARGES**
- **Location:** Line 74800
- **Data:** Per-orbital charge distribution
- **Size:** ~23 atoms × ~10 orbital shells
- **Estimated Effort:** 2-3 hours

**30. LOEWDIN ORBITAL CHARGES**
- **Location:** Line 90021
- **Data:** Loewdin version of orbital charges
- **Estimated Effort:** 2-3 hours

**31. LOEWDIN BOND ORDERS** ⭐⭐
- **Location:** Line 90585
- **Data:** Loewdin-based bond strengths (threshold 0.05)
- **Size:** ~66 bonds
- **Use Case:** Compare with Mayer bond orders
- **Estimated Effort:** 1-2 hours

**32-34. MO POPULATIONS PER MO** ⚠️ HUGE (3 sections)
- **Locations:** Lines 75406, 90611, 103939
- **Size:** ~11,609 MOs × contributions = **~580k entries**
- **Memory:** ~50-100 MB for full dataset
- **Strategy:** Use `detailed=True` flag, sparse storage
- **Estimated Effort:** 6-8 hours

#### NMR - Tensor Components (3 sections)

**35. J-COUPLING TENSOR COMPONENTS** 🌟 TOP PRIORITY
- **Location:** Line 3109+
- **Current:** Only parse isotropic J-value
- **Missing Components:**
  - DSO (Diamagnetic Spin-Orbit) - 3×3 tensor
  - PSO (Paramagnetic Spin-Orbit) - 3×3 tensor
  - FC (Fermi Contact) - 3×3 tensor
  - SD (Spin-Dipolar) - 3×3 tensor
  - SD/FC (Cross term) - 3×3 tensor
  - Total J tensor + eigenvalues
- **Size:** 15 atom pairs × 6 tensors = 810 values
- **Use Case:** Advanced NMR analysis, mechanism studies
- **Estimated Effort:** 4-6 hours

**36. CHEMICAL SHIELDING TENSORS** 🌟 TOP PRIORITY
- **Location:** Line 2491
- **Current:** Only parse isotropic shielding
- **Missing:** Full anisotropic tensor + principal components
- **Size:** 18 nuclei × 9 tensor elements = 162 values
- **Use Case:** Anisotropic NMR, solid-state NMR
- **Estimated Effort:** 3-4 hours

**37. CHEMICAL SHIELDING SUMMARY**
- **Location:** Line 3079
- **Note:** May be redundant with tensor parsing
- **Estimated Effort:** 1-2 hours

### Medium Priority Unparsed Sections (12)

#### Geometry & Coordinates (3 sections)

**38. INTERNAL COORDINATES** ⭐
- **Location:** Lines 375, 402
- **Data:** Bond lengths, angles, dihedrals (Z-matrix)
- **Size:** ~100 internal coordinates
- **Estimated Effort:** 2-3 hours

**39. CARTESIAN COORDINATES (A.U.)**
- **Location:** Line 347
- **Data:** Coordinates in atomic units (Bohr)
- **Estimated Effort:** 30 min

**40. GEOMETRIC PERTURBATIONS**
- **Location:** Line 111531
- **Data:** Numerical derivatives for gradients
- **Size:** 69 perturbations (23 atoms × 3 dims)
- **Estimated Effort:** 2 hours

#### Basis Set & Integrals (3 sections)

**41. BASIS SET INFORMATION (Enhanced)**
- **Location:** Line 429
- **Missing:** Contraction schemes, exponents, coefficients
- **Estimated Effort:** 4-5 hours

**42. SHARK INTEGRAL PACKAGE**
- **Location:** Line 526
- **Data:** Integral method, screening thresholds
- **Estimated Effort:** 1-2 hours

**43. COSX GRID GENERATION**
- **Location:** Line 618
- **Data:** Chain-of-spheres exchange grid parameters
- **Estimated Effort:** 1-2 hours

#### SCF Details (4 sections)

**44. INITIAL GUESS: MOREAD**
- **Location:** Line 795
- **Estimated Effort:** 30 min

**45. INITIAL GUESS ORBITALS**
- **Location:** Line 4673
- **Estimated Effort:** 1 hour

**46. DIIS DETAILS**
- **Location:** Line 854
- **Estimated Effort:** 2-3 hours

**47. SOSCF DETAILS**
- **Location:** Line 864
- **Estimated Effort:** 2-3 hours

#### Solvation & Grid (2 sections)

**48. CPCM SOLVATION (Enhanced)**
- **Location:** Line 828
- **Missing:** Cavity details, surface area, radii
- **Estimated Effort:** 2-3 hours

**49. DFT GRID GENERATION (Enhanced)**
- **Location:** Line 617
- **Missing:** Construction details, partition functions
- **Estimated Effort:** 2 hours

### Low Priority Unparsed Sections (7)

**50. NORMAL MODES (Complete)** - All 63 modes (3-4 hours)
**51. SCF HESSIAN** ⚠️ LARGE - 69×69 matrix (4-5 hours)
**52. MO COEFFICIENTS** ⚠️ MASSIVE - 15M values (6-8 hours)
**53. POPLE LINEAR SOLVER** - CPKS/CPHF details (2-3 hours)
**54. SCF SETTINGS (Detailed)** - All parameters (3-4 hours)
**55. CITATIONS** - Reference existing bibtex parser (1 hour)
**56. TOTAL RUN TIME** - Summary only (15 min)

### Implementation Roadmap

**Phase 1: High-Value NMR (8-10 hours)**
1. J-Coupling Tensor Components (4-6 hrs)
2. Chemical Shielding Tensors (3-4 hrs)

**Phase 2: Charge & Bond Analysis (5-7 hours)**
3. Loewdin Bond Orders (1-2 hrs)
4. Mulliken Overlap Charges (2-3 hrs)
5. Mulliken/Loewdin Orbital Charges (2-3 hrs)

**Phase 3: Geometry Details (3-4 hours)**
6. Internal Coordinates (2-3 hrs)
7. Cartesian Coordinates (a.u.) (30 min)

**Phase 4: Computational Details (6-8 hours)**
8. Enhanced Basis Set Info (4-5 hrs)
9. Enhanced CPCM Solvation (2-3 hrs)
10. SHARK Integral Package (1-2 hrs)

**Phase 5: Advanced/Optional (15-20 hours)**
11. MO Populations Per MO (6-8 hrs)
12. Complete Normal Modes (3-4 hrs)
13. SCF Hessian (4-5 hrs)

**Total Estimated Effort:** 37-49 hours for all high/medium priority sections

### Storage Considerations

- **Current parsing:** ~500 KB JSON
- **With HIGH priority:** ~2-3 MB JSON
- **With MO populations:** ~50-100 MB JSON
- **With MO coefficients:** ~150-200 MB JSON

### Implementation Notes

**Parser Architecture:**
- Standalone functions returning None for missing sections
- Optional `detailed=True` parameter for large datasets
- Maintain backwards compatibility

**Dataclass Design:**
- Use `@dataclass` with `to_dict()` methods
- `Optional[T]` for nullable fields
- `field(default_factory=list)` for mutable defaults

**Testing Strategy:**
- Manual verification (2-3 samples per section)
- Unit tests for each parser
- Integration tests with JSON export
- Performance monitoring (<10 sec for 113k lines)

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

