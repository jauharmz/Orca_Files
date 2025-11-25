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

**Currently Parsed (32 sections):**
1. **Job Info:** Method, basis set, charge, multiplicity
2. **Final Energy:** Total SCF energy
3. **SCF Energies:** Energy at each SCF iteration
4. **Optimization Energies:** Energy at each geometry step
5. **Coordinates:** Cartesian coordinates (Angstrom)
6. **Coordinates (a.u.):** Cartesian coordinates with atomic numbers and masses (Bohr)
7. **Internal Coordinates:** Z-matrix (bond lengths, angles, dihedrals) in Angstrom and a.u.
8. **Dipole Moment:** X, Y, Z components (a.u. and Debye)
9. **Polarizability:** Static polarizability tensor + eigenvalues
10. **Orbital Energies:** HOMO, LUMO, all occupied/virtual orbitals
11. **Vibrational Frequencies:** All normal modes (cm⁻¹)
12. **IR Spectrum:** Frequencies with intensities (km/mol)
13. **Raman Spectrum:** Frequencies with activities and depolarization
14. **Dispersion Correction:** DFTD3 (E6, E8 components)
15. **Mulliken Charges:** Atomic partial charges
16. **Mulliken Overlap Charges:** Charge overlap between atom pairs (105 pairs)
17. **Loewdin Charges:** Alternative charge analysis
18. **Mayer Bond Orders:** Covalent bond strengths
19. **Loewdin Bond Orders:** Alternative bond order analysis
20. **Thermochemistry:**
    - Zero-point energy (ZPE)
    - Thermal corrections
    - Enthalpy, Entropy (electronic/vibrational/rotational/translational)
    - Gibbs free energy
21. **NMR Chemical Shifts:** Isotropic shielding (ppm)
22. **NMR J-Couplings:** Isotropic coupling constants (Hz)
23. **Normal Modes:** Vibrational displacement vectors (partial)
24. **SCF Iterations:** Detailed convergence data per iteration
25. **Timing Data:** Computational timing breakdown + total run time
26. **DFT Grid Info:** Integration grid parameters
27. **Basis Set Info:** Basis set name, functions, primitives
28. **Energy Components:** Nuclear/electronic/kinetic/potential/virial/XC
29. **CPCM Solvation:** Surface charge, dielectric energy
30. **SCF Convergence:** Final convergence metrics
31. **Mulliken Orbital Populations:** s,p,d,f,g breakdown per atom
32. **Total Run Time:** Complete execution time (days, hours, minutes, seconds)

**Unparsed Sections Available (25 more):**

**HIGH PRIORITY (9 sections):**
- Mulliken/Loewdin Orbital Populations Per MO (~11k MOs - HUGE dataset)
- Mulliken/Loewdin Orbital Charges
- J-Coupling Tensor Components (DSO/PSO/FC/SD/SD-FC)
- Chemical Shielding Tensors (full anisotropic)
- Chemical Shielding Summary

**MEDIUM PRIORITY (10 sections):**
- Geometric Perturbations
- Basis Set Details (contractions, exponents)
- SHARK Integral Package
- COSX Grid Generation
- Initial Guess Orbitals
- DIIS/SOSCF Details
- CPCM Detailed Parameters

**LOW PRIORITY (6 sections):**
- Normal Modes (complete all 63 modes)
- SCF Hessian Matrix
- MO Coefficients (MASSIVE dataset)
- Pople Linear Equation Solver
- SCF Settings (detailed)
- Citations

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
**Current Coverage:** 34/57 sections parsed (60%)

### Status Summary

| Category | Parsed | Remaining | Total |
|----------|--------|-----------|-------|
| **HIGH Priority** | 9 | 7 | 16 |
| **MEDIUM Priority** | 17 | 10 | 27 |
| **LOW Priority** | 8 | 6 | 14 |
| **TOTAL** | **34** | **23** | **57** |

### High Priority Unparsed Sections (7)

#### Electronic Structure - Charge Analysis (1 section)

**29. MULLIKEN ORBITAL CHARGES** ✅ DONE
- **Location:** Line 74800
- **Data:** Per-MO charge distribution (371 MO entries)
- **Status:** Implemented with uncorrected charges

**30. LOEWDIN ORBITAL CHARGES** ✅ DONE
- **Location:** Line 90021
- **Data:** Loewdin per-MO charge distribution (371 MO entries)
- **Status:** Implemented

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

### Medium Priority Unparsed Sections (10)

#### Geometry & Coordinates (1 section)

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

### Low Priority Unparsed Sections (6)

**50. NORMAL MODES (Complete)** - All 63 modes (3-4 hours)
**51. SCF HESSIAN** ⚠️ LARGE - 69×69 matrix (4-5 hours)
**52. MO COEFFICIENTS** ⚠️ MASSIVE - 15M values (6-8 hours)
**53. POPLE LINEAR SOLVER** - CPKS/CPHF details (2-3 hours)
**54. SCF SETTINGS (Detailed)** - All parameters (3-4 hours)
**55. CITATIONS** - Reference existing bibtex parser (1 hour)

### Implementation Roadmap

**Phase 1: High-Value NMR (8-10 hours)**
1. J-Coupling Tensor Components (4-6 hrs)
2. Chemical Shielding Tensors (3-4 hrs)

**Phase 2: Charge & Bond Analysis (2-3 hours)** ✅ Partially Complete
3. ✅ Loewdin Bond Orders (completed)
4. ✅ Mulliken Overlap Charges (completed)
5. Mulliken/Loewdin Orbital Charges (2-3 hrs)

**Phase 3: Geometry Details** ✅ COMPLETED
6. ✅ Internal Coordinates (completed)
7. ✅ Cartesian Coordinates (a.u.) (completed)

**Phase 4: Computational Details (6-8 hours)**
8. Enhanced Basis Set Info (4-5 hrs)
9. Enhanced CPCM Solvation (2-3 hrs)
10. SHARK Integral Package (1-2 hrs)

**Phase 5: Advanced/Optional (15-20 hours)**
11. MO Populations Per MO (6-8 hrs)
12. Complete Normal Modes (3-4 hrs)
13. SCF Hessian (4-5 hrs)

**Total Estimated Effort:** 30-42 hours for all remaining high/medium priority sections (reduced from 37-49 after completing 4 sections)

---

## Quick Continuation Guide

**When context runs out, use this guide to continue development:**

### Current State (Session Checkpoint)
- **Branch:** `claude/parse-orca-output-014tZra7WsNVQwk3ZtzDJnCk`
- **Coverage:** 34/57 sections (60%)
- **Test File:** `p1xs0p.out` (113,234 lines, 23 atoms)
- **Last Commit:** 005e8b0 - "Add 4 new ORCA output parsers (32/57 sections, 56% coverage)"

### Recent Additions (Latest Session)
1. ✅ Cartesian coordinates (a.u.) - line 347, with atomic numbers/masses
2. ✅ Internal coordinates - lines 375, 402 (Z-matrix format)
3. ✅ Mulliken overlap charges - line 75366 (105 atom pairs)
4. ✅ Total run time - line 113339 (days/hours/min/sec/msec)
5. ✅ Mulliken orbital charges - line 74800 (371 MO entries, per-orbital distribution)
6. ✅ Loewdin orbital charges - line 90021 (371 MO entries, per-orbital distribution)

### Development Pattern (Follow This)
```bash
# 1. Start from designated branch
git checkout claude/parse-orca-output-014tZra7WsNVQwk3ZtzDJnCk

# 2. Read the parser file first
Read: /home/user/Orca_Files/parsers/out_parser.py

# 3. Find section in test file
grep -n "SECTION NAME" p1xs0p.out
sed -n 'LINE_START,LINE_END p' p1xs0p.out

# 4. Implement parser
# - Add @dataclass for complex data
# - Add parse_SECTION() function
# - Call in parse_out_content()
# - Update OrcaOutput dataclass
# - Update to_dict() for JSON
# - Update __main__ display

# 5. Test parser
python parsers/out_parser.py p1xs0p.out | grep "SECTION"

# 6. Commit and push
git add -A
git commit -m "Add [section name] parsing"
git push -u origin claude/parse-orca-output-014tZra7WsNVQwk3ZtzDJnCk
```

### Next Priority Sections (Recommended Order)
**Quick Wins (1-2 hours each):**
1. Chemical Shielding Summary (line 3079) - 18 nuclei with isotropic values
2. Mulliken Orbital Charges (line 74800) - per-orbital charge distribution
3. Loewdin Orbital Charges (line 75015) - alternative charges

**High Value (4-6 hours each):**
4. J-Coupling Tensor Components (lines 107219-109419) - DSO/PSO/FC/SD tensors
5. Chemical Shielding Tensors (lines 88583-103469) - full anisotropic NMR

**Important Implementation Notes:**
- All parsers return `Optional[T]` or `list[T]`
- Use `re.search()` with `re.DOTALL` for section extraction
- Use `re.findall()` with `re.MULTILINE` for line-by-line parsing
- Test with actual data before committing
- Update README coverage stats after each commit

### File Locations
- **Main Parser:** `/home/user/Orca_Files/parsers/out_parser.py` (1600+ lines)
- **Test File:** `/home/user/Orca_Files/p1xs0p.out`
- **README:** `/home/user/Orca_Files/README.md`
- **Preview:** `/home/user/Orca_Files/previews/unified_preview.py`

### Key Code Patterns
```python
# Pattern 1: Simple extraction
@dataclass
class SectionData:
    field1: float = 0.0
    field2: str = ""

def parse_section(content: str) -> Optional[SectionData]:
    section = re.search(r'SECTION NAME.*?-+\s*(.*?)(?:\n\n|$)', content, re.DOTALL)
    if section:
        # Parse data
        return SectionData(...)
    return None

# Pattern 2: List extraction
def parse_section_list(content: str) -> list[tuple]:
    data = []
    section = re.search(r'SECTION.*?-+\s*(.*?)(?:\n\n|$)', content, re.DOTALL)
    if section:
        matches = re.findall(r'PATTERN', section.group(1), re.MULTILINE)
        for match in matches:
            data.append((match[0], float(match[1])))
    return data

# Pattern 3: Line-by-line parsing (for complex formats)
def parse_complex_section(content: str) -> list[DataClass]:
    items = []
    section_match = re.search(r'SECTION.*?\n(.*?)(?:\n\n)', content, re.DOTALL)
    if section_match:
        lines = section_match.group(1).split('\n')
        for line in lines:
            # Parse each line
            items.append(DataClass(...))
    return items
```

---

## Test Plan

### Test Strategy Overview
**Goal:** Ensure all 32 parsed sections extract correct data with proper error handling.

### 1. Unit Tests (Per Parser Function)

**Location:** Create `/home/user/Orca_Files/tests/test_out_parser.py`

**Test Categories:**

#### A. Data Extraction Tests
```python
def test_parse_coordinates():
    """Test Cartesian coordinate parsing."""
    # Test normal case
    result = parse_out_file('p1xs0p.out')
    assert len(result.coordinates) == 23
    assert result.coordinates[0][0] == 'C'  # Element
    assert isinstance(result.coordinates[0][1], float)  # X

def test_parse_coordinates_au():
    """Test atomic unit coordinates with mass."""
    result = parse_out_file('p1xs0p.out')
    assert len(result.coordinates_au) == 23
    assert result.coordinates_au[0][4] == 6.0  # Atomic number
    assert result.coordinates_au[0][5] == 12.011  # Mass

def test_parse_internal_coordinates():
    """Test Z-matrix internal coordinates."""
    result = parse_out_file('p1xs0p.out')
    assert len(result.internal_coords) == 23
    # Second atom should have bond to atom 1
    assert result.internal_coords[1].bond_to == 1
    assert result.internal_coords[1].bond_length > 0
```

#### B. Edge Case Tests
```python
def test_missing_section():
    """Test graceful handling of missing sections."""
    # Create test file without certain section
    result = parse_out_content("MINIMAL CONTENT")
    assert result.mulliken_overlap_charges == []
    assert result.internal_coords == []

def test_empty_output():
    """Test with minimal ORCA output."""
    result = parse_out_content("")
    assert result.final_energy is None
    assert result.coordinates == []
```

#### C. Data Integrity Tests
```python
def test_energy_consistency():
    """Verify final energy matches last SCF energy."""
    result = parse_out_file('p1xs0p.out')
    assert result.final_energy is not None
    if result.scf_energies:
        # Within numerical precision
        assert abs(result.final_energy - result.scf_energies[-1]) < 1e-6

def test_coordinate_consistency():
    """Verify coordinate count consistency."""
    result = parse_out_file('p1xs0p.out')
    # All coordinate formats should have same atom count
    assert len(result.coordinates) == len(result.coordinates_au)
    assert len(result.coordinates) == len(result.internal_coords)
    if result.mulliken_charges:
        assert len(result.coordinates) == len(result.mulliken_charges)
```

#### D. Performance Tests
```python
def test_parsing_speed():
    """Ensure parsing completes within time limit."""
    import time
    start = time.time()
    result = parse_out_file('p1xs0p.out')  # 113k lines
    elapsed = time.time() - start
    assert elapsed < 10.0  # Should complete in <10 seconds

def test_memory_efficiency():
    """Monitor memory usage during parsing."""
    import tracemalloc
    tracemalloc.start()
    result = parse_out_file('p1xs0p.out')
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 500 * 1024 * 1024  # <500 MB peak memory
```

### 2. Integration Tests

#### A. JSON Export Test
```python
def test_json_serialization():
    """Test complete JSON export."""
    result = parse_out_file('p1xs0p.out')
    data = result.to_dict()

    # Verify all sections present
    assert 'job_info' in data
    assert 'coordinates' in data
    assert 'coordinates_au' in data
    assert 'internal_coords' in data
    assert 'mulliken_overlap_charges' in data

    # Test JSON serialization
    import json
    json_str = json.dumps(data, indent=2)
    assert len(json_str) > 0

    # Test deserialization
    parsed = json.loads(json_str)
    assert parsed['final_energy'] == data['final_energy']
```

#### B. Cross-Format Consistency
```python
def test_xyz_vs_out_coordinates():
    """Compare coordinates from .xyz and .out files."""
    xyz_result = parse_xyz_file('p1xs0p.xyz')
    out_result = parse_out_file('p1xs0p.out')

    assert len(xyz_result.atoms) == len(out_result.coordinates)
    # Coordinates should match within tolerance
    for i, (xyz_atom, out_coord) in enumerate(zip(xyz_result.atoms, out_result.coordinates)):
        assert xyz_atom.element == out_coord[0]
        assert abs(xyz_atom.x - out_coord[1]) < 0.001
```

### 3. Regression Tests

**Track Changes Across Versions:**
```python
def test_backward_compatibility():
    """Ensure new parsers don't break old functionality."""
    result = parse_out_file('p1xs0p.out')

    # Original sections should still work
    assert result.final_energy is not None
    assert len(result.frequencies) == 63
    assert len(result.mayer_bond_orders) == 66

    # New sections should work
    assert len(result.coordinates_au) == 23
    assert len(result.mulliken_overlap_charges) == 105
```

### 4. Manual Verification Checklist

**For Each New Parser:**
- [ ] View raw section format in test file
- [ ] Verify regex pattern matches all lines
- [ ] Check first/last entries are captured
- [ ] Test with minimal section (if available)
- [ ] Verify data types are correct
- [ ] Check for off-by-one errors in atom indexing
- [ ] Confirm units are documented
- [ ] Validate against known reference values

### 5. Test Execution

**Run All Tests:**
```bash
# Run unit tests
python -m pytest tests/test_out_parser.py -v

# Run with coverage
python -m pytest tests/test_out_parser.py --cov=parsers --cov-report=html

# Run specific test
python -m pytest tests/test_out_parser.py::test_parse_coordinates -v

# Run performance tests
python -m pytest tests/test_out_parser.py -k "performance" -v
```

**Quick Manual Test:**
```bash
# Test parser output
python parsers/out_parser.py p1xs0p.out

# Test JSON export
python parsers/out_parser.py p1xs0p.out > output.json

# Verify section counts
python parsers/out_parser.py p1xs0p.out | grep -E "(Coordinates|Bond Orders|Overlap)"
```

### 6. Test Data Requirements

**Test Files Needed:**
1. **p1xs0p.out** (primary, 113k lines) - B3LYP/NMR/freq calculation
2. **minimal.out** (create) - Minimal SCF calculation
3. **optimization.out** (optional) - Geometry optimization
4. **tddft.out** (optional) - TD-DFT calculation

**Create Minimal Test File:**
```bash
# Extract essential sections for quick testing
head -n 500 p1xs0p.out > tests/fixtures/minimal.out
```

---

## UI Master Plan: Unified ORCA File Reader

### Vision
**Single-window application with tabbed interface to view all ORCA output files in one place.**

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│  ORCA Quantum Chemistry File Viewer                         │
│  File: p1xs0p.out                                   [×]     │
├─────────────────────────────────────────────────────────────┤
│  📁 File  🔧 Tools  📊 Export  ❓ Help                      │
├──────────┬──────────────────────────────────────────────────┤
│          │  Tab Bar                                         │
│  File    │  [Summary] [Geometry] [Energy] [Spectrum]       │
│  Tree    │  [Orbitals] [NMR] [Population] [Advanced]       │
│          │                                                  │
│  📄 Info │──────────────────────────────────────────────────│
│  ├─ Geo  │                                                  │
│  ├─ SCF  │         Main Content Area                       │
│  ├─ Opt  │      (Dynamic based on selected tab)            │
│  ├─ Freq │                                                  │
│  ├─ NMR  │                                                  │
│  └─ Orb  │                                                  │
│          │                                                  │
└──────────┴──────────────────────────────────────────────────┘
```

### Technology Stack

**Option A: Web-Based (Recommended)**
- **Frontend:** React + TypeScript
- **UI Framework:** Material-UI or Chakra UI
- **3D Viewer:** 3Dmol.js
- **Plots:** Plotly.js or Recharts
- **Backend:** FastAPI (Python)
- **Deployment:** Electron wrapper for desktop, or pure web

**Option B: Native Python**
- **Framework:** PyQt6 or tkinter
- **3D Viewer:** PyMOL or py3Dmol
- **Plots:** Matplotlib embedded
- **Advantage:** Single codebase, no web dependency

### Tab Structure

#### Tab 1: Summary 📊
**Quick overview of calculation**
```
┌────────────────────────────────────────────────┐
│  Calculation Summary                           │
│  ────────────────────────────────────────────  │
│  Method: B3LYP         Basis: pcSseg-3        │
│  Charge: 0             Multiplicity: 1        │
│  Atoms: 23             Electrons: 102         │
│                                                │
│  Final Energy: -662.998375 Eh                 │
│  Gibbs Free Energy: -662.858808 Eh            │
│                                                │
│  Status: ✅ SCF Converged (16 iterations)     │
│          ✅ Geometry Optimized                │
│          ✅ Frequencies Computed (63 modes)   │
│          ✅ NMR Calculated (18 shifts)        │
│                                                │
│  Runtime: 1h 3min 0s                          │
│  Files: 15 output files                       │
└────────────────────────────────────────────────┘
```

#### Tab 2: Geometry 🧬
**3D structure + coordinates**
```
┌────────────────┬───────────────────────────────┐
│  3D Viewer     │  Coordinate Table             │
│                │  ───────────────────────────  │
│   [Rotate]     │  Atom  Element  X    Y    Z  │
│   [Zoom]       │  ────────────────────────────│
│   [Reset]      │  0     C        0.98  4.01   │
│                │  1     H        0.94  6.07   │
│  [Structure]   │  2     C        2.37  2.72   │
│                │  ...                          │
│  Cartoon  ○    │                               │
│  Ball+Stick ●  │  Format: [Angstrom ▼]        │
│  Surface  ○    │  [Export XYZ] [Copy]         │
│                │                               │
│  Show: ☑ Labels│  Internal Coordinates         │
│        ☑ Bonds │  ───────────────────────────  │
│        ☐ Axes  │  Bond: H(1)-C(0) = 1.083 Å  │
└────────────────┴───────────────────────────────┘
```

#### Tab 3: Energy ⚡
**SCF, optimization, components**
```
┌──────────────────────────────────────────────────┐
│  SCF Convergence Plot                            │
│  ──────────────────────────────────────────────  │
│  [Line chart: Energy vs Iteration]               │
│                                                   │
│  Energy Components                                │
│  ──────────────────────────────────────────────  │
│  Nuclear Repulsion:      839.662 Eh              │
│  Electronic Energy:   -1502.722 Eh              │
│  Kinetic Energy:        660.289 Eh              │
│  XC Energy:             -72.671 Eh              │
│  Virial Ratio:            2.0042                 │
│                                                   │
│  Solvation (CPCM)                                │
│  ──────────────────────────────────────────────  │
│  Dielectric Energy:      -0.020 Eh              │
│  Surface Charge:         -0.032                  │
└──────────────────────────────────────────────────┘
```

#### Tab 4: Spectrum 📈
**IR, Raman, UV-Vis**
```
┌───────────────────────────────────────────────┐
│  Spectrum Type: [IR ▼]                        │
│  ────────────────────────────────────────────  │
│  [Interactive plot with hover]                │
│                                                │
│  Peak List                                     │
│  Mode  Frequency  Intensity  Assignment       │
│  ────────────────────────────────────────────  │
│  7     44.9       1.2        C-C-C bend       │
│  8     95.3       0.5        O-C-O bend       │
│  ...                                           │
│                                                │
│  [Export CSV] [Export Image]                  │
└───────────────────────────────────────────────┘
```

#### Tab 5: Orbitals 🔬
**MO energies, densities**
```
┌────────────────────┬──────────────────────────┐
│  MO Energy Diagram │  Orbital Details         │
│  ────────────────  │  ──────────────────────  │
│  [Energy level     │  Orbital: HOMO           │
│   diagram with     │  Index: 393              │
│   HOMO-LUMO gap]   │  Energy: -0.000 eV      │
│                    │  Occupation: 2.000       │
│  LUMO: -0.814 eV  │                          │
│  ═══════════════   │  Contributions:          │
│  Gap: 0.814 eV    │  C(0): 34.2%            │
│  ───────────────   │  C(2): 28.1%            │
│  HOMO:  0.000 eV  │  N(10): 15.6%           │
│                    │                          │
│  [Show MO ▼]      │  [View Density]          │
└────────────────────┴──────────────────────────┘
```

#### Tab 6: NMR 🧲
**Chemical shifts, J-couplings**
```
┌──────────────────────────────────────────────┐
│  Chemical Shifts                             │
│  ──────────────────────────────────────────  │
│  Atom  Element  Isotropic  Anisotropy       │
│  ──────────────────────────────────────────  │
│  0     C        50.2       45.8             │
│  1     H        3.1        12.4             │
│  ...                                         │
│                                              │
│  J-Coupling Constants                        │
│  ──────────────────────────────────────────  │
│  [Matrix visualization: Atom pairs]          │
│                                              │
│  [Export to SpinWorks] [Export CSV]         │
└──────────────────────────────────────────────┘
```

#### Tab 7: Population 👥
**Mulliken, Loewdin, charges**
```
┌────────────────────┬─────────────────────────┐
│  Charge Analysis   │  Visualization          │
│  ────────────────  │  ─────────────────────  │
│  Type: [Mulliken▼]│  [3D structure colored  │
│                    │   by partial charge]    │
│  Atom  Charge  Pop │                         │
│  ────────────────  │  Red: δ+ (positive)    │
│  0 C   -0.145 6.14 │  Blue: δ- (negative)   │
│  1 H    0.215 0.78 │                         │
│  ...               │  [Export PNG]           │
│                    │                         │
│  Bond Orders       │  Overlap Charges        │
│  ────────────────  │  ─────────────────────  │
│  C(0)-H(1): 0.945 │  C(0)-H(1): 1.174      │
│  C(0)-C(2): 1.324 │  C(6)-C(8): -24.514    │
└────────────────────┴─────────────────────────┘
```

#### Tab 8: Advanced 🔧
**Raw data, logs, citations**
```
┌──────────────────────────────────────────────┐
│  Advanced Data                               │
│  ──────────────────────────────────────────  │
│  📋 Timing Data                              │
│     Total: 443.7 sec | Fock: 395.9 sec      │
│                                              │
│  ⚙️  DFT Grid Info                           │
│     291,858 points | Lebedev-590            │
│                                              │
│  📚 Basis Set                                │
│     1305 functions | 0 primitives           │
│                                              │
│  📄 Full Output Log                          │
│     [Scrollable text viewer]                │
│                                              │
│  📖 Citations                                │
│     [BibTeX entries]                         │
└──────────────────────────────────────────────┘
```

### Implementation Phases

#### Phase 1: Backend API (2-3 days)
```python
# /home/user/Orca_Files/api/main.py
from fastapi import FastAPI, UploadFile
from parsers.out_parser import parse_out_file

app = FastAPI()

@app.post("/parse/out")
async def parse_out(file: UploadFile):
    content = await file.read()
    result = parse_out_content(content.decode())
    return result.to_dict()

@app.get("/formats")
def get_supported_formats():
    return {"formats": [".out", ".xyz", ".hess", ...]}
```

#### Phase 2: Frontend Shell (2-3 days)
```typescript
// src/App.tsx
import { Tabs, Tab } from '@mui/material';
import SummaryTab from './tabs/SummaryTab';
import GeometryTab from './tabs/GeometryTab';
// ...

function App() {
  const [activeTab, setActiveTab] = useState(0);
  const [orcaData, setOrcaData] = useState(null);

  return (
    <div>
      <FileUploader onParse={setOrcaData} />
      <Tabs value={activeTab} onChange={setActiveTab}>
        <Tab label="Summary" />
        <Tab label="Geometry" />
        {/* ... */}
      </Tabs>
      {activeTab === 0 && <SummaryTab data={orcaData} />}
      {/* ... */}
    </div>
  );
}
```

#### Phase 3: Individual Tabs (1-2 days each)
- Implement each tab component
- Add visualizations
- Add export functionality

#### Phase 4: Polish & Deploy (2-3 days)
- Error handling
- Loading states
- Responsive design
- Electron packaging

**Total Estimated Time: 15-20 days**

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

