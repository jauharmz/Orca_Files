# ORCA Parser Tests

This directory contains tests for the ORCA quantum chemistry parser.

## Test Files

| Test File | Description | Output |
|-----------|-------------|--------|
| `test_comprehensive.py` | Full modular parser test | 16 CSV files + JSON |
| `test_original_parser.py` | Original `orca_praser.py` test | CSV + JSON |
| `test_comparison.py` | Side-by-side comparison | Console report |
| `test_visualizations.py` | Generate HTML plots | HTML files |
| `test_parsers.py` | Unit tests for parser modules | pytest |
| `test_analysis.py` | Analysis module tests | pytest |
| `test_real.py` | Integration tests | pytest |

## Quick Run

```bash
# Full parser test (recommended)
python tests/test_comprehensive.py

# Compare with original parser
python tests/test_original_parser.py
python tests/test_comparison.py

# Generate visualizations
python tests/test_visualizations.py

# Unit tests
pytest tests/test_parsers.py -v
pytest tests/test_analysis.py -v
```

## Output Files

### From `test_comprehensive.py`:

**Main data:**
- `parsed_molecules.csv` - Scalar fields + nested data counts

**Geometry:**
- `parsed_cart_coords.csv` - All atom coordinates

**Energy:**
- `parsed_orbitals.csv` - All orbital energies (OCC, Eh, eV, lvl)

**Spectroscopy:**
- `parsed_ir_spectra.csv` - IR spectrum
- `parsed_vibrations.csv` - Vibrational frequencies
- `parsed_raman.csv` - Raman spectrum
- `parsed_mulliken.csv` - Mulliken charges
- `parsed_nmr_shielding.csv` - NMR shielding
- `parsed_nmr_coupling.csv` - NMR J-coupling

**TD-DFT:**
- `parsed_tddft_states.csv` - Excited state transitions
- `parsed_electric_dipole.csv` - Electric dipole absorption
- `parsed_velocity_dipole.csv` - Velocity dipole absorption

**Internal Coordinates:**
- `parsed_internal_bonds.csv` - Bond data
- `parsed_internal_angles.csv` - Angle data
- `parsed_internal_dihedrals.csv` - Dihedral data

**Complete:**
- `parsed_data.json` - All data in JSON format

### From `test_visualizations.py`:

- `viz_energy.html` - Energy diagram
- `viz_pathway.html` - Pathway visualization
- `viz_orbitals.html` - Orbital energy levels
- `viz_ir.html` - IR spectrum
- `viz_uvvis.html` - UV-Vis spectrum

## Data Source

Tests automatically download sample data from HuggingFace:

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="JauharMz/Orca",
    repo_type="dataset",
    local_dir="./test_data_hf"
)
```

## 28 Data Components

The parser extracts 28 data fields:

| Category | Fields |
|----------|--------|
| Identity | molecule_id, smiles, charge, multiplicity |
| Energy | gibbs_Eh, single_point_Eh |
| Orbitals | homo_energy, lumo_energy, homo_lumo_gap, orbitals |
| Geometry | cart_coords, bonds, angles, dihedrals |
| Spectroscopy | ir, vibrations, raman, mulliken, nmr_shielding, nmr_coupling |
| TD-DFT | tddft_states, electric_dipole_abs/soc, velocity_dipole_abs/soc |
| Metadata | is_optimization, optimized_state, calc_class, esd_type |
