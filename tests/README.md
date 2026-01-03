# Tests

## Quick Start

```bash
cd d:\Antigravity\Orca_Files
```

### 1. Basic Parser Test
```bash
python tests/test_parser_csv.py
```
**Output:** `parsed_data.csv`

### 2. Comprehensive Test (Enhanced Export + JSON + Comparison)
```bash
python tests/test_comprehensive.py
```
**Output:**
- `parsed_data_full.csv` - with nested data counts
- `parsed_data.json` - ALL data
- `comparison_report.txt` - modular vs original comparison

### 3. Visualization Test
```bash
python tests/test_visualizations.py
```
**Output:**
- `viz_energy.html` - energy bar chart
- `viz_pathway.html` - reaction pathway
- `viz_orbitals.html` - HOMO/LUMO levels
- `viz_ir.html` - IR spectrum
- `viz_uvvis.html` - UV-Vis spectrum

### 4. Unit Tests
```bash
pytest tests/test_parsers.py tests/test_analysis.py -v -s
```

---

## Test Files

| File | Purpose |
|------|---------|
| `test_parser_csv.py` | Basic parsing → CSV |
| `test_comprehensive.py` | Full export + comparison |
| `test_visualizations.py` | Generate HTML plots |
| `test_parsers.py` | Parser unit tests |
| `test_analysis.py` | Analysis unit tests |

---

## HuggingFace Data

Cached in `./test_data_hf/`. To re-download:
```bash
rm -rf test_data_hf
```
