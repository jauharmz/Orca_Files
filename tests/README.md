# Tests

## Quick Start

```bash
cd d:\Antigravity\Orca_Files
```

### 1. Parse Real Data + Export CSV
```bash
python tests/test_parser_csv.py
```
**Output:** `parsed_data.csv`

### 2. Unit Tests
```bash
pytest tests/test_parsers.py tests/test_analysis.py -v -s
```

### 3. Full Integration Test
```bash
python tests/test_real.py
```
**Output:** `test_output.html`, `test_export.json`, `test_export.csv`

---

## Test Files

| File | Purpose | Data |
|------|---------|------|
| `test_parsers.py` | Unit test parsers | Inline |
| `test_analysis.py` | Unit test analysis | Inline |
| `test_parser_csv.py` | Parse real data → CSV | HuggingFace |
| `test_real.py` | Full integration | HuggingFace |

---

## HuggingFace Data

Data is cached in `./test_data_hf/` after first download.

To force re-download:
```bash
rm -rf test_data_hf
python tests/test_parser_csv.py
```
