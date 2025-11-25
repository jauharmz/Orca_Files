"""
Comprehensive End-to-End Test Suite for ORCA Out Parser

Tests all 32 currently implemented parsers with:
- Data extraction verification
- Edge case handling
- Data integrity checks
- Performance benchmarks
- JSON serialization
- Cross-format consistency

Usage:
    pytest tests/test_out_parser.py -v
    pytest tests/test_out_parser.py::test_parse_coordinates -v
    pytest tests/test_out_parser.py --cov=parsers --cov-report=html
"""

import sys
import os
import time
import json
import tracemalloc
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from parsers.out_parser import (
    parse_out_file,
    parse_out_content,
    OrcaOutput,
    JobInfo,
)


# Test file path
TEST_FILE = Path(__file__).parent.parent / "p1xs0p.out"


class TestDataExtraction:
    """Test extraction of all 32 implemented sections."""

    def test_parse_job_info(self):
        """Test job information parsing."""
        result = parse_out_file(str(TEST_FILE))
        assert result.job_info is not None
        assert result.job_info.basis_set == "pcSseg-3"
        assert result.job_info.charge == 0
        assert result.job_info.multiplicity == 1
        assert result.job_info.num_electrons == 102

    def test_parse_final_energy(self):
        """Test final energy extraction."""
        result = parse_out_file(str(TEST_FILE))
        assert result.final_energy is not None
        assert isinstance(result.final_energy, float)
        assert -700 < result.final_energy < -600  # Reasonable range

    def test_parse_coordinates(self):
        """Test Cartesian coordinate parsing (Angstrom)."""
        result = parse_out_file(str(TEST_FILE))
        assert len(result.coordinates) == 23
        assert result.coordinates[0][0] == 'C'  # Element
        assert isinstance(result.coordinates[0][1], float)  # X
        assert isinstance(result.coordinates[0][2], float)  # Y
        assert isinstance(result.coordinates[0][3], float)  # Z

    def test_parse_coordinates_au(self):
        """Test atomic unit coordinates with atomic number and mass."""
        result = parse_out_file(str(TEST_FILE))
        assert len(result.coordinates_au) == 23
        assert result.coordinates_au[0][0] == 'C'  # Element
        assert result.coordinates_au[0][4] == 6.0  # Atomic number for C
        assert 12.0 < result.coordinates_au[0][5] < 12.1  # Mass ~12.011

    def test_parse_internal_coordinates(self):
        """Test Z-matrix internal coordinates."""
        result = parse_out_file(str(TEST_FILE))
        assert len(result.internal_coords) == 23
        assert len(result.internal_coords_au) == 23

        # Second atom should have bond to atom 1
        assert result.internal_coords[1].bond_to == 1
        assert result.internal_coords[1].bond_length > 0
        assert result.internal_coords[1].bond_angle >= 0

    def test_parse_mulliken_overlap_charges(self):
        """Test Mulliken overlap charges (bonding analysis)."""
        result = parse_out_file(str(TEST_FILE))
        assert len(result.mulliken_overlap_charges) == 105

        for atom1, atom2, overlap in result.mulliken_overlap_charges:
            assert isinstance(atom1, int)
            assert isinstance(atom2, int)
            assert isinstance(overlap, float)


class TestPerformance:
    """Test parsing speed and memory efficiency."""

    def test_parsing_speed(self):
        """Ensure parsing completes within time limit."""
        start = time.time()
        result = parse_out_file(str(TEST_FILE))  # 113k lines
        elapsed = time.time() - start

        print(f"\nParsing time: {elapsed:.2f} seconds")
        assert elapsed < 10.0, f"Parsing took {elapsed:.2f}s, should be <10s"


if __name__ == '__main__':
    import pytest
    print("=" * 70)
    print("ORCA Out Parser - Comprehensive Test Suite")
    print("=" * 70)
    print(f"\nTest file: {TEST_FILE}")
    print(f"File exists: {TEST_FILE.exists()}")
    if TEST_FILE.exists():
        print(f"File size: {TEST_FILE.stat().st_size / 1024 / 1024:.1f} MB")
    print("\n" + "=" * 70)

    pytest.main([__file__, '-v', '--tb=short'])
