"""Tests for XYZ parser."""

import pytest
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.xyz_parser import (
    parse_xyz_file,
    parse_xyz_content,
    parse_trajectory_xyz,
    extract_energy,
    Atom,
    MoleculeGeometry
)


class TestXYZParser:
    """Test XYZ file parsing."""

    def test_parse_simple_xyz(self):
        """Test parsing a simple XYZ content."""
        content = """3
Water molecule
O   0.000   0.000   0.117
H   0.000   0.757  -0.469
H   0.000  -0.757  -0.469
"""
        geom = parse_xyz_content(content)

        assert geom.num_atoms == 3
        assert geom.comment == "Water molecule"
        assert len(geom.atoms) == 3
        assert geom.atoms[0].symbol == "O"
        assert geom.atoms[1].symbol == "H"

    def test_extract_energy_orca_format(self):
        """Test energy extraction from ORCA comment line."""
        comment = "Coordinates from ORCA-job p1xs0 E -662.998375154159"
        energy = extract_energy(comment)

        assert energy is not None
        assert abs(energy - (-662.998375154159)) < 1e-10

    def test_extract_energy_no_energy(self):
        """Test energy extraction when no energy present."""
        comment = "Water molecule"
        energy = extract_energy(comment)

        assert energy is None

    def test_parse_real_xyz_file(self):
        """Test parsing the actual p1xs0.xyz file."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.xyz'
        )

        if os.path.exists(filepath):
            geom = parse_xyz_file(filepath)

            assert geom.num_atoms == 23
            assert geom.energy is not None
            assert abs(geom.energy - (-662.998375154159)) < 1e-6
            assert geom.atoms[0].symbol == "C"

    def test_parse_trajectory_file(self):
        """Test parsing trajectory XYZ file with multiple frames."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0_trj.xyz'
        )

        if os.path.exists(filepath):
            frames = parse_trajectory_xyz(filepath)

            assert len(frames) > 1
            # All frames should have same number of atoms
            num_atoms = frames[0].num_atoms
            for frame in frames:
                assert frame.num_atoms == num_atoms

    def test_to_dict(self):
        """Test conversion to dictionary."""
        content = """2
Test E -100.5
C   0.0   0.0   0.0
H   1.0   0.0   0.0
"""
        geom = parse_xyz_content(content)
        d = geom.to_dict()

        assert d['num_atoms'] == 2
        assert d['energy'] == -100.5
        assert len(d['atoms']) == 2
        assert d['atoms'][0]['symbol'] == 'C'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
