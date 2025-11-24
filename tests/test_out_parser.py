"""Tests for ORCA output parser."""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.out_parser import (
    parse_out_file,
    parse_final_energy,
    parse_dipole_moment,
    parse_frequencies,
    parse_mulliken_charges,
    parse_thermochemistry
)


class TestOutParser:
    """Test ORCA output file parsing."""

    def test_parse_real_out_file(self):
        """Test parsing the actual p1xs0.out file."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.out'
        )

        if os.path.exists(filepath):
            result = parse_out_file(filepath)

            # Check final energy
            assert result.final_energy is not None
            assert abs(result.final_energy - (-662.998375)) < 0.001

            # Check job info
            assert result.job_info.charge == 0
            assert result.job_info.multiplicity == 1
            assert result.job_info.basis_set == '6-311++G(d,p)'

            # Check dipole moment
            assert result.dipole_moment is not None
            assert abs(result.dipole_moment.magnitude_debye - 4.686) < 0.01

            # Check frequencies
            assert len(result.frequencies) == 63
            assert min(result.frequencies) > 40
            assert max(result.frequencies) < 3600

            # Check thermochemistry
            assert result.thermochemistry is not None
            assert result.thermochemistry.gibbs_free_energy is not None

            # Check Mulliken charges
            assert len(result.mulliken_charges) > 0

    def test_parse_optimization_energies(self):
        """Test extraction of optimization cycle energies."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.out'
        )

        if os.path.exists(filepath):
            result = parse_out_file(filepath)

            # Should have multiple optimization cycles
            assert len(result.optimization_energies) > 1
            # Energy should decrease during optimization
            assert result.optimization_energies[-1] < result.optimization_energies[0]

    def test_to_dict(self):
        """Test conversion to dictionary."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.out'
        )

        if os.path.exists(filepath):
            result = parse_out_file(filepath)
            d = result.to_dict()

            assert 'final_energy' in d
            assert 'job_info' in d
            assert 'dipole_moment' in d
            assert 'frequencies' in d


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
