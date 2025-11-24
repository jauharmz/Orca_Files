"""Tests for Hessian parser."""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.hess_parser import (
    parse_hess_file,
    get_ir_spectrum_data,
    find_strongest_peaks
)


class TestHessParser:
    """Test Hessian file parsing."""

    def test_parse_real_hess_file(self):
        """Test parsing the actual p1xs0.hess file."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.hess'
        )

        if os.path.exists(filepath):
            hess = parse_hess_file(filepath)

            # Check basic data
            assert hess.num_modes == 69
            assert len(hess.frequencies) == 69
            assert len(hess.ir_spectrum) > 0

            # Check we have real modes
            real_freqs = [f for f in hess.frequencies if f > 0]
            assert len(real_freqs) == 63

    def test_get_ir_spectrum_data(self):
        """Test extraction of IR spectrum data."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.hess'
        )

        if os.path.exists(filepath):
            hess = parse_hess_file(filepath)
            freqs, intensities = get_ir_spectrum_data(hess)

            # Should have same length
            assert len(freqs) == len(intensities)
            # Should only have positive frequencies
            assert all(f > 0 for f in freqs)

    def test_find_strongest_peaks(self):
        """Test finding strongest IR peaks."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.hess'
        )

        if os.path.exists(filepath):
            hess = parse_hess_file(filepath)
            peaks = find_strongest_peaks(hess, 5)

            # Should return requested number
            assert len(peaks) == 5

            # Should be sorted by intensity
            for i in range(len(peaks) - 1):
                assert peaks[i].intensity >= peaks[i+1].intensity

    def test_to_dict(self):
        """Test conversion to dictionary."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0.hess'
        )

        if os.path.exists(filepath):
            hess = parse_hess_file(filepath)
            d = hess.to_dict()

            assert 'frequencies' in d
            assert 'ir_spectrum' in d
            assert 'num_real_modes' in d


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
