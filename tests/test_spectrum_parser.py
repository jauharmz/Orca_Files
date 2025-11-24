"""Tests for spectrum parser."""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.spectrum_parser import (
    parse_spectrum_file,
    parse_spectrum_content,
    find_peaks,
    SpectrumData
)


class TestSpectrumParser:
    """Test spectrum file parsing."""

    def test_parse_real_spectrum_file(self):
        """Test parsing the actual p1xs0vg.spectrum file."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0vg.spectrum'
        )

        if os.path.exists(filepath):
            spectrum = parse_spectrum_file(filepath)

            # Check we have data
            assert len(spectrum.energy) > 0
            assert len(spectrum.total_spectrum) > 0
            assert len(spectrum.energy) == len(spectrum.total_spectrum)

            # Check energy conversion
            assert len(spectrum.energy_nm) == len(spectrum.energy)
            assert len(spectrum.energy_ev) == len(spectrum.energy)

    def test_find_peaks(self):
        """Test peak finding."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0vg.spectrum'
        )

        if os.path.exists(filepath):
            spectrum = parse_spectrum_file(filepath)
            peaks = find_peaks(spectrum, threshold=0.1)

            # Should find some peaks
            assert len(peaks) > 0

            # Peaks should be sorted by intensity
            for i in range(len(peaks) - 1):
                assert peaks[i]['intensity'] >= peaks[i+1]['intensity']

    def test_to_dict(self):
        """Test conversion to dictionary."""
        filepath = os.path.join(
            os.path.dirname(os.path.dirname(__file__)),
            'p1xs0vg.spectrum'
        )

        if os.path.exists(filepath):
            spectrum = parse_spectrum_file(filepath)
            d = spectrum.to_dict()

            assert 'energy_cm' in d
            assert 'energy_nm' in d
            assert 'total_spectrum' in d
            assert 'num_points' in d

    def test_energy_conversion(self):
        """Test energy unit conversion."""
        # Simple test data
        content = """Energy\tTotalSpectrum
10000\t1.0
20000\t2.0
"""
        spectrum = parse_spectrum_content(content)

        # Check nm conversion: nm = 1e7 / cm^-1
        assert abs(spectrum.energy_nm[0] - 1000) < 0.1
        assert abs(spectrum.energy_nm[1] - 500) < 0.1


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
