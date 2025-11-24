"""Tests for remaining parsers: engrad, property, inp, cpcm, bibtex."""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.engrad_parser import parse_engrad_file
from parsers.property_parser import parse_property_file
from parsers.inp_parser import parse_inp_file
from parsers.cpcm_parser import parse_cpcm_file
from parsers.bibtex_parser import parse_bibtex_file


class TestEngradParser:
    def test_parse_real_file(self):
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.engrad')
        if os.path.exists(filepath):
            data = parse_engrad_file(filepath)
            assert data.num_atoms == 23
            assert abs(data.energy - (-662.998375)) < 0.001
            assert len(data.gradients) == 69
            assert data.max_gradient > 0


class TestPropertyParser:
    def test_parse_real_file(self):
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.property.txt')
        if os.path.exists(filepath):
            data = parse_property_file(filepath)
            assert data.num_atoms == 23
            assert data.status == "NORMAL TERMINATION"
            assert len(data.mulliken_charges) == 23
            assert len(data.loewdin_charges) == 23


class TestInpParser:
    def test_parse_real_file(self):
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.inp')
        if os.path.exists(filepath):
            data = parse_inp_file(filepath)
            assert data.method == "B3LYP"
            assert data.basis_set == "6-311++G(d,p)"
            assert data.charge == 0
            assert data.multiplicity == 1
            assert len(data.coordinates) == 23


class TestCPCMParser:
    def test_parse_real_file(self):
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.cpcm')
        if os.path.exists(filepath):
            data = parse_cpcm_file(filepath)
            assert data.num_atoms == 23
            assert data.num_surface_points > 0
            assert data.volume > 0
            assert data.area > 0


class TestBibtexParser:
    def test_parse_real_file(self):
        filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'p1xs0.bibtex')
        if os.path.exists(filepath):
            data = parse_bibtex_file(filepath)
            assert len(data.citations) > 0
            assert data.citations[0].author != ""
            assert data.citations[0].year != ""


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
