"""
Parser Unit Tests

Run: pytest tests/test_parsers.py -v -s
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestGeometryParser:
    """Test geometry parsing."""
    
    def test_parse_cartesian_coords(self):
        """Test cartesian coordinate parsing."""
        from src.parser.geometry import GeometryParser
        
        sample_text = """
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
C      0.000000    0.000000    0.000000
H      1.089000    0.000000    0.000000
H     -0.363000    1.027350    0.000000

        """
        
        parser = GeometryParser(sample_text)
        result = parser.parse()
        
        print("\n=== GEOMETRY OUTPUT ===")
        print(f"Cart coords:\n{result.cart_coords}")
        print(f"SMILES: {result.smiles}")
        
        assert result.cart_coords is not None
        assert len(result.cart_coords) == 3
    
    def test_parse_geometry_info(self):
        """Test geometry info extraction."""
        from src.parser.geometry import GeometryParser
        
        sample_text = """
================================================================================
INPUT FILE
================================================================================
! OPT B3LYP DEF2-SVP
* xyzfile 0 1 molecule.xyz
****END OF INPUT****
        """
        
        parser = GeometryParser(sample_text)
        result = parser.parse()
        
        print("\n=== GEOMETRY INFO ===")
        print(f"Filename: {result.filename}")
        print(f"Charge: {result.charge}")
        print(f"Multiplicity: {result.multiplicity}")
        
        assert result.charge == 0
        assert result.multiplicity == 1


class TestEnergyParser:
    """Test energy parsing."""
    
    def test_parse_gibbs_energy(self):
        """Test Gibbs energy parsing."""
        from src.parser.energy import EnergyParser
        
        sample_text = """
Final Gibbs free energy         ...   -854.123456789 Eh
        """
        
        parser = EnergyParser(sample_text)
        result = parser.parse()
        
        print("\n=== ENERGY OUTPUT ===")
        print(f"Gibbs: {result.gibbs_Eh} Eh")
        
        assert result.gibbs_Eh is not None
        assert abs(result.gibbs_Eh - (-854.123456789)) < 0.001
    
    def test_parse_single_point(self):
        """Test single point energy parsing."""
        from src.parser.energy import EnergyParser
        
        sample_text = """
FINAL SINGLE POINT ENERGY      -854.987654321
        """
        
        parser = EnergyParser(sample_text)
        result = parser.parse()
        
        print(f"Single point: {result.single_point_Eh} Eh")
        
        assert result.single_point_Eh is not None


class TestOrbitalParser:
    """Test orbital parsing."""
    
    def test_parse_orbitals(self):
        """Test orbital parsing with HOMO/LUMO."""
        from src.parser.orbitals import OrbitalParser
        
        sample_text = """
ORBITAL ENERGIES
-----------------
  NO   OCC          E(Eh)            E(eV)
   0   2.0000    -19.123456       -520.5
   1   2.0000    -10.234567       -278.5
   2   2.0000     -0.456789        -12.4
   3   0.0000      0.123456          3.4
   4   0.0000      0.234567          6.4

        """
        
        parser = OrbitalParser(sample_text)
        result = parser.parse()
        
        print("\n=== ORBITAL OUTPUT ===")
        print(f"Orbitals:\n{result.orbitals}")
        print(f"HOMO: {result.homo_energy} eV")
        print(f"LUMO: {result.lumo_energy} eV")
        print(f"Gap: {result.homo_lumo_gap} eV")
        
        assert result.orbitals is not None
        assert result.homo_energy is not None
        assert result.lumo_energy is not None


class TestParserFactory:
    """Test parser factory."""
    
    def test_parse_text(self):
        """Test parsing from text."""
        from src.parser.factory import ParserFactory
        
        sample_text = """
================================================================================
INPUT FILE
================================================================================
! OPT B3LYP
* xyzfile 0 1 test.xyz
****END OF INPUT****

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
C      0.000000    0.000000    0.000000
H      1.089000    0.000000    0.000000

Final Gibbs free energy         ...   -100.123456 Eh
        """
        
        factory = ParserFactory()
        result = factory.parse_text(sample_text, "test.out")
        
        print("\n=== FACTORY OUTPUT ===")
        print(f"Filename: {result.filename}")
        print(f"Geometry atoms: {len(result.geometry.cart_coords) if result.geometry.cart_coords is not None else 0}")
        print(f"Energy: {result.energy.gibbs_Eh}")
        
        assert result.geometry.cart_coords is not None
        assert result.energy.gibbs_Eh is not None


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v", "-s"])
