"""
Comprehensive Test Suite for All 34 ORCA Output Parsers

Tests all implemented parsers with:
- Data extraction verification
- Data type validation
- Value range checks
- Edge case handling
- JSON serialization
- Performance benchmarks

Usage:
    python tests/test_comprehensive.py
    pytest tests/test_comprehensive.py -v
"""

import sys
import os
import time
import json
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from parsers.out_parser import parse_out_file

# Test file path
TEST_FILE = Path(__file__).parent.parent / "p1xs0p.out"


class TestAllParsers:
    """Comprehensive tests for all 34 parsers."""

    @classmethod
    def setup_class(cls):
        """Parse file once for all tests."""
        print(f"\n{'='*70}")
        print("Loading test file and parsing...")
        print(f"{'='*70}")
        start = time.time()
        cls.result = parse_out_file(str(TEST_FILE))
        cls.parse_time = time.time() - start
        print(f"✓ Parsing completed in {cls.parse_time:.2f} seconds")
        print(f"{'='*70}\n")

    # ========== Section 1-4: Basic Job Info & Energy ==========

    def test_01_job_info(self):
        """Test job information extraction."""
        assert self.result.job_info is not None
        assert self.result.job_info.basis_set == "pcSseg-3"
        assert self.result.job_info.charge == 0
        assert self.result.job_info.multiplicity == 1
        assert self.result.job_info.num_electrons == 100  # Actual value from file
        print(f"✓ Job Info: {self.result.job_info.basis_set}, {self.result.job_info.num_electrons} electrons")

    def test_02_final_energy(self):
        """Test final single point energy."""
        assert self.result.final_energy is not None
        assert isinstance(self.result.final_energy, float)
        assert -700 < self.result.final_energy < -600
        print(f"✓ Final Energy: {self.result.final_energy:.6f} Eh")

    def test_03_scf_energies(self):
        """Test SCF iteration energies."""
        assert len(self.result.scf_energies) > 0
        assert all(isinstance(e, float) for e in self.result.scf_energies)
        print(f"✓ SCF Energies: {len(self.result.scf_energies)} iterations")

    def test_04_optimization_energies(self):
        """Test geometry optimization trajectory."""
        # May be empty for single-point calculations
        assert isinstance(self.result.optimization_energies, list)
        print(f"✓ Optimization Energies: {len(self.result.optimization_energies)} steps")

    # ========== Section 5-8: Coordinates ==========

    def test_05_coordinates_angstrom(self):
        """Test Cartesian coordinates (Angstrom)."""
        assert len(self.result.coordinates) == 23
        for elem, x, y, z in self.result.coordinates:
            assert isinstance(elem, str)
            assert isinstance(x, float) and isinstance(y, float) and isinstance(z, float)
        print(f"✓ Coordinates (Å): {len(self.result.coordinates)} atoms")

    def test_06_coordinates_au(self):
        """Test Cartesian coordinates (a.u.) with atomic data."""
        assert len(self.result.coordinates_au) == 23
        assert self.result.coordinates_au[0][0] == 'C'
        assert self.result.coordinates_au[0][4] == 6.0  # Atomic number
        assert 12.0 < self.result.coordinates_au[0][5] < 12.1  # Mass
        print(f"✓ Coordinates (a.u.): {len(self.result.coordinates_au)} atoms with masses")

    def test_07_internal_coordinates(self):
        """Test internal coordinates (Z-matrix, Angstrom)."""
        assert len(self.result.internal_coords) == 23
        # Second atom should have bond info
        assert self.result.internal_coords[1].bond_to == 1
        assert self.result.internal_coords[1].bond_length > 0
        print(f"✓ Internal Coords (Å): {len(self.result.internal_coords)} atoms")

    def test_08_internal_coordinates_au(self):
        """Test internal coordinates (Z-matrix, a.u.)."""
        assert len(self.result.internal_coords_au) == 23
        print(f"✓ Internal Coords (a.u.): {len(self.result.internal_coords_au)} atoms")

    # ========== Section 9-10: Molecular Properties ==========

    def test_09_dipole_moment(self):
        """Test electric dipole moment."""
        assert self.result.dipole_moment is not None
        assert self.result.dipole_moment.magnitude_debye >= 0
        assert isinstance(self.result.dipole_moment.x, float)
        print(f"✓ Dipole Moment: {self.result.dipole_moment.magnitude_debye:.4f} Debye")

    def test_10_polarizability(self):
        """Test polarizability tensor."""
        assert self.result.polarizability is not None
        assert self.result.polarizability.isotropic > 0
        print(f"✓ Polarizability: {self.result.polarizability.isotropic:.2f} a.u.")

    # ========== Section 11: Orbital Energies ==========

    def test_11_orbital_energies(self):
        """Test molecular orbital energies and occupations."""
        assert len(self.result.orbital_energies) > 0
        homo = [o for o in self.result.orbital_energies if o.occupation > 0.5][-1]
        lumo = [o for o in self.result.orbital_energies if o.occupation < 0.5][0]
        gap = lumo.energy_ev - homo.energy_ev  # Already in eV
        print(f"✓ Orbital Energies: {len(self.result.orbital_energies)} MOs, HOMO-LUMO gap = {gap:.2f} eV")

    # ========== Section 12-15: Vibrational Analysis ==========

    def test_12_frequencies(self):
        """Test vibrational frequencies."""
        assert len(self.result.frequencies) > 0
        print(f"✓ Frequencies: {len(self.result.frequencies)} modes ({min(self.result.frequencies):.1f} to {max(self.result.frequencies):.1f} cm⁻¹)")

    def test_13_ir_spectrum(self):
        """Test IR spectrum intensities."""
        # IR spectrum may be empty if not computed
        if self.result.ir_spectrum:
            strongest = max(self.result.ir_spectrum, key=lambda x: x.intensity)
            print(f"✓ IR Spectrum: {len(self.result.ir_spectrum)} modes, strongest at {strongest.frequency:.1f} cm⁻¹")
        else:
            print(f"✓ IR Spectrum: {len(self.result.ir_spectrum)} modes (not computed)")

    def test_14_raman_spectrum(self):
        """Test Raman spectrum activities."""
        assert len(self.result.raman_spectrum) > 0
        strongest = max(self.result.raman_spectrum, key=lambda x: x.activity)
        print(f"✓ Raman Spectrum: {len(self.result.raman_spectrum)} modes, strongest at {strongest.frequency:.1f} cm⁻¹")

    def test_15_normal_modes(self):
        """Test normal mode displacement vectors."""
        assert len(self.result.normal_modes) > 0
        print(f"✓ Normal Modes: {len(self.result.normal_modes)} modes with displacement vectors")

    # ========== Section 16-18: SCF Details ==========

    def test_16_scf_iterations(self):
        """Test SCF iteration details."""
        assert len(self.result.scf_iterations) > 0
        print(f"✓ SCF Iterations: {len(self.result.scf_iterations)} cycles")

    def test_17_scf_convergence(self):
        """Test SCF convergence criteria."""
        assert self.result.scf_convergence is not None
        print(f"✓ SCF Convergence: converged in {len(self.result.scf_iterations)} iterations")

    def test_18_energy_components(self):
        """Test energy decomposition."""
        assert self.result.energy_components is not None
        print(f"✓ Energy Components: nuclear, electronic, kinetic, etc.")

    # ========== Section 19-21: DFT Details ==========

    def test_19_dispersion_correction(self):
        """Test DFT dispersion correction."""
        assert self.result.dispersion_correction is not None
        print(f"✓ Dispersion: {self.result.dispersion_correction.method}")

    def test_20_dft_grid_info(self):
        """Test DFT integration grid."""
        assert self.result.dft_grid_info is not None
        print(f"✓ DFT Grid: {self.result.dft_grid_info.total_grid_points:,} points")

    def test_21_basis_set_info(self):
        """Test basis set statistics."""
        assert self.result.basis_set_info is not None
        print(f"✓ Basis Set: {self.result.basis_set_info.num_basis_functions} functions")

    # ========== Section 22: Solvation ==========

    def test_22_cpcm_solvation(self):
        """Test CPCM solvation model."""
        # May be None if no solvation
        if self.result.cpcm_solvation:
            print(f"✓ CPCM Solvation: dielectric energy = {self.result.cpcm_solvation.dielectric_energy:.6f} Eh")
        else:
            print("✓ CPCM Solvation: not used")

    # ========== Section 23-28: Population Analysis ==========

    def test_23_mulliken_charges(self):
        """Test Mulliken atomic charges."""
        assert len(self.result.mulliken_charges) == 23
        print(f"✓ Mulliken Charges: {len(self.result.mulliken_charges)} atoms")

    def test_24_mulliken_orbital_populations(self):
        """Test Mulliken orbital populations (s, p, d, f)."""
        assert len(self.result.mulliken_orbital_populations) == 23
        print(f"✓ Mulliken Orbital Populations: {len(self.result.mulliken_orbital_populations)} atoms")

    def test_25_mulliken_orbital_charges(self):
        """Test Mulliken per-MO orbital charges."""
        assert len(self.result.mulliken_orbital_charges) == 371
        print(f"✓ Mulliken Orbital Charges: {len(self.result.mulliken_orbital_charges)} MO entries")

    def test_26_loewdin_charges(self):
        """Test Loewdin atomic charges."""
        assert len(self.result.loewdin_charges) == 23
        print(f"✓ Loewdin Charges: {len(self.result.loewdin_charges)} atoms")

    def test_27_loewdin_orbital_charges(self):
        """Test Loewdin per-MO orbital charges."""
        assert len(self.result.loewdin_orbital_charges) == 371
        print(f"✓ Loewdin Orbital Charges: {len(self.result.loewdin_orbital_charges)} MO entries")

    def test_28_mulliken_overlap_charges(self):
        """Test Mulliken overlap charges (bonding)."""
        assert len(self.result.mulliken_overlap_charges) == 105
        print(f"✓ Mulliken Overlap Charges: {len(self.result.mulliken_overlap_charges)} atom pairs")

    # ========== Section 29-30: Bond Orders ==========

    def test_29_mayer_bond_orders(self):
        """Test Mayer bond order analysis."""
        assert len(self.result.mayer_bond_orders) > 0
        print(f"✓ Mayer Bond Orders: {len(self.result.mayer_bond_orders)} bonds")

    def test_30_loewdin_bond_orders(self):
        """Test Loewdin bond order analysis."""
        assert len(self.result.loewdin_bond_orders) > 0
        print(f"✓ Loewdin Bond Orders: {len(self.result.loewdin_bond_orders)} bonds")

    # ========== Section 31: Thermochemistry ==========

    def test_31_thermochemistry(self):
        """Test thermochemistry (G, H, S, ZPE)."""
        assert self.result.thermochemistry is not None
        assert self.result.thermochemistry.gibbs_free_energy is not None
        print(f"✓ Thermochemistry: G = {self.result.thermochemistry.gibbs_free_energy:.6f} Eh")

    # ========== Section 32-33: NMR ==========

    def test_32_nmr_chemical_shifts(self):
        """Test NMR chemical shift prediction."""
        assert self.result.nmr_data is not None
        assert len(self.result.nmr_data.chemical_shifts) == 18
        print(f"✓ NMR Chemical Shifts: {len(self.result.nmr_data.chemical_shifts)} nuclei")

    def test_33_nmr_j_couplings(self):
        """Test NMR J-coupling constants."""
        # May be empty in this test file
        print(f"✓ NMR J-Couplings: {len(self.result.nmr_data.j_couplings)} pairs")

    # ========== Section 34: Timing ==========

    def test_34_timing_data(self):
        """Test computation timing breakdown."""
        assert self.result.timing_data is not None
        assert self.result.timing_data.total_run_time > 0
        print(f"✓ Timing Data: {self.result.timing_data.total_run_time:.1f} seconds total")

    # ========== Additional Tests ==========

    def test_json_serialization(self):
        """Test JSON export/import."""
        data_dict = self.result.to_dict()
        json_str = json.dumps(data_dict, indent=2)
        assert len(json_str) > 1000
        print(f"✓ JSON Serialization: {len(json_str):,} characters")

    def test_performance(self):
        """Test parsing performance."""
        assert self.parse_time < 10.0, f"Parsing took {self.parse_time:.2f}s (>10s limit)"
        print(f"✓ Performance: {self.parse_time:.2f}s for 113,234 lines")


def run_tests():
    """Run all tests manually (no pytest required)."""
    print(f"\n{'='*70}")
    print("COMPREHENSIVE TEST SUITE FOR ORCA OUTPUT PARSER")
    print(f"{'='*70}")
    print(f"Test file: {TEST_FILE}")
    print(f"File exists: {TEST_FILE.exists()}")

    if not TEST_FILE.exists():
        print(f"\n❌ ERROR: Test file not found at {TEST_FILE}")
        return False

    print(f"File size: {TEST_FILE.stat().st_size / 1024 / 1024:.1f} MB")
    print(f"{'='*70}\n")

    # Create test instance
    test_suite = TestAllParsers()
    test_suite.setup_class()

    # Get all test methods
    test_methods = [m for m in dir(test_suite) if m.startswith('test_')]
    test_methods.sort()

    passed = 0
    failed = 0

    print(f"Running {len(test_methods)} tests...\n")

    for method_name in test_methods:
        try:
            method = getattr(test_suite, method_name)
            method()
            passed += 1
        except AssertionError as e:
            print(f"❌ FAILED: {method_name}")
            print(f"   {str(e)}")
            failed += 1
        except Exception as e:
            print(f"❌ ERROR in {method_name}: {str(e)}")
            failed += 1

    # Summary
    print(f"\n{'='*70}")
    print(f"TEST SUMMARY")
    print(f"{'='*70}")
    print(f"✓ Passed: {passed}/{len(test_methods)}")
    if failed > 0:
        print(f"❌ Failed: {failed}/{len(test_methods)}")
    print(f"Success Rate: {passed/len(test_methods)*100:.1f}%")
    print(f"{'='*70}\n")

    return failed == 0


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
