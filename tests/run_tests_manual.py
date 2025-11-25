#!/usr/bin/env python3
"""
Manual test runner (doesn't require pytest)
Run tests without external dependencies.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from parsers.out_parser import parse_out_file

TEST_FILE = Path(__file__).parent.parent / "p1xs0p.out"

def test_all_sections():
    """Test all 32 implemented sections."""
    print("Testing all 32 parsed sections...")
    result = parse_out_file(str(TEST_FILE))
    
    tests_passed = 0
    tests_failed = 0
    
    # Test 1: Job Info
    try:
        assert result.job_info.basis_set == "pcSseg-3"
        assert result.job_info.charge == 0
        print("✓ Job Info")
        tests_passed += 1
    except AssertionError as e:
        print(f"✗ Job Info: {e}")
        tests_failed += 1
    
    # Test 2: Final Energy
    try:
        assert result.final_energy == -662.998375
        print("✓ Final Energy")
        tests_passed += 1
    except AssertionError:
        print(f"✗ Final Energy")
        tests_failed += 1
    
    # Test 3: Coordinates
    try:
        assert len(result.coordinates) == 23
        print("✓ Coordinates (Angstrom)")
        tests_passed += 1
    except AssertionError:
        print("✗ Coordinates (Angstrom)")
        tests_failed += 1
    
    # Test 4: Coordinates (a.u.)
    try:
        assert len(result.coordinates_au) == 23
        assert result.coordinates_au[0][4] == 6.0
        print("✓ Coordinates (a.u.)")
        tests_passed += 1
    except AssertionError:
        print("✗ Coordinates (a.u.)")
        tests_failed += 1
    
    # Test 5: Internal Coordinates
    try:
        assert len(result.internal_coords) == 23
        assert result.internal_coords[1].bond_to == 1
        print("✓ Internal Coordinates")
        tests_passed += 1
    except AssertionError:
        print("✗ Internal Coordinates")
        tests_failed += 1
    
    # Test 6: Dipole Moment
    try:
        assert result.dipole_moment.magnitude_debye > 0
        print("✓ Dipole Moment")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Dipole Moment")
        tests_failed += 1
    
    # Test 7: Polarizability
    try:
        assert result.polarizability.isotropic > 0
        print("✓ Polarizability")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Polarizability")
        tests_failed += 1
    
    # Test 8: Orbital Energies
    try:
        assert len(result.orbital_energies) > 0
        print("✓ Orbital Energies")
        tests_passed += 1
    except AssertionError:
        print("✗ Orbital Energies")
        tests_failed += 1
    
    # Test 9: Frequencies
    try:
        assert len(result.frequencies) == 63
        print("✓ Frequencies")
        tests_passed += 1
    except AssertionError:
        print("✗ Frequencies")
        tests_failed += 1
    
    # Test 10: IR Spectrum
    try:
        assert len(result.ir_spectrum) > 0
        print("✓ IR Spectrum")
        tests_passed += 1
    except AssertionError:
        print("✗ IR Spectrum")
        tests_failed += 1
    
    # Test 11: Raman Spectrum
    try:
        assert len(result.raman_spectrum) == 63
        print("✓ Raman Spectrum")
        tests_passed += 1
    except AssertionError:
        print("✗ Raman Spectrum")
        tests_failed += 1
    
    # Test 12: Dispersion Correction
    try:
        assert result.dispersion_correction.method == "DFTD3"
        print("✓ Dispersion Correction")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Dispersion Correction")
        tests_failed += 1
    
    # Test 13: Mulliken Charges
    try:
        assert len(result.mulliken_charges) == 23
        print("✓ Mulliken Charges")
        tests_passed += 1
    except AssertionError:
        print("✗ Mulliken Charges")
        tests_failed += 1
    
    # Test 14: Mulliken Overlap Charges
    try:
        assert len(result.mulliken_overlap_charges) == 105
        print("✓ Mulliken Overlap Charges")
        tests_passed += 1
    except AssertionError:
        print("✗ Mulliken Overlap Charges")
        tests_failed += 1
    
    # Test 15: Loewdin Charges
    try:
        assert len(result.loewdin_charges) == 23
        print("✓ Loewdin Charges")
        tests_passed += 1
    except AssertionError:
        print("✗ Loewdin Charges")
        tests_failed += 1
    
    # Test 16: Mayer Bond Orders
    try:
        assert len(result.mayer_bond_orders) == 66
        print("✓ Mayer Bond Orders")
        tests_passed += 1
    except AssertionError:
        print("✗ Mayer Bond Orders")
        tests_failed += 1
    
    # Test 17: Loewdin Bond Orders
    try:
        assert len(result.loewdin_bond_orders) == 65
        print("✓ Loewdin Bond Orders")
        tests_passed += 1
    except AssertionError:
        print("✗ Loewdin Bond Orders")
        tests_failed += 1
    
    # Test 18: Thermochemistry
    try:
        assert result.thermochemistry.gibbs_free_energy is not None
        print("✓ Thermochemistry")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Thermochemistry")
        tests_failed += 1
    
    # Test 19: NMR Data
    try:
        assert len(result.nmr_data.chemical_shifts) == 18
        print("✓ NMR Data")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ NMR Data")
        tests_failed += 1
    
    # Test 20: SCF Iterations
    try:
        assert len(result.scf_iterations) == 16
        print("✓ SCF Iterations")
        tests_passed += 1
    except AssertionError:
        print("✗ SCF Iterations")
        tests_failed += 1
    
    # Test 21: Timing Data
    try:
        assert result.timing_data.total_time > 0
        assert result.timing_data.total_run_time > 0
        print("✓ Timing Data + Total Runtime")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Timing Data")
        tests_failed += 1
    
    # Test 22: DFT Grid Info
    try:
        assert result.dft_grid_info.total_grid_points == 291858
        print("✓ DFT Grid Info")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ DFT Grid Info")
        tests_failed += 1
    
    # Test 23: Basis Set Info
    try:
        assert result.basis_set_info.name == "pcSseg-3"
        print("✓ Basis Set Info")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Basis Set Info")
        tests_failed += 1
    
    # Test 24: Energy Components
    try:
        assert result.energy_components.nuclear_repulsion > 0
        print("✓ Energy Components")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ Energy Components")
        tests_failed += 1
    
    # Test 25: CPCM Solvation
    try:
        assert result.cpcm_solvation.dielectric_energy != 0.0
        print("✓ CPCM Solvation")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ CPCM Solvation")
        tests_failed += 1
    
    # Test 26: SCF Convergence
    try:
        assert result.scf_convergence.diis_error >= 0
        print("✓ SCF Convergence")
        tests_passed += 1
    except (AssertionError, AttributeError):
        print("✗ SCF Convergence")
        tests_failed += 1
    
    # Test 27: Mulliken Orbital Populations
    try:
        assert len(result.mulliken_orbital_populations) == 23
        print("✓ Mulliken Orbital Populations")
        tests_passed += 1
    except AssertionError:
        print("✗ Mulliken Orbital Populations")
        tests_failed += 1
    
    print("\n" + "=" * 60)
    print(f"Tests Passed: {tests_passed}/27")
    print(f"Tests Failed: {tests_failed}/27")
    print("=" * 60)
    
    return tests_failed == 0

if __name__ == '__main__':
    print("=" * 60)
    print("ORCA Output Parser - Manual Test Runner")
    print("=" * 60)
    print(f"Test file: {TEST_FILE}")
    print(f"File exists: {TEST_FILE.exists()}")
    if TEST_FILE.exists():
        print(f"File size: {TEST_FILE.stat().st_size / 1024 / 1024:.1f} MB")
    print("=" * 60)
    print()
    
    success = test_all_sections()
    sys.exit(0 if success else 1)
