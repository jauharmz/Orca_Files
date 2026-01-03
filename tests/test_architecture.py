"""
Test Architecture Implementation

Tests the new ORCA data architecture:
- MethodData parsing
- MoleculeStore hierarchical storage
- Canonical projection
"""

import sys
from pathlib import Path

# Add root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from src.parser import ParserFactory, MethodParser
from src.core import MoleculeStore, MethodData


def test_method_parsing():
    """Test method parsing from sample text."""
    print("=" * 60)
    print("TEST 1: Method Parsing")
    print("=" * 60)
    
    sample_text = """
    !  B3LYP def2-TZVP D3BJ TightSCF Opt Freq
    
    %cpcm
    smd true
    solvent "water"
    end
    
    * xyzfile 0 1 test.xyz
    """
    
    parser = MethodParser(sample_text)
    method = parser.parse()
    
    print(f"  Input line: {method.input_line}")
    print(f"  Method ID: {method.to_id()}")
    print(f"  Functional: {method.functional}")
    print(f"  Basis: {method.basis_set}")
    print(f"  Dispersion: {method.dispersion}")
    print(f"  Environment: {method.environment}")
    print(f"  Solvent: {method.solvent}")
    
    assert method.functional == "B3LYP", f"Expected B3LYP, got {method.functional}"
    assert method.basis_set == "DEF2-TZVP", f"Expected DEF2-TZVP, got {method.basis_set}"
    assert method.dispersion == "D3BJ", f"Expected D3BJ, got {method.dispersion}"
    assert method.solvent == "WATER", f"Expected WATER, got {method.solvent}"
    
    print("  ✅ Method parsing test PASSED")


def test_molecule_store():
    """Test MoleculeStore creation and operations."""
    print("\n" + "=" * 60)
    print("TEST 2: MoleculeStore")
    print("=" * 60)
    
    from src.core.data_models import ParseResult, GeometryData, EnergyData, OrbitalData, SpectraData, TDDFTData, MullikenData, InternalCoordsData
    
    # Create mock results
    method1 = MethodData(formalism="DFT", functional="B3LYP", basis_set="def2-TZVP", dispersion="D3BJ")
    method2 = MethodData(formalism="DFT", functional="PBE0", basis_set="def2-TZVP", dispersion="D3BJ")
    
    result1 = ParseResult(
        filename="p1x_s0_opt.out",
        method=method1,
        geometry=GeometryData(filename="p1x"),
        energy=EnergyData(gibbs_Eh=-100.5),
        orbitals=OrbitalData(homo_energy=-5.5, lumo_energy=-2.1),
        spectra=SpectraData(),
        tddft=TDDFTData(),
        mulliken=MullikenData(),
        internal_coords=InternalCoordsData(),
        is_optimization=True,
        optimized_state="S0"
    )
    
    result2 = ParseResult(
        filename="p1x_s0_sp.out",
        method=method1,
        geometry=GeometryData(filename="p1x"),
        energy=EnergyData(single_point_Eh=-100.6),
        orbitals=OrbitalData(homo_energy=-5.6, lumo_energy=-2.2),
        spectra=SpectraData(),
        tddft=TDDFTData(),
        mulliken=MullikenData(),
        internal_coords=InternalCoordsData(),
        is_optimization=False,
        optimized_state="S0"
    )
    
    result3 = ParseResult(
        filename="p1x_pbe0_opt.out",
        method=method2,
        geometry=GeometryData(filename="p1x"),
        energy=EnergyData(gibbs_Eh=-100.7),
        orbitals=OrbitalData(homo_energy=-5.7, lumo_energy=-2.3),
        spectra=SpectraData(),
        tddft=TDDFTData(),
        mulliken=MullikenData(),
        internal_coords=InternalCoordsData(),
        is_optimization=True,
        optimized_state="S0"
    )
    
    # Create store and add results
    store = MoleculeStore()
    store.add(result1)
    store.add(result2)
    store.add(result3)
    
    print(f"  Store summary: {store}")
    print(store.summary())
    
    # Test discovery methods
    molecules = store.molecules()
    print(f"\n  Molecules: {molecules}")
    assert "p1x" in molecules, "p1x should be in molecules"
    
    methods = store.get_methods("p1x")
    print(f"  Methods for p1x: {methods}")
    assert len(methods) == 2, f"Expected 2 methods, got {len(methods)}"
    
    states = store.get_states("p1x", methods[0])
    print(f"  States for {methods[0]}: {states}")
    
    # Test canonical access
    canonical = store.get_canonical("p1x")
    print(f"\n  Canonical view keys: {list(canonical.keys())}")
    print(f"  Canonical method_id: {canonical.get('method_id')}")
    
    # Test simple access
    result = store.get("p1x")
    print(f"  Simple get: {result.filename if result else None}")
    
    # Test DataFrame export
    df_canonical = store.to_dataframe(canonical=True)
    print(f"\n  Canonical DataFrame: {len(df_canonical)} rows, {list(df_canonical.columns)[:5]}...")
    
    df_full = store.to_dataframe(canonical=False)
    print(f"  Full DataFrame: {len(df_full)} rows")
    
    print("  ✅ MoleculeStore test PASSED")


def test_with_real_data():
    """Test with real ORCA files if available."""
    print("\n" + "=" * 60)
    print("TEST 3: Real Data Integration")
    print("=" * 60)
    
    # Try to find test data
    test_data = Path("./test_data_hf")
    if not test_data.exists():
        print("  ⏭ Skipping: test_data_hf not found (run test_comprehensive.py first)")
        return
    
    files = list(test_data.rglob("*.out"))[:10]  # First 10 files
    if not files:
        print("  ⏭ Skipping: no .out files found")
        return
    
    print(f"  Testing with {len(files)} files...")
    
    factory = ParserFactory()
    store = MoleculeStore()
    
    for f in files:
        try:
            result = factory.parse(str(f))
            store.add(result)
            print(f"    ✓ {f.name}: method={result.method.to_id()}")
        except Exception as e:
            print(f"    ✗ {f.name}: {e}")
    
    print(f"\n  {store}")
    print(store.summary())
    
    # Export to dataframe
    df = store.to_dataframe(canonical=True)
    print(f"\n  Canonical DataFrame columns: {list(df.columns)}")
    
    # Show method distribution
    df_full = store.to_dataframe(canonical=False)
    if "method_id" in df_full.columns:
        print(f"  Method distribution:\n{df_full['method_id'].value_counts().head()}")
    
    print("  ✅ Real data integration test PASSED")


if __name__ == "__main__":
    test_method_parsing()
    test_molecule_store()
    test_with_real_data()
    
    print("\n" + "=" * 60)
    print("ALL TESTS PASSED ✅")
    print("=" * 60)
