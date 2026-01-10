"""
Parser Comparison Test - Compare Original vs Modular Parser

Runs both parsers and shows side-by-side comparison of coverage.

Run: python tests/test_comparison.py
"""

import sys
from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Add root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))


def count_df(x):
    """Count items in DataFrame or list."""
    if isinstance(x, pd.DataFrame):
        return len(x) if not x.empty else 0
    if isinstance(x, list):
        return len(x)
    return 0


def main():
    print("=" * 80)
    print("PARSER COMPARISON: Original (orca_praser.py) vs Modular (src/parser/)")
    print("=" * 80)
    
    # 1. Download from HuggingFace
    print("\n[STEP 1] Loading data from HuggingFace...")
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    if not data_dir.exists():
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
    
    all_files = list(data_dir.rglob("*.out"))
    print(f"  Found: {len(all_files)} files")
    
    # 2. Parse with ORIGINAL parser
    print("\n[STEP 2] Parsing with ORIGINAL parser...")
    from orca_praser import ORCAParser
    import re
    
    original_results = []
    for file_path in tqdm(all_files, desc="Original"):
        try:
            parser = ORCAParser(str(file_path))
            result = parser.parse(as_df=True)
            
            # Extract molecule_id
            base_name = file_path.stem
            mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt)$', '', base_name, flags=re.I)
            result["molecule_id"] = mol_id.strip("_-")
            result["source_file"] = file_path.name
            
            original_results.append(result)
        except Exception as e:
            print(f"  Error: {file_path.name}: {e}")
    
    # 3. Parse with MODULAR parser
    print("\n[STEP 3] Parsing with MODULAR parser...")
    from src.parser.factory import ParserFactory
    
    factory = ParserFactory()
    modular_results = []
    
    for file_path in tqdm(all_files, desc="Modular"):
        try:
            result = factory.parse(str(file_path))
            row = result.to_dict()
            row["source_file"] = file_path.name
            modular_results.append(row)
        except Exception as e:
            print(f"  Error: {file_path.name}: {e}")
    
    print(f"\n  Original: {len(original_results)} parsed")
    print(f"  Modular:  {len(modular_results)} parsed")
    
    # 4. Compare coverage
    print("\n" + "=" * 80)
    print("COMPARISON RESULTS")
    print("=" * 80)
    
    # Build comparison table
    comparisons = []
    
    # --- SCALAR FIELDS ---
    scalar_checks = [
        ("smiles", "smiles", "smiles"),
        ("charge", lambda r: r.get("geometry", {}).get("charge") if isinstance(r.get("geometry"), dict) else None, "charge"),
        ("multiplicity", lambda r: r.get("geometry", {}).get("multiplicity") if isinstance(r.get("geometry"), dict) else None, "multiplicity"),
        ("gibbs_Eh", "gibbs_energy_Eh", "gibbs_Eh"),
        ("single_point_Eh", "single_point_energy_Eh", "single_point_Eh"),
        ("is_optimization", "is_optimization", "is_optimization"),
        ("optimized_state", "optimized_state", "optimized_state"),
        ("calc_class", "calc_class", "calc_class"),
        ("esd_type", "esd_type", "esd_type"),
    ]
    
    # HOMO/LUMO from orbitals
    def get_homo_lumo_original(r):
        orbitals = r.get("orbitals", {})
        if isinstance(orbitals, dict):
            for key, orb_df in orbitals.items():
                if isinstance(orb_df, pd.DataFrame) and not orb_df.empty:
                    occupied = orb_df[orb_df["OCC"] > 0]
                    if not occupied.empty:
                        return occupied["eV"].max()
        return None
    
    def get_lumo_original(r):
        orbitals = r.get("orbitals", {})
        if isinstance(orbitals, dict):
            for key, orb_df in orbitals.items():
                if isinstance(orb_df, pd.DataFrame) and not orb_df.empty:
                    virtual = orb_df[orb_df["OCC"] == 0]
                    if not virtual.empty:
                        return virtual["eV"].min()
        return None
    
    scalar_checks.extend([
        ("homo_eV", get_homo_lumo_original, "homo_energy"),
        ("lumo_eV", get_lumo_original, "lumo_energy"),
    ])
    
    print("\n--- SCALAR DATA ---")
    print(f"{'Field':<25} {'Original':>15} {'Modular':>15} {'Match':>10}")
    print("-" * 70)
    
    for field_name, orig_key, mod_key in scalar_checks:
        # Count original
        orig_count = 0
        for r in original_results:
            if callable(orig_key):
                val = orig_key(r)
            else:
                val = r.get(orig_key)
            if val is not None and (not isinstance(val, float) or not pd.isna(val)):
                orig_count += 1
        
        # Count modular
        mod_count = 0
        for r in modular_results:
            val = r.get(mod_key)
            if val is not None and (not isinstance(val, float) or not pd.isna(val)):
                mod_count += 1
        
        match = "✓" if orig_count == mod_count else "✗"
        print(f"{field_name:<25} {orig_count:>10}/{len(original_results):<4} {mod_count:>10}/{len(modular_results):<4} {match:>5}")
        comparisons.append((field_name, orig_count, mod_count, orig_count == mod_count))
    
    # --- NESTED DATA ---
    nested_checks = [
        ("n_atoms", lambda r: count_df(r.get("cart_coords")), "cart_coords"),
        ("n_orbitals", lambda r: sum(count_df(v) for v in r.get("orbitals", {}).values() if isinstance(r.get("orbitals"), dict)), "orbitals"),
        ("n_tddft_states", lambda r: count_df(r.get("tddft_states")), "tddft_states"),
        ("n_ir_peaks", lambda r: count_df(r.get("ir_spectrum")), "ir"),
        ("n_vibrations", lambda r: count_df(r.get("vibrations")), "vibrations"),
        ("n_bonds", lambda r: count_df(r.get("internal", {}).get("bonds")) if isinstance(r.get("internal"), dict) else 0, "bonds"),
        ("n_angles", lambda r: count_df(r.get("internal", {}).get("angles")) if isinstance(r.get("internal"), dict) else 0, "angles"),
        ("n_dihedrals", lambda r: count_df(r.get("internal", {}).get("dihedrals")) if isinstance(r.get("internal"), dict) else 0, "dihedrals"),
        ("n_elec_dipole", lambda r: count_df(r.get("electric_dipole_spectrum", {}).get("abs")) if isinstance(r.get("electric_dipole_spectrum"), dict) else 0, "electric_dipole_abs"),
        ("n_vel_dipole", lambda r: count_df(r.get("velocity_dipole_spectrum", {}).get("abs")) if isinstance(r.get("velocity_dipole_spectrum"), dict) else 0, "velocity_dipole_abs"),
        ("n_raman", lambda r: count_df(r.get("raman_spectrum")), "raman"),
        ("n_mulliken", lambda r: count_df(r.get("mulliken_charges")), "mulliken"),
    ]
    
    print("\n--- NESTED DATA (files with data / total items) ---")
    print(f"{'Field':<25} {'Original':>20} {'Modular':>20} {'Match':>10}")
    print("-" * 80)
    
    for field_name, orig_getter, mod_key in nested_checks:
        # Count original
        orig_files = 0
        orig_items = 0
        for r in original_results:
            count = orig_getter(r)
            if count > 0:
                orig_files += 1
                orig_items += count
        
        # Count modular
        mod_files = 0
        mod_items = 0
        for r in modular_results:
            data = r.get(mod_key)
            count = count_df(data)
            if count > 0:
                mod_files += 1
                mod_items += count
        
        match = "✓" if orig_files == mod_files else "✗"
        orig_str = f"{orig_files} files / {orig_items}"
        mod_str = f"{mod_files} files / {mod_items}"
        print(f"{field_name:<25} {orig_str:>20} {mod_str:>20} {match:>5}")
        comparisons.append((field_name, orig_files, mod_files, orig_files == mod_files))
    
    # Summary
    print("\n" + "=" * 80)
    matches = sum(1 for _, _, _, m in comparisons if m)
    total = len(comparisons)
    print(f"SUMMARY: {matches}/{total} fields match ({100*matches/total:.0f}%)")
    
    mismatches = [(f, o, m) for f, o, m, match in comparisons if not match]
    if mismatches:
        print("\nMISMATCHES:")
        for field, orig, mod in mismatches:
            print(f"  {field}: Original={orig}, Modular={mod}")
    
    print("=" * 80)
    
    return original_results, modular_results


if __name__ == "__main__":
    orig, mod = main()
