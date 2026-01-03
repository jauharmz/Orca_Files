"""
Original Parser Test - Uses orca_praser.py directly for comparison

Uses rglob to find .out files and parses each with ORCAParser.

Run: python tests/test_original_parser.py
"""

import sys
from pathlib import Path
import json
import pandas as pd
from tqdm import tqdm

# Add root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))


def df_to_serializable(obj):
    """Convert DataFrames and other objects to JSON-serializable format."""
    if obj is None:
        return None
    if isinstance(obj, pd.DataFrame):
        if obj.empty:
            return None
        return obj.to_dict(orient='records')
    if isinstance(obj, dict):
        return {k: df_to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [df_to_serializable(v) for v in obj]
    if pd.isna(obj):
        return None
    return obj


def main():
    print("=" * 70)
    print("ORIGINAL PARSER TEST (orca_praser.py)")
    print("=" * 70)
    
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
    
    # 2. Find all .out files using rglob
    print("\n[STEP 2] Finding .out files...")
    all_files = list(data_dir.rglob("*.out"))
    print(f"  Found: {len(all_files)} .out files")
    
    if not all_files:
        print("ERROR: No .out files found!")
        return None
    
    # Show first few files
    for f in all_files[:5]:
        print(f"    - {f}")
    if len(all_files) > 5:
        print(f"    ... and {len(all_files) - 5} more")
    
    # 3. Parse each file using ORCAParser
    print("\n[STEP 3] Parsing with ORCAParser...")
    from orca_praser import ORCAParser
    
    molecules = []
    errors = []
    
    for file_path in tqdm(all_files, desc="Parsing ORCA files"):
        file_path_str = str(file_path)
        
        try:
            parser = ORCAParser(file_path_str)
            result = parser.parse(as_df=True)
            
            # Add filename to result
            result["source_file"] = file_path.name
            result["source_path"] = file_path_str
            
            # Extract molecule_id from filename
            base_name = file_path.stem  # filename without extension
            # Remove common suffixes
            import re
            mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt)$', '', base_name, flags=re.I)
            result["molecule_id"] = mol_id.strip("_-")
            
            molecules.append(result)
            
        except Exception as e:
            errors.append({"file": file_path.name, "error": str(e)})
            print(f"\n  Error parsing {file_path.name}: {e}")
    
    print(f"\n  Successfully parsed: {len(molecules)} files")
    print(f"  Errors: {len(errors)} files")
    
    if not molecules:
        print("ERROR: No molecules parsed!")
        return None
    
    # 4. Convert to DataFrame
    print("\n[STEP 4] Converting to DataFrame...")
    df = pd.DataFrame(molecules)
    print(f"  DataFrame shape: {df.shape}")
    print(f"  Columns: {list(df.columns)}")
    
    # 5. Export scalar CSV
    print("\n[STEP 5] Exporting scalar data to CSV...")
    csv_path = "original_parsed_data.csv"
    
    export_df = pd.DataFrame()
    export_df["molecule_id"] = df["molecule_id"]
    export_df["source_file"] = df["source_file"]
    
    # SMILES
    if "smiles" in df.columns:
        export_df["smiles"] = df["smiles"]
    
    # Geometry info
    if "geometry" in df.columns:
        export_df["charge"] = df["geometry"].apply(lambda x: x.get("charge") if isinstance(x, dict) else None)
        export_df["multiplicity"] = df["geometry"].apply(lambda x: x.get("multiplicity") if isinstance(x, dict) else None)
        export_df["mol_filename"] = df["geometry"].apply(lambda x: x.get("filename") if isinstance(x, dict) else None)
    
    # Energies (direct columns from parse())
    if "gibbs_energy_Eh" in df.columns:
        export_df["gibbs_Eh"] = df["gibbs_energy_Eh"]
    if "single_point_energy_Eh" in df.columns:
        export_df["single_point_Eh"] = df["single_point_energy_Eh"]
    
    # Get HOMO/LUMO from orbitals
    if "orbitals" in df.columns:
        def get_homo_lumo(orbitals):
            if not isinstance(orbitals, dict):
                return None, None, None
            
            for key, orb_df in orbitals.items():
                if isinstance(orb_df, pd.DataFrame) and not orb_df.empty:
                    occupied = orb_df[orb_df["OCC"] > 0]
                    virtual = orb_df[orb_df["OCC"] == 0]
                    
                    homo = occupied["eV"].max() if not occupied.empty else None
                    lumo = virtual["eV"].min() if not virtual.empty else None
                    gap = lumo - homo if homo is not None and lumo is not None else None
                    return homo, lumo, gap
            return None, None, None
        
        hl = df["orbitals"].apply(get_homo_lumo)
        export_df["homo_eV"] = [x[0] for x in hl]
        export_df["lumo_eV"] = [x[1] for x in hl]
        export_df["homo_lumo_gap"] = [x[2] for x in hl]
    
    # Calc info
    if "is_optimization" in df.columns:
        export_df["is_optimization"] = df["is_optimization"]
    if "optimized_state" in df.columns:
        export_df["optimized_state"] = df["optimized_state"]
    if "calc_class" in df.columns:
        export_df["calc_class"] = df["calc_class"]
    if "esd_type" in df.columns:
        export_df["esd_type"] = df["esd_type"]
    
    # Count nested data
    def count_df(x):
        if isinstance(x, pd.DataFrame):
            return len(x) if not x.empty else 0
        if isinstance(x, list):
            return len(x)
        return 0
    
    if "cart_coords" in df.columns:
        export_df["n_atoms"] = df["cart_coords"].apply(count_df)
    if "tddft_states" in df.columns:
        export_df["n_tddft_states"] = df["tddft_states"].apply(count_df)
    if "ir_spectrum" in df.columns:
        export_df["n_ir_peaks"] = df["ir_spectrum"].apply(count_df)
    if "vibrations" in df.columns:
        export_df["n_vibrations"] = df["vibrations"].apply(count_df)
    
    # Internal coords
    if "internal" in df.columns:
        def count_internal(internal):
            if not isinstance(internal, dict):
                return 0, 0, 0
            bonds = internal.get("bonds")
            angles = internal.get("angles")
            dihedrals = internal.get("dihedrals")
            return count_df(bonds), count_df(angles), count_df(dihedrals)
        
        ic = df["internal"].apply(count_internal)
        export_df["n_bonds"] = [x[0] for x in ic]
        export_df["n_angles"] = [x[1] for x in ic]
        export_df["n_dihedrals"] = [x[2] for x in ic]
    
    export_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")
    
    # 6. Export full JSON
    print("\n[STEP 6] Exporting ALL data to JSON...")
    json_path = "original_parsed_data.json"
    
    json_data = []
    for mol in molecules:
        record = {}
        for key, val in mol.items():
            record[key] = df_to_serializable(val)
        json_data.append(record)
    
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"  Saved: {json_path}")
    
    # 7. Summary
    print("\n" + "=" * 70)
    print("ORIGINAL PARSER TEST COMPLETE!")
    print("=" * 70)
    print(f"  Files found: {len(all_files)}")
    print(f"  Files parsed: {len(molecules)}")
    print(f"  Errors: {len(errors)}")
    print(f"  CSV: {csv_path}")
    print(f"  JSON: {json_path}")
    print("=" * 70)
    
    # 8. Show sample
    print("\nSample data (first 5):")
    display_cols = [c for c in ["molecule_id", "smiles", "gibbs_Eh", "homo_eV", "lumo_eV", "n_atoms", "calc_class"] 
                   if c in export_df.columns]
    print(export_df[display_cols].head().to_string())
    
    # 9. Show column stats
    print("\n--- Column Statistics ---")
    for col in export_df.columns:
        count = export_df[col].notna().sum()
        pct = 100 * count / len(export_df)
        print(f"  {col}: {count}/{len(export_df)} ({pct:.0f}%)")
    
    return df


if __name__ == "__main__":
    df = main()
