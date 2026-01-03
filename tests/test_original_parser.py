"""
Original Parser Test - Uses orca_praser.py directly for comparison

Creates:
- original_parsed_data.csv - scalar columns
- original_parsed_data.json - ALL data
- original_parse_log.txt - detailed parsing info

Run: python tests/test_original_parser.py
"""

import sys
from pathlib import Path
import json
import pandas as pd

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
    
    out_files = list(data_dir.rglob("*.out"))
    print(f"  Found: {len(out_files)} .out files")
    
    # 2. Parse using ORIGINAL parser
    print("\n[STEP 2] Parsing with ORIGINAL orca_praser.py...")
    from orca_praser import ORCABatchParser
    
    pattern = str(data_dir / "**/*.out")
    batch = ORCABatchParser(pattern)
    df = batch.parse_all(verbose=True)
    
    print(f"\n  Parsed: {len(df)} molecule records")
    
    # 3. Analyze columns
    print("\n[STEP 3] Analyzing parsed data...")
    print(f"  Columns in DataFrame: {list(df.columns)}")
    
    scalar_cols = []
    nested_cols = []
    
    for col in df.columns:
        sample = df[col].dropna()
        if sample.empty:
            scalar_cols.append(col)
        elif isinstance(sample.iloc[0], (pd.DataFrame, dict, list)):
            count = sum(1 for x in df[col] if x is not None and (not isinstance(x, pd.DataFrame) or not x.empty))
            if count > 0:
                nested_cols.append((col, count))
        else:
            scalar_cols.append(col)
    
    print(f"\n  Scalar columns ({len(scalar_cols)}):")
    for col in scalar_cols:
        count = df[col].notna().sum()
        pct = 100 * count / len(df)
        print(f"    {col}: {count}/{len(df)} ({pct:.0f}%)")
    
    print(f"\n  Nested data ({len(nested_cols)}):")
    for col, count in nested_cols:
        pct = 100 * count / len(df)
        print(f"    {col}: {count}/{len(df)} ({pct:.0f}%)")
    
    # 4. Export scalar CSV
    print("\n[STEP 4] Exporting scalar data to CSV...")
    csv_path = "original_parsed_data.csv"
    
    # Flatten some common scalar fields
    export_df = pd.DataFrame()
    
    # Get molecule_id
    if "molecule_id" in df.columns:
        export_df["molecule_id"] = df["molecule_id"]
    else:
        # Try to extract from index or other sources
        export_df["molecule_id"] = df.index.astype(str)
    
    if "smiles" in df.columns:
        export_df["smiles"] = df["smiles"]
    
    # Get geometry info
    if "geometry" in df.columns:
        export_df["charge"] = df["geometry"].apply(lambda x: x.get("charge") if isinstance(x, dict) else None)
        export_df["multiplicity"] = df["geometry"].apply(lambda x: x.get("multiplicity") if isinstance(x, dict) else None)
        export_df["filename"] = df["geometry"].apply(lambda x: x.get("filename") if isinstance(x, dict) else None)
    
    # Get energies from energies dict
    if "energies" in df.columns:
        def get_gibbs(energies):
            if not isinstance(energies, dict):
                return None
            for key, val in energies.items():
                if isinstance(val, dict) and "gibbs_Eh" in val:
                    return val["gibbs_Eh"]
            return None
        
        def get_sp(energies):
            if not isinstance(energies, dict):
                return None
            for key, val in energies.items():
                if isinstance(val, dict) and "single_point_Eh" in val:
                    return val["single_point_Eh"]
            return None
        
        export_df["gibbs_Eh"] = df["energies"].apply(get_gibbs)
        export_df["single_point_Eh"] = df["energies"].apply(get_sp)
    
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
                    gap = lumo - homo if homo and lumo else None
                    return homo, lumo, gap
            return None, None, None
        
        hl = df["orbitals"].apply(get_homo_lumo)
        export_df["homo_eV"] = [x[0] for x in hl]
        export_df["lumo_eV"] = [x[1] for x in hl]
        export_df["homo_lumo_gap"] = [x[2] for x in hl]
    
    # Count nested data
    def count_df(x):
        if isinstance(x, pd.DataFrame):
            return len(x) if not x.empty else 0
        return 0
    
    if "cart_coords" in df.columns:
        export_df["n_atoms"] = df["cart_coords"].apply(count_df)
    if "tddft_states" in df.columns:
        export_df["n_tddft_states"] = df["tddft_states"].apply(count_df)
    if "ir_spectrum" in df.columns:
        export_df["n_ir_peaks"] = df["ir_spectrum"].apply(lambda x: len(x) if isinstance(x, list) else 0)
    if "vibrations" in df.columns:
        export_df["n_vibrations"] = df["vibrations"].apply(lambda x: len(x) if isinstance(x, list) else 0)
    
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
    
    # 5. Export full JSON
    print("\n[STEP 5] Exporting ALL data to JSON...")
    json_path = "original_parsed_data.json"
    
    json_data = []
    for idx, row in df.iterrows():
        record = {}
        for col in df.columns:
            record[col] = df_to_serializable(row[col])
        json_data.append(record)
    
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"  Saved: {json_path}")
    
    # 6. Summary
    print("\n" + "=" * 70)
    print("ORIGINAL PARSER TEST COMPLETE!")
    print("=" * 70)
    print(f"  Molecules: {len(df)}")
    print(f"  CSV: {csv_path}")
    print(f"  JSON: {json_path}")
    print("")
    print("Compare with modular parser:")
    print("  python tests/test_comprehensive.py")
    print("=" * 70)
    
    # 7. Show sample
    print("\nSample data (first 5):")
    display_cols = [c for c in ["molecule_id", "smiles", "gibbs_Eh", "homo_eV", "lumo_eV", "n_atoms"] 
                   if c in export_df.columns]
    print(export_df[display_cols].head().to_string())
    
    return df


if __name__ == "__main__":
    df = main()
