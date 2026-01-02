"""
Parser Test with CSV Export - Downloads from HuggingFace

Run: python tests/test_parser_csv.py
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def main():
    print("=" * 60)
    print("ORCA PARSER TEST - HuggingFace → CSV")
    print("=" * 60)
    
    # 1. Download from HuggingFace
    print("\n[1] Downloading from HuggingFace...")
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    if not data_dir.exists():
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
    print(f"    Data dir: {data_dir}")
    
    # 2. Find .out files
    print("\n[2] Finding .out files...")
    out_files = list(data_dir.rglob("*.out"))
    print(f"    Found: {len(out_files)} files")
    
    if not out_files:
        print("ERROR: No .out files found!")
        return
    
    # Show first 10 files
    print("    Files:")
    for f in out_files[:10]:
        print(f"      - {f.name}")
    if len(out_files) > 10:
        print(f"      ... and {len(out_files) - 10} more")
    
    # 3. Parse all files
    print("\n[3] Parsing files...")
    from src.parser.batch import BatchParser
    
    parser = BatchParser()
    df = parser.parse_files([str(f) for f in out_files])
    print(f"    Parsed: {len(df)} molecules")
    
    # 4. Show parsed data summary
    print("\n[4] Parsed Data Summary:")
    print(f"    Columns: {list(df.columns)}")
    
    # Count non-null values
    print("\n    Non-null values per column:")
    for col in ["molecule_id", "smiles", "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy"]:
        if col in df.columns:
            count = df[col].notna().sum()
            print(f"      {col}: {count}/{len(df)}")
    
    # 5. Export to CSV (scalar columns only)
    print("\n[5] Exporting to CSV...")
    scalar_cols = [
        "filename", "molecule_id", "smiles", 
        "charge", "multiplicity",
        "gibbs_Eh", "single_point_Eh",
        "homo_energy", "lumo_energy", "homo_lumo_gap",
        "is_optimization", "has_tddft", "optimized_state"
    ]
    available_cols = [c for c in scalar_cols if c in df.columns]
    
    output_csv = "parsed_data.csv"
    df[available_cols].to_csv(output_csv, index=False)
    print(f"    Saved: {output_csv}")
    
    # 6. Show sample data
    print("\n[6] Sample Data:")
    print(df[available_cols].head(10).to_string())
    
    print("\n" + "=" * 60)
    print(f"✅ DONE! CSV exported to: {output_csv}")
    print("=" * 60)
    
    return df


if __name__ == "__main__":
    df = main()
