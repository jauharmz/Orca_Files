"""
Parser Test with CSV Export - Downloads from HuggingFace

Run: python tests/test_parser_csv.py
"""

import sys
from pathlib import Path
import logging

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Setup logging
from src.logger import setup_logging, get_logger
setup_logging(level=logging.INFO)
logger = get_logger("test_parser_csv")


def main():
    logger.info("=" * 60)
    logger.info("ORCA PARSER TEST - HuggingFace → CSV")
    logger.info("=" * 60)
    
    # 1. Download from HuggingFace
    logger.info("[1] Downloading from HuggingFace...")
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    if not data_dir.exists():
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
    logger.info(f"    Data dir: {data_dir}")
    
    # 2. Find .out files
    logger.info("[2] Finding .out files...")
    out_files = list(data_dir.rglob("*.out"))
    logger.info(f"    Found: {len(out_files)} files")
    
    if not out_files:
        logger.error("No .out files found!")
        return
    
    # Show first 10 files
    logger.info("    Files:")
    for f in out_files[:10]:
        logger.info(f"      - {f.name}")
    if len(out_files) > 10:
        logger.info(f"      ... and {len(out_files) - 10} more")
    
    # 3. Parse all files
    logger.info("[3] Parsing files...")
    from src.parser.batch import BatchParser
    
    parser = BatchParser()
    df = parser.parse_files([str(f) for f in out_files])
    logger.info(f"    Parsed: {len(df)} molecules")
    
    # 4. Show parsed data summary
    logger.info("[4] Parsed Data Summary:")
    logger.info(f"    Columns: {list(df.columns)}")
    
    # Count non-null values
    logger.info("    Non-null values per column:")
    for col in ["molecule_id", "smiles", "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy"]:
        if col in df.columns:
            count = df[col].notna().sum()
            logger.info(f"      {col}: {count}/{len(df)}")
    
    # 5. Export to CSV (scalar columns only)
    logger.info("[5] Exporting to CSV...")
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
    logger.info(f"    Saved: {output_csv}")
    
    # 6. Show sample data
    logger.info("[6] Sample Data:")
    print(df[available_cols].head(10).to_string())
    
    logger.info("=" * 60)
    logger.info(f"✅ DONE! CSV exported to: {output_csv}")
    logger.info("=" * 60)
    
    return df


if __name__ == "__main__":
    df = main()
