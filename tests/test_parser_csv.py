"""
Parser Test with CSV Export - Uses MODULAR parsers with DETAILED LOGGING

Run: python tests/test_parser_csv.py

For DEBUG level output:
    import logging; logging.getLogger().setLevel(logging.DEBUG)
"""

import sys
from pathlib import Path
import logging
import pandas as pd

# Add root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

# Setup logging - INFO shows summary, DEBUG shows all details
LOG_LEVEL = logging.INFO  # Change to logging.DEBUG for maximum detail

logging.basicConfig(
    level=LOG_LEVEL,
    format='%(asctime)s | %(name)-20s | %(levelname)-7s | %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("orca_viz.log", mode='w')  # Fresh log each run
    ]
)
logger = logging.getLogger("test_parser_csv")


def main():
    logger.info("=" * 70)
    logger.info("ORCA MODULAR PARSER TEST")
    logger.info("=" * 70)
    
    # 1. Download from HuggingFace
    logger.info("")
    logger.info("[STEP 1] Downloading from HuggingFace...")
    from huggingface_hub import snapshot_download
    
    data_dir = Path("./test_data_hf")
    if not data_dir.exists():
        snapshot_download(
            repo_id="JauharMz/Orca",
            repo_type="dataset",
            local_dir=str(data_dir),
            local_dir_use_symlinks=False
        )
    logger.info(f"  Data directory: {data_dir}")
    
    # 2. Find .out files
    logger.info("")
    logger.info("[STEP 2] Finding ORCA output files...")
    out_files = list(data_dir.rglob("*.out"))
    logger.info(f"  Found: {len(out_files)} .out files")
    
    if not out_files:
        logger.error("No .out files found!")
        return None
    
    # Show files
    for f in out_files[:5]:
        logger.info(f"    - {f.name}")
    if len(out_files) > 5:
        logger.info(f"    ... and {len(out_files)-5} more")
    
    # 3. Parse using MODULAR parsers
    logger.info("")
    logger.info("[STEP 3] Parsing with MODULAR parsers...")
    logger.info("  (See detailed output below)")
    logger.info("")
    
    from src.parser.batch import BatchParser
    
    pattern = str(data_dir / "**/*.out")
    batch = BatchParser(pattern)
    
    df = batch.parse_all(verbose=True)
    
    if df.empty:
        logger.error("No data parsed!")
        return None
    
    # 4. Column summary
    logger.info("")
    logger.info("[STEP 4] Parsed data columns...")
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
    
    logger.info(f"  Scalar columns ({len(scalar_cols)}):")
    for col in scalar_cols:
        count = df[col].notna().sum()
        logger.info(f"    {col}: {count}/{len(df)}")
    
    logger.info(f"  Nested data ({len(nested_cols)}):")
    for col, count in nested_cols:
        logger.info(f"    {col}: {count}/{len(df)} with data")
    
    # 5. Export to CSV
    logger.info("")
    logger.info("[STEP 5] Exporting to CSV...")
    output_csv = "parsed_data.csv"
    df[scalar_cols].to_csv(output_csv, index=False)
    logger.info(f"  Saved: {output_csv}")
    
    # 6. Sample data display
    logger.info("")
    logger.info("[STEP 6] Sample data (first 5 rows):")
    display_cols = [c for c in ["molecule_id", "smiles", "gibbs_Eh", "single_point_Eh", 
                                 "homo_energy", "lumo_energy", "calc_class"] 
                   if c in scalar_cols]
    print("\n" + df[display_cols].head().to_string() + "\n")
    
    # 7. Final summary
    logger.info("=" * 70)
    logger.info("PARSING COMPLETE!")
    logger.info(f"  Molecules parsed: {len(df)}")
    logger.info(f"  CSV exported: {output_csv}")
    logger.info(f"  Log file: orca_viz.log")
    logger.info("")
    logger.info("TIP: For maximum detail, change LOG_LEVEL to logging.DEBUG")
    logger.info("=" * 70)
    
    return df


if __name__ == "__main__":
    df = main()
