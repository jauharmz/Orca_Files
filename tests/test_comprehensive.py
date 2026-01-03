"""
Comprehensive Parser Test - Enhanced Export + JSON + Comparison

Creates:
- parsed_data.csv - scalar columns only
- parsed_data_full.csv - with nested data counts
- parsed_data.json - ALL data including nested DataFrames
- comparison_report.txt - comparison with original parser

Run: python tests/test_comprehensive.py
"""

import sys
from pathlib import Path
import logging
import json
import pandas as pd

# Add root to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("orca_viz.log", mode='w')
    ]
)
logger = logging.getLogger("test_comprehensive")


def count_nested_data(val):
    """Count items in nested data."""
    if val is None:
        return 0
    if isinstance(val, pd.DataFrame):
        return len(val) if not val.empty else 0
    if isinstance(val, (list, dict)):
        return len(val) if val else 0
    return 0


def df_to_dict(df):
    """Convert DataFrame to serializable dict."""
    if df is None:
        return None
    if isinstance(df, pd.DataFrame):
        if df.empty:
            return None
        return df.to_dict(orient='records')
    return df


def main():
    logger.info("=" * 70)
    logger.info("COMPREHENSIVE PARSER TEST")
    logger.info("=" * 70)
    
    # 1. Download from HuggingFace
    logger.info("\n[STEP 1] Loading data from HuggingFace...")
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
    logger.info(f"  Found: {len(out_files)} files")
    
    # 2. Parse with MODULAR parser
    logger.info("\n[STEP 2] Parsing with MODULAR parsers...")
    from src.parser.batch import BatchParser
    
    pattern = str(data_dir / "**/*.out")
    batch = BatchParser(pattern)
    df_modular = batch.parse_all(verbose=True)
    
    # 3. Parse with ORIGINAL parser for comparison
    logger.info("\n[STEP 3] Parsing with ORIGINAL parser for comparison...")
    try:
        from orca_praser import ORCABatchParser
        original_batch = ORCABatchParser(pattern)
        df_original = original_batch.parse_all(verbose=False)
        has_original = True
        logger.info(f"  Original parser: {len(df_original)} rows")
    except Exception as e:
        logger.warning(f"  Original parser not available: {e}")
        df_original = None
        has_original = False
    
    # 4. Create ENHANCED CSV with nested data counts
    logger.info("\n[STEP 4] Creating enhanced CSV...")
    
    # Add count columns for nested data
    nested_cols = ['cart_coords', 'orbitals', 'vibrations', 'ir', 'raman', 
                   'nmr_shielding', 'nmr_coupling', 'tddft_states', 'mulliken',
                   'bonds', 'angles', 'dihedrals']
    
    for col in nested_cols:
        if col in df_modular.columns:
            df_modular[f'{col}_count'] = df_modular[col].apply(count_nested_data)
    
    # Identify all scalar columns
    scalar_cols = []
    for col in df_modular.columns:
        sample = df_modular[col].dropna()
        if sample.empty:
            scalar_cols.append(col)
        elif not isinstance(sample.iloc[0], (pd.DataFrame, dict, list)):
            scalar_cols.append(col)
    
    # Export enhanced CSV
    df_modular[scalar_cols].to_csv("parsed_data_full.csv", index=False)
    logger.info(f"  Saved: parsed_data_full.csv ({len(scalar_cols)} columns)")
    
    # 5. Create JSON export with ALL data
    logger.info("\n[STEP 5] Creating JSON export with ALL data...")
    
    json_data = []
    for idx, row in df_modular.iterrows():
        record = {}
        for col in df_modular.columns:
            val = row[col]
            if isinstance(val, pd.DataFrame):
                record[col] = df_to_dict(val)
            elif pd.isna(val):
                record[col] = None
            else:
                record[col] = val
        json_data.append(record)
    
    with open("parsed_data.json", "w") as f:
        json.dump(json_data, f, indent=2, default=str)
    logger.info(f"  Saved: parsed_data.json ({len(json_data)} records)")
    
    # 6. Comparison report
    logger.info("\n[STEP 6] Creating comparison report...")
    
    report_lines = [
        "=" * 70,
        "ORCA PARSER COMPARISON REPORT",
        "=" * 70,
        "",
        f"Files parsed: {len(df_modular)}",
        "",
        "--- MODULAR PARSER RESULTS ---",
    ]
    
    # Modular parser stats
    for col in ['smiles', 'gibbs_Eh', 'single_point_Eh', 'homo_energy', 'lumo_energy']:
        if col in df_modular.columns:
            count = df_modular[col].notna().sum()
            report_lines.append(f"  {col}: {count}/{len(df_modular)} ({100*count/len(df_modular):.1f}%)")
    
    # Nested data stats
    report_lines.append("")
    report_lines.append("--- NESTED DATA (from modular parser) ---")
    for col in nested_cols:
        count_col = f'{col}_count'
        if count_col in df_modular.columns:
            with_data = (df_modular[count_col] > 0).sum()
            total_items = df_modular[count_col].sum()
            report_lines.append(f"  {col}: {with_data}/{len(df_modular)} files, {int(total_items)} total items")
    
    # Comparison with original
    if has_original and df_original is not None:
        report_lines.append("")
        report_lines.append("--- ORIGINAL PARSER RESULTS ---")
        for col in ['smiles', 'gibbs_Eh', 'single_point_Eh', 'homo_eV', 'lumo_eV']:
            if col in df_original.columns:
                count = df_original[col].notna().sum()
                report_lines.append(f"  {col}: {count}/{len(df_original)} ({100*count/len(df_original):.1f}%)")
        
        report_lines.append("")
        report_lines.append("--- COMPARISON ---")
        
        # Compare key columns
        if 'homo_eV' in df_original.columns and 'homo_energy' in df_modular.columns:
            # Match by molecule_id
            orig_homos = df_original.set_index('molecule_id')['homo_eV'].dropna()
            mod_homos = df_modular.set_index('molecule_id')['homo_energy'].dropna()
            common = orig_homos.index.intersection(mod_homos.index)
            if len(common) > 0:
                diff = (orig_homos.loc[common] - mod_homos.loc[common]).abs()
                report_lines.append(f"  HOMO energy match: {len(common)} molecules, max diff: {diff.max():.6f} eV")
    
    report_lines.append("")
    report_lines.append("=" * 70)
    
    report = "\n".join(report_lines)
    with open("comparison_report.txt", "w") as f:
        f.write(report)
    print("\n" + report)
    
    # 7. Summary
    logger.info("\n" + "=" * 70)
    logger.info("COMPLETE! Files created:")
    logger.info("  - parsed_data_full.csv (enhanced with nested counts)")
    logger.info("  - parsed_data.json (ALL data)")
    logger.info("  - comparison_report.txt")
    logger.info("  - orca_viz.log")
    logger.info("=" * 70)
    
    return df_modular


if __name__ == "__main__":
    df = main()
