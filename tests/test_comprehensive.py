"""
Comprehensive Parser Test - Enhanced Export with Multiple CSVs

Creates:
- parsed_molecules.csv        - Main molecule data (scalar fields)
- parsed_cart_coords.csv      - All atom coordinates
- parsed_orbitals.csv         - All orbital energies
- parsed_tddft_states.csv     - All TD-DFT transitions
- parsed_ir_spectra.csv       - All IR peaks
- parsed_internal_bonds.csv   - All bond data
- parsed_internal_angles.csv  - All angle data
- parsed_data.json            - Complete data

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
    
    # 2. Parse with MODULAR parser (detailed mode - one row per file with state)
    logger.info("\n[STEP 2] Parsing with MODULAR parsers (detailed mode)...")
    from src.parser.batch import BatchParser
    
    pattern = str(data_dir / "**/*.out")
    batch = BatchParser(pattern)
    df_modular = batch.parse_all_detailed(verbose=True)  # Use detailed mode
    
    # ========================================
    # EXPORT ALL DATA TO MULTIPLE CSVs
    # ========================================
    logger.info("\n[STEP 3] Exporting to multiple CSVs...")
    
    # 3A. Main molecules CSV (scalar data)
    molecules_df = pd.DataFrame()
    molecules_df["molecule_id"] = df_modular["molecule_id"]
    molecules_df["optimized_state"] = df_modular["optimized_state"]  # Include state
    molecules_df["filename"] = df_modular["filename"]
    
    # Add all scalar columns
    for col in ["smiles", "charge", "multiplicity", "gibbs_Eh", "single_point_Eh",
                "homo_energy", "lumo_energy", "homo_lumo_gap", 
                "is_optimization", "has_tddft", "optimized_state", "esd_type", "calc_class"]:
        if col in df_modular.columns:
            molecules_df[col] = df_modular[col]
    
    # Add counts for nested data
    nested_cols = ['cart_coords', 'orbitals', 'vibrations', 'ir', 'raman', 
                   'nmr_shielding', 'nmr_coupling', 'tddft_states', 
                   'electric_dipole_abs', 'electric_dipole_soc',
                   'velocity_dipole_abs', 'velocity_dipole_soc',
                   'bonds', 'angles', 'dihedrals']
    
    for col in nested_cols:
        if col in df_modular.columns:
            molecules_df[f'n_{col}'] = df_modular[col].apply(count_nested_data)
    
    molecules_df.to_csv("parsed_molecules.csv", index=False)
    logger.info(f"  Saved: parsed_molecules.csv ({len(molecules_df)} molecules, {len(molecules_df.columns)} columns)")
    
    # Helper to extract nested data with molecule_id and state
    def export_nested(df, col_name, csv_name, label):
        all_data = []
        for _, row in df.iterrows():
            mol_id = row["molecule_id"]
            state = row.get("optimized_state", "unknown")
            data = row.get(col_name)
            if isinstance(data, pd.DataFrame) and not data.empty:
                data_copy = data.copy()
                data_copy.insert(0, "molecule_id", mol_id)
                data_copy.insert(1, "optimized_state", state)
                all_data.append(data_copy)
        if all_data:
            result_df = pd.concat(all_data, ignore_index=True)
            result_df.to_csv(csv_name, index=False)
            logger.info(f"  Saved: {csv_name} ({len(result_df)} {label})")
            return all_data, result_df
        return all_data, None
    
    # 3B. Cartesian coordinates CSV
    all_coords, coords_df = export_nested(df_modular, "cart_coords", "parsed_cart_coords.csv", "atoms")
    
    # 3C. Orbitals CSV
    all_orbitals, orbitals_df = export_nested(df_modular, "orbitals", "parsed_orbitals.csv", "orbitals")
    
    # 3D. TD-DFT states CSV
    all_tddft, tddft_df = export_nested(df_modular, "tddft_states", "parsed_tddft_states.csv", "transitions")
    
    # 3E. IR spectra CSV
    all_ir, ir_df = export_nested(df_modular, "ir", "parsed_ir_spectra.csv", "peaks")
    
    # 3E2. Vibrations CSV
    all_vibs, vibs_df = export_nested(df_modular, "vibrations", "parsed_vibrations.csv", "modes")
    
    # 3E3. Raman CSV
    all_raman, raman_df = export_nested(df_modular, "raman", "parsed_raman.csv", "peaks")
    
    # 3E4. Mulliken CSV
    all_mulliken, mull_df = export_nested(df_modular, "mulliken", "parsed_mulliken.csv", "atoms")
    
    # 3E5. NMR Shielding CSV
    all_nmr_shield, nmr_shield_df = export_nested(df_modular, "nmr_shielding", "parsed_nmr_shielding.csv", "nuclei")
    
    # 3E6. NMR Coupling CSV
    all_nmr_coup, nmr_coup_df = export_nested(df_modular, "nmr_coupling", "parsed_nmr_coupling.csv", "couplings")
    
    # 3F. Internal coordinates - Bonds CSV
    all_bonds, bonds_df = export_nested(df_modular, "bonds", "parsed_internal_bonds.csv", "bonds")
    
    # 3G. Internal coordinates - Angles CSV
    all_angles, angles_df = export_nested(df_modular, "angles", "parsed_internal_angles.csv", "angles")
    
    # 3H. Internal coordinates - Dihedrals CSV
    all_dihedrals, dihedrals_df = export_nested(df_modular, "dihedrals", "parsed_internal_dihedrals.csv", "dihedrals")
    
    # 3I. Electric dipole absorption CSV
    all_elec_abs, elec_abs_df = export_nested(df_modular, "electric_dipole_abs", "parsed_electric_dipole_abs.csv", "transitions")
    
    # 3J. Electric dipole SOC CSV
    all_elec_soc, elec_soc_df = export_nested(df_modular, "electric_dipole_soc", "parsed_electric_dipole_soc.csv", "transitions")
    
    # 3K. Velocity dipole absorption CSV
    all_vel_abs, vel_abs_df = export_nested(df_modular, "velocity_dipole_abs", "parsed_velocity_dipole_abs.csv", "transitions")
    
    # 3L. Velocity dipole SOC CSV
    all_vel_soc, vel_soc_df = export_nested(df_modular, "velocity_dipole_soc", "parsed_velocity_dipole_soc.csv", "transitions")
    
    
    # 4. Create JSON export with ALL data
    logger.info("\n[STEP 4] Creating JSON export with ALL data...")
    
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
    
    # 5. Summary report
    logger.info("\n" + "=" * 70)
    logger.info("PARSING COMPLETE!")
    logger.info("=" * 70)
    
    report_lines = [
        "=" * 70,
        "ORCA PARSER OUTPUT SUMMARY",
        "=" * 70,
        "",
        f"Total molecules parsed: {len(df_modular)}",
        "",
        "--- CSV FILES CREATED ---",
        f"  parsed_molecules.csv       : {len(molecules_df)} molecules, {len(molecules_df.columns)} columns",
    ]
    
    if all_coords:
        report_lines.append(f"  parsed_cart_coords.csv     : {len(coords_df)} atoms")
    if all_orbitals:
        report_lines.append(f"  parsed_orbitals.csv        : {len(orbitals_df)} orbitals")
    if all_tddft:
        report_lines.append(f"  parsed_tddft_states.csv    : {len(tddft_df)} transitions")
    if all_ir:
        report_lines.append(f"  parsed_ir_spectra.csv      : {len(ir_df)} peaks")
    if all_bonds:
        report_lines.append(f"  parsed_internal_bonds.csv  : {len(bonds_df)} bonds")
    if all_angles:
        report_lines.append(f"  parsed_internal_angles.csv : {len(angles_df)} angles")
    if all_dihedrals:
        report_lines.append(f"  parsed_internal_dihedrals.csv: {len(dihedrals_df)} dihedrals")
    
    report_lines.extend([
        "",
        "--- SCALAR DATA COVERAGE ---"
    ])
    
    for col in ["smiles", "gibbs_Eh", "single_point_Eh", "homo_energy", "lumo_energy"]:
        if col in molecules_df.columns:
            count = molecules_df[col].notna().sum()
            pct = 100 * count / len(molecules_df)
            report_lines.append(f"  {col}: {count}/{len(molecules_df)} ({pct:.0f}%)")
    
    report_lines.append("")
    report_lines.append("=" * 70)
    
    report = "\n".join(report_lines)
    print("\n" + report)
    
    with open("parsing_summary.txt", "w") as f:
        f.write(report)
    
    logger.info("  Saved: parsing_summary.txt")
    
    # 6. Show sample from main CSV
    print("\n--- Sample from parsed_molecules.csv ---")
    display_cols = [c for c in ["molecule_id", "optimized_state", "gibbs_Eh", "homo_energy", "lumo_energy", 
                                 "calc_class", "n_cart_coords", "n_tddft_states"] 
                   if c in molecules_df.columns]
    print(molecules_df[display_cols].head(10).to_string())
    
    return df_modular


if __name__ == "__main__":
    df = main()
