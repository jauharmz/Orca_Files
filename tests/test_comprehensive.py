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
    
    # 3B. Cartesian coordinates CSV
    all_coords = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        coords = row.get("cart_coords")
        if isinstance(coords, pd.DataFrame) and not coords.empty:
            coords_copy = coords.copy()
            coords_copy.insert(0, "molecule_id", mol_id)
            all_coords.append(coords_copy)
    
    if all_coords:
        coords_df = pd.concat(all_coords, ignore_index=True)
        coords_df.to_csv("parsed_cart_coords.csv", index=False)
        logger.info(f"  Saved: parsed_cart_coords.csv ({len(coords_df)} atoms)")
    
    # 3C. Orbitals CSV
    all_orbitals = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        orbitals = row.get("orbitals")
        if isinstance(orbitals, pd.DataFrame) and not orbitals.empty:
            orb_copy = orbitals.copy()
            orb_copy.insert(0, "molecule_id", mol_id)
            all_orbitals.append(orb_copy)
    
    if all_orbitals:
        orbitals_df = pd.concat(all_orbitals, ignore_index=True)
        orbitals_df.to_csv("parsed_orbitals.csv", index=False)
        logger.info(f"  Saved: parsed_orbitals.csv ({len(orbitals_df)} orbitals)")
    
    # 3D. TD-DFT states CSV
    all_tddft = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        tddft = row.get("tddft_states")
        if isinstance(tddft, pd.DataFrame) and not tddft.empty:
            tddft_copy = tddft.copy()
            tddft_copy.insert(0, "molecule_id", mol_id)
            all_tddft.append(tddft_copy)
    
    if all_tddft:
        tddft_df = pd.concat(all_tddft, ignore_index=True)
        tddft_df.to_csv("parsed_tddft_states.csv", index=False)
        logger.info(f"  Saved: parsed_tddft_states.csv ({len(tddft_df)} transitions)")
    
    # 3E. IR spectra CSV
    all_ir = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        ir = row.get("ir")
        if isinstance(ir, pd.DataFrame) and not ir.empty:
            ir_copy = ir.copy()
            ir_copy.insert(0, "molecule_id", mol_id)
            all_ir.append(ir_copy)
    
    if all_ir:
        ir_df = pd.concat(all_ir, ignore_index=True)
        ir_df.to_csv("parsed_ir_spectra.csv", index=False)
        logger.info(f"  Saved: parsed_ir_spectra.csv ({len(ir_df)} peaks)")
    
    # 3E2. Vibrations CSV
    all_vibs = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        vibs = row.get("vibrations")
        if isinstance(vibs, pd.DataFrame) and not vibs.empty:
            vibs_copy = vibs.copy()
            vibs_copy.insert(0, "molecule_id", mol_id)
            all_vibs.append(vibs_copy)
    
    if all_vibs:
        vibs_df = pd.concat(all_vibs, ignore_index=True)
        vibs_df.to_csv("parsed_vibrations.csv", index=False)
        logger.info(f"  Saved: parsed_vibrations.csv ({len(vibs_df)} modes)")
    
    # 3E3. Raman CSV
    all_raman = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        raman = row.get("raman")
        if isinstance(raman, pd.DataFrame) and not raman.empty:
            raman_copy = raman.copy()
            raman_copy.insert(0, "molecule_id", mol_id)
            all_raman.append(raman_copy)
    
    if all_raman:
        raman_df = pd.concat(all_raman, ignore_index=True)
        raman_df.to_csv("parsed_raman.csv", index=False)
        logger.info(f"  Saved: parsed_raman.csv ({len(raman_df)} peaks)")
    
    # 3E4. Mulliken CSV
    all_mulliken = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        mull = row.get("mulliken")
        if isinstance(mull, pd.DataFrame) and not mull.empty:
            mull_copy = mull.copy()
            mull_copy.insert(0, "molecule_id", mol_id)
            all_mulliken.append(mull_copy)
    
    if all_mulliken:
        mull_df = pd.concat(all_mulliken, ignore_index=True)
        mull_df.to_csv("parsed_mulliken.csv", index=False)
        logger.info(f"  Saved: parsed_mulliken.csv ({len(mull_df)} atoms)")
    
    # 3E5. NMR Shielding CSV
    all_nmr_shield = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        nmr = row.get("nmr_shielding")
        if isinstance(nmr, pd.DataFrame) and not nmr.empty:
            nmr_copy = nmr.copy()
            nmr_copy.insert(0, "molecule_id", mol_id)
            all_nmr_shield.append(nmr_copy)
    
    if all_nmr_shield:
        nmr_shield_df = pd.concat(all_nmr_shield, ignore_index=True)
        nmr_shield_df.to_csv("parsed_nmr_shielding.csv", index=False)
        logger.info(f"  Saved: parsed_nmr_shielding.csv ({len(nmr_shield_df)} nuclei)")
    
    # 3E6. NMR Coupling CSV
    all_nmr_coup = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        nmr = row.get("nmr_coupling")
        if isinstance(nmr, pd.DataFrame) and not nmr.empty:
            nmr_copy = nmr.copy()
            nmr_copy.insert(0, "molecule_id", mol_id)
            all_nmr_coup.append(nmr_copy)
    
    if all_nmr_coup:
        nmr_coup_df = pd.concat(all_nmr_coup, ignore_index=True)
        nmr_coup_df.to_csv("parsed_nmr_coupling.csv", index=False)
        logger.info(f"  Saved: parsed_nmr_coupling.csv ({len(nmr_coup_df)} couplings)")
    
    # 3F. Internal coordinates - Bonds CSV
    all_bonds = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        bonds = row.get("bonds")
        if isinstance(bonds, pd.DataFrame) and not bonds.empty:
            bonds_copy = bonds.copy()
            bonds_copy.insert(0, "molecule_id", mol_id)
            all_bonds.append(bonds_copy)
    
    if all_bonds:
        bonds_df = pd.concat(all_bonds, ignore_index=True)
        bonds_df.to_csv("parsed_internal_bonds.csv", index=False)
        logger.info(f"  Saved: parsed_internal_bonds.csv ({len(bonds_df)} bonds)")
    
    # 3G. Internal coordinates - Angles CSV
    all_angles = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        angles = row.get("angles")
        if isinstance(angles, pd.DataFrame) and not angles.empty:
            angles_copy = angles.copy()
            angles_copy.insert(0, "molecule_id", mol_id)
            all_angles.append(angles_copy)
    
    if all_angles:
        angles_df = pd.concat(all_angles, ignore_index=True)
        angles_df.to_csv("parsed_internal_angles.csv", index=False)
        logger.info(f"  Saved: parsed_internal_angles.csv ({len(angles_df)} angles)")
    
    # 3H. Internal coordinates - Dihedrals CSV
    all_dihedrals = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        dihedrals = row.get("dihedrals")
        if isinstance(dihedrals, pd.DataFrame) and not dihedrals.empty:
            dih_copy = dihedrals.copy()
            dih_copy.insert(0, "molecule_id", mol_id)
            all_dihedrals.append(dih_copy)
    
    if all_dihedrals:
        dihedrals_df = pd.concat(all_dihedrals, ignore_index=True)
        dihedrals_df.to_csv("parsed_internal_dihedrals.csv", index=False)
        logger.info(f"  Saved: parsed_internal_dihedrals.csv ({len(dihedrals_df)} dihedrals)")
    
    # 3I. Electric dipole absorption CSV
    all_elec_abs = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        elec = row.get("electric_dipole_abs")
        if isinstance(elec, pd.DataFrame) and not elec.empty:
            elec_copy = elec.copy()
            elec_copy.insert(0, "molecule_id", mol_id)
            all_elec_abs.append(elec_copy)
    
    if all_elec_abs:
        elec_abs_df = pd.concat(all_elec_abs, ignore_index=True)
        elec_abs_df.to_csv("parsed_electric_dipole.csv", index=False)
        logger.info(f"  Saved: parsed_electric_dipole.csv ({len(elec_abs_df)} transitions)")
    
    # 3J. Velocity dipole absorption CSV
    all_vel_abs = []
    for _, row in df_modular.iterrows():
        mol_id = row["molecule_id"]
        vel = row.get("velocity_dipole_abs")
        if isinstance(vel, pd.DataFrame) and not vel.empty:
            vel_copy = vel.copy()
            vel_copy.insert(0, "molecule_id", mol_id)
            all_vel_abs.append(vel_copy)
    
    if all_vel_abs:
        vel_abs_df = pd.concat(all_vel_abs, ignore_index=True)
        vel_abs_df.to_csv("parsed_velocity_dipole.csv", index=False)
        logger.info(f"  Saved: parsed_velocity_dipole.csv ({len(vel_abs_df)} transitions)")
    
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
    display_cols = [c for c in ["molecule_id", "gibbs_Eh", "homo_energy", "lumo_energy", 
                                 "calc_class", "n_cart_coords", "n_tddft_states"] 
                   if c in molecules_df.columns]
    print(molecules_df[display_cols].head(10).to_string())
    
    return df_modular


if __name__ == "__main__":
    df = main()
