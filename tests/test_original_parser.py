"""
Original Parser Test - Uses orca_praser.py directly for comparison

Creates multiple CSVs with ALL nested data:
- original_molecules.csv       - Main molecule data
- original_cart_coords.csv     - All atom coordinates  
- original_orbitals.csv        - All orbital energies
- original_tddft_states.csv    - All TD-DFT transitions
- original_ir_spectra.csv      - All IR peaks
- original_internal_bonds.csv  - All bonds
- original_data.json           - Complete data

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


def count_df(x):
    """Count items in DataFrame or list."""
    if isinstance(x, pd.DataFrame):
        return len(x) if not x.empty else 0
    if isinstance(x, list):
        return len(x)
    return 0


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
    
    # 3. Parse each file using ORCAParser
    print("\n[STEP 3] Parsing with ORCAParser...")
    from orca_praser import ORCAParser
    import re
    
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
            base_name = file_path.stem
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
    
    # Convert to DataFrame for easier processing
    df = pd.DataFrame(molecules)
    print(f"  DataFrame columns: {list(df.columns)}")
    
    # ========================================
    # EXPORT ALL DATA TO MULTIPLE CSVs
    # ========================================
    print("\n[STEP 4] Exporting to multiple CSVs...")
    
    # 4A. Main molecules CSV (scalar data)
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
    
    # Energies
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
    for col in ["is_optimization", "optimized_state", "calc_class", "esd_type"]:
        if col in df.columns:
            export_df[col] = df[col]
    
    # Count nested data
    if "cart_coords" in df.columns:
        export_df["n_atoms"] = df["cart_coords"].apply(count_df)
    if "tddft_states" in df.columns:
        export_df["n_tddft_states"] = df["tddft_states"].apply(count_df)
    if "ir_spectrum" in df.columns:
        export_df["n_ir_peaks"] = df["ir_spectrum"].apply(count_df)
    if "vibrations" in df.columns:
        export_df["n_vibrations"] = df["vibrations"].apply(count_df)
    if "internal" in df.columns:
        def count_internal(internal):
            if not isinstance(internal, dict):
                return 0, 0, 0
            return count_df(internal.get("bonds")), count_df(internal.get("angles")), count_df(internal.get("dihedrals"))
        ic = df["internal"].apply(count_internal)
        export_df["n_bonds"] = [x[0] for x in ic]
        export_df["n_angles"] = [x[1] for x in ic]
        export_df["n_dihedrals"] = [x[2] for x in ic]
    
    export_df.to_csv("original_molecules.csv", index=False)
    print(f"  Saved: original_molecules.csv ({len(export_df)} rows, {len(export_df.columns)} cols)")
    
    # 4B. Cartesian coordinates CSV
    all_coords = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        coords = mol.get("cart_coords")
        if isinstance(coords, pd.DataFrame) and not coords.empty:
            coords_copy = coords.copy()
            coords_copy.insert(0, "molecule_id", mol_id)
            all_coords.append(coords_copy)
    
    if all_coords:
        coords_df = pd.concat(all_coords, ignore_index=True)
        coords_df.to_csv("original_cart_coords.csv", index=False)
        print(f"  Saved: original_cart_coords.csv ({len(coords_df)} atoms)")
    
    # 4C. Orbitals CSV
    all_orbitals = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        orbitals = mol.get("orbitals")
        if isinstance(orbitals, dict):
            for state_key, orb_df in orbitals.items():
                if isinstance(orb_df, pd.DataFrame) and not orb_df.empty:
                    orb_copy = orb_df.copy()
                    orb_copy.insert(0, "molecule_id", mol_id)
                    orb_copy.insert(1, "state", state_key)
                    all_orbitals.append(orb_copy)
    
    if all_orbitals:
        orbitals_df = pd.concat(all_orbitals, ignore_index=True)
        orbitals_df.to_csv("original_orbitals.csv", index=False)
        print(f"  Saved: original_orbitals.csv ({len(orbitals_df)} orbitals)")
    
    # 4D. TD-DFT states CSV
    all_tddft = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        tddft = mol.get("tddft_states")
        if isinstance(tddft, pd.DataFrame) and not tddft.empty:
            tddft_copy = tddft.copy()
            tddft_copy.insert(0, "molecule_id", mol_id)
            all_tddft.append(tddft_copy)
    
    if all_tddft:
        tddft_df = pd.concat(all_tddft, ignore_index=True)
        tddft_df.to_csv("original_tddft_states.csv", index=False)
        print(f"  Saved: original_tddft_states.csv ({len(tddft_df)} transitions)")
    
    # 4E. IR spectra CSV
    all_ir = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        ir = mol.get("ir_spectrum")
        if isinstance(ir, list) and ir:
            ir_df = pd.DataFrame(ir)
            ir_df.insert(0, "molecule_id", mol_id)
            all_ir.append(ir_df)
    
    if all_ir:
        ir_df = pd.concat(all_ir, ignore_index=True)
        ir_df.to_csv("original_ir_spectra.csv", index=False)
        print(f"  Saved: original_ir_spectra.csv ({len(ir_df)} peaks)")
    
    # 4F. Internal coordinates - Bonds CSV
    all_bonds = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        internal = mol.get("internal")
        if isinstance(internal, dict):
            bonds = internal.get("bonds")
            if isinstance(bonds, pd.DataFrame) and not bonds.empty:
                bonds_copy = bonds.copy()
                bonds_copy.insert(0, "molecule_id", mol_id)
                all_bonds.append(bonds_copy)
    
    if all_bonds:
        bonds_df = pd.concat(all_bonds, ignore_index=True)
        bonds_df.to_csv("original_internal_bonds.csv", index=False)
        print(f"  Saved: original_internal_bonds.csv ({len(bonds_df)} bonds)")
    
    # 4G. Internal coordinates - Angles CSV
    all_angles = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        internal = mol.get("internal")
        if isinstance(internal, dict):
            angles = internal.get("angles")
            if isinstance(angles, pd.DataFrame) and not angles.empty:
                angles_copy = angles.copy()
                angles_copy.insert(0, "molecule_id", mol_id)
                all_angles.append(angles_copy)
    
    if all_angles:
        angles_df = pd.concat(all_angles, ignore_index=True)
        angles_df.to_csv("original_internal_angles.csv", index=False)
        print(f"  Saved: original_internal_angles.csv ({len(angles_df)} angles)")
    
    # 4H. Internal coordinates - Dihedrals CSV
    all_dihedrals = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        internal = mol.get("internal")
        if isinstance(internal, dict):
            dihedrals = internal.get("dihedrals")
            if isinstance(dihedrals, pd.DataFrame) and not dihedrals.empty:
                dih_copy = dihedrals.copy()
                dih_copy.insert(0, "molecule_id", mol_id)
                all_dihedrals.append(dih_copy)
    
    if all_dihedrals:
        dihedrals_df = pd.concat(all_dihedrals, ignore_index=True)
        dihedrals_df.to_csv("original_internal_dihedrals.csv", index=False)
        print(f"  Saved: original_internal_dihedrals.csv ({len(dihedrals_df)} dihedrals)")
    
    # 4I. Electric dipole spectrum CSV
    all_elec = []
    for mol in molecules:
        mol_id = mol["molecule_id"]
        elec = mol.get("electric_dipole_spectrum")
        if isinstance(elec, dict):
            for key in ["abs", "soc"]:
                if key in elec:
                    elec_df = elec[key]
                    if isinstance(elec_df, pd.DataFrame) and not elec_df.empty:
                        elec_copy = elec_df.copy()
                        elec_copy.insert(0, "molecule_id", mol_id)
                        elec_copy.insert(1, "type", key)
                        all_elec.append(elec_copy)
    
    if all_elec:
        elec_df = pd.concat(all_elec, ignore_index=True)
        elec_df.to_csv("original_electric_dipole.csv", index=False)
        print(f"  Saved: original_electric_dipole.csv ({len(elec_df)} transitions)")
    
    # 5. Export full JSON
    print("\n[STEP 5] Exporting ALL data to JSON...")
    json_data = []
    for mol in molecules:
        record = {}
        for key, val in mol.items():
            record[key] = df_to_serializable(val)
        json_data.append(record)
    
    with open("original_data.json", "w") as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"  Saved: original_data.json")
    
    # 6. Summary
    print("\n" + "=" * 70)
    print("ORIGINAL PARSER TEST COMPLETE!")
    print("=" * 70)
    
    report_lines = [
        f"  Files found: {len(all_files)}",
        f"  Files parsed: {len(molecules)}",
        f"  Errors: {len(errors)}",
        "",
        "  CSV FILES CREATED:",
        f"    original_molecules.csv      : {len(export_df)} rows",
    ]
    
    if all_coords:
        report_lines.append(f"    original_cart_coords.csv    : {len(coords_df)} atoms")
    if all_orbitals:
        report_lines.append(f"    original_orbitals.csv       : {len(orbitals_df)} orbitals")
    if all_tddft:
        report_lines.append(f"    original_tddft_states.csv   : {len(tddft_df)} transitions")
    if all_ir:
        report_lines.append(f"    original_ir_spectra.csv     : {len(ir_df)} peaks")
    if all_bonds:
        report_lines.append(f"    original_internal_bonds.csv : {len(bonds_df)} bonds")
    if all_angles:
        report_lines.append(f"    original_internal_angles.csv: {len(angles_df)} angles")
    if all_dihedrals:
        report_lines.append(f"    original_internal_dihedrals.csv: {len(dihedrals_df)} dihedrals")
    if all_elec:
        report_lines.append(f"    original_electric_dipole.csv: {len(elec_df)} transitions")
    
    for line in report_lines:
        print(line)
    
    print("=" * 70)
    
    # 7. Show sample from main CSV
    print("\n--- Sample from original_molecules.csv ---")
    display_cols = [c for c in ["molecule_id", "gibbs_Eh", "homo_eV", "lumo_eV", 
                                 "calc_class", "n_atoms", "n_tddft_states"] 
                   if c in export_df.columns]
    print(export_df[display_cols].head(10).to_string())
    
    return df


if __name__ == "__main__":
    df = main()
