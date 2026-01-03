"""
Batch Parser - Parse multiple ORCA files with molecule grouping

Produces a DataFrame where each row is a standardized molecule record,
grouping data from multiple calculations (OPT, SP, TDDFT) for the same molecule.

Refactored from orca_praser.py ORCABatchParser
"""

import os
import re
import glob
from pathlib import Path
from typing import List, Optional, Dict, Any
import pandas as pd

from .factory import ParserFactory
from .spectrum_file import SpectrumFileParser
from ..logger import get_logger


def _is_empty_obj(x: Any) -> bool:
    """Check if object is empty."""
    if x is None:
        return True
    if isinstance(x, pd.DataFrame):
        return x.empty
    if isinstance(x, dict):
        return all(_is_empty_obj(v) for v in x.values()) if x else True
    if isinstance(x, list):
        return len(x) == 0
    return False


class BatchParser:
    """
    Batch parser for multiple ORCA output files.
    
    Groups calculations by molecule_id and produces standardized records.
    """
    
    def __init__(self, *patterns: str):
        """
        Initialize with glob patterns.
        
        Args:
            *patterns: Glob patterns like "*.out" or "data/**/*.out"
        """
        self.logger = get_logger("BatchParser")
        self.patterns = patterns
        self.files: List[str] = []
        self.states: Dict[str, int] = {}
        
        # Expand patterns to file list
        if not patterns:
            patterns = ("*.out",)
        
        for pattern in patterns:
            if os.path.isfile(pattern):
                self.files.append(pattern)
            else:
                self.files.extend(glob.glob(pattern, recursive=True))
        
        self.files = sorted(set(self.files))
        self.logger.info(f"Initialized with {len(patterns)} pattern(s), found {len(self.files)} files")
    
    def parse_files(self, filepaths: List[str]) -> pd.DataFrame:
        """Parse a list of file paths."""
        self.files = filepaths
        self.logger.info(f"Parsing {len(filepaths)} files")
        return self.parse_all(verbose=False)
    
    def parse_folder(self, folder: str, pattern: str = "**/*.out") -> pd.DataFrame:
        """Parse all matching files in folder."""
        folder_path = Path(folder)
        glob_pattern = str(folder_path / pattern)
        self.files = sorted(glob.glob(glob_pattern, recursive=True))
        self.logger.info(f"Parsing folder {folder}: found {len(self.files)} files")
        return self.parse_all(verbose=False)
    
    def parse_all(self, verbose: bool = True, limit_files: Optional[int] = None) -> pd.DataFrame:
        """
        Parse all files and return DataFrame with standardized molecule records.
        
        Groups calculations by molecule_id and consolidates data.
        """
        files = self.files[:limit_files] if limit_files else self.files
        
        if not files:
            self.logger.warning("No files to parse")
            return pd.DataFrame()
        
        self.logger.info(f"Starting batch parse of {len(files)} files")
        
        # Step 1: Parse all files individually
        factory = ParserFactory()
        spec_file_parser = SpectrumFileParser()
        parsed = []
        
        for i, filepath in enumerate(files):
            if verbose and (i + 1) % 5 == 0:
                self.logger.info(f"Progress: {i+1}/{len(files)} ({100*(i+1)/len(files):.0f}%)")
            
            try:
                result = factory.parse(filepath)
                row = result.to_dict()
                
                # Extract molecule_id
                row["molecule_id"] = self._extract_molecule_id(filepath, result.geometry.filename)
                row["_filepath"] = filepath
                row["_is_optimization"] = result.is_optimization
                row["_optimized_state"] = result.optimized_state
                
                # Try to find spectrum files
                spec_files = spec_file_parser.find_spectrum_files(
                    filepath, 
                    result.esd_type,
                    result.geometry.filename
                )
                if spec_files:
                    for spec_type, spec_path in spec_files.items():
                        row["spectrum_file"] = spec_file_parser.parse_spectrum_file(spec_path, spec_type)
                
                parsed.append(row)
            except Exception as e:
                if verbose:
                    self.logger.warning(f"Failed to parse {filepath}: {e}")
                continue
        
        if not parsed:
            self.logger.warning("No files parsed successfully")
            return pd.DataFrame()
        
        raw = pd.DataFrame(parsed)
        
        # Step 2: Propagate geometry-dependent fields from OPT to SP
        geometry_fields = [
            "cart_coords", "smiles", "vibrations", "ir", "raman",
            "nmr_shielding", "nmr_coupling", "mulliken", "bonds", "angles", "dihedrals"
        ]
        
        filled = raw.copy()
        
        for mid in raw["molecule_id"].unique():
            mask_mid = raw["molecule_id"] == mid
            states = raw[mask_mid]["_optimized_state"].unique()
            
            for st in states:
                opt_mask = mask_mid & (raw["_is_optimization"] == True) & (raw["_optimized_state"] == st)
                if not raw[opt_mask].empty:
                    opt_row = raw[opt_mask].iloc[0]
                    non_opt_mask = mask_mid & (raw["_is_optimization"] == False)
                    
                    for idx in raw[non_opt_mask].index:
                        for col in geometry_fields:
                            if col not in filled.columns:
                                continue
                            try:
                                if _is_empty_obj(filled.at[idx, col]) and not _is_empty_obj(opt_row.get(col)):
                                    filled.at[idx, col] = opt_row[col]
                            except Exception:
                                pass
        
        # Step 3: Build standardized molecule records
        standardized = []
        
        for mid in filled["molecule_id"].unique():
            mol_rows = filled[filled["molecule_id"] == mid]
            
            row = self._build_molecule_record(mid, mol_rows)
            standardized.append(row)
        
        df = pd.DataFrame(standardized)
        self.states = {row["molecule_id"]: idx for idx, row in df.iterrows()}
        
        # Final summary
        self.logger.info("=" * 50)
        self.logger.info(f"BATCH COMPLETE: {len(df)} molecules from {len(files)} files")
        self._log_batch_summary(df)
        self.logger.info("=" * 50)
        
        return df
    
    def _build_molecule_record(self, molecule_id: str, mol_rows: pd.DataFrame) -> Dict[str, Any]:
        """Build standardized molecule record from multiple calculations."""
        row = {
            "molecule_id": molecule_id,
            "filename": mol_rows.iloc[0].get("filename"),
            "smiles": None,
            "charge": None,
            "multiplicity": None,
            "cart_coords": None,
            "gibbs_Eh": None,
            "single_point_Eh": None,
            "homo_energy": None,
            "lumo_energy": None,
            "homo_lumo_gap": None,
            "orbitals": None,
            "vibrations": None,
            "ir": None,
            "raman": None,
            "nmr_shielding": None,
            "nmr_coupling": None,
            "mulliken": None,
            "bonds": None,
            "angles": None,
            "dihedrals": None,
            "tddft_states": None,
            "is_optimization": False,
            "has_tddft": False,
            "optimized_state": "S0",
            "esd_type": None,
            "calc_class": "single_point",
            "spectrum_file": None
        }
        
        # Prefer S0 OPT, fallback to S0 SP
        s0_opt = mol_rows[(mol_rows["_is_optimization"] == True) & (mol_rows["_optimized_state"] == "S0")]
        s0_sp = mol_rows[(mol_rows["_is_optimization"] == False) & (mol_rows["_optimized_state"] == "S0")]
        s0_source = s0_opt if not s0_opt.empty else s0_sp
        
        if not s0_source.empty:
            s0 = s0_source.iloc[0]
            for col in ["smiles", "charge", "multiplicity", "cart_coords", 
                        "vibrations", "ir", "raman", "nmr_shielding", "nmr_coupling",
                        "mulliken", "bonds", "angles", "dihedrals"]:
                if not _is_empty_obj(s0.get(col)):
                    row[col] = s0[col]
        
        # Get energies (prefer OPT for Gibbs, any for SP)
        for _, calc in mol_rows.iterrows():
            if calc.get("gibbs_Eh") and row["gibbs_Eh"] is None:
                row["gibbs_Eh"] = calc["gibbs_Eh"]
            if calc.get("single_point_Eh") and row["single_point_Eh"] is None:
                row["single_point_Eh"] = calc["single_point_Eh"]
            if calc.get("homo_energy") and row["homo_energy"] is None:
                row["homo_energy"] = calc["homo_energy"]
            if calc.get("lumo_energy") and row["lumo_energy"] is None:
                row["lumo_energy"] = calc["lumo_energy"]
            if calc.get("homo_lumo_gap") and row["homo_lumo_gap"] is None:
                row["homo_lumo_gap"] = calc["homo_lumo_gap"]
            if not _is_empty_obj(calc.get("orbitals")) and row["orbitals"] is None:
                row["orbitals"] = calc["orbitals"]
            if not _is_empty_obj(calc.get("tddft_states")):
                if row["tddft_states"] is None:
                    row["tddft_states"] = calc["tddft_states"]
                else:
                    try:
                        row["tddft_states"] = pd.concat([row["tddft_states"], calc["tddft_states"]], ignore_index=True)
                    except Exception:
                        pass
            if calc.get("is_optimization"):
                row["is_optimization"] = True
            if calc.get("has_tddft"):
                row["has_tddft"] = True
            if calc.get("esd_type") and row["esd_type"] is None:
                row["esd_type"] = calc["esd_type"]
            if not _is_empty_obj(calc.get("spectrum_file")) and row["spectrum_file"] is None:
                row["spectrum_file"] = calc["spectrum_file"]
        
        # Determine calc_class
        if row["is_optimization"]:
            row["calc_class"] = "optimization"
        elif row["has_tddft"]:
            row["calc_class"] = "tddft"
        elif row["esd_type"]:
            row["calc_class"] = "spectrum"
        
        return row
    
    def _log_batch_summary(self, df: pd.DataFrame):
        """Log batch summary statistics."""
        self.logger.info(f"  Total molecules: {len(df)}")
        
        stats = {
            "smiles": df["smiles"].notna().sum() if "smiles" in df else 0,
            "gibbs_Eh": df["gibbs_Eh"].notna().sum() if "gibbs_Eh" in df else 0,
            "single_point_Eh": df["single_point_Eh"].notna().sum() if "single_point_Eh" in df else 0,
            "homo_energy": df["homo_energy"].notna().sum() if "homo_energy" in df else 0,
        }
        
        for key, count in stats.items():
            if count > 0:
                pct = 100 * count / len(df)
                self.logger.info(f"  {key}: {count}/{len(df)} ({pct:.0f}%)")
        
        if "calc_class" in df:
            classes = df["calc_class"].value_counts().to_dict()
            self.logger.info(f"  Calc classes: {classes}")
    
    @staticmethod
    def _extract_molecule_id(filename: str, geometry_filename: Optional[str]) -> str:
        """Extract molecule ID from filename, removing state suffixes."""
        if geometry_filename:
            geom_file = os.path.splitext(os.path.basename(geometry_filename))[0]
            mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt)$', '', geom_file, flags=re.I)
            return mol_id.strip("_-")
        
        base = os.path.splitext(os.path.basename(filename))[0]
        mol_id = re.sub(r'(_?s0p?|_?s0|_?s1|_?t1|_?vg|_?ah|_?ahas|_?p|_?opt)$', '', base, flags=re.I)
        return mol_id.strip("_-")
