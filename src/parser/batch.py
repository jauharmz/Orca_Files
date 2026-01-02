"""
Batch Parser - Parse multiple ORCA files

Produces a DataFrame where each row is a parsed calculation.
"""

import os
import re
import glob
from pathlib import Path
from typing import List, Optional, Dict, Any
import pandas as pd

from .factory import ParserFactory
from ..logger import get_logger


class BatchParser:
    """Batch parser for multiple ORCA output files."""
    
    def __init__(self, *patterns: str):
        """
        Initialize with glob patterns.
        
        Args:
            *patterns: Glob patterns like "*.out" or "data/**/*.out"
        """
        self.logger = get_logger("BatchParser")
        self.patterns = patterns
        self.files: List[str] = []
        
        # Expand patterns to file list
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
        """Parse all files and return DataFrame."""
        files = self.files[:limit_files] if limit_files else self.files
        
        if not files:
            self.logger.warning("No files to parse")
            return pd.DataFrame()
        
        self.logger.info(f"Starting batch parse of {len(files)} files")
        
        factory = ParserFactory()
        rows = []
        success_count = 0
        error_count = 0
        
        for i, filepath in enumerate(files):
            if verbose and (i + 1) % 5 == 0:
                self.logger.info(f"Progress: {i+1}/{len(files)} ({100*(i+1)/len(files):.0f}%)")
            
            try:
                result = factory.parse(filepath)
                row = result.to_dict()
                
                # Extract molecule_id
                row["molecule_id"] = self._extract_molecule_id(filepath, result.geometry.filename)
                
                rows.append(row)
                success_count += 1
            except Exception as e:
                self.logger.warning(f"Failed to parse {filepath}: {e}")
                error_count += 1
                continue
        
        if not rows:
            self.logger.warning("No files parsed successfully")
            return pd.DataFrame()
        
        df = pd.DataFrame(rows)
        
        # Final summary
        self.logger.info("=" * 50)
        self.logger.info(f"BATCH COMPLETE: {success_count}/{len(files)} parsed, {error_count} errors")
        self._log_batch_summary(df)
        self.logger.info("=" * 50)
        
        return df
    
    def _log_batch_summary(self, df: pd.DataFrame):
        """Log batch summary statistics."""
        self.logger.info(f"  Total molecules: {len(df)}")
        
        # Count data availability
        stats = {
            "smiles": df["smiles"].notna().sum() if "smiles" in df else 0,
            "gibbs_Eh": df["gibbs_Eh"].notna().sum() if "gibbs_Eh" in df else 0,
            "single_point_Eh": df["single_point_Eh"].notna().sum() if "single_point_Eh" in df else 0,
            "homo_energy": df["homo_energy"].notna().sum() if "homo_energy" in df else 0,
            "lumo_energy": df["lumo_energy"].notna().sum() if "lumo_energy" in df else 0,
        }
        
        for key, count in stats.items():
            if count > 0:
                pct = 100 * count / len(df)
                self.logger.info(f"  {key}: {count}/{len(df)} ({pct:.0f}%)")
        
        # Calc class distribution
        if "calc_class" in df:
            classes = df["calc_class"].value_counts().to_dict()
            self.logger.info(f"  Calc classes: {classes}")
    
    @staticmethod
    def _extract_molecule_id(filename: str, geometry_filename: Optional[str]) -> str:
        """Extract molecule ID from filename."""
        if geometry_filename:
            return geometry_filename
        
        base = os.path.splitext(os.path.basename(filename))[0]
        cleaned = re.sub(r'[-_](opt|sp|freq|td|esd|vg|ahas|ah)$', '', base, flags=re.I)
        return cleaned
