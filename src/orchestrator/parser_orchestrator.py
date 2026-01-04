"""
Parser Orchestrator - High-level API for parsing ORCA files

This module provides a clean interface for the UI layer to use,
following the architecture from README.md.

Usage:
    from src.orchestrator.parser_orchestrator import ParserOrchestrator
    
    orchestrator = ParserOrchestrator()
    df = orchestrator.parse_folder("./data")
    # df has molecule_id, optimized_state, and all parsed fields
"""

from pathlib import Path
from typing import Optional, List
import pandas as pd

from ..parser.batch import BatchParser
from ..logger import get_logger


class ParserOrchestrator:
    """
    High-level parser interface for UI/application layer.
    
    Provides simplified access to parsing functionality with
    proper state extraction and data formatting.
    """
    
    def __init__(self):
        self.logger = get_logger("ParserOrchestrator")
        self._df: Optional[pd.DataFrame] = None
    
    def parse_folder(self, folder: str, pattern: str = "**/*.out") -> pd.DataFrame:
        """
        Parse all ORCA files in folder.
        
        Uses detailed parsing mode - one row per file with state from path.
        
        Args:
            folder: Path to folder containing .out files
            pattern: Glob pattern for files (default: **/*.out)
            
        Returns:
            DataFrame with columns:
                - molecule_id: Base molecule ID (e.g., p1a)
                - optimized_state: State from path (S0, S0-SP, VG)
                - All parsed data fields
        """
        folder_path = Path(folder)
        if not folder_path.exists():
            self.logger.warning(f"Folder not found: {folder}")
            return pd.DataFrame()
        
        glob_pattern = str(folder_path / pattern)
        self.logger.info(f"Parsing folder: {folder}")
        
        batch = BatchParser(glob_pattern)
        self._df = batch.parse_all_detailed(verbose=True)
        
        self.logger.info(f"Parsed {len(self._df)} files")
        return self._df
    
    def parse_files(self, filepaths: List[str]) -> pd.DataFrame:
        """
        Parse a list of specific file paths.
        
        Args:
            filepaths: List of paths to .out files
            
        Returns:
            DataFrame with parsed data
        """
        if not filepaths:
            return pd.DataFrame()
        
        # Create temporary batch with file list
        batch = BatchParser()
        batch.files = filepaths
        self._df = batch.parse_all_detailed(verbose=False)
        
        return self._df
    
    def get_dataframe(self) -> Optional[pd.DataFrame]:
        """Get the currently loaded DataFrame."""
        return self._df
    
    def get_molecules(self) -> List[str]:
        """Get list of unique molecule IDs."""
        if self._df is None:
            return []
        return sorted(self._df["molecule_id"].dropna().unique().tolist())
    
    def get_states(self) -> List[str]:
        """Get list of unique states."""
        if self._df is None or "optimized_state" not in self._df.columns:
            return []
        return sorted(self._df["optimized_state"].dropna().unique().tolist())
    
    def get_molecule_state_pairs(self) -> List[tuple]:
        """Get list of (molecule_id, state) pairs."""
        if self._df is None:
            return []
        
        pairs = []
        for _, row in self._df.iterrows():
            mol = row.get("molecule_id")
            state = row.get("optimized_state", "unknown")
            if mol:
                pairs.append((mol, state))
        return sorted(set(pairs))
    
    def filter_by_molecule(self, molecule_ids: List[str]) -> pd.DataFrame:
        """Filter DataFrame by molecule IDs."""
        if self._df is None:
            return pd.DataFrame()
        return self._df[self._df["molecule_id"].isin(molecule_ids)]
    
    def filter_by_state(self, states: List[str]) -> pd.DataFrame:
        """Filter DataFrame by states."""
        if self._df is None or "optimized_state" not in self._df.columns:
            return pd.DataFrame()
        return self._df[self._df["optimized_state"].isin(states)]
    
    def filter_by_molecule_and_state(
        self, 
        molecule_ids: Optional[List[str]] = None,
        states: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """Filter DataFrame by molecule IDs and/or states."""
        if self._df is None:
            return pd.DataFrame()
        
        df = self._df.copy()
        
        if molecule_ids:
            df = df[df["molecule_id"].isin(molecule_ids)]
        
        if states and "optimized_state" in df.columns:
            df = df[df["optimized_state"].isin(states)]
        
        return df
    
    def get_summary(self) -> dict:
        """Get summary statistics of loaded data."""
        if self._df is None:
            return {"total": 0}
        
        summary = {
            "total": len(self._df),
            "molecules": self._df["molecule_id"].nunique() if "molecule_id" in self._df.columns else 0,
        }
        
        if "optimized_state" in self._df.columns:
            summary["states"] = dict(self._df["optimized_state"].value_counts())
        
        # Count data availability
        for col in ["cart_coords", "orbitals", "ir", "raman", "tddft_states"]:
            if col in self._df.columns:
                count = self._df[col].apply(lambda x: x is not None and hasattr(x, '__len__') and len(x) > 0).sum()
                summary[f"has_{col}"] = int(count)
        
        return summary
