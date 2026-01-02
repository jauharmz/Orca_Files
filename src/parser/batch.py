"""Batch parser for multiple ORCA files."""

from typing import List, Dict, Tuple, Optional
from pathlib import Path
import pandas as pd

from ..core.data_models import ParseResult
from ..logger import get_logger
from .factory import ParserFactory


class BatchParser:
    """Parse multiple ORCA output files."""
    
    def __init__(self):
        self.logger = get_logger("BatchParser")
        self.factory = ParserFactory()
    
    def parse_files(self, filepaths: List[str], verbose: bool = True) -> pd.DataFrame:
        """
        Parse multiple files and return DataFrame.
        
        Args:
            filepaths: List of file paths
            verbose: Whether to log progress
            
        Returns:
            DataFrame with parsed data
        """
        results = []
        
        for i, filepath in enumerate(filepaths):
            if verbose:
                self.logger.info(f"Parsing [{i+1}/{len(filepaths)}]: {Path(filepath).name}")
            
            try:
                result = self.factory.parse(filepath)
                row = self._result_to_row(result)
                results.append(row)
            except Exception as e:
                self.logger.error(f"Failed to parse {filepath}: {e}")
                continue
        
        if not results:
            return pd.DataFrame()
        
        df = pd.DataFrame(results)
        
        # Extract molecule_id from geometry filename
        if "geometry_filename" in df.columns:
            df["molecule_id"] = df["geometry_filename"].apply(
                lambda x: self._extract_molecule_id(x) if x else None
            )
        
        self.logger.info(f"Parsed {len(df)} files")
        return df
    
    def parse_folder(self, folder: str, pattern: str = "*.out") -> pd.DataFrame:
        """
        Parse all .out files in a folder.
        
        Args:
            folder: Folder path
            pattern: Glob pattern (default: *.out)
            
        Returns:
            DataFrame with parsed data
        """
        path = Path(folder)
        files = list(path.rglob(pattern))
        self.logger.info(f"Found {len(files)} files in {folder}")
        return self.parse_files([str(f) for f in files])
    
    def _result_to_row(self, result: ParseResult) -> Dict:
        """Convert ParseResult to DataFrame row."""
        return {
            "filename": result.filename,
            "geometry_filename": result.geometry.filename,
            "smiles": result.geometry.smiles,
            "charge": result.geometry.charge,
            "multiplicity": result.geometry.multiplicity,
            "cart_coords": result.geometry.cart_coords,
            "gibbs_Eh": result.energy.gibbs_Eh,
            "single_point_Eh": result.energy.single_point_Eh,
            "orbitals": result.orbitals.orbitals,
            "homo_energy": result.orbitals.homo_energy,
            "lumo_energy": result.orbitals.lumo_energy,
            "homo_lumo_gap": result.orbitals.homo_lumo_gap,
            "ir": result.spectra.ir,
            "raman": result.spectra.raman,
            "vibrations": result.spectra.vibrations,
            "tddft_states": result.tddft.states,
            "electric_dipole": result.tddft.electric_dipole,
            "is_optimization": result.is_optimization,
            "has_tddft": result.has_tddft,
            "optimized_state": result.optimized_state,
        }
    
    def _extract_molecule_id(self, filename: str) -> str:
        """Extract molecule ID from filename."""
        # Remove common suffixes to get base molecule ID
        import re
        
        # Remove extension
        name = Path(filename).stem
        
        # Remove state suffixes like _s0, _t1, _sp, _opt
        name = re.sub(r'[_-]?(s[0-9]+|t[0-9]+|sp|opt|freq|tddft)$', '', name, flags=re.IGNORECASE)
        
        return name
