"""Comparison engine for multi-molecule comparison."""

from typing import List, Dict, Optional, Any
import pandas as pd
import numpy as np

from ..logger import get_logger


class ComparisonEngine:
    """
    Compare multiple molecules across various properties.
    
    Supports comparison of:
        - Energies
        - Orbital levels
        - Spectra
        - Geometry
    """
    
    def __init__(self, df: pd.DataFrame, id_column: str = "molecule_id"):
        """
        Initialize comparison engine.
        
        Args:
            df: DataFrame with molecule data
            id_column: Column containing molecule IDs
        """
        self.df = df
        self.id_column = id_column
        self.logger = get_logger("ComparisonEngine")
    
    def compare(
        self,
        mol_ids: List[str],
        properties: Optional[List[str]] = None
    ) -> Dict[str, pd.DataFrame]:
        """
        Compare molecules across specified properties.
        
        Args:
            mol_ids: List of molecule IDs to compare
            properties: Properties to compare (default: all)
            
        Returns:
            Dict mapping property name to comparison DataFrame
        """
        if properties is None:
            properties = ["energy", "orbitals", "spectra"]
        
        # Filter to selected molecules
        mask = self.df[self.id_column].isin(mol_ids)
        subset = self.df[mask].copy()
        
        if subset.empty:
            self.logger.warning(f"No molecules found: {mol_ids}")
            return {}
        
        self.logger.info(f"Comparing {len(subset)} molecules: {mol_ids}")
        
        results = {}
        
        if "energy" in properties:
            results["energy"] = self._compare_energy(subset)
        
        if "orbitals" in properties:
            results["orbitals"] = self._compare_orbitals(subset)
        
        if "spectra" in properties:
            results["spectra"] = self._compare_spectra(subset)
        
        return results
    
    def _compare_energy(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compare energy values."""
        cols = [self.id_column]
        
        if "gibbs_Eh" in df.columns:
            cols.append("gibbs_Eh")
        if "single_point_Eh" in df.columns:
            cols.append("single_point_Eh")
        
        result = df[cols].copy()
        
        # Calculate relative energies
        if "gibbs_Eh" in result.columns:
            ref = result["gibbs_Eh"].min()
            result["gibbs_rel_kcal"] = (result["gibbs_Eh"] - ref) * 627.509  # Eh to kcal/mol
        
        return result
    
    def _compare_orbitals(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compare orbital energies."""
        cols = [self.id_column]
        
        for col in ["homo_energy", "lumo_energy", "homo_lumo_gap"]:
            if col in df.columns:
                cols.append(col)
        
        return df[cols].copy()
    
    def _compare_spectra(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Compare spectra (returns dict, not DataFrame)."""
        result = {}
        
        for spectrum_type in ["ir", "raman"]:
            if spectrum_type in df.columns:
                spectra = {}
                for _, row in df.iterrows():
                    mol_id = row[self.id_column]
                    spectrum = row[spectrum_type]
                    if isinstance(spectrum, pd.DataFrame) and not spectrum.empty:
                        spectra[mol_id] = spectrum
                if spectra:
                    result[spectrum_type] = spectra
        
        return result
    
    def energy_ranking(self, energy_col: str = "gibbs_Eh") -> pd.DataFrame:
        """
        Rank molecules by energy.
        
        Args:
            energy_col: Energy column to use
            
        Returns:
            DataFrame with ranked molecules
        """
        if energy_col not in self.df.columns:
            self.logger.warning(f"Column {energy_col} not found")
            return pd.DataFrame()
        
        result = self.df[[self.id_column, energy_col]].copy()
        result = result.sort_values(energy_col)
        
        ref = result[energy_col].min()
        result["relative_kcal"] = (result[energy_col] - ref) * 627.509
        result["rank"] = range(1, len(result) + 1)
        
        return result
