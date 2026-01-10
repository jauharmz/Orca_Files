"""Data exporter for parsed ORCA data."""

import json
from pathlib import Path
from typing import Any, Dict, Optional
import pandas as pd

from ..logger import get_logger


class DataExporter:
    """
    Export parsed ORCA data to various formats.
    
    Supported formats:
        - JSON
        - CSV
        - Parquet
        - Pickle
    """
    
    def __init__(self, df: pd.DataFrame, metadata: Optional[Dict] = None):
        """
        Initialize exporter.
        
        Args:
            df: DataFrame to export
            metadata: Optional metadata dict
        """
        self.df = df
        self.metadata = metadata or {}
        self.logger = get_logger("DataExporter")
    
    def to_json(self, filepath: str, include_nested: bool = False) -> str:
        """
        Export to JSON.
        
        Args:
            filepath: Output path
            include_nested: Whether to include nested DataFrames
            
        Returns:
            Output path
        """
        path = Path(filepath)
        
        if include_nested:
            # Convert all data
            data = self._df_to_dict(self.df)
        else:
            # Only scalar columns
            scalar_cols = [c for c in self.df.columns 
                          if not self.df[c].apply(lambda x: isinstance(x, pd.DataFrame)).any()]
            data = self.df[scalar_cols].to_dict(orient='records')
        
        # Add metadata
        output = {
            "metadata": self.metadata,
            "data": data
        }
        
        with open(path, 'w') as f:
            json.dump(output, f, indent=2, default=str)
        
        self.logger.info(f"Exported JSON: {path}")
        return str(path)
    
    def to_csv(self, filepath: str) -> str:
        """
        Export to CSV (scalar columns only).
        
        Args:
            filepath: Output path
            
        Returns:
            Output path
        """
        path = Path(filepath)
        
        # Only scalar columns
        scalar_cols = [c for c in self.df.columns 
                      if not self.df[c].apply(lambda x: isinstance(x, (pd.DataFrame, dict, list))).any()]
        
        self.df[scalar_cols].to_csv(path, index=False)
        self.logger.info(f"Exported CSV: {path}")
        return str(path)
    
    def to_parquet(self, filepath: str) -> str:
        """
        Export to Parquet.
        
        Args:
            filepath: Output path
            
        Returns:
            Output path
        """
        path = Path(filepath)
        
        # Only scalar columns for parquet
        scalar_cols = [c for c in self.df.columns 
                      if not self.df[c].apply(lambda x: isinstance(x, (pd.DataFrame, dict, list))).any()]
        
        self.df[scalar_cols].to_parquet(path, index=False)
        self.logger.info(f"Exported Parquet: {path}")
        return str(path)
    
    def to_pickle(self, filepath: str) -> str:
        """
        Export to Pickle (preserves all data types).
        
        Args:
            filepath: Output path
            
        Returns:
            Output path
        """
        path = Path(filepath)
        self.df.to_pickle(path)
        self.logger.info(f"Exported Pickle: {path}")
        return str(path)
    
    def export_bundle(self, folder: str, formats: Optional[list] = None) -> Dict[str, str]:
        """
        Export to multiple formats in a folder.
        
        Args:
            folder: Output folder
            formats: List of formats (default: all)
            
        Returns:
            Dict of format -> filepath
        """
        if formats is None:
            formats = ["json", "csv", "parquet", "pkl"]
        
        path = Path(folder)
        path.mkdir(parents=True, exist_ok=True)
        
        outputs = {}
        
        if "json" in formats:
            outputs["json"] = self.to_json(path / "data.json")
        if "csv" in formats:
            outputs["csv"] = self.to_csv(path / "data.csv")
        if "parquet" in formats:
            outputs["parquet"] = self.to_parquet(path / "data.parquet")
        if "pkl" in formats:
            outputs["pkl"] = self.to_pickle(path / "data.pkl")
        
        # Export metadata separately
        if self.metadata:
            meta_path = path / "metadata.json"
            with open(meta_path, 'w') as f:
                json.dump(self.metadata, f, indent=2, default=str)
            outputs["metadata"] = str(meta_path)
        
        self.logger.info(f"Exported bundle to: {folder}")
        return outputs
    
    def _df_to_dict(self, df: pd.DataFrame) -> list:
        """Convert DataFrame with nested DataFrames to dict."""
        records = []
        for _, row in df.iterrows():
            record = {}
            for col in df.columns:
                val = row[col]
                if isinstance(val, pd.DataFrame):
                    record[col] = val.to_dict(orient='records')
                elif isinstance(val, dict):
                    # Handle nested dicts with DataFrames
                    record[col] = {}
                    for k, v in val.items():
                        if isinstance(v, pd.DataFrame):
                            record[col][k] = v.to_dict(orient='records')
                        else:
                            record[col][k] = v
                else:
                    record[col] = val
            records.append(record)
        return records
