"""Hierarchy detector for molecule naming patterns."""

import re
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import pandas as pd

from ..core.data_models import Hierarchy
from ..logger import get_logger


class HierarchyDetector:
    """
    Detect molecule hierarchy from naming patterns.
    
    Examples:
        p1x, p1a, p1b → p1 (root) with variants x, a, b
        cbzs0, cbzs1, cbzt1 → cbz (root) with states s0, s1, t1
    """
    
    def __init__(self, df: pd.DataFrame, id_column: str = "molecule_id"):
        """
        Initialize detector.
        
        Args:
            df: DataFrame with molecule data
            id_column: Column containing molecule IDs
        """
        self.df = df
        self.id_column = id_column
        self.logger = get_logger("HierarchyDetector")
    
    def detect(self) -> Hierarchy:
        """
        Detect hierarchy from molecule IDs.
        
        Returns:
            Hierarchy object with roots and variants
        """
        if self.id_column not in self.df.columns:
            self.logger.warning(f"Column {self.id_column} not found")
            return Hierarchy(roots=[], variants={})
        
        ids = self.df[self.id_column].dropna().unique().tolist()
        self.logger.info(f"Detecting hierarchy from {len(ids)} IDs")
        
        # Extract base names and variants
        roots, variants = self._parse_ids(ids)
        
        hierarchy = Hierarchy(roots=roots, variants=variants)
        self.logger.info(f"Found {len(roots)} root groups")
        
        return hierarchy
    
    def _parse_ids(self, ids: List[str]) -> Tuple[List[str], Dict[str, List[str]]]:
        """Parse IDs into roots and variants."""
        groups = defaultdict(list)
        
        for mol_id in ids:
            base, variant = self._split_id(mol_id)
            groups[base].append(mol_id)
        
        roots = sorted(groups.keys())
        variants = {k: sorted(v) for k, v in groups.items()}
        
        return roots, variants
    
    def _split_id(self, mol_id: str) -> Tuple[str, str]:
        """
        Split molecule ID into base and variant.
        
        Examples:
            p1x → (p1, x)
            p2a → (p2, a)
            cbzs0_opt → (cbzs0, opt)
        """
        # Pattern: base name + single letter/digit variant
        # e.g., p1x → p1 + x
        match = re.match(r'^(.+?)([a-z])$', mol_id, re.IGNORECASE)
        if match:
            return match.group(1), match.group(2)
        
        # Pattern: base_suffix
        if '_' in mol_id:
            parts = mol_id.rsplit('_', 1)
            return parts[0], parts[1]
        
        # No variant detected
        return mol_id, ""
    
    def get_group(self, root: str) -> pd.DataFrame:
        """Get all molecules in a root group."""
        hierarchy = self.detect()
        if root not in hierarchy.variants:
            return pd.DataFrame()
        
        ids = hierarchy.variants[root]
        return self.df[self.df[self.id_column].isin(ids)]
    
    def to_tree_string(self) -> str:
        """Generate tree representation."""
        return self.detect().to_tree()
