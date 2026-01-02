"""Partition detector for calculation types and states."""

import re
from typing import Dict, List
from collections import defaultdict
import pandas as pd

from ..core.data_models import Partitions
from ..logger import get_logger


class PartitionDetector:
    """
    Detect data partitions by state, calculation type, and ESD type.
    
    Partitions:
        - by_state: S0, S1, T1
        - by_calc_type: OPT, SP
        - by_esd_type: VG, AH, AHAS, PHOSP, FLUOR
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
        self.logger = get_logger("PartitionDetector")
    
    def detect(self) -> Partitions:
        """
        Detect all partitions.
        
        Returns:
            Partitions object
        """
        partitions = Partitions()
        
        if self.id_column not in self.df.columns:
            self.logger.warning(f"Column {self.id_column} not found")
            return partitions
        
        ids = self.df[self.id_column].dropna().unique().tolist()
        
        partitions.by_state = self._detect_states(ids)
        partitions.by_calc_type = self._detect_calc_types(ids)
        partitions.by_esd_type = self._detect_esd_types(ids)
        
        self.logger.info(
            f"Partitions: states={len(partitions.by_state)}, "
            f"calc_types={len(partitions.by_calc_type)}, "
            f"esd_types={len(partitions.by_esd_type)}"
        )
        
        return partitions
    
    def _detect_states(self, ids: List[str]) -> Dict[str, List[str]]:
        """Detect electronic state partitions."""
        states = defaultdict(list)
        
        for mol_id in ids:
            mol_lower = mol_id.lower()
            
            if 's0' in mol_lower or mol_lower.endswith('s0'):
                states['S0'].append(mol_id)
            elif 's1' in mol_lower or mol_lower.endswith('s1'):
                states['S1'].append(mol_id)
            elif 's2' in mol_lower:
                states['S2'].append(mol_id)
            elif 't1' in mol_lower or mol_lower.endswith('t1'):
                states['T1'].append(mol_id)
            elif 't2' in mol_lower:
                states['T2'].append(mol_id)
            else:
                # Default to ground state
                states['S0'].append(mol_id)
        
        return dict(states)
    
    def _detect_calc_types(self, ids: List[str]) -> Dict[str, List[str]]:
        """Detect calculation type partitions."""
        calc_types = defaultdict(list)
        
        for mol_id in ids:
            mol_lower = mol_id.lower()
            
            if '_opt' in mol_lower or 'opt' in mol_lower:
                calc_types['OPT'].append(mol_id)
            elif '_sp' in mol_lower or mol_lower.endswith('sp'):
                calc_types['SP'].append(mol_id)
            elif '_freq' in mol_lower:
                calc_types['FREQ'].append(mol_id)
            else:
                # Check DataFrame columns if available
                if 'is_optimization' in self.df.columns:
                    row = self.df[self.df[self.id_column] == mol_id]
                    if not row.empty and row.iloc[0].get('is_optimization', False):
                        calc_types['OPT'].append(mol_id)
                    else:
                        calc_types['SP'].append(mol_id)
                else:
                    calc_types['UNKNOWN'].append(mol_id)
        
        return dict(calc_types)
    
    def _detect_esd_types(self, ids: List[str]) -> Dict[str, List[str]]:
        """Detect ESD (Excited State Dynamics) type partitions."""
        esd_types = defaultdict(list)
        
        for mol_id in ids:
            mol_lower = mol_id.lower()
            
            if 'ahas' in mol_lower:
                esd_types['AHAS'].append(mol_id)
            elif 'ah' in mol_lower:
                esd_types['AH'].append(mol_id)
            elif 'vg' in mol_lower:
                esd_types['VG'].append(mol_id)
            elif 'phosp' in mol_lower:
                esd_types['PHOSP'].append(mol_id)
            elif 'fluor' in mol_lower:
                esd_types['FLUOR'].append(mol_id)
        
        return dict(esd_types)
    
    def get_partition(self, partition_type: str, value: str) -> pd.DataFrame:
        """
        Get molecules in a specific partition.
        
        Args:
            partition_type: 'state', 'calc_type', or 'esd_type'
            value: Partition value (e.g., 'S0', 'OPT', 'VG')
        """
        partitions = self.detect()
        
        if partition_type == 'state':
            ids = partitions.by_state.get(value, [])
        elif partition_type == 'calc_type':
            ids = partitions.by_calc_type.get(value, [])
        elif partition_type == 'esd_type':
            ids = partitions.by_esd_type.get(value, [])
        else:
            return pd.DataFrame()
        
        return self.df[self.df[self.id_column].isin(ids)]
