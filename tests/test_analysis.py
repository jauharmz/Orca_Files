"""
Analysis Unit Tests

Run: pytest tests/test_analysis.py -v -s
"""

import sys
from pathlib import Path
import logging
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Setup logging
from src.logger import setup_logging, get_logger
setup_logging(level=logging.DEBUG)
logger = get_logger("test_analysis")


class TestHierarchyDetector:
    """Test hierarchy detection."""
    
    def test_detect_variants(self):
        """Test detecting molecule variants."""
        from src.analysis.hierarchy_detector import HierarchyDetector
        
        df = pd.DataFrame({
            "molecule_id": ["p1x", "p1a", "p1b", "p2x", "p2a"]
        })
        
        detector = HierarchyDetector(df)
        hierarchy = detector.detect()
        
        logger.info("=== HIERARCHY OUTPUT ===")
        logger.info(f"Roots: {hierarchy.roots}")
        logger.info(f"Variants: {hierarchy.variants}")
        logger.info(f"Tree:\n{hierarchy.to_tree()}")
        
        assert len(hierarchy.roots) > 0


class TestPartitionDetector:
    """Test partition detection."""
    
    def test_detect_states(self):
        """Test detecting electronic states."""
        from src.analysis.partition_detector import PartitionDetector
        
        df = pd.DataFrame({
            "molecule_id": ["mol_s0_opt", "mol_s1_opt", "mol_t1_opt", "mol_sp"]
        })
        
        detector = PartitionDetector(df)
        partitions = detector.detect()
        
        logger.info("=== PARTITION OUTPUT ===")
        logger.info(f"By state: {partitions.by_state}")
        logger.info(f"By calc type: {partitions.by_calc_type}")
        
        assert len(partitions.by_state) > 0 or len(partitions.by_calc_type) > 0


class TestPathwayDetector:
    """Test pathway detection."""
    
    def test_detect_edges(self):
        """Test edge detection."""
        from src.analysis.pathway_detector import PathwayDetector
        
        df = pd.DataFrame({
            "molecule_id": ["p1x", "p2x", "p3x", "p1a", "p2a"]
        })
        
        detector = PathwayDetector(df)
        pathways = detector.detect()
        
        logger.info("=== PATHWAY OUTPUT ===")
        for i, pw in enumerate(pathways):
            logger.info(f"Pathway {i}: nodes={pw.nodes}, edges={pw.edges}, color={pw.color}")


class TestSpectralScaler:
    """Test spectral scaling."""
    
    def test_linear_scale(self):
        """Test linear scaling."""
        from src.analysis.spectral_scaler import SpectralScaler
        
        spectrum = pd.DataFrame({
            "freq_cm-1": [1000, 1500, 2000],
            "intensity": [0.5, 1.0, 0.3]
        })
        
        scaler = SpectralScaler(spectrum)
        scaled = scaler.linear_scale(0.97)
        
        logger.info("=== SPECTRAL SCALING OUTPUT ===")
        logger.info(f"Original: {list(spectrum['freq_cm-1'])}")
        logger.info(f"Linear scaled (0.97): {list(scaled['freq_cm-1_scaled'])}")
        
        # 1000 * 0.97 = 970
        assert scaled["freq_cm-1_scaled"].iloc[0] == 970
    
    def test_relative_scale(self):
        """Test relative scaling."""
        from src.analysis.spectral_scaler import SpectralScaler
        
        spectrum = pd.DataFrame({
            "freq_cm-1": [100, 500, 1000],
            "intensity": [0.5, 1.0, 0.3]
        })
        
        scaler = SpectralScaler(spectrum)
        scaled = scaler.relative_scale(2.0)
        
        logger.info(f"Relative scaled (2.0): {list(scaled['freq_cm-1_scaled'])}")
        
        # ν_min + 2 * (ν - ν_min)
        # 100 + 2 * (100 - 100) = 100
        # 100 + 2 * (500 - 100) = 900
        assert scaled["freq_cm-1_scaled"].iloc[0] == 100
        assert scaled["freq_cm-1_scaled"].iloc[1] == 900


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v", "-s"])
