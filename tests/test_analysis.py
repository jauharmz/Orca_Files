"""Tests for analysis modules."""

import pytest
import pandas as pd


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
        
        print("\n=== HIERARCHY OUTPUT ===")
        print(f"Roots: {hierarchy.roots}")
        print(f"Variants: {hierarchy.variants}")
        print(f"Tree:\n{hierarchy.to_tree()}")
        
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
        
        print("\n=== PARTITION OUTPUT ===")
        print(f"By state: {partitions.by_state}")
        print(f"By calc type: {partitions.by_calc_type}")
        
        assert "S0" in partitions.by_state or "S1" in partitions.by_state


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
        
        print("\n=== PATHWAY OUTPUT ===")
        for i, pw in enumerate(pathways):
            print(f"Pathway {i}: {pw.nodes}")
            print(f"  Edges: {pw.edges}")
            print(f"  Color: {pw.color}")


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
        
        print("\n=== SPECTRAL SCALING OUTPUT ===")
        print(f"Original: {list(spectrum['freq_cm-1'])}")
        print(f"Linear scaled (0.97): {list(scaled['freq_cm-1_scaled'])}")
        
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
        
        print(f"Relative scaled (2.0): {list(scaled['freq_cm-1_scaled'])}")
        
        # ν_min + 2 * (ν - ν_min)
        # 100 + 2 * (100 - 100) = 100
        # 100 + 2 * (500 - 100) = 900
        # 100 + 2 * (1000 - 100) = 1900
        assert scaled["freq_cm-1_scaled"].iloc[0] == 100
        assert scaled["freq_cm-1_scaled"].iloc[1] == 900


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
