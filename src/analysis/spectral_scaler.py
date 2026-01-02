"""Spectral scaler for frequency scaling."""

from typing import Optional
import numpy as np
import pandas as pd

from ..logger import get_logger


class SpectralScaler:
    """
    Scale spectral frequencies using linear or relative methods.
    
    Linear Scale: ν_s = s × ν
        - Frequencies scaled directly from zero
        - Used for DFT frequency correction
    
    Relative Scale: ν_s = ν_min + s × (ν - ν_min)
        - Frequencies scaled relative to minimum
        - Used for visual comparison
    """
    
    def __init__(self, spectrum: pd.DataFrame, freq_column: str = "freq_cm-1"):
        """
        Initialize scaler.
        
        Args:
            spectrum: DataFrame with frequency data
            freq_column: Column containing frequencies
        """
        self.spectrum = spectrum.copy()
        self.freq_column = freq_column
        self.logger = get_logger("SpectralScaler")
    
    def linear_scale(self, factor: float) -> pd.DataFrame:
        """
        Apply linear scaling: ν_s = s × ν
        
        Args:
            factor: Scale factor (e.g., 0.97 for DFT correction)
            
        Returns:
            Scaled spectrum DataFrame
        """
        result = self.spectrum.copy()
        
        if self.freq_column not in result.columns:
            self.logger.warning(f"Column {self.freq_column} not found")
            return result
        
        original = result[self.freq_column].values
        scaled = factor * original
        
        result[f"{self.freq_column}_scaled"] = scaled
        result["scale_type"] = "linear"
        result["scale_factor"] = factor
        
        self.logger.info(f"Linear scale factor={factor}: {original.min():.2f}-{original.max():.2f} → {scaled.min():.2f}-{scaled.max():.2f} cm⁻¹")
        
        return result
    
    def relative_scale(self, factor: float) -> pd.DataFrame:
        """
        Apply relative scaling: ν_s = ν_min + s × (ν - ν_min)
        
        Args:
            factor: Scale factor
            
        Returns:
            Scaled spectrum DataFrame
        """
        result = self.spectrum.copy()
        
        if self.freq_column not in result.columns:
            self.logger.warning(f"Column {self.freq_column} not found")
            return result
        
        original = result[self.freq_column].values
        freq_min = original.min()
        
        # ν_s = ν_min + s × (ν - ν_min)
        scaled = freq_min + factor * (original - freq_min)
        
        result[f"{self.freq_column}_scaled"] = scaled
        result["scale_type"] = "relative"
        result["scale_factor"] = factor
        result["anchor_freq"] = freq_min
        
        self.logger.info(f"Relative scale factor={factor} anchor={freq_min:.2f}: {original.max():.2f} → {scaled.max():.2f} cm⁻¹")
        
        return result
    
    def normalize_intensity(self, intensity_column: str = "intensity") -> pd.DataFrame:
        """
        Normalize intensity to 0-1 range.
        
        Args:
            intensity_column: Column to normalize
            
        Returns:
            Spectrum with normalized intensity
        """
        result = self.spectrum.copy()
        
        if intensity_column not in result.columns:
            # Try common column names
            for col in ["intensity_km/mol", "activity", "eps"]:
                if col in result.columns:
                    intensity_column = col
                    break
            else:
                self.logger.warning("No intensity column found")
                return result
        
        values = result[intensity_column].values
        max_val = values.max()
        
        if max_val > 0:
            result[f"{intensity_column}_norm"] = values / max_val
        
        return result
    
    def compare_scales(self, factor: float) -> pd.DataFrame:
        """
        Compare linear and relative scaling side by side.
        
        Args:
            factor: Scale factor to apply
            
        Returns:
            DataFrame with both scaling results
        """
        linear = self.linear_scale(factor)
        relative = self.relative_scale(factor)
        
        result = self.spectrum.copy()
        result["linear_scaled"] = linear[f"{self.freq_column}_scaled"]
        result["relative_scaled"] = relative[f"{self.freq_column}_scaled"]
        
        return result
