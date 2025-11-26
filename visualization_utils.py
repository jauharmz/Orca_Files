"""
Visualization Utilities for ORCA Output Analysis

Provides reusable functions for creating publication-quality spectroscopy
visualizations, including multi-dataset stacking and advanced plotting features.

Based on implementations from 0cbz.ipynb notebook.
"""

import numpy as np
import logging
from typing import List, Tuple, Optional, Dict, Any

# Configure logging
logger = logging.getLogger(__name__)

# Physical constants for Raman intensity calculations
H_PLANCK = 6.62607015e-34      # Planck constant (J·s)
C_LIGHT = 2.99792458e10        # Speed of light (cm/s)
KB_BOLTZMANN = 1.380649e-23    # Boltzmann constant (J/K)
TEMPERATURE_K = 298.15         # Standard temperature (K)
LASER_WL_NM = 532             # Default laser wavelength (nm)


def gaussian_broadening(x_grid: np.ndarray, peaks: List[Tuple[float, float]],
                        fwhm: float) -> np.ndarray:
    """
    Apply Gaussian broadening to discrete peaks.

    Args:
        x_grid: X-axis grid for continuous spectrum
        peaks: List of (position, intensity) tuples
        fwhm: Full width at half maximum in same units as x_grid

    Returns:
        Continuous spectrum with Gaussian broadening
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    spectrum = np.zeros_like(x_grid)

    for position, intensity in peaks:
        spectrum += intensity * np.exp(-(x_grid - position)**2 / (2 * sigma**2))

    return spectrum


def calculate_stack_offsets(spectra: List[np.ndarray], y_space: float = 0.1) -> List[float]:
    """
    Calculate vertical offsets for stacked spectra with minimal spacing.

    Stacking logic: next_offset = prev_top - current_bottom + y_space

    Args:
        spectra: List of spectrum arrays to stack
        y_space: Minimum vertical spacing between spectra

    Returns:
        List of y-offsets for each spectrum
    """
    offsets = [0.0]  # First spectrum at baseline
    current_top = np.max(spectra[0]) if len(spectra[0]) > 0 else 0.0

    for spectrum in spectra[1:]:
        spectrum_bottom = np.min(spectrum) if len(spectrum) > 0 else 0.0
        next_offset = current_top - spectrum_bottom + y_space
        offsets.append(next_offset)
        current_top = next_offset + (np.max(spectrum) if len(spectrum) > 0 else 0.0)

    logger.debug(f"Calculated {len(offsets)} stack offsets: {offsets[:5]}...")
    return offsets


def raman_activity_to_intensity(frequencies: np.ndarray, activities: np.ndarray,
                                laser_wavelength_nm: float = LASER_WL_NM,
                                temperature_k: float = TEMPERATURE_K,
                                low_freq_cutoff: float = 0.0) -> np.ndarray:
    """
    Convert Raman activity to temperature-dependent intensity.

    Physical intensity formula:
    I ∝ (ν₀ - νᵢ)⁴ × S / [νᵢ × (1 - exp(-hcνᵢ/kBT))]

    Args:
        frequencies: Raman shift values (cm⁻¹)
        activities: Raman activities from DFT
        laser_wavelength_nm: Excitation laser wavelength (nm)
        temperature_k: Temperature (K)
        low_freq_cutoff: Suppress modes below this frequency (cm⁻¹)

    Returns:
        Temperature-corrected Raman intensities
    """
    laser_wavenumber = 1e7 / laser_wavelength_nm  # Convert nm → cm⁻¹

    # Apply low-frequency cutoff
    mask = frequencies >= low_freq_cutoff
    filtered_freqs = np.where(mask, frequencies, 1.0)  # Avoid division by zero
    filtered_acts = np.where(mask, activities, 0.0)

    # Calculate physical intensity
    intensity = (
        (laser_wavenumber - filtered_freqs)**4
        * filtered_acts
        / (
            filtered_freqs
            * (1 - np.exp(-H_PLANCK * C_LIGHT * filtered_freqs / (KB_BOLTZMANN * temperature_k)))
        )
    )

    # Set suppressed modes to zero
    intensity = np.where(mask, intensity, 0.0)

    logger.debug(f"Converted {len(frequencies)} Raman activities to intensities")
    return intensity


def calculate_absorbance_transmittance(ir_intensities: np.ndarray,
                                      concentration_mol_l: float = 0.01,
                                      path_length_cm: float = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate IR absorbance and transmittance from DFT intensities.

    Beer-Lambert law: A = ε × c × l
    Transmittance: T = 10^(-A) × 100%

    Args:
        ir_intensities: IR intensities (km/mol)
        concentration_mol_l: Sample concentration (mol/L)
        path_length_cm: Path length (cm)

    Returns:
        (absorbance, transmittance_percent) tuple
    """
    # Normalize intensities to represent molar absorptivity (ε)
    epsilon = ir_intensities / np.max(ir_intensities) if np.max(ir_intensities) > 0 else ir_intensities

    # Calculate absorbance
    absorbance = epsilon * concentration_mol_l * path_length_cm

    # Calculate transmittance
    transmittance = 10**(-absorbance) * 100

    return absorbance, transmittance


def apply_frequency_scaling(frequencies: np.ndarray,
                           scale_factor: float = 1.0,
                           shift_cm: float = 0.0) -> np.ndarray:
    """
    Apply scaling and shift corrections to vibrational frequencies.

    DFT calculations often overestimate frequencies. Common scale factors:
    - B3LYP/6-311++G(d,p): 0.96-0.97
    - ωB97X-D/6-311++G(d,p): 0.95-0.96

    Args:
        frequencies: Original frequencies (cm⁻¹)
        scale_factor: Multiplicative scaling factor
        shift_cm: Additive shift (cm⁻¹)

    Returns:
        Scaled and shifted frequencies
    """
    scaled = frequencies * scale_factor + shift_cm
    logger.debug(f"Scaled frequencies: factor={scale_factor}, shift={shift_cm}")
    return scaled


def get_color_palette(n_colors: int, palette: str = 'default') -> List[str]:
    """
    Generate color palette for multi-dataset plots.

    Args:
        n_colors: Number of colors needed
        palette: Color scheme ('default', 'vibrant', 'pastel')

    Returns:
        List of color codes
    """
    palettes = {
        'default': ['black', 'red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray'],
        'vibrant': ['#E63946', '#F1FAEE', '#A8DADC', '#457B9D', '#1D3557', '#F77F00', '#06FFA5'],
        'pastel': ['#FFB4A2', '#E5989B', '#B5838D', '#6D6875', '#FEC89A', '#F0EFEB']
    }

    colors = palettes.get(palette, palettes['default'])

    # Repeat colors if needed
    if n_colors > len(colors):
        colors = colors * ((n_colors // len(colors)) + 1)

    return colors[:n_colors]


def get_ir_region_boundaries() -> List[int]:
    """
    Get standard IR spectroscopy region boundaries.

    Returns:
        List of wavenumbers (cm⁻¹) marking region boundaries
    """
    return [600, 1000, 1400, 1600, 2000, 2400, 3000]


def get_ir_region_labels() -> List[Tuple[float, str]]:
    """
    Get IR functional group region labels.

    Returns:
        List of (wavenumber, label) tuples
    """
    return [
        (3300, 'O-H/N-H'),
        (2900, 'C-H'),
        (1700, 'C=O'),
        (1500, 'C=C'),
        (1200, 'C-O'),
        (800, 'Fingerprint')
    ]


def get_raman_region_boundaries() -> List[int]:
    """
    Get standard Raman spectroscopy region boundaries.

    Returns:
        List of wavenumbers (cm⁻¹) marking region boundaries
    """
    return [300, 550, 900, 1400]


def get_raman_region_labels() -> List[Tuple[float, str]]:
    """
    Get Raman region labels.

    Returns:
        List of (wavenumber, label) tuples
    """
    return [
        (150, 'Low freq'),
        (425, 'Fingerprint'),
        (725, 'C-C/Ring'),
        (1150, 'C-H bend'),
        (1700, 'C=O/C=C'),
        (2800, 'C-H stretch'),
        (3200, 'O-H/N-H')
    ]


def normalize_spectrum(spectrum: np.ndarray, method: str = 'max') -> np.ndarray:
    """
    Normalize spectrum.

    Args:
        spectrum: Spectrum array
        method: Normalization method ('max', 'minmax', 'area')

    Returns:
        Normalized spectrum
    """
    if method == 'max':
        max_val = np.max(spectrum)
        return spectrum / max_val if max_val > 0 else spectrum

    elif method == 'minmax':
        min_val, max_val = np.min(spectrum), np.max(spectrum)
        if max_val > min_val:
            return (spectrum - min_val) / (max_val - min_val)
        return spectrum

    elif method == 'area':
        area = np.trapz(spectrum)
        return spectrum / area if area > 0 else spectrum

    else:
        logger.warning(f"Unknown normalization method '{method}', using 'max'")
        return normalize_spectrum(spectrum, 'max')


def interpolate_spectrum(
    x_data: np.ndarray,
    y_data: np.ndarray,
    x_new: np.ndarray,
    method: str = 'linear',
    fill_value: float = 0.0
) -> np.ndarray:
    """
    Interpolate spectrum onto a new x-axis grid.

    Useful for comparing experimental and calculated spectra where
    data points don't align.

    Args:
        x_data: Original x-axis values
        y_data: Original y-axis values (spectrum)
        x_new: New x-axis grid for interpolation
        method: Interpolation method ('linear', 'cubic', 'nearest')
        fill_value: Value to use for points outside data range

    Returns:
        Interpolated spectrum on x_new grid
    """
    from scipy import interpolate

    if method == 'linear':
        f = interpolate.interp1d(x_data, y_data, kind='linear',
                                bounds_error=False, fill_value=fill_value)
    elif method == 'cubic':
        if len(x_data) < 4:
            logger.warning("Not enough points for cubic interpolation, using linear")
            f = interpolate.interp1d(x_data, y_data, kind='linear',
                                    bounds_error=False, fill_value=fill_value)
        else:
            f = interpolate.interp1d(x_data, y_data, kind='cubic',
                                    bounds_error=False, fill_value=fill_value)
    elif method == 'nearest':
        f = interpolate.interp1d(x_data, y_data, kind='nearest',
                                bounds_error=False, fill_value=fill_value)
    else:
        logger.warning(f"Unknown interpolation method '{method}', using linear")
        f = interpolate.interp1d(x_data, y_data, kind='linear',
                                bounds_error=False, fill_value=fill_value)

    y_new = f(x_new)
    logger.debug(f"Interpolated spectrum from {len(x_data)} to {len(x_new)} points")
    return y_new


def align_spectra_to_common_grid(
    spectra_list: List[Tuple[np.ndarray, np.ndarray]],
    x_min: Optional[float] = None,
    x_max: Optional[float] = None,
    num_points: int = 1000,
    method: str = 'linear'
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Align multiple spectra onto a common x-axis grid.

    Args:
        spectra_list: List of (x_data, y_data) tuples for each spectrum
        x_min: Minimum x value (None = use global minimum)
        x_max: Maximum x value (None = use global maximum)
        num_points: Number of points in common grid
        method: Interpolation method

    Returns:
        (common_x_grid, list_of_interpolated_spectra)
    """
    # Find global x range if not specified
    if x_min is None:
        x_min = min(np.min(x) for x, y in spectra_list)
    if x_max is None:
        x_max = max(np.max(x) for x, y in spectra_list)

    # Create common grid
    x_common = np.linspace(x_min, x_max, num_points)

    # Interpolate each spectrum onto common grid
    interpolated_spectra = []
    for i, (x_data, y_data) in enumerate(spectra_list):
        y_interp = interpolate_spectrum(x_data, y_data, x_common, method=method)
        interpolated_spectra.append(y_interp)
        logger.debug(f"Aligned spectrum {i+1}/{len(spectra_list)}")

    logger.info(f"Aligned {len(spectra_list)} spectra to common grid "
                f"({x_min:.1f} to {x_max:.1f}, {num_points} points)")

    return x_common, interpolated_spectra


def find_optimal_label_positions(
    spectra: List[np.ndarray],
    y_offsets: List[float],
    x_range: Tuple[float, float],
    x_grid: np.ndarray,
    method: str = 'peak'
) -> List[Tuple[float, float]]:
    """
    Find optimal (x, y) positions for dataset labels on stacked spectra.

    Args:
        spectra: List of spectrum arrays
        y_offsets: Vertical offsets for each spectrum
        x_range: (x_min, x_max) range to consider
        x_grid: X-axis grid corresponding to spectra
        method: 'peak' (at maximum), 'flat' (at flattest region), 'right' (right side)

    Returns:
        List of (x, y) positions for each label
    """
    positions = []

    for i, (spectrum, offset) in enumerate(zip(spectra, y_offsets)):
        offset_spectrum = spectrum + offset

        # Find indices within x_range
        mask = (x_grid >= x_range[0]) & (x_grid <= x_range[1])
        x_subset = x_grid[mask]
        y_subset = offset_spectrum[mask]

        if len(y_subset) == 0:
            # Fallback to right edge
            positions.append((x_range[1], offset))
            continue

        if method == 'peak':
            # Place at maximum intensity
            max_idx = np.argmax(y_subset)
            x_pos = x_subset[max_idx]
            y_pos = y_subset[max_idx]

        elif method == 'flat':
            # Place at flattest region (minimum gradient)
            gradient = np.abs(np.gradient(y_subset))
            flat_idx = np.argmin(gradient[10:-10]) + 10  # Avoid edges
            x_pos = x_subset[flat_idx]
            y_pos = y_subset[flat_idx]

        elif method == 'right':
            # Place at right edge
            x_pos = x_range[1]
            y_pos = offset

        else:
            # Default to right edge
            x_pos = x_range[1]
            y_pos = offset

        positions.append((x_pos, y_pos))

    logger.debug(f"Found optimal label positions using method '{method}'")
    return positions


def calculate_spectrum_similarity(
    spectrum1: np.ndarray,
    spectrum2: np.ndarray,
    method: str = 'correlation'
) -> float:
    """
    Calculate similarity between two spectra.

    Useful for quantifying agreement between experimental and calculated spectra.

    Args:
        spectrum1: First spectrum
        spectrum2: Second spectrum (must be same length as spectrum1)
        method: Similarity metric ('correlation', 'rmse', 'mae', 'r2')

    Returns:
        Similarity score
    """
    if len(spectrum1) != len(spectrum2):
        raise ValueError(f"Spectra must have same length: {len(spectrum1)} vs {len(spectrum2)}")

    if method == 'correlation':
        # Pearson correlation coefficient
        correlation = np.corrcoef(spectrum1, spectrum2)[0, 1]
        return correlation

    elif method == 'rmse':
        # Root mean squared error (lower is better)
        rmse = np.sqrt(np.mean((spectrum1 - spectrum2)**2))
        return rmse

    elif method == 'mae':
        # Mean absolute error (lower is better)
        mae = np.mean(np.abs(spectrum1 - spectrum2))
        return mae

    elif method == 'r2':
        # R-squared coefficient of determination
        ss_res = np.sum((spectrum1 - spectrum2)**2)
        ss_tot = np.sum((spectrum1 - np.mean(spectrum1))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        return r2

    else:
        logger.warning(f"Unknown similarity method '{method}', using correlation")
        return calculate_spectrum_similarity(spectrum1, spectrum2, 'correlation')


# Export all functions
__all__ = [
    'gaussian_broadening',
    'calculate_stack_offsets',
    'raman_activity_to_intensity',
    'calculate_absorbance_transmittance',
    'apply_frequency_scaling',
    'get_color_palette',
    'get_ir_region_boundaries',
    'get_ir_region_labels',
    'get_raman_region_boundaries',
    'get_raman_region_labels',
    'normalize_spectrum',
    'interpolate_spectrum',
    'align_spectra_to_common_grid',
    'calculate_spectrum_similarity',
    'find_optimal_label_positions'
]
