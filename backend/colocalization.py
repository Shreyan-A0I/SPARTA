"""
Co-localization Map Computation
================================
Pixel-wise co-localization index calculation with safety masking.
"""

import numpy as np
from typing import Literal


def compute_colocalization_map(
    intensity_a: np.ndarray,
    intensity_b: np.ndarray,
    min_b_intensity: float = 10.0,
    handle_zeros: Literal['nan', 'zero'] = 'nan'
) -> np.ndarray:
    """
    Compute pixel-wise co-localization index (A/B) with safety masking.
    
    Mathematical Foundation:
    Co-loc(x,y) = I_A(x,y) / I_B(x,y) if I_B(x,y) >= I_min and I_B(x,y) > 0
                 = NaN otherwise
    
    Args:
        intensity_a: 2D array of metabolite A intensities
        intensity_b: 2D array of metabolite B intensities (after alignment shift)
        min_b_intensity: Minimum threshold for denominator (default: 10.0)
        handle_zeros: 'nan' (default) or 'zero' for invalid pixels
    
    Returns:
        Co-localization map (2D array), same shape as inputs
    """
    # Validate inputs
    if intensity_a.shape != intensity_b.shape:
        raise ValueError(
            f"Intensity arrays must have same shape. "
            f"Got {intensity_a.shape} and {intensity_b.shape}"
        )
    
    # Create valid mask: B must be above threshold and non-zero
    valid_mask = (intensity_b >= min_b_intensity) & (intensity_b > 0)
    
    # Initialize result array
    colocalization = np.full_like(intensity_a, np.nan, dtype=float)
    
    # Compute ratio only for valid pixels
    colocalization[valid_mask] = intensity_a[valid_mask] / intensity_b[valid_mask]
    
    # Handle invalid pixels
    if handle_zeros == 'zero':
        colocalization[~valid_mask] = 0.0
    # else: keep as NaN (default)
    
    return colocalization

