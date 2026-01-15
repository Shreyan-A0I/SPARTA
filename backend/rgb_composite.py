"""
RGB Composite Overlay
=====================
Per-channel normalized RGB overlays for metabolite co-visualization.
"""

import numpy as np
from typing import Tuple


def create_rgb_composite(
    intensity_a: np.ndarray,
    intensity_b: np.ndarray,
    percentile_clip: float = 99.0
) -> np.ndarray:
    """
    Create an RGB composite overlay with per-channel normalization.
    
    Mathematical Foundation:
    Normalized_A = (I_A - I_A_min) / (I_A_max - I_A_min)
    Normalized_B = (I_B - I_B_min) / (I_B_max - I_B_min)
    Overlay = [Normalized_A, 0, Normalized_B] (Red=A, Green=0, Blue=B)
    
    Args:
        intensity_a: 2D array of metabolite A intensities
        intensity_b: 2D array of metabolite B intensities (after shift)
        percentile_clip: Percentile for intensity clipping (default: 99th)
    
    Returns:
        RGB composite array of shape (height, width, 3) with values in [0, 1]
    """
    # Validate inputs
    if intensity_a.shape != intensity_b.shape:
        raise ValueError(
            f"Intensity arrays must have same shape. "
            f"Got {intensity_a.shape} and {intensity_b.shape}"
        )
    
    # Clip extreme outliers using percentile
    a_nonzero = intensity_a[intensity_a > 0]
    b_nonzero = intensity_b[intensity_b > 0]
    
    if len(a_nonzero) > 0:
        a_clip = np.percentile(a_nonzero, percentile_clip)
    else:
        a_clip = 1.0
    
    if len(b_nonzero) > 0:
        b_clip = np.percentile(b_nonzero, percentile_clip)
    else:
        b_clip = 1.0
    
    intensity_a_clipped = np.clip(intensity_a, 0, a_clip)
    intensity_b_clipped = np.clip(intensity_b, 0, b_clip)
    
    # Normalize each channel independently
    a_max = np.max(intensity_a_clipped)
    b_max = np.max(intensity_b_clipped)
    
    # Avoid division by zero
    normalized_a = intensity_a_clipped / (a_max + 1e-9)
    normalized_b = intensity_b_clipped / (b_max + 1e-9)
    
    # Stack into RGB: Red = A, Green = 0, Blue = B
    rgb_composite = np.stack([
        normalized_a,
        np.zeros_like(normalized_a),
        normalized_b
    ], axis=2)
    
    # Ensure values are in [0, 1]
    rgb_composite = np.clip(rgb_composite, 0.0, 1.0)
    
    return rgb_composite


def create_rgb_composite_with_green_blend(
    intensity_a: np.ndarray,
    intensity_b: np.ndarray,
    percentile_clip: float = 99.0
) -> np.ndarray:
    """
    Alternative RGB composite with green channel showing co-enrichment.
    
    Red = A, Green = Min(A, B), Blue = B
    Result: Magenta (A+B), White (strong co-enrichment)
    
    Args:
        intensity_a: 2D array of metabolite A intensities
        intensity_b: 2D array of metabolite B intensities (after shift)
        percentile_clip: Percentile for intensity clipping (default: 99th)
    
    Returns:
        RGB composite array of shape (height, width, 3) with values in [0, 1]
    """
    # Validate inputs
    if intensity_a.shape != intensity_b.shape:
        raise ValueError(
            f"Intensity arrays must have same shape. "
            f"Got {intensity_a.shape} and {intensity_b.shape}"
        )
    
    # Clip extreme outliers
    a_nonzero = intensity_a[intensity_a > 0]
    b_nonzero = intensity_b[intensity_b > 0]
    
    if len(a_nonzero) > 0:
        a_clip = np.percentile(a_nonzero, percentile_clip)
    else:
        a_clip = 1.0
    
    if len(b_nonzero) > 0:
        b_clip = np.percentile(b_nonzero, percentile_clip)
    else:
        b_clip = 1.0
    
    intensity_a_clipped = np.clip(intensity_a, 0, a_clip)
    intensity_b_clipped = np.clip(intensity_b, 0, b_clip)
    
    # Normalize each channel independently
    a_max = np.max(intensity_a_clipped)
    b_max = np.max(intensity_b_clipped)
    
    normalized_a = intensity_a_clipped / (a_max + 1e-9)
    normalized_b = intensity_b_clipped / (b_max + 1e-9)
    
    # Green channel shows co-enrichment (minimum of both)
    green_channel = np.minimum(normalized_a, normalized_b)
    
    # Stack into RGB: Red = A, Green = Min(A,B), Blue = B
    rgb_composite = np.stack([
        normalized_a,
        green_channel,
        normalized_b
    ], axis=2)
    
    # Ensure values are in [0, 1]
    rgb_composite = np.clip(rgb_composite, 0.0, 1.0)
    
    return rgb_composite

