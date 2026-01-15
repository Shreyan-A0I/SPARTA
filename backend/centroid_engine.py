"""
Centroid Engine
===============
Dual-mode centroid calculation with SNR-based noise floor filtering.

Methods:
- Intensity-Weighted: Emphasizes high-intensity regions
- Binary-Mask: Geometric center, noise-independent
"""

import numpy as np
from typing import Tuple, Dict, Optional


def calculate_centroid(
    intensity_map: np.ndarray,
    method: str = 'weighted',
    snr_floor_multiplier: float = 3.0,
    binary_threshold: Optional[float] = None
) -> Tuple[Optional[float], Optional[float], Dict]:
    """
    Calculate centroid using intensity-weighted or binary-mask method.
    
    Mathematical Foundation:
    - SNR Floor: I_floor = μ_noise + k * σ_noise
    - Weighted: centroid_x = Σ(x * I(x,y)) / Σ(I(x,y))
    - Binary: centroid_x = mean(x) where I(x,y) > threshold
    
    Args:
        intensity_map: 2D numpy array of normalized intensities
        method: 'weighted' (intensity-weighted) or 'binary' (geometric center)
        snr_floor_multiplier: SNR multiplier for noise floor (default: 3.0)
        binary_threshold: Custom threshold for binary method (if None, use SNR floor)
    
    Returns:
        (centroid_x, centroid_y, metadata_dict)
        Returns (None, None, {'error': '...'}) if calculation fails
    """
    # Validate input
    if intensity_map.size == 0:
        return None, None, {'error': 'Empty intensity map'}
    
    if not isinstance(intensity_map, np.ndarray):
        intensity_map = np.array(intensity_map)
    
    # Compute SNR floor
    noise_percentile = np.percentile(intensity_map, 10)
    noise_mask = intensity_map <= noise_percentile
    
    if np.any(noise_mask):
        noise_values = intensity_map[noise_mask]
        mu_noise = np.mean(noise_values)
        sigma_noise = np.std(noise_values)
        
        # Handle edge case where std is zero
        if sigma_noise == 0:
            sigma_noise = 1e-9
        
        snr_floor = mu_noise + snr_floor_multiplier * sigma_noise
    else:
        # Fallback if no noise pixels found
        snr_floor = np.min(intensity_map)
    
    # Apply SNR floor
    signal_mask = intensity_map > snr_floor
    intensity_filtered = intensity_map.copy()
    intensity_filtered[~signal_mask] = 0.0
    
    # Calculate centroid
    y_coords, x_coords = np.indices(intensity_map.shape)
    
    if method == 'weighted':
        total_intensity = np.sum(intensity_filtered)
        
        if total_intensity <= 0:
            return None, None, {
                'error': 'No signal above SNR floor',
                'snr_floor': snr_floor,
                'signal_pixels': 0
            }
        
        centroid_x = np.sum(x_coords * intensity_filtered) / total_intensity
        centroid_y = np.sum(y_coords * intensity_filtered) / total_intensity
        method_name = "Intensity-Weighted"
        
    elif method == 'binary':
        if binary_threshold is None:
            binary_threshold = snr_floor
        
        tissue_mask = intensity_filtered > binary_threshold
        num_tissue_pixels = np.sum(tissue_mask)
        
        if num_tissue_pixels <= 0:
            return None, None, {
                'error': 'No tissue pixels above threshold',
                'snr_floor': snr_floor,
                'threshold': binary_threshold,
                'signal_pixels': 0
            }
        
        centroid_x = np.mean(x_coords[tissue_mask])
        centroid_y = np.mean(y_coords[tissue_mask])
        method_name = "Binary-Mask"
    else:
        raise ValueError(f"Unknown method: {method}. Must be 'weighted' or 'binary'")
    
    # Build metadata
    metadata = {
        'method': method_name,
        'snr_floor': float(snr_floor),
        'signal_pixels': int(np.sum(signal_mask)),
        'total_pixels': int(intensity_map.size),
        'coverage_pct': 100.0 * np.sum(signal_mask) / intensity_map.size,
        'centroid_x': float(centroid_x),
        'centroid_y': float(centroid_y)
    }
    
    return float(centroid_x), float(centroid_y), metadata

