"""
Heatmap Reconstruction from imzML
==================================
Reconstruct 2D intensity heatmaps for specific m/z values from imzML data.
"""

import numpy as np
from typing import Union

try:
    from pyimzml.ImzMLParser import ImzMLParser
except ImportError:
    ImzMLParser = None  # For testing without pyimzml installed


def reconstruct_heatmap(
    parser,
    target_mz: float,
    tolerance: float = 0.5
) -> np.ndarray:
    """
    Reconstruct a 2D intensity heatmap for a specific m/z from imzML data.
    
    Args:
        parser: PyImzMLParser object (already initialized) or mock object for testing
        target_mz: Target m/z value
        tolerance: m/z tolerance (±Da) for binning (default: 0.5)
    
    Returns:
        2D numpy array (heatmap) with shape (y_max+1, x_max+1)
    """
    if parser is None:
        raise ValueError("Parser is None. Please initialize ImzMLParser first.")
    
    if len(parser.coordinates) == 0:
        raise ValueError("No coordinates found in imzML file.")
    
    # Get coordinate ranges - handle different coordinate formats safely
    try:
        coords = np.array(parser.coordinates)
        # Handle both 2D and 3D coordinates (x, y) or (x, y, z)
        if coords.shape[1] >= 2:
            x_max = int(np.max(coords[:, 0]))
            y_max = int(np.max(coords[:, 1]))
        else:
            raise ValueError("Coordinates must have at least x and y values")
    except (IndexError, ValueError) as e:
        raise ValueError(f"Error parsing coordinates: {e}. Check imzML file format.")
    
    # Initialize heatmap
    heatmap = np.zeros((y_max + 1, x_max + 1), dtype=np.float32)
    
    # Iterate through all pixels
    for i, coord in enumerate(parser.coordinates):
        try:
            # Safely extract x, y coordinates
            if len(coord) >= 2:
                x, y = int(coord[0]), int(coord[1])
            else:
                continue  # Skip invalid coordinates
            
            mzs, intensities = parser.getspectrum(i)
            
            # Find peaks within tolerance
            mask = np.abs(mzs - target_mz) <= tolerance
            
            if np.any(mask):
                # Sum intensities within tolerance window
                heatmap[y, x] = np.sum(intensities[mask])
        except (IndexError, ValueError, TypeError) as e:
            # Skip problematic spectra
            continue
        except Exception as e:
            # Skip any other errors
            continue
    
    return heatmap

