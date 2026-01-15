"""
Linear Shifter
==============
Rigid shift with integer pixel offsets and zero-padding.
No wrap-around allowed.
"""

import numpy as np
from scipy.ndimage import shift as ndimage_shift
from typing import Union


def rigid_shift_image(
    image: np.ndarray,
    shift_x: int,
    shift_y: int,
    pad_value: float = 0.0,
    allow_wraparound: bool = False
) -> np.ndarray:
    """
    Rigidly shift a 2D image by integer pixel offsets using scipy.ndimage.shift.
    
    Mathematical Foundation:
    Shifted(x,y) = Original(x - Δx, y - Δy)
    Out-of-bounds pixels are zero-padded; circular wrap-around is forbidden.
    
    Args:
        image: 2D numpy array to shift
        shift_x: Integer shift in x direction (columns, positive = right)
        shift_y: Integer shift in y direction (rows, positive = down)
        pad_value: Value for out-of-bounds pixels (default: 0.0)
        allow_wraparound: If False (default), forbid circular shifts
    
    Returns:
        Shifted 2D array, same shape as input
    
    Raises:
        ValueError: If wrap-around would occur and allow_wraparound=False
    """
    # Validate input
    if not isinstance(image, np.ndarray):
        image = np.array(image)
    
    if image.ndim != 2:
        raise ValueError(f"Expected 2D array, got {image.ndim}D")
    
    # Convert to integers
    shift_x = int(shift_x)
    shift_y = int(shift_y)
    
    # Check for wrap-around
    if not allow_wraparound:
        if abs(shift_x) >= image.shape[1] or abs(shift_y) >= image.shape[0]:
            raise ValueError(
                f"Shift ({shift_x}, {shift_y}) would cause wrap-around for "
                f"image shape {image.shape}. Use smaller shifts or set allow_wraparound=True."
            )
    
    # scipy.ndimage.shift with order=0 (nearest-neighbor, no interpolation)
    # Note: scipy uses (row, col) ordering, so we pass (shift_y, shift_x)
    shifted = ndimage_shift(
        image,
        shift=(shift_y, shift_x),
        order=0,  # Nearest neighbor (integer shifts)
        cval=pad_value,
        mode='constant'  # Zero-padding
    )
    
    return shifted

