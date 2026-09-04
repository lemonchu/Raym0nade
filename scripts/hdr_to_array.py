"""Decode an HDR environment for the embedded C++ image bridge."""

import imageio.v2 as imageio
import numpy as np


def hdr_to_array(input_hdr_file):
    image = np.asarray(imageio.imread(input_hdr_file, format="HDR-FI"))
    if image.ndim not in (2, 3):
        raise ValueError(f"HDR image must have two or three dimensions, got {image.ndim}.")

    height, width = image.shape[:2]
    channels = image.shape[2] if image.ndim == 3 else 1
    if width <= 0 or height <= 0 or channels < 1 or channels > 4:
        raise ValueError(f"HDR image has an unsupported shape: {image.shape}.")

    contiguous_image = np.ascontiguousarray(image, dtype=np.float32)
    return width, height, channels, contiguous_image
