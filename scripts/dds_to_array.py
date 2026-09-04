"""Decode a DDS texture for the embedded C++ image bridge."""

import imageio.v2 as imageio
import numpy as np


def dds_to_array(input_dds_file):
    reader = imageio.get_reader(input_dds_file, format="DDS")
    try:
        image = np.asarray(reader.get_data(0))
    finally:
        reader.close()

    if image.ndim not in (2, 3):
        raise ValueError(f"DDS image must have two or three dimensions, got {image.ndim}.")

    height, width = image.shape[:2]
    channels = image.shape[2] if image.ndim == 3 else 1
    if width <= 0 or height <= 0 or channels < 1 or channels > 4:
        raise ValueError(f"DDS image has an unsupported shape: {image.shape}.")
    if image.dtype != np.uint8:
        raise TypeError(f"DDS image must decode to uint8, got {image.dtype}.")

    contiguous_image = np.ascontiguousarray(image)
    return width, height, channels, contiguous_image.tobytes(order="C")
