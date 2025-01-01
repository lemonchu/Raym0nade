import imageio.v2 as imageio
import numpy as np

def hdr_to_array(input_hdr_file):
    try:
        hdr_image = imageio.imread(input_hdr_file, format='HDR-FI')

        width = hdr_image.shape[1]
        height = hdr_image.shape[0]
        channels = hdr_image.shape[2] if len(hdr_image.shape) > 2 else 1

        # Convert the image to a float32 array
        hdr_array = hdr_image.astype(np.float32)

        return width, height, channels, hdr_array
    except Exception as e:
        print(f"Conversion failed: {e}")
        return None, None, None, None