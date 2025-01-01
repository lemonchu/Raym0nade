import os

def hdr_to_array(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    if not os.access(filename, os.R_OK):
        raise PermissionError(f"File is not readable: {filename}")

    try:
        import OpenEXR
        import Imath
        exr_file = OpenEXR.InputFile(filename)
        dw = exr_file.header()['dataWindow']
        width = dw.max.x - dw.min.x + 1
        height = dw.max.y - dw.min.y + 1

        FLOAT = Imath.PixelType(Imath.PixelType.FLOAT)
        redstr = exr_file.channel('R', FLOAT)
        greenstr = exr_file.channel('G', FLOAT)
        bluestr = exr_file.channel('B', FLOAT)

        import numpy as np
        red = np.frombuffer(redstr, dtype=np.float32)
        green = np.frombuffer(greenstr, dtype=np.float32)
        blue = np.frombuffer(bluestr, dtype=np.float32)

        data = np.stack([red, green, blue], axis=-1).reshape(height, width, 3)

        return width, height, data

    except Exception as e:
        raise RuntimeError(f"Failed to load HDR image: {e}")

if __name__ == "__main__":
    filename = "model/Bistro_v5_2/san_giuseppe_bridge_4k.hdr"
    try:
        width, height, data = hdr_to_array(filename)
        print(f"Loaded HDR image with dimensions: {width}x{height}")
    except Exception as e:
        print(e)