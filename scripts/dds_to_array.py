import imageio.v2 as imageio

def dds_to_array(input_dds_file):
    try:
        reader = imageio.get_reader(input_dds_file, format='DDS')
        dds_image = reader.get_data(0)

        width = dds_image.shape[1]
        height = dds_image.shape[0]
        channels = dds_image.shape[2] if len(dds_image.shape) > 2 else 1

        return width, height, channels, dds_image.tobytes()
    except Exception as e:
        print(f"Conversion failed: {e}")
        return None, None, None, None