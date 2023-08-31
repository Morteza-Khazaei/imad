import os
import re
import argparse
# import ray
from osgeo import gdal

from .core import IRMAD






def create_geotiff(base_dir, input_files):
    # Create a list to store the individual bands
    bands = []

    # Open each input file, extract the band, and append it to the bands list
    for file in input_files:
        ds = gdal.Open(file, gdal.GA_ReadOnly)
        if dataset is not None:
            band = dataset.GetRasterBand(1)  # Get the first (and only) band
            bands.append(band)
            dataset = None  # Close the dataset
    # Create NRGB name
    fname = re.sub(r'B\d+', 'NRGB', file)
    output_file = os.path.join(base_dir, fname)
    
    # Create a new dataset for the NRGB bands
    driver = gdal.GetDriverByName("GTiff")
    output_dataset = driver.Create(output_file, band[0].XSize, band[0].YSize, len(bands), bands[0].DataType)

    # Loop through the bands and write them to the output dataset
    for i, band in enumerate(bands):
        output_band = output_dataset.GetRasterBand(i + 1)  # Band numbers start from 1
        output_band.WriteArray(band.ReadAsArray())

    # Close the output dataset
    output_dataset = None

    return output_file



def main():

    parser = argparse.ArgumentParser(description="Perfrom IR-MAD change detection on bitemporal, multispectral imagery.")

    parser.add_argument("-i", "--input", type=str, help="Where did you store L3A products?")
    parser.add_argument("-o", "--output", type=str, help="Where is your preferred direction to store change maps?")
    parser.add_argument("-n", "--cups", type=int, help="Number of cores to be utilized?", default=12)

    args = parser.parse_args()

    # # Start Ray.
    # ray.init(num_cpus=args.cups, include_dashboard=False)
    # assert ray.is_initialized()

    input_base_dir = args.input
    output_base_dir = args.output

    perfix = 'CHMAP'
    master = None

    mad_instance_list = []

    tiles = os.listdir(input_base_dir)
    for tile in tiles:
        out_dir = os.path.join(output_base_dir, tile)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        l3a_tile = os.path.join(input_base_dir, tile)
        print('l3a_tile:', l3a_tile)
        # Keep products dirs
        l3a_dirs = [d for d in os.listdir(l3a_tile) if os.path.isdir(os.path.join(l3a_tile, d))]
        print('l3a_dirs:', l3a_dirs)
        # Using a lambda function within the sorted function
        extract_date = lambda name: re.search(r'_(\d{8})-', name).group(1) if re.search(r'_(\d{8})-', name) else ""
        l3a_products = sorted(l3a_dirs, key=extract_date)
        print('l3a_products:', l3a_products)
        for l3a in l3a_products:
            l3a_path = os.path.join(l3a_tile, l3a)

            # List all files in the directory
            all_files = os.listdir(l3a_path)

            # Filter files that have B2, B3, B4, B8
            NRGB_bands = [fname for fname in all_files if any(band in fname for band in ['B2', 'B3', 'B4', 'B8'])]

            NRGB_file = create_geotiff(l3a_path, NRGB_bands)
            
            if master:
                product_name = re.sub(r'NRGB\d+', perfix, NRGB_file)
                
                # Create an actor process.
                # imad = IRMAD.remote(master=master, slave=NRGB_file, output=output_base_dir, filename=product_name, penalization=0.001)
                imad = IRMAD(master=master, slave=NRGB_file, output=output_base_dir, filename=product_name, penalization=0.001)
                imad.MAD_iteration()

                # mad_instance_list.append(imad.MAD_iteration.remote())
            
            master = NRGB_file
        master = None

    # chunks = [mad_instance_list[x:x+5] for x in range(0, len(mad_instance_list), 5)]

    # for chunk in chunks:
    #     ray.get(chunk)
    #     ready, not_ready = ray.wait(chunk, num_returns=len(chunk))
    #     print(ready, not_ready)
    
    # ray.shutdown()
    # assert not ray.is_initialized()