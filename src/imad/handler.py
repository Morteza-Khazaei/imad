import os
import sys
import re
import ray
import argparse
import logging
from osgeo import gdal
import numpy as np

from .core import IRMAD



class  imadHandler:

    version = 2023.08


    def __init__(self, input, output, **kwargs) -> None:

        self.logger = kwargs.get("logger", logging.getLogger("root"))
        self.logger.info(f'======================= This is IRMAD v{self.version} ======================')

        self.input_base_dir = input
        self.output_base_dir = output
        
        return None
        

    @staticmethod
    def init_loggers(msg_level=logging.DEBUG):
        """
        Init a stdout logger
        :param msg_level: Standard msgLevel for both loggers. Default is DEBUG
        """

        logging.getLogger().addHandler(logging.NullHandler())
        # Create default path or get the pathname without the extension, if there is one
        imad_logger = logging.getLogger("root")
        imad_logger.handlers = []  # Remove the standard handler again - Bug in logging module

        imad_logger.setLevel(msg_level)
        formatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s] %(message)s")

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        imad_logger.addHandler(console_handler)
        return imad_logger


    def create_geotiff(self, base_dir, input_files):
        # Create a list to store the individual bands
        bands = []

        # Open each input file, extract the band, and append it to the bands list
        for file in input_files:
            file_path = os.path.join(base_dir, file)
            self.logger.info(f'Read raster band with Id: {file}')
            ds = gdal.Open(file_path, gdal.GA_ReadOnly)
            cols = ds.RasterXSize
            rows = ds.RasterYSize
            geotransform = ds.GetGeoTransform()
            projection = ds.GetProjection()

            if ds is not None:
                band = ds.GetRasterBand(1)  # Get the first (and only) band
                nodata = band.GetNoDataValue()
                band_array = band.ReadAsArray()
                band_array[band_array == np.nan] = -10000
                bands.append(band_array)
                ds = None  # Close the dataset
        # Create NRGB name
        fname = re.sub(r'B\d+', 'NRGB', file)
        output_file = os.path.join(base_dir, fname)

        if not os.path.exists(output_file):
            self.logger.info(f'******* Write raster file with id: {fname}.')
            # Create a new dataset for the NRGB bands
            driver = gdal.GetDriverByName("GTiff")
            ds_out = driver.Create(output_file, cols, rows, len(bands), band.DataType, options=['COMPRESS=LZW'])
            ds_out.SetGeoTransform(geotransform)
            if projection is not None:
                ds_out.SetProjection(projection)

            # Loop through the bands and write them to the output dataset
            for i, band in enumerate(bands):
                output_band = ds_out.GetRasterBand(i + 1)  # Band numbers start from 1
                output_band.WriteArray(band)
                output_band.SetNoDataValue(nodata)

            # Close the output dataset
            ds_out = None

        return output_file


    def execute(self):
        master = None
        mad_instance_list = []

        tiles = os.listdir(self.input_base_dir)
        for tile in tiles:
            out_dir = os.path.join(self.output_base_dir, tile)
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            
            l3a_tile = os.path.join(self.input_base_dir, tile)
            # Keep products dirs
            l3a_dirs = [d for d in os.listdir(l3a_tile) if os.path.isdir(os.path.join(l3a_tile, d))]
            # Using a lambda function within the sorted function
            extract_date = lambda name: re.search(r'_(\d{8})-', name).group(1) if re.search(r'_(\d{8})-', name) else ""
            l3a_products = sorted(l3a_dirs, key=extract_date)
            for l3a in l3a_products:
                l3a_path = os.path.join(l3a_tile, l3a)

                # List all files in the directory
                all_files = os.listdir(l3a_path)

                # Filter files that have B2, B3, B4, B8
                NRGB_bands = [fname for fname in all_files if any(band in fname for band in ['B2', 'B3', 'B4', 'B8.'])]

                NRGB_file = self.create_geotiff(l3a_path, NRGB_bands)
                
                if master is not None:
                    head, tail = os.path.split(NRGB_file)
                    pname = tail.replace('NRGB', 'CHMAP')
                    chmap_path = os.path.join(out_dir, pname)
                    self.logger.info(f'-----> Check if CHMAP raster with id: {pname} exists.')
                    if not os.path.isfile(chmap_path):
                        self.logger.info(f'******* Master image is: {master}.')
                        self.logger.info(f'******* Slave image is: {NRGB_file}.')
                        self.logger.info(f'******* Write CHMAP raster file with id: {pname}.')
                    
                        # Create an actor process.
                        imad = IRMAD.remote(master=master, slave=NRGB_file, output=out_dir, 
                                            filename=pname, penalization=0.001, logger=self.logger)

                        mad_instance_list.append(imad.MAD_iteration.remote())
                
                master = NRGB_file
            master = None
        
        return mad_instance_list


def main():

    parser = argparse.ArgumentParser(description="Perfrom IR-MAD change detection on bitemporal, multispectral imagery.")

    parser.add_argument("-i", "--input", type=str, help="Where did you store L3A products?")
    parser.add_argument("-o", "--output", type=str, help="Where is your preferred direction to store change maps?")
    parser.add_argument("-n", "--cups", type=int, help="Number of cores to be utilized?", default=12)
    parser.add_argument("-v", "--verbose", help="Provides detailed (DEBUG) logging for Maja. Default is false",
                        default=False, action="store_true")

    args = parser.parse_args()

    # TODO Add error skipping
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logger = imadHandler.init_loggers(msg_level=logging_level)

    # Start Ray.
    ray.init(num_cpus=args.cups, include_dashboard=False)
    logger.info(f'Ray is started: {ray.is_initialized()}')

    imad = imadHandler(args.input, args.output, logger=logger)
    mad_instances = imad.execute()
    chunks = [mad_instances[x:x+5] for x in range(0, len(mad_instances), 5)]

    for chunk in chunks:
        try:
            ray.get(chunk)
            ready, not_ready = ray.wait(chunk, num_returns=len(chunk))
            logger.info(ready, not_ready)
        except ValueError as e:
            logger.info(f'This {chunk} could not be executed because of this error: {e}')


    ray.shutdown()
    logger.info(f'Ray is shutdown: {ray.is_initialized()}')