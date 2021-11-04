#!/usr/bin/python3

# Taryn Black, June 2020

import argparse
from subprocess import call

#%% Get user-specified parameters
def get_args():
    """Usage statement."""
    parser = argparse.ArgumentParser(description='Create RGB composite of \
        Landsat band TIFs, and pansharpen if available.')
    parser.add_argument('--id', help='USGS Product ID for scene, e.g. \
        LC08_L1GT_068018_20170910_20170927_01_T2')
    parser.add_argument('--path', help='Path to directory to store image',
                        default='/Volumes/insar5/teblack/data/LANDSAT/')
    args = parser.parse_args()
    return args

args = get_args()
landsat_id = args.id
image_path = args.path

#%% Determine which bands to use
sensor = landsat_id[0:4]

rgb_bands = {'LC08': [4, 3, 2, 8], # Landsat-8: 4-red, 3-green, 2-blue, 8-pan
             'LE07': [3, 2, 1, 8], # Landsat-7: 3-red, 2-green, 1-blue, 8-pan
             'LT05': [3, 2, 1],    # Landsat-4,5 TM: 3-red, 2-green, 1-blue
             'LT04': [3, 2, 1],    # Landsat-4,5 TM: 3-red, 2-green, 1-blue
             'LM05': [3, 2, 1],    # Landsat-4,5 MSS: 3-NIR, 2-red, 1-green
             'LM04': [3, 2, 1],    # Landsat-4,5 MSS: 3-NIR, 2-red, 1-green
             'LM03': [6, 5, 4],    # Landsat-1,2,3 MSS: 6-NIR, 5-red, 4-green
             'LM02': [6, 5, 4],    # Landsat-1,2,3 MSS: 6-NIR, 5-red, 4-green
             'LM01': [6, 5, 4]     # Landsat-1,2,3 MSS: 6-NIR, 5-red, 4-green
}

sensor_bands = rgb_bands.get(sensor)

#%% Get RGB file paths
red = image_path+landsat_id+'_B{}.TIF'.format(sensor_bands[0])
green = image_path+landsat_id+'_B{}.tif'.format(sensor_bands[1])
blue = image_path+landsat_id+'_B{}.tif'.format(sensor_bands[2])

#%% Create composite image, and pansharpen if pan band exists
if len(sensor_bands) == 4:
    print('Panchromatic band present. Compositing bands {}, {}, and {} and pansharpening with band {}...'.format(
            sensor_bands[0], sensor_bands[1], sensor_bands[2], sensor_bands[3]))
    pan = image_path+landsat_id+'_B{}.tif'.format(sensor_bands[3])
    cmd = 'gdal_pansharpen.py {} {} {} {} {}_pan-composite.TIF'.format(
        pan, red, green, blue, image_path+landsat_id)
    call(cmd, shell=True)
    print('Created {}_pan-composite.TIF'.format(image_path+landsat_id))
else:
    print('No panchromatic band present. Compositing bands {}, {}, and {}...'.format(sensor_bands[0], sensor_bands[1], sensor_bands[2]))
    cmd = 'gdalbuildvrt -separate {}_composite.vrt {} {} {}'.format(
        image_path+landsat_id, red, green, blue)
    call(cmd, shell=True)
    cmd = 'gdal_translate {}_composite.vrt {}_composite.TIF'.format(
        image_path+landsat_id, image_path+landsat_id)
    call(cmd, shell=True)
    print('Created {}_composite.TIF'.format(image_path+landsat_id))
