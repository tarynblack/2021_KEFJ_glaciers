#!/usr/bin/python3

# Taryn Black, August 2018

import argparse
import os
from subprocess import call, check_output


def get_args():
    """Usage statement."""
    parser = argparse.ArgumentParser(description='Download and reproject \
                                     Landsat scenes from Google Cloud')
    parser.add_argument('--productID', help='USGS Product ID for scene')
    parser.add_argument('--basedir', help='Base directory to store image',
                        default='/Volumes/insar5/teblack/data/LANDSAT/')
    parser.add_argument('--epsg', help='Coordinate reference system',
                        default='3413')
    args = parser.parse_args()
    return args


def create_scene_path(args):
    """Create Google Cloud Platform path for the input scene.

    The GCP directory structure for Landsat Collection 1 scenes is:
    gs://gcp-public-data-landsat/[SENSOR]/01/[PATH]/[ROW]/[PRODUCT_ID]

    Note: as of August 2018, some USGS products (e.g. GloVis) still
    output a Scene ID, which is deprecated and has been replaced by
    Product ID as part of the USGS transition away from Pre-Collection
    datasets. GCP still stores imagery with a Scene ID under a different
    directory structure.
    """
    pid = args.productID
    sensor = pid[:4]
    path = pid[10:13]
    row = pid[13:16]
    GCP_path = 'gs://gcp-public-data-landsat/' + sensor + '/01/' + path + '/' \
               + row + '/' + pid
    return GCP_path


def dest_directory(args):
    """Create the image's destination directory, if it does not exist.

    The directory is named for the hydrological year, e.g. 1996-1997/
    where the hydrological year runs from September 1 through August 31.
    """
    pid = args.productID
    basedir = args.basedir
    year = pid[17:21]
    month = pid[21:23]
    if 9 <= int(month) <= 12:
        dirname = year + '-' + str(int(year)+1)
    elif 1 <= int(month) <= 8:
        dirname = str(int(year)-1) + '-' + year
    destdir = basedir #+ dirname
    if not os.path.isdir(destdir):
        os.mkdir(destdir)
    return destdir


def choose_bands(args):
    """Choose which bands to download, depending on the sensor.

    Landsat 7-8: band 8 (panchromatic), band 4 (red), band 3 (green), band 2 (blue)
    Landsat 4-5 TM: band 4 (NIR), band 3 (red), band 2 (green)
    Landsat 4-5 MSS: band 3? (NIR), band 2 (red), band 1 (green)
    Landsat 1-3 MSS: band 6? (NIR), band 5 (red), band 4 (green)
    """
    pid = args.productID
    sensor = pid[:4]
    if sensor in {'LC08'}:
        bands = ['B8', 'B4', 'B3', 'B2']
    elif sensor in {'LE07'}:
        bands = ['B8', 'B3', 'B2', 'B1']
    elif sensor in {'LT05', 'LT04'}:
        bands = ['B3', 'B2', 'B1']
    elif sensor in {'LM05', 'LM04'}:
        bands = ['B3', 'B2', 'B1']
    elif sensor in {'LM03', 'LM02', 'LM01'}:
        bands = ['B6', 'B5', 'B4']
    return bands


def reproject_scene(args, bandfile, destdir):
    """Use GDAL to reproject Landsat scene.

    These commands are copied from Ian Joughin (getgoogle.py).
    """
    proj = str(check_output('gdalsrsinfo -o proj4 "EPSG:' + args.epsg + '"', shell=True), 'utf-8'
               ).strip()
    cmdstr = 'gdalwarp -srcnodata 0 -tr 15.0 15.0 -dstnodata 0 ' \
             '-co "COMPRESS=DEFLATE" -co "TILED=YES" -co "ZLEVEL=9" ' \
             '-co "PREDICTOR=2" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" ' \
             '-r cubic -t_srs ' + proj + ' ' \
             + args.basedir + '/temp/' + bandfile + ' ' \
             + destdir + '/' + bandfile
    print(cmdstr)
    call(cmdstr, shell=True)
    smallfile = 'small_' + bandfile
    cmdstr = 'gdalwarp -srcnodata 0 -tr 100.0 100.0 -dstnodata 0 ' \
             '-co "COMPRESS=DEFLATE" -co "TILED=YES" -co "ZLEVEL=9" ' \
             '-co "PREDICTOR=2" -co "BLOCKXSIZE=512" -co "BLOCKYSIZE=512" ' \
             '-r near -t_srs ' + proj + ' ' \
             + args.basedir + '/temp/' + bandfile + ' ' \
             + destdir + '/' + smallfile + ' ; ' \
             + 'gdaladdo -r average ' + destdir + '/' + smallfile + \
             ' 2 4 8 16 ; '
    call(cmdstr, shell=True)


def download_scene(bands, args, GCP_path, destdir):
    """Download Landsat scene bands and metadata from Google Cloud Platform.

    The scene is downloaded to a temp directory and reprojected, then moved to
    a directory based on its hydrological year.
    """
    for band in bands:
        bandfile = args.productID + '_' + band + '.TIF'
        command = 'gsutil cp ' + GCP_path + '/' + bandfile + ' ' \
                  + args.basedir + '/temp/'
        call(command, shell=True)
        reproject_scene(args, bandfile, destdir)
    command = 'gsutil cp ' + GCP_path + '/' + args.productID + '_MTL.txt ' \
              + destdir + ' ; ' \
              + 'rm ' + args.basedir + '/temp/*'
    call(command, shell=True)


def main():
    args = get_args()
    GCP_path = create_scene_path(args)
    destdir = dest_directory(args)
    bands = choose_bands(args)
    download_scene(bands, args, GCP_path, destdir)


main()
