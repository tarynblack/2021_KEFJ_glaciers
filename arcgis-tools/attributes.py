# Functions to evaluate geodatabase attributes based on Landsat Product ID
# T Black, 24 July 2020

from datetime import datetime


def getSensor(id):
    sensor = id[0:4]
    return sensor

def getTileCoordinate(id):
    tile_coordinate = id[10:16]
    return tile_coordinate

def getReferenceSystem(id):
    sensor_number = int(id[2:4])
    if sensor_number > 3:
        reference_system = 'WRS-2'
    elif sensor_number <= 3:
        reference_system = 'WRS-1'
    return reference_system

def getSourceDate(id):
    year = id[17:21]
    month = id[21:23]
    day = id[23:25]
    source_date = year + '-' + month + '-' + day
    return source_date

def getCircadianDate(id):
    year = int(id[17:21])
    month = int(id[21:23])
    day = int(id[23:25])
    day_of_year = datetime(year, month, day).timetuple().tm_yday
    return day_of_year

def getYear(id):
    year = id[17:21]
    return year

def getSeason(id):
    month = int(id[21:23])
    day = int(id[23:25])
    if month in [4, 5, 6]:
        season = 'SPR'
    elif month in [7]:
        season = 'SUM'
    elif month in [8]:
        if day < 17:
            season = 'SUM'
        else:
            season = 'AUT'
    elif month in [9]:
        season = 'AUT'
    elif month in [10]:
        if day < 17:
            season = 'AUT'
        else:
            season = 'WIN'
    elif month in [11, 12, 1, 2, 3]:
        season = 'WIN'
    return season