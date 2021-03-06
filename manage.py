# Classes and functions for handling intake and management of glacier terminus data.

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString


class Glacier:
    def __init__(self, gid):
        self.gid = gid
        self.refline = LineString()
        self.refbox = LineString()
        self.centerline = LineString()
        self.officialname = ""
        self.unofficialname = ""
        self.fjordname = ""
        self.armname = ""
        self.obsseries = []
        self.areas = []
        self.termareas = []
        self.lengths = []
        self.dates = []

    def add_observation(self, observation):
        if observation.gid != self.gid:
            print('Cannot add glacier %s observation to glacier %s observation series' % (observation.gid, self.gid))
        self.obsseries.append(observation)

    def sort_by_date(self):
        self.obsseries = sorted(self.obsseries, key=lambda k: k.date)
    
    def extract(self, attr):
        """Extract list of a given attribute for each TerminusObservation in the Glacier observation series"""
        data_list = eval("[obs.%s for obs in self.obsseries]" % attr)
        data_list = pd.Series(data_list)
        return data_list


class TerminusObservation:
    def __init__(self, gid, qflag, termination, imageid, sensor, date, circadiandate, year, season, geometry):
        # Attributes that must be defined on instantiation
        self.gid = gid
        self.qflag = qflag
        self.termination = termination
        self.imageid = imageid
        self.sensor = sensor
        self.date = pd.to_datetime(date)
        self.circadiandate = circadiandate
        self.year = year
        self.season = season
        self.geometry = geometry
        # Attributes that are determined from initial instance attributes
        self.hydroyear = hydrologicalYear(self.date)
        self.dayofhydroyear = dayOfHydroyear(self.date)
        # Attributes that are calculated elsewhere...
        self.area = 0.0
        self.termarea = 0.0
        self.centerlineintersection = Point()
        self.length = 0.0


def shp2gdf(file, epsg=3574):
    """Reads a shapefile of glacier data and reprojects to a specified EPSG (default is EPSG:3574 - WGS 84 / North Pole LAEA Atlantic, unit=meters)
    Result is a geodataframe containing the shapefile data.
    Also reindex the gdf so that index=GlacierID for direct selection by ID."""
    gdf = gpd.read_file(file)
    gdf = gdf.to_crs(epsg=epsg)
    # TODO: check whether glacier has multiple entries
    # (then will have multiple indices of that value, which screws up indexing)
    gdf = gdf.sort_values(by='GlacierID').set_index('GlacierID', drop=False)
    return gdf


def glacierInfo(termini_gdf, box_gdf, gid):
    """Get terminus and metadata for a given glacier ID in the geodataframe of
    termini data, as well as the glacier's reference box."""
    terminus = termini_gdf.loc[gid]
    box = box_gdf.loc[gid]
    return terminus, box


def hydrologicalYear(date):
    """Determine hydrological year of a given date. Greenland hydrological year
    is defined as September 1 through August 31. Returns the starting year of
    the hydrological year (aligned with September)."""
    date = pd.to_datetime(date)
    if pd.notnull(date):
        if date.month >= 9:
            hydroyear = date.year
        elif date.month < 9:
            hydroyear = date.year - 1
        return hydroyear


def dayOfHydroyear(date):
    """Convert date to number of days since start of hydrological year,
    defined as September 1. Analagous to day-of-year."""
    date = pd.to_datetime(date)
    if pd.notnull(date):
        hydroyear = hydrologicalYear(date)
        start_hydroyear = pd.to_datetime('%s-09-01' % str(hydroyear))
        day_of_hydroyear = (date - start_hydroyear).days
        return day_of_hydroyear