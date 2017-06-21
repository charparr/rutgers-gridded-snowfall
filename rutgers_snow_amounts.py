#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:23:45 2017

@author: Charlie Parr
"""


import os
import pandas as pd
import geopandas as gpd
import numpy as np
from math import radians, sqrt, sin, cos, atan2
from shapely.geometry import Polygon

directory = '/home/cparr/workspace/rutgers_snow_lab/noram_1degv4_snow/'

datfnames = []
datfpaths = []

for root, dirs, files in os.walk(directory):
    for file in files:
        if file.endswith('.grid'):
            datfnames.append(file)
            datfpaths.append(root)

datfnames.sort()
datfpaths.sort()

fullpaths = []

for p, n in zip(datfpaths, datfnames):
    fullpaths.append(p + '/' + n)

fullpaths.sort()

lats = [i for i in np.arange(-125.0, -102.0)]
lons = [i for i in np.arange(31.0, 49.0)]
possible_corner_pts = [(x, y) for x in lons for y in lats]
possible_corner_pts = [str(y)+'      '+str(x) for x, y in possible_corner_pts]
corner_pts = []

with open(fullpaths[-1], 'r') as f:
    for line in f:
        for p in possible_corner_pts:    
            if p in line:
                corner_pts.append(p)

result_dct = dict()

del (datfnames, datfpaths, p, n, root, dirs, directory,
     file, files, line, lats, lons, possible_corner_pts)


def water_day(indate):
    """
    Function which turns calendar day of year into
    the water day of year. Water days sync with seasonal
    melt cycles in the Northern Hemisphere
    :param indate: (int) Calendar DOY
    :return: (int) Water DOY
    """
    doy = indate.timetuple().tm_yday
    if doy >= 274:
        outdate = doy - 274
    else:
        outdate = doy + 91
    return outdate


def geocalc(lat1, lon1, lat2, lon2):
    """
    Function to calculate greatest circle distance
    between two coordinates of latitude and longitude.
    Accounts for curve of Earth's surface.
    :param lat1:
    :param lon1:
    :param lat2:
    :param lon2:
    :return: Distance [km]
    """
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)

    dlon = lon1 - lon2
    earth_r = 6372.8

    y = sqrt(
        (cos(lat2) * sin(dlon)) ** 2
        + (cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)) ** 2
        )
    x = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon)
    c = atan2(y, x)
    return earth_r * c


def read_cell(loc):
    """
    Reading daily snowfall [mm] data for time series in
    1 deg. lat. by 1 deg. lon. grid cell. Each time step and value
    is stored in a DataFrame and then resampled by summing each winter.
    We then take the average of the winters.
    :param loc: (str) Lower left corner of grid cell e.g. '<-Lon>      <Lat>'
    :return: Mean Annual Winter Snowfall [cm]
    """
    vals = []
    keys = []
    # This part is pretty slow...
    for files_ in fullpaths:
        with open(files_, 'r') as f_:
            for line_ in f_:
                if loc in line_:
                    vals.append(float(line_[-6:-1]))
                    keys.append(files_[68:78])
            f_.close()

    ts_idx = pd.to_datetime(keys, format='%Y.%m.%d')

    df = pd.DataFrame(vals, index=ts_idx)
    df['Snow_cm'] = df[0] / 10.0  # data is snowfall in mm
    del df[0]
    df['day_of_year'] = df.index.strftime('%j')
    df['date'] = df.index
    df['water_doy'] = df['date'].apply(lambda x: water_day(x))
    winter = df.loc[df['water_doy'].isin(range(0, 200))]
    ann_winter = winter.resample('a').sum()
    recent_mean_ann_snow_cm = ann_winter['Snow_cm'].mean()
    return round(recent_mean_ann_snow_cm)


def calc_swe(mean_snow, ll_lat, ll_lon):
    """
    Calulate cell area, snow depths, SWE values, and water amounts.
    :param mean_snow: Mean Annual Winter Snowfall [cm]
    :param ll_lat: Lower Left Cell Corner Lat.
    :param ll_lon: Lower Left Cell Corner Lon.
    :return: Populates a results dictionary where each cell's LL coordinates are a key
    which corresponds to a dictionary with snow depth, SWE, Area, and Water Volume in
    both metric and English units.
    """
    rho = 0.1  # density value of snow

    grid_cell_area_sqkm = geocalc(ll_lat, ll_lon, ll_lat, ll_lon+1) * geocalc(ll_lat, ll_lon, ll_lat + 1, ll_lon)
    cell_acres_area = grid_cell_area_sqkm * 247.105
    mean_snow_d_ft = mean_snow * 0.0328084
    swe_cm = mean_snow * (rho / 1)
    swe_ft = swe_cm / 30.48
    swe_in = swe_ft * 12
    acre_feet_for_cell = round(swe_ft*cell_acres_area)

    cell_key = str(ll_lat)+' N ' + str(ll_lon) + ' W'
    result_dct[cell_key] = {}
    result_dct[cell_key]['Mean Snow Depth [ft]'] = round(mean_snow_d_ft, 1)
    result_dct[cell_key]['Mean Snow Depth [cm]'] = round(mean_snow)
    result_dct[cell_key]['Ave. SWE [in]'] = round(swe_in)
    result_dct[cell_key]['Ave. SWE [cm]'] = round(swe_cm)
    result_dct[cell_key]['Density [g/cm^3]'] = rho
    result_dct[cell_key]['Acre Feet'] = acre_feet_for_cell
    result_dct[cell_key]['Area [km^2]'] = round(grid_cell_area_sqkm)
    result_dct[cell_key]['Area [acres'] = round(cell_acres_area)
    result_dct[cell_key]['Water Volume [km^3]'] = acre_feet_for_cell*(1.23348 * 10**-6)

i = 0
for p in corner_pts:
    lllon = int(p[:4])
    lllat = int(p[-4:-2])
    calc_swe(read_cell(p), lllat, lllon)
    i += 1

df = pd.DataFrame.from_dict(result_dct)
df = df.T
df['ll'] = df.index
df['lly'] = df.ll.str[0:2].astype(int)
df['llx'] = df.ll.str[-6:-2].astype(int)
df['lry'] = df['lly']
df['lrx'] = df['llx']+1
df['ury'] = df['lly']+1
df['urx'] = df['lrx']
df['uly'] = df['lly']+1
df['ulx'] = df['llx']
df['ll'] = [xy for xy in zip(df['llx'], df['lly'])]
df['lr'] = [xy for xy in zip(df['lrx'], df['lry'])]
df['ur'] = [xy for xy in zip(df['urx'], df['ury'])]
df['ul'] = [xy for xy in zip(df['ulx'], df['uly'])]
df['close_geo'] = df['ll'].copy()

geo = []

for p in zip(df['ll'], df['ul'], df['ur'], df['lr'], df['close_geo']):
    poly = Polygon(p)
    geo.append(poly)

crs = {'init': 'epsg:4326'}
geo_df = gpd.GeoDataFrame(crs=crs, geometry=geo)
geo_df['Acre Feet'] = df['Acre Feet'].values
geo_df.to_file('results/rutgers_acrefeet.shp', driver='ESRI Shapefile')
geo_df.to_file('results/rutgers_acrefeet.geojson', driver='GeoJSON')
geo_df.to_csv('results/rutgers_acrefeet.csv')
