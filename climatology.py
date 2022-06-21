# -*- coding: utf-8 -*-
"""
Created on Tue May  3 16:36:04 2022

@author: kse18nru
"""

import iris
import numpy as np
import sys
import math
import xarray as xr
import count_storms

def subset_cube_by_coord(cube,vertexSW, vertexNE, lonname = 'grid_longitude', latname = 'grid_latitude', surfdiag = False, dimtime = True):
    
    ## first extract horizontal coordiate arrays
    try:
        gridlon = cube.coord(lonname).points
        gridlat = cube.coord(latname).points
        condlon = (gridlon < vertexNE[0]) * (gridlon > vertexSW[0])
        condlat = (gridlat < vertexNE[1]) * (gridlat > vertexSW[1])
    except:
        gridlon = cube.coord(lonname)[0].points
        gridlat = cube.coord(latname)[0].points
        
        condlon = (gridlon < vertexNE[0]) * (gridlon > vertexSW[0])
        condlat = (gridlat < vertexNE[1]) * (gridlat > vertexSW[1])
    
    # if surfdiag == True:
    #     if dimtime == False:
    #         subcube = cube[]
    
    dimnum = len(np.shape(cube.data))
    if dimnum == 2:
        subcube = cube[condlat,:]
        subcube = subcube[:,condlon]
    elif dimnum == 3:
        subcube = cube[:,condlat,:]
        subcube = subcube[:,:,condlon]
    elif dimnum == 4:
        subcube = cube[:,:,condlat,:]
        subcube = subcube[:,:,:,condlon]
    else:
        print('Too many or too few dimensions in cube')
        sys.exit()
    
    return subcube

def quadmag(X,Y):
    ## element wise magntiude by quadrature
    return (X**2 + Y**2)**0.5

def haversine(lons, lats):
    """
    Calculate the great circle distance in meters between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    ## first get horizontal step in distance grid
    lons = lons*(np.pi/180) # convert to radians
    lats = lats*(np.pi/180)

    # haversine formula 
    dlon = lons[1:] - lons[:-1]       #lon2 - lon1 
    dlat = lats[1:] - lats[:-1]                 #lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lats[:-1]) * np.cos(lats[1:]) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return (c * r)*1000

def divergence_horizontal(u,v,lon,lat, r = 6371e3):
    #print('u', np.shape(u), 'lon', np.shape(lon), 'v', np.shape(v), 'lat', np.shape(lat))
    ## first get horizontal step in distance grid
    rlon = lon*(np.pi/180) # convert to radians
    rlat = lat*(np.pi/180)
    
    # take angle difference; change in position from cell to cell vertex
    drlon = rlon[:,1:] - rlon[:,:-1]
    drlat = rlat[1:] - rlat[:-1]
    drlon, drlat
    
    dx = r*np.sin(drlon)*np.cos(rlat[:,:-1])
    dy = r*np.sin(drlat)
    dx = dx[:-1]
    dy = dy[:,:-1]
    ## get cell to cell velcity changes
    
    du = u[:,1:]-u[:,:-1]
    dv = v[1:]-v[:-1]
    du = du[:-1]
    dv = dv[:,:-1]
    
    ## compute divergence in 2D
    #print('du', np.shape(du), 'dx', np.shape(dx), 'dv', np.shape(dv), 'dy', np.shape(dy))
    div = du/dx + dv/dy
    #print('div',div)
    return div

def altitude_pressure_adjust(msp, alt):
    p = msp
    return p

def rh_to_q_calculate(RH, T, p, Tref = 273.15):
    """
    Function to convert relative humidity to specfic humidity. Numpy array compatible.

    Parameters
    ----------
    RH : Relative humidity in %.
    T : Air temperature in Kelvin.
    p : Air pressure in Pa.
    

    Returns
    -------
    q : specific humidity

    """
    es = 611.2*np.exp(17.67*(T-Tref)/(T-29.65)) # partial vapor pressure
    rvs = 0.622*es/(p - es) # saturation water vapor mixing ratio
    rv = (RH/100) * rvs # water vapor mixing ratio
    q = rv/(1 + rv) # specific humditity
    
    return q
    
def regrid_cubes(cube1, cube2, scheme, guess_bounds = False):
    if guess_bounds:
        try: # guess bounds for area weighting
            cube1.coord(axis = 'x').guess_bounds()
            cube1.coord(axis = 'y').guess_bounds()
            cube2.coord(axis = 'x').guess_bounds()
            cube2.coord(axis = 'y').guess_bounds()
            
            return cube1.regrid(cube2, scheme)
        except: # if bounds are already guessed once, error is thrown
            return cube1.regrid(cube2, scheme)
    else:
        return cube1.regrid(cube2, scheme)
    
    
def regrid_numpy(data1,lon1,lat1, lon2,lat2,  return_coords = False):
    
    
    xlen = len(lon2)
    ylen = len(lat2)
    histo = np.histogram2d(lon1,lat1,bins = [xlen,ylen], weights = data1)
    
    if return_coords:
        return histo, lon2, lat2
    else:
        return histo
    
def convert_glm_coordsystem(glmcube, regcube):
    polelat = regcube.coord(axis = 'y').coord_system.grid_north_pole_latitude
    polelon = regcube.coord(axis = 'x').coord_system.grid_north_pole_longitude
    # extract regional model coord system
    regcoordsys = regcube.coord_system(iris.coord_systems.CoordSystem)
    glmcoordsys = glmcube.coord_system(iris.coord_systems.CoordSystem)
    print(regcoordsys)
    print(glmcube.coord(axis = 'x'))
    
def regrid_glm_and_compute_q(regcube, glm_RH, glm_P, glm_T, scheme, Tref = 273.15):
    
    ## regrid cubes
    regrid_glm_RH = regrid_cubes(glm_RH, regcube, scheme)
    regrid_glm_P = regrid_cubes(glm_P, regcube, scheme)
    regrid_glm_T = regrid_cubes(glm_T, regcube, scheme)
    
    ## compute q
    glm_q = rh_to_q_calculate(regrid_glm_RH.data, regrid_glm_T.data, regrid_glm_P.data, Tref = Tref)
    return glm_q
        

def satellite_seasonal_grouping(variable_name, season, path_datafile = 'D:/Project/Climatology/Data/',
                                years = [10,11,12,13,14,15,16,17,18,19,20],
                                name_datafile = 'KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_12_ASC_20',
                                concatenate_data = True):
    DataArray_list = []
    for year in years:
        yearfilename = f'{name_datafile}{year}.nc'
        
        dataset = xr.open_dataset(f'{path_datafile}{yearfilename}')
        DataArray = dataset[variable_name]
        
        if season == 'MAM':
            season_DataArray = DataArray.sel(time=slice(f'20{year}-03-01', f'20{year}-05-31'))
        elif season == 'JJA':
            season_DataArray = DataArray.sel(time=slice(f'20{year}-06-01', f'20{year}-08-31'))
        elif season == 'SON':
            season_DataArray = DataArray.sel(time=slice(f'20{year}-09-01', f'20{year}-11-30'))
        elif season == 'DJF':
            try:
                JF_DataArray = DataArray.sel(time=slice(f'20{year}-01-01', f'20{year}-02-28'))
                yearfilename = f'{name_datafile}{year-1}.nc'
                dataset = xr.open_dataset(f'{path_datafile}{yearfilename}')
                DataArray = dataset[variable_name]
                Dec_DataArray = DataArray.sel(time=slice(f'20{year - 1}-12-01', f'20{year - 1}-12-31'))
                season_DataArray = xr.concat([Dec_DataArray, JF_DataArray], dim = 'time')
            except: 
                print(f'year 20{year} skipped')
                continue
        else:
            print('Incorrect season definition given, or specfific months not provided for. \nCompatible codes are: DJF, MAM, JJA, SON')
        DataArray_list.append(season_DataArray)
        
    if concatenate_data:
        ## concatenate DataArrays into single DataArray along axis "time"
        conc_data = xr.concat(DataArray_list, dim = 'time')
        return conc_data
    else:
        ## return non-concatenated list instead
        return DataArray_list

def sample_in_space(DataArray, sample_lons, sample_lats):
    indexers = {'lon':sample_lons, 'lat':sample_lats}
    da_sample = DataArray.sel(indexers, method = 'nearest')
    array_sample = []
    try: np.array(sample_lons) - np.array(sample_lats)
    except: 
        print('sample latitude and longitude arrays passed to climatology.sample_in_space() must be arrays or lists of floats of equal length.')
        sys.exit()
    for i in range(len(sample_lons)):
        timeseries = np.array(da_sample[:,i,i])
        array_sample.append(timeseries)
    return np.array(array_sample)

def filter_by_direction(DataArray, wind_dir, test_range,test_lons, test_lats):
    indexers = {'lon':test_lons, 'lat':test_lats}
    sample_dir = sample_in_space(wind_dir, test_lons, test_lats)
    condition = True
    for i in range(len(sample_dir)):
        condition_upper = sample_dir[i] >= test_range[0]
        condition_lower = sample_dir[i] <= test_range[1]
        condition_point = condition_upper*condition_lower
        condition = condition*condition_point
    da_filtered = DataArray[condition]
    return da_filtered

def get_wakes(DataArray, da_dir, threshold, interval = (0,90),method = 'ptp_method',
              point1 = (66.6,338.35), point2 = (67.3,338), point3 = (67.1,340),condition_type = 'smaller'):
    
    ## calc mean wrt time
    da_mean = DataArray.mean(dim = 'time')
    ## filter samples
    test_lons = np.array([point1[1], point2[1]])
    test_lats = np.array([point1[0], point2[0]]) # for my purposes I will just be using points1 and 2 anyway
    # filter by direction
    da_filtered = filter_by_direction(DataArray, da_dir, interval, test_lons,test_lats )
    print('Number of filtered days =', len(np.array(da_filtered.coords['time'])))
    ## count wakes and return dataarray
    if method == 'ptp_method':
        # select method to count by; new ones may be added in future
        da_counted = count_storms.two_point_difference_method(da_filtered, da_mean, threshold = threshold, condition_type = condition_type,
                                                              point1 = point1, point2 = point2, point3 = point3)
    print('Number of counted wakes =', len(np.array(da_counted.coords['time'])))
    
    return da_counted