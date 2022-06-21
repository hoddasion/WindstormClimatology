# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 16:32:54 2022

@author: kse18nru
"""
import iris
import numpy as np
import sys
import math
import xarray as xr


def get_2point_difference(DataArray, point1 = (66.6,338.35), point2 = (67.3,338)):
    
    datapoint1 = DataArray.sel({'lon':point1[1],'lat':point1[0]}, method = 'Nearest')
    datapoint2 = DataArray.sel({'lon':point2[1],'lat':point2[0]}, method = 'Nearest')
    
    #print(datapoint1)
    print('Datapoint 2 =',datapoint2)
    ptp_diff = np.array(datapoint1 - datapoint2)
    
    return ptp_diff

def two_point_difference_method(DataArray, da_mean, threshold, point1 = (66.6,338.35), point2 = (67.3,338), point3 = (67.1,340),condition_type = 'larger'):
    #da_diff = da_mean - da_wsp
    datapoint3 = DataArray.sel({'lon':point3[1],'lat':point3[0]}, method = 'Nearest')
    #print(datapoint3 )
    ptp_diff = get_2point_difference(DataArray, point1 = point1, point2 = point2)
    
    print('Number of input days =', len(ptp_diff))
    if condition_type == 'larger':
        condition = ptp_diff > threshold
    if condition_type == 'smaller':
        condition = ptp_diff < threshold
    
    ## point 3 data availability test
    # check for nans, then invert boolean array; finally combine with ptp condition
    condition = condition * np.invert(np.isnan(np.array(datapoint3))) 
    
    return DataArray[condition] 

def curve_method(da_wsp, threshold, ):
    da_wsp_mean = da_wsp.mean(dim = 'time')
    time_dim = np.array(da_wsp.coords['time'])
    for day in time_dim:
        wsp_case = np.array(da_wsp.sel(time=day)[0])
        wsp_mean =  np.array(da_wsp_mean)
        wsp_diff = wsp_mean - wsp_case
    return 0