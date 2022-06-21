# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:37:12 2022

@author: kse18nru
"""

#%% imports and globals
import xarray as xr ## to load and handle satellite nc file data
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import climatology as clima
from cartopy import config
import cartopy.crs as ccrs ## for coastlines and projections
import pandas as pd
from metpy.interpolate import cross_section
import datetime
import sys

matplotlib.rcParams.update({'font.size': 22})

#%%
xleft = 750
xright = 850
ybttm = 380
ytop = 480

#%%
ds_nc = xr.open_dataset('D:/Project/Climatology/Data/CARRA/CARRA_reanalysis_scalar_vars_qT_MAM_2018.nc')
ds_orog = xr.open_dataset('D:/Project/Climatology/Data/CARRA/CARRA_west_orog_and_land.nc').orog#[xleft:xright,ybttm:ytop]
ds_lsm =  xr.open_dataset('D:/Project/Climatology/Data/CARRA/CARRA_west_orog_and_land.nc').lsm[ybttm:ytop,xleft:xright]
print(ds_nc)
ds_t2m_nc = ds_nc.t2m[:,ybttm:ytop,xleft:xright]
#%% select case study
ds_t2m_case = ds_t2m_nc.sel(time = '2018-03-19')[0]
ds_t2m_mean = ds_t2m_nc.mean(dim = 'time')
ds_t2m_diff = ds_t2m_mean - ds_t2m_case 
print(ds_t2m_mean)
#%% extract arrays
x = np.array(ds_t2m_nc.coords['x'])
y = np.array(ds_t2m_nc.coords['y'])

lsm = np.array(ds_lsm)
t2m = np.array(ds_t2m_case)
t2m_noland = np.copy(t2m); t2m_noland[lsm == 1] = np.nan


#t2m_case = ds_t2m_case.data
#print(t2m_case)
#%% test plotting


fig, ax0 = plt.subplots(1,1, figsize = (12.69,10.69))

kwargs = {'shading':'auto', 'cmap':'bwr', 'vmin':-4,'vmax':4}

plt.pcolormesh(x,y,np.array(ds_t2m_diff), **kwargs)
plt.contour(x,y, np.array(ds_lsm))