# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 12:50:04 2022

@author: kse18nru
"""
#%%
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from cartopy import config
import cartopy.crs as ccrs
#%%
path_datafile = 'D:/Project/Climatology/Data/'
name_datafile = 'KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_12_ASC_2013.nc'

dataset = xr.open_dataset(f'{path_datafile}{name_datafile}')


matplotlib.rcParams.update({'font.size': 22})
 
#%%
print(dataset)
#

#%%

ds_wsp = dataset.wind_speed
#ds_wdir = dataset.wind_dir
ds_lat = dataset.lat
ds_lon = dataset.lon
print(ds_wsp)

#%%
points_wsp = np.array(ds_wsp)
points_lat = np.array(ds_lat)
points_lon = np.array(ds_lon)
points_v = np.array(dataset.northward_wind)
points_u = np.array(dataset.eastward_wind)
print(points_wsp[18])




#%%
tidx = 18
qstep = 6
fig = plt.figure(figsize = (18,18))

gs = fig.add_gridspec(1,1)

ax0 = fig.add_subplot(gs[0,0],projection=ccrs.PlateCarree() )

#ax0.set_ylim(bottom = 55.1, top = 75.8)
ax0.pcolormesh(points_lon, points_lat,points_wsp[tidx])
ax0.quiver(points_lon[::qstep], points_lat[::qstep], points_u[tidx,::qstep,::qstep], points_v[tidx,::qstep,::qstep],headlength = 4, headwidth = 2)
ax0.coastlines()
ax0.set_xlim(left = -30.9, right = -12)