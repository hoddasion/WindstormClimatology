# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 13:19:25 2022

@author: kse18nru
"""

#%% imports and globals
import xarray as xr ## to load and handle satellite nc file data
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import climatology as clima
from cartopy import config
import cartopy.crs as ccrs ## for coastlines and projections
import windrose
from windrose import WindroseAxes
import pandas as pd

matplotlib.rcParams.update({'font.size': 22})
#%%



wsp_MAM = clima.satellite_seasonal_grouping('wind_speed', 'MAM')
dir_MAM = clima.satellite_seasonal_grouping('wind_to_dir', 'MAM')
print(dir_MAM)


#%%
points_wsp = np.array(wsp_MAM.sel({'time':'2018-03-19'}, method = None))[0]
points_lat = np.array(wsp_MAM.coords['lat'])
points_lon = np.array(wsp_MAM.coords['lon'])

fig = plt.figure(figsize = (10,10))

gs = fig.add_gridspec(1,1)
blanks = np.copy(points_wsp)
blanks[:,:] = np.nan
ax0 = fig.add_subplot(gs[0,0],projection=ccrs.PlateCarree() )
ax0.pcolormesh(points_lon, points_lat, blanks, shading = 'auto')
ax0.coastlines()
ax0.set_xlim(left = -25, right = -19)
ax0.set_ylim(top = 68, bottom = 64)

sample_lons = np.array([338.35,336,335.8, 340])
sample_lats = np.array([66.6,67,65.2, 67.1])
sample_labels = [f'A: {sample_lats[0],sample_lons[0]}',
                 f'B: {sample_lats[1],sample_lons[1]}',
                 f'C: {sample_lats[2],sample_lons[2]}',
                 f'A2: {sample_lats[3],sample_lons[3]}']
ax0.scatter(sample_lons, sample_lats, s = 100, c = 'black')
for i, txt in enumerate(sample_labels):
    ax0.annotate(txt, (sample_lons[i]+0.1-360, sample_lats[i]))

#%%
dir_sample = clima.sample_in_space(dir_MAM, sample_lons, sample_lats)
wsp_sample = clima.sample_in_space(wsp_MAM, sample_lons, sample_lats)
#sample = sample.sel(lon = [340,341,342], method = 'nearest')
print(dir_sample)
print(wsp_sample)
df_A = pd.DataFrame({'speed':wsp_sample[0], 'direction':dir_sample[0]})
print(df_A)
df_A = df_A.dropna()
print(df_A)

#%%
# for year in years:
#     ## frst loa whole dataset
#     path_datafile = 'D:/Project/Climatology/Data/'
#     name_datafile = f'KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_12_ASC_20{year}.nc'
#     dataset = xr.open_dataset(f'{path_datafile}{name_datafile}')
    
#     ## extract windspeed, components, coordinates
#     points_wsp = np.array(dataset.wind_speed)
#     points_lat = np.array(dataset.lat)
#     points_lon = np.array(dataset.lon)
#     points_v = np.array(dataset.northward_wind)
#     points_u = np.array(dataset.eastward_wind)
#     points_dir = np.array(dataset.wind_to_dir)
#     print(dataset['wind_speed'])
#     #print(np.shape(points_dir))
#     #print(dataset.wind_speed.sel(time=slice(f'20{year}-03-01', f'20{year}-05-30')))
#     ## append to all lists
#     wsp_all_ds.append(dataset.wind_speed)
#     lat_all_ds.append(dataset.lat)
#     lon_all_ds.append(dataset.lon)
#     v_all_ds.append(dataset.northward_wind)
#     u_all_ds.append(dataset.eastward_wind)
#     dir_all_ds.append(dataset.wind_to_dir)
#%%
# =============================================================================
# u_all = xr.concat(u_all_ds, 'year')
# v_all = xr.concat
# print(u_all_ds[0])
# print(u_all[-1])
# =============================================================================
