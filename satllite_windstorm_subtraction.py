# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 17:47:53 2022

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

matplotlib.rcParams.update({'font.size': 22})

#%%
sample_lons = np.array([338.35,336,335.8, 340])
sample_lats = np.array([66.6,67,65.2, 67.1])
sample_labels = ['A1', 'B', 'C']#, 'A2']

#%%
seasons = ['MAM']#,'JJA','SON','DJF']
minlon, maxlon, minlat, maxlat = (-25+360, -19+360, 64, 69)
for season in seasons:
    ##load and subset wsp and wdir data by season
    da_wsp = clima.satellite_seasonal_grouping('wind_speed', season)#.sel({'lon':slice(minlon,maxlon),'lat':slice(minlat,maxlat)})#, method = 'Nearest')
    da_dir = clima.satellite_seasonal_grouping('wind_to_dir', season)#.sel({'lon':slice(minlon,maxlon),'lat':slice(minlat,maxlat)})#, method = 'Nearest')
    interval = (0,90)#(180,270)
    da_wsp_fil = clima.filter_by_direction(da_wsp,da_dir,interval,[sample_lons[0]], [sample_lats[0]] ) # filtered DataArray
    print(da_wsp)
    print(da_wsp_fil)
    #%% average the data
    da_wsp_mean = da_wsp_fil.mean(dim = 'time')
    #print(da_wsp_mean)
    
    #%% select 19th March 2018: Case Study
    da_wsp_case = da_wsp.sel(time='2018-03-19')[0]
    #print(da_wsp_case)

    #%% extract numpy arrays and take difference
    wsp_mean = np.array(da_wsp_mean)
    points_lat = np.array(da_wsp_mean.coords['lat'])
    points_lon = np.array(da_wsp_mean.coords['lon'])
    
    wsp_case = np.array(da_wsp_case)
    
    wsp_diff = wsp_mean - wsp_case
    
    #%%
    wspmin = 0; wspmax = 20
    
    print(np.nanmax(wsp_mean))
    print(np.nanmax(wsp_case))
    print(np.nanmin(wsp_diff),np.nanmax(wsp_diff))
    #%% plot fields
    plot_threepanels = False
    if plot_threepanels:
    
        fig, (ax0,ax1,ax2) = plt.subplots(1,3,figsize = (20,10))
        pcol_kwargs = {'shading':'auto', 'cmap':'Oranges', 'vmin':0,'vmax':wspmax}
        
        
        ax0.pcolormesh(points_lon, points_lat, wsp_mean, **pcol_kwargs)
        ax0.set_title('2010-2020 MAM Average ')
        
        ax1.pcolormesh(points_lon,points_lat, wsp_case, **pcol_kwargs)
        ax1.set_title('Case 2018-03-19')
        
        ax2.pcolormesh(points_lon, points_lat, wsp_diff, shading = 'auto', cmap = 'bwr', vmin = -6, vmax = 6)
        ax2.set_title('Seasonal avr. -  case field')
        
        wsp_norm = matplotlib.colors.Normalize(vmin =wspmin, vmax = wspmax)
        wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
        cbar1 = fig.colorbar(mappable=wsp_map,ax = (ax0), orientation = 'horizontal', fraction = 0.15, pad = 0.01, aspect = 10, ticks = [0,5,10,15,20] )
        cbar1.ax.set(xlabel = r'ms$^{-1}$')
        
        cbar12 = fig.colorbar(mappable=wsp_map,ax = (ax1), orientation = 'horizontal', fraction = 0.15, pad = 0.01, aspect = 10, ticks = [5,10,15,20] )
        cbar12.ax.set(xlabel = r'ms$^{-1}$')
        
        diff_norm = matplotlib.colors.Normalize(vmin =-6, vmax = 6)
        diff_map = matplotlib.cm.ScalarMappable(norm = diff_norm, cmap = 'bwr')
        cbar2 = fig.colorbar(mappable=diff_map,ax = (ax2), orientation = 'horizontal', fraction = 0.15, pad = 0.01, aspect = 10 )
        cbar2.ax.set(xlabel = r'ms$^{-1}$')
        
        for ax in (ax0,ax1,ax2):
            ax.set_xticks([])
            ax.set_yticks([])
        
        plt.subplots_adjust(left=0.1,
                            right=0.9, 
                            top = 0.87,
                            bottom = 0, 
                            wspace=0.05, 
                            hspace=0.35)
        
        fig.suptitle(f'METOP-A ASCAT-A Conditional Sampling: {interval}degrees; 10m windspeed')
        
        plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/Satellite/ASCATA_10m_wsp_MAMavr_Case2field_difference_plot_wholeregion_{interval[0]}_{interval[1]}.png')
        
        
    #%%
    plot_case = True
    if plot_case:
        ## load additional variables
        da_u = clima.satellite_seasonal_grouping('eastward_wind', season)
        da_v = clima.satellite_seasonal_grouping('northward_wind', season)
        da_u_case = da_u.sel(time='2018-03-19')[0]
        da_v_case = da_v.sel(time='2018-03-19')[0]
        
        A1_dir = da_dir.sel({'lon':sample_lons[0], 'lat':sample_lats[0], 'time':'2018-03-19'}, method = 'Nearest')
        print(A1_dir)
        
        #%% plot
        fig, ax0 = plt.subplots(1,1, figsize = (7,7))
        
        ax0.quiver(points_lon[::5], points_lat[::5],np.array(da_u_case)[::5,::5], np.array(da_v_case)[::5,::5])
        ax0.set_title('ASCAT-A wind direction for Case 2 from wind components', fontsize = 13)
        ax0.scatter(sample_lons[0], sample_lats[0])
        ax0.annotate('A1',(sample_lons[0]-0.9,sample_lats[0]-0.4))
        