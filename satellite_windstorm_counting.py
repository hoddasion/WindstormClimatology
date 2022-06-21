# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 16:24:16 2022

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
    interval = (0,90)
    da_wsp_fil = clima.filter_by_direction(da_wsp,da_dir,interval,[sample_lons[0]], [sample_lats[0]] ) # filtered DataArray
    
    #%% average the data
    da_wsp_mean = da_wsp_fil.mean(dim = 'time')
    
    #%% select 19th March 2018: Case Study
    da_wsp_case = da_wsp.sel(time='2018-03-19')[0]
    
    #%% get transect 
    da_diff = da_wsp_mean - da_wsp_case
    
    #%%
    #ptp_wsp_case = clima.count_storms.get_2point_difference(da_wsp_case)
    #print(ptp_wsp_case)
    #ptp_diff_case = clima.count_storms.get_2point_difference(da_diff)
    #print(ptp_diff_case)
   
    
    #%% extract data arrays
    wsp_case = np.array(da_wsp_case)
    diff_case = np.array(da_diff)
    
    points_lon = np.array(da_wsp_case.coords['lon'])
    points_lat = np.array(da_wsp_case.coords['lat'])
    
    #%%
    pthreshold = 2
    da_counted = clima.get_wakes(da_wsp, da_dir, -pthreshold)
    print(da_counted)
    #%% plot fields
    plot_point_locations = True
    if plot_point_locations:
        fig, (ax0,ax1) = plt.subplots(2,1,figsize = (10,12))#,gridspec_kw={'height_ratios': [3, 1]})
        pcol_kwargs = {'shading':'auto', 'cmap':'Oranges', 'vmin':0,'vmax':20}
        diff_kwargs = {'shading':'auto', 'cmap':'bwr', 'vmin':-6,'vmax':6}
        
        ## wsp field
        p0 = ax0.pcolormesh(points_lon, points_lat, wsp_case, **pcol_kwargs)
        ax0.set_xlim(left = minlon, right = maxlon)
        ax0.set_ylim(bottom = minlat, top = maxlat)
        ax0.set_title('windspeed')
        point1 = (66.6,338.5); point2 = (67.3,338)
        ax0.scatter((point1[1], point2[1]),(point1[0], point2[0]), color = 'k')
        ax0.scatter((338.5625, 338.0625), (66.5625, 67.3125), color = 'green')
        ax0.set_xticks([])
        ax0.set_yticks([])
        
        for i, point in enumerate([point1, point2]):
            ax0.annotate(i+1, (point[1]+0.1, point[0]))
        
        cbar1 = fig.colorbar(mappable=p0,ax = ax0, orientation = 'vertical', fraction = 0.15, pad = 0.04 )
        cbar1.ax.set(ylabel = r'ms$^{-1}$')
        
            
        ## difference field
        p1 = ax1.pcolormesh(points_lon, points_lat, diff_case, **diff_kwargs)
        ax1.set_xlim(left = minlon, right = maxlon)
        ax1.set_ylim(bottom = minlat, top = maxlat)
        ax1.set_title('mean field anomaly: mean - case')
        
        ax1.scatter((point1[1], point2[1]),(point1[0], point2[0]), color = 'k')
        ax1.scatter((338.5625, 338.0625), (66.5625, 67.3125), color = 'green')
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        for i, point in enumerate([point1, point2]):
            ax1.annotate(i+1, (point[1]+0.1, point[0]))
        
        fig.suptitle('Case 2 - 2018-03-19 - field comparison')
        
        cbar2 = fig.colorbar(mappable=p1,ax = ax1, orientation = 'vertical', fraction = 0.15, pad = 0.04 )
        cbar2.ax.set(ylabel = r'ms$^{-1}$')
        
    #%%
    plot_counted_fields = True
    if plot_counted_fields:
        da_wsp = clima.satellite_seasonal_grouping('wind_speed', season)
        da_u = clima.satellite_seasonal_grouping('eastward_wind', season)
        da_v = clima.satellite_seasonal_grouping('northward_wind', season)
        da_dir = clima.satellite_seasonal_grouping('wind_to_dir', season)
        interval = (0,90)
        
        pthreshold = 6
        da_wsp_counted = clima.get_wakes(da_wsp, da_dir, -pthreshold)
        
        
        ptp_time = np.array(da_wsp_counted.coords['time'])
        da_u_counted = da_u.sel(time = ptp_time)
        da_v_counted = da_v.sel(time = ptp_time)
        
        
        print(da_u_counted)
        
        #%%
        mass_plots = True
        if mass_plots:
            for i,time in enumerate(np.array(da_wsp_counted.coords['time'])):
                print(time)
                da_wsp_day = da_wsp_counted[i]
                da_u_day = da_u_counted[i]
                da_v_day = da_v_counted[i]
                quiv = 5
                
                points_wsp = np.array(da_wsp_day)
                points_u = np.array(da_u_day)
                points_v = np.array(da_v_day)
                points_lon = np.array(da_wsp_day.coords['lon'])
                points_lat = np.array(da_wsp_day.coords['lat'])
                
                
                fig, ax0 = plt.subplots(1,1, figsize = (14,8),subplot_kw=dict(projection=ccrs.PlateCarree()))
                pcol_kwargs = {'shading':'auto', 'cmap':'Oranges'}#, 'vmin':0,'vmax':20}
                p0 = ax0.pcolormesh(points_lon, points_lat, points_wsp, **pcol_kwargs)
                ax0.quiver(points_lon[::quiv], points_lat[::quiv], points_u[::quiv,::quiv], points_v[::quiv,::quiv])
                ax0.set_xticks([])
                ax0.set_yticks([])
                ax0.coastlines()
                cbar1 = fig.colorbar(mappable=p0,ax = ax0, orientation = 'vertical', fraction = 0.025, pad = 0.04 )
                cbar1.ax.set(ylabel = r'ms$^{-1}$')
                
                ax0.set_title(f'{str(time)[:10]} 10m windspeed')
                plt.tight_layout()
                plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/Counting/ptp_method/wsp/threshold_m{pthreshold}/ASCATA_10m_windspeed_{str(time)[:10]}.png')
                plt.close()
                
        #%%
        plot_histogram = False
        if plot_histogram:
           
            date_nparray = np.array(da_wsp_counted.coords['time'])
            print(date_nparray)
            months = np.array([])
            years = np.array([])
            for i, date in enumerate(date_nparray):
                years = np.concatenate((years, np.array([date.year])))
                months = np.concatenate((months,np.array([date.month])))
            print(months)
            print(np.histogram(years))
            ## plot
            hist_years = np.histogram(years)
            fig,ax0  = plt.subplots(1,1, figsize = (10,10))
            
            ax0.plot(hist_years)