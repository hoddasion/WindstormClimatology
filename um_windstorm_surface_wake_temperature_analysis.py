# -*- coding: utf-8 -*-
"""
Created on Fri May  6 15:03:20 2022

@author: kse18nru
"""

#%%
import climatology
import iris
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.analysis.cartography import rotate_pole
import warnings
import xarray as xr
from cartopy import config
import cartopy.crs as ccrs


#%% globals
um_suite = 'u-cc134'
um_datapath = f'D:/Project/Model_Data/{um_suite}/'
resolutions = ['0p5km','1p5km','4p4km']
matplotlib.rcParams.update({'font.size': 22})
#warnings.filterwarnings("ignore", message="warning")

#%%
wake_hoz_xsecs = True
if wake_hoz_xsecs:
    for chunk in [4]:
        
        ##
        if chunk == 4:
            titletime = '12.5h-15h'
            glmtitle = '1500hrs'
            glmidx = 0
            glm_file = 3
            maxlon_0p5km = [360.21249]
            maxlat_0p5km = [0.49]
            maxlon_1p5km = [360.28753]
            maxlat_1p5km = [0.37]
            maxlon_4p4km = [360.26]
            maxlat_4p4km = [0.44]
        elif chunk == 5:
            titletime = '15.5h-18h'
            glmtitle = '1800hrs'
            glmidx = 1
            glm_file = 3
            maxlon_0p5km = [360.21249]
            maxlat_0p5km = [0.49]
            maxlon_1p5km = [360.25753]
            maxlat_1p5km = [0.493]
            maxlon_4p4km = [360.26]
            maxlat_4p4km = [0.48]
        elif chunk == 6:
            titletime = '18.5h-21h'
            glmtitle = '2100hrs'
            glmidx = 0
            glm_file = 4
            maxlon_0p5km = [360.21249]
            maxlat_0p5km = [0.49]
            maxlon_1p5km = [360.21253]
            maxlat_1p5km = [0.49]
            maxlon_4p4km = [360.30002]
            maxlat_4p4km = [0.44]
        ## also perform 3 hour averages
        #chunk = 5
        ss = 6 # ss : sample size; for half hourly outputs, a sample of 6 represents 3 hours
        ## all of these below are surface diagnostics
        temp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        temp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        temp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        temp_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc','air_temperature')[glmidx]
        
        ## load binary land mask
        lbm_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_1_flt306.nc', 'land_binary_mask')
        
        ## load msp
        msp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_glm = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'air_pressure_at_sea_level')[glmidx]
        
        print(temp_0p5km)
        #%% subsetting to region
        glmNE = (341,68)
        glmSW = (335,65)
        
        polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
        rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        
        regNE = (rot_db_lon[1], rot_db_lat[1])
        regSW = (rot_db_lon[0], rot_db_lat[0])  
        
        print(regSW, regNE)
        
        ## subset wind
        sub_temp_0p5km = climatology.subset_cube_by_coord(temp_0p5km,regSW, regNE)
        sub_temp_1p5km = climatology.subset_cube_by_coord(temp_1p5km,regSW, regNE)
        sub_temp_4p4km = climatology.subset_cube_by_coord(temp_4p4km,regSW, regNE)
        sub_temp_glm   = climatology.subset_cube_by_coord(temp_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        ## subset land-binary-mask
        sub_lbm_0p5km = climatology.subset_cube_by_coord(lbm_0p5km,regSW, regNE)
        sub_lbm_1p5km = climatology.subset_cube_by_coord(lbm_1p5km,regSW, regNE)
        sub_lbm_4p4km = climatology.subset_cube_by_coord(lbm_4p4km,regSW, regNE)
        sub_lbm_glm   = climatology.subset_cube_by_coord(lbm_glm, glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        ## subset msp
        sub_msp_0p5km = climatology.subset_cube_by_coord(msp_0p5km,regSW, regNE)
        sub_msp_1p5km = climatology.subset_cube_by_coord(msp_1p5km,regSW, regNE)
        sub_msp_4p4km = climatology.subset_cube_by_coord(msp_4p4km,regSW, regNE)
        sub_msp_glm = climatology.subset_cube_by_coord(msp_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        lon_0p5km = sub_temp_0p5km.coord('grid_longitude').points
        lat_0p5km = sub_temp_0p5km.coord('grid_latitude').points
        
        lon_1p5km = sub_temp_1p5km.coord('grid_longitude').points
        lat_1p5km = sub_temp_1p5km.coord('grid_latitude').points
        
        lon_4p4km = sub_temp_4p4km.coord('grid_longitude').points
        lat_4p4km = sub_temp_4p4km.coord('grid_latitude').points
        
        lon_glm = sub_temp_glm.coord('longitude').points
        lat_glm = sub_temp_glm.coord('latitude').points
        
        #%%
        #print(sub_msp_0p5km.data)
        msp_levels = np.arange(90000,130000, 100)
        #print(msp_levels)
        
        #%% check bounds
        print(np.nanmin(sub_temp_0p5km.data), np.nanmax(sub_temp_0p5km.data))
        
        #%% cursory plotting
        
        
        plot_wholeregion = False
        if plot_wholeregion:
        
            fig = plt.figure(figsize = (12,12))
            gs = fig.add_gridspec(2,2)
            kwargs = {'cmap' : 'plasma', 'vmin' : 270, 'vmax' : 281, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            
            fig.suptitle(f'{titletime} avr. 1.5m temperature - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.pcolormesh(lon_0p5km, lat_0p5km, sub_temp_0p5km.data, **kwargs)
            ax0.contour(lon_0p5km, lat_0p5km, sub_lbm_0p5km.data, colors = 'green', levels = [1], linewidths = 4)
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim( right = regNE[0])
            ax0.set_ylim( top = regNE[1])
            
            ax1 = fig.add_subplot(gs[0,1])
            ax1.pcolormesh(lon_1p5km, lat_1p5km, sub_temp_1p5km.data, **kwargs)
            ax1.contour(lon_1p5km, lat_1p5km, sub_lbm_1p5km.data, colors = 'green', levels = [1], linewidths = 4)
            ax1.contour(lon_1p5km, lat_1p5km, sub_msp_1p5km.data, **mspkwargs)
            ax1.set_title('1p5km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            
            ax2 = fig.add_subplot(gs[1,0])
            ax2.pcolormesh(lon_4p4km, lat_4p4km, sub_temp_4p4km.data, **kwargs)
            ax2.contour(lon_4p4km, lat_4p4km, sub_lbm_4p4km.data, colors = 'green', levels = [1], linewidths = 4)
            ax2.contour(lon_4p4km, lat_4p4km, sub_msp_4p4km.data, **mspkwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            
            ax3 = fig.add_subplot(gs[1,1])#,projection=ccrs.PlateCarree())
            ax3.pcolormesh(lon_glm, lat_glm, sub_temp_glm.data[glmidx], **kwargs)
            ax3.contour(lon_glm, lat_glm, sub_msp_glm.data[glmidx], **mspkwargs)
            ax3.contour(sub_lbm_glm.coord('longitude').points, sub_lbm_glm.coord('latitude').points, sub_lbm_glm.data, colors = 'green', levels = [1], linewidths = 4)
            ax3.set_title(f'GLM - {glmtitle}')
            ax3.set_xticks([])
            ax3.set_yticks([])
            
            norm = matplotlib.colors.Normalize(vmin =270, vmax = 281)
            Tmap = matplotlib.cm.ScalarMappable(norm = norm, cmap = 'plasma')
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.87])
            cbar = fig.colorbar(mappable=Tmap,cax = cbar_ax)
            cbar.ax.set(ylabel = 'K') 
            
            plt.tight_layout()
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_air_temperature_3houravr_upper{glmtitle}.png')
            
        
        #%% remove land points
        ## nl means "no-land"
        nl_temp_0p5km = sub_temp_0p5km.data.copy()
        nl_temp_1p5km = sub_temp_1p5km.data.copy()
        nl_temp_4p4km = sub_temp_4p4km.data.copy()
        nl_temp_glm   =sub_temp_glm.data.copy()
        
        nl_temp_0p5km[sub_lbm_0p5km.data == 1] = np.nan
        nl_temp_1p5km[sub_lbm_1p5km.data == 1] = np.nan
        nl_temp_4p4km[sub_lbm_4p4km.data == 1] = np.nan
        nl_temp_glm[sub_lbm_glm.data == 1] = np.nan
        
        #%% check scale
        print(np.nanmax(nl_temp_0p5km), np.nanmin(nl_temp_0p5km))
        
        #%% set new subsetting coordinates
        ## using "p" to denote plotting purpose
        plot_wake0p5kmonly = False
        if plot_wake0p5kmonly:
            ## plot new
            
            fig, ax0 = plt.subplots(1,1,figsize = (8,10))
            
            kwargs = {'cmap' : 'plasma', 'vmin' : 275, 'vmax' : 279, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            
            fig.suptitle(f'{titletime} avr. 1.5m temperature\n2018-03-19')
            
            
            ax0.pcolormesh(lon_0p5km, lat_0p5km, nl_temp_0p5km.data, **kwargs)
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            
            div_norm = matplotlib.colors.Normalize(vmin = 275, vmax = 279)
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'plasma')
            cbar = fig.colorbar(mappable=div_map,ax = ax0, orientation = 'horizontal', fraction = 0.15, pad = 0.01 )
            cbar.ax.set(xlabel = 'K')
            plt.subplots_adjust(left=0.125,
                    bottom=0.25, 
                    right=0.9, 
                    top=0.85, 
                    wspace=0.05, 
                    hspace=0.35)
            
            plt.tight_layout()
            # ax3.pcolormesh(lon_glm, lat_glm, nl_temp_glm.data, **kwargs)
            # ax3.contour(lon_glm, lat_glm, sub_msp_glm.data, **mspkwargs)
            # ax3.set_title(f'GLM - {glmtitle}')
            # ax3.set_xticks([])
            # ax3.set_yticks([])
            
        #%%
        plot_wakeregion = False
        if plot_wakeregion:
            pglmNE = (340,67.4)
            pglmSW = (337,66)
            
            polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            pregNE = (rot_db_lon[1], rot_db_lat[1])
            pregSW = (rot_db_lon[0], rot_db_lat[0])  
            
            ## plot new
            
            fig = plt.figure(figsize = (12,12))
            gs = fig.add_gridspec(2,2)
            kwargs = {'cmap' : 'plasma', 'vmin' : 275, 'vmax' : 279, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            
            fig.suptitle(f'{titletime} avr. 1.5m temperature - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.pcolormesh(lon_0p5km, lat_0p5km, nl_temp_0p5km.data, **kwargs)
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.scatter([maxlon_0p5km], [maxlat_0p5km], s = 600, marker = 'x', linewidths = 5, color = 'blue')
            
            ax1 = fig.add_subplot(gs[0,1])
            ax1.pcolormesh(lon_1p5km, lat_1p5km, nl_temp_1p5km.data, **kwargs)
            ax1.contour(lon_1p5km, lat_1p5km, sub_msp_1p5km.data, **mspkwargs)
            ax1.set_title('1p5km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.scatter([maxlon_1p5km], [maxlat_1p5km], s = 600,marker = 'x', linewidths = 5, color = 'blue')
            
            ax2 = fig.add_subplot(gs[1,0])
            ax2.pcolormesh(lon_4p4km, lat_4p4km, nl_temp_4p4km.data, **kwargs)
            ax2.contour(lon_4p4km, lat_4p4km, sub_msp_4p4km.data, **mspkwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.scatter([maxlon_4p4km], [maxlat_4p4km], s = 600,marker = 'x', linewidths = 5, color = 'blue')
            
            ax3 = fig.add_subplot(gs[1,1])#,projection=ccrs.PlateCarree())
            ax3.pcolormesh(lon_glm, lat_glm, nl_temp_glm.data, **kwargs)
            ax3.contour(lon_glm, lat_glm, sub_msp_glm.data[glmidx], **mspkwargs)
            ax3.set_title(f'GLM - {glmtitle}')
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            
            norm = matplotlib.colors.Normalize(vmin =275, vmax = 279)
            Tmap = matplotlib.cm.ScalarMappable(norm = norm, cmap = 'plasma')
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.87])
            cbar = fig.colorbar(mappable=Tmap,cax = cbar_ax)
            cbar.ax.set(ylabel = 'K') 
            
            plt.tight_layout()
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_air_temperature_3houravr_upper{glmtitle}_noland.png')
            
        
        #%% section for surface altitude plots
        surf_alt_plot = False
        if surf_alt_plot:
            
            alt_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_surface_altitude_0_24hrs_pi_306.nc', 'surface_altitude_0')
            alt_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_surface_altitude_0_24hrs_pi_306.nc', 'surface_altitude_0')
            alt_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_surface_altitude_0_24hrs_pi_306.nc', 'surface_altitude_0')
            alt_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_1_flt306.nc', 'surface_altitude')
            
            lbm_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
            lbm_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
            lbm_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
            lbm_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_1_flt306.nc', 'land_binary_mask')
            #%% subsetting to region
            glmNE = (341,68)
            glmSW = (335,65)
            
            polelat = alt_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = alt_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            regNE = (rot_db_lon[1], rot_db_lat[1])
            regSW = (rot_db_lon[0], rot_db_lat[0])  
            
            print(regSW, regNE)
            ## subset land-binary-mask
            sub_lbm_0p5km = climatology.subset_cube_by_coord(lbm_0p5km,regSW, regNE)
            sub_lbm_1p5km = climatology.subset_cube_by_coord(lbm_1p5km,regSW, regNE)
            sub_lbm_4p4km = climatology.subset_cube_by_coord(lbm_4p4km,regSW, regNE)
            sub_lbm_glm   = climatology.subset_cube_by_coord(lbm_glm, glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
            
            ## subset orography
            sub_alt_0p5km = climatology.subset_cube_by_coord(alt_0p5km,regSW, regNE)
            sub_alt_1p5km = climatology.subset_cube_by_coord(alt_1p5km,regSW, regNE)
            sub_alt_4p4km = climatology.subset_cube_by_coord(alt_4p4km,regSW, regNE)
            sub_alt_glm   = climatology.subset_cube_by_coord(alt_glm, glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
            
            lon_0p5km = sub_alt_0p5km.coord('grid_longitude').points
            lat_0p5km = sub_alt_0p5km.coord('grid_latitude').points
            
            lon_1p5km = sub_alt_1p5km.coord('grid_longitude').points
            lat_1p5km = sub_alt_1p5km.coord('grid_latitude').points
            
            lon_4p4km = sub_alt_4p4km.coord('grid_longitude').points
            lat_4p4km = sub_alt_4p4km.coord('grid_latitude').points
            
            lon_glm = sub_alt_glm.coord('longitude').points
            lat_glm = sub_alt_glm.coord('latitude').points
            
            #%% remove sea points
            ## ns means "no-sea"
            ns_alt_0p5km = sub_alt_0p5km.data.copy()
            ns_alt_1p5km = sub_alt_1p5km.data.copy()
            ns_alt_4p4km = sub_alt_4p4km.data.copy()
            ns_alt_glm   =sub_alt_glm.data.copy()
            
            ns_alt_0p5km[sub_lbm_0p5km.data == 0] = np.nan
            ns_alt_1p5km[sub_lbm_1p5km.data == 0] = np.nan
            ns_alt_4p4km[sub_lbm_4p4km.data == 0] = np.nan
            ns_alt_glm[sub_lbm_glm.data == 0] =np.nan
            
            #%%
            print(np.nanmax(ns_alt_0p5km))
            #%% compute axis scalar distance 
            
            lon_dists_0p5km = climatology.haversine(lon_0p5km, np.ones(len(lon_0p5km))*lat_0p5km[0])
            #print(lon_dists_0p5km)
            #%%
             ## using "p" to denote plotting purpose
            pglmNE = (339,67.0)
            pglmSW = (335,65)
            
            polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            pregNE = (rot_db_lon[1], rot_db_lat[1])
            pregSW = (rot_db_lon[0], rot_db_lat[0])
            
            fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize= (12,12))
            fig.suptitle('UM smoothed orography')
            
            kwargs = {'cmap' : 'cividis', 'vmin' : 0, 'vmax' : 900, 'shading' : 'auto'}
            
            ax0.pcolormesh(lon_0p5km, lat_0p5km, ns_alt_0p5km, **kwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.scatter([maxlon_0p5km], [maxlat_0p5km], s = 600, marker = 'x', linewidths = 5, color = 'blue')
            
            ax1.pcolormesh(lon_1p5km, lat_1p5km, ns_alt_1p5km, **kwargs)
            ax1.set_title('1p5km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.scatter([maxlon_1p5km], [maxlat_1p5km], s = 600,marker = 'x', linewidths = 5, color = 'blue')
            
            ax2.pcolormesh(lon_4p4km, lat_4p4km, ns_alt_4p4km, **kwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.scatter([maxlon_4p4km], [maxlat_4p4km], s = 600,marker = 'x', linewidths = 5, color = 'blue')
            
            ax3.pcolormesh(lon_glm, lat_glm, ns_alt_glm.data, **kwargs)
            ax3.set_title('GLM')
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            
            norm = matplotlib.colors.Normalize(vmin =0, vmax = 900)
            altmap = matplotlib.cm.ScalarMappable(norm = norm, cmap = 'cividis')
            #fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.87])
            cbar = fig.colorbar(mappable=altmap,cax = cbar_ax)
            cbar.ax.set(ylabel = 'm') 
            
            plt.tight_layout()
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_orography_of_westfjords_with_storm_maxima_{glmtitle}.png')
            
        #%% difference plotting
        plot_differences =True
        if plot_differences:
            ## reload regional data
            temp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')[ss+chunk*ss-1]
            temp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')[ss+chunk*ss-1]
            temp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')[ss+chunk*ss-1]
            
            
            glmNE = (341,68)
            glmSW = (335,65)
            
            polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            regNE = (rot_db_lon[1], rot_db_lat[1])
            regSW = (rot_db_lon[0], rot_db_lat[0])  
            
            print(regSW, regNE)
            
            ## subset wind
            sub_temp_0p5km = climatology.subset_cube_by_coord(temp_0p5km,regSW, regNE)
            sub_temp_1p5km = climatology.subset_cube_by_coord(temp_1p5km,regSW, regNE)
            sub_temp_4p4km = climatology.subset_cube_by_coord(temp_4p4km,regSW, regNE)
            
            #%%
            print(sub_temp_0p5km)
            #%%
            ## preemptively remove glm and 4p4km landpoints
            sub_temp_glm.data[sub_lbm_glm.data == 1] = np.nan
            sub_temp_4p4km.data[sub_lbm_4p4km.data == 1] = np.nan
            ## regrid using iris (0p5km and 1p5km onto 4p4km grid and glm onto all reg)
            regrid_scheme = iris.analysis.Linear(extrapolation_mode = 'mask')
            nearest_scheme = iris.analysis.Linear(extrapolation_mode='mask') # actually use linear for now
            
            regrid_4p4_temp_0p5km = climatology.regrid_cubes( sub_temp_4p4km,sub_temp_0p5km, regrid_scheme)
            regrid_4p4_temp_1p5km = climatology.regrid_cubes( sub_temp_4p4km,sub_temp_1p5km, regrid_scheme)
            
            regrid_glm_temp_0p5km = climatology.regrid_cubes(sub_temp_glm, sub_temp_0p5km,nearest_scheme)
            regrid_glm_temp_1p5km = climatology.regrid_cubes(sub_temp_glm, sub_temp_1p5km,nearest_scheme)
            regrid_glm_temp_4p4km = climatology.regrid_cubes(sub_temp_glm, sub_temp_4p4km,nearest_scheme)
            
            #%%mask land points and calcualte differences
            nl_4p4_temp_0p5km = regrid_4p4_temp_0p5km.data.copy()
            nl_4p4_temp_1p5km = regrid_4p4_temp_1p5km.data.copy()
            #nl_4p4_temp_0p5km[sub_lbm_4p4km.data == 1] = np.nan
            #nl_4p4_temp_1p5km[sub_lbm_4p4km.data == 1] = np.nan
            diff_4p4_temp_0p5km = nl_temp_0p5km - nl_4p4_temp_0p5km
            diff_4p4_temp_1p5km = nl_temp_1p5km - nl_4p4_temp_1p5km
            
            nl_glm_temp_0p5km = regrid_glm_temp_0p5km.data.copy()
            nl_glm_temp_1p5km = regrid_glm_temp_1p5km.data.copy()
            nl_glm_temp_4p4km = regrid_glm_temp_4p4km.data.copy()
            nl_glm_temp_0p5km[sub_lbm_0p5km.data == 1] = np.nan
            nl_glm_temp_1p5km[sub_lbm_1p5km.data == 1] = np.nan
            nl_glm_temp_4p4km[sub_lbm_4p4km.data == 1] = np.nan
            
            diff_glm_temp_0p5km = nl_temp_0p5km - nl_glm_temp_0p5km
            diff_glm_temp_1p5km = nl_temp_1p5km - nl_glm_temp_1p5km
            diff_glm_temp_4p4km = nl_temp_4p4km - nl_glm_temp_4p4km
            #%%
            print(np.nanmax(diff_glm_temp_0p5km))
            print(np.nanmin(diff_glm_temp_0p5km))
            print(np.nanmax(diff_glm_temp_1p5km))
            print(np.nanmin(diff_glm_temp_1p5km))
            #%%
            plot_twopanels = True
            if plot_twopanels:
                matplotlib.rcParams.update({'font.size': 20})
                pglmNE = (340.45,67.4)
                pglmSW = (337,66)
                
                polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                pregNE = (rot_db_lon[1], rot_db_lat[1])
                pregSW = (rot_db_lon[0], rot_db_lat[0])  
                
                fig, (ax2,ax0) = plt.subplots(1,2, figsize = (18,10))
                gs = fig.add_gridspec(3,2)
                kwargs = {'vmin': -1.2,'vmax':1.2, 'cmap':'bwr', 'shading':'auto'}
                lb_kwargs = {'vmin': 0,'vmax':1, 'cmap':'Greys', 'shading':'auto'}
                fig.suptitle(f'{glmtitle} 1.5m temperature differences - 2018-03-19')
            
                 
                ax0.pcolormesh(lon_0p5km,lat_0p5km,diff_4p4_temp_0p5km, **kwargs )
                ax0.set_xticks([])
                ax0.set_yticks([])
                ax0.set_xlim(left = pregSW[0], right = pregNE[0])
                ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax0.set_title('0p5km - 4p4km')
                
                ax2.pcolormesh(lon_0p5km, lat_0p5km, diff_glm_temp_0p5km, **kwargs)
                ax2.set_xticks([])
                ax2.set_yticks([])
                ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax2.set_title('0p5km - GLM')
                
                div_norm = matplotlib.colors.Normalize(vmin = -1.2, vmax = 1.2)
                div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'bwr')
                cbar = fig.colorbar(mappable=div_map,ax = (ax2,ax0), orientation = 'horizontal', fraction = 0.15, pad = 0.01 )
                cbar.ax.set(xlabel = 'K')
                plt.subplots_adjust(left=0.125,
                        bottom=0.25, 
                        right=0.9, 
                        top=0.85, 
                        wspace=0.05, 
                        hspace=0.35)
            #%%
            
            plot_allpanels = False
            if plot_allpanels:
                matplotlib.rcParams.update({'font.size': 20})
                pglmNE = (340.8,67.6)
                pglmSW = (337,66)
                
                polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                pregNE = (rot_db_lon[1], rot_db_lat[1])
                pregSW = (rot_db_lon[0], rot_db_lat[0])  
                
                fig = plt.figure(figsize = (12,18))
                gs = fig.add_gridspec(3,2)
                kwargs = {'vmin': -1.2,'vmax':1.2, 'cmap':'bwr', 'shading':'auto'}
                lb_kwargs = {'vmin': 0,'vmax':1, 'cmap':'Greys', 'shading':'auto'}
                fig.suptitle(f'{glmtitle} avr. 1.5m temperature differences - 2018-03-19')
                
                ax0 = fig.add_subplot(gs[0,1])
                ax0.pcolormesh(lon_0p5km,lat_0p5km,diff_4p4_temp_0p5km.data, **kwargs )
                ax0.set_xticks([])
                ax0.set_yticks([])
                ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax0.set_title('0p5km - 4p4km')
                
                ax1 = fig.add_subplot(gs[1,1])
                ax1.pcolormesh(lon_1p5km,lat_1p5km,diff_4p4_temp_1p5km.data, **kwargs )
                ax1.set_xticks([])
                ax1.set_yticks([])
                ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax1.set_title('1p5km - 4p4km')
                
                ax2 = fig.add_subplot(gs[0,0])
                ax2.pcolormesh(lon_0p5km, lat_0p5km, diff_glm_temp_0p5km, **kwargs)
                ax2.set_xticks([])
                ax2.set_yticks([])
                ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax2.set_title('0p5km - GLM')
                
                ax3 = fig.add_subplot(gs[1,0])
                ax3.pcolormesh(lon_1p5km, lat_1p5km, diff_glm_temp_1p5km, **kwargs)
                ax3.set_xticks([])
                ax3.set_yticks([])
                ax3.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax3.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax3.set_title('1p5km - GLM')
                
                ax4 = fig.add_subplot(gs[2,0])
                ax4.pcolormesh(lon_4p4km, lat_4p4km, diff_glm_temp_4p4km, **kwargs)
                ax4.set_xticks([])
                ax4.set_yticks([])
                ax4.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax4.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax4.set_title('4p4km - GLM')
                
                div_norm = matplotlib.colors.Normalize(vmin = -1.2, vmax = 1.2)
                div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'bwr')
                
                #fig.subplots_adjust(right = 0.8)
                cbar_ax = fig.add_axes([0.99,0.03, 0.05,0.93])
                cbar = fig.colorbar(mappable=div_map,cax = cbar_ax)
                cbar.ax.set(ylabel = 'K')
                
                plt.tight_layout()
                
                plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_air_temperature_differences_3houravr_upper{glmtitle}_noland.png')