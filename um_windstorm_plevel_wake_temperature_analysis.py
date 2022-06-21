# -*- coding: utf-8 -*-
"""
Created on Mon May 16 18:09:10 2022

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
matplotlib.rcParams.update({'font.size': 20})
#warnings.filterwarnings("ignore", message="warning")
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s}" if plt.rcParams["text.usetex"] else f"{s}"
glmhidx = 4; reghidx = 12
#%%
wake_hoz_xsecs = True
if wake_hoz_xsecs:
    for chunk in [4]:
        if chunk == 4:
            titletime = '12.5h-15h'
            glmtitle = '1500hrs'
            glmidx = 0
            glm_file = 3
        elif chunk == 5:
            titletime = '15.5h-18h'
            glmtitle = '1800hrs'
            glmidx = 1
            glm_file = 3
        elif chunk == 6:
            titletime = '18.5h-21h'
            glmtitle = '2100hrs'
            glmidx = 0
            glm_file = 4
        ss = 3 # ss : sample size; for half hourly outputs, a sample of 6 represents 3 hours
        umdata_0p5km =  iris.load(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc')
        print(umdata_0p5km)
        umdata_0p5km = None
        #%%
        #glmcubes = iris.load(f'{um_datapath}RA1M_glm_3_flt306.nc')
        #%%
        #print(glmcubes)
        #%% load model data
        ## load temprature
        temp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'air_temperature')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)[reghidx]
        temp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'air_temperature')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)[reghidx]
        temp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'air_temperature')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)[reghidx]
        #wbp_glm = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'wet_bulb_potential_temperature')[glmidx, glmhidx]
        ## load binary land mask
        lbm_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        
        ## load geopotential height
        gph_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'geopotential_height')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)[reghidx]
        gph_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'geopotential_height')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)[reghidx]
        gph_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'geopotential_height')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)[reghidx]
        #gph_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'geopotential_height')
        
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
        #sub_wbp_glm   = climatology.subset_cube_by_coord(wbp_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        ## subset land-binary-mask
        sub_lbm_0p5km = climatology.subset_cube_by_coord(lbm_0p5km,regSW, regNE)
        sub_lbm_1p5km = climatology.subset_cube_by_coord(lbm_1p5km,regSW, regNE)
        sub_lbm_4p4km = climatology.subset_cube_by_coord(lbm_4p4km,regSW, regNE)
        #sub_lbm_glm   = climatology.subset_cube_by_coord(lbm_glm, glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        ## subset msp
        sub_gph_0p5km = climatology.subset_cube_by_coord(gph_0p5km,regSW, regNE)
        sub_gph_1p5km = climatology.subset_cube_by_coord(gph_1p5km,regSW, regNE)
        sub_gph_4p4km = climatology.subset_cube_by_coord(gph_4p4km,regSW, regNE)
        #sub_msp_glm = climatology.subset_cube_by_coord(msp_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        lon_0p5km = sub_temp_0p5km.coord('grid_longitude').points
        lat_0p5km = sub_temp_0p5km.coord('grid_latitude').points
        
        lon_1p5km = sub_temp_1p5km.coord('grid_longitude').points
        lat_1p5km = sub_temp_1p5km.coord('grid_latitude').points
        
        lon_4p4km = sub_temp_4p4km.coord('grid_longitude').points
        lat_4p4km = sub_temp_4p4km.coord('grid_latitude').points
        
        #lon_glm = sub_wbp_glm.coord('longitude').points
        #lat_glm = sub_wbp_glm.coord('latitude').points
        
    
        #%%
        plot_wholeregion = True
        if plot_wholeregion:
            
            
            
            ## set gph levels for contour plots
            gph_levels = np.arange(700,2500, 10)
            
            
                
                
            #%% check bounds 
            print(np.nanmin(sub_temp_0p5km.data), np.nanmax(sub_temp_0p5km.data))
            
            #%%
            plot_0p5kmpanel = True
            if plot_0p5kmpanel:
                fig, ax0=  plt.subplots(1,1, figsize = (8,10))    #plt.figure(figsize = (18,6))
                #gs = fig.add_gridspec(1,3)
                kwargs = {'cmap' : 'plasma', 'vmin' : 266, 'vmax' : 275, 'shading' : 'auto'}
                gphkwargs = {'levels':gph_levels, 'colors':'k'}
                
                fig.suptitle(f'{titletime} avr. 850hPa air temperature\n2018-03-19')
                
                ax0.pcolormesh(lon_0p5km, lat_0p5km, sub_temp_0p5km.data, **kwargs)
                ax0.contour(lon_0p5km, lat_0p5km, sub_gph_0p5km.data, **gphkwargs)
                ax0.contour(sub_lbm_0p5km.coords('grid_longitude')[0].points, sub_lbm_0p5km.coords('grid_latitude')[0].points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
                ax0.set_title('0p5km')
                ax0.set_xticks([])
                ax0.set_yticks([])
                
                div_norm = matplotlib.colors.Normalize(vmin = 266, vmax = 275)
                div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'plasma')
                
                #fig.subplots_adjust(right = 0.8)
                #cbar_ax = fig.add_axes([0.03,0.99, 0.93,0.05])
                cbar = fig.colorbar(mappable=div_map,ax = ax0, orientation = 'horizontal', fraction = 0.15, pad = 0.01 )
                cbar.ax.set(xlabel = 'K')
                
                plt.tight_layout()
            #%% plot
            plot_threepanels = False
            if plot_threepanels:
                pglmNE = (340,67.4)
                pglmSW = (337,66)
                
                polelat = temp_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = temp_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                pregNE = (rot_db_lon[1], rot_db_lat[1])
                pregSW = (rot_db_lon[0], rot_db_lat[0])  
                
                fig, (ax0,ax1,ax2)=  plt.subplots(1,3, figsize = (18,8))    #plt.figure(figsize = (18,6))
                #gs = fig.add_gridspec(1,3)
                kwargs = {'cmap' : 'plasma', 'vmin' : 266, 'vmax' : 275, 'shading' : 'auto'}
                gphkwargs = {'levels':gph_levels, 'colors':'k'}
                
                fig.suptitle(f'{titletime} avr. 850hPa temperature - 2018-03-19')
                
                #ax0 = fig.add_subplot(gs[0,0])
                ax0.pcolormesh(lon_0p5km, lat_0p5km, sub_temp_0p5km.data, **kwargs)
                ax0.contour(lon_0p5km, lat_0p5km, sub_gph_0p5km.data, **gphkwargs)
                ax0.contour(sub_lbm_0p5km.coords('grid_longitude')[0].points, sub_lbm_0p5km.coords('grid_latitude')[0].points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
                ax0.set_title('0p5km')
                ax0.set_xticks([])
                ax0.set_yticks([])
                #ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
                #ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
                
                #ax1 = fig.add_subplot(gs[0,1])
                ax1.pcolormesh(lon_1p5km, lat_1p5km, sub_temp_1p5km.data, **kwargs)
                ax1.contour(lon_1p5km, lat_1p5km, sub_gph_1p5km.data, **gphkwargs)
                ax1.contour(sub_lbm_1p5km.coords('grid_longitude')[0].points, sub_lbm_1p5km.coords('grid_latitude')[0].points, sub_lbm_1p5km.data, colors = 'green',levels = [1], linewidths = 4)
                ax1.set_title('1p5km')
                ax1.set_xticks([])
                ax1.set_yticks([])
                
                #ax2 = fig.add_subplot(gs[0,2])
                ax2.pcolormesh(lon_4p4km, lat_4p4km, sub_temp_4p4km.data, **kwargs)
                ax2.contour(lon_4p4km, lat_4p4km, sub_gph_4p4km.data, **gphkwargs)
                ax2.contour(sub_lbm_4p4km.coords('grid_longitude')[0].points, sub_lbm_4p4km.coords('grid_latitude')[0].points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
                ax2.set_title('4p4km')
                ax2.set_xticks([])
                ax2.set_yticks([])
                
                div_norm = matplotlib.colors.Normalize(vmin = 266, vmax = 275)
                div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'plasma')
                
                #fig.subplots_adjust(right = 0.8)
                #cbar_ax = fig.add_axes([0.03,0.99, 0.93,0.05])
                cbar = fig.colorbar(mappable=div_map,ax = (ax0,ax1,ax2), orientation = 'horizontal', fraction = 0.15, pad = 0.01 )
                cbar.ax.set(xlabel = 'K')
                plt.subplots_adjust(left=0.125,
                        bottom=0.25, 
                        right=0.9, 
                        top=0.85, 
                        wspace=0.05, 
                        hspace=0.35)
                
                
                
                #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/Plevel/RA1M_UM_850hPa_temperature_3houravr_upper{glmtitle}.png')
                print('Slay')
                ## plot end
            
        #%%
        plot_wholeregion_diff = False
        if plot_wholeregion_diff:
            
            ## reload data
            ## load temprature
            temp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'air_temperature')[ss+chunk*ss-1][reghidx]
            temp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'air_temperature')[ss+chunk*ss-1][reghidx]
            temp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'air_temperature')[ss+chunk*ss-1][reghidx]
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
            
            
            #%% regrid take difference
            ## preemptively remove glm and 4p4km landpoints
            #sub_temp_4p4km.data[sub_lbm_4p4km.data == 1] = np.nan
            ## regrid using iris (0p5km and 1p5km onto 4p4km grid and glm onto all reg)
            regrid_scheme = iris.analysis.Linear(extrapolation_mode = 'mask')
            regrid_4p4_temp_0p5km = climatology.regrid_cubes( sub_temp_4p4km,sub_temp_0p5km, regrid_scheme)
            regrid_4p4_temp_1p5km = climatology.regrid_cubes( sub_temp_4p4km,sub_temp_1p5km, regrid_scheme)
            
            diff_4p4_temp_0p5km = sub_temp_0p5km.data - regrid_4p4_temp_0p5km.data
            diff_4p4_temp_1p5km = sub_temp_1p5km.data - regrid_4p4_temp_1p5km.data
            
            #%%
            print(np.nanmin(diff_4p4_temp_0p5km), np.nanmax(diff_4p4_temp_0p5km))
            print(np.nanmin(diff_4p4_temp_1p5km), np.nanmax(diff_4p4_temp_1p5km))
            #%% plot
            fig, (ax0,ax1)=  plt.subplots(1,2, figsize = (12,8))    #plt.figure(figsize = (18,6))
            #gs = fig.add_gridspec(1,3)
            kwargs = {'cmap' : 'bwr', 'vmin' : -3.2, 'vmax' : 3.2, 'shading' : 'auto'}
            fig.suptitle(f'{glmtitle} 850hPa temperature differences - 2018-03-19')
            
            ax0.pcolormesh(lon_0p5km, lat_0p5km, diff_4p4_temp_0p5km, **kwargs)
            ax0.contour(sub_lbm_4p4km.coords('grid_longitude')[0].points, sub_lbm_4p4km.coords('grid_latitude')[0].points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            ax0.set_title('0p5km - 4p4km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            
            ax1.pcolormesh(lon_1p5km, lat_1p5km, diff_4p4_temp_1p5km, **kwargs)
            ax1.contour(sub_lbm_4p4km.coords('grid_longitude')[0].points, sub_lbm_4p4km.coords('grid_latitude')[0].points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            ax1.set_title('1p5km - 4p4km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            
            div_norm = matplotlib.colors.Normalize(vmin = -3.2, vmax = 3.2)
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'bwr')
            
            #fig.subplots_adjust(right = 0.8)
            #cbar_ax = fig.add_axes([0.03,0.99, 0.93,0.05])
            cbar = fig.colorbar(mappable=div_map,ax = (ax0,ax1), orientation = 'horizontal', fraction = 0.15, pad = 0.01 )
            cbar.ax.set(xlabel = 'K')
            
            plt.subplots_adjust(left=0.125,
                    bottom=0.25, 
                    right=0.9, 
                    top=0.85, 
                    wspace=0.05, 
                    hspace=0.35)
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/Plevel/RA1M_UM_850hPa_temperature_differences_upper{glmtitle}.png')