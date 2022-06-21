# -*- coding: utf-8 -*-
"""
Created on Fri May 13 14:43:55 2022

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
    for chunk in [4,5,6]:
        
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
        #glm_cubes = iris.load(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc')
        #%%
       # print(glm_cubes)
        #%%
        ss = 6 # ss : sample size; for half hourly outputs, a sample of 6 represents 3 hours
        ## all of these below are surface diagnostics
        q_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_specific_humidity_24hrs_pg_306.nc', 'specific_humidity')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        q_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_specific_humidity_24hrs_pg_306.nc', 'specific_humidity')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        q_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_specific_humidity_24hrs_pg_306.nc', 'specific_humidity')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        rh_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc','m01s03i245')[glmidx]
        temp_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc','air_temperature')[glmidx]
        msp_glm = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'air_pressure_at_sea_level')[glmidx]
        
        
        
        ## load binary land mask
        lbm_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_1_flt306.nc', 'land_binary_mask')
        
        ## load msp
        msp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        #%%
        print(rh_glm.units)
        print(temp_glm.units)
        print(msp_glm.units)
        
        
        #%% subsetting to region
        glmNE = (341,68)
        glmSW = (335,65)
        
        polelat = q_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = q_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
        rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        
        regNE = (rot_db_lon[1], rot_db_lat[1])
        regSW = (rot_db_lon[0], rot_db_lat[0])  
        
        print(regSW, regNE)
        
        ## subset wind
        sub_q_0p5km = climatology.subset_cube_by_coord(q_0p5km,regSW, regNE)
        sub_q_1p5km = climatology.subset_cube_by_coord(q_1p5km,regSW, regNE)
        sub_q_4p4km = climatology.subset_cube_by_coord(q_4p4km,regSW, regNE)
        sub_rh_glm   = climatology.subset_cube_by_coord(rh_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
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
        
        lon_0p5km = sub_q_0p5km.coord('grid_longitude').points
        lat_0p5km = sub_q_0p5km.coord('grid_latitude').points
        
        lon_1p5km = sub_q_1p5km.coord('grid_longitude').points
        lat_1p5km = sub_q_1p5km.coord('grid_latitude').points
        
        lon_4p4km = sub_q_4p4km.coord('grid_longitude').points
        lat_4p4km = sub_q_4p4km.coord('grid_latitude').points
        
        lon_glm = sub_rh_glm.coord('longitude').points
        lat_glm = sub_rh_glm.coord('latitude').points
        
        #%%
        #print(sub_msp_0p5km.data)
        msp_levels = np.arange(90000,130000, 100)
        #print(msp_levels)
        
        #%% check bounds
        print(np.nanmin(sub_q_0p5km.data), np.nanmax(sub_q_0p5km.data))
        
        
        #%% remove land points
        ## nl means "no-land"
        nl_q_0p5km = sub_q_0p5km.data.copy()
        nl_q_1p5km = sub_q_1p5km.data.copy()
        nl_q_4p4km = sub_q_4p4km.data.copy()
        
        nl_rh_glm   =sub_rh_glm.data.copy()
        
        
        nl_q_0p5km[sub_lbm_0p5km.data == 1] = np.nan
        nl_q_1p5km[sub_lbm_1p5km.data == 1] = np.nan
        nl_q_4p4km[sub_lbm_4p4km.data == 1] = np.nan
        nl_rh_glm[sub_lbm_glm.data == 1] = np.nan
        
        #%% check scale
        print(np.nanmax(nl_q_0p5km), np.nanmin(nl_q_0p5km))
        
        #%% calculate q for glm
        q_glm = climatology.rh_to_q_calculate(nl_rh_glm, sub_temp_glm.data , sub_msp_glm.data)
        print(np.shape(q_glm))
        
        #%% set new subsetting coordinates
        ## using "p" to denote plotting purpose
        
        plot_wakeregion = False
        if plot_wakeregion:
            pglmNE = (340,67.4)
            pglmSW = (337,66)
            
            polelat = q_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = q_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            pregNE = (rot_db_lon[1], rot_db_lat[1])
            pregSW = (rot_db_lon[0], rot_db_lat[0])  
            
            ## plot new
            
            fig = plt.figure(figsize = (12,12))
            gs = fig.add_gridspec(2,2)
            kwargs = {'cmap' : 'Blues', 'vmin' : 0.003*1000, 'vmax' : 0.005*1000, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            
            fig.suptitle(f'{titletime} avr. 1.5m specific humidity - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.pcolormesh(lon_0p5km, lat_0p5km, nl_q_0p5km.data*1000, **kwargs)
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.scatter([maxlon_0p5km], [maxlat_0p5km], s = 600, marker = 'x', linewidths = 5, color = 'green')
            
            ax1 = fig.add_subplot(gs[0,1])
            ax1.pcolormesh(lon_1p5km, lat_1p5km, nl_q_1p5km.data*1000, **kwargs)
            ax1.contour(lon_1p5km, lat_1p5km, sub_msp_1p5km.data, **mspkwargs)
            ax1.set_title('1p5km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.scatter([maxlon_1p5km], [maxlat_1p5km], s = 600,marker = 'x', linewidths = 5, color = 'green')
            
            ax2 = fig.add_subplot(gs[1,0])
            ax2.pcolormesh(lon_4p4km, lat_4p4km, nl_q_4p4km.data*1000, **kwargs)
            ax2.contour(lon_4p4km, lat_4p4km, sub_msp_4p4km.data, **mspkwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.scatter([maxlon_4p4km], [maxlat_4p4km], s = 600,marker = 'x', linewidths = 5, color = 'green')
            
            ax3 = fig.add_subplot(gs[1,1])#,projection=ccrs.PlateCarree())
            ax3.pcolormesh(lon_glm, lat_glm, q_glm*1000, **kwargs)
            ax3.contour(lon_glm, lat_glm, sub_msp_glm.data, **mspkwargs)
            ax3.set_title(f'GLM - {glmtitle}')
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            
            norm = matplotlib.colors.Normalize(vmin =0.003*1000, vmax = 0.005*1000)
            Tmap = matplotlib.cm.ScalarMappable(norm = norm, cmap = 'Blues')
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.87])
            cbar = fig.colorbar(mappable=Tmap,cax = cbar_ax)
            cbar.ax.set(ylabel = r'gkg$^{-1}$') 
            
            plt.tight_layout()
            
            #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_spec_humidity_3houravr_upper{glmtitle}_noland.png')
            
        
        #%% differences processing
        plot_differences = True
        if plot_differences:
            ## reload regional data
            q_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_specific_humidity_24hrs_pg_306.nc', 'specific_humidity')[ss+chunk*ss-1]
            q_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_specific_humidity_24hrs_pg_306.nc', 'specific_humidity')[ss+chunk*ss-1]
            q_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_specific_humidity_24hrs_pg_306.nc', 'specific_humidity')[ss+chunk*ss-1]
            
            glmNE = (341,68)
            glmSW = (335,65)
            
            polelat = q_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = q_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            regNE = (rot_db_lon[1], rot_db_lat[1])
            regSW = (rot_db_lon[0], rot_db_lat[0])  
            
            print(regSW, regNE)
            
            ## subset wind
            sub_q_0p5km = climatology.subset_cube_by_coord(q_0p5km,regSW, regNE)
            sub_q_1p5km = climatology.subset_cube_by_coord(q_1p5km,regSW, regNE)
            sub_q_4p4km = climatology.subset_cube_by_coord(q_4p4km,regSW, regNE)
            
            #%%
            ## preemptively remove glm and 4p4km land points
            sub_rh_glm.data[sub_lbm_glm.data == 1] = np.nan
            sub_q_4p4km.data[sub_lbm_4p4km.data == 1] = np.nan
            ## regrid using iris (0p5km and 1p5km onto 4p4km grid and glm onto all reg)
            area_scheme = iris.analysis.Linear(extrapolation_mode='mask')
            linear_scheme = iris.analysis.Linear(extrapolation_mode='mask')
            
            regrid_4p4_q_0p5km = climatology.regrid_cubes(sub_q_4p4km, sub_q_0p5km, area_scheme)
            regrid_4p4_q_1p5km = climatology.regrid_cubes( sub_q_4p4km,sub_q_1p5km, area_scheme)
            
            glm_q_0p5km = climatology.regrid_glm_and_compute_q(sub_lbm_0p5km, sub_rh_glm, sub_msp_glm, sub_temp_glm, linear_scheme)
            glm_q_1p5km = climatology.regrid_glm_and_compute_q(sub_lbm_1p5km, sub_rh_glm, sub_msp_glm, sub_temp_glm, linear_scheme)
            glm_q_4p4km = climatology.regrid_glm_and_compute_q(sub_lbm_4p4km, sub_rh_glm, sub_msp_glm, sub_temp_glm, linear_scheme)
            
            ## remove land points
            nl_4p4_q_0p5km = regrid_4p4_q_0p5km.data.copy()
            nl_4p4_q_1p5km = regrid_4p4_q_1p5km.data.copy()
            
            nl_glm_q_0p5km = glm_q_0p5km.copy()
            nl_glm_q_1p5km = glm_q_1p5km.copy()
            nl_glm_q_4p4km = glm_q_4p4km.copy()
            
            
            #nl_4p4_q_0p5km[sub_lbm_4p4km.data == 1] = np.nan
            #nl_4p4_q_1p5km[sub_lbm_4p4km.data == 1] = np.nan
            
            #nl_glm_q_0p5km[sub_lbm_0p5km.data == 1] = np.nan
            #nl_glm_q_1p5km[sub_lbm_1p5km.data == 1] = np.nan
            #nl_glm_q_4p4km[sub_lbm_4p4km.data == 1] = np.nan
            
            ## take differences
            diff_4p4_q_0p5km = nl_q_0p5km - nl_4p4_q_0p5km
            diff_4p4_q_1p5km = nl_q_1p5km - nl_4p4_q_1p5km
            
            diff_glm_q_0p5km = sub_q_0p5km.data - nl_glm_q_0p5km
            diff_glm_q_1p5km = sub_q_1p5km.data - nl_glm_q_1p5km
            diff_glm_q_4p4km = sub_q_4p4km.data - nl_glm_q_4p4km
            
            print('Slay')  
            #%% check ranges
            print(np.nanmin(diff_glm_q_0p5km), np.nanmax(diff_glm_q_0p5km))
            #%% plot processed fields
            matplotlib.rcParams.update({'font.size': 20})
            pglmNE = (340.8,67.6)
            pglmSW = (337,66)
            
            polelat = q_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = q_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            pregNE = (rot_db_lon[1], rot_db_lat[1])
            pregSW = (rot_db_lon[0], rot_db_lat[0])  
            
            fig = plt.figure(figsize = (12,18))
            gs = fig.add_gridspec(3,2)
            kwargs = {'vmin': -0.4,'vmax':0.4, 'cmap':'PiYG', 'shading':'auto'}
            lb_kwargs = {'vmin': 0,'vmax':1, 'cmap':'Greys', 'shading':'auto'}
            fig.suptitle(f'{glmtitle} avr. 1.5m spec. humidity differences - 2018-03-19')
            
            
            ax0 = fig.add_subplot(gs[0,1])
            #ax0.pcolormesh(lon_4p4km,lat_4p4km,sub_lbm_4p4km.data, **lb_kwargs )
            ax0.pcolormesh(lon_0p5km,lat_0p5km,diff_4p4_q_0p5km*1000, **kwargs )
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.set_title('0p5km - 4p4km')
            
            ax1 = fig.add_subplot(gs[1,1])
            #ax1.pcolormesh(lon_4p4km,lat_4p4km,sub_lbm_4p4km.data, **lb_kwargs )
            ax1.pcolormesh(lon_1p5km,lat_1p5km,diff_4p4_q_1p5km*1000, **kwargs )
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.set_title('1p5km - 4p4km')
            
            ax2 = fig.add_subplot(gs[0,0])
            ax2.pcolormesh(lon_0p5km, lat_0p5km, diff_glm_q_0p5km*1000, **kwargs)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.set_title('0p5km - GLM')
            
            ax3 = fig.add_subplot(gs[1,0])
            ax3.pcolormesh(lon_1p5km, lat_1p5km, diff_glm_q_1p5km*1000, **kwargs)
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax3.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax3.set_title('1p5km - GLM')
            
            ax4 = fig.add_subplot(gs[2,0])
            ax4.pcolormesh(lon_4p4km, lat_4p4km, diff_glm_q_4p4km*1000, **kwargs)
            ax4.set_xticks([])
            ax4.set_yticks([])
            ax4.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax4.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax4.set_title('4p4km - GLM')
            
            div_norm = matplotlib.colors.Normalize(vmin = -0.4, vmax = 0.4)
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'PiYG')
            
            #fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([0.99,0.03, 0.05,0.93])
            cbar = fig.colorbar(mappable=div_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'gkg$^{-1}$')
            
            plt.tight_layout()
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_specific_humidity_differences_3houravr_upper{glmtitle}_noland.png')