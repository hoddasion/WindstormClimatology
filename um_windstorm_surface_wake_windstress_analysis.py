# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:50:45 2022

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
matplotlib.rcParams.update({'font.size': 18})
#warnings.filterwarnings("ignore", message="warning")

#%%
wake_hoz_xsecs = True
if wake_hoz_xsecs:
    for chunk in [5]:
    
        #chunk = 5
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
        glm_cubes = iris.load(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc')
       
            #%%
        ss = 6 # ss : sample size; for half hourly outputs, a sample of 6 represents 3 hours
        ## all of these below are surface diagnostics
        ## also perform 3 hour averages
        estress_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_surface_downward_eastward_stress_24hrs_pg_306.nc', 'surface_downward_eastward_stress')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        nstress_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_surface_downward_northward_stress_24hrs_pg_306.nc', 'surface_downward_northward_stress')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        estress_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_surface_downward_eastward_stress_24hrs_pg_306.nc', 'surface_downward_eastward_stress')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        nstress_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_surface_downward_northward_stress_24hrs_pg_306.nc','surface_downward_northward_stress')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        estress_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_surface_downward_eastward_stress_24hrs_pg_306.nc', 'surface_downward_eastward_stress')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        nstress_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_surface_downward_northward_stress_24hrs_pg_306.nc', 'surface_downward_northward_stress')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        
        
        glm_xwind = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'm01s03i225')#glm_cubes[30]
        glm_ywind = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'm01s03i226')#glm_cubes[32]
        
        ## load binary land mask
        lbm_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
        lbm_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_1_flt306.nc', 'land_binary_mask')
        
        ## load msp
        msp_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        msp_glm = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'air_pressure_at_sea_level')
        #%% load satellite data
        path_sat = 'D:/Project/Climatology/Data/'
        name_sat = 'KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_12_ASC_2018.nc'

        satdata = xr.open_dataset(f'{path_sat}{name_sat}')
        print(satdata)
        ds_wsp = satdata.wind_speed
        #ds_wdir = dataset.wind_dir
        ds_lat = satdata.lat
        ds_lon = satdata.lon
        #print(ds_wsp[77])


        points_wsp = np.array(ds_wsp)
        points_lat = np.array(ds_lat)
        points_lon = np.array(ds_lon)
        points_v = np.array(satdata.northward_wind)
        points_u = np.array(satdata.eastward_wind)
        points_rho = np.array(satdata.air_density)
        measurement_time = np.array(satdata.measurement_time)
        points_wdiv = np.array(satdata.wind_divergence)
        #print(points_rho)