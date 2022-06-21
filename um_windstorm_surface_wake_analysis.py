# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:06:35 2022

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
    for chunk in [4,5,6]:
    
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
            #%%
        ss = 6 # ss : sample size; for half hourly outputs, a sample of 6 represents 3 hours
        ## all of these below are surface diagnostics
        ## also perform 3 hour averages
        xwind_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        ywind_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_y_wind_24hrs_pg_306.nc', 'y_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        xwind_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        ywind_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_y_wind_24hrs_pg_306.nc', 'y_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        xwind_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_x_wind_24hrs_pg_306.nc', 'x_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        ywind_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_y_wind_24hrs_pg_306.nc', 'y_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        
        
        glm_xwind = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'm01s03i225')[glmidx]#glm_cubes[30]
        glm_ywind = iris.load_cube(f'{um_datapath}RA1M_glm_{glm_file}_flt306.nc', 'm01s03i226')[glmidx]#glm_cubes[32]
        
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
        
        
        #%% subsetting to region
        glmNE = (341,68)
        glmSW = (335,65)
        
        polelat = xwind_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = xwind_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
        rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        
        regNE = (rot_db_lon[1], rot_db_lat[1])
        regSW = (rot_db_lon[0], rot_db_lat[0])  
        
        #print(regSW, regNE)
        
        
        ## subset wind
        sub_xwind_0p5km = climatology.subset_cube_by_coord(xwind_0p5km,regSW, regNE)
        sub_ywind_0p5km = climatology.subset_cube_by_coord(ywind_0p5km,regSW, regNE)
        
        sub_xwind_1p5km = climatology.subset_cube_by_coord(xwind_1p5km,regSW, regNE)
        sub_ywind_1p5km = climatology.subset_cube_by_coord(ywind_1p5km,regSW, regNE)
        
        sub_xwind_4p4km = climatology.subset_cube_by_coord(xwind_4p4km,regSW, regNE)
        sub_ywind_4p4km = climatology.subset_cube_by_coord(ywind_4p4km,regSW, regNE)
        
        sub_xwind_glm = climatology.subset_cube_by_coord(glm_xwind,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        sub_ywind_glm = climatology.subset_cube_by_coord(glm_ywind,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
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
        
        #print(sub_xwind_glm)
        
        lon_0p5km = sub_xwind_0p5km.coord('grid_longitude').points
        lat_0p5km = sub_xwind_0p5km.coord('grid_latitude').points
        
        lon_1p5km = sub_xwind_1p5km.coord('grid_longitude').points
        lat_1p5km = sub_xwind_1p5km.coord('grid_latitude').points
        
        lon_4p4km = sub_xwind_4p4km.coord('grid_longitude').points
        lat_4p4km = sub_xwind_4p4km.coord('grid_latitude').points
        
        lon_glm = sub_xwind_glm.coord('longitude').points
        lat_glm = sub_xwind_glm.coord('latitude').points
        
        
        
        #%% subset measurement time field from ASCAT data
        check_satellite_time = False
        if check_satellite_time:
            condlon = (points_lon < glmNE[0]) * (points_lon > glmSW[0])
            condlat = (points_lat < glmNE[1]) * (points_lat > glmSW[1])
            sub_mestime =measurement_time[77][condlat,:]
            sub_mestime =sub_mestime[:, condlon]
            print(sub_mestime)
            #print(np.nanmean(sub_mestime))
        
        #%%
        #print(sub_msp_0p5km.data)
        msp_levels = np.arange(90000,130000, 100)
        #print(msp_levels)
        
        
        
        #%% cursory plotting
        quiv0p5 = 50
        quiv1p5 = 10
        quiv4p4 = 8
        quivglm = 4
        quivsat = 3
        
        
            
        #%% compute magnitudes
        wsp_0p5km = climatology.quadmag(sub_xwind_0p5km.data, sub_ywind_0p5km.data[:-1,:])
        wsp_1p5km = climatology.quadmag(sub_xwind_1p5km.data, sub_ywind_1p5km.data[:-1,:])
        wsp_4p4km = climatology.quadmag(sub_xwind_4p4km.data, sub_ywind_4p4km.data)
        wsp_glm = climatology.quadmag(sub_xwind_glm.data, sub_ywind_glm.data)
        #%%
        plot_wholeregion_3panels = True
        if plot_wholeregion_3panels:
            
            print('blub')
            
            fig, (ax0,ax2,ax3) = plt.subplots(3,1, figsize = (6,15))
            kwargs = {'cmap' : 'Oranges', 'vmin' : 0, 'vmax' : 20, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            ax0.pcolormesh(lon_0p5km, lat_0p5km, climatology.quadmag(sub_xwind_0p5km.data, sub_ywind_0p5km.data[:-1,:]), **kwargs)
            ax0.quiver(lon_0p5km[::quiv0p5], lat_0p5km[::quiv0p5], sub_xwind_0p5km.data[::quiv0p5,::quiv0p5], sub_ywind_0p5km.data[::quiv0p5,::quiv0p5])
            ax0.contour(lon_0p5km, lat_0p5km, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim( right = regNE[0])
            ax0.set_ylim( top = regNE[1])
            
            ax2.pcolormesh(lon_4p4km, lat_4p4km, climatology.quadmag(sub_xwind_4p4km.data, sub_ywind_4p4km.data), **kwargs)
            ax2.quiver(lon_4p4km[::quiv4p4], lat_4p4km[::quiv4p4], sub_xwind_4p4km.data[::quiv4p4,::quiv4p4], sub_ywind_4p4km.data[::quiv4p4,::quiv4p4])
            ax2.contour(lon_4p4km, lat_4p4km, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            ax2.contour(lon_4p4km, lat_4p4km, sub_msp_4p4km.data, **mspkwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            
            ax3.pcolormesh(lon_glm, lat_glm, climatology.quadmag(sub_xwind_glm.data, sub_ywind_glm.data), **kwargs)
            ax3.quiver(lon_glm[::quivglm], lat_glm[::quivglm], sub_xwind_glm.data[::quivglm,::quivglm], sub_ywind_glm.data[::quivglm,::quivglm])
            ax3.contour(lon_glm, lat_glm, sub_msp_glm.data[:,:-1], **mspkwargs)
            ax3.contour(sub_lbm_glm.coord('longitude').points, sub_lbm_glm.coord('latitude').points, sub_lbm_glm.data, colors = 'green', levels = [1], linewidths = 4)
            ax3.set_title(f'GLM - {glmtitle}')
            ax3.set_xticks([])
            ax3.set_yticks([])
            
            diff_norm = matplotlib.colors.Normalize(vmin =0, vmax = 20)
            diff_map = matplotlib.cm.ScalarMappable(norm = diff_norm, cmap = 'Oranges')
            cbar2 = fig.colorbar(mappable=diff_map,ax = (ax0,ax2,ax3), orientation = 'vertical', fraction = 0.15, pad = 0.06)
            cbar2.ax.set(ylabel = r'ms$^{-1}$')
            
            plt.subplots_adjust(left=0,
                            right=0.75, 
                            top = 0.9,
                            bottom = 0.1, 
                            wspace=0.05, 
                            hspace=0.1)
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_wind_wASCAT_3houravr_upper{glmtitle}_3panel.png')
        
        #%%
        plot_wholeregion = False
        if plot_wholeregion:
            fig = plt.figure(figsize = (12,18))
            gs = fig.add_gridspec(3,2)
            kwargs = {'cmap' : 'Oranges', 'vmin' : 0, 'vmax' : 20, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            fig.suptitle(f'{titletime} avr. 10m horizontal windspeed - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.pcolormesh(lon_0p5km, lat_0p5km, climatology.quadmag(sub_xwind_0p5km.data, sub_ywind_0p5km.data[:-1,:]), **kwargs)
            ax0.quiver(lon_0p5km[::quiv0p5], lat_0p5km[::quiv0p5], sub_xwind_0p5km.data[::quiv0p5,::quiv0p5], sub_ywind_0p5km.data[::quiv0p5,::quiv0p5])
            ax0.contour(lon_0p5km, lat_0p5km, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim( right = regNE[0])
            ax0.set_ylim( top = regNE[1])
            #ax0.scatter([maxlon_0p5km], [maxlat_0p5km], s = 600)#, marker = 'x')
            
            ax1 = fig.add_subplot(gs[0,1])
            ax1.pcolormesh(lon_1p5km, lat_1p5km, climatology.quadmag(sub_xwind_1p5km.data, sub_ywind_1p5km.data[:-1,:]), **kwargs)
            ax1.quiver(lon_1p5km[::quiv1p5], lat_1p5km[::quiv1p5], sub_xwind_1p5km.data[::quiv1p5,::quiv1p5], sub_ywind_1p5km.data[::quiv1p5,::quiv1p5])
            ax1.contour(lon_1p5km, lat_1p5km, sub_lbm_1p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax1.contour(lon_1p5km, lat_1p5km, sub_msp_1p5km.data, **mspkwargs)
            ax1.set_title('1p5km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            #ax1.scatter([maxlon_1p5km], [maxlat_1p5km], s = 600)
            
            ax2 = fig.add_subplot(gs[1,0])
            ax2.pcolormesh(lon_4p4km, lat_4p4km, climatology.quadmag(sub_xwind_4p4km.data, sub_ywind_4p4km.data), **kwargs)
            ax2.quiver(lon_4p4km[::quiv4p4], lat_4p4km[::quiv4p4], sub_xwind_4p4km.data[::quiv4p4,::quiv4p4], sub_ywind_4p4km.data[::quiv4p4,::quiv4p4])
            ax2.contour(lon_4p4km, lat_4p4km, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            ax2.contour(lon_4p4km, lat_4p4km, sub_msp_4p4km.data, **mspkwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            
            ax3 = fig.add_subplot(gs[2,0])#,projection=ccrs.PlateCarree())
            ax3.pcolormesh(lon_glm, lat_glm, climatology.quadmag(sub_xwind_glm.data, sub_ywind_glm.data), **kwargs)
            ax3.quiver(lon_glm[::quivglm], lat_glm[::quivglm], sub_xwind_glm.data[::quivglm,::quivglm], sub_ywind_glm.data[::quivglm,::quivglm])
            ax3.contour(lon_glm, lat_glm, sub_msp_glm.data[glmidx,:,:-1], **mspkwargs)
            ax3.contour(sub_lbm_glm.coord('longitude').points, sub_lbm_glm.coord('latitude').points, sub_lbm_glm.data, colors = 'green', levels = [1], linewidths = 4)
            ax3.set_title(f'GLM - {glmtitle}')
            ax3.set_xticks([])
            ax3.set_yticks([])
            #ax3.coastlines()
            #ax3.contour(lon_4p4km, lat_4p4km, sub_lbm_4p4km.data)
            
            ax4 = fig.add_subplot(gs[2,1])#,projection=ccrs.PlateCarree())
            ax4.pcolormesh(points_lon, points_lat,points_wsp[77], **kwargs)
            ax4.quiver(points_lon[::quivsat], points_lat[::quivsat], points_u[77,::quivsat,::quivsat], points_v[77,::quivsat,::quivsat],headlength = 6, headwidth = 4)
            ax4.contour(points_lon, points_lat, points_rho[77], colors = 'k')
            #ax4.coastlines()
            ax4.set_xlim(left = glmSW[0], right = glmNE[0] )
            ax4.set_ylim(bottom = glmSW[1], top = glmNE[1])
            ax4.set_title('ASCAT 12.5km Coastal')
            ax4.set_xticks([])
            ax4.set_yticks([])
            
            wsp_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 20)
            wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.02, 0.05,0.91])
            cbar = fig.colorbar(mappable=wsp_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'$ms^{-1}$') 
            
            plt.tight_layout()
            
            #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_wind_wASCAT_3houravr_upper{glmtitle}.png')
        #%% take location of windspeed maxima (hopefully the windstorm)
        
        ## narrow search by subsetting data to area over storm only
        stormreg0p5 = wsp_0p5km[(lat_0p5km >0) * (lat_0p5km < 0.6),:]
        stormreg0p5 = stormreg0p5[:,(lon_0p5km >360) * (lon_0p5km < 360.4)]
        stormreg1p5 = wsp_1p5km[(lat_1p5km >0) * (lat_1p5km < 0.6),:]
        stormreg1p5 = stormreg1p5[:,(lon_1p5km >360) * (lon_1p5km < 360.4)]
        stormreg4p4 = wsp_4p4km[(lat_4p4km >0) * (lat_4p4km < 0.6),:]
        stormreg4p4 = stormreg4p4[:,(lon_4p4km >360) * (lon_4p4km < 360.4)]
        
        wspmax_0p5km = np.nanmax(stormreg0p5)
        wspmax_1p5km = np.nanmax(stormreg1p5)
        wspmax_4p4km = np.nanmax(stormreg4p4)
        #print(wspmax_0p5km, wspmax_1p5km, wspmax_4p4km)
        lonmesh_0p5km, latmesh_0p5km = np.meshgrid(lon_0p5km, lat_0p5km)
        lonmesh_1p5km, latmesh_1p5km = np.meshgrid(lon_1p5km, lat_1p5km)
        lonmesh_4p4km, latmesh_4p4km = np.meshgrid(lon_4p4km, lat_4p4km)
        #print(np.shape(wsp_0p5km), np.shape(lonmesh_0p5km))
        
        
        maxlon_0p5km = lonmesh_0p5km[wsp_0p5km == wspmax_0p5km]
        maxlat_0p5km = latmesh_0p5km[wsp_0p5km == wspmax_0p5km]
        maxlon_1p5km = lonmesh_1p5km[wsp_1p5km == wspmax_1p5km]
        maxlat_1p5km = latmesh_1p5km[wsp_1p5km == wspmax_1p5km]
        maxlon_4p4km = lonmesh_4p4km[wsp_4p4km == wspmax_4p4km]
        maxlat_4p4km = latmesh_4p4km[wsp_4p4km == wspmax_4p4km]
        
        print(titletime,(maxlon_0p5km, maxlat_0p5km), (maxlon_1p5km,maxlat_1p5km), (maxlon_4p4km, maxlat_4p4km))
        #%%
        print(sub_lbm_glm.data[:,:-1])
        #%% remove land points
        ## nl means "no-land"
        nl_wsp_0p5km = wsp_0p5km.copy()
        nl_wsp_1p5km = wsp_1p5km.copy()
        nl_wsp_4p4km = wsp_4p4km.copy()
        nl_wsp_glm   = wsp_glm.copy()
        nl_wsp_0p5km[sub_lbm_0p5km.data == 1] = np.nan
        nl_wsp_1p5km[sub_lbm_1p5km.data == 1] = np.nan
        nl_wsp_4p4km[sub_lbm_4p4km.data == 1] = np.nan
        nl_wsp_glm[sub_lbm_glm.data[:,:-1] == 1] = np.nan
        
        #%% set new subsetting coordinates
        ## using "p" to denote plotting purpose
        plot_wakeregion = True
        if plot_wakeregion:
            pglmNE = (340,67.4)
            pglmSW = (337,66)
            
            polelat = xwind_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = xwind_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            pregNE = (rot_db_lon[1], rot_db_lat[1])
            pregSW = (rot_db_lon[0], rot_db_lat[0])  
            
            ## plot new
            
            quiv0p5 = 20
            quiv1p5 = 7
            quiv4p4 = 3
            quivglm = 1
            quivsat = 2
            
            
            
            fig = plt.figure(figsize = (12,18))
            gs = fig.add_gridspec(3,2)
            kwargs = {'cmap' : 'Oranges', 'vmin' : 0, 'vmax' : 20, 'shading' : 'auto'}
            mspkwargs = {'levels':msp_levels, 'colors':'k'}
            fig.suptitle(f'{titletime} avr. 10m horizontal windspeed - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.pcolormesh(lon_0p5km, lat_0p5km, nl_wsp_0p5km, **kwargs)
            ax0.quiver(lon_0p5km[::quiv0p5], lat_0p5km[::quiv0p5], sub_xwind_0p5km.data[::quiv0p5,::quiv0p5], sub_ywind_0p5km.data[::quiv0p5,::quiv0p5])
            #ax0.contour(lon_0p5km, lat_0p5km, sub_lbm_0p5km.data, colors = 'green')
            ax0.set_title('0p5km')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.scatter([maxlon_0p5km], [maxlat_0p5km], s = 600, marker = 'x', linewidths = 5, color = 'blue')
            ax0.contour(lon_0p5km, lat_0p5km, sub_msp_0p5km.data, **mspkwargs)
            
            ax1 = fig.add_subplot(gs[0,1])
            ax1.pcolormesh(lon_1p5km, lat_1p5km, nl_wsp_1p5km, **kwargs)
            ax1.quiver(lon_1p5km[::quiv1p5], lat_1p5km[::quiv1p5], sub_xwind_1p5km.data[::quiv1p5,::quiv1p5], sub_ywind_1p5km.data[::quiv1p5,::quiv1p5])
            #ax1.contour(lon_1p5km, lat_1p5km, sub_lbm_1p5km.data, colors = 'green')
            ax1.contour(lon_1p5km, lat_1p5km, sub_msp_1p5km.data, **mspkwargs)
            ax1.set_title('1p5km')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.scatter([maxlon_1p5km], [maxlat_1p5km], s = 600,marker = 'x', linewidths = 5, color = 'blue')
            
            ax2 = fig.add_subplot(gs[1,0])
            ax2.pcolormesh(lon_4p4km, lat_4p4km, nl_wsp_4p4km, **kwargs)
            ax2.quiver(lon_4p4km[::quiv4p4], lat_4p4km[::quiv4p4], sub_xwind_4p4km.data[::quiv4p4,::quiv4p4], sub_ywind_4p4km.data[::quiv4p4,::quiv4p4])
            #ax2.contour(lon_4p4km, lat_4p4km, sub_lbm_4p4km.data, colors = 'green')
            ax2.contour(lon_4p4km, lat_4p4km, sub_msp_4p4km.data, **mspkwargs)
            ax2.set_title('4p4km')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.scatter([maxlon_4p4km], [maxlat_4p4km], s = 600,marker = 'x', linewidths = 5, color = 'blue')
            
            ax3 = fig.add_subplot(gs[2,0])#,projection=ccrs.PlateCarree())
            ax3.pcolormesh(lon_glm, lat_glm, nl_wsp_glm, **kwargs)
            ax3.quiver(lon_glm[::quivglm], lat_glm[::quivglm], sub_xwind_glm.data[::quivglm,::quivglm], sub_ywind_glm.data[::quivglm,::quivglm])
            ax3.contour(lon_glm, lat_glm, sub_msp_glm.data[:,:-1], **mspkwargs)
            ax3.set_title(f'GLM - {glmtitle}')
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            ax3.set_xticks([])
            ax3.set_yticks([])
            #ax3.coastlines()
            #ax3.contour(lon_4p4km, lat_4p4km, sub_lbm_4p4km.data)
            
            ax4 = fig.add_subplot(gs[2,1])#,projection=ccrs.PlateCarree())
            ax4.pcolormesh(points_lon, points_lat,points_wsp[77], **kwargs)
            ax4.quiver(points_lon[::quivsat], points_lat[::quivsat], points_u[77,::quivsat,::quivsat], points_v[77,::quivsat,::quivsat],headlength = 6, headwidth = 4)
            #ax4.coastlines()
            ax4.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax4.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            ax4.set_title('ASCAT 12.5km Coastal')
            ax4.set_xticks([])
            ax4.set_yticks([])
            
            wsp_norm = matplotlib.colors.Normalize(vmin = 0, vmax = 20)
            wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
            
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.02, 0.05,0.91])
            cbar = fig.colorbar(mappable=wsp_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'$ms^{-1}$') 
            
            plt.tight_layout()
            #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_wind_wASCAT_3houravr_upper{glmtitle}_noland.png')
            
        #%%
        horizontal_convergence = False
        if horizontal_convergence:
            
            import climatology
            lonmesh0p5km, latmesh0p5km = np.meshgrid(lon_0p5km, lat_0p5km)
            lonmesh1p5km, latmesh1p5km = np.meshgrid(lon_1p5km, lat_1p5km)
            lonmesh4p4km, latmesh4p4km = np.meshgrid(lon_4p4km, lat_4p4km)
            lonmeshglm, latmeshglm = np.meshgrid(lon_glm, lat_glm)
            div_0p5km = climatology.divergence_horizontal(sub_xwind_0p5km.data, sub_ywind_0p5km.data[:-1], lonmesh0p5km, latmesh0p5km)
            
            div_1p5km = climatology.divergence_horizontal(sub_xwind_1p5km.data, sub_ywind_1p5km.data[:-1], lonmesh1p5km, latmesh1p5km)
            div_4p4km = climatology.divergence_horizontal(sub_xwind_4p4km.data, sub_ywind_4p4km.data, lonmesh4p4km, latmesh4p4km)
            div_glm = climatology.divergence_horizontal(sub_xwind_glm.data[glmidx], sub_ywind_glm.data[glmidx], lonmeshglm, latmeshglm)
            
            ##%% remove land points from div fields
            div_0p5km[sub_lbm_0p5km.data[:-1,:-1] == 1] = np.nan
            div_1p5km[sub_lbm_1p5km.data[:-1,:-1] == 1] = np.nan
            div_4p4km[sub_lbm_4p4km.data[:-1,:-1] == 1] = np.nan
            div_glm[sub_lbm_glm.data[:-1,:-2] == 1]     = np.nan
            #print(div_glm)
            #%%
            print(np.shape(lon_glm), np.shape(lat_glm), np.shape(div_glm))
            #%%
            ## using "p" to denote plotting purpose
            pglmNE = (340,67.4)
            pglmSW = (337,66)
            
            polelat = xwind_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = xwind_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            pregNE = (rot_db_lon[1], rot_db_lat[1])
            pregSW = (rot_db_lon[0], rot_db_lat[0])  
            
            
            fig = plt.figure(figsize = (12,18))
            gs = fig.add_gridspec(3,2)
            kwargs = {'cmap' : 'RdYlBu_r',  'shading' : 'auto', 'vmin':-0.0025, 'vmax':0.0025}
            fig.suptitle(f'{titletime} avr. 10m horizontal wind divergence - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.pcolormesh(lon_0p5km[:-1], lat_0p5km[:-1], div_0p5km, **kwargs)
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.scatter([maxlon_0p5km], [maxlat_0p5km], s = 600, marker = 'x', linewidths = 5, color = 'green')
            ax0.set_title('0p5km')
            
            ax1 = fig.add_subplot(gs[0,1])
            ax1.pcolormesh(lon_1p5km[:-1], lat_1p5km[:-1], div_1p5km, **kwargs)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.scatter([maxlon_1p5km], [maxlat_1p5km], s = 600, marker = 'x', linewidths = 5, color = 'green')
            ax1.set_title('1p5km')
            
            ax2 = fig.add_subplot(gs[1,0])
            ax2.pcolormesh(lon_4p4km[:-1], lat_4p4km[:-1], div_4p4km, **kwargs)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.scatter([maxlon_4p4km], [maxlat_4p4km], s = 600, marker = 'x', linewidths = 5, color = 'green')
            ax2.set_title('4p4km')
            
            ax3 = fig.add_subplot(gs[2,0])
            ax3.pcolormesh(lon_glm[:-1], lat_glm[:-1], div_glm, **kwargs)
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            ax3.set_title(f'GLM - {glmtitle}')
            
            ax4 = fig.add_subplot(gs[2,1])#,projection=ccrs.PlateCarree())
            ax4.pcolormesh(points_lon, points_lat,points_wdiv[77], **kwargs)
            
            ax4.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax4.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            ax4.set_title('ASCAT 12.5km Coastal')
            ax4.set_xticks([])
            ax4.set_yticks([])
            
            div_norm = matplotlib.colors.Normalize(vmin = -0.0025, vmax = 0.0025)
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'RdYlBu_r')
            
            #fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.02, 0.05,0.91])
            cbar = fig.colorbar(mappable=div_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'$s^{-1}$') 
            
            plt.tight_layout()
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_2Dwind_divergence_wASCAT_3houravr_upper{glmtitle}_noland.png')
            
        #%%
        difference_plotting = False
        if difference_plotting:
            ## reload data
            xwind_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')[ss+chunk*ss-1]
            ywind_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_um_y_wind_24hrs_pg_306.nc', 'y_wind')[ss+chunk*ss-1]
            
            xwind_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')[ss+chunk*ss-1]
            ywind_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_um_y_wind_24hrs_pg_306.nc', 'y_wind')[ss+chunk*ss-1]
            
            xwind_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_x_wind_24hrs_pg_306.nc', 'x_wind')[ss+chunk*ss-1]
            ywind_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_um_y_wind_24hrs_pg_306.nc', 'y_wind')[ss+chunk*ss-1]
            
            
            #%% subsetting to region
            glmNE = (341,68)
            glmSW = (335,65)
            
            polelat = xwind_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = xwind_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(np.array([glmSW[0],glmNE[0]]), np.array([glmSW[1], glmNE[1]]), polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            regNE = (rot_db_lon[1], rot_db_lat[1])
            regSW = (rot_db_lon[0], rot_db_lat[0])  
            
            #print(regSW, regNE)
            
            
            ## subset wind
            sub_xwind_0p5km = climatology.subset_cube_by_coord(xwind_0p5km,regSW, regNE)
            sub_ywind_0p5km = climatology.subset_cube_by_coord(ywind_0p5km,regSW, regNE)
            
            sub_xwind_1p5km = climatology.subset_cube_by_coord(xwind_1p5km,regSW, regNE)
            sub_ywind_1p5km = climatology.subset_cube_by_coord(ywind_1p5km,regSW, regNE)
            
            sub_xwind_4p4km = climatology.subset_cube_by_coord(xwind_4p4km,regSW, regNE)
            sub_ywind_4p4km = climatology.subset_cube_by_coord(ywind_4p4km,regSW, regNE)
            
            
            
            
            #%% regrid using iris (0p5km and 1p5km onto 4p4km grid and glm onto all reg)
            
            ## preemptively remove land points from glm data
            sub_xwind_glm.data[sub_lbm_glm.data[:,:-1]==1] = np.nan
            sub_ywind_glm.data[sub_lbm_glm.data[:,:-1]==1] = np.nan
            ## preemptively remove land points from 4p4km data
            #sub_xwind_0p5km.data[sub_lbm_0p5km.data==1] = np.nan
            #sub_ywind_0p5km.data[sub_lbm_0p5km.data==1] = np.nan
            #sub_xwind_1p5km.data[sub_lbm_1p5km.data==1] = np.nan
            #sub_ywind_1p5km.data[sub_lbm_1p5km.data==1] = np.nan
            sub_xwind_4p4km.data[sub_lbm_4p4km.data==1] = np.nan
            sub_ywind_4p4km.data[sub_lbm_4p4km.data==1] = np.nan
            
            regrid_scheme = iris.analysis.Linear(extrapolation_mode='mask')         #AreaWeighted(mdtol=0.5)
            regrid_4p4_xwind_0p5km = climatology.regrid_cubes(sub_xwind_4p4km, sub_lbm_0p5km, regrid_scheme)
            regrid_4p4_ywind_0p5km = climatology.regrid_cubes(sub_ywind_4p4km, sub_ywind_0p5km, regrid_scheme)
            regrid_4p4_xwind_1p5km = climatology.regrid_cubes(sub_xwind_4p4km, sub_xwind_1p5km, regrid_scheme)
            regrid_4p4_ywind_1p5km = climatology.regrid_cubes(sub_ywind_4p4km, sub_ywind_1p5km, regrid_scheme)
            nearest_scheme = iris.analysis.Linear(extrapolation_mode='mask')
            regrid_glm_xwind_0p5km = climatology.regrid_cubes(sub_xwind_glm, sub_xwind_0p5km,nearest_scheme)
            regrid_glm_xwind_1p5km = climatology.regrid_cubes(sub_xwind_glm, sub_xwind_1p5km,nearest_scheme)
            regrid_glm_xwind_4p4km = climatology.regrid_cubes(sub_xwind_glm, sub_xwind_4p4km,nearest_scheme)
            regrid_glm_ywind_0p5km = climatology.regrid_cubes(sub_ywind_glm, sub_ywind_0p5km,nearest_scheme)
            regrid_glm_ywind_1p5km = climatology.regrid_cubes(sub_ywind_glm, sub_ywind_1p5km,nearest_scheme)
            regrid_glm_ywind_4p4km = climatology.regrid_cubes(sub_ywind_glm, sub_ywind_4p4km,nearest_scheme)
            
            #%% calculate magnitudes, remove landpoints, perform subtraction
            ## regrid regional onto regional via area weighted
            regrid_4p4_wsp_0p5km = (regrid_4p4_xwind_0p5km.data**2 + regrid_4p4_ywind_0p5km.data[:-1]**2)**0.5
            regrid_4p4_wsp_1p5km = (regrid_4p4_xwind_1p5km.data**2 + regrid_4p4_ywind_1p5km.data[:-1]**2)**0.5
            nl_4p4_wsp_0p5km = regrid_4p4_wsp_0p5km.copy()
            nl_4p4_wsp_1p5km = regrid_4p4_wsp_1p5km.copy()
            nl_4p4_wsp_0p5km[sub_lbm_0p5km.data == 1] = np.nan
            nl_4p4_wsp_1p5km[sub_lbm_1p5km.data == 1] = np.nan
            diff_4p4_wsp_0p5km = nl_wsp_0p5km - nl_4p4_wsp_0p5km
            diff_4p4_wsp_1p5km = nl_wsp_1p5km - nl_4p4_wsp_1p5km
            
            
            regrid_glm_wsp_0p5km = (regrid_glm_xwind_0p5km.data**2+regrid_glm_xwind_0p5km.data**2)**0.5
            regrid_glm_wsp_1p5km = (regrid_glm_xwind_1p5km.data**2+regrid_glm_xwind_1p5km.data**2)**0.5
            regrid_glm_wsp_4p4km = (regrid_glm_xwind_4p4km.data**2+regrid_glm_xwind_4p4km.data**2)**0.5
            
            
            nl_glm_wsp_0p5km = regrid_glm_wsp_0p5km.copy()
            nl_glm_wsp_1p5km = regrid_glm_wsp_1p5km.copy()
            nl_glm_wsp_4p4km = regrid_glm_wsp_4p4km.copy()
            nl_glm_wsp_0p5km[sub_lbm_0p5km.data == 1] = np.nan
            nl_glm_wsp_1p5km[sub_lbm_1p5km.data == 1] = np.nan
            nl_glm_wsp_4p4km[sub_lbm_4p4km.data == 1] = np.nan
            
            diff_glm_wsp_0p5km = nl_wsp_0p5km - nl_glm_wsp_0p5km
            diff_glm_wsp_1p5km = nl_wsp_1p5km - nl_glm_wsp_1p5km
            diff_glm_wsp_4p4km = nl_wsp_4p4km - nl_glm_wsp_4p4km
            #%%
            print('MEANS',np.nanmean(nl_4p4_wsp_0p5km), np.nanmean(nl_wsp_0p5km))
            #%%
            print(np.nanmax(diff_4p4_wsp_0p5km))
            print(np.nanmin(diff_4p4_wsp_0p5km))
            print(np.nanmax(diff_4p4_wsp_1p5km))
            print(np.nanmin(diff_4p4_wsp_1p5km))
            #%%
            plot_twopanels = True
            if plot_twopanels:
                matplotlib.rcParams.update({'font.size': 20})
                pglmNE = (340.45,67.4)
                pglmSW = (337,66)
                
                polelat = xwind_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = xwind_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                pregNE = (rot_db_lon[1], rot_db_lat[1])
                pregSW = (rot_db_lon[0], rot_db_lat[0])  
                
                fig, (ax2,ax0) = plt.subplots(1,2, figsize = (18,10))
                gs = fig.add_gridspec(3,2)
                kwargs = {'vmin': -8.5,'vmax':8.5, 'cmap':'BrBG', 'shading':'auto'}
                lb_kwargs = {'vmin': 0,'vmax':1, 'cmap':'Greys', 'shading':'auto'}
                fig.suptitle(f'{glmtitle} 10m horizontal wind differences - 2018-03-19')
            
                 
                ax0.pcolormesh(lon_0p5km,lat_0p5km,diff_4p4_wsp_0p5km, **kwargs )
                ax0.set_xticks([])
                ax0.set_yticks([])
                ax0.set_xlim(left = pregSW[0], right = pregNE[0])
                ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax0.set_title('0p5km - 4p4km')
                
                ax2.pcolormesh(lon_0p5km, lat_0p5km, diff_glm_wsp_0p5km, **kwargs)
                ax2.set_xticks([])
                ax2.set_yticks([])
                ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax2.set_title('0p5km - GLM')
                
                div_norm = matplotlib.colors.Normalize(vmin = -8.5, vmax = 8.5)
                div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'BrBG')
                cbar = fig.colorbar(mappable=div_map,ax = (ax2,ax0), orientation = 'horizontal', fraction = 0.15, pad = 0.01 )
                cbar.ax.set(xlabel = 'K')
                plt.subplots_adjust(left=0.125,
                        bottom=0.25, 
                        right=0.9, 
                        top=0.85, 
                        wspace=0.05, 
                        hspace=0.35)
                
                #plt.tight_layout()
            
            #%%
            plot_allpanels = False
            if plot_allpanels:
                matplotlib.rcParams.update({'font.size': 20})
                pglmNE = (340.8,67.6)
                pglmSW = (337,66)
                
                polelat = xwind_0p5km.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = xwind_0p5km.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(np.array([pglmSW[0],pglmNE[0]]), np.array([pglmSW[1], pglmNE[1]]), polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                pregNE = (rot_db_lon[1], rot_db_lat[1])
                pregSW = (rot_db_lon[0], rot_db_lat[0])  
                
                fig = plt.figure(figsize = (12,18))
                gs = fig.add_gridspec(3,2)
                kwargs = {'vmin': -8.5,'vmax':8.5, 'cmap':'BrBG', 'shading':'auto'}
                lb_kwargs = {'vmin': 0,'vmax':1, 'cmap':'Greys', 'shading':'auto'}
                fig.suptitle(f'{glmtitle} avr. 10m horizontal wind differences - 2018-03-19')
                
                ax0 = fig.add_subplot(gs[0,1])
                ax0.pcolormesh(lon_0p5km,lat_0p5km,diff_4p4_wsp_0p5km, **kwargs )
                ax0.set_xticks([])
                ax0.set_yticks([])
                ax0.set_xlim(left = pregSW[0], right = pregNE[0])
                ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax0.set_title('0p5km - 4p4km')
                
                ax1 = fig.add_subplot(gs[1,1])
                ax1.pcolormesh(lon_1p5km,lat_1p5km,diff_4p4_wsp_1p5km, **kwargs )
                ax1.set_xticks([])
                ax1.set_yticks([])
                ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax1.set_title('1p5km - 4p4km')
                
                ax2 = fig.add_subplot(gs[0,0])
                ax2.pcolormesh(lon_0p5km, lat_0p5km, diff_glm_wsp_0p5km, **kwargs)
                ax2.set_xticks([])
                ax2.set_yticks([])
                ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax2.set_title('0p5km - GLM')
                
                ax3 = fig.add_subplot(gs[1,0])
                ax3.pcolormesh(lon_1p5km, lat_1p5km, diff_glm_wsp_1p5km, **kwargs)
                ax3.set_xticks([])
                ax3.set_yticks([])
                ax3.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax3.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax3.set_title('1p5km - GLM')
                
                ax4 = fig.add_subplot(gs[2,0])
                ax4.pcolormesh(lon_4p4km, lat_4p4km, diff_glm_wsp_4p4km, **kwargs)
                ax4.set_xticks([])
                ax4.set_yticks([])
                ax4.set_xlim(left = pregSW[0], right = pregNE[0] )
                ax4.set_ylim(bottom = pregSW[1], top = pregNE[1])
                ax4.set_title('4p4km - GLM')
                
                #ax5 = fig.add_subplot(gs[2,1])
                #ax5.pcolormesh(lon_0p5km, lat_0p5km, nl_wsp_0p5km-nl_4p4_wsp_0p5km, vmin = -15, vmax = 15, shading = 'auto', cmap = 'BrBG')
                
                div_norm = matplotlib.colors.Normalize(vmin = -8.5, vmax = 8.5)
                div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'BrBG')
                
                #fig.subplots_adjust(right = 0.8)
                cbar_ax = fig.add_axes([0.99,0.03, 0.05,0.93])
                cbar = fig.colorbar(mappable=div_map,cax = cbar_ax)
                cbar.ax.set(ylabel = r'ms$^{-1}$')
                
                plt.tight_layout()
                
                plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_surface_differences_wASCAT_3houravr_upper{glmtitle}_noland.png')
                
                plt.show()