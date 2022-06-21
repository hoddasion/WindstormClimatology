# -*- coding: utf-8 -*-
"""
Created on Tue May 10 17:23:19 2022

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
        #%%
        
        xwind_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'x_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        ywind_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'y_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        xwind_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'x_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        ywind_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'y_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        xwind_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'x_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        ywind_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'y_wind')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        
        xwind_glm = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'm01s15i201')#glm_cubes[30]
        ywind_glm = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'm01s15i202')#glm_cubes[32]
        
        gph_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'geopotential_height')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        gph_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'geopotential_height')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        gph_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'geopotential_height')[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
        gph_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'geopotential_height')
        
        lbm_0p5km = iris.load_cube(f'{um_datapath}RA1M_0p5km_cc134_24hrs_pressurevars_306.nc', 'land_binary_mask')
        lbm_1p5km = iris.load_cube(f'{um_datapath}RA1M_1p5km_cc134_24hrs_pressurevars_306.nc', 'land_binary_mask')
        lbm_4p4km = iris.load_cube(f'{um_datapath}RA1M_4p4km_cc134_24hrs_pressurevars_306.nc', 'land_binary_mask')
        lbm_glm   = iris.load_cube(f'{um_datapath}RA1M_glm_1_flt306.nc', 'land_binary_mask')
        #%%
        #glm_cubes = iris.load(f'{um_datapath}RA1M_glm_3_flt306.nc')
        #print(glm_cubes)
        #%%
        
        print(xwind_glm.coord('pressure').points)
        print(xwind_0p5km.coord('pressure').points)
        reg_plevels = xwind_0p5km.coord('pressure').points
        
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
        
        sub_xwind_glm = climatology.subset_cube_by_coord(xwind_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')[glmidx, glmhidx]
        sub_ywind_glm = climatology.subset_cube_by_coord(ywind_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')[glmidx, glmhidx]
        
        ## subset land-binary-mask
        sub_lbm_0p5km = climatology.subset_cube_by_coord(lbm_0p5km,regSW, regNE)
        sub_lbm_1p5km = climatology.subset_cube_by_coord(lbm_1p5km,regSW, regNE)
        sub_lbm_4p4km = climatology.subset_cube_by_coord(lbm_4p4km,regSW, regNE)
        sub_lbm_glm   = climatology.subset_cube_by_coord(lbm_glm, glmSW, glmNE,lonname = 'longitude', latname = 'latitude')
        
        ## subset msp
        sub_gph_0p5km = climatology.subset_cube_by_coord(gph_0p5km,regSW, regNE)
        sub_gph_1p5km = climatology.subset_cube_by_coord(gph_1p5km,regSW, regNE)
        sub_gph_4p4km = climatology.subset_cube_by_coord(gph_4p4km,regSW, regNE)
        sub_gph_glm = climatology.subset_cube_by_coord(gph_glm,glmSW, glmNE,lonname = 'longitude', latname = 'latitude')[glmidx, glmhidx]
        
        #print(sub_xwind_glm)
        
        lon_0p5km = sub_xwind_0p5km.coord('grid_longitude').points
        lat_0p5km = sub_xwind_0p5km.coord('grid_latitude').points
        
        lon_1p5km = sub_xwind_1p5km.coord('grid_longitude').points
        lat_1p5km = sub_xwind_1p5km.coord('grid_latitude').points
        
        lon_4p4km = sub_xwind_4p4km.coord('grid_longitude').points
        lat_4p4km = sub_xwind_4p4km.coord('grid_latitude').points
        
        lon_glm = sub_xwind_glm.coord('longitude').points
        lat_glm = sub_xwind_glm.coord('latitude').points
        ## compute magnitudes
        wsp_0p5km = climatology.quadmag(sub_xwind_0p5km.data[reghidx], sub_ywind_0p5km.data[reghidx])
        wsp_1p5km = climatology.quadmag(sub_xwind_1p5km.data[reghidx], sub_ywind_1p5km.data[reghidx])
        wsp_4p4km = climatology.quadmag(sub_xwind_4p4km.data[reghidx], sub_ywind_4p4km.data[reghidx])
        wsp_glm = climatology.quadmag(sub_xwind_glm.data, sub_ywind_glm.data)
            
        
        print(np.shape(lon_glm), np.shape(lat_glm), np.shape(wsp_glm), np.shape(sub_xwind_glm))
        print(lon_glm)
        #%% cursory plotting
        plot_wholeregion_egu = True
        if plot_wholeregion_egu:
            quiv0p5 = 30
            quiv1p5 = 10
            quiv4p4 = 4
            quivglm = 2
            quivsat = 3
            
            ## set gph levels for contour plots
            gph_levels = np.arange(700,2500, 10)
            
            ## compute magnitudes
            
            
            fig, (ax0,ax1,ax2) = plt.subplots(3,1, figsize = (8,20))
            kwargs = {'cmap' : 'Oranges', 'vmin' : 0, 'vmax' : np.nanmax(wsp_0p5km), 'shading' : 'auto'}
            mspkwargs = {'levels':gph_levels, 'colors':'k'}
            fig.suptitle(f'{titletime} avr. 850hPa \nhorizontal windspeed - 2018-03-19')
            
            ax0.pcolormesh(lon_0p5km, lat_0p5km, wsp_0p5km, **kwargs)
            ax0.quiver(lon_0p5km[::quiv0p5], lat_0p5km[::quiv0p5], sub_xwind_0p5km.data[reghidx,::quiv0p5,::quiv0p5], sub_ywind_0p5km.data[reghidx,::quiv0p5,::quiv0p5])
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.contour(sub_lbm_0p5km.coords('grid_longitude')[0].points, sub_lbm_0p5km.coords('grid_latitude')[0].points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax0.set_xlim( right = regNE[0])
            ax0.set_ylim( top = regNE[1])
            ax0.set_title('0p5km')
            c0p5= ax0.contour(lon_0p5km, lat_0p5km, sub_gph_0p5km.data[reghidx], **mspkwargs)
            
            ax1.pcolormesh(lon_4p4km, lat_4p4km, wsp_4p4km, **kwargs)
            ax1.quiver(lon_4p4km[::quiv4p4], lat_4p4km[::quiv4p4], sub_xwind_4p4km.data[reghidx,::quiv4p4,::quiv4p4], sub_ywind_4p4km.data[reghidx,::quiv4p4,::quiv4p4])
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.contour(sub_lbm_4p4km.coords('grid_longitude')[0].points, sub_lbm_4p4km.coords('grid_latitude')[0].points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            ax1.set_xlim( right = regNE[0])
            ax1.set_ylim( top = regNE[1])
            ax1.set_title('4p4km')
            c4p4= ax1.contour(lon_4p4km, lat_4p4km, sub_gph_4p4km.data[reghidx], **mspkwargs)
            
            ax2.pcolormesh(lon_glm, lat_glm, wsp_glm, **kwargs)
            ax2.quiver(lon_glm[::quivglm], lat_glm[::quivglm], sub_xwind_glm.data[::quivglm,::quivglm], sub_ywind_glm.data[::quivglm,::quivglm])
            ax2.contour(sub_lbm_glm.coord('longitude').points, sub_lbm_glm.coord('latitude').points, sub_lbm_glm.data, colors = 'green', levels = [1], linewidths = 4)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim( right = glmNE[0])
            ax2.set_ylim( top = glmNE[1])
            ax2.set_title(f'GLM - {glmtitle}')
            cglm = ax2.contour(sub_gph_glm.coords('longitude')[0].points, sub_gph_glm.coords('latitude')[0].points, sub_gph_glm.data, **mspkwargs)
            
            div_norm = matplotlib.colors.Normalize(vmin = 0, vmax = np.nanmax(wsp_0p5km))
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'Oranges')
            
            #fig.subplots_adjust(right = 0.8)
            #cbar_ax = fig.add_axes([0.03,0.99, 0.93,0.05])
            cbar = fig.colorbar(mappable=div_map,ax = (ax0,ax1,ax2), orientation = 'vertical', fraction = 0.15, pad = 0.04 )
            cbar.ax.set(ylabel = r'ms$^{-1}$')
            plt.subplots_adjust(left=0.125,
                    bottom=0.12, 
                    right=0.75, 
                    top=0.9, 
                    wspace=0.05, 
                    hspace=0.1)
            
            print('slay')
        
        #%%
        plot_wholeregion = False
        if plot_wholeregion:
            quiv0p5 = 30
            quiv1p5 = 10
            quiv4p4 = 4
            quivglm = 2
            quivsat = 3
            
            ## set gph levels for contour plots
            gph_levels = np.arange(700,2500, 10)
            
            ## compute magnitudes
            wsp_0p5km = climatology.quadmag(sub_xwind_0p5km.data[reghidx], sub_ywind_0p5km.data[reghidx])
            wsp_1p5km = climatology.quadmag(sub_xwind_1p5km.data[reghidx], sub_ywind_1p5km.data[reghidx])
            wsp_4p4km = climatology.quadmag(sub_xwind_4p4km.data[reghidx], sub_ywind_4p4km.data[reghidx])
            wsp_glm = climatology.quadmag(sub_xwind_glm.data[glmidx, glmhidx], sub_ywind_glm.data[glmidx, glmhidx])
            
            fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (12,12))
            kwargs = {'cmap' : 'Oranges', 'vmin' : 0, 'vmax' : np.nanmax(wsp_0p5km), 'shading' : 'auto'}
            mspkwargs = {'levels':gph_levels, 'colors':'k'}
            fig.suptitle(f'{titletime} avr. horizontal windspeed - 2018-03-19')
            
            ax0.pcolormesh(lon_0p5km, lat_0p5km, wsp_0p5km, **kwargs)
            ax0.quiver(lon_0p5km[::quiv0p5], lat_0p5km[::quiv0p5], sub_xwind_0p5km.data[reghidx,::quiv0p5,::quiv0p5], sub_ywind_0p5km.data[reghidx,::quiv0p5,::quiv0p5])
            ax0.contour(sub_lbm_0p5km.coords('grid_longitude')[0].points, sub_lbm_0p5km.coords('grid_latitude')[0].points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            c0p5= ax0.contour(lon_0p5km, lat_0p5km, sub_gph_0p5km.data[reghidx], **mspkwargs)
            ax0.clabel(c0p5, c0p5.levels, fontsize = 'smaller', fmt = fmt)
            ax0.set_title(f'0p5km - {reg_plevels[reghidx]:.0f}hPa')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim( right = regNE[0])
            ax0.set_ylim( top = regNE[1])
            
            ax1.pcolormesh(lon_1p5km, lat_1p5km, wsp_1p5km, **kwargs)
            ax1.quiver(lon_1p5km[::quiv1p5], lat_1p5km[::quiv1p5], sub_xwind_1p5km.data[reghidx,::quiv1p5,::quiv1p5], sub_ywind_1p5km.data[reghidx,::quiv1p5,::quiv1p5])
            ax1.contour(sub_lbm_1p5km.coords('grid_longitude')[0].points, sub_lbm_1p5km.coords('grid_latitude')[0].points, sub_lbm_1p5km.data, colors = 'green',levels = [1], linewidths = 4)
            c1p5= ax1.contour(lon_1p5km, lat_1p5km, sub_gph_1p5km.data[reghidx], **mspkwargs)
            ax1.clabel(c1p5, c1p5.levels, fontsize = 'smaller', fmt = fmt)
            ax1.set_title(f'1p5km - {reg_plevels[reghidx]:.0f}hPa')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim( right = regNE[0])
            ax1.set_ylim( top = regNE[1])
            
            ax2.pcolormesh(lon_4p4km, lat_4p4km, wsp_4p4km, **kwargs)
            ax2.quiver(lon_4p4km[::quiv4p4], lat_4p4km[::quiv4p4], sub_xwind_4p4km.data[reghidx,::quiv4p4,::quiv4p4], sub_ywind_4p4km.data[reghidx,::quiv4p4,::quiv4p4])
            ax2.contour(sub_lbm_4p4km.coords('grid_longitude')[0].points, sub_lbm_4p4km.coords('grid_latitude')[0].points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            c4p4= ax2.contour(lon_4p4km, lat_4p4km, sub_gph_4p4km.data[reghidx], **mspkwargs)
            ax2.clabel(c4p4, c4p4.levels, fontsize = 'smaller', fmt = fmt)
            ax2.set_title(f'4p4km - {reg_plevels[reghidx]:.0f}hPa')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim( right = regNE[0])
            ax2.set_ylim( top = regNE[1])
            
            ax3.pcolormesh(lon_glm, lat_glm, wsp_glm, **kwargs)
            ax3.quiver(lon_glm[::quivglm], lat_glm[::quivglm], sub_xwind_glm.data[glmidx,glmhidx,::quivglm,::quivglm], sub_ywind_glm.data[glmidx,glmhidx,::quivglm,::quivglm])
            cglm = ax3.contour(sub_gph_glm.coords('longitude')[0].points, sub_gph_glm.coords('latitude')[0].points, sub_gph_glm.data[glmidx,glmhidx], **mspkwargs)
            ax3.clabel(cglm, cglm.levels, fontsize = 'smaller', fmt = fmt)
            ax3.contour(sub_lbm_glm.coord('longitude').points, sub_lbm_glm.coord('latitude').points, sub_lbm_glm.data, colors = 'green', levels = [1], linewidths = 4)
            ax3.set_title(f'GLM - {glmtitle} - 850hPa')
            ax3.set_xticks([])
            ax3.set_yticks([])
            
            wsp_norm = matplotlib.colors.Normalize(vmin = 0, vmax = np.nanmax(wsp_0p5km))
            wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.86])
            cbar = fig.colorbar(mappable=wsp_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'$ms^{-1}$') 
            
            plt.tight_layout()
            
            #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_{reg_plevels[reghidx]:.0f}hPa_wind_3houravr_upper{glmtitle}.png')
        
        #%% set new subsetting coordinates
        ## using "p" to denote plotting purpose
        plot_wakeregion = False
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
            glmhidx = 4; reghidx = 12
            ## set gph levels for contour plots
            gph_levels = np.arange(700,2500, 5)
            
            ## compute magnitudes
            wsp_0p5km = climatology.quadmag(sub_xwind_0p5km.data[reghidx], sub_ywind_0p5km.data[reghidx])
            wsp_1p5km = climatology.quadmag(sub_xwind_1p5km.data[reghidx], sub_ywind_1p5km.data[reghidx])
            wsp_4p4km = climatology.quadmag(sub_xwind_4p4km.data[reghidx], sub_ywind_4p4km.data[reghidx])
            wsp_glm = climatology.quadmag(sub_xwind_glm.data[glmidx, glmhidx], sub_ywind_glm.data[glmidx, glmhidx])
            
            fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (12,12))
            kwargs = {'cmap' : 'Oranges', 'vmin' : 0, 'vmax' : np.nanmax(wsp_0p5km), 'shading' : 'auto'}
            mspkwargs = {'levels':gph_levels, 'colors':'k'}
            fig.suptitle(f'{titletime} avr. horizontal windspeed - 2018-03-19')
            
            ax0.pcolormesh(lon_0p5km, lat_0p5km, wsp_0p5km, **kwargs)
            ax0.quiver(lon_0p5km[::quiv0p5], lat_0p5km[::quiv0p5], sub_xwind_0p5km.data[reghidx,::quiv0p5,::quiv0p5], sub_ywind_0p5km.data[reghidx,::quiv0p5,::quiv0p5])
            ax0.contour(sub_lbm_0p5km.coords('grid_longitude')[0].points, sub_lbm_0p5km.coords('grid_latitude')[0].points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            c0p5= ax0.contour(lon_0p5km, lat_0p5km, sub_gph_0p5km.data[reghidx], **mspkwargs)
            ax0.clabel(c0p5, c0p5.levels, fontsize = 'smaller', fmt = fmt)
            ax0.set_title(f'0p5km - {reg_plevels[reghidx]:.0f}hPa')
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            
            ax1.pcolormesh(lon_1p5km, lat_1p5km, wsp_1p5km, **kwargs)
            ax1.quiver(lon_1p5km[::quiv1p5], lat_1p5km[::quiv1p5], sub_xwind_1p5km.data[reghidx,::quiv1p5,::quiv1p5], sub_ywind_1p5km.data[reghidx,::quiv1p5,::quiv1p5])
            ax1.contour(sub_lbm_1p5km.coords('grid_longitude')[0].points, sub_lbm_1p5km.coords('grid_latitude')[0].points, sub_lbm_1p5km.data, colors = 'green',levels = [1], linewidths = 4)
            c1p5= ax1.contour(lon_1p5km, lat_1p5km, sub_gph_1p5km.data[reghidx], **mspkwargs)
            ax1.clabel(c1p5, c1p5.levels, fontsize = 'smaller', fmt = fmt)
            ax1.set_title(f'1p5km - {reg_plevels[reghidx]:.0f}hPa')
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            
            ax2.pcolormesh(lon_4p4km, lat_4p4km, wsp_4p4km, **kwargs)
            ax2.quiver(lon_4p4km[::quiv4p4], lat_4p4km[::quiv4p4], sub_xwind_4p4km.data[reghidx,::quiv4p4,::quiv4p4], sub_ywind_4p4km.data[reghidx,::quiv4p4,::quiv4p4])
            ax2.contour(sub_lbm_4p4km.coords('grid_longitude')[0].points, sub_lbm_4p4km.coords('grid_latitude')[0].points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            c4p4= ax2.contour(lon_4p4km, lat_4p4km, sub_gph_4p4km.data[reghidx], **mspkwargs)
            ax2.clabel(c4p4, c4p4.levels, fontsize = 'smaller', fmt = fmt)
            ax2.set_title(f'4p4km - {reg_plevels[reghidx]:.0f}hPa')
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            
            ax3.pcolormesh(lon_glm, lat_glm, wsp_glm, **kwargs)
            ax3.quiver(lon_glm[::quivglm], lat_glm[::quivglm], sub_xwind_glm.data[glmidx,glmhidx,::quivglm,::quivglm], sub_ywind_glm.data[glmidx,glmhidx,::quivglm,::quivglm])
            cglm = ax3.contour(sub_gph_glm.coords('longitude')[0].points, sub_gph_glm.coords('latitude')[0].points, sub_gph_glm.data[glmidx,glmhidx], **mspkwargs)
            ax3.clabel(cglm, cglm.levels, fontsize = 'smaller', fmt = fmt)
            ax3.contour(sub_lbm_glm.coord('longitude').points, sub_lbm_glm.coord('latitude').points, sub_lbm_glm.data, colors = 'green', levels = [1], linewidths = 4)
            ax3.set_title(f'GLM - {glmtitle} - 850hPa')
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            
            wsp_norm = matplotlib.colors.Normalize(vmin = 0, vmax = np.nanmax(wsp_0p5km))
            wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
            fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.86])
            cbar = fig.colorbar(mappable=wsp_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'$ms^{-1}$') 
            
            plt.tight_layout()
            
            #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_{reg_plevels[reghidx]:.0f}hPa_win_3houravr_upper{glmtitle}_wakeregion.png')
        
        #%%
        horizontal_convergence = False
        if horizontal_convergence:
            import climatology
            lonmesh0p5km, latmesh0p5km = np.meshgrid(lon_0p5km, lat_0p5km)
            lonmesh1p5km, latmesh1p5km = np.meshgrid(lon_1p5km, lat_1p5km)
            lonmesh4p4km, latmesh4p4km = np.meshgrid(lon_4p4km, lat_4p4km)
            lonmeshglm, latmeshglm = np.meshgrid(lon_glm, lat_glm)
            div_0p5km = climatology.divergence_horizontal(sub_xwind_0p5km.data[reghidx], sub_ywind_0p5km.data[reghidx], lonmesh0p5km, latmesh0p5km)
            
            div_1p5km = climatology.divergence_horizontal(sub_xwind_1p5km.data[reghidx], sub_ywind_1p5km.data[reghidx], lonmesh1p5km, latmesh1p5km)
            div_4p4km = climatology.divergence_horizontal(sub_xwind_4p4km.data[reghidx], sub_ywind_4p4km.data[reghidx], lonmesh4p4km, latmesh4p4km)
            div_glm = climatology.divergence_horizontal(sub_xwind_glm.data[glmidx, glmhidx], sub_ywind_glm.data[glmidx, glmhidx], lonmeshglm, latmeshglm)
            
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
            
            
            fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (12,12))
            gs = fig.add_gridspec(3,2)
            kwargs = {'cmap' : 'RdYlBu_r',  'shading' : 'auto', 'vmin':-0.0025, 'vmax':0.0025}
            fig.suptitle(f'{titletime} avr. horizontal wind divergence - 2018-03-19')
            
            
            ax0.pcolormesh(lon_0p5km[:-1], lat_0p5km[:-1], div_0p5km, **kwargs)
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            
            ax0.set_title(f'0p5km - {reg_plevels[reghidx]:.0f}hPa')
            
            
            ax1.pcolormesh(lon_1p5km[:-1], lat_1p5km[:-1], div_1p5km, **kwargs)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            
            ax1.set_title(f'1p5km - {reg_plevels[reghidx]:.0f}hPa')
            
            
            ax2.pcolormesh(lon_4p4km[:-1], lat_4p4km[:-1], div_4p4km, **kwargs)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0] )
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            
            ax2.set_title(f'4p4km - {reg_plevels[reghidx]:.0f}hPa')
            
            
            ax3.pcolormesh(lon_glm[:-1], lat_glm[:-1], div_glm, **kwargs)
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pglmSW[0], right = pglmNE[0] )
            ax3.set_ylim(bottom = pglmSW[1], top = pglmNE[1])
            ax3.set_title(f'GLM - {glmtitle} - {reg_plevels[reghidx]:.0f}hPa')
            
            div_norm = matplotlib.colors.Normalize(vmin = -0.0025, vmax = 0.0025)
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'RdYlBu_r')
            
            #fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([1,0.03, 0.05,0.87])
            cbar = fig.colorbar(mappable=div_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'$s^{-1}$') 
            
            plt.tight_layout()
            
            #plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_850hPa_2Dwind_divergence_3houravr_upper{glmtitle}_wakeregion.png')
            
        #%%
        plot_differences = False
        if plot_differences:
            ## remove land points from glm cube before interpolation
            
            ## regrid using iris (0p5km and 1p5km onto 4p4km grid and glm onto all reg)
            regrid_scheme = iris.analysis.Linear(extrapolation_mode = 'mask')
            regrid_4p4_xwind_0p5km = climatology.regrid_cubes( sub_xwind_4p4km, sub_xwind_0p5km, regrid_scheme)[reghidx]
            regrid_4p4_ywind_0p5km = climatology.regrid_cubes( sub_ywind_4p4km, sub_ywind_0p5km, regrid_scheme)[reghidx]
            regrid_4p4_xwind_1p5km = climatology.regrid_cubes(sub_xwind_4p4km,sub_xwind_1p5km,  regrid_scheme)[reghidx]
            regrid_4p4_ywind_1p5km = climatology.regrid_cubes( sub_ywind_4p4km,sub_ywind_1p5km, regrid_scheme)[reghidx]
            nearest_scheme = iris.analysis.Linear(extrapolation_mode='mask')
            regrid_glm_xwind_0p5km = climatology.regrid_cubes(sub_xwind_glm, sub_xwind_0p5km,nearest_scheme)
            regrid_glm_xwind_1p5km = climatology.regrid_cubes(sub_xwind_glm, sub_xwind_1p5km,nearest_scheme)
            regrid_glm_xwind_4p4km = climatology.regrid_cubes(sub_xwind_glm, sub_xwind_4p4km,nearest_scheme)
            regrid_glm_ywind_0p5km = climatology.regrid_cubes(sub_ywind_glm, sub_ywind_0p5km,nearest_scheme)
            regrid_glm_ywind_1p5km = climatology.regrid_cubes(sub_ywind_glm, sub_ywind_1p5km,nearest_scheme)
            regrid_glm_ywind_4p4km = climatology.regrid_cubes(sub_ywind_glm, sub_ywind_4p4km,nearest_scheme)
           
            
            ## compute windspeed magntitudes
            regrid_4p4_wsp_0p5km = (regrid_4p4_xwind_0p5km.data**2 + regrid_4p4_ywind_0p5km.data**2)**0.5
            regrid_4p4_wsp_1p5km = (regrid_4p4_xwind_1p5km.data**2 + regrid_4p4_ywind_1p5km.data**2)**0.5
            
            glm_wsp_0p5km = (regrid_glm_xwind_0p5km.data**2+regrid_glm_xwind_0p5km.data**2)**0.5
            glm_wsp_1p5km = (regrid_glm_xwind_1p5km.data**2+regrid_glm_xwind_1p5km.data**2)**0.5
            glm_wsp_4p4km = (regrid_glm_xwind_4p4km.data**2+regrid_glm_xwind_4p4km.data**2)**0.5
            
            ## compute differences
            diff_4p4_wsp_0p5km = wsp_0p5km-  regrid_4p4_wsp_0p5km
            diff_4p4_wsp_1p5km = wsp_1p5km - regrid_4p4_wsp_1p5km
            
            diff_glm_wsp_0p5km = wsp_0p5km - glm_wsp_0p5km
            diff_glm_wsp_1p5km = wsp_1p5km - glm_wsp_1p5km
            diff_glm_wsp_4p4km = wsp_4p4km - glm_wsp_4p4km
            print('Slay')
            
            #%%
            print(np.nanmax(diff_glm_wsp_0p5km))
            print(np.nanmin(diff_glm_wsp_0p5km))
            print(np.nanmax(diff_glm_wsp_1p5km))
            print(np.nanmin(diff_glm_wsp_1p5km))
            
            #%%
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
            kwargs = {'vmin': -16,'vmax':16, 'cmap':'BrBG', 'shading':'auto'}
            lb_kwargs = {'vmin': 0,'vmax':1, 'cmap':'Greys', 'shading':'auto'}
            fig.suptitle(f'{titletime} avr. 850hPa horizontal wind differences - 2018-03-19')
            
            ax0 = fig.add_subplot(gs[0,1])
            ax0.pcolormesh(lon_0p5km,lat_0p5km,diff_4p4_wsp_0p5km.data, **kwargs )
            ax0.contour(sub_lbm_0p5km.coord(axis='x').points, sub_lbm_0p5km.coord(axis='y').points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax0.set_xticks([])
            ax0.set_yticks([])
            ax0.set_xlim(left = pregSW[0], right = pregNE[0])
            ax0.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax0.set_title('0p5km - 4p4km')
            
            ax1 = fig.add_subplot(gs[1,1])
            ax1.pcolormesh(lon_1p5km,lat_1p5km,diff_4p4_wsp_1p5km.data, **kwargs )
            ax1.contour(sub_lbm_1p5km.coord(axis='x').points, sub_lbm_1p5km.coord(axis='y').points,sub_lbm_1p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax1.set_xticks([])
            ax1.set_yticks([])
            ax1.set_xlim(left = pregSW[0], right = pregNE[0])
            ax1.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax1.set_title('1p5km - 4p4km')
            
            ax2 = fig.add_subplot(gs[0,0])
            ax2.pcolormesh(lon_0p5km,lat_0p5km,diff_glm_wsp_0p5km.data, **kwargs )
            ax2.contour(sub_lbm_0p5km.coord(axis='x').points, sub_lbm_0p5km.coord(axis='y').points, sub_lbm_0p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.set_xlim(left = pregSW[0], right = pregNE[0])
            ax2.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax2.set_title('0p5km - GLM')
            
            ax3 = fig.add_subplot(gs[1,0])
            ax3.pcolormesh(lon_1p5km,lat_1p5km,diff_glm_wsp_1p5km.data, **kwargs )
            ax3.contour(sub_lbm_1p5km.coord(axis='x').points, sub_lbm_1p5km.coord(axis='y').points, sub_lbm_1p5km.data, colors = 'green',levels = [1], linewidths = 4)
            ax3.set_xticks([])
            ax3.set_yticks([])
            ax3.set_xlim(left = pregSW[0], right = pregNE[0])
            ax3.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax3.set_title('1p5km - GLM')
            
            ax4 = fig.add_subplot(gs[2,0])
            ax4.pcolormesh(lon_4p4km,lat_4p4km,diff_glm_wsp_4p4km.data, **kwargs )
            ax4.contour(sub_lbm_4p4km.coord(axis='x').points, sub_lbm_4p4km.coord(axis='y').points, sub_lbm_4p4km.data, colors = 'green',levels = [1], linewidths = 4)
            ax4.set_xticks([])
            ax4.set_yticks([])
            ax4.set_xlim(left = pregSW[0], right = pregNE[0])
            ax4.set_ylim(bottom = pregSW[1], top = pregNE[1])
            ax4.set_title('4p4km - GLM')
            
            div_norm = matplotlib.colors.Normalize(vmin = -16, vmax = 16)
            div_map = matplotlib.cm.ScalarMappable(norm = div_norm, cmap = 'BrBG')
            
            #fig.subplots_adjust(right = 0.8)
            cbar_ax = fig.add_axes([0.99,0.03, 0.05,0.93])
            cbar = fig.colorbar(mappable=div_map,cax = cbar_ax)
            cbar.ax.set(ylabel = r'ms$^{-1}$')
            
            plt.tight_layout()
            
            plt.savefig(f'D:/Project/Climatology/Figures/Windstorm/RA1M_UMGLM_850hPa_differences_windspeed_3houravr_upper{glmtitle}.png')