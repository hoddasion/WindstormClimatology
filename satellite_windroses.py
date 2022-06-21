# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 12:20:04 2022

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
import windrose
from windrose import WindroseAxes
import pandas as pd
import matplotlib.cm as cm

matplotlib.rcParams.update({'font.size': 22})



#%%
sample_lons = np.array([338.35,336,335.8, 340])
sample_lats = np.array([66.6,67,65.2, 67.1])
sample_labels = ['A1', 'B', 'C', 'A2']





#%%
seasons = ['MAM','JJA','SON','DJF']

for season in seasons:
    ##load and subset wsp and wdir data by season
    da_wsp = clima.satellite_seasonal_grouping('wind_speed', season)
    da_dir = clima.satellite_seasonal_grouping('wind_to_dir', season)
    ## use nearest neighbour method for taking sample time series
    dir_sample = clima.sample_in_space(da_dir, sample_lons, sample_lats)
    wsp_sample = clima.sample_in_space(da_wsp, sample_lons, sample_lats)
    ## convert to pandas dataframes and remove rows containing nans
    df_A = pd.DataFrame({'speed':wsp_sample[0], 'direction':dir_sample[0]}).dropna()
    df_B = pd.DataFrame({'speed':wsp_sample[1], 'direction':dir_sample[1]}).dropna()
    df_C = pd.DataFrame({'speed':wsp_sample[2], 'direction':dir_sample[2]}).dropna()
    df_A2 = pd.DataFrame({'speed':wsp_sample[3], 'direction':dir_sample[3]}).dropna()
    
    #%%
    df_A.direction = df_A.direction - 180
    condA = np.array(df_A.direction < 0)
    df_A.direction[condA] = df_A.direction[condA] + 360
    
    df_A2.direction = df_A2.direction - 180
    condA2 = np.array(df_A2.direction < 0)
    df_A2.direction[condA2] = df_A2.direction[condA2] + 360
    
    df_B.direction = df_B.direction - 180
    condB = np.array(df_B.direction < 0)
    df_B.direction[condB] = df_B.direction[condB] + 360
    
    df_C.direction = df_C.direction - 180
    condC = np.array(df_C.direction < 0)
    df_C.direction[condC] = df_C.direction[condC] + 360
    
    #%% plot
    minlon, maxlon, minlat, maxlat = (-25, -19, 64, 69)

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(12, 10))
    # Draw main ax on top of which we will add windroses
    main_ax = fig.add_subplot(1, 1, 1, projection=proj)
    main_ax.set_extent([minlon, maxlon, minlat, maxlat], crs=proj)
    main_ax.gridlines(draw_labels=True)
    main_ax.coastlines()
    
    # Inset axe it with a fixed size
    axes_dim = 2.8
    
    wrax_A = inset_axes(main_ax,
            width=axes_dim,                             # size in inches
            height=axes_dim,                            # size in inches
            loc='center',                        # center bbox at given position
            bbox_to_anchor=(sample_lons[0]-360, sample_lats[0]), # position of the axe
            bbox_transform=main_ax.transData,    # use data coordinate (not axe coordinate)
            axes_class=windrose.WindroseAxes,    # specify the class of the axe
            )
    wrax_B = inset_axes(main_ax,
            width=axes_dim,                             # size in inches
            height=axes_dim,                            # size in inches
            loc='center',                        # center bbox at given position
            bbox_to_anchor=(sample_lons[1]-360, sample_lats[1]), # position of the axe
            bbox_transform=main_ax.transData,    # use data coordinate (not axe coordinate)
            axes_class=windrose.WindroseAxes,    # specify the class of the axe
            )
    wrax_C = inset_axes(main_ax,
            width=axes_dim,                             # size in inches
            height=axes_dim,                            # size in inches
            loc='center',                        # center bbox at given position
            bbox_to_anchor=(sample_lons[2]-360, sample_lats[2]), # position of the axe
            bbox_transform=main_ax.transData,    # use data coordinate (not axe coordinate)
            axes_class=windrose.WindroseAxes,    # specify the class of the axe
            )
    wrax_A2 = inset_axes(main_ax,
            width=axes_dim,                             # size in inches
            height=axes_dim,                            # size in inches
            loc='center',                        # center bbox at given position
            bbox_to_anchor=(sample_lons[3]-360, sample_lats[3]), # position of the axe
            bbox_transform=main_ax.transData,    # use data coordinate (not axe coordinate)
            axes_class=windrose.WindroseAxes,    # specify the class of the axe
            )
    for ax in [wrax_A, wrax_B, wrax_C, wrax_A2]:
        ax.tick_params(labelleft=False, labelbottom=False)
    ## now add date to windrose axes
    kwargs = {'bins':np.arange(0,24,4), 'cmap':cm.Oranges} # use this to have equal scaling/binning across all
    wrax_A.bar(np.array(df_A.direction),np.array(df_A.speed), **kwargs)
    wrax_B.bar(np.array(df_B.direction),np.array(df_B.speed), **kwargs)
    wrax_C.bar(np.array(df_C.direction),np.array(df_C.speed), **kwargs)
    wrax_A2.bar(np.array(df_A2.direction),np.array(df_A2.speed), **kwargs)
    ## insert and position legend
    wrax_A.legend(bbox_to_anchor=(0.8, -0.9 ))
    
    ## title
    fig.suptitle(f'Windroses for {season} 2010-2020')
    ## labelling
    
    for i, txt in enumerate(sample_labels):
        if txt == 'C':
            main_ax.annotate(txt, (sample_lons[i]+0.8 -360, sample_lats[i]-0.7))
        else:
            main_ax.annotate(txt, (sample_lons[i] -360, sample_lats[i]+1))
    
    ## save figure
    plt.savefig(f'D:/Project/Climatology/Figures/Windroses/ASCATA_windroses_on_map_season{season}_20102020_convention.png')