# -*- coding: utf-8 -*-
"""
Created on Thu May 19 13:10:27 2022

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
glmhidx = 4

#%%

wbp_glm_15h = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'wet_bulb_potential_temperature')[0, glmhidx]
wbp_glm_18h = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'wet_bulb_potential_temperature')[1, glmhidx]
wbp_glm_21h = iris.load_cube(f'{um_datapath}RA1M_glm_4_flt306.nc', 'wet_bulb_potential_temperature')[0, glmhidx]

wbp_glm_15h = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'geopotential_height')[0, glmhidx]
wbp_glm_18h = iris.load_cube(f'{um_datapath}RA1M_glm_3_flt306.nc', 'geopotential_height')[1, glmhidx]
wbp_glm_21h = iris.load_cube(f'{um_datapath}RA1M_glm_4_flt306.nc', 'geopotential_height')[0, glmhidx
                                                                                          
#%%
                                                                                          ]
