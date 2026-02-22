#!/usr/bin/env python

import getopt, sys, os
import re

import datetime

import netCDF4 as nc4

import pyart
import numpy as np
import numpy.ma as ma
import shutil
import glob
import gzip

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import cmocean
import getpass, socket

import pandas as pd

import cftime

version = 0.1

from pathlib import Path
homepath = Path.home()

matplotlib.use('Agg')


import cartopy

import geopandas as gpd
import fiona
import pyart
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import contextily as ctx
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import pyproj

from pyproj import Geod

from contextily.tile import warp_img_transform, warp_tiles, _warper


# Enable KML reading
fiona.drvsupport.supported_drivers['KML'] = 'rw'

try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:bk", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')

tracking_tag = 'AMOF_20220922221548';

campaign = 'woest';


yaml_project_file = os.path.join(homepath,'amof_campaigns',f'{campaign}_project.yml')
yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

blflag = False;
darkmode = False;

for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    elif o == "-b":
        blflag = True;
    elif o == "-k":
        darkmode = True;
    else:
        assert False, "unhandled option"

yr = datestr[0:4];
mo = datestr[4:6];
dy = datestr[6:8];

ncas_radar_vol1_path = '/gws/nopw/j04/ncas_radar_vol1';
ncas_obs_vol2_path = '/gws/pw/j07/ncas_obs_vol2/'
#inpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1',campaign,'L1',datestr,'sop');
inpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1','20230602_woest','v1.0.1','level1','sop',yr,mo,dy);

#inpath = os.path.join(ncas_radar_vol1_path,'cjw','projects',campaign,'camra','L1b');

if darkmode:
    qlstr = 'quicklooks_dark'
else:
    qlstr = 'quicklooks'

figpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1',f'woest_{qlstr}','v1.0.1','sop');
#figpath = os.path.join(inpath,'quicklooks','iop')

def make_woest_sop_rhi_plot(ncfile,figpath,blflag=False,darkmode=False):

    if darkmode:
        plt.style.use('dark_background')

    if blflag:
        hmax = 4;
        xmin = 0;
        xmax = 20;
    else:
        hmax = 12;
        xmin = 0;
        xmax = 90;
    
    dbz_cmap = 'pyart_HomeyerRainbow';
    vel_cmap = 'pyart_balance';
    ldr_cmap = 'pyart_SpectralExtended';
    #ldr_cmap = 'pyart_ChaseSpectral';
    spw_cmap = 'pyart_SpectralExtended';
    #spw_cmap = 'pyart_ChaseSpectral';
    #spw_cmap = 'pyart_NWS_SPW';

    from matplotlib import colors

    RadarDS = pyart.io.read_cfradial(ncfile,delay_field_loading=True);

    dtime0 = cftime.num2pydate(RadarDS.time['data'][0],RadarDS.time['units']);
    dtime0_str = dtime0.strftime("%Y%m%d-%H%M%S");
    nsweeps = RadarDS.nsweeps;

    vel_field = RadarDS.fields['VEL_HV']

    vel_limit_lower = vel_field['field_limit_lower']
    vel_limit_upper = vel_field['field_limit_upper']

    #vel_limit_lower = -15.0
    #vel_limit_upper = 15.0
#    fig, ax = plt.subplots(nsweeps,4,figsize=(24,nsweeps*4),constrained_layout=True)
#    fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)

    figpath = os.path.join(figpath,'rhi',datestr);
    if not os.path.isdir(figpath): 
        os.makedirs(figpath);

    #figpath_region = os.path.join(figpath,f'region_{region}');
    #if not os.path.isdir(figpath_region):
    #    os.makedirs(figpath_region);


    if nsweeps>0:

        for s in range(nsweeps):

            Radar = RadarDS.extract_sweeps([s])

            display = pyart.graph.RadarDisplay(Radar);

            gatefilter = pyart.correct.GateFilter(Radar)
            #gatefilter.exclude_below('SNR', -8.5)
            #gatefilter = pyart.correct.despeckle_field(Radar, "SNR", threshold=-50, size=4, gatefilter=gatefilter)

            print(f"sweep {s}/{nsweeps}");
            rhi_az = Radar.get_azimuth(0)[0];
            #if 180.<rhi_az<360.:
            #    ax_reverse = True
            #else:
            #    ax_reverse = False
            ax_reverse = False

            fig, ax = plt.subplots(4,1,figsize=(15,18),constrained_layout=True)
    
            fig.set_constrained_layout_pads(w_pad=1 / 72, h_pad=1 / 72, hspace=0.2,wspace=0.2)
            display.plot_rhi("DBZ_H", ax=ax[0], sweep=0, vmin=-20, vmax=50, gatefilter=gatefilter,
                             cmap=dbz_cmap, colorbar_orient='horizontal'); #,reverse_xaxis=ax_reverse);         
            ax[0].set_ylim(0,hmax)
            ax[0].set_xlim(xmin,xmax)
            ax[0].grid(True)
            #if ax_reverse:
            #    ax[0,0].invert_xaxis();
            ax[0].set_aspect('equal');
#            display.plot_colorbar(ax=ax[0],shrink=0.8);
            display.plot_rhi("VEL_HV", ax=ax[1], sweep=0, vmin=vel_limit_lower, vmax=vel_limit_upper, gatefilter=gatefilter,
                             cmap=vel_cmap, colorbar_orient='horizontal'); #,reverse_xaxis=ax_reverse)
            ax[1].set_ylim(0,hmax)
            ax[1].set_xlim(xmin,xmax)
            ax[1].grid(True)
            #if ax_reverse:
            #    ax[1,0].invert_xaxis();
            ax[1].set_aspect('equal','box');
            display.plot_rhi("SPW_HV", ax=ax[2], sweep=0, 
                             norm=colors.LogNorm(vmin=1e-1*np.sqrt(1e-1),vmax=np.sqrt(1e1)), gatefilter=gatefilter,
                             cmap=spw_cmap, colorbar_orient='horizontal'); #,reverse_xaxis=ax_reverse)
            ax[2].set_ylim(0,hmax)
            ax[2].set_xlim(xmin,xmax)
            ax[2].grid(True)
            #if ax_reverse:
            #    ax[1,1].invert_xaxis();
            
            ax[2].set_aspect('equal','box');
            display.plot_rhi("LDR", ax=ax[3], sweep=0, vmin=-35, vmax=5, gatefilter=gatefilter,
                             cmap=ldr_cmap, colorbar_orient='horizontal'); #,reverse_xaxis=ax_reverse)
            ax[3].set_ylim(0,hmax)
            ax[3].set_xlim(xmin,xmax)
            ax[3].grid(True)
            #if ax_reverse:
            #    ax[0,1].invert_xaxis();
            ax[3].set_aspect('equal','box');
            #sweep_start_index = Radar.get_start(0);
            dtime_sweep = cftime.num2pydate(Radar.time['data'][0],Radar.time['units']);
            dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
            figname = f'ncas-radar-camra-1_cao_{dtime_sweep_str}_rhi_az{rhi_az:0.2f}_l1_v1.0.1.png';
            plt.savefig(os.path.join(figpath,figname),dpi=300);
            plt.close();


def make_woest_sop_ppi_plot(radar, save_path, sweep=0, vmin_vmax=None):

    ppi_el = radar.get_elevation(0)[0];

    #if radar.azimuth['data'][0]>radar.azimuth['data'][1]:
    #    print(f'Sweep {sweep} Anticlockwise PPI');
    #    radar.azimuth['data'] += (-1.0)*(0.07); 
    #else:
    #    print(f'Sweep {sweep} Clockwise PPI');
    #    radar.azimuth['data'] += 1.0*(0.07); 

    gatefilter = pyart.filters.GateFilter(radar)

    display = pyart.graph.RadarDisplay(radar)

    # Define color maps for different radar variables
    colormaps = {
        "DBZ_H": "pyart_HomeyerRainbow",
        "VEL_HV": "pyart_balance",
        "ZDR": "pyart_HomeyerRainbow",
        "LDR": "viridis"
    }

    VEL_HV_min = radar.fields['VEL_HV']['field_limit_lower'];
    VEL_HV_max = radar.fields['VEL_HV']['field_limit_upper'];
    
    # Default vmin and vmax values if not provided
    default_vmin_vmax = {
        "DBZ_H": (-20, 50),
        "VEL_HV": (VEL_HV_min, VEL_HV_max),
        "ZDR": (-5, 5),
        "LDR": (-35, 5)
    }

    # Create a 2x2 subplot with correct projection
    fig, ax = plt.subplots(2, 2, figsize=(30, 30))

    # Plot PPI maps for different radar fields
    variables = ["DBZ_H", "VEL_HV", "ZDR", "LDR"]
    positions = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Define RHI azimuths for radial lines
    rhi_azimuths = [0.0, 25.0, 44.0, 60.0, 80.0, 99.0, 120.0, 143.0, 168.0, 
                    189.0, 215.0, 232.0, 244.0, 255.0, 270.0, 283.0, 293.0, 
                    310.0, 331.0, 348.0]

    # Maximum radar range in km
    #max_range_km = np.max(radar.range['data']) / 1000.0  
    max_range_km = 200.0;

    for var, pos in zip(variables, positions):

        ax[pos].set_aspect('equal','box')
        ax[pos].set_xlim(-200,200)
        ax[pos].set_ylim(-200,200)

        display.plot_ppi(
            var, 0, 
            vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1],
            fig=fig,
            ax=ax[pos],
            cmap=colormaps[var],
            colorbar_flag=False,
            alpha=1.0,
            gatefilter=gatefilter,
            edges=True
        )


        display.plot_colorbar(ax=ax[pos], orient='horizontal', shrink=0.8, cmap=colormaps[var])
        # Plot range rings
        #display.plot_range_rings([50, 100, 150, 200], ax=ax[pos], col='0.5', lw=0.5)
        ax[pos].plot([-1,1],[0,0],'k', lw=0.5)
        ax[pos].plot([0,0],[-1,1],'k', lw=0.5)

        # Define RHI azimuths for radial lines
        rhi_azimuths = [0.0, 25.0, 44.0, 60.0, 80.0, 99.0, 120.0, 143.0, 168.0, 189.0, 
                        215.0, 232.0, 244.0, 255.0, 270.0, 283.0, 293.0, 310.0, 331.0, 348.0]


        range_km = 193.0;

        for azimuth in rhi_azimuths:
            azimuth_rad = np.deg2rad(azimuth)
            xend = 200.0 * np.sin(azimuth_rad)
            yend = 200.0 * np.cos(azimuth_rad)
            ax[pos].plot([0,xend],[0,yend],color='grey', lw=0.5)

            # Determine text rotation (aligning with the azimuth angle)
            rotation_angle = -azimuth if azimuth <= 180 else 360 - azimuth  # Ensures text is aligned correctly

            # Add labels at the end of each RHI line
            ax[pos].text(xend*range_km/200., yend*range_km/200., f"{int(azimuth)}°", 
                    fontsize=10, color='black', fontweight='bold', 
                    ha='center', va='center', rotation=rotation_angle,  # Apply rotation
                    bbox=dict(facecolor='white', edgecolor='grey', alpha=1.0, boxstyle='round,pad=0.3'))

    print(save_path);
    # Save figure to file if save_path is provided
    dtime_sweep = cftime.num2pydate(radar.time['data'][0],radar.time['units']);
    dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
    figname = f'ncas-radar-camra-1_cao_{dtime_sweep_str}_ppi_el{ppi_el:0.2f}_l1_v1.0.1.png';
    plt.savefig(os.path.join(save_path,figname),dpi=300, bbox_inches='tight');
    plt.close();


def reverse_radar_azimuths(radar):
    """
    Creates a copy of a Py-ART radar object with azimuth angles reversed.
    
    Parameters:
    radar (pyart.core.Radar): The original Py-ART radar object.

    Returns:
    pyart.core.Radar: A new radar object with reversed azimuths.
    """
    # Deep copy to avoid modifying the original radar object
    radar_reversed = copy.deepcopy(radar)

    # Reverse azimuth array
    radar_reversed.azimuth["data"] = radar.azimuth["data"][::-1]

    # Reverse data fields to maintain alignment with the azimuths
    for field_name in radar.fields:
        radar_reversed.fields[field_name]["data"] = radar.fields[field_name]["data"][::-1]

    return radar_reversed

import pyart
import numpy as np
import copy

def ensure_decreasing_azimuths(radar):
    """
    Checks if azimuths are increasing and reverses them if necessary.

    Parameters:
    radar (pyart.core.Radar): The original Py-ART radar object.

    Returns:
    pyart.core.Radar: A new radar object with azimuths in decreasing order.
    """
    # Deep copy to avoid modifying the original radar object
    radar_adjusted = copy.deepcopy(radar)

    # Check if azimuths are increasing
    #if np.all(np.diff(radar.azimuth["data"]) > 0):  
    if radar.azimuth["data"][0] < radar.azimuth["data"][1]:
        print("Azimuths are increasing, reversing them...")

        # Reverse azimuth array
        radar_adjusted.azimuth["data"] = radar.azimuth["data"][::-1]

        # Reverse all data fields to maintain alignment
        for field_name in radar.fields:
            radar_adjusted.fields[field_name]["data"] = radar.fields[field_name]["data"][::-1]

    else:
        print("Azimuths are already in decreasing order. No changes made.")

    return radar_adjusted

def ensure_increasing_azimuths(radar):
    """
    Checks if azimuths are increasing and reverses them if necessary.

    Parameters:
    radar (pyart.core.Radar): The original Py-ART radar object.

    Returns:
    pyart.core.Radar: A new radar object with azimuths in decreasing order.
    """
    # Deep copy to avoid modifying the original radar object
    radar_adjusted = copy.deepcopy(radar)

    # Check if azimuths are decreasing
    #if np.all(np.diff(radar.azimuth["data"]) > 0):  
    if radar.azimuth["data"][1] < radar.azimuth["data"][0]:
        print("Azimuths are decreasing, reversing them...")

        # Reverse azimuth array
        radar_adjusted.azimuth["data"] = radar.azimuth["data"][::-1]

        # Reverse all data fields to maintain alignment
        for field_name in radar.fields:
            radar_adjusted.fields[field_name]["data"] = radar.fields[field_name]["data"][::-1]

    else:
        print("Azimuths are already in increasing order. No changes made.")

    return radar_adjusted

def make_woest_sop_ppi_map_plot(radar, os_key, kml_paths, save_path, sweep=0, vmin_vmax=None):
    """
    Plots a 2x2 PPI map of radar variables using OpenStreetMap background.

    Parameters:
    radar_file (str): Path to radar file (e.g., .nc or .h5 format).
    os_key (str): API key for Ordnance Survey (optional, if using OS maps).
    kml_paths (dict): Dictionary containing KML file paths for geospatial overlays.
                      Expected keys: 'verticals', 'horizontals', 'outline', 'surface_sites'
    azimuth (int): Azimuth angle for RHI analysis (default is 0).
    """

    # Load Radar Data
    #radar = pyart.io.read(radar_file)

    # Load KML Files into GeoDataFrames
    #gdf_list = {}
    #for key, path in kml_paths.items():
    #    gdf_list[key] = gpd.read_file(path, driver='LIBKML')

    ppi_el = radar.get_elevation(0)[0];

    # Setup OpenStreetMap tiles
    #osm_tiles = cimgt.OSM()
    #projection = osm_tiles.crs  # Use correct map projection
    tiles = cimgt.GoogleTiles(style="satellite")  # Use "roadmap" for standard map

    lat_0=radar.latitude["data"][0];
    lon_0=radar.longitude["data"][0];

    # Setting projection and ploting the second tilt
    #lambert = ccrs.LambertConformal(
    #    central_latitude=lat_0,
    #    central_longitude=lon_0
    #)

    #projection = ccrs.Mercator()
    #projection = cimgt.OSM().crs
    projection = ccrs.AzimuthalEquidistant(central_longitude=lon_0, central_latitude=lat_0)


    # Setup PyART gate filtering
    gatefilter = pyart.filters.GateFilter(radar)
    #gatefilter.exclude_below("SNR", -4)  # Adjust filtering criteria

    # Radar display setup
    display = pyart.graph.RadarMapDisplay(radar)

    # Define color maps for different radar variables
    colormaps = {
        "DBZ_H": "pyart_HomeyerRainbow",
        "VEL_HV": "pyart_balance",
        "ZDR": "pyart_HomeyerRainbow",
        "LDR": "viridis"
    }

    VEL_HV_min = radar.fields['VEL_HV']['field_limit_lower'];
    VEL_HV_max = radar.fields['VEL_HV']['field_limit_upper'];
    
    # Default vmin and vmax values if not provided
    default_vmin_vmax = {
        "DBZ_H": (-20, 50),
        "VEL_HV": (VEL_HV_min, VEL_HV_max),
        "ZDR": (-5, 5),
        "LDR": (-35, 5)
    }

    # Create a 2x2 subplot with correct projection
    fig, ax = plt.subplots(2, 2, figsize=(30, 30), subplot_kw=dict(projection=projection))

    # Define radar bounds
    #lat_min = np.min(radar.gate_latitude['data'])
    #lon_min = np.min(radar.gate_longitude['data'])
    #lat_max = np.max(radar.gate_latitude['data'])
    #lon_max = np.max(radar.gate_longitude['data'])

    lat_min = 49.376
    lon_min = -4.2875
    lat_max = 52.9337 
    lon_max = 1.3892

    #lat_min = 51.14
    #lon_min = -1.45
    #lat_max = 51.152
    #lon_max = -1.43

    #lat_min = 51.115
    #lon_min = -1.55
    #lat_max = 51.145
    #lon_max = -1.48

    # Plot PPI maps for different radar fields
    variables = ["DBZ_H", "VEL_HV", "ZDR", "LDR"]
    positions = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Define RHI azimuths for radial lines
    rhi_azimuths = [0.0, 25.0, 44.0, 60.0, 80.0, 99.0, 120.0, 143.0, 168.0, 
                    189.0, 215.0, 232.0, 244.0, 255.0, 270.0, 283.0, 293.0, 
                    310.0, 331.0, 348.0]

    # Maximum radar range in km
    #max_range_km = np.max(radar.range['data']) / 1000.0  
    max_range_km = 200.0;
    
    # Get radar location
    radar_lon, radar_lat = radar.longitude['data'][0], radar.latitude['data'][0]

    # Define the geodetic projection (WGS84 Ellipsoid)
    geod = Geod(ellps="WGS84")
    
    for var, pos in zip(variables, positions):
        display.plot_ppi_map(
            var, 0, 
            vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1],
            min_lat=40.0, max_lat=53.0,
            min_lon=-4.0, max_lon=2.0,
            lon_lines=None,
            #np.arange(-4.0, 2.0, 0.5),
            lat_lines=None,
            #np.arange(40.0, 53.0, 0.5),
            projection=projection,
            fig=fig,
            ax=ax[pos],
            lat_0=lat_0,
            lon_0=lon_0,
            cmap=colormaps[var],
            colorbar_flag=False,
            alpha=0.6,
            gatefilter=gatefilter,
            edges=True,
            embellish=False,
            resolution="10m",
        )
        display.plot_colorbar(ax=ax[pos], orient='horizontal', shrink=0.8, cmap=colormaps[var])
        # Plot range rings
        display.plot_range_rings([50, 100, 150, 200], '0.5', lw=0.5)
        #display.plot_line_xy([-1000,1000],[0,0],'k', lw=0.5)
        #display.plot_line_xy([0,0],[-1000,1000],'k', lw=0.5)

        # Define RHI azimuths for radial lines
        rhi_azimuths = [0.0, 25.0, 44.0, 60.0, 80.0, 99.0, 120.0, 143.0, 168.0, 189.0, 
                        215.0, 232.0, 244.0, 255.0, 270.0, 283.0, 293.0, 310.0, 331.0, 348.0]


        range_m = 193000.0;


        for azimuth in rhi_azimuths:
            azimuth_rad = np.deg2rad(azimuth)
            xend = 200000.0 * np.sin(azimuth_rad)
            yend = 200000.0 * np.cos(azimuth_rad)
            display.plot_line_xy([0,xend],[0,yend],color='grey', lw=0.5)

            # Compute new lat/lon using forward geodetic transformation
            end_lon, end_lat = geod.fwd(lon_0, lat_0, azimuth, range_m)[:2]

            print(end_lat,end_lon)
   
            # Determine text rotation (aligning with the azimuth angle)
            rotation_angle = -azimuth if azimuth <= 180 else 360 - azimuth  # Ensures text is aligned correctly

            # Add labels at the end of each RHI line
            ax[pos].text(end_lon, end_lat, f"{int(azimuth)}°",  # Convert azimuth to an integer for no decimals
                    transform=ccrs.PlateCarree(), 
                    fontsize=10, color='black', fontweight='bold', 
                    ha='center', va='center', rotation=rotation_angle,  # Apply rotation
                    bbox=dict(facecolor='white', edgecolor='grey', alpha=1.0, boxstyle='round,pad=0.3'))

    # Set map extents
    for pos in positions:
        ax[pos].set_extent([lon_min, lon_max, lat_min, lat_max], ccrs.PlateCarree())

    # Add basemap from contextily
    for pos in positions:
        # Use Google Tiles in satellite mode
        google_tiles = cimgt.GoogleTiles(style="satellite")
        # Add Google Tiles as a basemap
        #ax[pos].add_image(google_tiles, 14, alpha=0.6, zorder=6)
        ctx.add_basemap(ax[pos], zoom=8,source=ctx.providers.CartoDB.Positron, crs=projection);

    # Add latitude and longitude gridlines

    for pos in positions:
        gl = ax[pos].gridlines(draw_labels=False, linestyle="None")
        gridliner = ax[pos].gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color='gray')
        gridliner.xformatter = LONGITUDE_FORMATTER
        gridliner.yformatter = LATITUDE_FORMATTER
        gridliner.xlabel_style = {'size': 12, 'color': 'black'}
        gridliner.ylabel_style = {'size': 12, 'color': 'black'}

    # Plot reference points (e.g., radar sites)
    #nxpol1_lat, nxpol1_lon = 51.507125, -2.005448333
    #kepler_lat, kepler_lon = 51.50888, -1.99907

    #display.plot_point(nxpol1_lon, nxpol1_lat, symbol='k+', zorder=12)
    #display.plot_point(kepler_lon, kepler_lat, symbol='k+', zorder=12)

    # Save figure to file if save_path is provided


    print(save_path);
    # Save figure to file if save_path is provided
    dtime_sweep = cftime.num2pydate(radar.time['data'][0],radar.time['units']);
    dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
    figname = f'ncas-radar-camra-1_cao_{dtime_sweep_str}_ppi_el{ppi_el:0.2f}_map_l1_v1.0.1.png';
    plt.savefig(os.path.join(save_path,figname),dpi=300, bbox_inches='tight');
    plt.close();


import pyart
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib.pyplot as plt
import numpy as np
import os
from pyproj import Geod


def make_woest_sop_ppi_map_plot_new(radar, os_key, kml_paths, sweep=0, vmin_vmax=None, save_path=None, separate_plots=False):
    """
    Plots PPI maps of radar variables using OpenStreetMap background.

    Parameters:
    radar: Py-ART Radar object
    os_key (str): API key for Ordnance Survey (optional, if using OS maps).
    kml_paths (dict): Dictionary containing KML file paths for geospatial overlays.
    sweep (int): Radar sweep to plot (default = 0).
    vmin_vmax (dict): Dictionary specifying vmin and vmax for each radar variable.
    save_path (str): Directory path to save the figure.
    separate_plots (bool): If True, creates 4 individual plots instead of a 2x2 matrix.
    """

    ppi_el = radar.get_elevation(0)[0];

    # Radar location and projection setup
    lat_0, lon_0 = radar.latitude["data"][0], radar.longitude["data"][0]
    lambert = ccrs.LambertConformal(central_latitude=lat_0, central_longitude=lon_0)
    
    projection = ccrs.Mercator()

    # Setup PyART gate filtering
    gatefilter = pyart.filters.GateFilter(radar)
    
    # Radar display setup
    display = pyart.graph.RadarMapDisplay(radar)


    # Define radar bounds
    lat_min = np.min(radar.gate_latitude['data'])
    lon_min = np.min(radar.gate_longitude['data'])
    lat_max = np.max(radar.gate_latitude['data'])
    lon_max = np.max(radar.gate_longitude['data'])
    
    # Define color maps for different radar variables
    colormaps = {
        "DBZ_H": "pyart_HomeyerRainbow",
        "VEL_HV": "pyart_balance",
        "ZDR": "pyart_HomeyerRainbow",
        "LDR": "viridis"
    }
    
    # Default vmin and vmax values if not provided
    VEL_HV_min = radar.fields['VEL_HV']['field_limit_lower']
    VEL_HV_max = radar.fields['VEL_HV']['field_limit_upper']
    
    default_vmin_vmax = {
        "DBZ_H": (-20, 50),
        "VEL_HV": (VEL_HV_min, VEL_HV_max),
        "ZDR": (-5, 5),
        "LDR": (-35, 5)
    }
    
    # Radar variables to plot
    variables = ["DBZ_H", "VEL_HV", "ZDR", "LDR"]
    
    if separate_plots:
        fig_size = (12, 10)
        plot_layout = [(1, 1)] * 4  # Each variable gets its own figure
    else:
        fig_size = (30, 30)
        plot_layout = [(2, 2)]  # 2x2 subplot layout
    
    for idx, var in enumerate(variables):
        if separate_plots:
            fig, ax = plt.subplots(figsize=fig_size, subplot_kw={'projection': lambert})
        else:
            if idx == 0:
                fig, axarr = plt.subplots(2, 2, figsize=fig_size, subplot_kw={'projection': lambert})
            ax = axarr[idx // 2, idx % 2]
        
        display.plot_ppi_map(
            var, sweep,
            vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1],
            min_lat=40.0, max_lat=53.0,
            min_lon=-4.0, max_lon=2.0,
            projection=lambert,
            fig=fig,
            ax=ax,
            lat_0=lat_0,
            lon_0=lon_0,
            cmap=colormaps[var],
            colorbar_flag=True,
            alpha=1.0,
            gatefilter=gatefilter,
            edges=True,
            embellish=False,
            resolution="10m",
        )

        ax.set_extent([lon_min, lon_max, lat_min, lat_max], ccrs.PlateCarree())

        ctx.add_basemap(ax, zoom=8,source=ctx.providers.CartoDB.Positron, crs=lambert);

        
        # Plot range rings and cross lines
        display.plot_range_rings([50, 100, 150, 200], '0.5', lw=0.5)
        display.plot_line_xy([-1000, 1000], [0, 0], 'k', lw=0.5)
        display.plot_line_xy([0, 0], [-1000, 1000], 'k', lw=0.5)
        
        print(save_path);
        dtime_sweep = cftime.num2pydate(radar.time['data'][0],radar.time['units']);
        dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
        figname = f'ncas-radar-camra-1_cao_{dtime_sweep_str}_{var}_ppi_el{ppi_el:0.2f}_map_l1_v1.0.1.png';
        plt.savefig(os.path.join(save_path, figname), dpi=300, bbox_inches='tight')
        
        if separate_plots:
            plt.close(fig)
    
    if not separate_plots:
        figname = f'ncas-radar-camra-1_cao_{dtime_sweep_str}_ppi_el{ppi_el:0.2f}_map_l1_v1.0.1.png';
        plt.savefig(os.path.join(save_path, figname), dpi=300, bbox_inches='tight')
        plt.close(fig)



def make_woest_sop_ppi_map_plot_old(radar, os_key, kml_paths, sweep=0, vmin_vmax=None,save_path=None):
    """
    Plots a 2x2 PPI map of radar variables using OpenStreetMap background.

    Parameters:
    radar_file (str): Path to radar file (e.g., .nc or .h5 format).
    os_key (str): API key for Ordnance Survey (optional, if using OS maps).
    kml_paths (dict): Dictionary containing KML file paths for geospatial overlays.
                      Expected keys: 'verticals', 'horizontals', 'outline', 'surface_sites'
    azimuth (int): Azimuth angle for RHI analysis (default is 0).
    """

    # Load Radar Data
    #radar = pyart.io.read(radar_file)

    # Load KML Files into GeoDataFrames
    #gdf_list = {}
    #for key, path in kml_paths.items():
    #    gdf_list[key] = gpd.read_file(path, driver='LIBKML')

    ppi_el = radar.get_elevation(0)[0];

    # Setup OpenStreetMap tiles
    #osm_tiles = cimgt.OSM()
    osm_tiles = GT(desired_tile_form='RGB', style='satellite', url='https://mts0.google.com/vt/lyrs={style}@177000000&hl=en&src=api&x={x}&y={y}&z={z}&s=G')

    projection = ccrs.Mercator()  # Use correct map projection

    # Setup PyART gate filtering
    gatefilter = pyart.filters.GateFilter(radar)
    #gatefilter.exclude_below("SNR", -4)  # Adjust filtering criteria

    # Radar display setup
    display = pyart.graph.RadarMapDisplay(radar)

    # Define color maps for different radar variables
    colormaps = {
        "DBZ_H": "pyart_HomeyerRainbow",
        "VEL_HV": "pyart_balance",
        "ZDR": "pyart_HomeyerRainbow",
        "LDR": "viridis"
    }

    VEL_HV_min = radar.fields['VEL_HV']['field_limit_lower'];
    VEL_HV_max = radar.fields['VEL_HV']['field_limit_upper'];
    
    # Default vmin and vmax values if not provided
    default_vmin_vmax = {
        "DBZ_H": (-20, 50),
        "VEL_HV": (VEL_HV_min, VEL_HV_max),
        "ZDR": (-5, 5),
        "LDR": (-35, 5)
    }

    # Create a 2x2 subplot with correct projection
    fig, ax = plt.subplots(2, 2, figsize=(20, 20), subplot_kw=dict(projection=projection))

    # Define radar bounds
    lat_min = np.min(radar.gate_latitude['data'])
    lon_min = np.min(radar.gate_longitude['data'])
    lat_max = np.max(radar.gate_latitude['data'])
    lon_max = np.max(radar.gate_longitude['data'])

    # Plot PPI maps for different radar fields
    variables = ["DBZ_H", "VEL_HV", "ZDR", "LDR"]
    positions = [(0, 0), (1, 0), (0, 1), (1, 1)]

    # Define RHI azimuths for radial lines
    rhi_azimuths = [0.0, 25.0, 44.0, 60.0, 80.0, 99.0, 120.0, 143.0, 168.0,189.0, 215.0, 232.0, 244.0, 255.0, 270.0, 283.0, 293.0, 310.0, 331.0, 348.0]

    # Maximum radar range in km
    #max_range_km = np.max(radar.range['data']) / 1000.0  
    max_range_km = 200.0;
    
    # Get radar location
    radar_lon, radar_lat = radar.longitude['data'][0], radar.latitude['data'][0]
    
    for var, pos in zip(variables, positions):
        display.plot_ppi_map(
            var, sweep, 
            vmin=default_vmin_vmax[var][0], vmax=default_vmin_vmax[var][1],
            min_lat=50.0, max_lat=53.0,
            min_lon=-4.0, max_lon=2.0,
            lon_lines=np.arange(-4.0, 2.0, 0.5),
            lat_lines=np.arange(40.0, 53.0, 0.5),
            projection=projection,
            fig=fig,
            ax=ax[pos],
            lat_0=radar.latitude["data"][0],
            lon_0=radar.longitude["data"][0],
            cmap=colormaps[var],
            colorbar_flag=False,
            alpha=1.0,
            gatefilter=gatefilter,
            edges=True,
            embellish=True,
            resolution="10m",

        )
        display.plot_colorbar(ax=ax[pos], orient='horizontal', shrink=0.8, cmap=colormaps[var])
        # Plot range rings
        display.plot_range_rings([50, 100, 150, 200], 'k', lw=1)

        # Overlay RHI azimuth lines with rotated labels
        for azimuth in rhi_azimuths:
            azimuth_rad = np.deg2rad(azimuth)
            end_lon = radar_lon + (max_range_km / 111) * np.sin(azimuth_rad)
            end_lat = radar_lat + (max_range_km / 111) * np.cos(azimuth_rad)

            ax[pos].plot([radar_lon, end_lon], [radar_lat, end_lat], transform = projection, #transform=ccrs.PlateCarree(),
                     linestyle='--', color='black', linewidth=0.5, alpha=0.8)

            # Determine text rotation (aligning with the azimuth angle)
            rotation_angle = -azimuth if azimuth <= 180 else 360 - azimuth  # Ensures text is aligned correctly

            # Add labels at the end of each RHI line
            ax[pos].text(end_lon, end_lat, f"{int(azimuth)}°",  # Convert azimuth to an integer for no decimals
                         transform=ccrs.PlateCarree(), 
                         fontsize=10, color='black', fontweight='bold', 
                         ha='center', va='center', rotation=rotation_angle,  # Apply rotation
                         bbox=dict(facecolor='white', edgecolor='black', alpha=0.7, boxstyle='round,pad=0.3'))

    # Set map extents
    for pos in positions:
        ax[pos].set_extent([lon_min, lon_max, lat_min, lat_max], ccrs.PlateCarree())

    # Add basemap from contextily
    for pos in positions:
        ctx.add_basemap(ax[pos], zoom=7, source=ctx.providers.CartoDB.Positron)

    # Add latitude and longitude gridlines
    for pos in positions:
        gridliner = ax[pos].gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color='gray')
        gridliner.xformatter = LONGITUDE_FORMATTER
        gridliner.yformatter = LATITUDE_FORMATTER
        gridliner.xlabel_style = {'size': 12, 'color': 'black'}
        gridliner.ylabel_style = {'size': 12, 'color': 'black'}

    # Plot reference points (e.g., radar sites)
    #nxpol1_lat, nxpol1_lon = 51.507125, -2.005448333
    #kepler_lat, kepler_lon = 51.50888, -1.99907

    #display.plot_point(nxpol1_lon, nxpol1_lat, symbol='k+', zorder=12)
    #display.plot_point(kepler_lon, kepler_lat, symbol='k+', zorder=12)

    print(save_path);
    # Save figure to file if save_path is provided
    dtime_sweep = cftime.num2pydate(radar.time['data'][0],radar.time['units']);
    dtime_sweep_str = dtime_sweep.strftime("%Y%m%d-%H%M%S");
    figname = f'ncas-radar-camra-1_cao_{dtime_sweep_str}_ppi_el{ppi_el:0.2f}_map_l1_v1.0.1.png';
    plt.savefig(os.path.join(save_path,figname),dpi=300);
    plt.close();

    #if save_path:
    #    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    #    print(f"Figure saved to {save_path}")

def make_woest_sop_rhi_plots_day(datestr,inpath,figpath,blflag=False,darkmode=False):
    #inpath_date = os.path.join(inpath,datestr);
    inpath_date = inpath;

    os.chdir(inpath_date);
    print(inpath_date)
    woest_rhi_files = [os.path.join(inpath_date,f) for f in glob.glob('*{}*rhi*.nc'.format(datestr))]

    print(f'srhi files = ',woest_rhi_files);

    for f in woest_rhi_files:
        DS = nc4.Dataset(f, 'r')
        #comment_attribute = DS.getncattr('comment')

        #print(comment_attribute);
        #region_number = re.search(r'\b\d+\b', comment_attribute).group()

        #region_number = re.search(r'\b\w+\b$', comment_attribute).group()
        #print(region_number);
        make_woest_sop_rhi_plot(f,figpath,blflag=blflag,darkmode=darkmode);

    return


os_key = "ngaMRdeOtnGTnXCb378JLN29j3H8AWAo"

from cartopy.io.img_tiles import GoogleTiles as GT


# Define KML file paths
kml_paths = {
    "verticals": "/home/users/cjwalden/WesConGrid_Verticals.kml",
    "horizontals": "/home/users/cjwalden/WesConGrid_Horizontals.kml",
    "outline": "/home/users/cjwalden/WesConGrid_Outline.kml",
    "surface_sites": "/home/users/cjwalden/surface_sites.kml"
}

def make_woest_sop_ppi_map_plots_day(datestr,inpath,figpath):
    #inpath_date = os.path.join(inpath,datestr);
    inpath_date = inpath;

    os.chdir(inpath_date);
    print(inpath_date)
    woest_ppi_files = [f for f in glob.glob('*{}*ppi*.nc'.format(datestr))]

    print(f'ppi files = ',woest_ppi_files);

    save_path = os.path.join(figpath,'ppi',datestr);
    if not os.path.isdir(save_path): 
        os.makedirs(save_path);

    for f in woest_ppi_files:

        radar = pyart.io.read(os.path.join(inpath_date,f))



        for sweep in range(radar.nsweeps):
            radar_sweep = radar.extract_sweeps([sweep])
            radar_adjusted = ensure_increasing_azimuths(radar_sweep)
            #make_woest_sop_ppi_plot(radar_adjusted, save_path);
            make_woest_sop_ppi_map_plot(radar_adjusted, os_key, kml_paths, save_path);
    return


make_woest_sop_rhi_plots_day(datestr,inpath,figpath,blflag=False,darkmode=darkmode);

print(inpath)
print(figpath)
make_woest_sop_ppi_map_plots_day(datestr,inpath,figpath);
