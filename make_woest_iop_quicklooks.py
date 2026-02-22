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
#inpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1',campaign,'L1',datestr,'iop');
inpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1','20230602_woest','v1.0.1','level1','iop',yr,mo,dy);

if darkmode:
    qlstr = 'quicklooks_dark'
else:
    qlstr = 'quicklooks'

#figpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1',campaign,'L1',qlstr,'iop');
figpath = os.path.join(ncas_obs_vol2_path,'cao','processing','ncas-radar-camra-1',f'woest_{qlstr}','v1.0.1','iop');

#figpath = os.path.join(inpath,'quicklooks','iop')

def make_woest_iop_rhi_plot(ncfile,figpath,region,blflag=False,darkmode=False):

    if darkmode:
        plt.style.use('dark_background')

    if blflag:
        hmax = 4;
        xmin = 0;
        xmax = 20;
    else:
        hmax = 12;
        xmin = 0;
        xmax = 125;
    
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

    figpath_region = os.path.join(figpath,f'region_{region}');
    if not os.path.isdir(figpath_region):
        os.makedirs(figpath_region);


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
            plt.savefig(os.path.join(figpath_region,figname),dpi=300);
            plt.close();



def make_woest_iop_rhi_plots_day(datestr,inpath,figpath,blflag=False,darkmode=False):
    #inpath_date = os.path.join(inpath,datestr);
    inpath_date = inpath;

    os.chdir(inpath_date);
    print(inpath_date)
    #woest_srhi_files = [os.path.join(inpath_date,f) for f in glob.glob('*{}*srhi*.nc'.format(datestr))]
    woest_rhi_files = [os.path.join(inpath_date,f) for f in glob.glob('*{}*_rhi*.nc'.format(datestr))]

    #print(f'srhi files = ',woest_srhi_files);
    print(f'rhi files = ',woest_rhi_files);

    #for f in woest_srhi_files:
    for f in woest_rhi_files:
        DS = nc4.Dataset(f, 'r')
        comment_attribute = DS.getncattr('comment')

        print(comment_attribute);
        #region_number = re.search(r'\b\d+\b', comment_attribute).group()

        region_number = re.search(r'\b\w+\b$', comment_attribute).group()
        print(region_number);
        make_woest_iop_rhi_plot(f,figpath,f'{region_number}',blflag=blflag,darkmode=darkmode);

    return


make_woest_iop_rhi_plots_day(datestr,inpath,figpath,blflag=False,darkmode=darkmode);

