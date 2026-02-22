#!/usr/bin/env python

import getopt, sys, os

import datetime

import netCDF4 as nc4

import pyart
import numpy as np
import numpy.ma as ma
import shutil
import glob
import gzip

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import cmocean
import getpass, socket

import pandas as pd

import cftime

version = 0.1

sys.path.append('/home/users/cjwalden/git/camra-radar-utils')

import campaign_processing

from pathlib import Path
homepath = Path.home()


try:
    opts, args = getopt.getopt(sys.argv[1:], "d:i:o:", ["date=","inpath=","outpath="])
except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized
        sys.exit(2)

data_date = datetime.datetime.now()
datestr = data_date.strftime('%Y%m%d')
tracking_tag = 'AMOF_20230201132601';

campaign = 'ccrest-m';

data_version = "1.0.0"

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
yaml_project_file = os.path.join(script_dir, 'campaigns', f'{campaign}_project.yml')
yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

print(yaml_project_file);

for o, a in opts:
    if o == "-d":
        datestr = a;
    elif o == "-i":
        inpath = a;
    elif o == "-o":
        outpath = a;
    else:
        assert False, "unhandled option"


camrapath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-camra-1';
ncas_radar_path = '/gws/nopw/j04/ncas_radar_vol1';
amof_proc_path = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1'
inpath = os.path.join(camrapath,'data','campaign',campaign,'ts',datestr);

outpath = os.path.join(amof_proc_path,campaign); #,'L0b-ts',datestr);


l0bpath = os.path.join(outpath,'ts','L0b',datestr)
l1path = os.path.join(outpath,'ts','L1',datestr)
mompath = os.path.join(amof_proc_path,campaign,'L1a');

mom_l1b_path = mompath.replace('L1a','L1b')

print(mom_l1b_path);
if not os.path.isdir(mom_l1b_path):
    os.makedirs(mom_l1b_path);

print(f'l0bpath = {l0bpath}')
if not os.path.isdir(l0bpath): 
    os.makedirs(l0bpath);

if not os.path.isdir(l1path): 
    os.makedirs(l1path);


campaign_processing.process_camra_ccrest_vpt_day_ts(datestr,inpath,mompath,outpath,yaml_project_file,yaml_instrument_file,data_version=data_version);

