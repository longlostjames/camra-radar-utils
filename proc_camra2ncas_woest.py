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

import camra_utils

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
tracking_tag = ' AMOF_20220922221548';

campaign = 'woest';

data_version = "1.0.0"


yaml_project_file = os.path.join(homepath,'amof_campaigns',f'{campaign}_project.yml')
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
amof_proc_path = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1'
inpath = os.path.join(camrapath,'data','campaign',campaign,'raw',datestr);

outpath = os.path.join(amof_proc_path,campaign,'L1a',datestr);

if not os.path.isdir(outpath): 
    os.makedirs(outpath);

camra_utils.process_camra_woest_day_step1(datestr,inpath,outpath,yaml_project_file,yaml_instrument_file,data_version=data_version);
