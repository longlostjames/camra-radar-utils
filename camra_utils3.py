#!/usr/bin/python

# ==============================================================================
# CAMRa POSTPROCESSING UTILITIES MODULE
# Tools related to postprocessing NetCDF data from CAMRa
#
# Author:        Chris Walden, NCAS
# History:
# Version:	 0.2
# Last modified: 30/01/20
# ==============================================================================
module_version = 0.2

# ------------------------------------------------------------------------------
# Import required tools
# ------------------------------------------------------------------------------
import numpy as np
import numpy.ma as ma;
import os, re, sys, getopt, shutil, zipfile, string, pwd, getpass
import netCDF4 as nc4
import socket

import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt  #import plotting package
import pyart


from datetime import tzinfo, datetime, time
#from pylab import *


# ----------------------
# Define useful function
# ----------------------
def in_interval(seq,xmin,xmax):
    for i, x in enumerate(seq):
        if x>=xmin and x<xmax:
            yield i

# ------------------------------------------------------------------------------
# Define the netcdf recalibration function
# ------------------------------------------------------------------------------
def recalibrate(ncfile,zedcal,zdrcal):

    user = getpass.getuser()
    print(user)

    zedcal = float(zedcal)
    zdrcal = float(zdrcal)
    #pdpcal = float(pdpcal)
    #ldrcal = float(ldrcal)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    # -------------------------------
    # Get existing calibration values
    # -------------------------------

    zed = dataset.variables['ZED_H']
    zdr = dataset.variables['ZDR']
    pdp = dataset.variables['PDP']
    ldr = dataset.variables['LDR']

    # -------------------------------
    # Get existing calibration values
    # -------------------------------
    zed_oldcal = zed.getncattr('applied_calibration_offset')
    zdr_oldcal = zdr.getncattr('applied_calibration_offset')
    #pdp_oldcal = pdp.getncattr('applied_calibration_offset')
    #ldr_oldcal = ldr.getncattr('applied_calibration_offset')

    zed_recal = np.round(zedcal - zed_oldcal,2);
    zdr_recal = np.round(zdrcal - zdr_oldcal,2);
    #pdp_recal = pdpcal - pdp_oldcal
    #ldr_recal = ldrcal - ldr_oldcal

    print('zed_oldcal = ', zed_oldcal)
    print('zdr_oldcal = ', zdr_oldcal)
    #print('pdp_oldcal = ', pdp_oldcal)
    #print('ldr_oldcal = ', ldr_oldcal)

    print('zed_recal = ', zed_recal)
    print('zdr_recal = ', zdr_recal)
    #print('pdp_recal = ', pdp_recal)
    #print('ldr_recal = ', ldr_recal)

    caltime = datetime.utcnow()
    caltimestr = caltime.ctime()

    zed[:] += zed_recal;
    zdr[:] += zdr_recal;

    #np.where(zed>-999.0,zed+zed_recal,zed)
    #np.where(zdr>-999.0,zdr+zdr_recal,zdr)
    #np.where(pdp>-999.0,pdp+pdp_recal,pdp)
    #np.where(ldr>-999.0,ldr+ldr_recal,ldr)

    zed.applied_calibration_offset = zedcal
    zdr.applied_calibration_offset = zdrcal
    #pdp.applied_calibration_offset = pdpcal
    #ldr.applied_calibration_offset = ldrcal

    history = caltimestr + (" - user:" + user
        + " machine: " + socket.gethostname()
        + " program: camra_utils.py recalibrate_netcdf_camra"
        + " version:" + str(module_version) + " arguments: "
        + ncfile + " " + str(zedcal) + " " + str(zdrcal) + "\n")
#        + str(pdpcal) + " " + str(ldrcal) + "\n")

    dataset.history += history
    print(dataset.history)
    dataset.close()

# ------------------------------------------------------------------------------
# Define function to produce clutter classification based on LDR threshold
# ------------------------------------------------------------------------------
def classify_clutter(ncfile,ldrmax):

    user = getpass.getuser()
    print(user)

    # -----------------
    # Get LDR threshold
    # -----------------
    ldrmax = float(ldrmax)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    nc = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    # ------------------------------
    # Create classification variable
    # ------------------------------
    classification = nc.createVariable('classification','i1',
        ('time','range'),'fill_value',-99)

    classification.long_name = 'Classification'
    classification.short_name = 'Classification'
    classification.units = ''
    classification.units_html = ''
    classification.comment = ('Clutter classification using LDR threshold.\n'
        + 'Data with LDR greater than LDRmax are classed as'
        + ' clutter and set to clutter_val')
    classification.LDRmax = ldrmax
    missing_val = -99
    clutter_val = 1

    classification.missing_value = missing_val
    classification.clutter       = clutter_val

    # ----------------------------
    # Initialise to missing values
    # ----------------------------
    classification[:]            = missing_val

    # -----------
    # Read in LDR
    # -----------
    ldr = nc.variables['LDR']

    modtime = datetime.utcnow()
    modtimestr = modtime.ctime()

    history = modtimestr + (" - user:" + user
        + " machine: " + socket.gethostname()
        + " program: camra_utils.py classify_clutter" + " version:"
        + str(module_version) + " arguments: " + ncfile
        + " " + str(ldrmax) + "\n")

    iclutter = nonzero(ldr[:,:]>ldrmax)

    print(iclutter)

    classification[iclutter] = clutter_val

    nc.history += "\n" + history
    print(nc.history)
    nc.close

# ------------------------------------------------------------------------------
# Set as "missing value" gates 591 onwards for phase-derived parameters
# Needed for historic data where the transmit I & Q were not properly stored
# ------------------------------------------------------------------------------
def clean_distant_gates(ncfile):

    user = getpass.getuser()
    print(user)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    nc = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    rng = nc.variables['range']
    ngates = size(rng[:])
    gate_width = rng.getncattr('gate_width')
    gates_per_cell = int(rng.getncattr('gates_range_averaged'))

    r0 = rng[0]
    print("r0 = ", r0)
    ifargates = nonzero(rng-r0>172600)[0]
    g1 = ifargates[0]

    print("gate_width =" , gate_width, " rounded ", int(gate_width+0.5))
    print("gates_per_cell = ", gates_per_cell)
    #if (int(gate_width+0.5)==75)  & (gates_per_cell==1):
    #   g1 = 2355
    #else:
    #   g1 = 591


    # -------------------------------
    # Get phase-derived parameters
    # -------------------------------
    varlist = nc.variables.keys()

    print(len(ifargates))
    print(ifargates)
    print(ngates)

    if len(ifargates)>0:

       if 'PDP' in varlist:
          print('PDP present')
          pdp                  = nc.variables['PDP']
          pdp_missing          = pdp.getncattr('missing_value')
          print("missing_value = ", pdp_missing)
          #print g1,", ",ngates
          print(pdp[:,g1:ngates])
          pdp[:,ifargates]     = pdp_missing
          print("shape ", pdp.shape)
          print(pdp[:,g1:ngates])

       if 'VEL_H' in varlist:
          print('VEL_H present')
          vel_h                = nc.variables['VEL_H']
          vel_h_missing        = vel_h.getncattr('missing_value')
          vel_h[:,g1:ngates]   = vel_h_missing

       if 'VEL_HV' in varlist:
          print('VEL_HV present')
          vel_hv               = nc.variables['VEL_HV']
          vel_hv_missing       = vel_hv.getncattr('missing_value')
          vel_hv[:,g1:ngates]  = vel_hv_missing

       if 'VEL_UHV' in varlist:
          print('VEL_UHV present')
          vel_uhv              = nc.variables['VEL_UHV']
          vel_uhv_missing      = vel_uhv.getncattr('missing_value')
          vel_uhv[:,g1:ngates] = vel_uhv_missing

       if 'DDV' in varlist:
          print('DDV present')
          ddv                  = nc.variables['DDV']
          ddv_missing          = ddv.getncattr('missing_value')
          ddv[:,g1:ngates]     = ddv_missing

       if 'SPW_H' in varlist:
          print('SPW_H present')
          spw_h                = nc.variables['SPW_H']
          spw_h_missing        = spw_h.getncattr('missing_value')
          spw_h[:,g1:ngates]   = spw_h_missing

       if 'SPW_HV' in varlist:
          print('SPW_HV present')
          spw_hv               = nc.variables['SPW_HV']
          spw_hv_missing       = spw_hv.getncattr('missing_value')
          spw_hv[:,g1:ngates]  = spw_hv_missing

    reproctime = datetime.utcnow()
    reproctimestr = reproctime.ctime()

    history = reproctimestr + (" - user:" + user
        + " machine: " + socket.gethostname()
        + " program: camra_utils.py clean_distant_gates"
        + " version:" + str(module_version))

    print(size(pdp[:,0]))
    #print history
    nc.history += "\n" + history
    print(nc.history)
    nc.close

# ------------------------------------------------------------------------------
# Get history from NetCDF file
# ------------------------------------------------------------------------------
def get_history(ncfile):

    user = getpass.getuser()
    print(user)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    nc = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    print(nc.history)
    nc.close

# ------------------------------------------------------------------------------
# Update global attributes - institution and references
# ------------------------------------------------------------------------------
def update_metadata(ncfile):

    user = getpass.getuser()
    print(user)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    nc = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    institution  = ("National Centre for Atmospheric Science,"
                    + " UK: https://www.ncas.ac.uk\n")
    institution += ("Data held at the Centre for Environmental Data Analysis,"
                    + " UK: https://www.ceda.ac.uk")
    nc.institution = institution

    references = "doi: 10.1049/ecej:19940205"
    nc.references = references

    updttime = datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: camra_utils.py update_metadata"
    + " version:" + str(module_version))

    nc.history += "\n" + history
    print(nc.history)
    nc.close


# ------------------------------------------------------------------------------
# Perform noise subtraction
# ------------------------------------------------------------------------------
def subtract_noise_polar(ncfile,nezh,nezv,nezhv,threshold):

    user = getpass.getuser()
    print(user)

    thresh       = np.float(threshold);
    nezh         = np.float(nezh);
    nezv         = np.float(nezv);
    nezx         = np.float(nezhv);

    print(nezh,nezv,nezx,thresh);

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    nray    = len(dataset.dimensions['time']);
    ngate   = len(dataset.dimensions['range']);

    elv = np.transpose(np.tile(dataset.variables['elevation'][:],(ngate,1)));
    rng = np.tile(dataset.variables['range'][:],(nray,1));

    height = rng*np.sin(elv*np.pi/180.)

    zh = np.copy(dataset.variables['ZED_H']);
    zed = ma.masked_where(height<14000, zh[:]);

    rngkm = ma.masked_where(rng<=0.0, rng/1000.);

    range2 = 20.*ma.log10(rngkm);

    zh[:] = zed - range2;
    zv = zh.copy();
    zv[:] = zh[:] - dataset.variables['ZDR'][:]

    zx = zh.copy();
    zx[:] = zh[:] + dataset.variables['LDR'][:]

    rangesq = rng*rng/(1000.*1000.0);

    Pnezh = np.power(10,0.1*nezh);
    Pnezv = np.power(10,0.1*nezv);
    Pnezx = np.power(10,0.1*nezx);

    PZh = dataset.variables['ZED_H'][:];
    PZh = np.power(10,0.1*PZh)/rangesq;

    dataset.variables['ZED_H'][:] = ma.masked_where(PZh<Pnezh*thresh,
                        10.0*ma.log10(PZh-Pnezh)+range2);

    PZv = dataset.variables['ZED_H'][:]-dataset.variables['ZDR'][:];
    PZv = np.power(10,0.1*PZv)/rangesq;
    dataset.variables['ZDR'][:] = ma.masked_where(PZv<Pnezv*thresh,
                        10.0*ma.log10((PZh-Pnezh)/(PZv-Pnezv)));

    PZx = dataset.variables['ZED_H'][:]+dataset.variables['LDR'][:];
    PZx = np.power(10,0.1*PZx)/rangesq;
    dataset.variables['LDR'][:] = ma.masked_where(PZx<Pnezx*thresh,
                        10.0*ma.log10((PZx-Pnezx)/(PZh-Pnezh)));

    reproctime = datetime.utcnow()
    reproctimestr = reproctime.ctime()

    history = reproctimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: camra_utils.py subtract_noise_polar"
    + " {} {} {} {}".format(nezh, nezv, nezhv,  thresh)
    + " version:" + str(module_version))

    dataset.history += "\n" + history
    print(dataset.history)


    dataset.close()

    return

# ------------------------------------------------------------------------------
# Apply PHI_DP offset
# ------------------------------------------------------------------------------
def apply_phidp_offset(ncfile,phidp_offset):

    user = getpass.getuser()
    print(user)

    phidp_offset = np.float(phidp_offset);

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    dataset.variables['PDP'][:] -= phidp_offset;

    reproctime = datetime.utcnow()
    reproctimestr = reproctime.ctime()

    history = reproctimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: camra_utils.py apply_phidp_offset"
    + " {}".format(phidp_offset)
    + " version:" + str(module_version))

    dataset.history += "\n" + history
    print(dataset.history)


    dataset.close()

    return


# ------------------------------------------------------------------------------
# Perform clean to remove isolated pixels
# ------------------------------------------------------------------------------
def clean(ncfile):

    user = getpass.getuser()
    print(user)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    nray    = len(dataset.dimensions['time']);
    ngate   = len(dataset.dimensions['range']);

    zh  = dataset.variables['ZED_H'][:].copy();
    zdr = dataset.variables['ZDR'][:].copy();

    for iray in np.arange(nray):
        #print("iray = ",iray, "\n");
        for igate in np.arange(ngate-1):
            if (ma.is_masked(zh[iray,igate-1]) and
                ma.is_masked(zh[iray,igate+1])):
                zh[iray,igate] = ma.masked;

    dataset.variables['ZED_H'][:] = zh;

    zh_mask = ma.getmask(zh)

    for var in dataset.variables:
        if dataset.variables[var].dimensions == (u'time', u'range'):
            val = dataset.variables[var][:]
            val = ma.masked_where(zh_mask, val);
            dataset.variables[var][:]= val;

    reproctime = datetime.utcnow()
    reproctimestr = reproctime.ctime()

    history = reproctimestr + (" - user:" + user
        + " machine: " + socket.gethostname()
        + " program: camra_utils.py clean"
        + " version:" + str(module_version))

    dataset.history += "\n" + history
    print(dataset.history)

    dataset.close()

def doubleclean(ncfile):

    user = getpass.getuser()
    print(user)

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    nray    = len(dataset.dimensions['time']);
    ngate   = len(dataset.dimensions['range']);

    zh  = dataset.variables['ZED_H'][:].copy();

    for iray in np.arange(nray):
        for igate in np.arange(ngate-2):
            if (ma.is_masked(zh[iray,igate-1])
                and ma.is_masked(zh[iray,igate+1])
                or ma.is_masked(zh[iray,igate+2])):
                    zh[iray,igate] = ma.masked

    dataset.variables['ZED_H'][:] = zh;

    zh_mask = ma.getmask(zh)

    for var in dataset.variables:
        if dataset.variables[var].dimensions == (u'time', u'range'):
            val = dataset.variables[var][:]
            val = ma.masked_where(zh_mask, val);
            dataset.variables[var][:]= val;

    reproctime = datetime.utcnow()
    reproctimestr = reproctime.ctime()

    history = reproctimestr + (" - user:" + user
        + " machine: " + socket.gethostname()
        + " program: camra_utils.py doubleclean"
        + " version:" + str(module_version))

    dataset.history += "\n" + history
    print(dataset.history)

    dataset.close()


# ------------------------------------------------------------------------------
# Derive noise from gates above 14km
# ------------------------------------------------------------------------------
MAX_ERR = 0.5;

def get_noise_levels(ncfile):

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    nray    = len(dataset.dimensions['time']);
    ngate   = len(dataset.dimensions['range']);

    elv = np.transpose(np.tile(dataset.variables['elevation'][:],(ngate,1)));
    rng = np.tile(dataset.variables['range'][:],(nray,1))

    height = rng*np.sin(elv*np.pi/180.)

    zh = dataset.variables['ZED_H'][:];
    zed = ma.masked_where(height<14000, zh);

    rngkm = ma.masked_where(rng<=0.0, rng/1000.);

    range2 = 20.*ma.log10(rngkm);

    zh[:] = zed - range2;
    zv = zh.copy();
    zv[:] = zh[:] - dataset.variables['ZDR'][:]

    zx = zh.copy();
    zx[:] = zh[:] + dataset.variables['LDR'][:]

    nezharr = ma.mean(zh,axis=1)
    nezherr = ma.std(zh,axis=1)
    nezvarr = ma.mean(zv,axis=1)
    nezverr = ma.std(zv,axis=1)
    nezxarr = ma.mean(zx,axis=1)
    nezxerr = ma.std(zx,axis=1)

    nezharr = ma.masked_where(nezherr>MAX_ERR,nezharr)
    nezvarr = ma.masked_where(nezverr>MAX_ERR,nezvarr)
    nezxarr = ma.masked_where(nezxerr>MAX_ERR,nezxarr)

    nezh = ma.median(nezharr)
    nezv = ma.median(nezvarr)
    nezx = ma.median(nezxarr)

    dataset.close()


    return np.round(nezh,2), np.round(nezv,2), np.round(nezx,2)

# ------------------------------------------------------------------------------
# Derive PHI_DP offset
# ------------------------------------------------------------------------------
def get_phidp_offset(ncfile):

    # ----------------
    # Open NetCDF file
    # ----------------
    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    nray    = len(dataset.dimensions['time']);
    ngate   = len(dataset.dimensions['range']);

    phidp = dataset.variables['PDP'][:];

    phi0arr = ma.mean(phidp,axis=1)

    for iray in np.arange(nray):
        pdpsum = 0
        nphipixel = 0
        phidp_use = phidp[iray,:].copy()
        for igate in np.arange(6,ngate):
            tmp_mean = ma.mean(phidp_use[igate-4:igate])
            tmp_var  = ma.var(phidp_use[igate-4:igate])
            if tmp_var<9:
                pdpsum = pdpsum+tmp_mean;
                nphipixel = nphipixel+1
                break

        if nphipixel>0:
            phidp_floor = pdpsum/nphipixel
            phi0arr[iray] = phidp_floor

    phi0 = ma.median(phi0arr)

    dataset.close()

    return np.round(phi0,2)


# ------------------------------------------------------------------------
# Define function to read noise files and extract median values
# ------------------------------------------------------------------------
def get_median_calnoise(noisefile):

    import codecs
    # -------------------------------
    # Read in required data from file
    # -------------------------------
    m = 0
    f = codecs.open(noisefile,'r', encoding='utf-8',errors='ignore')
    print(noisefile)

    for k in range(2):
        line = f.readline()

    lines = f.readlines()

    f.close()

    pdpcal = []
    nh = []
    nv = []
    nhv = []
    for elem in lines:
        cols = elem.split()
        if len(cols)==5:
           pdpcal.append(float(cols[4]))
        if len(cols)==7:
           nhstr  = cols[4]
           nvstr  = cols[5]
           nhvstr = cols[6]
           l = nhstr.split(chr(177))
           nh.append(float(l[0]))
           l = nvstr.split(chr(177))
           nv.append(float(l[0]))
           l = nhvstr.split(chr(177))
           nhv.append(float(l[0]))
           #nh.append(float(nhstr[0:-6]))
           #nv.append(float(nvstr[0:-6]))
           #nhv.append(float(nhvstr[0:-6]))

    noiseh   = np.array(nh)
    noisev   = np.array(nv)
    noisehv  = np.array(nhv)
    phidpcal = np.array(pdpcal)

    return np.median(noiseh),np.median(noisev),np.median(noisehv),np.median(phidpcal)


# ------------------------------------------------------------------------
# Define function to produce processed files from raw files
# ------------------------------------------------------------------------
def create_processed_file(infile):

    user = getpass.getuser()
    print(user)

    proctime = datetime.utcnow()
    proctimestr = proctime.ctime()

    history = proctimestr + (" - user:" + user
        + " machine: " + socket.gethostname()
        + " program: camra_utils.py create_processed_file"
        + " version:" + str(module_version))

    # ---------------------------
    # Open netCDF file for input
    # ---------------------------
    nc0 = nc4.Dataset(infile, 'r',format='NETCDF3_CLASSIC')

    if len(nc0.dimensions['time'])==0:
        nc0.close()
        return

    nray   = len(nc0.dimensions['time']);
    ngate0 = len(nc0.dimensions['range']);
    ngate1 = ngate0;

    time0       = nc0.variables['time'];
    range0      = nc0.variables['range'];

    elevation0  = nc0.variables['elevation'];
    azimuth0    = nc0.variables['azimuth'];

    latitude0   = nc0.variables['latitude'];
    longitude0  = nc0.variables['longitude'];
    height0     = nc0.variables['height'];

    frequency0  = nc0.variables['frequency'];
    prf0        = nc0.variables['prf'];
    beamwidthH0 = nc0.variables['beamwidthH'];
    beamwidthV0 = nc0.variables['beamwidthV'];

    antenna_diameter0 = nc0.variables['antenna_diameter'];
    pulse_period0     = nc0.variables['pulse_period'];
    transmit_power0   = nc0.variables['transmit_power'];
    clock0            = nc0.variables['clock'];

    # ---------------------------
    # Open netCDF file for output
    # ---------------------------
    outfile=infile.replace("raw.nc",".nc");

    nc1 = nc4.Dataset(outfile,'w',format="NETCDF4")

    print('Creating dimensions in output file')

    #Copy dimensions
    nc1.createDimension("time", None)
    nc1.createDimension("range",ngate1)

    try:
        ZED_H0 = nc0.variables['ZED_H'];
        print('Creating ZED_H variable in output file');
        ZED_H1 = nc1.createVariable("ZED_H",ZED_H0.datatype,('time','range'));
        print('Copying ZED_H attributes');
        ZED_H1.setncatts({k: ZED_H0.getncattr(k) for k in ZED_H0.ncattrs()});
        ZED_H1.setncattr_string('proposed_standard_name',"radar_equivalent_reflectivity_factor_h");
        ZED_H1[:,:] = ZED_H0[:,:];
    except KeyError:
        print("ZED_H doesn't exist");

    try:
        ZDR0 = nc0.variables['ZDR'];
        print('Creating ZDR variable in output file');
        ZDR1 = nc1.createVariable("ZDR",ZDR0.datatype,('time','range'));
        print('Copying ZDR attributes');
        ZDR1.setncatts({k: ZDR0.getncattr(k) for k in ZDR0.ncattrs()});
        ZDR1.setncattr_string('proposed_standard_name',"radar_differential_reflectivity_hv");
        ZDR1[:,:] = ZDR0[:,:];
    except KeyError:
        print("ZDR doesn't exist");

    try:
        LDR0 = nc0.variables['LDR'];
        print('Creating LDR variable in output file');
        LDR1 = nc1.createVariable("LDR",LDR0.datatype,('time','range'));
        print('Copying LDR attributes');
        LDR1.setncatts({k: LDR0.getncattr(k) for k in LDR0.ncattrs()});
        LDR1.setncattr_string('proposed_standard_name',"radar_linear_depolarization_ratio_h");
        LDR1[:,:] = LDR0[:,:];
    except KeyError:
        print("LDR doesn't exist");

    try:
        CXC0 = nc0.variables['CXC'];
        print('Creating CXC variable in output file');
        CXC1 = nc1.createVariable("CXC",CXC0.datatype,('time','range'));
        print('Copying CXC attributes');
        CXC1.setncatts({k: CXC0.getncattr(k) for k in CXC0.ncattrs()});
        CXC1[:,:] = CXC0[:,:];
    except KeyError:
        print("CXC doesn't exist");

    try:
        PDP0 = nc0.variables['PDP'];
        print('Creating PDP variable in output file');
        PDP1 = nc1.createVariable("PDP",PDP0.datatype,('time','range'));
        print('Copying PDP attributes');
        PDP1.setncatts({k: PDP0.getncattr(k) for k in PDP0.ncattrs()});
        PDP1.setncattr_string('proposed_standard_name',"radar_differential_phase_hv");
        PDP1[:,:] = PDP0[:,:];
    except KeyError:
        print("PDP doesn't exist");

    try:
        VEL_HV0 = nc0.variables['VEL_HV'];
        print('Creating VEL_HV variable in output file');
        VEL_HV1 = nc1.createVariable("VEL_HV",VEL_HV0.datatype,('time','range'));
        print('Copying VEL_HV attributes');
        VEL_HV1.setncatts({k: VEL_HV0.getncattr(k) for k in VEL_HV0.ncattrs()});
        VEL_HV1.setncattr_string('standard_name',"radial_velocity_of_scatterers_away_from_instrument");
        VEL_HV1[:,:] = VEL_HV0[:,:];
    except KeyError:
        print("VEL_HV doesn't exist");

    try:
        VEL_H0 = nc0.variables['VEL_H'];
        print('Creating VEL_H variable in output file');
        VEL_H1 = nc1.createVariable("VEL_H",VEL_H0.datatype,('time','range'));
        print('Copying VEL_H attributes');
        VEL_H1.setncatts({k: VEL_H0.getncattr(k) for k in VEL_H0.ncattrs()});
        VEL_H1.setncattr_string('proposed_standard_name',"radial_velocity_of_scatterers_away_from_instrument_h");
        VEL_H1[:,:] = VEL_H0[:,:];
    except KeyError:
        print("VEL_HV doesn't exist");

    try:
        VEL_V0 = nc0.variables['VEL_V'];
        print('Creating VEL_V variable in output file');
        VEL_V1 = nc1.createVariable("VEL_V",VEL_V0.datatype,('time','range'));
        print('Copying VEL_V attributes');
        VEL_V1.setncatts({k: VEL_V0.getncattr(k) for k in VEL_V0.ncattrs()});
        VEL_V1.setncattr_string('proposed_standard_name',"radial_velocity_of_scatterers_away_from_instrument_v");
        VEL_V1[:,:] = VEL_V0[:,:];
    except KeyError:
        print("VEL_V doesn't exist");

    try:
        DDV0 = nc0.variables['DDV'];
        print('Creating DDV variable in output file');
        DDV1 = nc1.createVariable("DDV",DDV0.datatype,('time','range'));
        print('Copying DDV attributes');
        DDV1.setncatts({k: DDV0.getncattr(k) for k in DDV0.ncattrs()});
        DDV1[:,:] = DDV0[:,:];
    except KeyError:
        print("DDV doesn't exist");

    try:
        SPW_HV0     = nc0.variables['SPW_HV'];
        print('Creating SPW_HV variable in output file');
        SPW_HV1 = nc1.createVariable("SPW_HV",SPW_HV0.datatype,('time','range'));
        print('Copying SPW_HV attributes');
        SPW_HV1.setncatts({k: SPW_HV0.getncattr(k) for k in SPW_HV0.ncattrs()});
        SPW_HV1.setncattr_string('proposed_standard_name',"radar_doppler_spectrum_width");
        SPW_HV1[:,:] = SPW_HV0[:,:];
    except Keyerror:
        print("SPW_HV doesn't exist");

    try:
        SPW_H0      = nc0.variables['SPW_H'];
        print('Creating SPW_H variable in output file');
        SPW_H1 = nc1.createVariable("SPW_H",SPW_H0.datatype,('time','range'));
        print('Copying SPW_H attributes');
        SPW_H1.setncatts({k: SPW_H0.getncattr(k) for k in SPW_H0.ncattrs()});
        SPW_H1.setncattr_string('proposed_standard_name',"radar_doppler_spectrum_width_h");
        SPW_H1[:,:] = SPW_H0[:,:];
    except KeyError:
        print("SPW_H doesn't exist")

    try:
        PHI_HV0     = nc0.variables['PHI_HV'];
        print('Creating PHI_HV variable in output file');
        PHI_HV1 = nc1.createVariable("PHI_HV",PHI_HV0.datatype,('time','range'));
        print('Copying PHI_HV attributes');
        PHI_HV1.setncatts({k: PHI_HV0.getncattr(k) for k in PHI_HV0.ncattrs()});
        PHI_HV1[:,:] = PHI_HV0[:,:];
    except KeyError:
        print("PHI_HV doesn't exist");

    try:
        PHI_H0      = nc0.variables['PHI_H'];
        print('Creating PHI_H variable in output file');
        PHI_H1 = nc1.createVariable("PHI_H",PHI_H0.datatype,('time','range'));
        print('Copying PHI_H attributes');
        PHI_H1.setncatts({k: PHI_H0.getncattr(k) for k in PHI_H0.ncattrs()});
        PHI_H1[:,:] = PHI_H0[:,:];
    except KeyError:
        print("PHI_H doesn't exist");

    try:
        PHI_HVD0    = nc0.variables['PHI_HVD'];
        print('Creating PHI_HVD variable in output file');
        PHI_HVD1 = nc1.createVariable("PHI_HVD",PHI_HVD0.datatype,('time','range'));
        print('Copying PHI_HVD attributes');
        PHI_HVD1.setncatts({k: PHI_HVD0.getncattr(k) for k in PHI_HVD0.ncattrs()});
        PHI_HVD1[:,:] = PHI_HVD0[:,:];
    except KeyError:
        print("PHI_HVD doesn't exist");


    print('Creating scalar variables in output file');
    latitude1         = nc1.createVariable("latitude"  ,latitude0.datatype,[])
    longitude1        = nc1.createVariable("longitude" ,longitude0.datatype,[])
    height1           = nc1.createVariable("height"    ,height0.datatype,[])
    frequency1        = nc1.createVariable("frequency" ,frequency0.datatype,[])
    prf1              = nc1.createVariable("prf"       ,prf0.datatype,[])
    beamwidthH1       = nc1.createVariable("beamwidthH",beamwidthH0.datatype,[])
    beamwidthV1       = nc1.createVariable("beamwidthV",beamwidthV0.datatype,[])
    antenna_diameter1 = nc1.createVariable("antenna_diameter",antenna_diameter0.datatype,[])
    pulse_period1     = nc1.createVariable("pulse_period"  ,pulse_period0.datatype,[])
    transmit_power1   = nc1.createVariable("transmit_power",transmit_power0.datatype,[])
    clock1            = nc1.createVariable("clock"    ,clock0.datatype,[])

    print('Copying attributes of scalar variables');

    latitude1.setncatts({k: latitude0.getncattr(k) for k in latitude0.ncattrs()})
    longitude1.setncatts({k: longitude0.getncattr(k) for k in longitude0.ncattrs()})

    height1.setncatts({k: height0.getncattr(k) for k in height0.ncattrs()})

    frequency1.setncatts({k: frequency0.getncattr(k) for k in frequency0.ncattrs()})
    prf1.setncatts({k: prf0.getncattr(k) for k in prf0.ncattrs()})

    beamwidthH1.setncatts({k: beamwidthH0.getncattr(k) for k in beamwidthH0.ncattrs()})
    beamwidthV1.setncatts({k: beamwidthV0.getncattr(k) for k in beamwidthV0.ncattrs()})
    antenna_diameter1.setncatts({k: antenna_diameter0.getncattr(k) for k in antenna_diameter0.ncattrs()})

    pulse_period1.setncatts({k: pulse_period0.getncattr(k) for k in pulse_period0.ncattrs()})
    transmit_power1.setncatts({k: transmit_power0.getncattr(k) for k in transmit_power0.ncattrs()})
    clock1.setncatts({k: clock0.getncattr(k) for k in clock0.ncattrs()})

    print('Copying scalar variables');
    latitude1[:]         = latitude0[:];
    longitude1[:]        = longitude0[:];
    height1[:]           = height0[:];
    frequency1[:]        = frequency0[:];
    prf1[:]              = prf0[:];
    beamwidthH1[:]       = beamwidthH0[:];
    beamwidthV1[:]       = beamwidthV0[:];
    antenna_diameter1[:] = antenna_diameter0[:];
    pulse_period1[:]     = pulse_period0[:];
    transmit_power1[:]   = transmit_power0[:];
    clock1[:]            = clock0[:];

    print('Copying global attributes');
    nc1.setncatts({k: nc0.getncattr(k) for k in nc0.ncattrs()});

    print('Creating time variable in output file');
    time1  = nc1.createVariable("time",time0.datatype,time0.dimensions);
    print('Copying time attributes');
    time1.setncatts({k: time0.getncattr(k) for k in time0.ncattrs()});
    print('Copying time variable');
    time1[:] = time0[:];

    print('Creating range variable in output file');
    range1 = nc1.createVariable("range",range0.datatype,('range'));
    print('Copying range attributes');
    range1.setncatts({k: range0.getncattr(k) for k in range0.ncattrs()});

    print('Creating azimuth variable in output file');
    azimuth1 = nc1.createVariable("azimuth",azimuth0.datatype,('time'));
    print('Copying azimuth attributes');
    azimuth1.setncatts({k: azimuth0.getncattr(k) for k in azimuth0.ncattrs()});

    print('Creating elevation variable in output file');
    elevation1 = nc1.createVariable("elevation",elevation0.datatype,('time'));
    print('Copying azimuth attributes');
    elevation1.setncatts({k: elevation0.getncattr(k) for k in elevation0.ncattrs()});

    range1[:]     = range0[:];
    elevation1[:] = elevation0[:];
    azimuth1[:]  = azimuth0[:];

    nc0.close();

    oldhistory = nc1.history;

    if oldhistory.endswith('\n'):
        nc1.history += history;
    else:
        nc1.history += '\n' + history;


    nc1.close();

    return

def create_cfradial_file(ncfile):

    user = getpass.getuser()

    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    scantype = dataset.getncattr('scantype');

    dataset.close();

    # --------------------------------
    # Open NetCDF file using arm_pyart
    # --------------------------------
    if (scantype=='RHI'):
        radar = pyart.aux_io.read_camra_rhi(ncfile);
    elif (scantype=='PPI'):
        radar = pyart.aux_io.read_camra_ppi(ncfile);

    # -----------------------------------
    # Write cfradial file using arm_pyart
    # -----------------------------------
    cfradfile=ncfile.replace("raw.nc","cfrad.nc");

    pyart.io.cfradial.write_cfradial(cfradfile, radar, format='NETCDF4',
        time_reference=False)


# ------------------------------------------------------------------------------
# Quicklook generation
# ------------------------------------------------------------------------------
def make_quicklooks(ncfile,figpath):

    user = getpass.getuser()

    Figurename=ncfile.replace(".nc",".png");

    print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    scantype = dataset.getncattr('scantype');

    dataset.close();

    # --------------------------------
    # Open NetCDF file using arm_pyart
    # --------------------------------
    if (scantype=='RHI'):
        radar = pyart.aux_io.read_camra_rhi(ncfile);
    elif (scantype=='PPI'):
        radar = pyart.aux_io.read_camra_ppi(ncfile);

    dealias_data = pyart.correct.dealias_region_based(radar,
                          ref_vel_field=None, interval_splits=3,
                          interval_limits=None, skip_between_rays=100,
                          skip_along_ray=100, centered=True, nyquist_vel=14.90923,
                          check_nyquist_uniform=True, gatefilter=False,
                          rays_wrap_around=None, keep_original=False, set_limits=True,
                          vel_field='VEL_HV')

    radar.add_field('VEL_UHV', dealias_data)

    from matplotlib import rcParams

    # Define paramters from package
    rcParams['axes.labelsize'] = 16
    rcParams['axes.titlesize'] = 16
    rcParams['xtick.labelsize'] = 14
    rcParams['ytick.labelsize'] = 14

    # create a plot of the first and sixth sweeps
    fig = plt.figure(figsize=(25, 35))
    display = pyart.graph.RadarDisplay(radar)

    ax1 = fig.add_subplot(421)
    display.plot('ZED_H', 0, vmin=-10, vmax=60, ax=ax1,cmap='pyart_HomeyerRainbow',colorbar_orient='horizontal')
    plt.grid()
    ax1.set_ylim(0,12)

    ax2 = fig.add_subplot(423)
    display.plot('ZDR', 0, vmin=-5, vmax=5,ax=ax2,cmap='pyart_HomeyerRainbow',colorbar_orient='horizontal')
    plt.grid()
    ax2.set_ylim(0,12)

    ax3 = fig.add_subplot(425)
    display.plot('LDR', 0, vmin=-35, vmax=5,ax=ax3,cmap='pyart_HomeyerRainbow',colorbar_orient='horizontal')
    plt.grid()
    ax3.set_ylim(0,12)

    ax4 = fig.add_subplot(427)
    display.plot('PDP', 0, vmin=-5, vmax=60, ax=ax4,cmap='pyart_HomeyerRainbow',colorbar_orient='horizontal')
    plt.grid()
    ax4.set_ylim(0,12)

    ax5 = fig.add_subplot(422)
    display.plot('VEL_HV', 0, vmin=-15, vmax=15, ax=ax5,cmap='RdYlBu_r',colorbar_orient='horizontal')
    plt.grid()
    ax5.set_ylim(0,12)

    ax6 = fig.add_subplot(424)
    display.plot('VEL_UHV', 0, vmin=-30, vmax=30, ax=ax6,cmap='RdYlBu_r',colorbar_orient='horizontal')
    plt.grid()
    ax6.set_ylim(0,12)

    ax7 = fig.add_subplot(426)
    display.plot('SPW_HV', 0, vmin=0, vmax=5, ax=ax7,cmap='pyart_HomeyerRainbow',colorbar_orient='horizontal')
    plt.grid()
    ax7.set_ylim(0,12)

    ## Playing with Data and Masks
    CXC=radar.fields['CXC']['data']

    new_CXC=np.ma.masked_where(CXC.data<0.0, CXC)

    L = -np.log10(1-ma.sqrt(new_CXC))

    radar.add_field_like('CXC', 'L', L,replace_existing=True)


    ax8 = fig.add_subplot(428)
    display.plot('L', 0, vmin=0, vmax=3, ax=ax8,cmap='pyart_HomeyerRainbow',colorbar_orient='horizontal')
    plt.grid()
    ax8.set_ylim(0,12)

    plt.savefig(os.path.join(figpath,Figurename),dpi=200)

    plt.close()

if __name__ == "__main__":
    import sys, PythonCall
    PythonCall.PythonCall(sys.argv).execute()
