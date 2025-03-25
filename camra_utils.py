#!/usr/bin/env python
# coding: utf-8

# ==========================================================================
# Module for processing raw radar files from Chilbolton Advanced 
# Meteorological Radar (CAMRa)
# Author: Chris Walden, UK Research & Innovation and
#                       National Centre for Atmospheric Science
# Last modified: 20-05-2023
# ==========================================================================

"""Module for processing raw moment data from Chilbolton Advance Meteorological Radar (CAMRa)."""

module_version = 0.3;

import datetime, cftime

import netCDF4 as nc4
import numpy as np

import getpass, socket
import pyart

from pyart.config import FileMetadata, get_fillvalue
from pyart.core.radar import Radar
from pyart.io.common import _test_arguments, make_time_unit_str

import yaml

import os, fnmatch, glob

from io import StringIO

import pandas as pd


# ------------------------------------------------------------------------------
# Import required tools
# ------------------------------------------------------------------------------
import numpy.ma as ma;
import re, sys, getopt, shutil, zipfile, string, pwd, getpass

import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt  #import plotting package


from datetime import tzinfo, time
#from pylab import *


ZhOFFSET=0 #+7 %dB
ZDROFFSET=0 #+0.6 %+0.6 %dB 

SNR_threshold_HH=2.0;  	# for copolar H - default
SNR_threshold_VV=2.0; #3.5; 	# for copolar V - default
SNR_threshold_X=2.0; #3.5; 	# cross polar SNR threshold

SNR_threshold_CXC=100; 	# rho_hv and L : 200 corresponds to bias of 0.995
SNR_threshold_SPW=6; 	# spectral width thresholding


def read_camra_raw(filename, **kwargs):
    """
    Read a netCDF raw file from CAMRa.

    Parameters
    ----------
    filename : str
        Name of raw netCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.
    """
    
    # ---------------------------------------------------
    # The following are not required for a fixed platform
    # ---------------------------------------------------
    rotation = None;
    tilt = None;
    roll = None;
    drift = None;
    heading = None;
    pitch = None;
    georefs_applied = None;
    antenna_transition = None;

    # -------------------------
    # test for non empty kwargs
    # -------------------------
    _test_arguments(kwargs)

    # --------------------------------
    # create metadata retrieval object
    # --------------------------------
    filemetadata = FileMetadata('camra')

    # -----------------
    # Open netCDF4 file
    # -----------------
    DS = nc4.Dataset(filename)

    nrays = len(DS.dimensions["time"]);
    ngates = len(DS.dimensions["range"]);
    nsweeps = 1; # We only have single sweep files 
    
    ncvars = DS.variables;

    # --------------------------------
    # latitude, longitude and altitude
    # --------------------------------
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    latitude["data"] = np.array([ncvars["latitude"][0]], "f8");
    longitude["data"] = np.array([ncvars["longitude"][0]], "f8");
    altitude["data"]  = np.array([ncvars["height"][0]], "f8")

    
    metadata_keymap = {
        "location": "platform",
        "Longitude": "longitude",
        "Latitude": "latitude",
        "Altitude": "altitude",
        "system": "instrument_serial_number",
        "title": "title",
        "institution": "institution",
        "reference": "reference",
    }

 # time, range, fields, metadata, scan_type, latitude, longitude, altitude, altitude_agl,
    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index, sweep_end_ray_index, rays_per_sweep,
    # target_scan_rate, rays_are_indexed, ray_angle_res,
    # azimuth, elevation, gate_x, gate_y, gate_z, gate_longitude, gate_latitude, projection, gate_altitude,
    # scan_rate, antenna_transition, 
    # None rotation, tilt, roll, drift, heading, pitch
    # ?? georefs_applied
    # instrument_parameters
    # radar_calibration

    variables_keymap = {
        "elevation": "elevation",
        "azimuth": "azimuth_angle",
        "prf": "prf",
    }

    fields_keymap = {
        "ZED_H": "DBZ_H",
        "ZED_V": "DBZ_V",
        "ZDR": "ZDR",
        "LDR": "LDR",
        "CXC": "CXC",
        "VEL_H": "VEL_H",
        "VEL_V": "VEL_V",
        "VEL_HV": "VEL_HV",
        "SPW_H": "SPW_H",
        "SPW_V": "SPW_V",
        "SPW_HV": "SPW_HV",
        "DDV": "DDV",
        "PDP": "PHIDP",
        "PHI_H": "PHI_H",
        "PHI_V": "PHI_V",
        "PHI_HV": "PHI_HV",
        "PHI_HD": "PHI_HD",
        "PHI_VD": "PHI_VD",
        "PHI_HVD": "PHI_HVD",
    }


#    float frequency;
#      :chilbolton_standard_name = "frequency";
#      :long_name = "frequency of transmitted radiation";
#      :C_format = "%.4f";
#      :units = "GHz";

#    float antenna_diameter;
#      :chilbolton_standard_name = "antenna_diameter";
#      :long_name = "antenna diameter";
#      :units = "m";

#    float transmit_power;
#      :chilbolton_standard_name = "peak_transmitted_power";
#      :long_name = "peak transmitted power";
#      :units = "W";

#    float clock;
#      :chilbolton_standard_name = "clock";
#      :long_name = "clock input to ISACTRL";
#      :units = "Hz";
#      :clock_divfactor = 2; // int

#  // global attributes:
#  :operator = "rad";
#  :pulses_per_daq_cycle = 128; // int
#  :ADC_channels = 8; // int
#  :delay_clocks = 8; // int
#  :pulses_per_ray = 128; // int
#  :extra_attenuation = 0.0f; // float
#  :radar_constant = 64.7f; // float
#  :receiver_gain = 45.5f; // float
#  :cable_losses = 4.8f; // float

    variables = list(variables_keymap.keys())

    # --------
    # metadata
    # --------
    metadata = filemetadata('metadata')
    for k in ['history', 'title']:
        if k in DS.ncattrs(): 
            metadata[k] = DS.getncattr(k)

    metadata['instrument_name']='ncas-radar-camra-1'

    # ------------------------------------------
    # sweep_start_ray_index, sweep_end_ray_index
    # ------------------------------------------
    sweep_start_ray_index = filemetadata("sweep_start_ray_index")
    sweep_end_ray_index = filemetadata("sweep_end_ray_index")
    sweep_start_ray_index["data"] = np.array([0], dtype="int32")
    sweep_end_ray_index["data"] = np.array([nrays - 1], dtype="int32")

    # ------------
    # sweep number
    # ------------
    sweep_number = filemetadata("sweep_number")
    sweep_number["data"] = np.array([0], dtype="int32")

    # -----------------------
    # sweep_mode, fixed_angle
    # -----------------------
    sweep_modes = {'ppi' : 'ppi', 'rhi' : 'rhi', 'fixed' : 'pointing', 'csp' : 'coplane', 'man' : 'manual_rhi'}

    sweep_mode = filemetadata("sweep_mode")

    scan_type = DS.getncattr('scantype');

    print(f'scantype = {scan_type}')

    sweep_mode["data"] = np.array(1 * [None]);

    for key, value in sweep_modes.items():
        if key in scan_type.lower(): 
            scan_name = value;
            break;


    fixed_angles = {'ppi' : ncvars['elevation'][0], 'rhi' : ncvars['azimuth'][0] % 360, 'pointing' : ncvars['elevation'][0], "manual_rhi" : ncvars['azimuth'][0] % 360}

    fixed_angle = filemetadata("fixed_angle")


    if scan_name is not None:
        fixed_angle["data"] = np.array(1 * [fixed_angles[scan_name]],dtype='f')
    else:
        fixed_angle["data"] = np.array(1 * [None],dtype='f')


    # time
    # ----
    time = filemetadata('time')
    
    dtime = cftime.num2pydate(ncvars['time'][:],ncvars['time'].units)

    base_time = dtime[0].replace(hour=0, minute=0, second=0, microsecond=0)
    
    time['units'] = make_time_unit_str(base_time)  
    time['data']  = nc4.date2num(dtime,time['units']);

    metadata['time_coverage_start'] = datetime.datetime.strftime(dtime[0],'%Y-%m-%dT%H:%M:%SZ');
    metadata['time_coverage_end'] = datetime.datetime.strftime(dtime[-1],'%Y-%m-%dT%H:%M:%SZ');

    # range
    # -----
    _range = filemetadata('range');
    _range['data'] = ncvars['range'][:];
    _range['units'] = 'metres';
    #_range['metres_to_centre_of_first_gate'] = _range['data'][0];
    _range['proposed_standard_name'] = "projection_range_coordinate";
    _range['long_name'] = "distance_to_centre_of_each_range_gate";
    # assuming the distance between all gates is constant, may not
    # always be true.
    #_range['metres_between_gates'] = (_range['data'][1] - _range['data'][0])

    # azimuth, elevation
    # ------------------
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    azimuth['data'] = ncvars['azimuth'][:];
    azimuth['units'] = "degrees";
    azimuth['proposed_standard_name'] = "ray_azimuth_angle";
    azimuth['long_name'] = "azimuth_angle_from_true_north";
        
    elevation['data'] = ncvars['elevation'][:];
    elevation['proposed_standard_name'] = "ray_elevation_angle";
    elevation['long_name'] = "elevation_angle_from_horizontal_plane";
    elevation['units'] = "degrees";

    azimuth_span = max(azimuth['data'])-min(azimuth['data']);

    ray_duration = time['data'][1:]-time['data'][:-1];
    target_ray_duration = np.round(np.mean(ray_duration),3);
    ray_duration = np.insert(ray_duration,0,target_ray_duration);
    print(f'Target ray duration = {target_ray_duration}')
    print(f'Ray duration = {ray_duration[:10]}')
    print(f'Ray duration = {ray_duration[-10:]}')

    if scan_name == 'rhi':
        scan_rate = (elevation['data'][1:]-elevation['data'][:-1])/(time['data'][1:]-time['data'][:-1]);
        target_scan_rate = (elevation['data'][-1]-elevation['data'][0])/(time['data'][-1]-time['data'][0]);
        scan_rate = np.insert(scan_rate,0,target_scan_rate);
        print(f'Target scan rate = {target_scan_rate}')


    elif scan_name == 'ppi':
        scan_rate = (azimuth['data'][1:]-azimuth['data'][:-1])/(time['data'][1:]-time['data'][:-1]);
        target_scan_rate = (azimuth['data'][-1]-azimuth['data'][0])/(time['data'][-1]-time['data'][0]);
        scan_rate = np.insert(scan_rate,0,target_scan_rate);
        print(f'Target scan rate = {target_scan_rate}')
        print(f'Scan rate = {scan_rate[0:10]}')
        print(f'Scan rate = {scan_rate[-10:]}')

    print(f"Az = {azimuth['data'][:10]}")
    print(f"Az = {azimuth['data'][-10:]}")

    # Assume angles are at the end of each ray.  

    duration_offset = -0.05;
    duration_offset_ppi = -0.12;
    if scan_name == 'rhi':
        if elevation['data'][-1]>elevation['data'][0]:
            elevation['data'] -= (0.5*ray_duration+duration_offset) * scan_rate ;
        else:
            elevation['data'] -= (0.5*ray_duration+duration_offset) * scan_rate ;
    elif scan_name == 'ppi':
        azimuth['data'] -= (0.5*ray_duration+duration_offset_ppi) * scan_rate ;
   
    print(f"Az new = {azimuth['data'][:10]}")
    print(f"Az new = {azimuth['data'][-10:]}")

    if scan_name == 'ppi':
        print(azimuth_span);

        if azimuth_span > 350.0:
            sweep_mode["data"] = np.array(1 * ["azimuth_surveillance"]);
            sweep_mode["data"][0] = "azimuth_surveillance";
        else:
            print("It is a sector scan")
            sweep_mode["data"] = np.array(1 * ["sector"]);
            sweep_mode["data"][0] = "sector";
    else:
        sweep_mode["data"] = np.array(1 * [scan_type]);
        sweep_mode["data"][0] = scan_type;
    
    print(sweep_mode['data'][:])

    # Handle azimuths beyond 360
    azimuth['data'] = azimuth['data'] % 360;


    fields = {}

    if "ZED_H" in ncvars:
        field_name = fields_keymap['ZED_H']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dBZ'
        field_dic['data'] = ncvars['ZED_H'][:];
        field_dic['long_name'] =  "radar_equivalent_reflectivity_factor_for_copolar_horizontal_receive_signal";
        field_dic['proposed_standard_name'] =  "radar_equivalent_reflectivity_factor_h";   
        fields[field_name] = field_dic
    else:
        print("ZED_H does not exist")

    if "ZED_V" in ncvars:
        field_name = fields_keymap['ZED_V']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dBZ'
        field_dic['data'] = ncvars['ZED_V'][:];
        field_dic['long_name'] =  "radar_equivalent_reflectivity_factor_for_copolar_vertical_receive_signal";
        field_dic['proposed_standard_name'] =  "radar_equivalent_reflectivity_factor_v";   
        fields[field_name] = field_dic
    else:
        print("ZED_V does not exist")

    lightspeed = 299792458;
    wavelength = lightspeed/(ncvars['frequency'][:]*1.0e9);
    vfold_hv = ncvars['prf'][:]*wavelength/4.0;

    if "VEL_H" in ncvars:
        field_name = fields_keymap['VEL_H']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['VEL_H'][:];
        field_dic['long_name'] =  "radial_velocity_of_scatterers_away_from_instrument_for_copolar_horizontal_receive_signal";
        field_dic['standard_name'] = "radial_velocity_of_scatterers_away_from_instrument_h";
        field_dic['field_folds'] = 'true'
        field_dic['field_limit_lower'] = -vfold_hv/2.0;
        field_dic['field_limit_upper'] = vfold_hv/2.0;
        fields[field_name] = field_dic
    else:
        print("VEL_H does not exist")

    if "VEL_V" in ncvars:
        field_name = fields_keymap['VEL_V']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['VEL_V'][:];
        field_dic['long_name'] =  "radial_velocity_of_scatterers_away_from_instrument_for_copolar_vertical_receive_signal";
        field_dic['standard_name'] = "radial_velocity_of_scatterers_away_from_instrument_v";
        field_dic['field_folds'] = 'true'
        field_dic['field_limit_lower'] = -vfold_hv/2.0;
        field_dic['field_limit_upper'] = vfold_hv/2.0;
        fields[field_name] = field_dic
    else:
        print("VEL_V does not exist")

    if "VEL_HV" in ncvars:
        field_name = fields_keymap['VEL_HV']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['VEL_HV'][:];
        field_dic['long_name'] =  "radial_velocity_of_scatterers_away_from_instrument";
        field_dic['standard_name'] = "radial_velocity_of_scatterers_away_from_instrument";
        field_dic['field_folds'] = 'true'
        field_dic['field_limit_lower'] = -vfold_hv;
        field_dic['field_limit_upper'] = vfold_hv;
        fields[field_name] = field_dic
    else:
        print("VEL_HV does not exist")

    if "ZDR" in ncvars:
        field_name = fields_keymap['ZDR']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dB'
        field_dic['data'] = ncvars['ZDR'][:];
        field_dic['long_name'] =  "radar_differential_reflectivity_hv";
        field_dic['proposed_standard_name'] = "radar_differential_reflectivity_hv";
        fields[field_name] = field_dic
    else:
        print("ZDR does not exist")

    if "LDR" in ncvars:
        field_name = fields_keymap['LDR']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'dB'
        field_dic['data'] = ncvars['LDR'][:];
        field_dic['long_name'] =  "radar_linear_depolarization_ratio";
        field_dic['proposed_standard_name'] = "radar_linear_depolarization_ratio";
        fields[field_name] = field_dic
    else:
        print("LDR does not exist")

    if "SPW_HV" in ncvars:
        field_name = fields_keymap['SPW_HV']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['SPW_HV'][:];
        field_dic['long_name'] =  "radar_doppler_spectrum_width";
        field_dic['proposed_standard_name'] = "radar_doppler_spectrum_width";
        fields[field_name] = field_dic
    else:
        print("SPW_HV does not exist")

    if "SPW_H" in ncvars:
        field_name = fields_keymap['SPW_H']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['SPW_H'][:];
        field_dic['long_name'] =  "radar_doppler_spectrum_width_for_copolar_horizontal_receive_signal";
        field_dic['proposed_standard_name'] = "radar_doppler_spectrum_width_h";
        fields[field_name] = field_dic
    else:
        print("SPW_H does not exist")

    if "SPW_V" in ncvars:
        field_name = fields_keymap['SPW_V']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'm s-1'
        field_dic['data'] = ncvars['SPW_V'][:];
        field_dic['long_name'] =  "radar_doppler_spectrum_width_for_copolar_vertical_receive_signal";
        field_dic['proposed_standard_name'] = "radar_doppler_spectrum_width_v";
        fields[field_name] = field_dic
    else:
        print("SPW_V does not exist")

    if "PDP" in ncvars:
        field_name = fields_keymap['PDP']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PDP'][:];
        field_dic['long_name'] =  "radar_differential_phase_hv";
        field_dic['proposed_standard_name'] = "radar_differential_phase_hv";
        fields[field_name] = field_dic
    else:
        print("PDP does not exist")

    if "CXC" in ncvars:
        field_name = fields_keymap['CXC']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = ''
        field_dic['data'] = ncvars['CXC'][:];
        field_dic['long_name'] =  "radar_copolar_cross_correlation";
        fields[field_name] = field_dic
    else:
        print("CXC does not exist")
    
    if "DDV" in ncvars:
        field_name = fields_keymap['DDV']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = ''
        field_dic['data'] = ncvars['DDV'][:];
        field_dic['long_name'] =  "radar_differential_doppler_velocity";
        fields[field_name] = field_dic
    else:
        print("DDV does not exist")
   
    if "PHI_H" in ncvars:
        field_name = fields_keymap['PHI_H']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PHI_H'][:];
        field_dic['long_name'] =  "radar_absolute_phase_for_copolar_horizontal_receive_signal";
        fields[field_name] = field_dic
    else:
        print("PHI_H does not exist")

    if "PHI_V" in ncvars:
        field_name = fields_keymap['PHI_V']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PHI_V'][:];
        field_dic['long_name'] =  "radar_absolute_phase_for_copolar_vertical_receive_signal";
        fields[field_name] = field_dic
    else:
        print("PHI_V does not exist")

    if "PHI_HV" in ncvars:
        field_name = fields_keymap['PHI_HV']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PHI_HV'][:];
        field_dic['long_name'] =  "radar_absolute_phase";
        fields[field_name] = field_dic
    else:
        print("PHI_HV does not exist")

    if "PHI_HD" in ncvars:
        field_name = fields_keymap['PHI_HD']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PHI_HD'][:];
        field_dic['long_name'] =  "radar_standard_deviation_of_absolute_phase_for_copolar_horizontal_receive_signal";
        fields[field_name] = field_dic
    else:
        print("PHI_HD does not exist")

    if "PHI_VD" in ncvars:
        field_name = fields_keymap['PHI_VD']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PHI_VD'][:];
        field_dic['long_name'] =  "radar_standard_deviation_of_absolute_phase_for_copolar_vertical_receive_signal";
        fields[field_name] = field_dic
    else:
        print("PHI_VD does not exist")

    if "PHI_HVD" in ncvars:
        field_name = fields_keymap['PHI_HVD']
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue();
        field_dic['units'] = 'degree'
        field_dic['data'] = ncvars['PHI_HVD'][:];
        field_dic['long_name'] =  "radar_standard_deviation_of_absolute_phase";
        fields[field_name] = field_dic
    else:
        print("PHI_HVD does not exist")

    # instrument_parameters
    instrument_parameters = {}
    radar_parameters = {}

    sweep_start_ray_index["data"] = np.array([0], dtype="int32")


    if "prf" in ncvars:
        valarray = np.empty_like(elevation['data'])  # Initializes prt with zeros, same shape as elevation['data']  
        valarray[:] = 1.0 / ncvars["prf"][0]
        dic = filemetadata("prt")
        dic["data"] = valarray;
        instrument_parameters["prt"] = dic

    if "pulse_period" in ncvars:
        valarray = np.empty_like(elevation['data'])  # Initializes prt with zeros, same shape as elevation['data']  
        valarray[:] = ncvars["pulse_period"][0]
        dic = filemetadata("pulse_width")
        dic["data"] = valarray
        instrument_parameters["pulse_width"] = dic

    if "beamwidthH" in ncvars:
        dic = filemetadata("radar_beam_width_h")
        dic["data"] = ncvars["beamwidthH"][:]
        instrument_parameters["radar_beam_width_h"] = dic

    if "beamwidthV" in ncvars:
        dic = filemetadata("radar_beam_width_v")
        dic["data"] = ncvars["beamwidthV"][:]
        instrument_parameters["radar_beam_width_v"] = dic
    
    radar_calibration = {}

    DS.close()

    return Radar(
        time,
        _range,
        fields,
        metadata,
        scan_name,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        instrument_parameters=instrument_parameters,
        radar_calibration=radar_calibration,
    )

def multi_camraraw2cfrad(
    rawfiles,
    output_dir,
    scan_name="RHI",
    azimuth_offset=0.0,
    tracking_tag="AMOF_20220922221548",
    campaign="woest",
    data_version="1.0.0",
):
    """
    Aggregates single-sweep CAMRa raw data to a cfradial1 data.
    output_dir(str): Enter the path for output data,
    scan_name(str): "RHI"
    """
    
    from pathlib import Path
    homepath = Path.home()

    yaml_project_file = os.path.join(homepath,'amof_campaigns','{}_project.yml'.format(campaign))
    yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

    out_dir = output_dir
    files = sorted(rawfiles)
    print(files)
    print("Number of files: ", len(files))

    print('Start to read raw CAMRa file')

    RadarDS = read_camra_raw(files[0]);

    print('Done reading raw CAMRa file')

    print("Merging all scans into one Volume")
    for i in range(1, len(files)):

        newRadarDS = read_camra_raw(files[i])

        if 'RHI' in scan_name or 'rhi' in scan_name:
            if np.max(newRadarDS.elevation['data'])-np.min(newRadarDS.elevation['data'])!=0:
                print(f'sweep = {i}');
                RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
        elif 'PPI' in scan_name or 'ppi' in scan_name or 'VAD' in scan_name or 'vad' in scan_name:
            if np.max(newRadarDS.azimuth['data'])-np.min(newRadarDS.azimuth['data'])!=0:
                print(f'sweep = {i}');
                RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
        elif 'pointing' in scan_name:        
            RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
 
    
    fname = os.path.basename(files[0]).split(".")[0]

    out_file = f"{fname}_{scan_name.lower()}.nc"

    print(out_file);
    out_path = os.path.join(out_dir, out_file)
    print(out_path);
    pyart.io.write_cfradial(out_path, RadarDS, format="NETCDF4")

    DS = nc4.Dataset(out_path,'r+');
    DS.scan_name = scan_name.lower();
    DS.close();

    amend_unitless(out_path)
    lowercase_long_names(out_path)
    time_long_name(out_path)

    # Update history
    update_string = 'Merge single sweep files into cfradial file'
    update_history_attribute(out_path,update_string)

    print(rawfiles[0]);
    cfradial_add_instrument_parameters(rawfiles[0],out_path,yaml_project_file,yaml_instrument_file,tracking_tag, data_version)

    #cfradial_add_geometry_correction(out_path,revised_northangle);


    print(out_path);
    print('Going to add NCAS metadata')
    outfile = cfradial_add_ncas_metadata(out_path,yaml_project_file,yaml_instrument_file,tracking_tag,data_version);

    return outfile

def amend_unitless(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');

    # Loop through each variable in the NetCDF file
    for var_name in DS.variables:
        var = DS.variables[var_name]
        if hasattr(var, 'units') and var.units == 'unitless':
            var.units = ""
        if hasattr(var, 'units') and var.units == 'count':
            var.units = ""

    DS.close()

def lowercase_long_names(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');

    # Loop through each variable in the NetCDF file
    for var_name in DS.variables:
        var = DS.variables[var_name]
        if hasattr(var, 'long_name'):
            var.long_name = var.long_name.lower();
            if "utc" in var.long_name:
                var.long_name = var.long_name.replace("utc", "UTC")


    DS.close()


def update_long_name_attributes(cfradfile):  
    # Open the NetCDF file in read/write mode  
    with nc4.Dataset(cfradfile, 'r+') as dataset:  
        # Loop through all variables in the dataset  
        for var_name in dataset.variables:  
            var = dataset.variables[var_name]  
            # Check if the variable has a 'long_name' attribute  
            if 'long_name' in var.ncattrs():  
                # Get the current long_name  
                long_name = var.long_name  
                # Replace white spaces with underscores  
                updated_long_name = long_name.replace(' ', '_')  
                #updated_long_name = re.sub(r'(?<!\S) (?!\S)', '_', long_name)  
                # Convert to lower case  
                updated_long_name = updated_long_name.lower()  
                # Assign the updated long_name back to the variable  
                var.long_name = updated_long_name  
                print(f'Updated long_name for variable "{var_name}": {updated_long_name}')  
        dataset['azimuth'].long_name = 'azimuth_angle_from_true_north'
        dataset['elevation'].long_name = 'elevation_angle_from_horizontal_plane'



def update_geospatial_bounds(cfradfile):
    # Open the NetCDF file in read-write mode
    with nc4.Dataset(cfradfile, 'r+') as DS:

        updated = False
        # Check if the 'geospatial_bounds' attribute exists
        if 'geospatial_bounds' in DS.ncattrs():
            # Get the attribute value
            geospatial_bounds = DS.getncattr('geospatial_bounds')
            # Look for "Bounding Box:" in the attribute
            if geospatial_bounds.startswith("Bounding box:"):
                # Extract the coordinate pairs
                pairs = geospatial_bounds.replace("Bounding box:", "").strip().split(",")
                if len(pairs) == 2 and pairs[0].strip() == pairs[1].strip():
                    # If the first and second pairs are the same, update the attribute
                    new_value = pairs[0].strip()
                    DS.setncattr('geospatial_bounds', new_value)
                    print(f"Updated 'geospatial_bounds' to: {new_value}")
                    DS.featureType = "timeSeriesProfile"
                    DS.scan_name = 'vpt'
                    updated = True

    if updated:
        update_history_attribute(cfradfile,'Updated geospatial bounds')


def update_field_names(cfradfile):
    """                                                                                                                                                                      
    Replace the field_names global attribute with the complete list of fields                                                                                                
    found in the NetCDF file.                                                                                                                                                
                                                                                                                                                                             
    :param file_path: Path to the CfRadial NetCDF file.                                                                                                                      
    """
    with nc4.Dataset(cfradfile, 'r+') as DS:
        # Identify field variables based on their dimensions (time and range)                                                                                                
        field_variables = [
            var_name for var_name, var in DS.variables.items()
            if 'time' in var.dimensions and 'range' in var.dimensions
        ]

        # Create the new field_names value as a comma-separated string                                                                                                       
        new_field_names = ','.join(field_variables)

        # Update the field_names global attribute                                                                                                                            
        if 'field_names' in DS.ncattrs():
            print(f"Updating 'field_names' from: {DS.getncattr('field_names')}")
        else:
            print("Adding new 'field_names' attribute.")

        DS.setncattr('field_names', new_field_names)
        print(f"'field_names' updated to: {new_field_names}")

    update_history_attribute(cfradfile,'Update field_names attribute')


def lowercase_platform(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');

    val = DS.platform;
    DS.platform= val.lower()

    DS.close()

def fix_platform_altitude(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');
    DS.platform_altitude = "84 m (orthometric height above EGM2008 geoid)"

def time_long_name(cfradfile):
    DS = nc4.Dataset(cfradfile,'r+');
    time_var = DS.variables['time'];
    if 'time_reference' in DS.variables:
        time_var.long_name = "time_since_time_reference"

    DS.close()


def cfradial_get_bbox(cfradfile):
    print(cfradfile)
    Radar = pyart.io.read_cfradial(cfradfile);
    latmin = np.min(Radar.gate_latitude['data']);
    lonmin = np.min(Radar.gate_longitude['data']);
    latmax = np.max(Radar.gate_latitude['data']);
    lonmax = np.max(Radar.gate_longitude['data']); 
    print(latmin,latmax,lonmin,lonmax)
    boundingbox = f"Bounding box: {latmin:.4f}N {lonmin:.4f}E, {latmax:.4f}N {lonmax:.4f}E"
    return boundingbox

def cfradial_add_ncas_metadata(cfradfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):
    
    instrument_tagname = "ncas-radar-camra-1"

    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    with open(yaml_instrument_file, "r") as stream:
        try:
            instruments = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname];

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    print(tracking_tag);
    print(yaml_project_file);
    print(yaml_instrument_file);

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];

    radar_name = instrument["instrument_name"].lower();

    print(radar_name);

    for n in project["ncas_instruments"]:
        if radar_name in n:
            project_instrument = n[radar_name];

    print(project_instrument);

    location = instrument['platform']['location'].lower();
    
    RadarDataset = nc4.Dataset(cfradfile);

    scan_name = RadarDataset.scan_name;

    print(scan_name);
    print('HERE')

    str_start = RadarDataset.variables['time_coverage_start'][:].tobytes().decode('utf-8')

    time_coverage_start = str_start[0:20];

    print('and HERE')

    str_end = RadarDataset.variables['time_coverage_end'][:].tobytes().decode('utf-8').strip()
    time_coverage_end = str_end[0:20];

    print('and now HERE')
    print(time_coverage_start)
    print(type(time_coverage_start))   
    print(len(time_coverage_start))   
    file_timestamp = datetime.datetime.strptime(time_coverage_start,'%Y-%m-%dT%H:%M:%SZ');

    print('and then HERE')

    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')

    
    import pathlib
    outpath = pathlib.Path(cfradfile).parent.resolve();

    outfile = os.path.join(outpath,'{}_{}_{}_{}_l1_v{}.nc'.format(radar_name,location,dtstr,scan_name.replace('_','-',1),data_version));

    if os.path.isfile(outfile):
        print("The file already exists")
    else:
        # Rename the file
        os.rename(cfradfile,outfile);
    
    RadarDataset.close();
    
    # -------------------------------------------------------
    # Read cfradial file to add NCAS metadata
    # -------------------------------------------------------
    DS = nc4.Dataset(outfile,'r+');

    DS.Conventions = "NCAS-Radar-1.0 CfRadial-1.4 instrument_parameters radar_parameters geometry_correction"

    if 'version' in DS.ncattrs():
        DS.delncattr('version');

    DS.product_version = f"v{data_version}";
    DS.processing_level = "1" ;

    DS.licence = project_instrument["data_licence"];
    DS.acknowledgement = project_instrument["acknowledgement"];

    DS.platform = instrument["platform"]["location"];
    DS.platform_type = instrument["platform"]["type"];
    DS.location_keywords = instrument["platform"]["location_keywords"];
    DS.platform_is_mobile = "false";
    DS.deployment_mode = instrument["platform"]["deployment_mode"];

    DS.platform_altitude = instrument["platform"]["altitude"];

    DS.title = project_instrument["title"];
    DS.source = project_instrument["source"];

    DS.creator_name = project_instrument["data_creator"]["name"];
    DS.creator_email = project_instrument["data_creator"]["email"];
    DS.creator_url = project_instrument["data_creator"]["pid"];
    DS.institution = project_instrument["data_creator"]["institution"];
    DS.instrument_name = instrument["instrument_name"];
    DS.instrument_software = project_instrument["instrument_software"]["name"];
    DS.instrument_software_version = project_instrument["instrument_software"]["version"];
    DS.instrument_manufacturer = instrument['instrument_manufacturer'];
    DS.instrument_model = instrument['instrument_model'];
    DS.instrument_serial_number = instrument['instrument_serial_number'];
    DS.instrument_pid = instrument['instrument_pid']

    DS.references = instrument['references'];
    DS.comment = " ";
    DS.project = project["project_name"];
    DS.project_principal_investigator = project["principal_investigator"]["name"];
    DS.project_principal_investigator_email = project["principal_investigator"]["email"];
    DS.project_principal_investigator_url = project["principal_investigator"]["pid"];

    DS.processing_software_url = "https://github.com/longlostjames/camra-radar-utils/releases/tag/v1.0.0";
    DS.processing_software_version = "v1.0.0";

    DS.time_coverage_start = time_coverage_start;
    DS.time_coverage_end = time_coverage_end;

    print('getting bbox')
    print(cfradfile)
    str = cfradial_get_bbox(outfile)
    print(str)
    DS.geospatial_bounds = str;


    if "vpt" in DS.scan_name or "VPT" in DS.scan_name or "vertical_pointing" in DS.scan_name:
        DS.featureType = 'timeSeriesProfile';

    # -------------------------------------------------------
    # Now clean up some variable attributes
    # -------------------------------------------------------
    DS['time'].comment = "";
    try:
        DS['range'].delncattr('standard_name');
    except: 
        pass
    
    DS['range'].comment = 'Range to centre of each bin';
    DS['range'].meters_to_center_of_first_gate = DS['range'][0];
    try:
        DS['azimuth'].delncattr('standard_name');
    except:
        pass
    try:
        DS['elevation'].delncattr('standard_name');
    except:
        pass
    #DS['DBZ'].standard_name = 'equivalent_reflectivity_factor';
    try:
        DS['sweep_number'].delncattr('standard_name');
    except:
        pass
    try:
        DS['sweep_mode'].delncattr('standard_name');
    except:
        pass
    try:
        DS['fixed_angle'].delncattr('standard_name');
    except:
        pass


    # Set the variable to a sequence of integers starting from 0  
    DS['sweep_number'][:] = np.arange(len(DS['sweep_number'][:]))  

    DS['latitude'].long_name = 'latitude';
    DS['latitude'].standard_name = 'latitude';
    DS['latitude'][:] = float(instrument["latitude"])
    DS['longitude'].long_name = 'longitude';
    DS['longitude'].standard_name = 'longitude'; 
    DS['longitude'][:] = float(instrument["longitude"])

    if DS['longitude'][0]<0:
        DS['longitude'][:]+=360.0;

    DS['altitude'].standard_name = 'altitude';
    DS['altitude'].comment = 'Altitude of the centre of rotation of the antenna above the geoid using the WGS84 ellipsoid and EGM2008 geoid model' 
    DS['altitude'].long_name = 'altitude';
    DS['altitude'].units = 'metres';
    DS['altitude'][:] = float(instrument["altitude"]["value"])
    try:
        DS['altitude'].delncattr('positive'); 
    except:
        pass
    DS['volume_number'].long_name = 'data_volume_index_number';
    DS['volume_number'].units = "" ;
    #DS['volume_number']._FillValue = -9999 ;

    altitude_agl = DS.createVariable('altitude_agl', 'f8')
    altitude_agl.assignValue(float(instrument["altitude_agl"]["value"]))
    altitude_agl.long_name = 'altitude_above_ground_level';
    altitude_agl.units = instrument['altitude_agl']['units'];

    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['latitude'];
    #varout = DSout.createVariable('latitude',varin.datatype);
    #varout.standard_name = 'latitude';
    #varout.long_name = 'latitude of the antenna';
    #varout.units = 'degree_north';
    #varout[:]=51.1450;

    #varin = DSin['longitude'];
    #varout = DSout.createVariable('longitude',varin.datatype);
    #varout.standard_name = 'longitude';
    #varout.long_name = 'longitude of the antenna';
    #varout.units = 'degree_east';
    #varout[:]=-1.4384;

    #varin = DSin['height'];
    #varout = DSout.createVariable('altitude',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above the geoid (WGS84)';
    #varout.units = 'm';
    #varout[:]=146.7;

    #varout = DSout.createVariable('altitude_agl',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above ground';
    #varout.units = 'm';
    #varout[:]=16.0;

    #varin = DSin['frequency'];
    #varout = DSout.createVariable('frequency',varin.datatype);
    #varout.standard_name = 'radiation_frequency';
    #varout.long_name = 'frequency of transmitted radiation';
    #varout.units = 'GHz';
    #varout[:]=varin[:];

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];


    # ---------------
    # Field variables
    # ---------------
    # These are:
    # SNR, VEL, RMS, LDR, NPK, SNRg, VELg, NPKg, RHO, DPS, RHOwav, DPSwav, 
    # HSDco, HSDcx, Ze, Zg, ISDRco, ISDRcx
    # The ones to use are:
    # SNR, VEL, RMS, LDR, NPK, RHO, DPS, HSDco, HSDcx, Ze, ISDRco, ISDRcx

    DS.close();

    # -----------------------
    # Update history metadata
    # -----------------------
    updatestr = "Add NCAS metadata"
    update_history_attribute(outfile,updatestr)

    return outfile

# ===================
# CONVERSION ROUTINES
# ===================


def convert_camra_raw2l0b(infile,outpath,yaml_project_file,yaml_instrument_file,tracking_tag):

    """This routine converts L0 raw moment data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 0b (cfradial) data.
    Metadata are added using information in two YAML files the yaml_project_file, and yaml_instrument_file.

    :param infile: Full path of NetCDF Level 0 raw data file
    :type infile: str

    :param outfile: Full path of NetCDF Level 1 output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20220907-071502_ppi_l1_v1.0.nc`
    :type outfile: str
    """

    data_version = '1.0';

    instrument_tagname = "ncas-radar-camra-1"

    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    with open(yaml_instrument_file, "r") as stream:
        try:
            instruments = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname];

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];

    radar_name = instrument["instrument_name"].lower();

    print(radar_name);

    for n in project["ncas_instruments"]:
        if radar_name in n:
            project_instrument = n[radar_name];

    print(project_instrument);

    location = project_instrument['platform']['location'].lower();
    
    RadarDataset = read_camra_raw(infile);

    scan_type = RadarDataset.scan_type

    file_timestamp = datetime.datetime.strptime(RadarDataset.metadata["time_coverage_start"],'%Y-%m-%dT%H:%M:%SZ');

    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')

    outfile = os.path.join(outpath,'{}_{}_{}_{}_l0b_v{}.nc'.format(radar_name,location,dtstr,scan_type.replace('_','-',1),data_version));

    # Use PyART to create CfRadial file
    pyart.io.write_cfradial(outfile, RadarDataset, format='NETCDF4', time_reference=True)

    # -------------------------------------------------------
    # Read freshly created cfradial file to add NCAS metadata
    # -------------------------------------------------------
    DS = nc4.Dataset(outfile,'r+');

    DS.product_version = "v{}".format(data_version) ;
    DS.processing_level = "1" ;

    DS.licence = project_instrument["data_licence"];
    DS.acknowledgement = project_instrument["acknowledgement"];

    DS.platform = instrument["platform"]["location"];
    DS.platform_type = instrument["platform"]["type"];
    DS.location_keywords = instrument["platform"]["location_keywords"];
    DS.platform_altitude = instrument["platform"]["altitude"];

    DS.title = project_instrument["title"];

    DS.creator_name = project_instrument["data_creator"]["name"];
    DS.creator_email = project_instrument["data_creator"]["email"];
    DS.creator_url = project_instrument["data_creator"]["pid"];
    DS.institution = project_instrument["data_creator"]["institution"];
    DS.instrument_name = instrument["instrument_name"];
    DS.instrument_software = project_instrument["instrument_software"]["name"];
    DS.instrument_software_version = project_instrument["instrument_software"]["version"];
    DS.instrument_manufacturer = instrument['instrument_manufacturer'];
    DS.instrument_model = instrument['instrument_model'];
    DS.instrument_serial_number = instrument['instrument_serial_number'];

    DS.references = instrument['references'];
    DS.source = "CAMRa";
    DS.comment = " ";
    DS.project = project["project_name"];
    DS.project_principal_investigator = project["principal_investigator"]["name"];
    DS.project_principal_investigator_email = project["principal_investigator"]["email"];
    DS.project_principal_investigator_url = project["principal_investigator"]["pid"];

    DS.processing_software_url = " ";
    DS.processing_software_version = " ";

    #DS.time_coverage_start = datetime.datetime.strftime(dt_start,'%Y-%m-%dT%H:%M:%SZ');
    #DS.time_coverage_end = datetime.datetime.strftime(dt_end,'%Y-%m-%dT%H:%M:%SZ');
    #DS.geospatial_bounds = "51.1450N -1.4384E";

    
    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['latitude'];
    #varout = DSout.createVariable('latitude',varin.datatype);
    #varout.standard_name = 'latitude';
    #varout.long_name = 'latitude of the antenna';
    #varout.units = 'degree_north';
    #varout[:]=51.1450;

    #varin = DSin['longitude'];
    #varout = DSout.createVariable('longitude',varin.datatype);
    #varout.standard_name = 'longitude';
    #varout.long_name = 'longitude of the antenna';
    #varout.units = 'degree_east';
    #varout[:]=-1.4384;

    #varin = DSin['height'];
    #varout = DSout.createVariable('altitude',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above the geoid (WGS84)';
    #varout.units = 'm';
    #varout[:]=146.7;

    #varout = DSout.createVariable('altitude_agl',varin.datatype);
    #varout.standard_name = 'altitude';
    #varout.long_name = 'altitude of the elevation axis above ground';
    #varout.units = 'm';
    #varout[:]=16.0;

    #varin = DSin['frequency'];
    #varout = DSout.createVariable('frequency',varin.datatype);
    #varout.standard_name = 'radiation_frequency';
    #varout.long_name = 'frequency of transmitted radiation';
    #varout.units = 'GHz';
    #varout[:]=varin[:];

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];



    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()

    updttime = datetime.datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: camra_utils.convert_camra_raw2l1"
    + " version:" + str(module_version));

    DS.history = history + "\n" + DS.history;

    DS.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DS.close();

    return


def convert_camra_ts_l0a2l0b(infile,output_dir,
    scan_name="VPT",
    tracking_tag="AMOF_20230201132601",
    campaign="ccrest-m",
    data_version="1.0.0"):

    """This routine converts raw (Level 0a) time series data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 0b data.
    Processing involves removing redundant dimensions, and removing bias from the ADC samples of transmit and receive I and Q.
    Single estimates of I and Q for each transmit pulse are produced and stored in the output file.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of NetCDF Level 0a raw data file, e.g. `<path-to-file>/radar-camra_20201210212823_fix-ts.nc`
    :type infile: str

    :param outfile: Full path of NetCDF Level 0b output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l0b_v1.0.nc`
    :type outfile: str
    """

    from pathlib import Path
    homepath = Path.home()

    yaml_project_file = os.path.join(homepath,'amof_campaigns','{}_project.yml'.format(campaign))
    yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

    instrument_tagname = "ncas-radar-camra-1"

    print(yaml_instrument_file)
    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    with open(yaml_instrument_file, "r") as stream:
        try:
            instruments = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    
    for elem in instruments:
        if instrument_tagname in elem:
            instrument = elem[instrument_tagname];
    
    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag];
    
    radar_name = instrument["instrument_name"].lower();

    for n in project["ncas_instruments"]:
        print(n)
        if radar_name in n:
            project_instrument = n[radar_name];


    location = project_instrument['platform']['location'].lower();
    
    scan_type = 'vpt';

    print(infile) 

    DSin = nc4.Dataset(infile);

    dt_start = cftime.num2pydate(DSin['time'][0],DSin['time'].units)
    dt_end   = cftime.num2pydate(DSin['time'][-1],DSin['time'].units)

    print(dt_start)
    print(dt_end)
    dtstr = dt_start.strftime('%Y%m%d-%H%M%S')

    print(dtstr);
    print(output_dir)

    outfile = os.path.join(output_dir,f"{radar_name}_{location}_{dtstr}_{scan_type}-ts_l0b_v{data_version}.nc");

    print(outfile)

    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass

    print('Going to create');

    DSout = nc4.Dataset(outfile,mode='w',format='NETCDF4')
    print("Creating {}".format(outfile));

    # ------------------------
    # Set up global attributes
    # ------------------------

    DSout.product_version = f"v{data_version}";
    DSout.processing_level = "0b" ;

    DSout.licence = project_instrument["data_licence"];
    DSout.acknowledgement = project_instrument["acknowledgement"];

    DSout.platform = instrument["platform"]["location"];
    DSout.platform_type = instrument["platform"]["type"];
    DSout.platform_is_mobile = "false";
    DSout.deployment_mode = instrument["platform"]["deployment_mode"];
    DSout.location_keywords = instrument["platform"]["location_keywords"];
    DSout.platform_altitude = instrument["platform"]["altitude"];

    DSout.title = project_instrument["title"];
    DSout.source = project_instrument["source"];


    DSout.creator_name = project_instrument["data_creator"]["name"];
    DSout.creator_email = project_instrument["data_creator"]["email"];
    DSout.creator_url = project_instrument["data_creator"]["pid"];
    DSout.institution = project_instrument["data_creator"]["institution"];
    DSout.instrument_name = instrument["instrument_name"];
    DSout.instrument_software = project_instrument["instrument_software"]["name"];
    DSout.instrument_software_version = project_instrument["instrument_software"]["version"];
    DSout.instrument_manufacturer = instrument['instrument_manufacturer'];
    DSout.instrument_model = instrument['instrument_model'];
    DSout.instrument_serial_number = instrument['instrument_serial_number'];
    DSout.instrument_pid = instrument['instrument_pid']


    DSout.references = instrument['references'];
    DSout.comment = " ";
    DSout.project = project["project_name"];
    DSout.project_principal_investigator = project["principal_investigator"]["name"];
    DSout.project_principal_investigator_email = project["principal_investigator"]["email"];
    DSout.project_principal_investigator_url = project["principal_investigator"]["pid"];

    DSout.processing_software_url = "https://github.com/longlostjames/camra-radar-utils/releases/tag/v1.0.2";
    DSout.processing_software_version = "v1.0.2";

    DSout.scantype = "vertical_pointing";

    DSout.time_coverage_start = dt_start.strftime('%Y-%m-%dT%H:%M:%SZ');
    DSout.time_coverage_end = dt_end.strftime('%Y-%m-%dT%H:%M:%SZ');

    print("done time coverage")
    DSout.geospatial_bounds = "51.1450N -1.4384E";


    # ----------------
    # Scalar variables
    # ----------------

    varin = DSin['latitude'];
    varout = DSout.createVariable('latitude',"f8");
    varout.standard_name = 'latitude';
    varout.long_name = 'latitude of the antenna';
    varout.units = 'degrees_north';
    varout[:]=51.1450;

    print('done lat')
    varin = DSin['longitude'];
    varout = DSout.createVariable('longitude',"f8");
    varout.standard_name = 'longitude';
    varout.long_name = 'longitude of the antenna';
    varout.units = 'degrees_east';
    varout[:]=-1.4384;


    varout = DSout.createVariable('altitude_agl','float');
    varout.standard_name = 'altitude';
    varout.long_name = 'altitude of the elevation axis above ground';
    varout.units = 'm';
    varout[:]=16.0;

    varin = DSin['frequency'];
    varout = DSout.createVariable('frequency',varin.datatype);
    varout.standard_name = 'radiation_frequency';
    varout.long_name = 'frequency of transmitted radiation';
    varout.units = 'GHz';
    varout[:]=varin[:];

    varin = DSin['prf'];
    varout = DSout.createVariable('prf',varin.datatype);
    varout.long_name = 'pulse repetition frequency';
    varout.units = 'Hz';
    varout[:]=varin[:];

    varin = DSin['beamwidthH'];
    varout = DSout.createVariable('beamwidthH',varin.datatype);
    varout.long_name = 'horizontal angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];

    varin = DSin['beamwidthV'];
    varout = DSout.createVariable('beamwidthV',varin.datatype);
    varout.long_name = 'vertical angular beamwidth';
    varout.units = 'degree';
    varout[:]=varin[:];

    varin = DSin['antenna_diameter'];
    varout = DSout.createVariable('antenna_diameter',varin.datatype);
    varout.long_name = 'antenna diameter';
    varout.units = 'm';
    varout[:]=varin[:];

    varout = DSout.createVariable('antenna_focal_length','f4');
    varout.long_name = 'focal length of antenna';
    varout.units = 'm';
    varout[:] = 9.0;

    varout = DSout.createVariable('antenna_focus_radial_location','f4');
    varout.long_name = 'distance along boresight from elevation axis to antenna focus';
    varout.units = 'm';
    varout[:] = 14.34;

    varin = DSin['pulse_period'];
    varout = DSout.createVariable('pulse_width',varin.datatype);
    varout.long_name = 'pulse width';
    varout.units = 'us';
    varout[:]=varin[:];

    varin = DSin['transmit_power'];
    varout = DSout.createVariable('transmit_power',varin.datatype);
    varout.long_name = 'peak transmitted power';
    varout.units = 'W';
    varout[:]=varin[:];

    varin = DSin['clock'];
    varout = DSout.createVariable('clock',varin.datatype);
    varout.long_name = 'clock input to timer card';
    varout.units = 'Hz';
    varout[:]=varin[:];

    varout = DSout.createVariable('clock_divfactor','i4');
    varout.long_name = 'clock divide factor';
    varout.units = '1';
    varout[:]=DSin['clock'].clock_divfactor;

    varout = DSout.createVariable('delay_clocks','i4');
    varout.long_name = 'clock cycles before sampling is initiated';
    varout.units = '1';
    varout[:]=DSin.delay_clocks;

    varout = DSout.createVariable('samples_per_pulse','i4');
    varout.long_name = 'number of samples per pulse';
    varout.units = '1';
    varout[:]=DSin.samples_per_pulse;

    varout = DSout.createVariable('pulses_per_daq_cycle','i4');
    varout.long_name = 'number of pulses per data acquisition cycle';
    varout.units = '1';
    varout[:]=DSin.pulses_per_daq_cycle;

    varout = DSout.createVariable('pulses_per_ray','i4');
    varout.long_name = 'number of pulses per ray';
    varout.units = '1';
    varout[:]=DSin.pulses_per_ray;

    varout = DSout.createVariable('radar_constant','f4');
    varout.long_name = 'radar constant';
    varout.units = 'dB';
    varout[:]=DSin.radar_constant;

    varout = DSout.createVariable('receiver_gain','f4');
    varout.long_name = 'receiver gain';
    varout.units = 'dB';
    varout[:]=DSin.receiver_gain;

    varout = DSout.createVariable('cable_losses','f4');
    varout.long_name = 'cable losses';
    varout.units = 'dB';
    varout[:]=DSin.cable_losses;

    varout = DSout.createVariable('extra_attenuation','f4');
    varout.long_name = 'extra attenuation';
    varout.units = 'dB';
    varout[:]=DSin.extra_attenuation;


    # ---------------
    # Copy dimensions
    # ---------------
    the_dim = DSin.dimensions['time'];
    DSout.createDimension('time', len(the_dim) if not the_dim.isunlimited() else None)

    the_dim = DSin.dimensions['pulses'];
    DSout.createDimension('pulse', len(the_dim) if not the_dim.isunlimited() else None)

    the_dim = DSin.dimensions['samples'];
    DSout.createDimension('range', len(the_dim) if not the_dim.isunlimited() else None)

    # --------------------
    # Coordinate variables
    # --------------------
    varin = DSin['time'];
    varout = DSout.createVariable('time',varin.datatype,('time'));
    varout.standard_name = 'time';
    varout.long_name = 'time at the end of each recorded ray';
    varout.units = varin.units;
    varout[:]=varin[:];

    varin = DSin['range'];
    varout = DSout.createVariable('range',varin.datatype,('range'));
    varout.long_name = varin.long_name;
    varout.range_offset_applied = np.float32(varin.range_offset);
    varout.units = varin.units;
    varout[:]=varin[:];

    # --------------------------
    # Antenna pointing variables
    # --------------------------
    varin = DSin['elevation'];
    varout = DSout.createVariable('elevation',varin.datatype,'time');
    varout.long_name = 'elevation angle of the antenna boresight above the horizon';
    varout.units = 'degree';
    varout.elevation_offset_applied = np.float32(0.);
    varout[:] = varin[:];

    varin = DSin['azimuth'];
    varout = DSout.createVariable('azimuth',varin.datatype,'time');
    varout.long_name = 'azimuth angle of the antenna boresight clockwise from the grid north';
    varout.comment = 'More generally this is the azimuth angle of the plane perpendicular to the elevation axis, which remains defined even when the elevation is 90 degree';
    varout.units = 'degree';
    varout.azimuth_offset_applied = np.float32(0.);
    varout[:] = varin[:];

    # ---------------
    # Field variables
    # ---------------
    varin = DSin['ZLO'];
    varout = DSout.createVariable('ZLO',varin.datatype,('time','pulse','range'),zlib=True);
    varout.long_name = 'ZLO log amplifier output (channel with +12dB gain)';
    varout.units = 'dB';
    varout[:]=varin[:];
    comment_string  = "This is an estimator for co-polar radar equivalent reflectivity factor.\n"
    comment_string += "It does not take into account range correction, the radar constant, receiver gain or cable losses.\n"
    comment_string += "The data are packed and only values in the range [0,3840] (equivalent to the actual_range when the data are unpacked) should be used."
    varout.comment = comment_string;
    varout.scale_factor = np.float32(0.015625);
    varout.add_offset = np.float32(-70.);
    varout.actual_range = [np.float32(-70.),np.float32(-10.)];

    varin = DSin['ZHI'];
    varout = DSout.createVariable('ZHI',varin.datatype,('time','pulse','range'),zlib=True);
    varout.long_name = 'ZHI log amplifier output (channel with -20dB attenuation)';
    varout.units = 'dB';
    varout[:] = varin[:];
    comment_string  = "This is an estimator for co-polar radar equivalent reflectivity factor.\n"
    comment_string += "It does not take into account range correction, the radar constant, receiver gain or cable losses.\n"
    comment_string += "The data are packed and only values in the range [1793,4095] (equivalent to the actual_range when the data are unpacked) should be used."
    varout.comment = comment_string;
    varout.scale_factor = np.float32(0.015625);
    varout.add_offset = np.float32(-38.);
    varout.actual_range = [np.float32(-9.984375), np.float32(25.984375)];

    varin = DSin['ZCX'];
    varout = DSout.createVariable('ZCX',varin.datatype,('time','pulse','range'),zlib=True);
    varout.long_name = 'cross-polar log amplifier output';
    varout.units = 'dB';
    varout[:]=varin[:];
    comment_string  = "This is an estimator for cross-polar radar equivalent reflectivity factor.\n"
    comment_string += "It does not take into account range correction, the radar constant, receiver gain or cable losses."
    varout.comment = comment_string;
    varout.scale_factor = np.float32(0.03125);
    varout.add_offset = np.float32(-77.);

    # -------------------------------------------------------
    # Determine bias-corrected I and Q of each transmit pulse
    # -------------------------------------------------------
    delay_clocks = DSout['delay_clocks'][:];
    clock_divfactor = DSout['clock_divfactor'][:];

    pre_tx    = (18   - delay_clocks) // clock_divfactor;
    tx_pulse  = (24   - delay_clocks) // clock_divfactor;
    post_tx   = (68   - delay_clocks) // clock_divfactor;
    hold_end  = (4708 - delay_clocks) // clock_divfactor;
    post_hold = (4748 - delay_clocks) // clock_divfactor;

    ITXin = DSin['ITX'][:,:,:];
    QTXin = DSin['QTX'][:,:,:];

    ITX_bias_by_gate = np.mean(ITXin,axis=1);
    QTX_bias_by_gate = np.mean(QTXin,axis=1);

    ITXnew = ITXin - ITX_bias_by_gate[:,None,:];
    QTXnew = QTXin - QTX_bias_by_gate[:,None,:];

    # Only use data while sample and hold is active
    # ---------------------------------------------
    sample_end = min([hold_end,len(DSin.dimensions['samples'])]);

    ITXout = np.mean(ITXnew[:,:,post_tx:sample_end],axis=2);
    QTXout = np.mean(QTXnew[:,:,post_tx:sample_end],axis=2);

    varout = DSout.createVariable('ITX','f4',('time','pulse'),zlib=True);
#    add_offset = np.min(ITXout[:]);
#    scale_factor = (np.max(ITXout[:])-np.min(ITXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((ITXout[:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of in-phase video signal for each transmitted pulse';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=ITXout;

    varout = DSout.createVariable('QTX','f4',('time','pulse'),zlib=True);
#    add_offset = np.min(QTXout[:]);
#    scale_factor = (np.max(QTXout[:])-np.min(QTXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((QTXout[:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of quadrature video signal for each transmitted pulse';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=QTXout;

    # ----------------------------------------
    # Determine bias-corrected receive I and Q
    # ----------------------------------------
    IRXin = DSin['IRX'][:,:,:];
    QRXin = DSin['QRX'][:,:,:];
    IRX_bias_by_gate = np.mean(IRXin,axis=1);
    QRX_bias_by_gate = np.mean(QRXin,axis=1);
    IRXout = IRXin - IRX_bias_by_gate[:,None,:];
    QRXout = QRXin - QRX_bias_by_gate[:,None,:];

    varout = DSout.createVariable('IRX','f4',('time','pulse','range'),zlib=True);
#    add_offset = np.min(IRXout[:]);
#    scale_factor = (np.max(IRXout[:])-np.min(IRXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((IRXout[:,:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of received in-phase video signal at output of IF limiting amplifier';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=IRXout;

    varin = DSin['QRX'];
    varout = DSout.createVariable('QRX','f4',('time','pulse','range'),zlib=True);
#    add_offset = np.min(QRXout[:]);
#    scale_factor = (np.max(QRXout[:])-np.min(QRXout[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((QRXout[:,:,:] - add_offset)/scale_factor));
#    varout.scale_factor = np.float32(scale_factor);
#    varout.add_offset = np.float32(add_offset);
    varout.long_name = 'bias-corrected samples of received quadrature video signal at output of IF limiting amplifier';
    varout.units = '1';
#    varout[:] = packed_data;
    varout[:]=QRXout;

    print('updating history')
    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()

    print('OK')
    updttime = datetime.datetime.utcnow()
    print('OK1')
    updttimestr = updttime.ctime()
    print('OK2')
    print(updttimestr)

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: camra_utils.convert_camra_ts_l0a2l0b"
    + " version:" + str(module_version));


    DSout.history = history + "\n" + DSin.history;

    DSout.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DSin.close();
    DSout.close();

    return





def convert_camra_tsl1tomoments(ts_l1path,mom_l1a_path,mom_l1b_path):

    ts_l1path = os.path.join(ts_l1path);
    ts_l1files = glob.glob(f"{ts_l1path}/*.nc");

    sweep = 0

    vpt_ts_file = ts_l1files[sweep];

    print(mom_l1a_path)
    vptfile = glob.glob(f"{mom_l1a_path}/*.nc")[0];

    radar_vpt = pyart.io.read(vptfile);


    radar_vpt.add_field_like("DBZ_H", "DBZ_HS",radar_vpt.fields["DBZ_H"]['data'].copy())
    radar_vpt.fields['DBZ_HS']['long_name'] = 'radar_equivalent_reflectivity_factor_for_copolar_horizontal_receive_signal_from_spectral_processing'
    
    radar_vpt.add_field_like("SNR_H", "SNR_HS",radar_vpt.fields["SNR_H"]['data'].copy())
    radar_vpt.fields['SNR_HS']['long_name'] = 'radar_signal_to_noise_ratio_copolar_h_from_spectral_processing'

    radar_vpt.add_field_like("ZDR", "ZDR_S",radar_vpt.fields["ZDR"]['data'].copy())
    radar_vpt.fields['ZDR_S']['long_name'] = 'radar_differential_reflectivity_hv_from_spectral_processing'


    radar_vpt.add_field_like("VEL_HV", "VEL_HVS",radar_vpt.fields["VEL_HV"]['data'].copy())
    radar_vpt.fields['VEL_HVS']['long_name'] = 'radial_velocity_of_scatterers_away_from_instrument_from_spectral_processing'

    radar_vpt.add_field_like("SPW_HV", "SPW_HVS",radar_vpt.fields["SPW_HV"]['data'].copy())
    radar_vpt.fields['SPW_HVS']['long_name'] = 'radar_doppler_spectrum_width_from_spectral_processing'

    radar_vpt.add_field_like("VEL_H", "VEL_HS",radar_vpt.fields["VEL_H"]['data'].copy())
    radar_vpt.fields['VEL_HS']['long_name'] = 'radial_velocity_of_scatterers_away_from_instrument_h_from_spectral_processing'

    radar_vpt.add_field_like("VEL_V", "VEL_VS",radar_vpt.fields["VEL_V"]['data'].copy())
    radar_vpt.fields['VEL_VS']['long_name'] = 'radial_velocity_of_scatterers_away_from_instrument_v_from_spectral_processing'

    nsweep = len(radar_vpt.sweep_number['data'])

    numspec = 5;
    
    for s in range(nsweep):
        vptDS = radar_vpt.extract_sweeps([s]);
        tsDS = nc4.Dataset(ts_l1files[s]);

        rng_km,dBZ_all,VEL_HVS,SPW_HVS = get_camra_spectral_moments(tsDS,numspec,-1,2);
        rng_km,dBZ_HS,VEL_HS,SPW_HS = get_camra_spectral_moments(tsDS,numspec,1,2);
        rng_km,dBZ_VS,VEL_VS,SPW_VS = get_camra_spectral_moments(tsDS,numspec,2,2);

        # Get the start and end indices for the specified sweep
        start_idx = radar_vpt.sweep_start_ray_index['data'][s]
        end_idx = radar_vpt.sweep_end_ray_index['data'][s] + 1
    
        field_data = radar_vpt.fields['DBZ_HS']['data']
        field_data[start_idx:end_idx, :] = dBZ_HS
    
        field_data = radar_vpt.fields['ZDR_S']['data']
        field_data[start_idx:end_idx, :] = dBZ_HS-dBZ_VS;

        field_data = radar_vpt.fields['VEL_HVS']['data']
        field_data[start_idx:end_idx, :] = VEL_HVS;

        field_data = radar_vpt.fields['SPW_HVS']['data']
        field_data[start_idx:end_idx, :] = SPW_HVS;

        field_data = radar_vpt.fields['VEL_HS']['data']
        field_data[start_idx:end_idx, :] = VEL_HS;

        field_data = radar_vpt.fields['VEL_VS']['data']
        field_data[start_idx:end_idx, :] = VEL_VS;

    radar_vpt.fields['DBZ_HS']['data'] += 11.5;

    radar_vpt.fields['SNR_HS']['data'][:,:] = radar_vpt.fields['DBZ_HS']['data'][:,:] - 20.0*np.log10(radar_vpt.range['data'][None,:]);
 

    SNR_HS_linear = 10.0**(0.1*radar_vpt.fields['SNR_HS']['data'][:,:]);

    noise = ma.median(SNR_HS_linear[:,-20:],axis=1);
    noise = ma.median(noise)

    SNR_HS_linear /= noise;

    radar_vpt.fields['SNR_HS']['data'][:,:] = 10.0*ma.log10(SNR_HS_linear);

    # Save the modified radar object to a new file
    mom_l1b_file = os.path.basename(vptfile).replace('fix','vpt');
    mom_l1b_file = os.path.join(mom_l1b_path,mom_l1b_file)
    
    pyart.io.write_cfradial(mom_l1b_file, radar_vpt)

    update_string = 'Add moments from spectral processing'

    update_history_attribute(mom_l1b_file,update_string);

def convert_camra_tsl1tomoments_v2(ts_l1path,mom_l1a_path,mom_l1b_path):

    ts_l1path = os.path.join(ts_l1path);
    ts_l1files = glob.glob(f"{ts_l1path}/*.nc");

    sweep = 0

    vpt_ts_file = ts_l1files[sweep];

    print(mom_l1a_path)
    vptfile = glob.glob(f"{mom_l1a_path}/*.nc")[0];

    radar_vpt = pyart.io.read(vptfile);


    radar_vpt.add_field_like("DBZ_H", "DBZ_HS",radar_vpt.fields["DBZ_H"]['data'].copy())
    radar_vpt.fields['DBZ_HS']['long_name'] = 'radar_equivalent_reflectivity_factor_for_copolar_horizontal_receive_signal_from_spectral_processing'
    
    radar_vpt.add_field_like("SNR_H", "SNR_HS",radar_vpt.fields["SNR_H"]['data'].copy())
    radar_vpt.fields['SNR_HS']['long_name'] = 'radar_signal_to_noise_ratio_copolar_h_from_spectral_processing'

    radar_vpt.add_field_like("ZDR", "ZDR_S",radar_vpt.fields["ZDR"]['data'].copy())
    radar_vpt.fields['ZDR_S']['long_name'] = 'radar_differential_reflectivity_hv_from_spectral_processing'


    radar_vpt.add_field_like("VEL_HV", "VEL_HVS",radar_vpt.fields["VEL_HV"]['data'].copy())
    radar_vpt.fields['VEL_HVS']['long_name'] = 'radial_velocity_of_scatterers_away_from_instrument_from_spectral_processing'

    radar_vpt.add_field_like("SPW_HV", "SPW_HVS",radar_vpt.fields["SPW_HV"]['data'].copy())
    radar_vpt.fields['SPW_HVS']['long_name'] = 'radar_doppler_spectrum_width_from_spectral_processing'

    radar_vpt.add_field_like("VEL_H", "VEL_HS",radar_vpt.fields["VEL_H"]['data'].copy())
    radar_vpt.fields['VEL_HS']['long_name'] = 'radial_velocity_of_scatterers_away_from_instrument_h_from_spectral_processing'

    radar_vpt.add_field_like("VEL_V", "VEL_VS",radar_vpt.fields["VEL_V"]['data'].copy())
    radar_vpt.fields['VEL_VS']['long_name'] = 'radial_velocity_of_scatterers_away_from_instrument_v_from_spectral_processing'

    nsweep = len(radar_vpt.sweep_number['data'])

    numspec = 5;
    
    for s in range(nsweep):
        vptDS = radar_vpt.extract_sweeps([s]);
        tsDS = nc4.Dataset(ts_l1files[s]);

        rng_km,dBZ_all,VEL_HVS,SPW_HVS = get_camra_spectral_moments_v2(tsDS,numspec,-1,2);
        rng_km,dBZ_HS,VEL_HS,SPW_HS = get_camra_spectral_moments_v2(tsDS,numspec,1,2);
        rng_km,dBZ_VS,VEL_VS,SPW_VS = get_camra_spectral_moments_v2(tsDS,numspec,2,2);

        # Get the start and end indices for the specified sweep
        start_idx = radar_vpt.sweep_start_ray_index['data'][s]
        end_idx = radar_vpt.sweep_end_ray_index['data'][s] + 1
    
        field_data = radar_vpt.fields['DBZ_HS']['data']
        field_data[start_idx:end_idx, :] = dBZ_HS
    
        field_data = radar_vpt.fields['ZDR_S']['data']
        field_data[start_idx:end_idx, :] = dBZ_HS-dBZ_VS;

        field_data = radar_vpt.fields['VEL_HVS']['data']
        field_data[start_idx:end_idx, :] = VEL_HVS;

        field_data = radar_vpt.fields['SPW_HVS']['data']
        field_data[start_idx:end_idx, :] = SPW_HVS;

        field_data = radar_vpt.fields['VEL_HS']['data']
        field_data[start_idx:end_idx, :] = VEL_HS;

        field_data = radar_vpt.fields['VEL_VS']['data']
        field_data[start_idx:end_idx, :] = VEL_VS;

    radar_vpt.fields['DBZ_HS']['data'] += 11.5;

    radar_vpt.fields['SNR_HS']['data'][:,:] = radar_vpt.fields['DBZ_HS']['data'][:,:] - 20.0*np.log10(radar_vpt.range['data'][None,:]);
 

    SNR_HS_linear = 10.0**(0.1*radar_vpt.fields['SNR_HS']['data'][:,:]);

    noise = ma.median(SNR_HS_linear[:,-20:],axis=1);
    noise = ma.median(noise)

    SNR_HS_linear /= noise;

    radar_vpt.fields['SNR_HS']['data'][:,:] = 10.0*ma.log10(SNR_HS_linear);

    # Save the modified radar object to a new file
    mom_l1b_file = os.path.basename(vptfile).replace('fix','vpt');
    mom_l1b_file = os.path.join(mom_l1b_path,mom_l1b_file)
    
    pyart.io.write_cfradial(mom_l1b_file, radar_vpt)

    update_string = 'Add moments from spectral processing'

    update_history_attribute(mom_l1b_file,update_string);




def convert_camra_ts_l0b2l1(infile,outpath,tracking_tag="AMOF_20230201132601",data_version="1.0.0"):

    """This routine converts Level 0b time series data from the Chilbolton Advanced Meteorological Radar (CAMRa) to Level 1 data.

    Processing is applied to produce I and Q time series for received data at each range gate.
    I and Q values from the limiting amplifiers in the receiver are compared with transmit I and Q for each pulse.
    These are then scaled to the unit circle and multiplied the square root of the sampled linear reflectivity.
    Processing includes separate calibration offsets for even (H-polarized) and odd (V-polarized) pulses.
    Metadata are added specific to the WIVERN-2 project ground-based observations collected in 2020-2021.

    :param infile: Full path of NetCDF Level 0b input file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l0b_v1.0.nc`
    :type infile: str
    :param outfile: Full path of NetCDF Level 1 output file, e.g. `<path-to-file>/ncas-radar-camra-1_cao_20201210-212823_fix-ts_l1_v1.0.nc`
    :type outfile: str
    :param dBZh_offset: Calibration offset to be applied to H polarized reflectivity factor
    :type dBZh_offset: float
    :param ZDR_offset: Calibration offset that would be applied to differential reflectivity (used to calculate the calibration offset to apply to V polarized reflectivity factor)
    :type ZDR_offset: float
    :param range_offset: Additional range offset in metres to be applied
    :type range_offset: float
    :param version: Version of data product in the format `N.m`, where `N` denotes the major verion and `m` a minor revision.
    :type version: str
    """

    dBZh_offset = 0;
    ZDR_offset = 0;
    range_offset=0;
    
    instrument_tagname = "ncas-radar-camra-1"

    print(f'instrument_tagname={instrument_tagname}')

    # ---------------------------------------
    # Read metadata from YAML instrument file
    # ---------------------------------------  
    #with open(yaml_instrument_file, "r") as stream:
    #    try:
    #        instruments = yaml.safe_load(stream)
    #    except yaml.YAMLError as exc:
    #        print(exc)

    #for elem in instruments:
    #    if instrument_tagname in elem:
    #        instrument = elem[instrument_tagname];

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    #with open(yaml_project_file, "r") as stream:
    #    try:
    #        projects = yaml.safe_load(stream)
    #    except yaml.YAMLError as exc:
    #        print(exc)

    #for p in projects:
    #    if tracking_tag in p:
    #        project = p[tracking_tag];

    #radar_name = instrument["instrument_name"].lower();

    #for n in project["ncas_instruments"]:
    #    if radar_name in n:
    #        project_instrument = n[radar_name];


    #location = project_instrument['platform']['location'].lower();
    
    #RadarDataset = read_camra_raw(infile);

    #scan_type = RadarDataset.scan_type

    #file_timestamp = datetime.datetime.strptime(RadarDataset.metadata["time_coverage_start"],'%Y-%m-%dT%H:%M:%SZ');

    #dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S')

    outfile = os.path.basename(infile)
    print(outfile)
    outfile = outfile.replace('l0b','l1')
    outfile = os.path.join(outpath,outfile);
    #os.path.join(outpath,'{}_{}_{}_{}_l0b_v{}.nc'.format(radar_name,location,dtstr,scan_type.replace('_','-',1),data_version));

    print(f'outfile = {outfile}')

    toexclude = ['ZLO', 'ZHI', 'ZCX', 'ITX', 'QTX', 'IRX', 'QRX'];

    with nc4.Dataset(infile) as src, nc4.Dataset(outfile,mode='w',format='NETCDF4') as dst:

        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)

        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))

        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name not in toexclude:
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]


    try: outfile.close()  # just to be safe, make sure dataset is not already open.
    except: pass
    DSout = nc4.Dataset(outfile,mode='r+',format='NETCDF4')
    print(f'reopened {outfile}')

    DSin = nc4.Dataset(infile);
    dt = cftime.num2pydate(DSin['time'][:],DSin['time'].units);

    DSout.product_version = "v{}".format(data_version) ;
    DSout.processing_level = "1" ;

    #DSout.title = "Calibrated time series from 3 GHz CAMRa radar collected for ESA WIVERN-2 campaign at Chilbolton Observatory";

    comment_string = "Correction to account for inverse square power loss with range has not been applied";
    if len(DSin.comment)>0:
        DSout.comment = DSin.comment + "\n " + comment_string;
    else:
        DSout.comment = comment_string;

    varout = DSout.createVariable('dBZ_offsets_applied','f4',('pulse'),zlib=True);
    varout.long_name = 'dBZ calibration offsets applied for even and odd pulses';
    varout.units = 'dB';
    varout.comment = 'dBZ offsets for even pulses (H-polarized) and odd pulses (V-polarized)';

    varout = DSout.createVariable('I','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar in-phase video signal';
    varout.units = '1';
    varout.comment = 'Values are derived from I/Q on unit circle multiplied by square root of linear reflectivity factor';

    varout = DSout.createVariable('Q','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'co-polar quadrature video signal';
    varout.units = '1';
    varout.comment = 'Values are derived from I/Q on unit circle multiplied by square root of linear reflectivity factor';

    varout = DSout.createVariable('ZCX','f4',('time','pulse','range'),zlib=True);
    varout.ancillary_variables = 'qc_flag';
    varout.long_name = 'cross-polar radar equivalent reflectivity factor';
    varout.units = 'dBZ';
    varout.comment = '';

    varout = DSout.createVariable('qc_flag','u1',('time','pulse','range'),fill_value=255);
    varout.is_quality = 'true';
    varout.qualified_variables = 'I Q ZCX';
    varout.long_name = 'Quality control flag';
    varout.flag_values = np.uint8(0),np.uint8(1), np.uint8(2), np.uint8(3), np.uint8(4), np.uint8(255);
    varout.flag_meanings = 'not_used good_data probably_good_data bad_data data_in_blind_range no_qc_performed';
    varout[:] = 2;

    cable_losses = DSout['cable_losses'][:];
    radar_const  = DSout['radar_constant'][:];
    rec_gain     = DSout['receiver_gain'][:];
    freq         = DSout['frequency'][:];
    prf          = DSout['prf'][:];

    dBZv_offset = dBZh_offset-ZDR_offset;

    dBZcal = radar_const-rec_gain+cable_losses;

    DSout['dBZ_offsets_applied'][::2]  = dBZh_offset;
    DSout['dBZ_offsets_applied'][1::2] = dBZv_offset;

    Zh_cal = 10**((dBZcal+dBZh_offset)/10);
    Zv_cal = 10**((dBZcal+dBZv_offset)/10);

    range_m  = DSin['range'][:];
    range_km = (range_m+range_offset)/1000.; # range in km

    ITX = DSin['ITX'][:,:]; #*DSin['ITX'].scale_factor+DSin['ITX'].add_offset;
    QTX = DSin['QTX'][:,:]; #*DSin['QTX'].scale_factor+DSin['QTX'].add_offset;
    IRX = DSin['IRX'][:,:,:]; #*DSin['IRX'].scale_factor+DSin['IRX'].add_offset;
    QRX = DSin['QRX'][:,:,:]; #*DSin['QRX'].scale_factor+DSin['QRX'].add_offset;

    Vtx = ITX - 1j*QTX;
    Vrx = IRX + 1j*QRX;

    V = np.multiply(Vrx[:,:,:], Vtx[:,:,None]);
    V = ma.masked_where(V==0,V);
    V = ma.divide(np.conj(V),np.abs(V));

    V[:,1::2,:] = V[:,1::2,:]*-1.;

    ZLO  = DSin['ZLO'][:,:,:];
    ZHI  = DSin['ZHI'][:,:,:];

    threshold          = DSin['ZLO'].actual_range[1];

    ZED                = ZLO.copy();
    ZED[ZLO>threshold] = ZHI[ZLO>threshold];

    # Convert to linear ZED
    ZED = np.power(10,ZED/10.);

    ZED[:, ::2,:] = ZED[:, ::2,:]*Zh_cal;
    ZED[:,1::2,:] = ZED[:,1::2,:]*Zv_cal;

    ZCX = DSin['ZCX'][:,:,:];
    ZCX[:, ::2,:] = ZCX[:, ::2,:] + dBZcal + dBZh_offset;
    ZCX[:,1::2,:] = ZCX[:,1::2,:] + dBZcal + dBZv_offset;
 #   add_offset = np.min(ZCX[:]);
 #   scale_factor = (np.max(ZCX[:])-np.min(ZCX[:])) / (2**16-1)
 #   packed_data = np.int16(np.rint((ZCX[:,:,:] - add_offset)/scale_factor));
 #   DSout['ZCX'].scale_factor = np.float32(scale_factor);
 #   DSout['ZCX'].add_offset = np.float32(add_offset);
 #   DSout['ZCX'][:] = packed_data;
    DSout['ZCX'][:] = ZCX;

    V = V*np.sqrt(ZED);

    I = V.real;
    Q = V.imag;

#    add_offset = np.min(I[:]);
#    scale_factor = (np.max(I[:])-np.min(I[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((I[:,:,:] - add_offset)/scale_factor));
#    DSout['I'].scale_factor = np.float32(scale_factor);
#    DSout['I'].add_offset = np.float32(add_offset);
    DSout['I'][:] = I;

#    add_offset = np.min(Q[:]);
#    scale_factor = (np.max(Q[:])-np.min(Q[:])) / (2**16-1)
#    packed_data = np.int16(np.rint((Q[:,:,:] - add_offset)/scale_factor));
#    DSout['Q'].scale_factor = np.float32(scale_factor);
#    DSout['Q'].add_offset = np.float32(add_offset);
    DSout['Q'][:] = Q

    blind_range = np.arange(15);
    DSout['qc_flag'][:,:,blind_range] = 4;

    DSout['range'][:] = range_m;
    DSout['range'].range_offset_applied += range_offset;
    DSout['range'].comment = "range_offset_applied includes offset applied by the data acquisition program";

    # -----------------------
    # Update history metadata
    # -----------------------
    user = getpass.getuser()

    updttime = datetime.datetime.now(datetime.timezone.utc)
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " program: wivern_chilbolton_utils.py convert_camra_ts_l0b2l1"
    + " version:" + str(module_version));

    DSout.history = history + "\n" + DSin.history;

    DSout.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    DSin.close();
    DSout.close();

    return


def camra_ts2spectra(tsDS,numspec,pol,win):
    
    if pol==1:
        I = tsDS.variables['I'][:,::2,:];
        Q = tsDS.variables['Q'][:,::2,:];
    elif pol==2:  
        I = tsDS.variables['I'][:,1::2,:];
        Q = tsDS.variables['Q'][:,1::2,:];
    else:
        I = tsDS.variables['I'][:,:,:];
        Q = tsDS.variables['Q'][:,:,:];
    
    rng = tsDS['range'][:];
    start_time = tsDS['time'][0]*100.;

    nray    = I.shape[0]
    npulse  = I.shape[1];
    nsample = I.shape[2];
    
    num = npulse//numspec; 
        
    P = np.empty((nray,num,nsample));
    
    for ispec in np.arange(numspec):
        ind = np.arange(num)+ispec*num;
        [tmpS, norm] = calc_spectra(I[:,ind,:],Q[:,ind,:],win);
        S = np.absolute(tmpS);
        P = P+np.power(S,2)*norm;
        
    P = P/numspec;
    P = np.fft.fftshift(P,axes=1);
    
    return P


def calc_spectra(I,Q,win):
    # Calculates spectra from IQ data
    # Inputs: IQ, win (windowing switch)
    # win = -1 No window; 0 Blackmann; 1 Hanning; 2 Hamming
    
    #from numpy.fft import fft, fftshift

    nray = I.shape[0]
    nfft = I.shape[1];
    nsample = I.shape[2];
    data = I + 1j*Q;
    
    S = np.empty((nray,nfft,nsample)).astype(complex);
    
    if win==-1:
        wind = np.ones((nfft));
    elif win==0:
        wind = np.blackman(nfft);
    elif win==1:
        wind = np.hanning(nfft);
    elif win==2:
        wind = np.hamming(nfft);
        
    norm = 1/np.mean(np.power(wind,2));
 
    dmean = np.mean(data,axis=1);
    S = np.fft.fft(data*wind[None,:,None],axis=1);
    
    return [S, norm]

def get_camra_spectral_moments(tsDS,numspec,pol,win):
    # inputs: raynum (ray number to read)
    #         numspec (subdivide 6100 pulses into numspec spectra)
    #         pol (H=1, V=2, HV=0, H+V=-1)
        
    if pol>-1:
        P = camra_ts2spectra(tsDS,numspec,pol,win);
    else:
        P1 = camra_ts2spectra(tsDS,numspec,1,win);
        P2 = camra_ts2spectra(tsDS,numspec,2,win);
        P = (P1+P2)/2.;

    print(P.shape);
   
    #Range gate spacing (74.94811m = Sampling at 2.0 MHz)
    gate_len     = 0.07494811;
    range_offset = -0.864; #-0.915; #-0.840;       # Range offset (0.0km)
    freq         = 3.0765e9;  # Hz
    PRF = 612;
    if pol!=0:
        PRF=PRF/2.0;
      
    noise_nom = 2e-6; # Nominal noise level
#    noise_nom = 4e-6; # Nominal noise level

            
    Zcal=-60; #dB
    blindrange = 12; #13; # Blind range in 75-m range gates       
    P=P*np.power(10,(Zcal/10));
    
     
    #Remove blind range
    maP = ma.array(P);
    maP[:,:,0:blindrange] = ma.masked;

    [nrays,nbins,ngates]=P.shape;
      
    c = 299792458; # m/s
    lamda    = c/freq;
    vel_fold  = lamda*PRF/4.0;
    bin_width = 2*vel_fold/nbins;
    delay_clocks = 4;

    # Bins
    bins=np.arange(-nbins/2,nbins/2)*bin_width;  

    rng_km = (0.5*299792458/2000000*(np.arange(ngates)+delay_clocks-11.5)-2.1)/1000.;


    # Blank central bins
    cent=round(nbins/2); #+1;

    maP[:,cent-1:cent+1,:]=ma.masked;

    farP = maP[:,:,-20::];
    med_noise = np.ma.median(farP,axis=2);
    farPH = farP[:,::2,:];
    farPV = farP[:,1::2,:];
    
    med_noiseH = np.ma.median(farPH,axis=2);
    med_noiseV = np.ma.median(farPV,axis=2);
    
    noise = med_noise*nbins;

    maP[:,::2,:] -= med_noiseH[:,:,None];
    maP[:,1::2,:]-= med_noiseV[:,:,None];
    
    maP[:,cent-1:cent+1,:]=ma.masked;

    # Remove spikes from far range
    #for ibin in np.arange(slen):
    #    tmp = maP[ibin,120::];
    #    bin_noise = np.ma.median(tmp);
    #    if bin_noise>1e-4:
    #        maP[ibin,:]=maP[ibin,:]-bin_noise;
        
    # Remove clutter (basic approach: triangular blanking zone)
    for ih in np.arange(ngates):
        thresh=(6.0-rng_km[ih])*0.15;
        ind=np.where(np.abs(bins)<thresh);
        maP[:,ind,ih]=ma.masked;

    maP[:,:,0:blindrange] = ma.masked;
    
    Z   = ma.sum(maP,axis=1);

    #vel = ma.array((nrays,ngates));
    #spw = ma.array((nrays,ngates));
    
    #for ih in np.arange(ngates):
    #    vkern = maP[:,:,ih]*bins[None,:];
    #    vel[:,ih] = ma.sum(vkern,axis=1)/Z[:,ih];
    #    spwkern = maP[:,:,ih]*np.power(bins[None,:]-vel[:,ih],2);
    #    spw[:,ih] = np.sqrt(spwkern.sum(axis=1)/Z[:,ih]);

    
    velkern = maP*bins[None,:,None]
    
    vel = ma.sum(velkern,axis=1)/Z;
    
    spwkern = maP*np.power(bins[None,:,None]-vel[:,None,:],2);
    spw = np.sqrt(ma.sum(spwkern,axis=1)/Z)
    
    marng = ma.array(rng_km);
    marng[marng<=0.] = ma.masked;
        
    dBR = 20*np.log10(marng);
    
    dBZ = 10*np.log10(Z)+dBR;

    return rng_km, dBZ, vel, spw

def get_camra_spectral_moments_v2(tsDS,numspec,pol,win):
    # inputs: raynum (ray number to read)
    #         numspec (subdivide 6100 pulses into numspec spectra)
    #         pol (H=1, V=2, HV=0, H+V=-1)
        
    if pol>-1:
        P = camra_ts2spectra(tsDS,numspec,pol,win);
    else:
        P1 = camra_ts2spectra(tsDS,numspec,1,win);
        P2 = camra_ts2spectra(tsDS,numspec,2,win);
        P = (P1+P2)/2.;

    print(P.shape);

   
    #Range gate spacing (74.94811m = Sampling at 2.0 MHz)
    gate_len     = 0.07494811;
    range_offset = -0.864; #-0.915; #-0.840;       # Range offset (0.0km)
    freq         = 3.0765e9;  # Hz
    PRF = 612;
    if pol!=0:
        PRF=PRF/2.0;
      
    noise_nom = 2e-6; # Nominal noise level
#    noise_nom = 4e-6; # Nominal noise level

            
    Zcal=-60; #dB
    blindrange = 12; #13; # Blind range in 75-m range gates       
    P=P*np.power(10,(Zcal/10));
    
    [nrays,nbins,ngates]=P.shape;

    Pmean_noise = np.empty([nrays,ngates],dtype=float)
    Pthresh = np.empty([nrays,ngates],dtype=float)
    Pvar = np.empty([nrays,ngates],dtype=float)

    for iray in range(nrays):
        for igate in range(ngates):
            Pmean, Pthreshold, Pvariance, nnoise = pyart.util.estimate_noise_hs74(P[iray,:,igate], navg=1, nnoise_min=1)   
            Pthresh[iray,igate] = Pthreshold
            Pvar[iray,igate] = Pvariance
            Pmean_noise[iray,igate] = Pmean
        #print(f'Ray = {iray}')
        #print(Pthresh[iray,:])
    

    #Remove blind range
    maP = ma.array(P);
    maP[:,:,0:blindrange] = ma.masked;

      
    c = 299792458; # m/s
    lamda    = c/freq;
    vel_fold  = lamda*PRF/4.0;
    bin_width = 2*vel_fold/nbins;
    delay_clocks = 4;

    # Bins
    bins=np.arange(-nbins/2,nbins/2)*bin_width;  

    rng_km = (0.5*299792458/2000000*(np.arange(ngates)+delay_clocks-11.5)-2.1)/1000.;


    # Blank central bins
    cent=round(nbins/2); #+1;

    maP[:,cent-1:cent+1,:]=ma.masked;

    farP = maP[:,:,-20::];
    med_noise = np.ma.median(farP,axis=2);
    farPH = farP[:,::2,:];
    farPV = farP[:,1::2,:];
    
    med_noiseH = np.ma.median(farPH,axis=2);
    med_noiseV = np.ma.median(farPV,axis=2);
    
    noise = med_noise*nbins;

    maP[:,::2,:] -= med_noiseH[:,:,None];
    maP[:,1::2,:]-= med_noiseV[:,:,None];
    
    maP[:,cent-1:cent+1,:]=ma.masked;

    # Remove spikes from far range
    #for ibin in np.arange(slen):
    #    tmp = maP[ibin,120::];
    #    bin_noise = np.ma.median(tmp);
    #    if bin_noise>1e-4:
    #        maP[ibin,:]=maP[ibin,:]-bin_noise;
        
    # Remove clutter (basic approach: triangular blanking zone)
    for ih in np.arange(ngates):
        thresh=(6.0-rng_km[ih])*0.15;
        ind=np.where(np.abs(bins)<thresh);
        maP[:,ind,ih]=ma.masked;

    maP[:,:,0:blindrange] = ma.masked;
    
    Z   = ma.sum(maP,axis=1);

    #vel = ma.array((nrays,ngates));
    #spw = ma.array((nrays,ngates));
    
    #for ih in np.arange(ngates):
    #    vkern = maP[:,:,ih]*bins[None,:];
    #    vel[:,ih] = ma.sum(vkern,axis=1)/Z[:,ih];
    #    spwkern = maP[:,:,ih]*np.power(bins[None,:]-vel[:,ih],2);
    #    spw[:,ih] = np.sqrt(spwkern.sum(axis=1)/Z[:,ih]);

    
    velkern = maP*bins[None,:,None]
    
    vel = ma.sum(velkern,axis=1)/Z;
    
    spwkern = maP*np.power(bins[None,:,None]-vel[:,None,:],2);
    spw = np.sqrt(ma.sum(spwkern,axis=1)/Z)
    
    marng = ma.array(rng_km);
    marng[marng<=0.] = ma.masked;
        
    dBR = 20*np.log10(marng);
    
    dBZ = 10*np.log10(Z)+dBR;

    return rng_km, dBZ, vel, spw


#def convert_camra_ts_l12mom(infile,outpath):
    #

def near_field_correction_camra_westbrook(rng):
#
#    db_correction = near_field_correction_camra_chriswestbrook(range)
#
#    Calculates near field correction to be applied to the CAMRa radar
#    Based on Sekelsky method (ref below) with empirical adjustment 
#    determined by comparison with 35GHz profiles in low Z ice clouds
#
#    INPUTS: range (optional)   range in m
# 
#    Output: corrections in dB units 
#  
#    Accuracy of correction is not well established at very short ranges (<3km or so)
#
#    Reference: Sekelsky, 2002: Near-Field Reflectivity and Antenna Boresight 
#               Gain Corrections for Millimeter-Wave Atmospheric Radars. 
#               JTECH, 19, 468-477. 
    
    antenna_diameter=25; # metres for CAMRa
    lamda=0.0975; # metres for CAMRa
    
    far_field_distance = 2 * np.power(antenna_diameter,2) / lamda;
    
    a1 = 5.25e-5 + np.power(rng/far_field_distance,2.5) 
    a2 = 0.0117 + np.power(rng/far_field_distance,2.5)
                                                            
    linear_correction = np.real(a2/a1);
    
    linear_correction[linear_correction < 1] = 1;

    db_correction = 10.*np.log10(linear_correction);

    # EMPIRICAL TWEAK TO SEKELSKY CORRECTION - BASED ON 3 vs 35GHz PROFILE COMPARISONS
    #CHRIS_EXTRA_NEAR_FIELD_CORRECTION=(rng-6750)*0.37e-3;
    #db_correction[rng<7250]=db_correction[rng<7250]-CHRIS_EXTRA_NEAR_FIELD_CORRECTION[rng<7250];
    #db_correction[rng>7250]=0;
    
    return db_correction


def get_noise_levels_cfradial(RadarDS):

    # ----------------
    # Open NetCDF file
    # ----------------
    #RadarDS = pyart.io.read_cfradial(ncfile);

    nsweeps = RadarDS.nsweeps;

    noiseh_arr = [];
    noisev_arr = [];
    noisex_arr = [];
    times_arr = [];
    sdnoiseh_arr = [];
    sdnoisev_arr = [];
    sdnoisex_arr = [];
    
    for sweep in range(nsweeps):

        SweepDS = RadarDS.extract_sweeps([sweep]);
    
        x,y,height = SweepDS.get_gate_x_y_z(0)

        dtime = cftime.num2pydate(SweepDS.time['data'],SweepDS.time['units']);

        if len(dtime)>0:
            timestamp = datetime.datetime.strftime(dtime[0],'%Y-%m-%dT%H:%M:%SZ');
        else:
            timestamp = ""

        Zh = SweepDS.fields['DBZ_H']['data'];
        Zv = Zh - SweepDS.fields['ZDR']['data']; # Vertically polarised reflectivity [dBZ] (copolar)
        Zx = Zh + SweepDS.fields['LDR']['data']; 

        linZh = 10.0**(Zh/10.); 
        linZv = 10.0**(Zv/10.); 
        linZx = 10.0**(Zx/10.); 

        rng = SweepDS.range['data'];
        
        rangem = ma.masked_where(rng<=0.0, rng);
        rangekm = rangem/1000.;

        signalpowerh = linZh/rangem[None,:]**2;
        signalpowerv = linZv/rangem[None,:]**2;
        signalpowerx = linZx/rangem[None,:]**2;

        masked_signalpowerh = np.ma.masked_where(height < 14000, signalpowerh);

        noiseh = np.ma.median(masked_signalpowerh);

        unmasked_indices = np.argwhere(~masked_signalpowerh.mask) 

        noisev = np.ma.median(np.ma.masked_where(height < 14000, signalpowerv));
        noisex = np.ma.median(np.ma.masked_where(height < 14000, signalpowerx));

        signalpowerh-=noiseh;	# subtract noise from signal
        signalpowerv-=noisev;
        signalpowerx-=noisex;

        sdnoiseh = np.ma.std(np.ma.masked_where(height < 14000, signalpowerh));
        sdnoisev = np.ma.std(np.ma.masked_where(height < 14000, signalpowerv));
        sdnoisex = np.ma.std(np.ma.masked_where(height < 14000, signalpowerx));

        signalpowerh = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh, signalpowerh); 
        signalpowerv = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev, signalpowerv); 
        signalpowerx = ma.masked_where(signalpowerx<SNR_threshold_X*sdnoisex, signalpowerx);

        SNRperpulseH=signalpowerh/noiseh; 
        SNRperpulseV=signalpowerv/noisev;
        SNRperpulseX=signalpowerx/noisex;


        noiseh_arr.append(10.0*np.log10(noiseh)+60.);
        noisev_arr.append(10.0*np.log10(noisev)+60.);
        noisex_arr.append(10.0*np.log10(noisex)+60.);
        sdnoiseh_arr.append(10.0*np.log10(sdnoiseh));
        sdnoisev_arr.append(10.0*np.log10(sdnoisev));
        sdnoisex_arr.append(10.0*np.log10(sdnoisex));

        times_arr.append(timestamp);

    data_list = []
    for ts, nh, nv, nx, sdh, sdv, sdx in zip(times_arr, noiseh_arr, noisev_arr, noisex_arr, sdnoiseh_arr, sdnoisev_arr, sdnoisex_arr):
        data_list.append({
            'timestamp': ts,'noiseh': nh,'noisev': nv,'noisex': nx, 'sdnoiseh': sdh, 'sdnoisev': sdv, 'sdnoisex': sdx})

    df = pd.DataFrame(data_list)

    valid_entries= df.map(lambda x: np.nan if isinstance(x, np.ma.core.MaskedConstant) else x).dropna()

    def has_positive_noise(row):
        return row['noiseh'] > 0 or row['noisev'] > 0 or row['noisex'] > 0

    # Drop rows with any positive values
    df_filtered = valid_entries[~valid_entries.apply(has_positive_noise, axis=1)]

    def remove_outliers(df):
        Q1 = df[['noiseh', 'noisev', 'noisex', 'sdnoiseh', 'sdnoisev', 'sdnoisex']].quantile(0.25)
        Q3 = df[['noiseh', 'noisev', 'noisex', 'sdnoiseh', 'sdnoisev', 'sdnoisex']].quantile(0.75)
        IQR = Q3 - Q1
        return df[~((df[['noiseh', 'noisev', 'noisex', 'sdnoiseh', 'sdnoisev', 'sdnoisex']] < (Q1 - 1.5 * IQR)) | (df[['noiseh', 'noisev', 'noisex', 'sdnoiseh', 'sdnoisev', 'sdnoisex']] > (Q3 + 1.5 * IQR))).any(axis=1)]

    def remove_sd_outliers(df):
        Q1 = df[['sdnoiseh', 'sdnoisev', 'sdnoisex']].quantile(0.25)
        Q3 = df[['sdnoiseh', 'sdnoisev', 'sdnoisex']].quantile(0.75)
        IQR = Q3 - Q1
        return df[~((df[['sdnoiseh', 'sdnoisev', 'sdnoisex']] < (Q1 - 1.5 * IQR)) | (df[['sdnoiseh', 'sdnoisev', 'sdnoisex']] > (Q3 + 1.5 * IQR))).any(axis=1)]
    
    df_new = remove_sd_outliers(valid_entries);

    #noise_list = df_new.to_dict('records')

    return df_new

def do_noise_subtraction_cfradial(RadarDS,noisefile):

    df_noise = pd.read_csv(noisefile);

    dtimes_noise = [datetime.datetime.strptime(dt, '%Y-%m-%dT%H:%M:%SZ') for dt in df_noise.loc[:,'timestamp']];

    # ----------------
    # Open NetCDF file
    # ----------------
    # RadarDS = pyart.io.read_cfradial(ncfile);
    
    # Create a new field dictionary for SNR
    snrh_field = {
        'data': RadarDS.fields['DBZ_H']['data'].copy(),
        'units': 'dB',
        'long_name': 'radar_signal_to_noise_ratio_copolar_h',
        'proposed_standard_name': 'radar_signal_to_noise_ratio_copolar_h'
    }

    snrv_field = {
        'data': RadarDS.fields['DBZ_H']['data'].copy(),
        'units': 'dB',
        'long_name': 'radar_signal_to_noise_ratio_copolar_v',
        'proposed_standard_name': 'radar_signal_to_noise_ratio_copolar_v'
    }

    snrx_field = {
        'data': RadarDS.fields['DBZ_H']['data'].copy(),
        'units': 'dB',
        'long_name': 'radar_signal_to_noise_ratio_crosspolar_v',
        'proposed_standard_name': 'radar_signal_to_noise_ratio_crosspolar_v'
    }

    # Add the SNR field to the radar object
    RadarDS.add_field('SNR_H', snrh_field, replace_existing=True)
    RadarDS.add_field('SNR_V', snrv_field, replace_existing=True)
    RadarDS.add_field('SNR_X', snrx_field, replace_existing=True)
    
    nsweeps = RadarDS.nsweeps;

    close_idx = [];
    time_diff = [];
    
    for sweep in range(nsweeps):

        SweepDS = RadarDS.extract_sweeps([sweep]);
    
        dtime = cftime.num2pydate(SweepDS.time['data'][0],SweepDS.time['units']);

        time_diffs = [abs((dt-dtime).total_seconds()) for dt in dtimes_noise];

        closest = np.argmin(time_diffs);

        print(f'Scan:{dtime} Noise:{dtimes_noise[closest]}');
        
        close_idx.append(closest);
        time_diff.append(time_diffs[closest]);

        df_use = df_noise.loc[closest];
        noiseh = 10.0**((df_use['noiseh']-60.0)/10.)
        noisev = 10.0**((df_use['noisev']-60.0)/10.)
        noisex = 10.0**((df_use['noisex']-60.0)/10.);
        sdnoiseh = 10.0**(df_use['sdnoiseh']/10.);
        sdnoisev = 10.0**(df_use['sdnoisev']/10.);
        sdnoisex = 10.0**(df_use['sdnoisex']/10.);

        print(noiseh,noisev,noisex);
        print(sdnoiseh,sdnoisev,sdnoisex);
        
        Zh = SweepDS.fields['DBZ_H']['data'].copy();
        Zv = Zh - SweepDS.fields['ZDR']['data'].copy(); # Vertically polarised reflectivity [dBZ] (copolar)
        Zx = Zh + SweepDS.fields['LDR']['data'].copy(); 

        linZh = 10.0**(Zh/10.); 
        linZv = 10.0**(Zv/10.); 
        linZx = 10.0**(Zx/10.); 

        rng = SweepDS.range['data'].copy();
        
        rangem = ma.masked_where(rng<=0.0, rng);
        rangekm = rangem/1000.;

        signalpowerh = linZh/rangem[None,:]**2;
        signalpowerv = linZv/rangem[None,:]**2;
        signalpowerx = linZx/rangem[None,:]**2;

        signalpowerh-=noiseh;	# subtract noise from signal
        signalpowerv-=noisev;
        signalpowerx-=noisex;

        #signalpowerh = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh, signalpowerh); 
        #signalpowerv = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev, signalpowerv); 
        #signalpowerx = ma.masked_where(signalpowerx<SNR_threshold_X*sdnoisex, signalpowerx);

        SNRperpulseH=signalpowerh/noiseh; 
        SNRperpulseV=signalpowerv/noisev;
        SNRperpulseX=signalpowerx/noisex;

        sweep_slice = RadarDS.get_slice(sweep);

        print(sweep_slice);
        
        linZh = signalpowerh * rangem[None,:]**2;
        linZv = signalpowerv * rangem[None,:]**2;
        linZx = signalpowerx * rangem[None,:]**2;

        Zh = 10.*ma.log10(linZh);
        Zv = 10.*ma.log10(linZv);
        Zx = 10.*ma.log10(linZx);

        Zh = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh, Zh); 
        Zv = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev, Zv); 
        Zx = ma.masked_where(signalpowerx<SNR_threshold_X*sdnoisex, Zx);

        RadarDS.fields['DBZ_H']['data'][sweep_slice] = Zh;
        RadarDS.fields['ZDR']['data'][sweep_slice] = Zh-Zv;
        RadarDS.fields['LDR']['data'][sweep_slice] = Zx-Zh;

        RadarDS.fields['SNR_H']['data'][sweep_slice] = 10.0*np.log10(SNRperpulseH.copy());
        RadarDS.fields['SNR_V']['data'][sweep_slice] = 10.0*np.log10(SNRperpulseV.copy());
        RadarDS.fields['SNR_X']['data'][sweep_slice] = 10.0*np.log10(SNRperpulseX.copy());

        PDP    = RadarDS.fields['PHIDP']['data'][sweep_slice]; #-5.8; % Differential Phase Shift (deg)
        SPW_HV = RadarDS.fields['SPW_HV']['data'][sweep_slice]; # spectral width using H & V pulses
        VEL_HV = RadarDS.fields['VEL_HV']['data'][sweep_slice]; # Velocity from both H&V pulses 
        VEL_H  = RadarDS.fields['VEL_H']['data'][sweep_slice]; # Velocity from H pulses
        VEL_V  = RadarDS.fields['VEL_V']['data'][sweep_slice]; # Velocity from H pulses
        CXC    = RadarDS.fields['CXC']['data'][sweep_slice];    # copolar cross correlation rho_hv^2

        VEL_HV = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh,VEL_HV);
        VEL_HV = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev,VEL_HV);
        VEL_H = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh,VEL_H);
        VEL_V = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev,VEL_V);
        #CXC = ma.masked_where(SNRperpulseH<SNR_threshold_CXC, CXC); 
        CXC = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh,CXC);
        CXC = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev,CXC);
        SPW_HV = ma.masked_where(signalpowerh<SNR_threshold_SPW*noiseh,SPW_HV); 

        PDP = ma.masked_where(signalpowerh<SNR_threshold_HH*sdnoiseh,PDP);
        PDP = ma.masked_where(signalpowerv<SNR_threshold_VV*sdnoisev,PDP);

        RadarDS.fields['VEL_HV']['data'][sweep_slice] = VEL_HV;
        RadarDS.fields['VEL_H']['data'][sweep_slice] = VEL_H;
        RadarDS.fields['VEL_V']['data'][sweep_slice] = VEL_V;
        RadarDS.fields['PHIDP']['data'][sweep_slice] = PDP;
        RadarDS.fields['CXC']['data'][sweep_slice] = CXC;
        RadarDS.fields['SPW_HV']['data'][sweep_slice] = SPW_HV;
        RadarDS.fields['PHIDP']['data'][sweep_slice] = PDP;

    #CXC = ma.masked_where(CXC<=0, CXC);
    #RHO_HV=ma.power(CXC,0.5);
    #L=-ma.log10(1-RHO_HV);

    clean(RadarDS)
    doubleclean(RadarDS);

    return 


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
        + " program: camra_utils.recalibrate_netcdf_camra"
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
        + " program: camra_utils.classify_clutter" + " version:"
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
        + " program: camra_utils.clean_distant_gates"
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
    + " program: camra_utils.update_metadata"
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
    + " program: camra_utils.subtract_noise_polar"
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
    + " program: camra_utils.apply_phidp_offset"
    + " {}".format(phidp_offset)
    + " version:" + str(module_version))

    dataset.history += "\n" + history
    print(dataset.history)


    dataset.close()

    return


# ------------------------------------------------------------------------------
# Perform clean to remove isolated pixels
# ------------------------------------------------------------------------------
def clean(RadarDS):

    nray    = len(RadarDS.time['data']);
    ngate   = len(RadarDS.range['data']);

    zh  = RadarDS.fields['DBZ_H']['data'][:].copy();
    zdr = RadarDS.fields['ZDR']['data'][:].copy();

    for iray in np.arange(nray):
        #print("iray = ",iray, "\n");
        for igate in np.arange(ngate-1):
            if (ma.is_masked(zh[iray,igate-1]) and
                ma.is_masked(zh[iray,igate+1])):
                zh[iray,igate] = ma.masked;

    RadarDS.fields['DBZ_H']['data'][:] = zh;

    zh_mask = ma.getmask(zh)

    for var in RadarDS.fields:
        val = RadarDS.fields[var]['data'][:]
        val = ma.masked_where(zh_mask, val);
        RadarDS.fields[var]['data'][:]= val;


def clean_old(ncfile):

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
        + " program: camra_utils.clean"
        + " version:" + str(module_version))

    dataset.history += "\n" + history
    print(dataset.history)

    dataset.close()

def doubleclean(RadarDS):


    nray    = len(RadarDS.time['data']);
    ngate   = len(RadarDS.range['data']);

    zh  = RadarDS.fields['DBZ_H']['data'][:].copy();

    for iray in np.arange(nray):
        for igate in np.arange(ngate-2):
            if (ma.is_masked(zh[iray,igate-1])
                and ma.is_masked(zh[iray,igate+1])
                or ma.is_masked(zh[iray,igate+2])):
                    zh[iray,igate] = ma.masked

    RadarDS.fields['DBZ_H']['data'][:] = zh;

    zh_mask = ma.getmask(zh)

    for var in RadarDS.fields:
        val = RadarDS.fields[var]['data'][:]
        val = ma.masked_where(zh_mask, val);
        RadarDS.fields[var]['data'][:]= val;


def doubleclean_old(ncfile):

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
        + " program: camra_utils.doubleclean"
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
    #print('Opening NetCDF file ' + ncfile)
    dataset = nc4.Dataset(ncfile,'r+',format='NETCDF3_CLASSIC')

    dtime = cftime.num2pydate(dataset['time'][:],dataset['time'].units);

    if len(dtime)>0:
        timestamp = datetime.datetime.strftime(dtime[0],'%Y-%m-%dT%H:%M:%SZ');
    else:
        timestamp = ""

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


    return timestamp, np.round(nezh,2), np.round(nezv,2), np.round(nezx,2)


def get_noise_levels_cfradial_old(ncfile):

    # ----------------
    # Open NetCDF file
    # ----------------
    RadarDS = pyart.io.read_cfradial(ncfile);

    x,y,z = RadarDS.get_gate_x_y_z()

    dtime = cftime.num2pydate(RadarDS.time[:],RadarDS.dataset['time'].units);

    if len(dtime)>0:
        timestamp = datetime.datetime.strftime(dtime[0],'%Y-%m-%dT%H:%M:%SZ');
    else:
        timestamp = ""

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

    return timestamp, np.round(nezh,2), np.round(nezv,2), np.round(nezx,2)


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
        + " program: camra_utils.create_processed_file"
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


def process_camra(datestr,inpath,outpath,yaml_project_file,yaml_instrument_file,tracking_tag):

    pattern = '*{}*raw.nc'.format(datestr);

    print(datestr);
    print(inpath);
    datepath = os.path.join(inpath,datestr);

    print(datepath);

    rawfiles = [];
    rawdirs = [];

    for root,dirs,files in os.walk(datepath):
        rawfiles += [os.path.join(root,f) for f in fnmatch.filter(files, pattern)];
        rawdirs += dirs;
    
    print(rawdirs);

    data_version = "1.0";

    #dBZ_offset = 9.0;
    #range_offset = -865.56+864.0;

    l0bpath = os.path.join(outpath,'L0b',datestr);

    os.makedirs(l0bpath,exist_ok=True);

    for dir in rawdirs:
        print("I am Here!");
        os.makedirs(os.path.join(l0bpath,dir),exist_ok=True);

    for f in rawfiles:
        convert_camra_raw2l0b(f,l0bpath,yaml_project_file,yaml_instrument_file,tracking_tag);

    return

def multi_camra2cfrad(
    rawfiles,
    output_dir,
    scan_name="HSRHI",
    tracking_tag="AMOF_20220922221548",
    campaign="woest",
    data_version="1.0.0",
):
    """
    Aggregates single-sweep CAMRa data to a cfradial1 data.
    output_dir(str): Enter the path for output data,
    scan_name(str): "HSRHI"
    """

    from pathlib import Path
    homepath = Path.home()

    yaml_project_file = os.path.join(homepath,'amof_campaigns','{}_project.yml'.format(campaign))
    yaml_instrument_file = os.path.join(homepath,'amof_instruments','amof_instruments.yml')

    out_dir = output_dir
    files = sorted(rawfiles)
    print(files)
    print("Number of files: ", len(files))
    #print(f"gzip_flag={gzip_flag}");

    print('Start to read raw file')

    RadarDS = read_camra_raw(files[0]);

    print('Done reading raw file')
    # Read time and microsec directly from mmclx file

    nc = nc4.Dataset(files[0])
    dtsec = cftime.num2pydate(nc['time'][:],nc['time'].units)
    nc.close();

    dt_ref = dtsec[0].replace(hour=0,minute=0,second=0,microsecond=0);
    time_reference = dt_ref.strftime('%Y-%m-%dT%H:%M:%SZ');
    time_units = f'seconds since {time_reference}'; 

    tsec = cftime.date2num(dtsec,time_units);
    print(time_units);
    
    print(RadarDS.latitude['data']);
    print("Merging all scans into one Volume")
    for i in range(1, len(files)):

        newRadarDS = read_camra_raw(files[i])
        nc = nc4.Dataset(files[i])
        dtsec_new = cftime.num2pydate(nc['time'][:],nc['time'].units)
        nc.close();

        if 'RHI' in scan_name or 'rhi' in scan_name:
            if np.max(newRadarDS.elevation['data'])-np.min(newRadarDS.elevation['data'])!=0:
                print(f'sweep = {i}');
                RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
                tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));
        #        RadarDS.scan_rate['data'] = np.append(RadarDS.scan_rate['data'],newRadarDS.scan_rate['data']);
        #        RadarDS.antenna_transition['data'] = np.append(RadarDS.antenna_transition['data'],newRadarDS.antenna_transition['data']);
        elif 'PPI' in scan_name or 'ppi' in scan_name or 'VAD' in scan_name or 'vad' in scan_name:
            if np.max(newRadarDS.azimuth['data'])-np.min(newRadarDS.azimuth['data'])!=0:
                print(f'sweep = {i}');
                RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
                tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));

        #        RadarDS.scan_rate['data'] = np.append(RadarDS.scan_rate['data'],newRadarDS.scan_rate['data']);
        #        RadarDS.antenna_transition['data'] = np.append(RadarDS.antenna_transition['data'],newRadarDS.antenna_transition['data']);
        elif 'FIX' in scan_name or 'fix' in scan_name:        
            RadarDS = pyart.util.join_radar(RadarDS,newRadarDS);
            tsec = np.append(tsec,cftime.date2num(dtsec_new,time_units));
 
        
    RadarDS.time['units'] = time_units;
    RadarDS.time['data'][:] = tsec;


    fname = os.path.basename(files[0]).split(".")[0]

    out_file = f"{fname}_{scan_name.lower()}.nc"
    print(out_file);
    out_path = os.path.join(out_dir, out_file)
    print(out_path);

    #noisefile = out_path.replace('.nc', '.noise.csv')
    noisefile = os.path.join(out_dir,'ppi.noise.csv')

    #df_noise = get_noise_levels_cfradial(RadarDS)

    #df_noise.to_csv(noisefile, index=False);

    do_noise_subtraction_cfradial(RadarDS,noisefile);

    pyart.io.write_cfradial(out_path, RadarDS, format="NETCDF4")


    DS = nc4.Dataset(out_path,'r+');
    DS.scan_name = scan_name.lower();

    DS.close();

    amend_unitless(out_path)
    lowercase_long_names(out_path)
    time_long_name(out_path)

    # Update history
    update_string = 'Merge single sweep files into cfradial file'
    update_history_attribute(out_path,update_string)

    print(rawfiles[0]);
    cfradial_add_instrument_parameters(rawfiles[0],out_path,yaml_project_file,yaml_instrument_file,tracking_tag, data_version)


    print(out_path);
    print(f'Going to add NCAS metadata - outpath {out_path}')
    outfile = cfradial_add_ncas_metadata(out_path,yaml_project_file,yaml_instrument_file,tracking_tag,data_version);


    return outfile

def convert_angle(angle,offset):
    print(angle,offset)
    if angle >= 360+round(offset):    
        angle -= 360
    return angle

def find_camra_rhi_files(start_time, end_time,azim_min,azim_max,inpath):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S') 
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    print(start_datetime);
    print(end_datetime);
    hrstr = start_datetime.strftime("%H");
    # Define a list to store the found files
    matching_files = []

    #az_search_offset = -38.0; # July ? onwards
    az_search_offset = -8.0; # Before ? July


    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            # Check if the file matches the criteria
            # Example: check if the file name contains the sweep_type and falls within the time range

            if "rhi" in file and file.endswith('rhi-raw.nc'):
                nc = nc4.Dataset(os.path.join(root, file))

                if len(nc.dimensions['time'])>0:

                    file_time = cftime.num2pydate(nc['time'][0],nc['time'].units)
                    azim = (nc['azimuth'][0]) % 360;

                    if start_datetime <= file_time <= end_datetime:
                        if azim_min <= convert_angle(azim,az_search_offset) < azim_max:
                            print(f'{file}: {file_time} {azim_min} {convert_angle(azim,az_search_offset)} {azim_max}');
                            #print(f'{file_time} {convert_angle(azim)}');
                            matching_files.append(os.path.join(root, file))
                nc.close()         
    return sorted(matching_files)


def find_camra_vpt_files(start_time, end_time,inpath):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S') 
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    print(start_datetime);
    print(end_datetime);
    hrstr = start_datetime.strftime("%H");
    # Define a list to store the found files
    matching_files = []

    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(inpath):  # Replace 'path_to_directory' with the actual directory path
        for file in files:
            # Check if the file matches the criteria
            # Example: check if the file name contains the sweep_type and falls within the time range

            if "fix" in file and file.endswith('fix-raw.nc'):
                nc = nc4.Dataset(os.path.join(root, file))

                if len(nc.dimensions['time'])>0:

                    file_time = cftime.num2pydate(nc['time'][0],nc['time'].units)
                    elev = nc['elevation'][0];

                    if start_datetime <= file_time <= end_datetime:
                        if elev > 88.0:
                            print(f'{file}: {file_time} {elev}');
                            matching_files.append(os.path.join(root, file))
                nc.close()         
    return sorted(matching_files)

def find_camra_vpt_ts_files(start_time, end_time,searchpath,level='l0a'):
    # Convert the input times to datetime objects
    start_datetime = datetime.datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S') 
    end_datetime = datetime.datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
    print(start_datetime);
    print(end_datetime);
    hrstr = start_datetime.strftime("%H");
    # Define a list to store the found files
    matching_files = []

    if level=='l0b':
        searchstring = 'vpt'
    else:
        searchstring = 'fix-ts'

    print(f'searchpath={searchpath}')

    # Iterate through files in a specific directory
    for root, dirs, files in os.walk(searchpath):  # Replace 'path_to_directory' with the actual directory path
        print(files)
        for file in files:
            # Check if the file matches the criteria
            # Example: check if the file name contains the sweep_type and falls within the time range

            if searchstring in file and file.endswith('.nc'):
                nc = nc4.Dataset(os.path.join(root, file))

                if len(nc.dimensions['time'])>0:

                    file_time = cftime.num2pydate(nc['time'][0],nc['time'].units)
                    elev = nc['elevation'][0];

                    if start_datetime <= file_time <= end_datetime:
                        if elev > 88.0:
                            print(f'{file}: {file_time} {elev}');
                            matching_files.append(os.path.join(root, file))
                nc.close()    
    return sorted(matching_files)


def process_camra_woest_sop_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,data_version="1.0.0"):


    # Define the list of directories to search  
    vols = [os.path.basename(elem) for elem in glob.glob(os.path.join(indir,'*'))]

    print(vols)

    # Loop through each directory  
    for vol in vols:  
        rhi_files = [];
        ppi_files = [];

        # Construct the pattern to search for  
        pattern = os.path.join(indir,vol, '*rhi*.nc')  
        ppi_pattern = os.path.join(indir,vol,'*ppi*.nc')
    
        # Use glob to find all matching files  
        rhi_files = glob.glob(pattern)  
        ppi_files = glob.glob(ppi_pattern)

        print(rhi_files);
        print(ppi_files);

        if len(rhi_files)>0:
            RadarDS_SOP = multi_camra2cfrad(rhi_files,outdir,scan_name='RHI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");
            
        if len(ppi_files)>0:
            RadarDS_SOP = multi_camra2cfrad(ppi_files,outdir,scan_name='PPI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");


def process_camra_woest_iop_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,data_version="1.0.0"):

    ioppath = os.path.join(indir,'iop');

    # Define the list of directories to search  
    regions = [os.path.basename(elem) for elem in glob.glob(os.path.join(ioppath,'region_*'))]

    print(regions)

    # Initialize the result dictionary  
    result = {}  

    # Loop through each directory  
    for region in regions:  
        # Construct the pattern to search for  
        pattern = os.path.join(ioppath,region, '*.nc')  
    
        # Use glob to find all matching files  
        found_files = glob.glob(pattern)  
    
        # Check if any files were found  
        if found_files:  
            # Extract the region key from the directory name  
            region_key = region.split('_')[1]  
            # Store the found files in the result dictionary  
            result[region_key] = found_files  

    for key in result:

        iop_files = result[key];
        RadarDS_IOP = multi_camra2cfrad(iop_files,outdir,scan_name='SRHI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");

        for elem in iop_files:
            print(key, os.path.basename(elem))
            
        DS = nc4.Dataset(RadarDS_IOP,'r+');
        DS.comment = f"WOEST IOP tracking storm region {key}";
        DS.close();


def process_camra_ccrest_day_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,data_version="1.0.0"):

    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d');
    end_date = start_date + datetime.timedelta(days=1); # - datetime.timedelta(minutes=30);


    rhi_files_270 = find_camra_rhi_files(start_date.strftime('%Y-%m-%d %H:%M:%S'), end_date.strftime('%Y-%m-%d %H:%M:%S'), 260, 280, indir)
    rhi_files_246 = find_camra_rhi_files(start_date.strftime('%Y-%m-%d %H:%M:%S'), end_date.strftime('%Y-%m-%d %H:%M:%S'), 236, 256, indir)

    
    print(rhi_files_270);
    print(rhi_files_246);
    print(len(rhi_files_270));
    print(len(rhi_files_246));

    if (len(rhi_files_270)>0):
        RadarDS_RHI = multi_camra2cfrad(rhi_files_270,outdir,scan_name='RHI-CCREST1',data_version=data_version,tracking_tag="AMOF_20230201132601",campaign="ccrest-m");

    if (len(rhi_files_246)>0):
        RadarDS_RHI = multi_camra2cfrad(rhi_files_246,outdir,scan_name='RHI-CCREST2',data_version=data_version,tracking_tag="AMOF_20230201132601",campaign="ccrest-m");
    
    # Vertically pointing files for whole day
    #try:
    #    vpt_files = find_camra_vpt_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),'vert', indir,gzip_flag=True);
    #    print(vpt_files)
    #    #vpt_files_unzipped = find_mmclxfiles(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),'vert', indir,gzip_flag=False);
    #    if (len(vpt_files)>0):
    #        RadarDS_VPT = multi_mmclx2cfrad(vpt_files,outdir,scan_name='VPT',gzip_flag=True,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);
    #    #elif (len(vpt_files_unzipped)>0):
    #    #    RadarDS_VPT = multi_mmclx2cfrad(vpt_files_unzipped,outdir,scan_name='VPT',gzip_flag=False,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=55.7);
    #except:
    #    print("VPT problem")
    #    pass
    # VAD files for whole day
    #vad_dt_start = datetime.datetime.strptime(datestr,"%Y%m%d");
    #vad_dt_end = vad_dt_start+datetime.timedelta(days=1);
    #try:
    #    vad_files = find_mmclx_vad_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),80,90, indir,gzip_flag=True);
    #    if (len(vad_files)>0):
    #        RadarDS_VAD = multi_mmclx2cfrad(vad_files,outdir,scan_name='VAD',gzip_flag=True,azimuth_offset=azimuth_offset,tracking_tag='AMOF_20230201132601',campaign='ccrest-m',revised_northangle=revised_northangle,data_version=data_version);
    #except:
    #    pass


def process_camra_ccrest_vpt_day_step1(datestr,indir,outdir,yaml_project_file,yaml_instrument_file,data_version="1.0.0"):

    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d');
    end_date = start_date + datetime.timedelta(days=1); # - datetime.timedelta(minutes=30);

 # Vertically pointing files for whole day
    try:
        vpt_files = find_camra_vpt_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),indir);
        print(vpt_files)
        if (len(vpt_files)>0):
            RadarDS_VPT = multi_camra2cfrad(vpt_files,outdir,scan_name='fix',data_version=data_version,tracking_tag="AMOF_20230201132601",campaign="ccrest-m");
    except:
        print("VPT problem")
        pass


def process_camra_ccrest_vpt_day_ts(datestr,ts_indir,mom_indir,outdir,yaml_project_file,yaml_instrument_file,data_version="1.0.0"):

    # Define the start and end times for the loop
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d');
    end_date = start_date + datetime.timedelta(days=1); # - datetime.timedelta(minutes=30);

    l0bpath = os.path.join(outdir,'ts','L0b',datestr);
    l1path = os.path.join(outdir,'ts','L1',datestr)
    mom_l1a_path = os.path.join(mom_indir,datestr)
    mom_l1b_path = os.path.join(outdir,'L1b',datestr)


 # Vertically pointing files for whole day
    try:
        vpt_ts_files_l0a = find_camra_vpt_ts_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),ts_indir,level='l0a');
        print(vpt_ts_files_l0a)
        if (len(vpt_ts_files_l0a)>0):
            for f in vpt_ts_files_l0a:
                convert_camra_ts_l0a2l0b(f,l0bpath,tracking_tag="AMOF_20230201132601",data_version="1.0.0");
    except:
        print("VPT problem l0a2l0b")
        pass
     # Vertically pointing files for whole day
    try:
        vpt_ts_files_l0b = find_camra_vpt_ts_files(start_date.strftime('%Y-%m-%d %H:%M:%S'),end_date.strftime('%Y-%m-%d %H:%M:%S'),l0bpath,level='l0b');
        print(vpt_ts_files_l0b)
        if (len(vpt_ts_files_l0b)>0):
            for f in vpt_ts_files_l0b:
                print(f);
                convert_camra_ts_l0b2l1(f,l1path,tracking_tag="AMOF_20230201132601",data_version="1.0.0");

    except:
        print("VPT problem l0b2l1")
        pass


    convert_camra_tsl1tomoments(l1path,mom_l1a_path,mom_l1b_path);



def update_history_attribute(ncfile,update):

    dataset = nc4.Dataset(ncfile,'r+');

    user = getpass.getuser();

    updttime = datetime.datetime.utcnow();
    updttimestr = updttime.ctime();

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " " + update);

    dataset.history = history + "\n" + dataset.history;

    dataset.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%SZ')

    dataset.close();

def cfradial_add_instrument_parameters(rawfile,cfradfile,yaml_project_file,yaml_instrument_file,tracking_tag,data_version):
    # -------------------------------------------------------
    # Read cfradial file to add instrument parameters
    # -------------------------------------------------------
    DS = nc4.Dataset(cfradfile,'r+');

    print(cfradfile)

    DSin = nc4.Dataset(rawfile,'r');

    print('Creating frequency dimension')
    frequency = DS.createDimension("frequency", 1);

    #varin = DSin['lambda'];
    #lightspeed = 299792458
    #tx_freq = lightspeed/varin[:];
    #print(tx_freq)
    #varout = DS.createVariable('frequency',varin.datatype,("frequency"));
    #varout.standard_name = 'radiation_frequency';
    #varout.long_name = 'frequency_of_transmitted_radiation';
    #varout.units = 's-1';
    #varout[:]=tx_freq;
    #print('Creating meta_group attribute')
    #varout.meta_group = "instrument_parameters";

    #varout = DS['radar_measured_transmit_power_h'];
    #varout.long_name = "radar_measured_transmit_power_h";
    #varout.units = 'dBm'
    #varout.meta_group = "radar_parameters";


    # ----------------
    # Scalar variables
    # ----------------

    #varin = DSin['prf'];
    #varout = DSout.createVariable('prf',varin.datatype);
    #varout.long_name = 'pulse repetition frequency';
    #varout.units = 'Hz';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthH'];
    #varout = DSout.createVariable('beamwidthH',varin.datatype);
    #varout.long_name = 'horizontal angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['beamwidthV'];
    #varout = DSout.createVariable('beamwidthV',varin.datatype);
    #varout.long_name = 'vertical angular beamwidth';
    #varout.units = 'degree';
    #varout[:]=varin[:];

    #varin = DSin['antenna_diameter'];
    #varout = DSout.createVariable('antenna_diameter',varin.datatype);
    #varout.long_name = 'antenna diameter';
    #varout.units = 'm';
    #varout[:]=varin[:];

    #varin = DSin['pulse_period'];
    #varout = DSout.createVariable('pulse_width',varin.datatype);
    #varout.long_name = 'pulse width';
    #varout.units = 'us';
    #varout[:]=varin[:];

    #varin = DSin['transmit_power'];
    #varout = DSout.createVariable('transmit_power',varin.datatype);
    #varout.long_name = 'peak transmitted power';
    #varout.units = 'W';
    #varout[:]=varin[:];


    DSin.close();
    DS.close();

    return

if __name__ == "__main__":
    import sys, PythonCall
    PythonCall.PythonCall(sys.argv).execute()
