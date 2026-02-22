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
SNR_threshold_SPW=10.0; 	# spectral width thresholding


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

    #duration_offset = -0.05; # used for woest
    #duration_offset_ppi = -0.12; # used for woest
    duration_offset = -0.0; # used for kasbex
    duration_offset_ppi = -0.0; # used for kasbex
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

def multi_camra2cfrad(
    rawfiles,
    output_dir,
    scan_name="RHI",
    tracking_tag="AMOF_20220922221548",
    campaign="woest",
    data_version="1.0.1",
    yaml_project_file=None,
    yaml_instrument_file=None,
    noise_file=None
):
    """
    Aggregates single-sweep CAMRa data to a cfradial1 data.
    
    Args:
        rawfiles: List of raw CAMRa files
        output_dir(str): Enter the path for output data,
        scan_name(str): "RHI", "PPI", "MAN", "FIX"
        tracking_tag(str): AMOF tracking tag
        campaign(str): Campaign name
        data_version(str): Data version string
        yaml_project_file(str): Path to project YAML file (optional)
        yaml_instrument_file(str): Path to instrument YAML file (optional)
        noise_file: Path to noise file (optional)
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
    #noisefile = os.path.join(out_dir,'ppi.noise.csv')
    noisefile = os.path.join(out_dir,'camra_noise_20250731.csv')

    #df_noise = get_noise_levels_cfradial(RadarDS)

    #df_noise.to_csv(noisefile, index=False);

    # Use the provided noise file with proper path handling
    if noise_file:
        # Check if the provided noise file exists as specified
        if os.path.exists(noise_file):
            noisefile = noise_file
            print(f"Using specified noise file: {noisefile}")
        else:
            # If the specified file doesn't exist, check if it's meant to be in the output directory
            noise_basename = os.path.basename(noise_file)
            noisefile_in_outdir = os.path.join(output_dir, noise_basename)
            
            if os.path.exists(noisefile_in_outdir):
                noisefile = noisefile_in_outdir
                print(f"Using noise file in output directory: {noisefile}")
            else:
                print(f"Specified noise file not found at: {noise_file}")
                print(f"Also checked in output directory: {noisefile_in_outdir}")
                noisefile = None
    else:
        # Fall back to default noise file location in output directory
        noisefile = os.path.join(output_dir, 'camra_noise_20250731.csv')
        if not os.path.exists(noisefile):
            print(f"Default noise file not found: {noisefile}")
            noisefile = None
    
    # Apply noise processing if file exists
    if noisefile and os.path.exists(noisefile):
        try:
            print(f"Found noise file: {noisefile}")
            df_noise = get_noise_levels_cfradial(RadarDS)
            if not df_noise.empty:
                df_noise.to_csv(noisefile.replace('.csv', '_calculated.csv'), index=False)
            do_noise_subtraction_cfradial(RadarDS, noisefile)
            print(f"Applied noise subtraction using: {noisefile}")
        except Exception as e:
            print(f"Error in noise processing: {e}")
    else:
        print(f"Noise file not found. Skipping noise subtraction.")
    
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

def amend_unitless(cfradfile):
    """
    Amend unitless fields in NetCDF file by removing 'unitless' and 'count' units.
    
    Args:
        cfradfile: Path to CF-Radial NetCDF file
    """
    DS = nc4.Dataset(cfradfile, 'r+')

    # Loop through each variable in the NetCDF file
    for var_name in DS.variables:
        var = DS.variables[var_name]
        if hasattr(var, 'units') and var.units == 'unitless':
            var.units = ""
        if hasattr(var, 'units') and var.units == 'count':
            var.units = ""

    DS.close()

def lowercase_long_names(cfradfile):
    """
    Convert long_name attributes to lowercase.
    
    Args:
        cfradfile: Path to CF-Radial NetCDF file
    """
    DS = nc4.Dataset(cfradfile, 'r+')

    # Loop through each variable in the NetCDF file
    for var_name in DS.variables:
        var = DS.variables[var_name]
        if hasattr(var, 'long_name'):
            var.long_name = var.long_name.lower()
            if "utc" in var.long_name:
                var.long_name = var.long_name.replace("utc", "UTC")

    DS.close()

def time_long_name(cfradfile):
    """
    Update time variable long_name attribute.
    
    Args:
        cfradfile: Path to CF-Radial NetCDF file
    """
    DS = nc4.Dataset(cfradfile, 'r+')
    time_var = DS.variables['time']
    if 'time_reference' in DS.variables:
        time_var.long_name = "time_since_time_reference"

    DS.close()

def update_history_attribute(ncfile, update):
    """
    Update history attribute in NetCDF file.
    
    Args:
        ncfile: Path to NetCDF file
        update: String describing the update
    """
    dataset = nc4.Dataset(ncfile, 'r+')

    user = getpass.getuser()

    updttime = datetime.datetime.utcnow()
    updttimestr = updttime.ctime()

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " " + update)

    dataset.history = history + "\n" + dataset.history

    dataset.last_revised_date = datetime.datetime.strftime(updttime, '%Y-%m-%dT%H:%M:%S')

    dataset.close()

def cfradial_add_instrument_parameters(rawfile, cfradfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version):
    """
    Add instrument parameters to CF-Radial file.
    
    Args:
        rawfile: Path to raw input file
        cfradfile: Path to CF-Radial file
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        tracking_tag: AMOF tracking tag
        data_version: Data version string
    """
    # Read cfradial file to add instrument parameters
    DS = nc4.Dataset(cfradfile, 'r+')

    print(cfradfile)

    DSin = nc4.Dataset(rawfile, 'r')

    print('Creating frequency dimension')
    frequency = DS.createDimension("frequency", 1)

    # Add other instrument parameters as needed
    # This is a placeholder - implement based on your requirements

    DSin.close()
    DS.close()

    return

def single_camra2cfrad(
    rawfile,
    output_dir,
    scan_name="RHI",
    tracking_tag="AMOF_20220922221548",
    campaign="woest",
    data_version="1.0.1",
    yaml_project_file=None,
    yaml_instrument_file=None,
    noise_file=None
):
    """Process a single CAMRa file to CF-Radial format."""
    return multi_camraraw2cfrad(
        [rawfile],
        output_dir,
        scan_name=scan_name,
        tracking_tag=tracking_tag,
        campaign=campaign,
        data_version=data_version,
        yaml_project_file=yaml_project_file,
        yaml_instrument_file=yaml_instrument_file,
        noise_file=noise_file
    )

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


def process_camra_woest_sop_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.1"):


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


def process_camra_woest_iop_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.1"):

    ioppath = os.path.join(indir,'iop');

    # Define the list of directories to search  
    regions = [os.path.basename(elem) for elem in glob.glob(os.path.join(ioppath,'region_*'))]

    print(regions)

    SINGLE_SWEEP_FILES = True

    # Initialize the result dictionary  
    result = {}  

    # Loop through each directory  
    for region in regions:  
        # Construct the pattern to search for  
        pattern = os.path.join(ioppath,region, '*rhi*.nc')  
    
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
    
        if SINGLE_SWEEP_FILES:

            for elem in iop_files:
                RadarDS_IOP = multi_camra2cfrad([elem],outdir,scan_name='RHI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");
                print(key, os.path.basename(elem))
            
                DS = nc4.Dataset(RadarDS_IOP,'r+');
                DS.comment = f"WOEST IOP tracking storm region {key}";
                DS.close();
        
        else:

            RadarDS_IOP = multi_camra2cfrad(iop_files,outdir,scan_name='SRHI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");

            for elem in iop_files:
                print(key, os.path.basename(elem))
            
            DS = nc4.Dataset(RadarDS_IOP,'r+');
            DS.comment = f"WOEST IOP tracking storm region {key}";
            DS.close();

def process_camra_woest_other_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.1"):

    other_path = os.path.join(indir,'other_sop');

    SINGLE_SWEEP_FILES = True

    rhi_files = [];
    ppi_files = [];

    # Construct the pattern to search for  
    pattern = os.path.join(other_path,'*rhi*.nc')  
    ppi_pattern = os.path.join(other_path,'*ppi*.nc')
    
    # Use glob to find all matching files  
    rhi_files = glob.glob(pattern)  
    ppi_files = glob.glob(ppi_pattern)

    print(rhi_files);
    print(ppi_files);
    
    if SINGLE_SWEEP_FILES:

        for elem in rhi_files:
            RadarDS = multi_camra2cfrad([elem],outdir,scan_name='RHI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");
            print(os.path.basename(elem))

        for elem in ppi_files:
            RadarDS = multi_camra2cfrad([elem],outdir,scan_name='PPI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");
            print(os.path.basename(elem))

    else:

        RadarDS = multi_camra2cfrad(rhi_files,outdir,scan_name='SRHI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");
        RadarDS = multi_camra2cfrad(ppi_files,outdir,scan_name='PPI',data_version=data_version,tracking_tag="AMOF_20220922221548",campaign="woest");


def process_camra_ccrest_day_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.0"):

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


def process_camra_ccrest_vpt_day_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.0"):

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


def process_camra_ccrest_vpt_day_ts(datestr, ts_indir, mom_indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.0"):

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

def process_camra_dymecs_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.0"):

    SINGLE_SWEEP_FILES = True

    rhi_files = [];
    ppi_files = [];

    # Construct the pattern to search for  
    pattern = os.path.join(indir,'*rhi*.nc')  
    ppi_pattern = os.path.join(indir,'*ppi*.nc')
    
    # Use glob to find all matching files  
    rhi_files = glob.glob(pattern)  
    ppi_files = glob.glob(ppi_pattern)

    print(rhi_files);
    print(ppi_files);
    
    if SINGLE_SWEEP_FILES:

        for elem in rhi_files:
            RadarDS = multi_camra2cfrad([elem],outdir,scan_name='RHI',data_version=data_version,tracking_tag="CRF_85",campaign="dymecs");
            print(os.path.basename(elem))

        for elem in ppi_files:
            RadarDS = multi_camra2cfrad([elem],outdir,scan_name='PPI',data_version=data_version,tracking_tag="CRF_85",campaign="dymecs");
            print(os.path.basename(elem))

    else:

        RadarDS = multi_camra2cfrad(rhi_files,outdir,scan_name='SRHI',data_version=data_version,tracking_tag="CRF_85",campaign="dymecs");
        RadarDS = multi_camra2cfrad(ppi_files,outdir,scan_name='PPI',data_version=data_version,tracking_tag="CRF_85",campaign="dymecs");


def process_camra_step1(datestr, indir, outdir, yaml_project_file, yaml_instrument_file, data_version="1.0.0", tracking_tag="CRF_85", campaign="dymecs"):

    SINGLE_SWEEP_FILES = True

    rhi_files = [];
    ppi_files = [];

    # Construct the pattern to search for  
    pattern = os.path.join(indir,'*rhi*.nc')  
    ppi_pattern = os.path.join(indir,'*ppi*.nc')
    
    # Use glob to find all matching files  
    rhi_files = glob.glob(pattern)  
    ppi_files = glob.glob(ppi_pattern)

    print(rhi_files);
    print(ppi_files);
    
    if SINGLE_SWEEP_FILES:

        for elem in rhi_files:
            RadarDS = multi_camra2cfrad([elem],outdir,scan_name='RHI',data_version=data_version,tracking_tag=tracking_tag,campaign=campaign);
            print(os.path.basename(elem))

        for elem in ppi_files:
            RadarDS = multi_camra2cfrad([elem],outdir,scan_name='PPI',data_version=data_version,tracking_tag=tracking_tag,campaign=campaign);
            print(os.path.basename(elem))

    else:

        RadarDS = multi_camra2cfrad(rhi_files,outdir,scan_name='SRHI',data_version=data_version,tracking_tag=tracking_tag,campaign=campaign);
        RadarDS = multi_camra2cfrad(ppi_files,outdir,scan_name='PPI',data_version=data_version,tracking_tag=tracking_tag,campaign=campaign);



def update_history_attribute(ncfile,update):

    dataset = nc4.Dataset(ncfile,'r+');

    user = getpass.getuser();

    updttime = datetime.datetime.utcnow();
    updttimestr = updttime.ctime();

    history = updttimestr + (" - user:" + user
    + " machine: " + socket.gethostname()
    + " " + update);

    dataset.history = history + "\n" + dataset.history;

    dataset.last_revised_date = datetime.datetime.strftime(updttime,'%Y-%m-%dT%H:%M:%S')

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



def do_noise_subtraction_cfradial(radar, noisefile):
    """
    Perform noise subtraction on CF-Radial radar data and add SNR and QC fields.
    Uses time-matching to find closest noise measurements for each sweep.
    
    Args:
        radar: PyART radar object
        noisefile: Path to noise CSV file
    """
    import pandas as pd
    from datetime import datetime
    import pytz
    
    try:
        if not os.path.exists(noisefile):
            print(f"Noise file not found: {noisefile}. Skipping noise subtraction.")
            return
            
        print(f"Applying noise subtraction using: {noisefile}")
        
        # Read noise data from CSV file
        df_noise = pd.read_csv(noisefile)
        
        if df_noise.empty:
            print("Noise file is empty. Skipping noise subtraction.")
            return
        
        print(f"Noise file columns: {list(df_noise.columns)}")
        print(f"First few rows:\n{df_noise.head()}")
        
        # Check which columns exist (handle both naming conventions)
        if 'noise_hh' in df_noise.columns:
            hh_col, vv_col, xx_col = 'noise_hh', 'noise_vv', 'noise_xx'
            sd_hh_col, sd_vv_col, sd_xx_col = 'sd_noise_hh', 'sd_noise_vv', 'sd_noise_xx'
        elif 'noiseh' in df_noise.columns:
            hh_col, vv_col, xx_col = 'noiseh', 'noisev', 'noisex'
            sd_hh_col, sd_vv_col, sd_xx_col = 'sdnoiseh', 'sdnoisev', 'sdnoisex'
        else:
            print(f"Error: Could not find noise columns in CSV file. Available columns: {list(df_noise.columns)}")
            return
        
        # Filter out non-numeric values and convert to numeric
        noise_cols = [hh_col, vv_col, xx_col, sd_hh_col, sd_vv_col, sd_xx_col]
        for col in noise_cols:
            if col in df_noise.columns:
                # Replace '--' and other non-numeric values with NaN
                df_noise[col] = pd.to_numeric(df_noise[col], errors='coerce')
        
        # Remove rows with NaN values
        df_noise_clean = df_noise.dropna(subset=noise_cols)
        
        if df_noise_clean.empty:
            print("No valid numeric noise data found after cleaning. Skipping noise subtraction.")
            return
        
        print(f"Using {len(df_noise_clean)} valid noise measurements out of {len(df_noise)} total")
        
        # Convert timestamp column to datetime for time matching
        # Make sure all datetimes are timezone-naive (UTC)
        df_noise_clean['datetime'] = pd.to_datetime(df_noise_clean['timestamp'], utc=True).dt.tz_localize(None)
        
        # Create SNR field dictionaries
        snrh_field = {
            'data': radar.fields['DBZ_H']['data'].copy(),
            'units': 'dB',
            'long_name': 'radar_signal_to_noise_ratio_copolar_h',
            'proposed_standard_name': 'radar_signal_to_noise_ratio_copolar_h'
        }

        snrv_field = {
            'data': radar.fields['DBZ_H']['data'].copy(),
            'units': 'dB',
            'long_name': 'radar_signal_to_noise_ratio_copolar_v',
            'proposed_standard_name': 'radar_signal_to_noise_ratio_copolar_v'
        }

        snrx_field = {
            'data': radar.fields['DBZ_H']['data'].copy(),
            'units': 'dB',
            'long_name': 'radar_signal_to_noise_ratio_crosspolar_v',
            'proposed_standard_name': 'radar_signal_to_noise_ratio_crosspolar_v'
        }

        # Create QC field dictionaries using integer flags
        qc_dbz_h_field = {
            'data': np.zeros_like(radar.fields['DBZ_H']['data'], dtype=np.int8),
            'units': '',
            'long_name': 'quality_control_flag_for_radar_equivalent_reflectivity_factor_h',
            'flag_values': np.array([0, 1, 2, 4, 8], dtype=np.int8),
            'flag_meanings': 'good low_snr_h low_snr_v low_correlation high_spectrum_width',
            'comment': 'QC flags are additive: 0=good, 1=low SNR H-pol, 2=low SNR V-pol, 4=low correlation, 8=high spectrum width'
        }

        qc_zdr_field = {
            'data': np.zeros_like(radar.fields['ZDR']['data'], dtype=np.int8),
            'units': '',
            'long_name': 'quality_control_flag_for_radar_differential_reflectivity_hv',
            'flag_values': np.array([0, 1, 2, 4], dtype=np.int8),
            'flag_meanings': 'good low_snr_h low_snr_v low_correlation',
            'comment': 'QC flags are additive: 0=good, 1=low SNR H-pol, 2=low SNR V-pol, 4=low correlation'
        }

        qc_ldr_field = {
            'data': np.zeros_like(radar.fields['LDR']['data'], dtype=np.int8),
            'units': '',
            'long_name': 'quality_control_flag_for_radar_linear_depolarization_ratio',
            'flag_values': np.array([0, 1, 16], dtype=np.int8),
            'flag_meanings': 'good low_snr_h low_snr_x',
            'comment': 'QC flags are additive: 0=good, 1=low SNR H-pol, 16=low SNR X-pol'
        }

        # Initialize QC fields for velocity, correlation, and other parameters
        qc_fields = {}
        if 'VEL_HV' in radar.fields:
            qc_fields['VEL_HV'] = {
                'data': np.zeros_like(radar.fields['VEL_HV']['data'], dtype=np.int8),
                'units': '',
                'long_name': 'quality_control_flag_for_radial_velocity_of_scatterers_away_from_instrument',
                'flag_values': np.array([0, 1, 2, 8], dtype=np.int8),
                'flag_meanings': 'good low_snr_h low_snr_v high_spectrum_width',
                'comment': 'QC flags are additive'
            }

        if 'CXC' in radar.fields:
            qc_fields['CXC'] = {
                'data': np.zeros_like(radar.fields['CXC']['data'], dtype=np.int8),
                'units': '',
                'long_name': 'quality_control_flag_for_radar_copolar_cross_correlation',
                'flag_values': np.array([0, 1, 2], dtype=np.int8),
                'flag_meanings': 'good low_snr_h low_snr_v',
                'comment': 'QC flags are additive'
            }

        if 'PHIDP' in radar.fields:
            qc_fields['PHIDP'] = {
                'data': np.zeros_like(radar.fields['PHIDP']['data'], dtype=np.int8),
                'units': '',
                'long_name': 'quality_control_flag_for_radar_differential_phase_hv',
                'flag_values': np.array([0, 1, 2, 4], dtype=np.int8),
                'flag_meanings': 'good low_snr_h low_snr_v low_correlation',
                'comment': 'QC flags are additive'
            }

        if 'SPW_HV' in radar.fields:
            qc_fields['SPW_HV'] = {
                'data': np.zeros_like(radar.fields['SPW_HV']['data'], dtype=np.int8),
                'units': '',
                'long_name': 'quality_control_flag_for_radar_doppler_spectrum_width',
                'flag_values': np.array([0, 1], dtype=np.int8),
                'flag_meanings': 'good low_snr_h',
                'comment': 'QC flags: 0=good, 1=low SNR H-pol'
            }

        # Add the fields to the radar object
        radar.add_field('SNR_H', snrh_field, replace_existing=True)
        radar.add_field('SNR_V', snrv_field, replace_existing=True)
        radar.add_field('SNR_X', snrx_field, replace_existing=True)
        radar.add_field('QC_DBZ_H', qc_dbz_h_field, replace_existing=True)
        radar.add_field('QC_ZDR', qc_zdr_field, replace_existing=True)
        radar.add_field('QC_LDR', qc_ldr_field, replace_existing=True)
        
        for field_name, qc_field in qc_fields.items():
            radar.add_field(f'QC_{field_name}', qc_field, replace_existing=True)
        
        nsweeps = radar.nsweeps
        
        # Process each sweep individually with time matching (following camra_utils_old approach)
        for sweep in range(nsweeps):
            print(f"Processing sweep {sweep + 1}/{nsweeps}")
            
            # Extract this sweep
            SweepDS = radar.extract_sweeps([sweep])
            
            # Get time for this sweep
            dtime = cftime.num2pydate(SweepDS.time['data'], SweepDS.time['units'])
            if len(dtime) > 0:
                sweep_time = dtime[0]
                print(f"  Sweep time: {sweep_time}")
                
                # Ensure sweep_time is timezone-naive
                if hasattr(sweep_time, 'tzinfo') and sweep_time.tzinfo is not None:
                    if sweep_time.tzinfo != pytz.UTC:
                        sweep_time = sweep_time.astimezone(pytz.UTC)
                    sweep_time = sweep_time.replace(tzinfo=None)
                
                # Find closest noise measurement in time
                time_diffs = abs(df_noise_clean['datetime'] - sweep_time)
                closest_idx = time_diffs.idxmin()
                
                closest_noise = df_noise_clean.loc[closest_idx]
                time_diff_minutes = time_diffs.loc[closest_idx].total_seconds() / 60.0
                
                print(f"  Using noise from {closest_noise['timestamp']} (diff: {time_diff_minutes:.1f} min)")
                
                # Get noise levels for this sweep - convert from dBm+60 to linear
                noiseh = 10**((closest_noise[hh_col] - 60.0) / 10.)
                noisev = 10**((closest_noise[vv_col] - 60.0) / 10.)
                noisex = 10**((closest_noise[xx_col] - 60.0) / 10.)
                
                # Get standard deviations - convert from dB to linear
                sdnoiseh = 10**(closest_noise[sd_hh_col] / 10.)
                sdnoisev = 10**(closest_noise[sd_vv_col] / 10.)
                sdnoisex = 10**(closest_noise[sd_xx_col] / 10.)
                
                print(f"  Noise levels: H={noiseh:.2e}, V={noisev:.2e}, X={noisex:.2e}")
                
            else:
                print(f"  Warning: No time data for sweep {sweep}, using median noise")
                # Fall back to median if no time data
                noiseh = 10**((df_noise_clean[hh_col].median() - 60.0) / 10.)
                noisev = 10**((df_noise_clean[vv_col].median() - 60.0) / 10.)
                noisex = 10**((df_noise_clean[xx_col].median() - 60.0) / 10.)
                sdnoiseh = 10**(df_noise_clean[sd_hh_col].median() / 10.)
                sdnoisev = 10**(df_noise_clean[sd_vv_col].median() / 10.)
                sdnoisex = 10**(df_noise_clean[sd_xx_col].median() / 10.)
            
            # Get the reflectivity data (following camra_utils_old approach exactly)
            Zh = SweepDS.fields['DBZ_H']['data'].copy()
            Zv = Zh - SweepDS.fields['ZDR']['data'].copy()  # Vertically polarised reflectivity [dBZ] (copolar)
            Zx = Zh + SweepDS.fields['LDR']['data'].copy()

            linZh = 10.0**(Zh/10.)
            linZv = 10.0**(Zv/10.)
            linZx = 10.0**(Zx/10.)

            rng = SweepDS.range['data'].copy()
            
            rangem = ma.masked_where(rng <= 0.0, rng)
            rangekm = rangem / 1000.

            signalpowerh = linZh / rangem[None, :]**2
            signalpowerv = linZv / rangem[None, :]**2
            signalpowerx = linZx / rangem[None, :]**2

            signalpowerh -= noiseh  # subtract noise from signal
            signalpowerv -= noisev
            signalpowerx -= noisex

            # Calculate SNR
            SNRperpulseH = signalpowerh / noiseh
            SNRperpulseV = signalpowerv / noisev
            SNRperpulseX = signalpowerx / noisex

            # Get sweep slice indices
            sweep_slice = radar.get_slice(sweep)
            
            print(f"  Sweep slice: {sweep_slice}")
            
            # Convert back to reflectivity (no masking applied here)
            linZh_corrected = signalpowerh * rangem[None, :]**2
            linZv_corrected = signalpowerv * rangem[None, :]**2
            linZx_corrected = signalpowerx * rangem[None, :]**2

            # Only mask where signal power is negative (below noise floor)
            linZh_corrected = ma.masked_where(signalpowerh <= 0, linZh_corrected)
            linZv_corrected = ma.masked_where(signalpowerv <= 0, linZv_corrected)
            linZx_corrected = ma.masked_where(signalpowerx <= 0, linZx_corrected)

            Zh_corrected = 10. * ma.log10(linZh_corrected)
            Zv_corrected = 10. * ma.log10(linZv_corrected)
            Zx_corrected = 10. * ma.log10(linZx_corrected)

            # Update radar fields (noise-corrected but not quality-masked)
            radar.fields['DBZ_H']['data'][sweep_slice] = Zh_corrected
            radar.fields['ZDR']['data'][sweep_slice] = Zh_corrected - Zv_corrected
            radar.fields['LDR']['data'][sweep_slice] = Zx_corrected - Zh_corrected

            # Store SNR values
            radar.fields['SNR_H']['data'][sweep_slice] = 10.0 * ma.log10(ma.masked_where(SNRperpulseH <= 0, SNRperpulseH))
            radar.fields['SNR_V']['data'][sweep_slice] = 10.0 * ma.log10(ma.masked_where(SNRperpulseV <= 0, SNRperpulseV))
            radar.fields['SNR_X']['data'][sweep_slice] = 10.0 * ma.log10(ma.masked_where(SNRperpulseX <= 0, SNRperpulseX))

            # ============================================
            # SET QC FLAGS instead of masking data
            # ============================================
            
            # Create QC flag arrays for this sweep
            qc_dbz_h = np.zeros(Zh_corrected.shape, dtype=np.int8)
            qc_zdr = np.zeros(Zh_corrected.shape, dtype=np.int8)
            qc_ldr = np.zeros(Zh_corrected.shape, dtype=np.int8)
            
            # Flag 1: Low SNR H-pol
            low_snr_h = signalpowerh < SNR_threshold_HH * sdnoiseh
            qc_dbz_h[low_snr_h] += 1
            qc_zdr[low_snr_h] += 1
            qc_ldr[low_snr_h] += 1
            
            # Flag 2: Low SNR V-pol
            low_snr_v = signalpowerv < SNR_threshold_VV * sdnoisev
            qc_dbz_h[low_snr_v] += 2
            qc_zdr[low_snr_v] += 2
            
            # Flag 16: Low SNR X-pol (for LDR)
            low_snr_x = signalpowerx < SNR_threshold_X * sdnoisex
            qc_ldr[low_snr_x] += 16
            
            # Flag 4: Low correlation (if CXC field exists)
            if 'CXC' in radar.fields:
                cxc_data = radar.fields['CXC']['data'][sweep_slice]
                low_corr = cxc_data < 0.7  # correlation threshold
                qc_dbz_h[low_corr] += 4
                qc_zdr[low_corr] += 4
            
            # Flag 8: High spectrum width (if SPW_HV field exists)
            if 'SPW_HV' in radar.fields:
                spw_data = radar.fields['SPW_HV']['data'][sweep_slice]
                high_spw = spw_data > 5.0  # spectrum width threshold
                qc_dbz_h[high_spw] += 8
            
            # Update QC fields
            radar.fields['QC_DBZ_H']['data'][sweep_slice] = qc_dbz_h
            radar.fields['QC_ZDR']['data'][sweep_slice] = qc_zdr
            radar.fields['QC_LDR']['data'][sweep_slice] = qc_ldr
            
            # Set QC flags for other fields
            if 'VEL_HV' in radar.fields and 'QC_VEL_HV' in radar.fields:
                qc_vel_hv = np.zeros(Zh_corrected.shape, dtype=np.int8)
                qc_vel_hv[low_snr_h] += 1
                qc_vel_hv[low_snr_v] += 2
                if 'SPW_HV' in radar.fields:
                    qc_vel_hv[high_spw] += 8
                radar.fields['QC_VEL_HV']['data'][sweep_slice] = qc_vel_hv
            
            if 'CXC' in radar.fields and 'QC_CXC' in radar.fields:
                qc_cxc = np.zeros(Zh_corrected.shape, dtype=np.int8)
                qc_cxc[low_snr_h] += 1
                qc_cxc[low_snr_v] += 2
                radar.fields['QC_CXC']['data'][sweep_slice] = qc_cxc
            
            if 'PHIDP' in radar.fields and 'QC_PHIDP' in radar.fields:
                qc_phidp = np.zeros(Zh_corrected.shape, dtype=np.int8)
                qc_phidp[low_snr_h] += 1
                qc_phidp[low_snr_v] += 2
                if 'CXC' in radar.fields:
                    qc_phidp[low_corr] += 4
                radar.fields['QC_PHIDP']['data'][sweep_slice] = qc_phidp
            
            if 'SPW_HV' in radar.fields and 'QC_SPW_HV' in radar.fields:
                qc_spw = np.zeros(Zh_corrected.shape, dtype=np.int8)
                qc_spw[low_snr_h] += 1
                radar.fields['QC_SPW_HV']['data'][sweep_slice] = qc_spw

        # Apply only the basic cleaning functions (not doubleclean which masks more data)
        clean(radar)
        
        print("Noise subtraction with QC flags completed successfully")
        
    except Exception as e:
        print(f"Error in noise subtraction: {e}")
        import traceback
        traceback.print_exc()

def get_noise_levels_cfradial(RadarDS):
    """
    Get noise levels from CF-Radial radar data.
    
    Args:
        RadarDS: PyART radar object
        
    Returns:
        pandas.DataFrame: Noise levels data
    """
    import pandas as pd
    
    # ----------------
    # Process radar data
    # ----------------
    nsweeps = RadarDS.nsweeps

    noiseh_arr = []
    noisev_arr = []
    noisex_arr = []
    times_arr = []
    sdnoiseh_arr = []
    sdnoisev_arr = []
    sdnoisex_arr = []
    
    for sweep in range(nsweeps):

        SweepDS = RadarDS.extract_sweeps([sweep])
    
        x, y, height = SweepDS.get_gate_x_y_z(0)

        dtime = cftime.num2pydate(SweepDS.time['data'], SweepDS.time['units'])

        if len(dtime) > 0:
            timestamp = datetime.datetime.strftime(dtime[0], '%Y-%m-%dT%H:%M:%SZ')
        else:
            timestamp = ""

        Zh = SweepDS.fields['DBZ_H']['data']
        Zv = Zh - SweepDS.fields['ZDR']['data']  # Vertically polarised reflectivity [dBZ] (copolar)
        Zx = Zh + SweepDS.fields['LDR']['data'] 

        linZh = 10.0**(Zh/10.) 
        linZv = 10.0**(Zv/10.) 
        linZx = 10.0**(Zx/10.) 

        rng = SweepDS.range['data']
        
        rangem = ma.masked_where(rng <= 0.0, rng)
        rangekm = rangem/1000.;

        signalpowerh = linZh/rangem[None, :]**2
        signalpowerv = linZv/rangem[None, :]**2
        signalpowerx = linZx/rangem[None, :]**2

        masked_signalpowerh = np.ma.masked_where(height < 14000, signalpowerh)

        noiseh = np.ma.median(masked_signalpowerh)

        unmasked_indices = np.argwhere(~masked_signalpowerh.mask) 

        noisev = np.ma.median(np.ma.masked_where(height < 14000, signalpowerv))
        noisex = np.ma.median(np.ma.masked_where(height < 14000, signalpowerx))

        signalpowerh -= noiseh  # subtract noise from signal
        signalpowerv -= noisev
        signalpowerx -= noisex

        sdnoiseh = np.ma.std(np.ma.masked_where(height < 14000, signalpowerh))
        sdnoisev = np.ma.std(np.ma.masked_where(height < 14000, signalpowerv))
        sdnoisex = np.ma.std(np.ma.masked_where(height < 14000, signalpowerx))

        signalpowerh = ma.masked_where(signalpowerh < SNR_threshold_HH*sdnoiseh, signalpowerh) 
        signalpowerv = ma.masked_where(signalpowerv < SNR_threshold_VV*sdnoisev, signalpowerv) 
        signalpowerx = ma.masked_where(signalpowerx < SNR_threshold_X*sdnoisex, signalpowerx)

        SNRperpulseH = signalpowerh/noiseh 
        SNRperpulseV = signalpowerv/noisev
        SNRperpulseX = signalpowerx/noisex

        noiseh_arr.append(10.0*np.log10(noiseh)+60.)
        noisev_arr.append(10.0*np.log10(noisev)+60.)
        noisex_arr.append(10.0*np.log10(noisex)+60.)
        times_arr.append(timestamp)
        sdnoiseh_arr.append(10.0*np.log10(sdnoiseh))
        sdnoisev_arr.append(10.0*np.log10(sdnoisev))
        sdnoisex_arr.append(10.0*np.log10(sdnoisex))

    df = pd.DataFrame({
        'timestamp': times_arr,
        'noiseh': noiseh_arr,      # Changed from 'noise_hh'
        'noisev': noisev_arr,      # Changed from 'noise_vv' 
        'noisex': noisex_arr,      # Changed from 'noise_xx'
        'sdnoiseh': sdnoiseh_arr,  # Changed from 'sd_noise_hh'
        'sdnoisev': sdnoisev_arr,  # Changed from 'sd_noise_vv'
        'sdnoisex': sdnoisex_arr,  # Changed from 'sd_noise_xx'
    })
    
    return df

def calculate_noise_levels_from_files(
    rawfiles,
    noisefile_output,
    append_mode=True
):
    """
    Calculate noise levels from a list of CAMRa files and save to CSV.
    
    Args:
        rawfiles: List of raw CAMRa files
        noisefile_output: Path to output noise CSV file
        append_mode: If True, append to existing file; if False, overwrite
    
    Returns:
        pandas.DataFrame: Combined noise levels data
    """
    import pandas as pd
    
    all_noise_data = []
    
    print(f"Calculating noise levels from {len(rawfiles)} files...")
    
    for i, rawfile in enumerate(rawfiles):
        try:
            print(f"Processing file {i+1}/{len(rawfiles)}: {os.path.basename(rawfile)}")
            
            # Read the raw file
            RadarDS = read_camra_raw(rawfile)
            
            # Calculate noise levels for this file
            df_noise = get_noise_levels_cfradial(RadarDS)
            
            if not df_noise.empty:
                # Add file information
                df_noise['source_file'] = os.path.basename(rawfile)
                all_noise_data.append(df_noise)
                print(f"  Added {len(df_noise)} noise measurements")
            else:
                print(f"  No valid noise measurements found")
                
        except Exception as e:
            print(f"  Error processing {rawfile}: {e}")
            continue
    
    if all_noise_data:
        # Combine all noise data
        combined_df = pd.concat(all_noise_data, ignore_index=True)
        
        # Save to file
        if append_mode and os.path.exists(noisefile_output):
            # Load existing data and append
            try:
                existing_df = pd.read_csv(noisefile_output)
                combined_df = pd.concat([existing_df, combined_df], ignore_index=True)
                print(f"Appended to existing noise file with {len(existing_df)} previous measurements")
            except Exception as e:
                print(f"Warning: Could not load existing noise file: {e}")
        
        # Remove duplicates based on timestamp and source_file
        combined_df = combined_df.drop_duplicates(subset=['timestamp', 'source_file'], keep='last')
        
        # Save the combined data
        combined_df.to_csv(noisefile_output, index=False)
        print(f"Saved {len(combined_df)} total noise measurements to: {noisefile_output}")
        
        return combined_df
    else:
        print("No noise data calculated from any files")
        return pd.DataFrame()

def process_camra_day_noise_calculation(
    datestr: str,
    inpath: str,
    noisefile_output: str,
    file_pattern: str = "*rhi*.nc"
) -> str:
    """
    Process all RHI files for a day to calculate noise levels.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        noisefile_output: Path to output noise CSV file
        file_pattern: Glob pattern for files to process
    
    Returns:
        Path to noise file
    """
    print(f"=== NOISE CALCULATION for {datestr} ===")
    
    # Find all RHI files for the day
    search_pattern = os.path.join(inpath, file_pattern)
    rhi_files = glob.glob(search_pattern)
    
    if not rhi_files:
        print(f"No RHI files found matching pattern: {search_pattern}")
        return None
    
    print(f"Found {len(rhi_files)} RHI files for noise calculation")
    
    # Calculate noise levels
    calculate_noise_levels_from_files(
        rhi_files,
        noisefile_output,
        append_mode=True
    )
    
    return noisefile_output

def cfradial_add_ncas_metadata(cfradfile, yaml_project_file, yaml_instrument_file, tracking_tag, data_version):
    """
    Add NCAS metadata to CF-Radial file using YAML configuration files.
    
    Args:
        cfradfile: Path to CF-Radial NetCDF file
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        tracking_tag: AMOF tracking tag
        data_version: Data version string
        
    Returns:
        Path to renamed output file
    """
    import yaml
    import pathlib
    
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
            instrument = elem[instrument_tagname]

    # -------------------------------------
    # Read metadata from YAML projects file
    # -------------------------------------  
    with open(yaml_project_file, "r") as stream:
        try:
            projects = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    print(tracking_tag)
    print(yaml_project_file)
    print(yaml_instrument_file)

    for p in projects:
        if tracking_tag in p:
            project = p[tracking_tag]

    radar_name = instrument["instrument_name"].lower()
    print(radar_name)

    for n in project["ncas_instruments"]:
        if radar_name in n:
            project_instrument = n[radar_name]

    print(project_instrument)

    location = instrument['platform']['location'].lower();
    
    RadarDataset = nc4.Dataset(cfradfile)
    scan_name = RadarDataset.scan_name
    print(scan_name)

    str_start = RadarDataset.variables['time_coverage_start'][:].tobytes().decode('utf-8')
    time_coverage_start = str_start[0:20]

    str_end = RadarDataset.variables['time_coverage_end'][:].tobytes().decode('utf-8').strip()
    time_coverage_end = str_end[0:20]

    file_timestamp = datetime.datetime.strptime(time_coverage_start,'%Y-%m-%dT%H:%M:%SZ')

    dtstr = file_timestamp.strftime('%Y%m%d-%H%M%S');
    
    outpath = pathlib.Path(cfradfile).parent.resolve()
    outfile = os.path.join(outpath, f'{radar_name}_{location}_{dtstr}_{scan_name.replace("_", "-", 1)}_l1_v{data_version}.nc')

    if os.path.isfile(outfile):
        print("The file already exists")
    else:
        # Rename the file
        os.rename(cfradfile, outfile)
    
    RadarDataset.close()
    
    # -------------------------------------------------------
    # Read cfradial file to add NCAS metadata
    # -------------------------------------------------------
    DS = nc4.Dataset(outfile, 'r+')

    DS.Conventions = "NCAS-Radar-1.0 CfRadial-1.4 instrument_parameters radar_parameters geometry_correction"

    if 'version' in DS.ncattrs():
        DS.delncattr('version')

    DS.product_version = f"v{data_version}"
    DS.processing_level = "1"

    DS.licence = project_instrument["data_licence"]
    DS.acknowledgement = project_instrument["acknowledgement"]

    DS.platform = instrument["platform"]["location"]
    DS.platform_type = instrument["platform"]["type"]
    DS.location_keywords = instrument["platform"]["location_keywords"]
    DS.platform_is_mobile = "false"
    DS.deployment_mode = instrument["platform"]["deployment_mode"]

    DS.platform_altitude = instrument["platform"]["altitude"]

    DS.title = project_instrument["title"]
    DS.source = project_instrument["source"]

    DS.creator_name = project_instrument["data_creator"]["name"]
    DS.creator_email = project_instrument["data_creator"]["email"]
    DS.creator_url = project_instrument["data_creator"]["pid"]
    DS.institution = project_instrument["data_creator"]["institution"]
    DS.instrument_name = instrument["instrument_name"]
    DS.instrument_software = project_instrument["instrument_software"]["name"]
    DS.instrument_software_version = project_instrument["instrument_software"]["version"]
    DS.instrument_manufacturer = instrument['instrument_manufacturer']
    DS.instrument_model = instrument['instrument_model']
    DS.instrument_serial_number = instrument['instrument_serial_number']
    DS.instrument_pid = instrument['instrument_pid']

    DS.references = instrument['references']
    DS.comment = " "
    DS.project = project["project_name"]
    DS.project_principal_investigator = project["principal_investigator"]["name"]
    DS.project_principal_investigator_email = project["principal_investigator"]["email"]
    DS.project_principal_investigator_url = project["principal_investigator"]["pid"]

    DS.processing_software_url = "https://github.com/longlostjames/camra-radar-utils/releases/tag/v1.0.0"
    DS.processing_software_version = "v1.0.0"

    DS.time_coverage_start = time_coverage_start
    DS.time_coverage_end = time_coverage_end

    print('getting bbox')
    bbox_str = cfradial_get_bbox(outfile)
    print(bbox_str)
    DS.geospatial_bounds = bbox_str

    if "vpt" in DS.scan_name or "VPT" in DS.scan_name or "vertical_pointing" in DS.scan_name:
        DS.featureType = 'timeSeriesProfile'

    # -------------------------------------------------------
    # Now clean up some variable attributes
    # -------------------------------------------------------
    DS['time'].comment = ""
    try:
        DS['range'].delncattr('standard_name')
    except: 
        pass
    
    DS['range'].comment = 'Range to centre of each bin'
    DS['range'].meters_to_center_of_first_gate = DS['range'][0]
    try:
        DS['azimuth'].delncattr('standard_name')
    except:
        pass
    try:
        DS['elevation'].delncattr('standard_name')
    except:
        pass
    try:
        DS['sweep_number'].delncattr('standard_name')
    except:
        pass
    try:
        DS['sweep_mode'].delncattr('standard_name')
    except:
        pass
    try:
        DS['fixed_angle'].delncattr('standard_name')
    except:
        pass

    # Set the variable to a sequence of integers starting from 0  
    DS['sweep_number'][:] = np.arange(len(DS['sweep_number'][:]))  

    DS['latitude'].long_name = 'latitude'
    DS['latitude'].standard_name = 'latitude'
    DS['latitude'][:] = float(instrument["latitude"])
    DS['longitude'].long_name = 'longitude'
    DS['longitude'].standard_name = 'longitude'
    DS['longitude'][:] = float(instrument["longitude"])

    if DS['longitude'][0] < 0:
        DS['longitude'][:] += 360.0

    DS['altitude'].standard_name = 'altitude'
    DS['altitude'].comment = 'Altitude of the centre of rotation of the antenna above the geoid using the WGS84 ellipsoid and EGM2008 geoid model'
    DS['altitude'].long_name = 'altitude'
    DS['altitude'].units = 'metres'
    DS['altitude'][:] = float(instrument["altitude"]["value"])
    try:
        DS['altitude'].delncattr('positive')
    except:
        pass
    DS['volume_number'].long_name = 'data_volume_index_number'
    DS['volume_number'].units = ""

    altitude_agl = DS.createVariable('altitude_agl', 'f8')
    altitude_agl.assignValue(float(instrument["altitude_agl"]["value"]))
    altitude_agl.long_name = 'altitude_above_ground_level'
    altitude_agl.units = instrument['altitude_agl']['units']

    DS.close()

    # -----------------------
    # Update history metadata
    # -----------------------
    updatestr = "Add NCAS metadata"
    update_history_attribute(outfile, updatestr)

    return outfile

def cfradial_get_bbox(cfradfile):
    """
    Get bounding box from CF-Radial file.
    
    Args:
        cfradfile: Path to CF-Radial file
        
    Returns:
        String describing the bounding box
    """
    print(cfradfile)
    Radar = pyart.io.read_cfradial(cfradfile)
    latmin = np.min(Radar.gate_latitude['data'])
    lonmin = np.min(Radar.gate_longitude['data'])
    latmax = np.max(Radar.gate_latitude['data'])
    lonmax = np.max(Radar.gate_longitude['data'])
    print(latmin, latmax, lonmin, lonmax)
    boundingbox = f"Bounding box: {latmin:.4f}N {lonmin:.4f}E, {latmax:.4f}N {lonmax:.4f}E"
    return boundingbox
