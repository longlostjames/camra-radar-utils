#!/usr/bin/env python
# coding: utf-8

"""
campaign_processing.py

Campaign-specific processing functions for CAMRa radar data

This module provides processing functions for different CAMRa campaigns:
- WOEST (WESCON Observing the Evolving Structures of Turbulence) 
  - IOP (Intensive Observation Period)
  - SOP (Standard Observation Period) 
  - Other (Other operations)
- KASBEX (Ka- and S-Band EXperiment)
- CCREST-M (Characterising CiRrus and icE cloud acrosS the specTrum - Microwave)
  - Standard RHI processing at specific azimuths
  - VPT (Vertically Pointing) processing
  - VPT time series processing
- DYMECS (DYnamics, Microphysics and Entrainment in Convective Systems)
  - RHI and PPI processing

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
Last modified: 22-02-2026
Version: 1.3.0
"""

from typing import List, Dict, Tuple, Optional, Any
import datetime
import os
import glob
import yaml
import netCDF4 as nc4

# Import CAMRa utilities
import camra_utils

# ==============================================================================
# WOEST-SPECIFIC PROCESSING FUNCTIONS
# ==============================================================================

def process_camra_woest_iop_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.1",
    tracking_tag: str = "AMOF_20220922221548",
    campaign: str = "woest",
    single_sweep: bool = True  # Default to single sweep for IOP
) -> None:
    """
    Process CAMRa WOEST IOP (Intensive Observation Period) data.
    
    IOP processing handles regional storm tracking files organized by region.
    Files are typically located in subdirectories like region_01, region_02, etc.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string
        tracking_tag: AMOF tracking tag for WOEST CAMRa data
        campaign: Campaign name
        single_sweep: If True, process each file separately (recommended for IOP)
    """
    print(f"Processing CAMRa WOEST IOP day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # IOP files are organized in region directories
    ioppath = os.path.join(inpath, 'iop')
    
    if not os.path.exists(ioppath):
        print(f"IOP directory does not exist: {ioppath}")
        return
    
    # Find regional directories
    regions = [os.path.basename(elem) for elem in glob.glob(os.path.join(ioppath, 'region_*'))]
    print(f"Found regions: {regions}")
    
    if not regions:
        print("No regional directories found for IOP processing")
        return
    
    # Process each region
    result = {}
    for region in regions:
        # Look for RHI files in this region
        pattern = os.path.join(ioppath, region, '*rhi*.nc')
        found_files = glob.glob(pattern)
        
        if found_files:
            region_key = region.split('_')[1]  # Extract region number
            result[region_key] = found_files
            print(f"Region {region_key}: Found {len(found_files)} RHI files")
    
    # Process files for each region
    for region_key, iop_files in result.items():
        print(f"Processing region {region_key} with {len(iop_files)} files")
        
        if single_sweep:
            # Process each file individually (recommended for IOP)
            for file_path in iop_files:
                try:
                    print(f"Processing IOP file: {os.path.basename(file_path)}")
                    radar_ds_path = camra_utils.multi_camra2cfrad(
                        [file_path],  # Single file in list
                        outdir,
                        scan_name='RHI',
                        data_version=data_version,
                        tracking_tag=tracking_tag,
                        campaign=campaign,
                        yaml_project_file=yaml_project_file,
                        yaml_instrument_file=yaml_instrument_file
                    )
                    
                    # Add region-specific comment to the file
                    if radar_ds_path and os.path.exists(radar_ds_path):
                        with nc4.Dataset(radar_ds_path, 'r+') as ds:
                            ds.comment = f"WOEST IOP tracking storm region {region_key}"
                        print(f"Added region comment to {os.path.basename(radar_ds_path)}")
                        
                except Exception as e:
                    print(f"Error processing IOP file {file_path}: {e}")
        else:
            # Process all files in region concatenated
            try:
                radar_ds_path = camra_utils.multi_camra2cfrad(
                    iop_files,
                    outdir,
                    scan_name='SRHI',  # Sequential RHI
                    data_version=data_version,
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file
                )
                
                # Add region-specific comment
                if radar_ds_path and os.path.exists(radar_ds_path):
                    with nc4.Dataset(radar_ds_path, 'r+') as ds:
                        ds.comment = f"WOEST IOP tracking storm region {region_key}"
                    print(f"Added region comment to concatenated file")
                    
            except Exception as e:
                print(f"Error processing concatenated IOP files for region {region_key}: {e}")

def process_camra_woest_sop_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.1",
    tracking_tag: str = "AMOF_20220922221548",
    campaign: str = "woest",
    single_sweep: bool = False
) -> None:
    """
    Process CAMRa WOEST SOP (Standard Observation Period) data.
    
    SOP processing handles special observation files, typically with
    more intensive scanning strategies.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string
        tracking_tag: AMOF tracking tag for WOEST CAMRa data
        campaign: Campaign name
        single_sweep: If True, process each file separately
    """
    print(f"Processing CAMRa WOEST SOP day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # SOP files may be in a specific subdirectory
    sop_paths = [
        os.path.join(inpath, 'sop'),  # Dedicated SOP directory
        inpath  # Or in main directory with SOP naming
    ]
    
    files_by_scan = {}
    
    # Search for SOP files in potential locations
    for search_path in sop_paths:
        if os.path.exists(search_path):
            print(f"Searching for SOP files in: {search_path}")
            
            # Look for files with SOP indicators or enhanced scan patterns
            sop_files = find_camra_files(search_path, datestr, file_filter='sop')
            
            # Merge results
            for scan_type, files in sop_files.items():
                if scan_type not in files_by_scan:
                    files_by_scan[scan_type] = []
                files_by_scan[scan_type].extend(files)
    
    if not files_by_scan:
        print(f"No CAMRa SOP files found for {datestr}")
        return
    
    # Process each scan type with SOP-specific handling
    scan_types = ['rhi', 'ppi', 'man', 'fix']
    
    for scan_type in scan_types:
        if scan_type in files_by_scan:
            try:
                process_camra_scan_type(
                    files_by_scan[scan_type],
                    scan_type,
                    outdir,
                    yaml_project_file,
                    yaml_instrument_file,
                    data_version,
                    tracking_tag,
                    campaign,
                    single_sweep,
                    comment_suffix="SOP"  # Add SOP identifier
                )
            except Exception as e:
                print(f"Error processing CAMRa SOP {scan_type.upper()}: {e}")

def process_camra_woest_other_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.1",
    tracking_tag: str = "AMOF_20220922221548",
    campaign: str = "woest",
    single_sweep: bool = False
) -> None:
    """
    Process CAMRa WOEST Other (regular operations) data.
    
    Other processing handles routine operational scanning data
    that doesn't fall into IOP or SOP categories.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string
        tracking_tag: AMOF tracking tag for WOEST CAMRa data
        campaign: Campaign name
        single_sweep: If True, process each file separately
    """
    print(f"Processing CAMRa WOEST Other day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Find regular operational files (excluding IOP and SOP)
    files_by_scan = find_camra_files(inpath, datestr, exclude_filters=['iop', 'sop'])
    
    if not files_by_scan:
        print(f"No CAMRa Other files found for {datestr}")
        return
    
    # Process each scan type with standard handling
    scan_types = ['rhi', 'ppi', 'man', 'fix']
    
    for scan_type in scan_types:
        if scan_type in files_by_scan:
            try:
                process_camra_scan_type(
                    files_by_scan[scan_type],
                    scan_type,
                    outdir,
                    yaml_project_file,
                    yaml_instrument_file,
                    data_version,
                    tracking_tag,
                    campaign,
                    single_sweep
                )
            except Exception as e:
                print(f"Error processing CAMRa Other {scan_type.upper()}: {e}")

def process_camra_woest_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.1",
    tracking_tag: str = "AMOF_20220922221548",
    campaign: str = "woest",
    single_sweep: bool = False,
    woest_mode: str = "auto"
) -> None:
    """
    Process CAMRa WOEST campaign data with automatic mode detection.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string
        tracking_tag: AMOF tracking tag for WOEST CAMRa data
        campaign: Campaign name
        single_sweep: If True, process each file separately
        woest_mode: Processing mode ('auto', 'iop', 'sop', 'other')
    """
    print(f"Processing CAMRa WOEST day: {datestr} (mode: {woest_mode})")
    
    if woest_mode == "auto":
        # Auto-detect processing mode based on directory structure
        iop_path = os.path.join(inpath, 'iop')
        sop_path = os.path.join(inpath, 'sop')
        
        if os.path.exists(iop_path) and glob.glob(os.path.join(iop_path, 'region_*')):
            print("Auto-detected IOP mode")
            process_camra_woest_iop_step1(
                datestr, inpath, outpath, yaml_project_file, yaml_instrument_file,
                data_version, tracking_tag, campaign, single_sweep
            )
        elif os.path.exists(sop_path) or any('sop' in f.lower() for f in glob.glob(os.path.join(inpath, '*.nc'))):
            print("Auto-detected SOP mode")
            process_camra_woest_sop_step1(
                datestr, inpath, outpath, yaml_project_file, yaml_instrument_file,
                data_version, tracking_tag, campaign, single_sweep
            )
        else:
            print("Auto-detected Other mode")
            process_camra_woest_other_step1(
                datestr, inpath, outpath, yaml_project_file, yaml_instrument_file,
                data_version, tracking_tag, campaign, single_sweep
            )
    elif woest_mode == "iop":
        process_camra_woest_iop_step1(
            datestr, inpath, outpath, yaml_project_file, yaml_instrument_file,
            data_version, tracking_tag, campaign, single_sweep
        )
    elif woest_mode == "sop":
        process_camra_woest_sop_step1(
            datestr, inpath, outpath, yaml_project_file, yaml_instrument_file,
            data_version, tracking_tag, campaign, single_sweep
        )
    elif woest_mode == "other":
        process_camra_woest_other_step1(
            datestr, inpath, outpath, yaml_project_file, yaml_instrument_file,
            data_version, tracking_tag, campaign, single_sweep
        )
    else:
        raise ValueError(f"Unknown WOEST mode: {woest_mode}. Use 'auto', 'iop', 'sop', or 'other'")

# ==============================================================================
# KASBEX PROCESSING FUNCTIONS
# ==============================================================================

def process_camra_kasbex_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.1",
    tracking_tag: str = "AMOF_20250508133639",
    campaign: str = "kasbex",
    single_sweep: bool = False,
    noise_file: str = None
) -> None:
    """
    Process CAMRa KASBEX campaign data for a single day - Step 1.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string
        tracking_tag: AMOF tracking tag for KASBEX CAMRa data
        campaign: Campaign name
        single_sweep: If True, process each file separately
    """
    print(f"Processing CAMRa KASBEX day: {datestr}")
    print(f"Data version: {data_version}")
    print(f"Tracking tag: {tracking_tag}")
    print(f"Single sweep mode: {single_sweep}")
    
    # Debug noise file parameter
    if noise_file:
        print(f"Noise file parameter received: {noise_file}")
        print(f"Noise file exists: {os.path.exists(noise_file)}")
    else:
        print("No noise file parameter received")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Check if inpath already includes the date (i.e., when using --latest)
    path_parts = inpath.rstrip('/').split('/')
    if datestr in path_parts or 'latest' in path_parts:
        search_path = inpath
        print(f"Using direct path for file search: {search_path}")
    else:
        search_path = os.path.join(inpath, datestr)
        print(f"Using traditional date-based path: {search_path}")
    
    # Find CAMRa files for this date
    files_by_scan = find_camra_files(search_path, datestr)
    
    if not files_by_scan:
        print(f"No CAMRa files found for KASBEX {datestr}")
        return
    
    # Process each scan type
    scan_types = ['rhi', 'ppi', 'man', 'fix']
    
    for scan_type in scan_types:
        if scan_type in files_by_scan:
            try:
                print(f"Processing {scan_type.upper()} with noise_file: {noise_file}")
                process_camra_scan_type(
                    files_by_scan[scan_type],
                    scan_type,
                    outdir,
                    yaml_project_file,
                    yaml_instrument_file,
                    data_version,
                    tracking_tag,
                    campaign,
                    single_sweep,
                    noise_file=noise_file
                )
            except Exception as e:
                print(f"Error processing CAMRa {scan_type.upper()}: {e}")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def find_camra_files(
    search_path: str, 
    datestr: str, 
    file_filter: Optional[str] = None,
    exclude_filters: Optional[List[str]] = None
) -> Dict[str, List[str]]:
    """
    Find CAMRa files for the given date and path with optional filtering.
    
    Args:
        search_path: Directory to search for files
        datestr: Date string in YYYYMMDD format
        file_filter: Include only files containing this string
        exclude_filters: Exclude files containing any of these strings
        
    Returns:
        Dictionary mapping scan types to lists of files
    """
    if not os.path.exists(search_path):
        print(f"Search path does not exist: {search_path}")
        return {}
    
    # Search patterns for different scan types
    scan_patterns = {
        'rhi': [f'*{datestr}*rhi*.nc', f'*{datestr}*RHI*.nc'],
        'ppi': [f'*{datestr}*ppi*.nc', f'*{datestr}*PPI*.nc'],
        'man': [f'*{datestr}*man*.nc', f'*{datestr}*MAN*.nc', f'*{datestr}*manual*.nc'],
        'fix': [f'*{datestr}*fix*.nc', f'*{datestr}*FIX*.nc']
    }
    
    files_by_scan = {}
    
    for scan_type, patterns in scan_patterns.items():
        files = []
        for pattern in patterns:
            files.extend(glob.glob(os.path.join(search_path, pattern)))
        
        # Apply filters
        if file_filter:
            files = [f for f in files if file_filter.lower() in f.lower()]
        
        if exclude_filters:
            for exclude_filter in exclude_filters:
                files = [f for f in files if exclude_filter.lower() not in f.lower()]
        
        if files:
            files.sort()
            files_by_scan[scan_type] = files
            print(f"Found {len(files)} {scan_type.upper()} files")
        else:
            print(f"No {scan_type.upper()} files found")
    
    return files_by_scan


# CCREST-M PROCESSING FUNCTIONS
# ===============================

def process_camra_ccrest_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.0",
    tracking_tag: str = "AMOF_20230201132601",
    campaign: str = "ccrest-m"
) -> None:
    """
    Process CAMRa CCREST-M campaign data for a single day - Step 1.
    
    Processes RHI scans at two different azimuth angles for the CCREST-M campaign:
    - RHI-CCREST1: 270° azimuth (range 260-280°)
    - RHI-CCREST2: 246° azimuth (range 236-256°)
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path containing raw CAMRa files
        outpath: Output directory path for processed files
        yaml_project_file: Path to YAML project configuration file
        yaml_instrument_file: Path to YAML instrument configuration file
        data_version: Data version string (default: "1.0.0")
        tracking_tag: AMOF tracking tag (default: CCREST-M tag)
        campaign: Campaign name (default: "ccrest-m")
    """
    print(f"Processing CAMRa CCREST-M day: {datestr}")
    
    # Define the start and end times for the day
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    # Find RHI files at 270° azimuth (CCREST1)
    rhi_files_270 = camra_utils.find_camra_rhi_files(
        start_date.strftime('%Y-%m-%d %H:%M:%S'),
        end_date.strftime('%Y-%m-%d %H:%M:%S'),
        260, 280, inpath
    )
    
    # Find RHI files at 246° azimuth (CCREST2)
    rhi_files_246 = camra_utils.find_camra_rhi_files(
        start_date.strftime('%Y-%m-%d %H:%M:%S'),
        end_date.strftime('%Y-%m-%d %H:%M:%S'),
        236, 256, inpath
    )
    
    print(f"RHI files at 270°: {len(rhi_files_270)}")
    print(f"RHI files at 246°: {len(rhi_files_246)}")
    
    # Process RHI-CCREST1 (270° azimuth)
    if len(rhi_files_270) > 0:
        print("Processing RHI-CCREST1 (270° azimuth)")
        RadarDS_RHI1 = camra_utils.multi_camra2cfrad(
            rhi_files_270,
            outpath,
            scan_name='RHI-CCREST1',
            data_version=data_version,
            tracking_tag=tracking_tag,
            campaign=campaign
        )
    else:
        print("No RHI files found for CCREST1 (270° azimuth)")
    
    # Process RHI-CCREST2 (246° azimuth)
    if len(rhi_files_246) > 0:
        print("Processing RHI-CCREST2 (246° azimuth)")
        RadarDS_RHI2 = camra_utils.multi_camra2cfrad(
            rhi_files_246,
            outpath,
            scan_name='RHI-CCREST2',
            data_version=data_version,
            tracking_tag=tracking_tag,
            campaign=campaign
        )
    else:
        print("No RHI files found for CCREST2 (246° azimuth)")
    
    print(f"Completed CCREST-M processing for {datestr}")


def process_camra_ccrest_vpt_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.0",
    tracking_tag: str = "AMOF_20230201132601",
    campaign: str = "ccrest-m"
) -> None:
    """
    Process CAMRa CCREST-M VPT (Vertically Pointing) data for a single day.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path containing raw CAMRa VPT files
        outpath: Output directory path for processed files
        yaml_project_file: Path to YAML project configuration file
        yaml_instrument_file: Path to YAML instrument configuration file
        data_version: Data version string (default: "1.0.0")
        tracking_tag: AMOF tracking tag (default: CCREST-M tag)
        campaign: Campaign name (default: "ccrest-m")
    """
    print(f"Processing CAMRa CCREST-M VPT day: {datestr}")
    
    # Define the start and end times for the day
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    # Find VPT files for the day
    try:
        vpt_files = camra_utils.find_camra_vpt_files(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            inpath
        )
        print(f"Found {len(vpt_files)} VPT files")
        
        if len(vpt_files) > 0:
            RadarDS_VPT = camra_utils.multi_camra2cfrad(
                vpt_files,
                outpath,
                scan_name='fix',
                data_version=data_version,
                tracking_tag=tracking_tag,
                campaign=campaign
            )
            print(f"Successfully processed VPT data for {datestr}")
        else:
            print(f"No VPT files found for {datestr}")
            
    except Exception as e:
        print(f"VPT processing error for {datestr}: {e}")


def process_camra_ccrest_vpt_day_ts(
    datestr: str,
    ts_indir: str,
    mom_indir: str,
    outdir: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.0",
    tracking_tag: str = "AMOF_20230201132601"
) -> None:
    """
    Process CAMRa CCREST-M VPT time series data for a single day.
    
    This function processes VPT time series data through the following pipeline:
    1. Convert L0a to L0b
    2. Convert L0b to L1
    3. Convert L1 time series to moments
    
    Args:
        datestr: Date string in YYYYMMDD format
        ts_indir: Input directory path for time series L0a files
        mom_indir: Input directory path for moment files
        outdir: Output directory path for processed files
        yaml_project_file: Path to YAML project configuration file
        yaml_instrument_file: Path to YAML instrument configuration file
        data_version: Data version string (default: "1.0.0")
        tracking_tag: AMOF tracking tag (default: CCREST-M tag)
    """
    print(f"Processing CAMRa CCREST-M VPT time series day: {datestr}")
    
    # Define the start and end times for the day
    start_date = datetime.datetime.strptime(datestr, '%Y%m%d')
    end_date = start_date + datetime.timedelta(days=1)
    
    # Set up output paths
    l0bpath = os.path.join(outdir, 'ts', 'L0b', datestr)
    l1path = os.path.join(outdir, 'ts', 'L1', datestr)
    mom_l1a_path = os.path.join(mom_indir, datestr)
    mom_l1b_path = os.path.join(outdir, 'L1b', datestr)
    
    # Step 1: Convert L0a to L0b
    try:
        vpt_ts_files_l0a = camra_utils.find_camra_vpt_ts_files(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            ts_indir,
            level='l0a'
        )
        print(f"Found {len(vpt_ts_files_l0a)} L0a files")
        
        if len(vpt_ts_files_l0a) > 0:
            for f in vpt_ts_files_l0a:
                camra_utils.convert_camra_ts_l0a2l0b(
                    f, l0bpath,
                    tracking_tag=tracking_tag,
                    data_version=data_version
                )
            print("Completed L0a to L0b conversion")
    except Exception as e:
        print(f"VPT L0a to L0b conversion error: {e}")
    
    # Step 2: Convert L0b to L1
    try:
        vpt_ts_files_l0b = camra_utils.find_camra_vpt_ts_files(
            start_date.strftime('%Y-%m-%d %H:%M:%S'),
            end_date.strftime('%Y-%m-%d %H:%M:%S'),
            l0bpath,
            level='l0b'
        )
        print(f"Found {len(vpt_ts_files_l0b)} L0b files")
        
        if len(vpt_ts_files_l0b) > 0:
            for f in vpt_ts_files_l0b:
                print(f"Converting: {f}")
                camra_utils.convert_camra_ts_l0b2l1(
                    f, l1path,
                    tracking_tag=tracking_tag,
                    data_version=data_version
                )
            print("Completed L0b to L1 conversion")
    except Exception as e:
        print(f"VPT L0b to L1 conversion error: {e}")
    
    # Step 3: Convert L1 time series to moments
    try:
        camra_utils.convert_camra_tsl1tomoments(l1path, mom_l1a_path, mom_l1b_path)
        print("Completed L1 to moments conversion")
    except Exception as e:
        print(f"VPT L1 to moments conversion error: {e}")


# ==============================================================================
# DYMECS PROCESSING FUNCTIONS
# ==============================================================================

def process_camra_dymecs_day_step1(
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = "1.0.0",
    tracking_tag: str = "CRF_85",
    campaign: str = "dymecs",
    single_sweep: bool = True
) -> None:
    """
    Process CAMRa DYMECS campaign data for a single day - Step 1.
    
    DYMECS (DYnamics, Microphysics and Entrainment in Convective Systems)
    processing handles RHI and PPI scans.
    
    Args:
        datestr: Date string in YYYYMMDD format
        inpath: Input directory containing CAMRa files
        outpath: Output directory for processed files
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string (default: "1.0.0")
        tracking_tag: AMOF tracking tag (default: "CRF_85")
        campaign: Campaign name (default: "dymecs")
        single_sweep: If True, process each file separately (default: True)
    """
    print(f"Processing CAMRa DYMECS day: {datestr}")
    
    # Create output directory
    outdir = os.path.join(outpath, datestr)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Find CAMRa files for this date
    files_by_scan = find_camra_files(inpath, datestr)
    
    if not files_by_scan:
        print(f"No CAMRa files found for DYMECS {datestr}")
        return
    
    # Process RHI and PPI scan types
    scan_types = ['rhi', 'ppi']
    
    for scan_type in scan_types:
        if scan_type in files_by_scan:
            try:
                process_camra_scan_type(
                    files_by_scan[scan_type],
                    scan_type,
                    outdir,
                    yaml_project_file,
                    yaml_instrument_file,
                    data_version,
                    tracking_tag,
                    campaign,
                    single_sweep
                )
            except Exception as e:
                print(f"Error processing CAMRa DYMECS {scan_type.upper()}: {e}")


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def process_camra_scan_type(
    files: List[str],
    scan_type: str,
    outdir: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str,
    tracking_tag: str,
    campaign: str,
    single_sweep: bool,
    comment_suffix: Optional[str] = None,
    noise_file: Optional[str] = None
) -> None:
    """
    Process CAMRa files for a specific scan type.
    """
    if not files:
        print(f"No {scan_type.upper()} files to process")
        return
    
    print(f"Processing {len(files)} {scan_type.upper()} files")
    
    # Debug noise file parameter
    if noise_file:
        print(f"Scan processing received noise_file: {noise_file}")
    else:
        print("Scan processing: no noise file specified")
    
    if single_sweep:
        # Process each file individually
        for file_path in files:
            print(f"Processing single file: {os.path.basename(file_path)}")
            
            try:
                print(f"Calling multi_camra2cfrad with noise_file: {noise_file}")
                radar_ds_path = camra_utils.multi_camra2cfrad(
                    [file_path],
                    outdir,
                    scan_name=scan_type.upper(),
                    tracking_tag=tracking_tag,
                    campaign=campaign,
                    data_version=data_version,
                    yaml_project_file=yaml_project_file,
                    yaml_instrument_file=yaml_instrument_file,
                    noise_file=noise_file
                )
                
                # Add comment suffix if provided
                if comment_suffix and radar_ds_path and os.path.exists(radar_ds_path):
                    with nc4.Dataset(radar_ds_path, 'r+') as ds:
                        original_comment = getattr(ds, 'comment', '')
                        ds.comment = f"{original_comment} {comment_suffix}".strip()
                
                print(f"Successfully processed single {scan_type.upper()} file: {os.path.basename(file_path)}")
            except Exception as e:
                print(f"Error processing single {scan_type.upper()} file {file_path}: {e}")
                import traceback
                traceback.print_exc()
    else:
        # Process multiple files concatenated
        print(f"Processing {len(files)} {scan_type.upper()} files concatenated")
        
        try:
            radar_ds_path = camra_utils.multi_camra2cfrad(
                files,
                outdir,
                scan_name=scan_type.upper(),
                tracking_tag=tracking_tag,
                campaign=campaign,
                data_version=data_version,
                yaml_project_file=yaml_project_file,
                yaml_instrument_file=yaml_instrument_file,
                noise_file=noise_file  # Pass noise file to processing
            )
            
            # Add comment suffix if provided
            if comment_suffix and radar_ds_path and os.path.exists(radar_ds_path):
                with nc4.Dataset(radar_ds_path, 'r+') as ds:
                    original_comment = getattr(ds, 'comment', '')
                    ds.comment = f"{original_comment} {comment_suffix}".strip()
            
            print(f"Successfully processed {len(files)} {scan_type.upper()} files concatenated")
        except Exception as e:
            print(f"Error processing concatenated {scan_type.upper()} files: {e}")

# ==============================================================================
# UNIFIED CAMPAIGN PROCESSING FUNCTIONS
# ==============================================================================

def get_campaign_info(campaign: str) -> Dict[str, Any]:
    """
    Get campaign-specific configuration parameters for CAMRa radar.
    
    Args:
        campaign: Campaign name
        
    Returns:
        Dictionary of campaign-specific parameters
    """
    campaign_configs = {
        'woest': {
            'tracking_tag': 'AMOF_20220922221548',
            'data_version': '1.0.1',
            'modes': ['auto', 'iop', 'sop', 'other']
        },
        'kasbex': {
            'tracking_tag': 'AMOF_20250508133639',
            'data_version': '1.0.1'
        }
    }
    
    return campaign_configs.get(campaign, {
        'tracking_tag': f'AMOF_{campaign.upper()}',
        'data_version': '1.0.1'
    })

def load_yaml_config(yaml_project_file: str, yaml_instrument_file: str) -> Dict[str, Any]:
    """
    Load configuration from YAML files for CAMRa radar.
    
    Args:
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        
    Returns:
        Dictionary of configuration parameters
    """
    config = {
        'data_version': '1.0.1',
        'location': 'cao'
    }
    
    # Load project YAML
    try:
        with open(yaml_project_file, 'r') as f:
            project_data = yaml.safe_load(f)
        
        print(f"Loading configuration from: {yaml_project_file}")
        
        # Navigate through YAML structure for CAMRa radar
        if 'ncas_instruments' in project_data:
            for instrument in project_data['ncas_instruments']:
                if 'ncas-radar-camra-1' in instrument:
                    radar_config = instrument['ncas-radar-camra-1']
                    
                    # Extract processing software version
                    if 'processing_software' in radar_config and 'version' in radar_config['processing_software']:
                        config['data_version'] = radar_config['processing_software']['version']
                    
                    # Extract location
                    if 'platform' in radar_config and 'location' in radar_config['platform']:
                        config['location'] = radar_config['platform']['location'].lower()
                    
                    break
        
    except Exception as e:
        print(f"Warning: Could not load project YAML config: {e}")
        print(f"Using default data_version: {config['data_version']}")
    
    return config

def process_campaign_day(
    campaign: str,
    datestr: str,
    inpath: str,
    outpath: str,
    yaml_project_file: str,
    yaml_instrument_file: str,
    data_version: str = None,
    single_sweep: bool = False,
    noise_file: str = None,  # Add noise file parameter
    **kwargs
) -> None:
    """
    Process a single day of campaign data for CAMRa radar.
    
    Args:
        campaign: Campaign name ('woest', 'kasbex')
        datestr: Date string in YYYYMMDD format
        inpath: Input directory path
        outpath: Output directory path
        yaml_project_file: Path to project YAML file
        yaml_instrument_file: Path to instrument YAML file
        data_version: Data version string (if None, uses campaign default)
        single_sweep: If True, create separate files for each sweep
        noise_file: Path to noise CSV file (optional)
        **kwargs: Additional campaign-specific arguments (e.g., woest_mode)
    """
    campaign_lower = campaign.lower()
    
    # Get campaign-specific configuration
    campaign_info = get_campaign_info(campaign_lower)
    
    # Set default data version if not provided
    if data_version is None:
        data_version = campaign_info.get('data_version', '1.0.1')
    
    # Define available processors
    processors = {
        'woest': process_camra_woest_day_step1,
        'kasbex': process_camra_kasbex_day_step1
    }
    
    if campaign_lower not in processors:
        available_campaigns = list(processors.keys())
        raise ValueError(f"Unknown campaign '{campaign}' for CAMRa radar. Available: {available_campaigns}")
    
    # Set up processing arguments
    processing_args = {
        'datestr': datestr,
        'inpath': inpath,
        'outpath': outpath,
        'yaml_project_file': yaml_project_file,
        'yaml_instrument_file': yaml_instrument_file,
        'data_version': data_version,
        'single_sweep': single_sweep,
        'tracking_tag': campaign_info.get('tracking_tag', f'AMOF_{campaign.upper()}'),
        'campaign': campaign_lower,
        'noise_file': noise_file,  # Pass noise file to processing functions
        **kwargs
    }
    
    print(f"Processing {campaign.upper()} campaign data for CAMRa radar - {datestr}")
    print(f"Configuration:")
    print(f"  Data version: {data_version}")
    print(f"  Single sweep: {single_sweep}")
    print(f"  Tracking tag: {processing_args['tracking_tag']}")
    
    # Add campaign-specific parameter info
    if campaign_lower == 'woest' and 'woest_mode' in kwargs:
        print(f"  WOEST mode: {kwargs['woest_mode']}")
    
    # Call the appropriate processor
    processor = processors[campaign_lower]
    processor(**processing_args)
    
    print(f"Completed {campaign.upper()} processing for CAMRa radar - {datestr}")