#!/usr/bin/env python3
"""
proc_camra_kasbex_campaign_batch.py

Batch process KASBEX campaign data over a date range using campaign_processing module for CAMRA radar.

Usage:
    python proc_camra_kasbex_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD
    python proc_camra_kasbex_campaign_batch.py -d YYYYMMDD
    python proc_camra_kasbex_campaign_batch.py -d YYYYMMDD --latest

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
"""

import sys
import os
import argparse
import datetime
from pathlib import Path
import glob

# Add the script directory to Python path
script_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(script_dir))

from campaign_processing import process_campaign_day, get_campaign_info

def setup_kasbex_paths(use_latest=False):
    """Set up file and directory paths for CAMRA KASBEX campaign."""
    
    # Base KASBEX campaign paths for CAMRA radar
    base_inpath = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-camra-1/data/campaign/kasbex/raw'
    
    # Modify input path if using latest subdirectory
    if use_latest:
        print("Using 'latest' subdirectory for input data")
    
    # Use YAML files from the repository campaigns/ directory
    paths = {
        'inpath': base_inpath,
        'outpath': '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1/kasbex/L1',
        'yaml_project_file': script_dir / 'campaigns' / 'kasbex_project.yml',
        'yaml_instrument_file': script_dir / 'instrument_metadata.yml',
        'use_latest': use_latest
    }
    
    # Ensure output directory exists
    Path(paths['outpath']).mkdir(parents=True, exist_ok=True)
    
    return paths

def generate_date_list(start_date_str, end_date_str):
    """
    Generate a list of date strings between start and end dates (inclusive).
    
    Args:
        start_date_str: Start date in YYYYMMDD format
        end_date_str: End date in YYYYMMDD format
        
    Returns:
        List of date strings in YYYYMMDD format
    """
    try:
        start_date = datetime.datetime.strptime(start_date_str, '%Y%m%d')
        end_date = datetime.datetime.strptime(end_date_str, '%Y%m%d')
    except ValueError as e:
        raise ValueError(f"Invalid date format. Use YYYYMMDD. Error: {e}")
    
    if start_date > end_date:
        raise ValueError("Start date must be before or equal to end date")
    
    date_list = []
    current_date = start_date
    
    while current_date <= end_date:
        date_list.append(current_date.strftime('%Y%m%d'))
        current_date += datetime.timedelta(days=1)
    
    return date_list

def check_input_data(inpath, datestr, use_latest=False):
    """
    Check if input data exists for the specified date.
    
    Args:
        inpath: Input directory path
        datestr: Date string in YYYYMMDD format
        use_latest: If True, look in datestr/latest subdirectory
        
    Returns:
        bool: True if data exists, False otherwise
    """
    if use_latest:
        date_path = Path(inpath) / datestr / 'latest'
        print(f"Looking for data in LATEST subdirectory: {date_path}")
    else:
        date_path = Path(inpath) / datestr
        print(f"Looking for data in: {date_path}")
    
    if not date_path.exists():
        print(f"Date directory does not exist: {date_path}")
        return False
    
    # Check for CAMRa .nc files
    nc_files = list(date_path.glob('*.nc'))
    print(f"Found {len(nc_files)} nc files in {date_path}")
    
    if nc_files:
        print(f"Sample files: {[f.name for f in nc_files[:3]]}")
    else:
        # Debug: show what files are actually there
        all_files = list(date_path.glob('*'))
        print(f"All files in directory: {[f.name for f in all_files[:10]]}")
    
    return len(nc_files) > 0

def get_input_path(base_inpath, datestr, use_latest=False):
    """
    Get the appropriate input path for the given date.
    
    Args:
        base_inpath: Base input directory path
        datestr: Date string in YYYYMMDD format
        use_latest: If True, use datestr/latest subdirectory
        
    Returns:
        str: Full input path for the date
    """
    if use_latest:
        return os.path.join(base_inpath, datestr, 'latest')
    else:
        return os.path.join(base_inpath, datestr)

def main():
    """Main batch processing function for CAMRA KASBEX."""
    parser = argparse.ArgumentParser(
        description="Batch process CAMRA KASBEX campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--start-date', 
        required=False,
        help='Start date in YYYYMMDD format'
    )
    
    parser.add_argument(
        '--end-date', 
        required=False,
        help='End date in YYYYMMDD format'
    )
    
    parser.add_argument(
        '-d', '--date',
        help='Process a single date in YYYYMMDD format (alternative to --start-date/--end-date)'
    )
    
    parser.add_argument(
        '--latest', 
        action='store_true',
        help='Process files from the "latest" subdirectory within each date directory'
    )
    
    parser.add_argument(
        '--data-version', 
        default='1.0.1',
        help='Data version string'
    )
    
    parser.add_argument(
        '--noise-file',
        type=str,
        help='Path to noise CSV file for noise subtraction (optional)'
    )
    
    parser.add_argument(
        '--force', 
        action='store_true',
        help='Force processing even if output files exist'
    )
    
    parser.add_argument(
        '--dry-run', 
        action='store_true',
        help='Show what would be processed without actually processing'
    )
    
    parser.add_argument(
        '--skip-missing', 
        action='store_true',
        help='Skip dates with no input data instead of stopping'
    )

    parser.add_argument(
        '--single-sweep', 
        action='store_true',
        help='Create separate files for each sweep instead of multi-sweep files'
    )

    parser.add_argument(
        '--noise-only', 
        action='store_true',
        help='Only calculate noise levels, do not process data'
    )
    
    parser.add_argument(
        '--skip-noise-calc', 
        action='store_true',
        help='Skip noise calculation step and proceed directly to processing'
    )
    
    args = parser.parse_args()
    
    # Set up paths
    try:
        paths = setup_kasbex_paths(use_latest=args.latest)
        print(f"CAMRA KASBEX Campaign Processing")
        if args.latest:
            print(f"Mode: Processing from 'latest' subdirectories")
        print(f"Base input path: {paths['inpath']}")
        print(f"Output path: {paths['outpath']}")
        print(f"Project YAML: {paths['yaml_project_file']}")
        print(f"Instrument YAML: {paths['yaml_instrument_file']}")
        
        # Add noise file info
        if args.noise_file:
            print(f"Noise file: {args.noise_file}")
            if not Path(args.noise_file).exists():
                print(f"Warning: Noise file not found at {args.noise_file}")
        else:
            print("No noise file specified - noise subtraction will be skipped")
            
    except Exception as e:
        print(f"Error setting up paths: {e}")
        sys.exit(1)
    
    # Validate that required files exist
    for key, path in paths.items():
        if key.endswith('_file'):  # YAML files
            if not Path(path).exists():
                print(f"Warning: {key} not found at {path}")
                print("You may need to create this YAML file or update the path")
    
    # Generate date list
    try:
        if args.date:
            # Single date processing
            date_list = [args.date]
            # Validate the date format
            datetime.datetime.strptime(args.date, '%Y%m%d')
            print(f"Processing single date: {args.date}")
        elif args.start_date and args.end_date:
            # Date range processing
            date_list = generate_date_list(args.start_date, args.end_date)
            print(f"Processing {len(date_list)} dates from {args.start_date} to {args.end_date}")
        else:
            parser.error("Either --date (-d) or both --start-date and --end-date must be provided")
        
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Get campaign configuration
    try:
        campaign_info = get_campaign_info('kasbex')
        print(f"Campaign info: {campaign_info}")
    except Exception as e:
        print(f"Warning: Could not get campaign info for KASBEX: {e}")
        # Use default values for CAMRA
        campaign_info = {
            'tracking_tag': 'AMOF_CAMRA_KASBEX',
            'location': 'cao'
        }
    
    if args.dry_run:
        print("\nDRY RUN: Would process the following dates:")
        for datestr in date_list:
            has_data = check_input_data(paths['inpath'], datestr, use_latest=args.latest)
            status = "HAS DATA" if has_data else "NO DATA"
            input_path = get_input_path(paths['inpath'], datestr, use_latest=args.latest)
            print(f"  {datestr} - {status} (from: {input_path})")
        print(f"\nConfiguration:")
        print(f"  Data version: {args.data_version}")
        print(f"  Single sweep mode: {args.single_sweep}")
        print(f"  Use latest subdirectory: {args.latest}")
        print(f"  Noise file: {args.noise_file if args.noise_file else 'None'}")
        sys.exit(0)
    
    # Process each date
    total_processed = 0
    total_skipped = 0
    total_errors = 0
    
    overall_start_time = datetime.datetime.now()
    
    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing CAMRA KASBEX date {i}/{len(date_list)}: {datestr}")
        print(f"{'='*60}")
        
        # Check for input data
        if not check_input_data(paths['inpath'], datestr, use_latest=args.latest):
            print(f"No input data found for {datestr}")
            if args.skip_missing:
                print("Skipping (--skip-missing enabled)")
                total_skipped += 1
                continue
            else:
                print("Stopping (use --skip-missing to skip missing dates)")
                break
        
        # Get the appropriate input path for this date
        date_input_path = get_input_path(paths['inpath'], datestr, use_latest=args.latest)
        
        # Set up noise file path
        if args.noise_file:
            noise_file_path = args.noise_file
        else:
            # Default noise file in output directory
            output_date_dir = Path(paths['outpath']) / datestr
            output_date_dir.mkdir(parents=True, exist_ok=True)
            noise_file_path = str(output_date_dir / f'camra_noise_{datestr}.csv')
        
        try:
            # STEP 1: Calculate noise levels (unless skipped)
            if not args.skip_noise_calc:
                print(f"\n--- STEP 1: Calculating noise levels for {datestr} ---")
                from camra_utils import process_camra_day_noise_calculation
                
                calculated_noise_file = process_camra_day_noise_calculation(
                    datestr,
                    date_input_path,
                    noise_file_path
                )
                
                if calculated_noise_file:
                    print(f"✓ Noise calculation completed: {calculated_noise_file}")
                else:
                    print("⚠️  No noise data calculated")
                
                # If noise-only mode, skip to next date
                if args.noise_only:
                    print(f"Noise-only mode: Skipping data processing for {datestr}")
                    total_processed += 1
                    continue
            
            # STEP 2: Full data processing
            print(f"\n--- STEP 2: Processing radar data for {datestr} ---")
            start_time = datetime.datetime.now()
            
            # Process using the campaign processing module
            process_campaign_day(
                campaign='kasbex',
                datestr=datestr,
                inpath=date_input_path,
                outpath=paths['outpath'],
                yaml_project_file=str(paths['yaml_project_file']),
                yaml_instrument_file=str(paths['yaml_instrument_file']),
                data_version=args.data_version,
                single_sweep=args.single_sweep,
                noise_file=noise_file_path  # Use the calculated noise file
            )
            
            end_time = datetime.datetime.now()
            duration = end_time - start_time
            
            print(f"✓ Successfully completed CAMRA KASBEX processing for {datestr}")
            print(f"Processing time: {duration}")
            total_processed += 1
            
        except Exception as e:
            print(f"✗ Error during processing of {datestr}: {e}")
            import traceback
            traceback.print_exc()
            total_errors += 1
            
            if not args.skip_missing:
                print("Stopping due to error (use --skip-missing to continue)")
                break
    
    # Summary
    overall_end_time = datetime.datetime.now()
    total_duration = overall_end_time - overall_start_time
    
    print(f"\n{'='*60}")
    print("CAMRA KASBEX BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    if args.date:
        print(f"Single date: {args.date}")
    else:
        print(f"Date range: {args.start_date} to {args.end_date}")
    print(f"Total dates processed: {len(date_list)}")
    print(f"Successfully processed: {total_processed}")
    print(f"Skipped: {total_skipped}")
    print(f"Errors: {total_errors}")
    print(f"Latest subdirectory mode: {args.latest}")
    print(f"Total processing time: {total_duration}")
    print(f"Average per date: {total_duration / len(date_list) if date_list else 'N/A'}")
    print(f"Output directory: {paths['outpath']}")
    print(f"{'='*60}")
    
    # Exit with error code if there were any errors
    if total_errors > 0:
        print(f"⚠️  {total_errors} dates had processing errors")
        sys.exit(1)
    else:
        print("✓ All processing completed successfully!")

if __name__ == "__main__":
    main()