#!/usr/bin/env python3
"""
proc_camra_dymecs_campaign_batch.py

Batch process DYMECS campaign data over a date range using campaign_processing
module for CAMRA radar.

Usage:
    python proc_camra_dymecs_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD
    python proc_camra_dymecs_campaign_batch.py -d YYYYMMDD
    python proc_camra_dymecs_campaign_batch.py -d YYYYMMDD --single-sweep

Author: Chris Walden, UK Research & Innovation and
        National Centre for Atmospheric Science
"""

import sys
import os
import argparse
import datetime
import traceback
from pathlib import Path

# Add the script directory to Python path
script_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(script_dir))

from campaign_processing import process_camra_dymecs_day_step1, get_campaign_info

# ---------------------------------------------------------------------------
# Campaign constants
# ---------------------------------------------------------------------------
CAMPAIGN     = 'dymecs'
TRACKING_TAG = 'CRF_85'
CAMRA_BASE   = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-camra-1'
PROC_BASE    = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1'

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def setup_dymecs_paths():
    """Return a dict of all relevant paths for DYMECS processing."""
    paths = {
        'inpath'               : os.path.join(CAMRA_BASE, 'data', 'campaign', CAMPAIGN, 'raw'),
        'outpath'              : os.path.join(PROC_BASE, CAMPAIGN, 'L1'),
        'yaml_project_file'    : str(script_dir / 'campaigns' / f'{CAMPAIGN}_project.yml'),
        'yaml_instrument_file' : str(script_dir / 'instrument_metadata.yml'),
    }

    Path(paths['outpath']).mkdir(parents=True, exist_ok=True)
    return paths


def generate_date_list(start_date_str, end_date_str):
    """Generate a list of YYYYMMDD strings between start and end (inclusive)."""
    try:
        start = datetime.datetime.strptime(start_date_str, '%Y%m%d')
        end   = datetime.datetime.strptime(end_date_str,   '%Y%m%d')
    except ValueError as e:
        raise ValueError(f"Invalid date format. Use YYYYMMDD. Error: {e}")

    if start > end:
        raise ValueError("Start date must be before or equal to end date")

    dates = []
    current = start
    while current <= end:
        dates.append(current.strftime('%Y%m%d'))
        current += datetime.timedelta(days=1)
    return dates


def check_input_data(inpath, datestr):
    """Return True if any .nc files exist for datestr in the input path."""
    date_path = Path(inpath) / datestr
    if not date_path.exists():
        print(f"  Directory does not exist: {date_path}")
        return False
    nc_files = list(date_path.glob('*.nc'))
    print(f"  Found {len(nc_files)} nc files in {date_path}")
    if nc_files:
        print(f"  Sample: {[f.name for f in nc_files[:3]]}")
    return len(nc_files) > 0

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch process CAMRA DYMECS campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--start-date', help='Start date in YYYYMMDD format')
    parser.add_argument('--end-date',   help='End date in YYYYMMDD format')
    parser.add_argument('-d', '--date', help='Single date in YYYYMMDD format')

    parser.add_argument('--data-version',  default='1.0.0', help='Data version string')
    parser.add_argument('--single-sweep',  action='store_true',
                        help='Create a separate output file for each sweep (default: True for DYMECS)')
    parser.add_argument('--no-single-sweep', action='store_true',
                        help='Concatenate all sweeps of each scan type into one file')
    parser.add_argument('--dry-run',       action='store_true',
                        help='Show what would be processed without processing')
    parser.add_argument('--skip-missing',  action='store_true',
                        help='Skip dates with no input data instead of stopping')

    args = parser.parse_args()

    # DYMECS default is single_sweep=True; --no-single-sweep overrides
    single_sweep = not args.no_single_sweep

    # --- Set up paths -------------------------------------------------------
    try:
        paths = setup_dymecs_paths()
    except Exception as e:
        print(f"Error setting up paths: {e}")
        sys.exit(1)

    print("CAMRA DYMECS Campaign Processing")
    print(f"Input path     : {paths['inpath']}")
    print(f"Output path    : {paths['outpath']}")
    print(f"Project YAML   : {paths['yaml_project_file']}")
    print(f"Instrument YAML: {paths['yaml_instrument_file']}")
    print(f"Single sweep   : {single_sweep}")

    for key in ('yaml_project_file', 'yaml_instrument_file'):
        if not Path(paths[key]).exists():
            print(f"Warning: {key} not found at {paths[key]}")

    # --- Build date list ----------------------------------------------------
    try:
        if args.date:
            datetime.datetime.strptime(args.date, '%Y%m%d')   # validate
            date_list = [args.date]
            print(f"Processing single date: {args.date}")
        elif args.start_date and args.end_date:
            date_list = generate_date_list(args.start_date, args.end_date)
            print(f"Processing {len(date_list)} dates: {args.start_date} to {args.end_date}")
        else:
            parser.error("Either --date (-d) or both --start-date and --end-date must be provided")
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # --- Dry run ------------------------------------------------------------
    if args.dry_run:
        print("\nDRY RUN - would process the following dates:")
        for datestr in date_list:
            has_data = check_input_data(paths['inpath'], datestr)
            print(f"  {datestr}  {'HAS DATA' if has_data else 'NO DATA'}")
        print(f"\nData version : {args.data_version}")
        print(f"Single sweep : {single_sweep}")
        sys.exit(0)

    # --- Processing loop ----------------------------------------------------
    total_processed = 0
    total_skipped   = 0
    total_errors    = 0
    overall_start   = datetime.datetime.now()

    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing CAMRA DYMECS date {i}/{len(date_list)}: {datestr}")
        print(f"{'='*60}")

        if not check_input_data(paths['inpath'], datestr):
            print(f"No input data found for {datestr}")
            if args.skip_missing:
                print("Skipping (--skip-missing enabled)")
                total_skipped += 1
                continue
            else:
                print("Stopping. Use --skip-missing to skip missing dates.")
                break

        date_inpath = os.path.join(paths['inpath'], datestr)
        start_time  = datetime.datetime.now()

        try:
            process_camra_dymecs_day_step1(
                datestr,
                date_inpath,
                paths['outpath'],
                paths['yaml_project_file'],
                paths['yaml_instrument_file'],
                data_version=args.data_version,
                tracking_tag=TRACKING_TAG,
                campaign=CAMPAIGN,
                single_sweep=single_sweep,
            )

            duration = datetime.datetime.now() - start_time
            print(f"✓ Successfully completed DYMECS processing for {datestr} ({duration})")
            total_processed += 1

        except Exception as e:
            duration = datetime.datetime.now() - start_time
            print(f"✗ Error processing {datestr} ({duration}): {e}")
            traceback.print_exc()
            total_errors += 1

            if not args.skip_missing:
                print("Stopping due to error. Use --skip-missing to continue.")
                break

    # --- Summary ------------------------------------------------------------
    total_duration = datetime.datetime.now() - overall_start

    print(f"\n{'='*60}")
    print("CAMRA DYMECS BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    if args.date:
        print(f"Date           : {args.date}")
    else:
        print(f"Date range     : {args.start_date} to {args.end_date}")
    print(f"Total dates    : {len(date_list)}")
    print(f"Processed OK   : {total_processed}")
    print(f"Skipped        : {total_skipped}")
    print(f"Errors         : {total_errors}")
    print(f"Single sweep   : {single_sweep}")
    print(f"Total time     : {total_duration}")
    print(f"Output path    : {paths['outpath']}")
    print(f"{'='*60}")

    if total_errors > 0:
        print(f"⚠  {total_errors} date(s) had processing errors")
        sys.exit(1)
    else:
        print("✓ All processing completed successfully!")


if __name__ == "__main__":
    main()
