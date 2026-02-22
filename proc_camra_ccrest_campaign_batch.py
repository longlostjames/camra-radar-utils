#!/usr/bin/env python3
"""
proc_camra_ccrest_campaign_batch.py

Batch process CCREST-M campaign data over a date range using campaign_processing
module for CAMRA radar.

Three processing modes are available:
  rhi      - RHI scans at 270° and 246° azimuths -> L1a
  vpt      - Vertically-pointing moment data -> L1a
  ts_vpt   - Vertically-pointing time series: L0a -> L0b -> L1 -> moments
  all      - Run all three modes in sequence

Usage:
    python proc_camra_ccrest_campaign_batch.py --start-date YYYYMMDD --end-date YYYYMMDD [--mode MODE]
    python proc_camra_ccrest_campaign_batch.py -d YYYYMMDD [--mode MODE]

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

import campaign_processing
from campaign_processing import get_campaign_info

# ---------------------------------------------------------------------------
# Campaign constants
# ---------------------------------------------------------------------------
CAMPAIGN         = 'ccrest-m'
TRACKING_TAG     = 'AMOF_20230201132601'
CAMRA_BASE       = '/gws/pw/j07/ncas_obs_vol2/cao/raw_data/ncas-radar-camra-1'
PROC_BASE        = '/gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1'

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def setup_ccrest_paths():
    """Return a dict of all relevant paths for CCREST-M processing."""
    raw_base  = os.path.join(CAMRA_BASE, 'data', 'campaign', CAMPAIGN, 'raw')
    ts_base   = os.path.join(CAMRA_BASE, 'data', 'campaign', CAMPAIGN, 'ts')
    proc_camp = os.path.join(PROC_BASE, CAMPAIGN)

    paths = {
        'raw_inpath'           : raw_base,
        'ts_inpath'            : ts_base,
        'outpath_l1a'          : os.path.join(proc_camp, 'L1a'),
        'outpath_ts'           : proc_camp,          # ts sub-structure built inside function
        'yaml_project_file'    : str(script_dir / 'campaigns' / f'{CAMPAIGN}_project.yml'),
        'yaml_instrument_file' : str(script_dir / 'instrument_metadata.yml'),
    }

    # Ensure top-level output directories exist
    Path(paths['outpath_l1a']).mkdir(parents=True, exist_ok=True)

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


def check_raw_data(raw_inpath, datestr):
    """Return True if any .nc files exist for datestr in the raw input path."""
    date_path = Path(raw_inpath) / datestr
    if not date_path.exists():
        print(f"  Raw directory does not exist: {date_path}")
        return False
    nc_files = list(date_path.glob('*.nc'))
    print(f"  Found {len(nc_files)} nc files in {date_path}")
    return len(nc_files) > 0


def check_ts_data(ts_inpath, datestr):
    """Return True if any time-series files exist for datestr."""
    date_path = Path(ts_inpath) / datestr
    if not date_path.exists():
        print(f"  TS directory does not exist: {date_path}")
        return False
    files = list(date_path.glob('*'))
    print(f"  Found {len(files)} files in {date_path}")
    return len(files) > 0


# ---------------------------------------------------------------------------
# Per-date processing functions
# ---------------------------------------------------------------------------

def process_rhi(datestr, paths, data_version):
    """Process RHI scans (270° and 246° azimuths) -> L1a."""
    inpath  = os.path.join(paths['raw_inpath'], datestr)
    outpath = paths['outpath_l1a']

    print(f"  RHI: {inpath} -> {outpath}")
    campaign_processing.process_camra_ccrest_day_step1(
        datestr,
        inpath,
        outpath,
        paths['yaml_project_file'],
        paths['yaml_instrument_file'],
        data_version=data_version,
        tracking_tag=TRACKING_TAG,
        campaign=CAMPAIGN,
    )


def process_vpt(datestr, paths, data_version):
    """Process VPT moment files -> L1a."""
    inpath  = os.path.join(paths['raw_inpath'], datestr)
    outpath = paths['outpath_l1a']

    print(f"  VPT: {inpath} -> {outpath}")
    campaign_processing.process_camra_ccrest_vpt_day_step1(
        datestr,
        inpath,
        outpath,
        paths['yaml_project_file'],
        paths['yaml_instrument_file'],
        data_version=data_version,
        tracking_tag=TRACKING_TAG,
        campaign=CAMPAIGN,
    )


def process_ts_vpt(datestr, paths, data_version):
    """Process VPT time series: L0a -> L0b -> L1 -> moments."""
    ts_indir  = os.path.join(paths['ts_inpath'], datestr)
    mom_indir = paths['outpath_l1a']
    outdir    = paths['outpath_ts']

    # Ensure intermediate directories exist
    Path(os.path.join(outdir, 'ts', 'L0b', datestr)).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(outdir, 'ts', 'L1',  datestr)).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(outdir, 'L1b', datestr)).mkdir(parents=True, exist_ok=True)

    print(f"  TS-VPT: {ts_indir} -> {outdir}")
    campaign_processing.process_camra_ccrest_vpt_day_ts(
        datestr,
        ts_indir,
        mom_indir,
        outdir,
        paths['yaml_project_file'],
        paths['yaml_instrument_file'],
        data_version=data_version,
        tracking_tag=TRACKING_TAG,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch process CAMRA CCREST-M campaign radar data over a date range",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument('--start-date', help='Start date in YYYYMMDD format')
    parser.add_argument('--end-date',   help='End date in YYYYMMDD format')
    parser.add_argument('-d', '--date', help='Single date in YYYYMMDD format')

    parser.add_argument(
        '--mode',
        choices=['rhi', 'vpt', 'ts_vpt', 'all'],
        default='all',
        help='Processing mode: rhi (RHI scans), vpt (VPT moments), ts_vpt (VPT time series), all',
    )
    parser.add_argument('--data-version', default='1.0.0', help='Data version string')
    parser.add_argument('--dry-run',      action='store_true', help='Show what would be processed without processing')
    parser.add_argument('--skip-missing', action='store_true', help='Skip dates with no input data instead of stopping')

    args = parser.parse_args()

    # --- Set up paths -------------------------------------------------------
    try:
        paths = setup_ccrest_paths()
    except Exception as e:
        print(f"Error setting up paths: {e}")
        sys.exit(1)

    print("CAMRA CCREST-M Campaign Processing")
    print(f"Mode           : {args.mode}")
    print(f"Raw input path : {paths['raw_inpath']}")
    print(f"TS input path  : {paths['ts_inpath']}")
    print(f"Output (L1a)   : {paths['outpath_l1a']}")
    print(f"Output (ts)    : {paths['outpath_ts']}")
    print(f"Project YAML   : {paths['yaml_project_file']}")
    print(f"Instrument YAML: {paths['yaml_instrument_file']}")

    # Warn if YAML files are missing
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
            raw_ok = check_raw_data(paths['raw_inpath'], datestr)
            ts_ok  = check_ts_data(paths['ts_inpath'],   datestr)
            print(f"  {datestr}  raw={'OK' if raw_ok else 'MISSING'}  ts={'OK' if ts_ok else 'MISSING'}")
        print(f"\nData version: {args.data_version}")
        sys.exit(0)

    # --- Processing loop ----------------------------------------------------
    modes_to_run = ['rhi', 'vpt', 'ts_vpt'] if args.mode == 'all' else [args.mode]

    total_processed = 0
    total_skipped   = 0
    total_errors    = 0
    overall_start   = datetime.datetime.now()

    for i, datestr in enumerate(date_list, 1):
        print(f"\n{'='*60}")
        print(f"Date {i}/{len(date_list)}: {datestr}  (mode: {args.mode})")
        print(f"{'='*60}")

        # Determine which data sources are needed
        need_raw = any(m in modes_to_run for m in ('rhi', 'vpt'))
        need_ts  = 'ts_vpt' in modes_to_run

        raw_ok = check_raw_data(paths['raw_inpath'], datestr) if need_raw else True
        ts_ok  = check_ts_data( paths['ts_inpath'],  datestr) if need_ts  else True

        if (need_raw and not raw_ok) or (need_ts and not ts_ok):
            print(f"Required input data missing for {datestr}")
            if args.skip_missing:
                print("Skipping (--skip-missing enabled)")
                total_skipped += 1
                continue
            else:
                print("Stopping. Use --skip-missing to skip missing dates.")
                break

        date_start = datetime.datetime.now()
        date_errors = 0

        for mode in modes_to_run:
            print(f"\n--- {mode.upper()} ---")
            try:
                if mode == 'rhi':
                    process_rhi(datestr, paths, args.data_version)
                elif mode == 'vpt':
                    process_vpt(datestr, paths, args.data_version)
                elif mode == 'ts_vpt':
                    process_ts_vpt(datestr, paths, args.data_version)
                print(f"  {mode.upper()} completed successfully")
            except Exception as e:
                print(f"  Error in {mode.upper()} for {datestr}: {e}")
                traceback.print_exc()
                date_errors += 1

        duration = datetime.datetime.now() - date_start

        if date_errors == 0:
            print(f"\n✓ {datestr} completed in {duration}")
            total_processed += 1
        else:
            print(f"\n✗ {datestr} finished with {date_errors} error(s) in {duration}")
            total_errors += 1
            if not args.skip_missing:
                print("Stopping due to error. Use --skip-missing to continue.")
                break

    # --- Summary ------------------------------------------------------------
    total_duration = datetime.datetime.now() - overall_start

    print(f"\n{'='*60}")
    print("CAMRA CCREST-M BATCH PROCESSING SUMMARY")
    print(f"{'='*60}")
    if args.date:
        print(f"Date           : {args.date}")
    else:
        print(f"Date range     : {args.start_date} to {args.end_date}")
    print(f"Mode           : {args.mode}")
    print(f"Total dates    : {len(date_list)}")
    print(f"Processed OK   : {total_processed}")
    print(f"Skipped        : {total_skipped}")
    print(f"Errors         : {total_errors}")
    print(f"Total time     : {total_duration}")
    print(f"Output (L1a)   : {paths['outpath_l1a']}")
    print(f"{'='*60}")

    if total_errors > 0:
        print(f"⚠  {total_errors} date(s) had processing errors")
        sys.exit(1)
    else:
        print("✓ All processing completed successfully!")


if __name__ == "__main__":
    main()
