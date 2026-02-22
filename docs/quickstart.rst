Quick Start
===========

This page shows how to run the most common processing tasks for each campaign.

Single-date processing
----------------------

All single-date scripts follow the same interface::

    python proc_camra2ncas_<campaign>.py -d YYYYMMDD

For example, to process KASBEX data for 1 June 2024::

    python proc_camra2ncas_kasbex.py -d 20240601

Batch processing over a date range
-----------------------------------

Each campaign has a dedicated batch script that accepts ``--start-date`` and
``--end-date`` (or ``-d`` for a single date).

KASBEX
~~~~~~
::

    # Single date
    python proc_camra_kasbex_campaign_batch.py -d 20240601

    # Date range, creating one file per sweep, skipping missing days
    python proc_camra_kasbex_campaign_batch.py \
        --start-date 20240601 --end-date 20240630 \
        --single-sweep --skip-missing

    # Preview without processing
    python proc_camra_kasbex_campaign_batch.py \
        --start-date 20240601 --end-date 20240630 --dry-run

CCREST-M
~~~~~~~~
::

    # All modes (RHI + VPT + VPT time series) for a single date
    python proc_camra_ccrest_campaign_batch.py -d 20230601

    # RHI only over a date range
    python proc_camra_ccrest_campaign_batch.py \
        --start-date 20230601 --end-date 20230630 --mode rhi

    # VPT time series only
    python proc_camra_ccrest_campaign_batch.py -d 20230601 --mode ts_vpt

DYMECS
~~~~~~
::

    # Single date (single-sweep mode by default)
    python proc_camra_dymecs_campaign_batch.py -d 20120820

    # Date range, concatenating sweeps
    python proc_camra_dymecs_campaign_batch.py \
        --start-date 20120801 --end-date 20120831 --no-single-sweep

WOEST
~~~~~
::

    # Auto-detect IOP/SOP/Other from directory structure
    python proc_camra2ncas_woest.py -d 20230818

    # Force IOP mode explicitly
    python proc_camra2ncas_woest_iop.py -d 20230818

Data version and output paths
------------------------------

Default output paths are set inside each batch script (pointing to the
``/gws/...`` project workspace on JASMIN).  To override, edit
``setup_<campaign>_paths()`` in the relevant batch script, or pass
``-o <outpath>`` to the single-date scripts.

The data version defaults to ``"1.0.0"`` and can be overridden in the batch
scripts with ``--data-version 1.0.1``.
