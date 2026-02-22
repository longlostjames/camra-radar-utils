DYMECS Campaign
===============

**DYnamics, Microphysics and Entrainment in Convective Systems**

Overview
--------

DYMECS was a campaign at the Chilbolton Observatory studying the dynamics and
microphysics of convective systems.  CAMRa collected RHI and PPI scans, each
stored as a single-sweep NetCDF file by default.

Configuration
-------------

.. list-table::
   :widths: 30 70

   * - Tracking tag
     - ``CRF_85``
   * - Default data version
     - ``1.0.0``
   * - Project YAML
     - ``campaigns/dymecs_project.yml``
   * - Instrument YAML
     - ``instrument_metadata.yml``

Processing scripts
------------------

``proc_camra2ncas_dymecs.py``
    Single-date processing.  Accepts ``-d YYYYMMDD``.

``proc_camra_dymecs_campaign_batch.py``
    Batch processing over a date range.  Recommended for production use.

Batch script options
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    -d / --date YYYYMMDD          Process a single date
    --start-date / --end-date     Process a date range
    --data-version VERSION        Default: 1.0.0
    --single-sweep                One output file per sweep (default behaviour)
    --no-single-sweep             Concatenate all sweeps per scan type
    --dry-run                     Preview without processing
    --skip-missing                Continue past dates with no input data

Usage
-----

.. code-block:: bash

    # Single date (default: single-sweep output)
    python proc_camra_dymecs_campaign_batch.py -d 20120820

    # Date range, skip missing days
    python proc_camra_dymecs_campaign_batch.py \
        --start-date 20120801 --end-date 20120831 --skip-missing

    # Concatenate sweeps into one file per scan type
    python proc_camra_dymecs_campaign_batch.py \
        -d 20120820 --no-single-sweep

    # Dry run
    python proc_camra_dymecs_campaign_batch.py \
        --start-date 20120801 --end-date 20120831 --dry-run

Output
------

Processed files are written to::

    /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1/dymecs/L1/<YYYYMMDD>/

API reference
-------------

.. seealso::

   :func:`campaign_processing.process_camra_dymecs_day_step1`
