KASBEX Campaign
===============

**Ka- and S-band Experiment (KASBEX)**

Overview
--------

KASBEX was a dual-wavelength radar campaign using the CAMRa S-band radar
alongside the NCAS Ka-band radar (Kepler).  The campaign focused on
characterising radar reflectivity differences between the two wavelengths
to study precipitation microphysics.

A noise-calibration step precedes data conversion: for each day, noise
levels are estimated from the raw files and stored in a per-day CSV file,
which is then passed to the processing function for noise subtraction.

Configuration
-------------

.. list-table::
   :widths: 30 70

   * - Tracking tag
     - ``AMOF_20250508133639``
   * - Default data version
     - ``1.0.1``
   * - Project YAML
     - ``campaigns/kasbex_project.yml``
   * - Instrument YAML
     - ``instrument_metadata.yml``

Processing scripts
------------------

``proc_camra2ncas_kasbex.py``
    Single-date processing.  Accepts ``-d YYYYMMDD``.

``proc_camra_kasbex_campaign_batch.py``
    Batch processing over a date range with full noise-subtraction pipeline.
    Recommended for production use.

Batch script options
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    -d / --date YYYYMMDD          Process a single date
    --start-date / --end-date     Process a date range
    --data-version VERSION        Default: 1.0.1
    --noise-file PATH             Path to a pre-computed noise CSV
    --single-sweep                One output file per sweep
    --latest                      Read from <date>/latest/ sub-directory
    --dry-run                     Preview without processing
    --skip-missing                Continue past dates with no input data
    --noise-only                  Calculate noise levels only
    --skip-noise-calc             Skip noise calculation, process directly

Usage
-----

.. code-block:: bash

    # Single date, full pipeline
    python proc_camra_kasbex_campaign_batch.py -d 20240601

    # Date range, single-sweep files, skip missing days
    python proc_camra_kasbex_campaign_batch.py \
        --start-date 20240601 --end-date 20240630 \
        --single-sweep --skip-missing

    # Dry run to check data availability
    python proc_camra_kasbex_campaign_batch.py \
        --start-date 20240601 --end-date 20240630 --dry-run

Output
------

Processed files are written to::

    /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1/kasbex/L1/<YYYYMMDD>/

API reference
-------------

.. seealso::

   :func:`campaign_processing.process_camra_kasbex_day_step1`
