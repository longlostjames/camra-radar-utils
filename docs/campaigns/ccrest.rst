CCREST-M Campaign
=================

**Cloud, Convection and Rain ExperimentS in the Tropics – Maritime**

Overview
--------

CCREST-M is a tropical field campaign studying cloud, convection and rainfall
processes.  CAMRa data were collected at the Chilbolton Observatory in 2023.
Three complementary processing streams are used:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Mode
     - Description
   * - **RHI**
     - Range-Height Indicator scans at two fixed azimuths: 270° (CCREST1)
       and 246° (CCREST2).  Output scan names are ``RHI-CCREST1`` and
       ``RHI-CCREST2``.
   * - **VPT**
     - Vertically-Pointing mode.  Moment data converted directly to L1a.
   * - **TS-VPT**
     - Vertically-Pointing time series.  Full pipeline:
       L0a → L0b → L1 → moments (L1a / L1b).

Configuration
-------------

.. list-table::
   :widths: 30 70

   * - Tracking tag
     - ``AMOF_20230201132601``
   * - Default data version
     - ``1.0.0``
   * - Project YAML
     - ``campaigns/ccrest-m_project.yml``
   * - Instrument YAML
     - ``instrument_metadata.yml``

Processing scripts
------------------

Single-date scripts
~~~~~~~~~~~~~~~~~~~

``proc_camra2ncas_ccrest.py``
    Processes RHI scans (both azimuths) for a single date.

``proc_camra2ncas_ccrest_vpt.py``
    Processes VPT moment data for a single date.

``proc_camra2ncas_ccrest_ts_vpt.py``
    Processes VPT time-series data through the full L0a → moments pipeline
    for a single date.

Batch script
~~~~~~~~~~~~

``proc_camra_ccrest_campaign_batch.py``
    Processes one or more dates with a selectable ``--mode``.

Batch script options
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    -d / --date YYYYMMDD          Process a single date
    --start-date / --end-date     Process a date range
    --mode rhi|vpt|ts_vpt|all     Processing mode (default: all)
    --data-version VERSION        Default: 1.0.0
    --dry-run                     Preview without processing
    --skip-missing                Continue past dates with no input data

Usage
-----

.. code-block:: bash

    # All modes for a single date
    python proc_camra_ccrest_campaign_batch.py -d 20230601

    # RHI only over a date range
    python proc_camra_ccrest_campaign_batch.py \
        --start-date 20230601 --end-date 20230630 --mode rhi

    # VPT time series only
    python proc_camra_ccrest_campaign_batch.py -d 20230601 --mode ts_vpt

    # Dry run
    python proc_camra_ccrest_campaign_batch.py \
        --start-date 20230601 --end-date 20230630 --dry-run

Output
------

RHI and VPT moment files are written to::

    /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1/ccrest-m/L1a/<YYYYMMDD>/

VPT time series intermediate and final files::

    .../ccrest-m/ts/L0b/<YYYYMMDD>/
    .../ccrest-m/ts/L1/<YYYYMMDD>/
    .../ccrest-m/L1b/<YYYYMMDD>/

API reference
-------------

.. seealso::

   :func:`campaign_processing.process_camra_ccrest_day_step1`,
   :func:`campaign_processing.process_camra_ccrest_vpt_day_step1`,
   :func:`campaign_processing.process_camra_ccrest_vpt_day_ts`
