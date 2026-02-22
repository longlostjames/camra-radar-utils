WOEST Campaign
==============

**WesCon - Observing the Evolving Structures of Turbulence (WOEST)**

Overview
--------

WOEST was a CAMRa campaign focused on observing the evolving structures of
turbulence in convective and stratiform precipitation.  Data were collected at
the Chilbolton Observatory (CAO) during summer 2023.

Three observation modes were used during the campaign:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Mode
     - Description
   * - **IOP** (Intensive Observation Period)
     - Targeted storm-tracking RHI scans, organised into regional
       sub-directories (``region_01``, ``region_02``, …).
   * - **SOP** (Standard Operating Period)
     - Enhanced scanning strategy, files identified by ``sop`` in filename
       or stored in a dedicated ``sop/`` sub-directory.
   * - **Other**
     - Routine operational scanning outside IOP and SOP periods.

Configuration
-------------

.. list-table::
   :widths: 30 70

   * - Tracking tag
     - ``AMOF_20220922221548``
   * - Default data version
     - ``1.0.1``
   * - Project YAML
     - ``campaigns/woest_project.yml``
   * - Instrument YAML
     - ``instrument_metadata.yml``

Processing scripts
------------------

Single-date scripts
~~~~~~~~~~~~~~~~~~~

``proc_camra2ncas_woest.py``
    Auto-detects the observation mode (IOP / SOP / Other) from the
    input directory structure.  Accepts ``-d YYYYMMDD``.

``proc_camra2ncas_woest_iop.py``
    Processes IOP (storm-tracking) data explicitly.

``proc_camra2ncas_woest_sop.py``
    Processes SOP data explicitly.

``proc_camra2ncas_woest_other.py``
    Processes routine (Other) data explicitly.

Usage
-----

.. code-block:: bash

    # Auto-detect mode
    python proc_camra2ncas_woest.py -d 20230818

    # Force IOP mode
    python proc_camra2ncas_woest_iop.py -d 20230818

    # Force SOP mode
    python proc_camra2ncas_woest_sop.py -d 20230818

Output
------

Processed files are written to::

    /gws/pw/j07/ncas_obs_vol2/cao/processing/ncas-radar-camra-1/woest/L1/<YYYYMMDD>/

File naming follows the NCAS standard::

    ncas-radar-camra-1_cao_<YYYYMMDD>-<HHMMSS>_<scan>_l1_v<version>.nc

API reference
-------------

.. seealso::

   :func:`campaign_processing.process_camra_woest_day_step1`,
   :func:`campaign_processing.process_camra_woest_iop_step1`,
   :func:`campaign_processing.process_camra_woest_sop_step1`,
   :func:`campaign_processing.process_camra_woest_other_step1`
