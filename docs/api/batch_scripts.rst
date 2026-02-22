Batch Processing Scripts
========================

These scripts are the recommended way to process a full campaign or a range
of dates.  They wrap the functions in :mod:`campaign_processing` and add
command-line argument handling, progress reporting, and error recovery.

.. contents:: Contents
   :local:
   :depth: 1

proc\_camra\_kasbex\_campaign\_batch
-------------------------------------

.. automodule:: proc_camra_kasbex_campaign_batch
   :members:
   :undoc-members:

proc\_camra\_ccrest\_campaign\_batch
--------------------------------------

.. automodule:: proc_camra_ccrest_campaign_batch
   :members:
   :undoc-members:

proc\_camra\_dymecs\_campaign\_batch
--------------------------------------

.. automodule:: proc_camra_dymecs_campaign_batch
   :members:
   :undoc-members:

Single-date scripts
-------------------

These simpler scripts process one date at a time using ``getopt``-style
``-d``/``-i``/``-o`` flags.  They are useful for quick one-off runs or for
submission as individual SLURM array jobs.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Script
     - Campaign / mode
   * - ``proc_camra2ncas_woest.py``
     - WOEST (auto-detects IOP/SOP/Other)
   * - ``proc_camra2ncas_woest_iop.py``
     - WOEST IOP
   * - ``proc_camra2ncas_woest_sop.py``
     - WOEST SOP
   * - ``proc_camra2ncas_woest_other.py``
     - WOEST Other
   * - ``proc_camra2ncas_kasbex.py``
     - KASBEX
   * - ``proc_camra2ncas_ccrest.py``
     - CCREST-M RHI
   * - ``proc_camra2ncas_ccrest_vpt.py``
     - CCREST-M VPT moments
   * - ``proc_camra2ncas_ccrest_ts_vpt.py``
     - CCREST-M VPT time series
   * - ``proc_camra2ncas_dymecs.py``
     - DYMECS
