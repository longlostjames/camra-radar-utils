campaign\_processing module
===========================

This is the core module of **CAMRa Radar Utils**.  It provides
campaign-specific processing functions for all supported campaigns and a set
of shared helper utilities.

.. contents:: Contents
   :local:
   :depth: 2

Module overview
---------------

.. automodule:: campaign_processing
   :no-members:

WOEST functions
---------------

.. autofunction:: campaign_processing.process_camra_woest_day_step1
.. autofunction:: campaign_processing.process_camra_woest_iop_step1
.. autofunction:: campaign_processing.process_camra_woest_sop_step1
.. autofunction:: campaign_processing.process_camra_woest_other_step1

KASBEX functions
----------------

.. autofunction:: campaign_processing.process_camra_kasbex_day_step1

CCREST-M functions
------------------

.. autofunction:: campaign_processing.process_camra_ccrest_day_step1
.. autofunction:: campaign_processing.process_camra_ccrest_vpt_day_step1
.. autofunction:: campaign_processing.process_camra_ccrest_vpt_day_ts

DYMECS functions
----------------

.. autofunction:: campaign_processing.process_camra_dymecs_day_step1

Helper functions
----------------

.. autofunction:: campaign_processing.find_camra_files
.. autofunction:: campaign_processing.process_camra_scan_type
.. autofunction:: campaign_processing.process_campaign_day
.. autofunction:: campaign_processing.get_campaign_info
.. autofunction:: campaign_processing.load_yaml_config
