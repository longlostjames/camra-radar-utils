.. CAMRa Radar Utils documentation master file

CAMRa Radar Utils
=================

**CAMRa Radar Utils** is a Python toolkit for processing raw CAMRa
(Chilbolton Advanced Meteorological Radar) data into CF/RADIAL NetCDF files
compliant with the NCAS Observation Data standard.

The library currently supports four field campaigns:

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Campaign
     - Full name
     - Tracking tag
   * - WOEST
     - WESCON Observing the Evolving Structures of Turbulence
     - AMOF_20220922221548
   * - KASBEX
     - Ka- and S-Band EXperiment
     - AMOF_20250508133639
   * - CCREST-M
     - Cloud, Convection and Rain ExperimentS in the Tropics – Maritime
     - AMOF_20230201132601
   * - DYMECS
     - DYnamics, Microphysics and Entrainment in Convective Systems
     - CRF_85

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Campaigns

   campaigns/woest
   campaigns/kasbex
   campaigns/ccrest
   campaigns/dymecs

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/campaign_processing
   api/batch_scripts

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
