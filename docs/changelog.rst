Changelog
=========

v1.3.0 (2026-02-22)
--------------------
* Added DYMECS campaign processing (``process_camra_dymecs_day_step1``).
* Added ``proc_camra_dymecs_campaign_batch.py``.
* Added ``instrument_metadata.yml`` — CAMRa-only extract of instrument metadata,
  making the repository self-contained.
* Created ``campaigns/`` directory; moved all ``*_project.yml`` files into it.
* Updated all ``proc_camra2ncas_*.py`` scripts to resolve YAML files from the
  repository ``campaigns/`` directory rather than ``~/amof_campaigns/``.
* Fixed ``proc_camra_kasbex_campaign_batch.py``: use repository YAML paths;
  corrected output level ``L1c`` → ``L1``.
* Added ``proc_camra_ccrest_campaign_batch.py`` with ``--mode`` flag
  (``rhi`` | ``vpt`` | ``ts_vpt`` | ``all``).
* Added Sphinx documentation.

v1.2.0
------
* Added CCREST-M processing capability.
* Added ``process_camra_ccrest_day_step1``, ``process_camra_ccrest_vpt_day_step1``,
  ``process_camra_ccrest_vpt_day_ts`` to ``campaign_processing.py``.
* Added ``proc_camra2ncas_ccrest.py``, ``proc_camra2ncas_ccrest_vpt.py``,
  ``proc_camra2ncas_ccrest_ts_vpt.py``.

v1.1.0
------
* Refactored all processing scripts to use the unified ``campaign_processing``
  module.
* Added ``process_camra_woest_day_step1`` with auto-detection of IOP/SOP/Other
  modes.
* Added ``process_camra_kasbex_day_step1``.
* Added ``proc_camra_kasbex_campaign_batch.py``.

v1.0.0
------
* Initial release — WOEST CAMRa processing scripts.
