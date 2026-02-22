Installation
============

Requirements
------------

* Python 3.9+
* `PyART <https://arm-doe.github.io/pyart/>`_
* `netCDF4 <https://unidata.github.io/netcdf4-python/>`_
* `NumPy <https://numpy.org/>`_
* `PyYAML <https://pyyaml.org/>`_

These are runtime dependencies already present in most atmospheric science
Python environments (e.g. the JASMIN Sci Analysis Server).

Getting the code
----------------

Clone the repository from GitHub::

    git clone https://github.com/longlostjames/camra-radar-utils.git
    cd camra-radar-utils

No package installation is required — the scripts add the repository
directory to ``sys.path`` at runtime.  If you want to import
``campaign_processing`` from outside the repository, add the clone path to
your ``PYTHONPATH``::

    export PYTHONPATH=/path/to/camra-radar-utils:$PYTHONPATH

Configuration files
-------------------

Two YAML configuration files are needed at runtime:

1. **Campaign project file** — one per campaign, located in ``campaigns/``
   within the repository (e.g. ``campaigns/woest_project.yml``).

2. **Instrument metadata file** — ``instrument_metadata.yml`` in the
   repository root, containing CAMRa-specific hardware metadata.

These files are included in the repository and are discovered automatically
by the processing scripts — no path configuration is needed.

Building the documentation
--------------------------

Documentation dependencies::

    pip install sphinx sphinx-rtd-theme

Build HTML locally::

    cd docs
    make html

The built site is written to ``docs/_build/html/``.
