CDMS Analysis Notebooks and Tools
=================================

GitHub: kelseymh/CDMS_Analysis

This directory contains my personal analysis notebooks and supporting
Python tools.  It is used exclusively on the Texas A&M Grace cluster,
as /scratch/user/kelsey/jupyter/.


Data Directories
----------------

For the notebooks to be functional, a couple of symbolic links need to
be set, pointing to top-level data areas.

    CDMSdata -> /scratch/group/mitchcomp/CDMS/data
    data -> /scratch/user/kelsey/data


Jupyter Notebooks
-----------------

The repository contaings all of my Jupyter notebooks.  They are not
well documented (yet?), although I try to include comments and
docstrings to remind myself how to use them.  If this repository ends
up being advertised, I will try to provide guidaces on notebooks which
might be more generally useful.


Jupyter Utilities
-----------------

imports.py
: The Python file `imports.py` is used in many of my notebooks to set
up a consistent environment with pre-loaded functionality.

phonon_eff.py
: Utility functions for computing and plotting phonon collection
efficiency from DMC ROOT files.  Use `dir(phonon_eff)` after import
for a list of available functions.

traces_rdf.py
: Collection of functions for loading and plotting TES and FET traces
from DMC ROOT files.  Uses ROOT's RDataFrame directly, with very poor
memory performance.  Superseded by `traceReader.py`.

traces_numpy.py
: This is the original root_numpy version of traces_rdf.py.  It should
not be used, but is kept in the repository for reference, and because
it is used in fit-traces_snolab.ipynb.


Python Classes
--------------

ChannelShape.py
: This defines a class which emulates the `CDMSChannelShape` class in
SuperSim.  See the docstrings for details.

traceReader.py
: Driver class for loading and plotting TES and FET traces from DMC ROOT
files.  Uses CATs' CDataFrame utility.  See the docstrings for details.


Python Programs
---------------

trace_fitter.py
: This is an executable program which processes a DMC ROOT file and
can produce simple fits for TES and FET traces. Use `./trace_fitter.py -h`
for more information.
