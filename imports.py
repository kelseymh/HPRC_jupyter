## imnports.py -- does all the import statements for my Jupyter analysis
## Usage: from imports import *
##
## 20241202  Add environment variable handling for CDMS and Internet stuff
## 20250204  Add CDataFrame
## 20250912  Add debugging output to isolate lags in Grace processing

from datetime import datetime
print(datetime.now(),"imports.py")

from rawio.IO import *
print(datetime.now(),"rawio done")

import math
import matplotlib
global plt
import matplotlib.pyplot as plt
print(datetime.now(),"math, matplotlib (x2) done")

global pd
import pandas as pd
pd.set_option("display.max_row", 100)
print(datetime.now(),"pandas done")

###from pylab import *
print(datetime.now(),"pylab skipped (overwrites datetime?)")

import scipy.io as sio
from scipy import interpolate
from scipy.optimize import curve_fit
print(datetime.now(),"scipy (x3) done")

import ROOT
print(datetime.now(),"ROOT done")

global np
import numpy as np
print(datetime.now(),"numpy done")

from cats.cdataframe import CDataFrame
print(datetime.now(),"cats done")

import sys
import os, os.path 
import fnmatch
import glob
from sklearn.metrics import auc
print(datetime.now(),"system stuff and sklearn done")

# Additional settings needed for CDMS work
global CDMS_SUPERSIM
CDMS_SUPERSIM = os.environ['CDMS_SUPERSIM']

os.environ['HTTP_PROXY'] = '10.76.5.24:8080'
os.environ['HTTPS_PROXY'] = '10.76.5.24:8080'
os.environ['FTP_PROXY'] = '10.76.5.24:8080'
