## imnports.py -- does all the import statements for my Jupyter analysis
## Usage: from imports import *
##
## 20241202  Add environment variable handling for CDMS and Internet stuff
## 20250204  Add CDataFrame

from rawio.IO import *
import math
import matplotlib
global plt
import matplotlib.pyplot as plt
global pd
import pandas as pd
pd.set_option("display.max_row", 100)
from pylab import *
import scipy.io as sio
from scipy import interpolate
from scipy.optimize import curve_fit
import ROOT
global np
import numpy as np
from cats.cdataframe import CDataFrame

import sys
import os, os.path 
import fnmatch
import glob
from sklearn.metrics import auc

# Additional settings needed for CDMS work
global CDMS_SUPERSIM
CDMS_SUPERSIM = os.environ['CDMS_SUPERSIM']

os.environ['HTTP_PROXY'] = '10.76.5.24:8080'
os.environ['HTTPS_PROXY'] = '10.76.5.24:8080'
os.environ['FTP_PROXY'] = '10.76.5.24:8080'
