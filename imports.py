## imnports.py -- does all the import statements for my Jupyter analysis
## Usage: from imports import * in Jupyter cell

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

import sys
import os, os.path 
import fnmatch
import glob
from sklearn.metrics import auc