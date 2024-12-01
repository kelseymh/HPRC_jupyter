#!/usr/bin/env python3

"""Load TES or FET traces from input file, for plotting, fitting, etc.
   Requires: Numpy, Pandas, ROOT
   
   Michael Kelsey <kelsey@tamu.edu>, Texas A&M University 2024
   
   20240120  Extracted functions from Jupyter notebooks
   20240202  Add support for non-signal traces (DataType)
   20241026  Add function to get dictionary of all traces in detector/event
"""

import numpy as np
import pandas as pd
import ROOT
import os

global CDMS_SUPERSIM
CDMS_SUPERSIM = os.environ['CDMS_SUPERSIM']

### DIAGNOSTIC OUTPUT FUNCTIONS ###

verbose=False                   # Global variable, value set in getargs()
def printVerbose(string):
    if verbose: print(string)

def setVerbose(vb=True):        # To set the global from outside
    global verbose
    verbose = vb

### LOAD CHANNEL NAMES FOR CONVENIENT PLOTTING ###

def getChannels(file, sensor, event):
    """Get list of channel names for TES or FET traces"""
    printVerbose(f"getChannels(file='{file}', event={event}, sensor={sensor})")

    branches = ['ChanName']
    chanNames = ROOT.RDataFrame(f"G4SimDir/g4dmc{sensor}", file, branches)\
    .Filter(f"DataType==0 & EventNum=={event}")\
    .AsNumpy(branches)
    convertRDFBranches(chanNames)

    return chanNames['ChanName']

def getBins(file, sensor, event):
    """Return just the bins array for TES or FET channels in a given event"""
    printVerbose(f"getBins(file='{file}', sensor='{sensor}', event={event})")
    
    numBins, T0, BinWidth = loadBinning(file, sensor, event)
    bins = np.arange(numBins)*BinWidth+T0    # More accurate than floating arange
    return bins


### LOAD ALL TRACES FOR A GIVEN DETECTOR/EVENT ###

def getAllTES(file, event):
    """Extract normalized (flipped and zero-baseline) TES traces from input file
       for specified event
       Output: bins    = Numpy array of bin edges, for fitting and plotting
               traces  = Numpy dict of bin value arrays, keyed on channel name
               I0      = Dict of baseline current of raw trace (subtracted away)
               PhononE = Dict of energy collected on channels, keyed on name
    """
    printVerbose(f"getAllTES(file='{file}', event={event}")

    # Get list of channels to use as keys and for looping
    channels = getChannels(file, "TES", event)

    traces = {}
    I0 = {}
    PhononE = {}
    for ich in range(len(channels)):
        chan = channels[ich]
        bins, traces[chan], I0[chan], PhononE[chan] = getTES(file, event, ich)

    return bins, traces, I0, PhononE

def getAllFET(file, event):
    """Extract normalized FET traces from input file for specified event.
       Output: bins    = Numpy array of bin edges, for fitting and plotting
               trace   = Numpy dict of bin value arrays, keyed on channel name
               ChargeQ = Dict of Ramo charge collected on channels, keyed on name
    """
    printVerbose(f"getAllTES(file='{file}', event={event}")
    
    # Get list of channels to use as keys and for looping
    channels = getChannels(file, "FET", event)

    traces = {}
    ChargeQ = {}
    for ich in range(len(channels)):
        chan = channels[ich]
        bins, traces[chan], ChargeQ[chan] = getFET(file, event, ich)

    return bins, traces, ChargeQ

def getPlottableTraces(file, sensor, event):
    """Extract all <sensor> traces from specified file.  Returns three
       arrays/lists: channel names, time bins, and a transposed array
       of the traces (i.e., instead of eight arrays of 4096 values, the
       transpose is one array of 4096 entries, each entry is a list of
       the eight values for that time bin)."""
    chans = getChannels(file, sensor, event)
    traces = [[] for i in range(len(chans))]  # Make slot for each trace
    for i in range(len(chans)):
        if (sensor == "TES"):
            bins, traces[i] = getRawTES(file, event, i)
            I0, _,_ = bestI0(traces[i])
            traces[i] = I0 - traces[i]
        if (sensor == "FET"):
            bins, traces[i] = getRawFET(file, event, i)

    plottable = list(zip(*traces))      # This is direct Python transpose

    return chans, bins, plottable


### LOAD TRACES AS RAW OR NORMALIZED ###

def getTES(file, event, chan):
    """Extract normalized (flipped and zero-baseline) TES trace from input file
       for specified event and channel index
       Output: bins    = Numpy array of bin edges, for fitting and plotting
               trace   = Numpy array of bin values, baseline subtract and flip
               I0      = Baseline current of raw trace (subtracted away)
               PhononE = Energy collected on channel, for normalization
    """
    printVerbose(f"getTES(file='{file}', event={event}, chan={chan})")

    # Get the TES signal trace and binning; TES trace is flipped and rebaselined
    bins, trace = getRawTES(file, event, chan, 0)
    start   = np.nonzero(bins>=0.)[0][0]
    I0, _,_ = bestI0(trace)
    PhononE = getEventSums(file, event, 'PhononE')[chan]   # Single channel collected energy
    printVerbose(f" I0={I0:.4e} uA, PhononE={PhononE:.4e} eV")
    
    trace = I0-trace        # Normalize: upward peak, zero baseline
    
    return bins, trace, I0, PhononE
    
def getRawTES(file, event, chan, datatype=0):
    """Extract raw TES trace information from input file, for specified event
       and channel index
       Output: bins    = Numpy array of bin edges, for fitting and plotting
               trace   = Numpy array of bin values, baseline subtract and flip
    """
    printVerbose(f"getRawTES(file='{file}', event={event}, chan={chan}, datatype={datatype})")
    
    # Get the TES trace and binning information (trace may be a list of lists)
    bins, trace = loadTrace(file, event, chan, "TES", datatype);
    
    return bins, trace


def getFET(file, event, chan):
    """Extract normalized FET trace information from input file, for specified event
       and channel index
       Output: bins    = Numpy array of bin edges, for fitting and plotting
               trace   = Numpy array of bin values, baseline subtract and flip
               ChargeQ = Value of Ramo charge collected on channel (units of e+)
    """
    printVerbose(f"getFET(file='{file}', event={event}, chan={chan})")
    
    # User gives absolute path of ROOT file and the g4dmcFET tree is opened
    bins, trace = getRawFET(file, event, chan, 0)
    ChargeQ = getEventSums(file, event, 'ChargeQ')[chan]   # Single channel collected charge
    printVerbose(f" ChargeQ={ChargeQ:.4e} eV")

    if abs(trace.max())<abs(trace.min()): trace = -trace    # Flip downward going traces
        
    return bins, trace, ChargeQ

def getRawFET(file, event, chan, datatype=0):
    """Extract raw FET trace (electrons upward, holes downward) from input file, for
       for specified event and channel index
       Output: bins    = Numpy array of bin edges, for fitting and plotting
               trace   = Numpy array of bin values, may be list of lists
    """
    printVerbose(f"getRawFET(file='{file}', event={event}, chan={chan}, datatype={datatype})")

    # User gives absolute path of ROOT file and the g4dmcFET tree is opened
    bins, trace = loadTrace(file, event, chan, "FET", datatype)
    
    return bins, trace

def loadTrace(file, event, chan, sensor, datatype=0):
    """Extract requested trace (TES or FET, signal or diagnostic data, from event in file
       NOTE: channel may include a decimal subchannel, or reference all subchannels.
    """
    printVerbose(f"loadTrace(file='{file}', event={event}, chan={chan},"
                 f" sensor={sensor}, datatype={datatype})")

    bins = getBins(file, sensor, event)
  
    if type(datatype) is str: datatype = getDatatype(datatype)

    filter = f"DataType=={datatype} & EventNum=={event}"
    if hasSubchans(datatype) and type(chan) is not int:
        filter += f" & ChanNum=={chan}"
    else:
        filter += f" & int(ChanNum)=={chan}"

    printVerbose(f"Applying filter '{filter}'")

    branches = ['Trace']
    g4dmcTrace = ROOT.RDataFrame(f"G4SimDir/g4dmc{sensor}", file, branches)\
    .Filter(filter).AsNumpy(branches)
    convertRDFBranches(g4dmcTrace)

    trace = g4dmcTrace['Trace']
    if not hasSubchans(datatype):      # Unroll single-entry list for trace
        trace = g4dmcTrace['Trace'][0]

    return bins, trace

def loadBinning(file, sensor, event, chan=0):
    """Retrieve binning information for TES or FET traces
       Output: NBins    = Number of time bins per trace
               T0       = Start time of trace (us)
               BinWidth = Spacing between time bins (us)"""
    printVerbose(f"loadBinning(file='{file}', sensor='{sensor}', event={event},"
                 f" chan={chan})")

    branches = ['Trace','BinWidth','T0']          # Need trace to get bin count
    g4dmcTrace = ROOT.RDataFrame(f"G4SimDir/g4dmc{sensor}", file, branches)\
    .Filter(f"DataType==0 & EventNum=={event} & ChanNum=={chan}")\
    .AsNumpy(branches)
    convertRDFBranches(g4dmcTrace)
    
    numBins  = len(g4dmcTrace['Trace'][0])        # number of bins
    T0       = g4dmcTrace['T0'][0]*1e-3           # start time of trace (ns -> us)
    BinWidth = g4dmcTrace['BinWidth'][0]*1e-3     # width of bin (ns -> us)

    return numBins, T0, BinWidth

### UTILITY FUNCTIONS ###

def bestI0(trace):
    """Compute 'adaptive I0' for trace, scanning the supposed pre-trigger
       baseline for the range of values with the smallest RMS.  This
       should exclude the region where the real trace starts.  Returns
       computed I0 value, along with index of average and RMS value."""
    start = 5                   # Need some bins to compute RMS

    eom = [np.std(trace[:i])/np.sqrt(i) for i in range(start,len(trace))]
    ibest = eom.index(min(eom))+start

    return np.mean(trace[:ibest]), ibest, eom[ibest]

### G4DMC READOUT ENUMERATORS (See CDMSinfo/CDMSReadoutTypes.hh) ###

typeCodes = { "Signal": 0, "Data": 0, "IbIs": 1, "Vb": 2, "TES_A": 3,
              "PowerIn": 4, "FET_Vdet": 5, "Template": 6, "Crosstalk": 7,
              "RamoScale": 8, "PowerOut": 9, "TES_G": 10, "TES_R": 11,
              "TES_T": 12, "TES_Tc": 13 }

typeNames = {value: key for key, value in typeCodes.items()}

typeSubchan = { 0: False, 1: False, 2: False, 3: True, 4: True, 5: False,
                6: True, 7: False, 8: False, 9: True, 10: True, 11: True,
                12: True, 13: True }
               
def getDatatype(dtype):
    """Convert between data type name strings and integer codes"""
    return typeCodes[dtype] if type(dtype) is str else typeNames[dtype]

def hasSubchans(dtype):
    """Flag whether the specified data type has decimal subchannels"""
    return typeSubchan[dtype] if type(dtype) is int else typeSubchan(typeCodes[dtype])
    
    
### G4DMC SUMMARY DATA ###

def getEventSums(file, event, branch):
    """Extract PhononE or ChargeQ branch from g4dmcEvent for specified event"""
    printVerbose(f"getEventSums(file='{file}', event={event}, branch='{branch}')")

    # Unwrap Numpy array-of-arrays to return just the list of channel values
    eventSums = ROOT.RDataFrame('G4SimDir/g4dmcEvent', file, [branch])\
    .Filter(f"EventNum=={event}").AsNumpy([branch])
    convertRDFBranches(eventSums)
    
    return eventSums[branch][0]     # Only one event in array of arrays


def getTemplate(detname, chan, sensor):
    """Extract channel template for specified detector, as Numpy array"""
    printVerbose(f"getTemplate(detname='{detname}', chan={chan}, sensor={sensor})")
    
    if detname is None or detname == "": return None
    
    templatePath = f"{CDMS_SUPERSIM}/CDMSgeometry/data/dmc/{detname}/{sensor}Templates"
    if not os.path.isfile(templatePath): return None

    colname = {"TES":"Traces", "FET":"Trace"}[sensor]
    templates = pd.read_csv(templatePath,sep="\t")
    template = templates.loc[chan,colname].split()
    template = np.array([float(i) for i in template])

    printVerbose(" got template from {} with {} bins".format(templatePath, len(template)))
    return template


### SPECIAL: Data converter copied out of CATs/cdataframe/methods.py AsNumpy()

def convertRDFBranches(rdf):
    """Convert all of the branches in the RDataFrame to Numpy arrays if possible."""
    for br,arr in rdf.items():
        rdf[br] = convertRDFArray(arr)
    return rdf      # Changes above are applied in situ; this is just for convenience

def convertRDFArray(arr):
    """Convert a NumPy array of C++ objects to Python objects if possible."""
    def isinstancecppyy(obj, classname):
        """Tests whether "obj" belongs to the C++ class named "classname".

           cppyy unfortunately provides no way of checking if an object comes from a
           class template, so the best we can do is a string comparison of its name.

           Use "." in place of "::" for namespacing (e.g. "std.vector" instead of
           "std::vector").
        """
        return str(type(obj)).startswith(f"<class cppyy.gbl.{classname}<")

    if len(arr) == 0: return arr # No type information and nothing to convert
    elif isinstancecppyy(arr[0], "ROOT.VecOps.RVec") or\
         isinstancecppyy(arr[0], "std.vector"):
        # RVec<char> needs special treatment because .data() returns a str
        converted = np.vectorize(*
            [lambda rvec: rvec.data()] if isinstance(arr[0].data(), str) else
            [np.array, [object]]
        )(arr)

        # Attempt to collapse nested NumPy arrays to a single multidimensional array
        try              : return np.array(converted.tolist())
        except ValueError: return converted # Arrays are different lengths, not possible

    # Why do bool columns return arrays of dtype object?
    elif isinstance     (arr[0], bool     ): return np.array(arr, dtype=bool)
    elif isinstancecppyy(arr[0], "std.set"): return np.vectorize(set)(arr)
    else                                   : return arr # Unrecognized type; leave it alone
    
### END OF FILE ###
