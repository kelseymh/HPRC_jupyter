# Utility functions to support getting traces out of DMC output files
#
# Usage: import traces_numpy

import root_numpy
import numpy as np

# Extract desired readout-trce branches from ROOT file for processing
def getBranches(file, sensor, branches, selection=None):
    data = root_numpy.root2array(file, treename="G4SimDir/g4dmc"+sensor, branches=branches,
                                 selection=selection)
    return data

# Extract TES or FET channel names from file
def getChannels(file, sensor):
    data = getBranches(file, sensor, branches=['ChanNum','ChanName'],selection="DataType==0")
    chans = np.unique(data)['ChanName']             # Preserves ChanNum sorting order
    chans = [ch.decode("ascii") for ch in chans]    # Converts b'...' to normal string
    return chans

# Extract TES or FET channel binning from file
def getBinning(file, sensor):
    data = getBranches(file, sensor, branches=['T0','BinWidth','Trace'])
    nbin = len(data['Trace'][0])
    t0 = data['T0'][0]
    dt = data['BinWidth'][0]   

    params = { 'NBins': nbin, 'T0': t0, 'BinWidth': dt }
    return params
    
def getBins(file, sensor):
    params = getBinning(file, sensor)
    
    # Create bin edges to match length of trace
    nbin = params['NBins']
    t0 = params['T0']
    dt = params['BinWidth']
       
    # Pre-create dictionary of lists of traces for each channel
    bins = { chan: list() for chan in getChannels(file, sensor) }
    for ch in bins.keys():
        bins[ch] = np.arange(t0, t0+dt*nbin, dt)
        
    return bins

# Generic function used by getTES() and getFET() below
# sensor="TES" or "FET"; function will take care of directory
def getTraces(file, sensor, dtype=0):
    data = getBranches(file, sensor, branches=['EventNum','Trace','ChanName'],
                       selection=f"DataType=={dtype}")
    
    # Pre-create dictionary of lists of traces for each channel
    traces = { chan: list() for chan in getChannels(file, sensor) }

    # Copy traces for each event into new array indexed by event number
    for i in range(len(data)):
        for ch in traces.keys():
            if data['ChanName'][i] == bytes(ch, "ascii"):
                traces[ch].append(data['Trace'][i])
    
    for ch in traces.keys():
        traces[ch] = np.array(traces[ch])

    return traces

def getRawTES(file, dtype=0):
    tes = getTraces(file, "TES", dtype)
    return tes

def getTES(file, dtype=0):
    tes = getRawTES(file, dtype=0)
    
    # Convert TES traces from downward to upward going, remove baseline offset
    for ch in tes.keys():
        tes[ch] = np.array([max(tr)-tr for tr in tes[ch]])

    return tes

def getFET(file, dtype=0):
    fet = getTraces(file, "FET", dtype)
    return fet