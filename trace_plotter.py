#!/usr/bin/env python3
"""
Usage: trace_plotter.py [-h] [-v] [-t] [-e #] [-s TES|FET] -o <output> <rootfile>

Generates PNG file (-o) with TES or FET (-s) traces from specified event (-e)
as either overlay or tiled (-t) display.  For diagnostic information,
use the verbose (-v) option.

Requires: dmc module from SuperSim
          CDataFrame from CATs
          ROOT
          Numpy
          Matplotlib
"""
# 20250314  Michael Kelsey (TAMU)

### CONFIGURATION ###

import dmc
import os, sys
import matplotlib.pyplot as plt

### UTILITY FUNCTIONS ###

def niceGrid(nchan):
    """Factor value 'nchan' into rectangular grid, either NxN or Nx(N+1).
       There may be empty slots in the bottom row."""
    sqrtn = int(ceil(sqrt(nchan)))
    lo,hi = sqrtn-1,sqrtn+1
    if (lo*sqrtn >= nchan): return lo,sqrtn
    if (sqrtn*sqrtn >= nchan): return sqrtn,sqrtn
    return sqrtn,hi

def loadTrace(reader,sensor,event,channel=None):
    """Maps sensor type string to correct call to traceReader."""
    if sensor=="TES":
        return reader.plottableTES(event) if channel is None else reader.TES(event,channel)

    if sensor=="FET":
        return reader.plottableFET(event) if channel is None else reader.FET(event,channel)

    return None
                                        
def plotTiles(file,outname,title=None,sensor="TES",event=0,verbose=False):
    """Make rectangular grid of single-channel plots."""
    reader = dmc.traceReader(file,verbose)
    chans = reader.channels(sensor)
    bins = reader.timeBins(sensor)

    rows,cols = niceGrid(len(chans))
    fig, axes = subplots(rows,cols,figsize=(9,6),dpi=200)
    fig.set_tight_layout(True)

    for chan in range(len(chans)):
        currentAxis = axes.flatten()[chan]
    
        # Get TES trace for plotting
        trace = loadTrace(reader,sensor,event,chan)
        currentAxis.plot(bins, trace, lw=1, color = 'C0', label=chans[chan])
    
        currentAxis.set_xlim([-100.,2000.])
        currentAxis.set_xlabel("Time [\u03BCs]")
        currentAxis.set_ylabel("Trace [\u03BCA]")  ## FIXME for FET
        currentAxis.legend()     # Puts channel name inside box

    if title is None: title = os.path.splitext(os.path.basename(file))[0]
    fig.suptitle(title)

    if '.png' in outname:
        outname = os.path.splitext(os.path.basename(outname))[0]
    fig.savefig(outname+".png")

def plotOverlay(file,outname,title=None,sensor="TES",event=0,verbose=False):
    """Make overlay plot of readout channels."""
    fig = plt.figure(figsize=(6,4),dpi=200)

    reader = dmc.traceReader(file,verbose)
    bins = reader.timeBins(sensor)
    chans,traces = loadTrace(reader,sensor,event)

    currentAxis = plt.gca()
    currentAxis.plot(bins,traces,lw=1,label=chans)
    
    currentAxis.set_xlim([-100.,2000.])
    currentAxis.set_xlabel("Time [\u03BCs]")
    currentAxis.set_ylabel("Trace [\u03BCA]")  ## FIXME for FET
    currentAxis.legend()     # Puts channel name inside box

    if title is None: title = os.path.splitext(os.path.basename(file))[0]
    fig.suptitle(title)
    fig.savefig(outname+".png")

### MAIN PROGRAM ###

def main():
    args = getargs()

    # Separate "overlay" from list to avoid passing unrecognized parameters
    goodargs = { k:args[k] for k in args.keys()-['overlay'] }

    if args['overlay']:
        plotOverlay(**goodargs)
    else:
        plotTiles(**goodargs)

def getargs():
    """Returns arguments from the command line as a dictionary, for easier use.
       Output: settings = { 'file':    <name-of-DMC-file>,
                            'outname': <name-of-output-plot>,
                            'event':   <event number>, from -e>,
                            'sensor':  <sensor type>, from -s>,
                            'verbose': <True|False>, from -v>,
                            'overlay': <True|False>, from -t,
                           }
    """
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'e:o:s:thv')
    except getopt.GetoptError as err:
        sys.exit(2)

    if args is None or len(args)==0:
        help(__name__)
        sys.exit(2)
        
    settings = {'file':    args[0],    # Filename is required
                'outname': None,       # Ouptut name is requred
                'event':   0,          # First event, first channel
                'sensor':  "TES",      # TES or FET
                'overlay': True,       # False to get one channel per tile
                'verbose': False,
               }

    for o,a in opts:
        if   o == '-e':
            settings['event'] = int(a)
        elif o == '-h':
            help(__main__)
            sys.exit(0)
        elif o == '-o':
            settings['outname'] = a
        elif o == '-s':
            settings['sensor'] = a
        elif o == '-t':
            settings['overlay'] = False
        elif o == '-v':
            settings['verbose'] = True

    if settings['outname'] is None:
        help(__name__)
        sys.exit(1)
        
    if settings['verbose']:
        print(f"settings = {settings}")
        
    return settings


### COMMAND-LINE CALL ###

if __name__ == "__main__":
    main()
    
### END OF FILE ###
