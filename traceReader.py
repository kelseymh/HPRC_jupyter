"""Load TES or FET traces from input file, for plotting, fitting, etc.
   Requires: Numpy, Pandas, CATs from SuperCDMS

   Usage: import traceReader
          myReader = traceReader("path/to/DMC-file.root")
          bins = myReader.timeBins("TES")
          tracePAS1 = myReader.TES(event=2,chan="PAS1",det=0)

   Michael Kelsey <kelsey@tamu.edu>, Texas A&M University 2024
   
   20241028  Adapted from traces_rdf.py to reduce RDF memory overhead
   20241201  Swap det and chan arguments in binning functions.
   20241202  Use try/except to support non-existent trace TTrees.
   20241204  Drop channel argument for time binning (all channels match).
	     Add support for channel sets and event sets; new TESsum()
             function to return 'PTOF' sum of traces.
"""

import numpy as np
import pandas as pd
from cats.cdataframe import CDataFrame
import os

global CDMS_SUPERSIM
CDMS_SUPERSIM = os.environ['CDMS_SUPERSIM']

### MAIN DRIVER CLASS -- ONE PER FILE ###

class traceReader:
    """Load DMC traces from DMCintermediate ROOT file.  One reader should be
       created per input file, or set of input files.  Reader loads and returns
       Numpy arrays or dictionaries of traces, selected by event number,
       detector number, or/and sensor type ("TES" or "FET").

       Constructor arguments:
       files = Name of DMC ROOT file, or list of names (to use as TChain)
       verbose = (optional) True reports function calls, progress

       Most functions will take three or four arguments:
       sensor   = String name for trace type, "TES" or "FET"
       event    = (default=0) Event number to read (-1 means array of all)
       channel  = (default=0) Channel index or name to read (-1 means all)
       detector = (default=0) Detector number to read

       If the sensor type is in the function name, it is not passed as an
       argument.
    """

     # Constructor
    def __init__(self, files, verbose=False):
        """Constructor: Creates RDataFrames for TES and FET traces."""
        self.files = files
        self.verbose = verbose

        try:		# Better if CDataFrame had a flag for this
            self.tesDF = CDataFrame("G4SimDir/g4dmcTES", files)
        except:
            self.tesDF = None

        try:		# Better if CDataFrame had a flag for this
            self.fetDF = CDataFrame("G4SimDir/g4dmcFET", files)
        except:
            self.fetDF = None

        self.dfs = { 1: self.tesDF, "TES": self.tesDF,
                     2: self.fetDF, "FET": self.fetDF }
        self.evtDF = CDataFrame("G4SimDir/g4dmcEvent", files)

        return

    # Useful diagnostic functions

    def setVerbose(self, verbose=True):
        self.verbose = verbose

    def printVerbose(self, string):
        if self.verbose: print(string)

    def __repr__(self):
        return f"<traceReader('{self.files}','{self.verbose}')"


    # Get channel information for specified detector

    def channels(self, sensor, det=0):
        """List of channel names (in index order) for specified sensor."""
        self.printVerbose(f"channels('{sensor}', det={det})")

        # Find first event containing requested detector, then get all channels
        filt = f"DetNum=={det} & DataType==0"
        event = self.dfs[sensor].Filter(filt).Range(0,1).AsNumpy(["EventNum"])["EventNum"][0]

        filt += f" & EventNum=={event}"
        chanNames = self.dfs[sensor].Filter(filt).AsNumpy(["ChanName"])

        self.printVerbose(f" returning {chanNames['ChanName']}")
        return chanNames['ChanName']

    def channum(self, chan, det=0):
        """Converts named channel (or floating point subchannel) to index
           corresponding to channel; always returns an int, chopping off
           any subchannel info.  Sensor type is inferred from channel
           name, if that was the input."""
        self.printVerbose(f"channum('{chan}', det={det})")

        if isinstance(chan, str):
            idot = chan.find(".")
            nosub = chan[:idot] if idot>=0 else chan	# Strip optional subchannel
            self.printVerbose(f"chan is string; nosub={nosub}")
            
            teschans = self.channels("TES", det)
            if nosub in teschans: return list(teschans).index(nosub)

            fetchans = self.channels("FET", det)
            if nosub in fetchans: return list(fetchans).index(nosub)

            return -1
        else:
            self.printVerbose(f"chan is number: {chan}")
            
        return int(chan)       			# Strip optional subchannel

    def subchannum(self, chan, det=0):
        """For 'NAME.index' subchannels, converts the name to a number and
           returns the resulting float.  Different from channum(), which
           only returns the integer part."""
        if not isinstance(chan, str): return chan

        idot = chan.find(".")
        subchan = float(chan[chan.find("."):]) if idot>=0 else 0.
        
        return self.channum(chan)+subchan

    # Get time binning for specified detector

    def timeBins(self, sensor, det=0):
        """Return array of time bin edges for specified sensor type."""
        self.printVerbose(f"timeBins('{sensor}', det={det})")

        nBins, T0, dT = self.binning(sensor, det)
        bins = np.arange(nBins)*dT+T0      # More accurate than floating arange
        return bins

    def binning(self, sensor, det=0):
        """Return time binning parameters for specified sensor type.
           Output: numBins  = Number of time bins per trace
                   T0       = Start time of trace (us)
                   BinWidth = Spacing between time bins (us)"""
        self.printVerbose(f"binning('{sensor}', det={det})")

        branches = ['Trace','BinWidth','T0']    # Need trace to get bin count
        filt = f"DetNum=={det} & DataType==0 & ChanNum==0"
        binData = self.dfs[sensor].Filter(filt).Range(0,1).AsNumpy(branches)

        numBins  = len(binData['Trace'][0])     # number of bins
        T0       = binData['T0'][0]*1e-3        # start time of trace (ns -> us)
        BinWidth = binData['BinWidth'][0]*1e-3  # width of bin (ns -> us)

        self.printVerbose(f"numBins {numBins}, T0 {T0}, binWidth {BinWidth}")

        return numBins, T0, BinWidth

    
    # Load baseline-normalized traces, or plottable arrays

    def TES(self, event, chan=0, det=0):
        """Read single-channel TES trace for specified event, normalize to
           upward going with zero baseline (analysis style).

           If chan=-1, returns all traces as Numpy dictionary with allTES()
           If event=-1, returns array of traces
           If both are -1, returns dictionary of array of traces
        """
        if (self.channum(chan,det)<0): return self.allTES(event,det)

        self.printVerbose(f"TES({event}, chan={chan}, det={det})")

        if (event<0):		# All events, loop for normalization
            events = self.getEvents(det)
            return [ self.TES(ev,chan,det) for ev in events ]

        raw = self.rawTES(event, chan, det)
        I0, _,_ = self.bestI0(raw)		# Discard RMS and index
        return I0-raw

    def TESsum(self, event, det=0):
        """Sum all TES traces for specified event, after normalizing to
           upward going with zero baseline ('PTOF' version).

           If event=-1, returns array of traces
        """
        self.printVerbose(f"TESsum({event}, det={det})")

        if (event<0):		# All events, loop for normalization
            events = self.getEvents(det)
            return [ self.TESsum(ev,det) for ev in events ]

        tesDict = self.allTES(event, det)	# Dictionary keyed by channel
        traces = np.array(list(tesDict.values()))
        return np.sum(traces, axis=0)		# Sums traces bin-by-bin

    def FET(self, event, chan=0, det=0):
        """Read single-channel FET trace for specified event, normalize to
           upward going with zero baseline (analysis style).

           If chan=-1, returns all traces as Numpy dictionary with allFET()
           If event=-1, returns array of traces
           If both are -1, returns dictionary of array of traces
        """
        if (self.channum(chan,det)<0): return self.allFET(event,det)

        self.printVerbose(f"FET({event}, chan={chan}, det={det})")

        trace = self.rawFET(event, chan, det)

        if (event<0):		# All events: normalize each trace in array
            for i in range(0,len(trace)):
                if abs(trace[i].max())>abs(trace[i].min()): trace[i] = -trace[i]
        else:
            if abs(trace.max())>abs(trace.min()): trace = -trace

        return trace

    def allTES(self, event, det=0):
        """Load baseline-normalized TES traces for specified event into
           Numpy dictionary, keyed with channel names.

           If event<0, returns dictionary of lists of traces for all events.
        """
        self.printVerbose(f"allTES({event}, det={det})")

        chans = self.channels("TES", det)
        traceSet = { ch: self.TES(event,ch,det) for ch in chans }
        return traceSet

    def allFET(self, event, det=0):
        """Load baseline-normalized TES traces for specified event into
           Numpy dictionary, keyed with channel names.

           If event<0, returns dictionary of lists of traces for all events.
        """
        self.printVerbose(f"allFET({event}, det={det})")

        chans = self.channels("FET", det)
        traceSet = { ch: self.FET(event,ch,det) for ch in chans }
        return traceSet


    # Convert dictionary of channel traces to plottable stack

    def plottableTES(self, event, det=0):
        """Load all TES traces for specified event/detector, convert to
           "transpose" array suitable for passing to matplotlib.
           Output: channels = List of channel names to use with "label="
                   plottable = Array of timebin entries, each with array
                               of values from each channel for that bin."""
        self.printVerbose(f"plottableTES({event}, det={det}")
        return self.channels("TES",det), self.plottable(self.allTES(event,det))

    def plottableFET(self, event, det=0):
        """Load all FET traces for specified event/detector, convert to
           "transpose" array suitable for passing to matplotlib.
           Output: channels  = List of channel names to use with "label="
                   plottable = Array of timebin entries, each with array
                               of values from each channel for that bin."""
        self.printVerbose(f"plottableFET({event}, det={det}")
        return self.channels("FET",det), self.plottable(self.allFET(event,det))

    def plottable(self, traceSet):
        """Convert either dictionary or 2D array of traces into a "transpose"
           array (stack of N values for all channels x time bins) for passing
           into matplotlib."""
        self.printVerbose(f"plottable(traceSet={type(traceSet)})")

        if isinstance(traceSet, dict):		# Dictionary keyed on channel
            return self.plottable(list(traceSet.values()))

        return list(zip(*traceSet))

    # Related information about traces

    def TES_I0(self, event, chan=0, det=0):
        """Returns adaptive average (minimum uncertainty) of TES baseline."""
        self.printVerbose(f"TES_I0({event}, chan={chan}, det={det}")

        I0, _,_ = self.bestI0(self.rawTES(event,chan,det))    # Ignore err,index
        return I0


    # Load raw traces (called from high-level functions)

    def rawTES(self, event, chan=0, det=0):
        """Read single-channel TES trace in original DMC format (I0 baseline,
           downward going pulse).  If event<0, all traces are returned as
           array of arrays."""
        self.printVerbose(f"rawTES({event}, chan={chan}, det={det})")
        return self.loadTrace(event, "TES", chan, det, dtype=0)

    def rawFET(self, event, chan=0, det=0):
        """Read single-channel FET trace in original DMC format (zero baseline,
           upward or downward going pulse, depending on charge carrier and
           detector bias).  If event<0, all traces are returned as array of
           arrays."""
        self.printVerbose(f"rawFET({event}, chan={chan}, det={det})")
        return self.loadTrace(event, "FET", chan, det, dtype=0)

    def loadTrace(self, event, sensor, chan=0, det=0, dtype=0):
        """Extract requested trace (TES or FET, signal or diagnostic data,
           from event in file.  If event<0, all traces are returned as array
           of arrays.
           NOTE: channel may reference a decimal subchannel, or request all
                 subchannels as an array of traces.
        """
        self.printVerbose(f"loadTrace({event}, '{sensor}', chan={chan},"
                          f" det={det}, dtype={dtype})")

        filt = f"EventNum=={event} & " if event>=0 else ""
        filt += f"DetNum=={det} & DataType=={dtype} & "
        filt += f"ChanNum=={self.subchannum(chan)}"

        self.printVerbose(f"Applying filter '{filt}'")

        # NOTE: Final bracket unwraps dictionary, to get value array of arrays
        trace = self.dfs[sensor].Filter(filt).AsNumpy(["Trace"])["Trace"]

	# Strip outer array if single-value result expected
        return trace if (self.hasSubchans(dtype) or event<0) else trace[0]

    def bestI0(self, trace):
        """Compute 'adaptive I0' for trace, scanning the supposed pre-trigger
           baseline for the range of values with the smallest RMS.  This
           should exclude the region where the real trace starts.
           Returns computed I0 value, along with uncertainty, and end
           index of range.

           If a list of traces is provided, then three lists, of the three
           computed values, will be returned."""
        self.printVerbose(f"bestI0({type(trace)})")

        if np.asarray(trace).ndim == 2:		# List of traces
            # Builds array of [[I0,eom,index], [I0,eom,index], ...]
            I0transpose = [self.bestI0(trace[i]) for i in range(len(trace))]

            # Convert to three arrays [I0,...], [eom,...], [index,...]
            return list(zip(*IOtranspose))

        start = 5                   # Need some bins to compute RMS
        eom = [np.std(trace[:i])/np.sqrt(i) for i in range(start,len(trace))]
        ibest = eom.index(min(eom))+start
        I0 = np.mean(trace[:ibest])
        self.printVerbose(f"I0 {I0} +- {eom[ibest]}, [0:{ibest}]")

        return I0, eom[ibest], ibest


    # Event level summary information, related to traces

    def getEvents(self, det=0):
        """Returns list of event numbers for specified detector."""
        self.printVerbose(f"getEvents(det={det})")

        filt = f"DetNum=={det}"
        return self.evtDF.Filter(filt).AsNumpy(["EventNum"])["EventNum"]

    def getPhononE(self, event, det=0):
        """Returns array of PhononE values for all channels in detector.
           User can use `traceReader.channels("TES", event, det)` to get
           list of channels and unpack the array."""
        filt = f"EventNum=={event} & DetNum=={det}"
        phE = self.evtDF.Filter(filt).AsNumpy(["PhononE"])["PhononE"][0]
        self.printVerbose(f"getPhononE({event}, det={det}) = {phE}")

        return phE

    def getChargeQ(self, event, det=0):
        """Returns array of ChargeQ values for all channels in detector.
           User can use `traceReader.channels("FET", event, det)` to get
           list of channels and unpack the array."""
        filt = f"EventNum=={event} & DetNum=={det}"
        chgQ = self.evtDF.Filter(filt).AsNumpy(["ChargeQ"])["ChargeQ"][0]
        self.printVerbose(f"getChargeQ({event}, det={det}) = {chgQ}")

        return chgQ

    # G4DMC enumerator/string mappings

    _typeCodes = {
        "Signal": 0, "Data": 0, "IbIs": 1, "Vb": 2, "TES_A": 3, "PowerIn": 4, 
        "FET_Vdet": 5, "Template": 6, "Crosstalk": 7, "RamoScale": 8,
        "PowerOut": 9, "TES_G": 10, "TES_R": 11,  "TES_T": 12, "TES_Tc": 13 }

    _typeNames = {value: key for key, value in _typeCodes.items()}

    _hasSubchan = { 0: False, 1: False, 2: False, 3: True, 4: True,
                    5: False, 6: True, 7: False, 8: False, 9: True,
                    10: True, 11: True, 12: True, 13: True }

    def getDatatype(self, dtype):
        """Convert between data type name strings and integer codes"""
        return (self._typeCodes[dtype] if type(dtype) is str
                else self._typeNames[dtype])

    def hasSubchans(self, dtype):
        """Flag whether the specified data type has decimal subchannels"""
        return (self._hasSubchan[dtype] if type(dtype) is int
                else self._hasSubchan[self._typeCodes[dtype]])
