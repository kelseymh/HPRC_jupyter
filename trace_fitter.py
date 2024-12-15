#!/usr/bin/env python3

"""Processes TES or FET trace from input file, fits for shape parameters

   Command: trace_fitter.py [options]

   Use trace_fitter.py -h for details about command line arguments.

   Jupyter: from trace_fitter import traceFitter
            myfit = traceFitter("ROOT-FILE","TES",[True|False]) 
            myfit.doFit("iZIP5",<event>,<channum>,[True|False])

   Constructor True/False controls verbosity (diagnostic messages).
   doFit True/False controls whether overplay plots are created.

   Michael Kelsey <kelsey@tamu.edu>, Texas A&M University 2023
"""

# 20231210  Adapted from Warren Perry's trace_fitter.ipynb notebook
# 20240110  Adding support for FET traces, RDF input instead of root_numpy
# 20240115  Add diagnostic output, improve TES vs. FET configurations
# 20241201  Replace traces_rdf with new traceReader class; compute Ceff for FET
# 20241204  For event==-1, average traces from all events.  For channel==-1,
#	      loop over channels and report all fit results.
# 20241212  Add separate reporting functions; return curve_fit result as
#	      np.array(), rather than as tuple, append computed values for
#	      reporting or multi-event averaging.
# 20241212  Convert most functionality to class, to better share settings.

def usage():
    print("""
trace_fitter.py [-h] [-p] [-v] [-d <det>] [-e <evt>] [-s <type>] [-c <chan>] <file>
    
Reads in TES or FET traces from single-detector file, fits for shape parameters
(I0, Ipeak/E, rise and fall).
             
Argument: <file>    DMC ROOT file from single-detector simulation
Options: -d <det>   Detector type name (iZIP5, NF-C, etc.)
         -e <event> Event number from file to be processed (-1 does average)
         -c <chan>  Index (0 to N-1) of channel to be processed (-1 iterates)
         -s <type>  Sensor type (TES or FET)
         -h         [Optional] Display this usage information
         -p         [Optional] Generate plots of fit results
         -v         [Optional] Verbose output

NOTE: When "-e -1" is specified, each event is fit individually, and resulting
      parameter values are averages to get a final result.  This does not
      cleanly suppress noise/fluctuations in the traces, but avoids problems
      of trace normalization and time shifting (c.f. template generation).

Requires: Numpy, Matplotlib, SciPy, ROOT""")

### CONFIGURATION ###

from traceReader import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os, sys

global CDMS_SUPERSIM
CDMS_SUPERSIM = os.environ['CDMS_SUPERSIM']


### DRIVER CLASS TO PERFORM FITS ###

class traceFitter:
    """Driver class to fit DMC-generated traces (TES or FET) to simple shapes
       with single rise and fall time.  Uses traceReader to access events.
       Pulse shapes are available via static methods TESshape(t,tR,tF) and
       FETshape(t,invTd,invTr).

       Constructor arguments: (see also command line arguments in usage())
       file = Name of DMC ROOT file, or list of names (to use as TChain)
       sensor = 'TES' or 'FET' for kind of traces to be fit
       verbose = (optional) True reports function calls, progress

       The sensor argument may be overridden after construction; the value
       is used for by internal functions to know how to extract bin ranges and
       initial guesses for fit parameters.
    """

    # Constructor
    def __init__(self, file, sensor, verbose=False):
        """Constructor: Creates traceReader to share data access across
           functions."""
        self.reader  = traceReader(file, verbose)
        self.sensor  = sensor
        self.verbose = verbose

        return

    # Useful diagnostic functions

    def setVerbose(self, verbose=True):
        self.verbose = verbose

    def printVerbose(self, string):
        if self.verbose: print(string)

    def __repr__(self):
        return f"<traceFitter('{self.reader.files}','{self.sensor}',{self.verbose})"

    # Main fitting functions

    def doFit(self, detname="", event=0, channel=0, doplot=False):
        """Get specified TES or FET trace (event and channel) from DMC file,
           fit for shape and make overlay plots if requested.

           Arguments: detname	Detector name string used for output.
                      event	Event to process; -1 averages all events.
                      channel	Channel index; -1 iterates all channels.
                      doplot	True to generate overlay of trace and fit."""
        self.printVerbose(f"doFit(detname='{detname}', event={event},"
                          f" channel={channel}, doplot={doplot})")

        if (channel<0):
            allChans = self.reader.channels(self.sensor)
            self.printVerbose(f"Processing all {len(allChans)} channels...")
            for ichan in range(len(allChans)):
                self.doFit(detname, event, ichan, doplot)
            return

        # Following is for a single channel
        if self.sensor == "TES":
            result = self.fitTES(detname, event, channel, doplot)
            cname = self.reader.channels("TES")[channel]
            self.reportTES(detname, cname, result)
        elif self.sensor == "FET":
            result = self.fitFET(detname, event, channel, doplot)
            cname = self.reader.channels("FET")[channel]
            self.reportFET(detname, cname, result)
        else:
            print(f"Invalid sensor type '{self.sensor}' specified.")

        return

    def fitTES(self, detname="", event=0, channel=0, doplot=False):
        """Get specified TES trace (event and channel) from DMC file,
           fit for shape.  Arguments same as doFit().
           Returns [Amplitude, tR, tF, Toffset, I0, IvsE]."""
        self.printVerbose(f"fitTES(detname='{detname}', event={event},"
                          f" channel={channel}, doplot={doplot})")

        if (event<0):		# Average results from all events in file
            events = self.reader.getEvents(0)	# FIXME: Need detector number!
            resultSet = np.array([ self.fitTES(detname,ev,channel,False)
                                   for ev in np.int32(events) ])
            results = np.average(resultSet, axis=0)	# Average each column
            return results

        # Below is single-event processing
        bins  = self.reader.timeBins("TES")
        trace = self.reader.TES(event, channel)

        #### Obtain figures of merit measurements for trace and template ####
        results = self.fitTrace(bins, trace, self.TESshape, self.guessTES)

        I0       = self.reader.TES_I0(event, channel)
        PhononE  = self.reader.getPhononE(event)[channel]
        IversusE = results[0]/PhononE	# (a, tR, tF, Toffset)

        results = np.append(results, I0);
        results = np.append(results, IversusE);

        if doplot:
            self.plots(detname, "TES", channel, bins, trace,
                       self.TESshape(bins,*results))

        return results

    def fitFET(self, detname="", event=0, channel=0, doplot=False):
        """Get specified FET trace (event and channel) from DMC file,
           fit for shape.  Arguments same as doFit().
           Returns [Amplitude, invTd, invTr, Toffset, Ceff]."""
        self.printVerbose("fitFET same args as traceFit()")

        bins  = self.reader.timeBins("FET")
        trace = self.reader.FET(event, channel)

        #### Obtain figures of merit measurements for trace and template ####
        results = self.fitTrace(bins, trace, self.FETshape, self.guessFET, False)

        ChargeQ = reader.getChargeQ(event)[channel]
        Ceff    = ChargeQ*1.60218e-4 / results[0]  # (a, invTd, invTr, timeShift)
        # = Q/V, to get pF need 1e12 * coulomb/volt

        results = np.append(results, Ceff)

        if doplot:
            self.plots(detname, "FET", channel, bins, trace,
                       self.FETshape(bins,*results))

        return results

    # Report fit results for appropriate sensor type

    def reportTES(self, detname, cname, results):
        """Print TES fit results in format suitable for use in TESConstants."""
        self.printVerbose(f"reportTES(detname='{detname}', cname='{cname}',"
                          f" results={results}")
        
        a, tR, tF, offset, I0, IvsE = results     # Unroll results for reporting
        print(f'# {detname} {cname} shape parameters (to generate templates)')
        print(f'I0\t\t{I0:.4e} microampere')
        print(f'IversusE\t{IvsE:.4e} microampere/eV')
        print(f'riseTime\t{tR:.4e} us')
        print(f'fallTime\t{tF:.4e} us')
        print(f'Offset  \t{offset:.4e} us')

    def reportFET(self, detname, cname, results):
        """Print FET fit results in format suitable for use in FETConstants."""
        self.printVerbose(f"reportFET(detname='{detname}', cname='{cname}',"
                          f" results={results}")
            
        a, invTd, invTr, offset, Ceff = results      # Unroll results for reporting
        print(f'# {detname} {cname} shape parameters (to generate templates)')
        print(f'templateCeff\t{Ceff:.4e} pF')
        print(f'decayTime   \t{1./invTd:.4e} us\t# decayRate {invTd:.4e}/us')
        print(f'recoveryTime\t{1./invTr:.4e} us\t# recoveryRate {invTr:.4e}/us')
        print(f'Offset      \t{offset:.4e} us')
            

    # Idealized pulse shapes for fitting or plotting

    @classmethod
    def TESshape(cls, x, a, t_r, t_f, offset):
        """Shape of flipped TES trace above baseline, with simple rise and
           fall times.  'a' is observed peak value."""
        # Peak is where d(TESshape)/dt = 0;
        # Interestingly, this is where FETshape crosses zero
        tpeak = np.log(t_f/t_r) * t_f*t_r / (t_f-t_r)
        peak = cls._TESshape(tpeak, t_f, t_r)
        return (a/peak)*cls._TESshape(x-offset, t_f, t_r)

    @classmethod
    def _TESshape(cls, t, t_r, t_f):
        """Normalized, zero-aligned TES shape.  For internal use only."""
        if np.isscalar(t):
            return np.exp(-t/t_f)-np.exp(-t/t_r) if (t>=0) else 0.
        else:
            return np.array([cls._TESshape(ti,t_r,t_f) for ti in t])

    @classmethod
    def FETshape(cls, x, a, invTd, invTr, offset):
        """Shape of FET trace above baseline, with simple decay and recovery
           rates. 'a' is observed peak value at t=0."""
        peak = invTr**2 - invTd**2
        return (a/peak)*cls._FETshape(x-offset, invTd, invTr)

    @classmethod
    def _FETshape(cls, t, invTd, invTr):
        """Normalized, zero-aligned FET shape.  For internal use only."""
        # This is derivative of TESshape; discontinuous @ t=0
        if np.isscalar(t):
            shape = np.exp(-t*invTd)*invTd - np.exp(-t*invTr)*invTr
            return shape if (t>=0) else 0.
        else:
            return np.array([cls._FETshape(ti,invTd,invTr) for ti in t])


    # Fitting bounds and initial value estimates

    def guessTES(self, bins, trace):
        """Returns initial guesses for TES fit rise and fall time for curve_fit."""
        peak = trace.max()
        ipeak = trace.argmax()
        self.printVerbose(f"guessTES: peak {peak} @ bin {ipeak} (t {bins[ipeak]})")
    
        # Rise time: look for two e-foldings on rising side
        rlo = np.nonzero(trace[:ipeak]<=0.1*peak)[0][-1]    # End of rising edge
        rhi = np.nonzero(trace[:ipeak]<=0.2*peak*np.e)[0][-1]
        riseGuess = (bins[rhi]-bins[rlo])
    
        # Fall time: look for two e-foldings on falling side
        flo = np.nonzero(trace[ipeak:]<=0.8*peak)[0][0]     # Start of falling tail
        fhi = np.nonzero(trace[ipeak:]<=0.4*peak/np.e)[0][0]
        fallGuess = (bins[fhi]-bins[flo])/2

        # Analytic peak position is where d/dt of pulse shape is zero
        # ==> t_peak = tR * ln[(tF+tR)/tR]
        tpeak = riseGuess * np.log((fallGuess+riseGuess)/riseGuess)
        offsetGuess = bins[ipeak] - tpeak
    
        # Scale factor should be max of shape scaled by actual peak value
        pmax = self.TESshape(tpeak, 1., riseGuess, fallGuess, 0.)
        scaleGuess = peak / pmax

        self.printVerbose(f"guessTES: scale {scaleGuess:.4e} rise {riseGuess:.4e},"
                          f" fall {fallGuess:.4e}, offset {offsetGuess:.4e} us")
    
        return scaleGuess, riseGuess, fallGuess, offsetGuess

    def guessFET(self, bins, trace):
        """Returns initial guesses for FET fit inverse decay and recovery times."""
    
        peak = trace.max()
        ipeak = trace.argmax()
        self.printVerbose(f"guessFET: peak {peak} @ bin {ipeak} (t {bins[ipeak]})")

        # Peak should be at t=+binWidth (first bin after trigger)
        istart = np.nonzero(bins>=0.)[0][0]+1
        offsetGuess = ipeak - istart
    
        # Decay time: look for second e-folding after the peak
        dhi = np.nonzero(trace[ipeak:]>=peak/(2.*np.e))[0][-1]
        decayGuess = 2./(bins[dhi]-bins[0])
    
        # Recovery time; look for second e-folding beyond minimum
        tmin = trace.min()
        imin = trace.argmin()
        self.printVerbose(f" bottom {tmin} @ bin {imin} (t {bins[imin]})")

        recoveryGuess = 0.
        if tmin < 0:
            tlast = trace[imin:].max()
            rlo = np.nonzero(trace[imin:]>=tmin*0.8)[0][0]
            rhi = np.nonzero(trace[imin:]>=min(tmin*0.4/np.e,tlast))[0][0]
            recoveryGuess = 2./(bins[rhi]-bins[rlo])
            if recoveryGuess < 0.1*decayGuess:
                self.printVerbose(f" recoveryGuess {recoveryGuess} not physical.")
                recoveryGuess = 0.

        # FET function is [A/(D-R)]*(D*exp(-t*D) - R*exp(-t*R))
        scaleGuess = peak / (decayGuess-recoveryGuess)

        self.printVerbose(f"guessTES: scale {scaleGuess:.4e} decay {decayGuess:.4e} /us,"
                          f" recovery {recoveryGuess:.4e} /us, offset {offsetGuess:.4e} us")
    
        return scaleGuess, decayGuess, recoveryGuess, offsetGuess


    def guessRange(self, guess=None):
        """Compute allowed parameter ranges for fit based on initial guess
           values."""
        self.printVerbose(f"guessRange(guess={guess})")

        if guess is None or guess==0.:
            return (-np.inf, np.inf)

        lower = 0.1*np.array(guess)
        upper = 5.*np.array(guess)
        bounds = (lower, upper)

        for i,g in enumerate(guess):
            if guess[i] == 0.:       	# Bounds can't both be zero!
                bounds[0][i] = -np.inf
                bounds[1][i] = np.inf
            if guess[i] < 0.:		# Negative bounds need to be swapped
                bounds[0][i],bounds[1][i] = bounds[1][i],bounds[0][i]

        return bounds


    def fittingRange(self, trace, cut=0.2):
        """Return starting and ending points for pulse fit, corresponding to
           'cut' height on either side of peak.  Assumes TES trace is flipped
           and baseline-subtracted, or FET trace is charge-flipped."""
        peak = max(trace)          # Peak Height
        ipeak = trace.tolist().index(peak)
        self.printVerbose(f"fittingRange: peak {peak} @ bin {ipeak}")

        ilo = 0
        ihi = len(trace)
    
        if self.sensor=="TES":
            ilo = np.nonzero(trace[:ipeak]<=cut*peak)[0][-1]          # End of rising edge
            ihi = ipeak+np.nonzero(trace[ipeak:]<=cut*peak)[0][0]     # Start of falling tail
            self.printVerbose(f"fittingRange: TES I>{cut}*peak [{ilo}:{ihi}]")
        elif self.sensor=="FET":
            ilo = trace.argmax()+1	# Start fit just after peak bin
            ihi = ilo+2000		# TODO: Find far end after recovery
            self.printVerbose(f"fittingRange: FET [{ilo}:{ihi}]")

        return ilo, ihi


    ### General fitting function: sensor-specific info is in 'pulseShape' and 'guessFunc'

    def fitTrace(self, bins, trace, pulseShape, guessFunc=None, dobounds=True):
        """Fits input trace with binning to TES or FET shape; using function
           for initial values
           Output: a      = scale factor from fit
                   t1     = rise time for TES, or decay rate for FET
                   t2     = fall time for TES, or recovery rate for FET
                   offset = t0 of best fit relative to t=0 bin
        """
        self.printVerbose(f"fitTrace(bins, trace, pulseShape={pulseShape},"
                          f" guessFunc={guessFunc}, dobounds={dobounds})")
    
        start, end = self.fittingRange(trace)

        guess = guessFunc(bins, trace) if guessFunc else None
        bounds = self.guessRange(guess) if dobounds else (-np.inf,np.inf)
    
        self.printVerbose(f" range [{start}:{end}]\n guess {guess}")
        if (dobounds): self.printVerbose(f"bounds {bounds}")
                
        params, _ = curve_fit(pulseShape, bins[start:end], trace[start:end],
                              p0=guess, bounds=bounds)
        self.printVerbose(f" final result {params}")

        return params


    def plot(self, detname, channel, bins, trace, fitshape):
        """Generate linear and log overlays of trace and fitted function."""
        self.printVerbose(f"plot(detname='{detname}', channel={channel},"
                          f" bins, trace, fitshape)")
    
        titleName = detname if detname else "Trace"

        template = self.template(detname, channel)
    
        self.overlay(titleName, bins, trace, fitshape, template)
        plt.savefig(f"{titleName}-{self.sensor}_traceFit.eps", format='eps')
        plt.savefig(f"{titleName}-{self.sensor}_traceFit.png")

    
    def overlay(self, detname, bins, trace, fitshape, template):
        """Plots TES or FET trace (log and linear) with specified binning, overlaid
           with fitted shape and template detname argument used for plot title."""
        self.printVerbose(f"overlay(detname='{detname}', bins, trace,"
                          f" fitshape, template)")
    
        units  = { "TES": "\mu A",
                   "FET": "mV" }
        xlim   = { "TES": [ [max(-100,bins.min()),min(3000,bins.max())],
                            [max(-10,bins.min()),min(300,bins.max())] ],
                   "FET": [ [max(-100,bins.min()),min(1000,bins.max())],
                            [max(-10,bins.min()),min(300,bins.max())] ] }
        yscale = { "TES": ["log","linear"],
                   "FET": ["linear","linear"] }

        start, end = self.fittingRange(trace)
        fig, axes = plt.subplots(1, 2, figsize=(12*0.7, 4*0.7), dpi=200)
        for plot in range(2):
            currentAxis = axes.flatten()[plot]

            if template is not None:
                currentAxis.plot(bins,template*max(trace),lw=1,ls='--',
                                 color='black', label='Template')

            currentAxis.plot(bins, trace,lw=1,ls='-',color='red',label='Simulation')
            currentAxis.plot(bins[start:end], fitshape[start:end], label='Fit')
        
            currentAxis.set_xlabel(r"Time [$\mathrm{\mu s}$]")
            currentAxis.set_ylabel(r"Amplitude [$\mathrm{"+units[self.sensor]+"}$]")
            currentAxis.legend()
            currentAxis.set_xlim(xlim[self.sensor][plot])
            currentAxis.set_yscale(yscale[self.sensor][plot])

        if detname: plt.title(detname)
        plt.tight_layout()


    def template(self, detname, chan):
        """Extract channel template for specified detector, as Numpy array."""
        self.printVerbose(f"template(detname='{detname}', chan={chan})")
    
        if detname is None or detname == "": return None
    
        templatePath = f"{CDMS_SUPERSIM}/CDMSgeometry/data/dmc/{detname}/{self.sensor}Templates"
        if not os.path.isfile(templatePath): return None

        colname = {"TES":"Traces", "FET":"Trace"}[self.sensor]
        templates = pd.read_csv(templatePath,sep="\t")
        template = templates.loc[chan,colname].split()
        template = np.array([float(i) for i in template])

        self.printVerbose(f" got template from {templatePath} with"
                          f" {len(template)} bins")
        return template


### MAIN PROGRAM ###

def main():
    args = getargs()

    # Split full argument list to avoid passing unrecognized parameters
    ctor = { k:args[k] for k in ['file','sensor','verbose'] }
    fitter = traceFitter(**ctor)

    fitargs = { k:args[k] for k in args.keys()-ctor.keys() }
    fitter.doFit(**fitargs)

def getargs():
    """Returns arguments from the command line as a dictionary, for easier use.
       Output: settings = { 'file':    <name-of-DMC-file>,
                             'detname': <detector type>, from -d>,
                             'event':   <event number>, from -e>,
                             'channel': <channel number>, from -c>,
                             'sensor':  <sensor type>, from -s>,
                             'doplot':  <True|False>, from -p>,
                             'verbose': <True|False>, from -v>,
                           }
    """
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'c:d:e:s:hpv')
    except getopt.GetoptError as err:
        sys.exit(2)

    if args is None or len(args)==0:
        usage()
        sys.exit(2)
        
    settings = {'file':    args[0],   # Filename is required
                'detname': "",        # Detector name not required
                'event':   0,         # First event, first channel
                'channel': 0,
                'sensor':  "TES",     # TES or FET
                'doplot':  False,     # Results only, no figures
                'verbose': False,
               }

    for o,a in opts:
        if   o == '-c':
            settings['channel'] = int(a)
        elif o == '-d':
            settings['detname'] = a
        elif o == '-e':
            settings['event'] = int(a)
        elif o == '-h':
            usage()
            sys.exit(0)
        elif o == '-p':
            settings['doplot'] = True
        elif o == '-s':
            settings['sensor'] = a
        elif o == '-v':
            settings['verbose'] = True

    if settings['verbose']: print(f"settings = {settings}")
        
    return settings


### COMMAND-LINE CALL ###

if __name__ == "__main__":
    main()
    
### END OF FILE ###
