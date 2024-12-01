#!/bin/env python3

"""Usage: import phonon_eff

   Defines useful functions for calculating and plotting phonon collection
   efficiency from single detector or full geometry DMC simulations.  See
   `dir(phonon_eff)` for a list of available functions.
"""

# 20240719  Michael Kelsey
# 20240723  Add function to draw pre-defined SimFiducial regions

from imports import *
from cats.cdataframe import CDataFrame

def load(files):
    """Loads geometry, events, and hits from specified file or group of files.
       Returns three Numpy dictionaries (i.e., RDF.AsNumpy())."""
    geomBranches = ['DetNum','DetType','Crystal','Voltage']
    geom = CDataFrame("G4SettingsInfoDir/Geometry", files).AsNumpy(geomBranches)

    evBranches = ['EventNum','DetNum','DetType','PhononE']
    events = CDataFrame("G4SimDir/g4dmcEvent", files).Filter("!Empty").AsNumpy(evBranches)

    # NOTE: Don't "lukeEnergyEst" branch, appears to come out zero (despite fixing voltage registration)
    hitBranches = ['EventNum','DetNum','trueEnergy','truedEdx','trueNIEL','truePairs','phononEnergy']
    hitBranches += ['Xdet','Ydet','Zdet']
    hits = CDataFrame("G4SimDir/g4dmcHits", files).AsNumpy(hitBranches)
    hits["Rdet"] = np.sqrt(hits["Xdet"]**2 + hits["Ydet"]**2)
    
    # Add derived quantities using functions defined below
    hitExpected(geom, hits)
    eventExpected(events, hits)
    eventHitPos(events, hits)
    efficiency(events)
    
    return geom, events, hits

def voltage(geom, detNum):
    """Retrieve detector voltage from Geometry data (must be loaded and passed)."""
    return geom["Voltage"][geom["DetNum"]==detNum][0]

def hitExpected(geom, hits):
    """Compute expected energy for each individual hit ('lukeEnergyEst' is broken)."""
    hits["NTLenergy"] = hits["truePairs"]*np.vectorize(voltage)(geom, hits["DetNum"])
    hits["Eexpected"] = hits["trueEnergy"]+hits["NTLenergy"]
    return

def eventExpected(events, hits):
    """Compute expected energy for each event, summing associated hits.
       hitExpected() must have already been called."""

    # Duplicate DetNum structure, which has one entry per event+detector
    events["Eexpected"] = np.zeros_like(events["DetNum"])

    # NOTE: Separate masks for different N-tuples, to identify correct entries
    for ev in np.unique(hits["EventNum"]):
        for det in np.unique(hits["DetNum"][hits["EventNum"]==ev]):
            mask = (hits["EventNum"]==ev)&(hits["DetNum"]==det)
            dest = (events["EventNum"]==ev)&(events["DetNum"]==det)
            events["Eexpected"][dest] = sum(hits["Eexpected"][mask])

    return

def eventHitPos(events, hits):
    """Compute 'average' hit position and spread. summing associated hits."""

    # Duplicate DetNum structure, which has one entry per event+detector
    newcols = ["Xavg","Yavg","Zavg"]
    for col in newcols:
        events[col] = np.zeros_like(events["DetNum"])

    # NOTE: Separate masks for different N-tuples, to identify correct entries
    for ev in np.unique(hits["EventNum"]):
        for det in np.unique(hits["DetNum"][hits["EventNum"]==ev]):
            mask = (hits["EventNum"]==ev)&(hits["DetNum"]==det)
            dest = (events["EventNum"]==ev)&(events["DetNum"]==det)

            events["Xavg"][dest] = np.average(hits["Xdet"][mask])
            events["Yavg"][dest] = np.average(hits["Ydet"][mask])
            events["Zavg"][dest] = np.average(hits["Zdet"][mask])

    events["Ravg"] = np.sqrt(events["Xavg"]**2+events["Yavg"]**2)
    return

def efficiency(events):
    """Compute phonon collection efficiency per event.  eventExpected() must
       have already been called."""
    events["sumPhE"] = events["PhononE"].sum(axis=1)
    events["PhEff"] = events["sumPhE"]/events["Eexpected"]
    return

def plotEfficiency(files, title="Phonon Collection Efficiency", figfile=None, detType=None,
                   xlim=None):
    """Generate histogram of phonon collection efficiency."""
    geom,events,hits = load(files)

    # Generate selection mask if detector type was specified
    mask = (events["DetType"]==detType) if detType else np.full_like(events["PhEff"],True,dtype=bool)
    
    plt.hist(events["PhEff"][mask],bins=100,range=xlim)
    plt.title(title)
    plt.xlabel("Phonon Efficiency")
    
    if figfile: plt.savefig(figfile)
    return

def stackEfficiency(files, title="SNOLAB Phonon Collection Efficiency", figfile=None,
                    xlim=None):
    """Generate stacked histogram of phonon collection efficiency for SNOLAB
       detector types."""
    geom,events,hits = load(files)

    detCode = {"iZIP7":700, "iZIP7Si":701, "HV100mm":710, "HV100mmSi":711}
    effByDet = {}
    for name,det in detCode.items():
        effByDet[name] = events["PhEff"][events["DetType"]==det]

    plt.hist(effByDet.values(),label=effByDet.keys(),bins=100,range=xlim,histtype='bar',stacked=True)
    plt.legend()
    plt.title(title)
    plt.xlabel("Phonon Efficiency")
        
    if figfile: plt.savefig(figfile)
    return

def mapEfficiency(files, detType, title=None, usefig=None, izipcut=0.99, hvcut=0.99):
    """Generate R-Z heatmaps showing the bad efficiency region for specified
       SNOLAB detector type.  If 'title' is specified, plot will include that
       title along with axis labels.  Otherwise, it's assumed to be part of a
       larger presentation."""
    detCodes = {700:"iZIP7", 701:"iZIP7Si", 710:"HV100mm", 711:"HV100mmSi"}
    detName = detCodes[detType]
    
    geom,events,hits = load(files)
    lowEff = hvcut if (detType>=710) else izipcut
    sel = (events["DetType"]==detType)&(events["PhEff"]<lowEff)

    if usefig is None: usefig = plt.gca()

    usefig.scatter(events["Ravg"][sel]*1e3,events["Zavg"][sel]*1e3,
                   c=events["PhEff"][sel],vmin=0.,vmax=1.1,cmap='viridis',
                   label=detName)
    usefig.set_xlim(0., 52.)
    usefig.set_ylim(-18., 18.)
    usefig.annotate(detName, xy=(3.,0.))
    showSimFid(detType, usefig)
    
    if title:
        usefig.set_title(title)
        usefig.set_xlabel("R [mm]")
        usefig.set_ylabel("Z [mm]")
    return

def mapAllEfficiencies(files, title="SNOLAB SimFiducial Regions", figfile=None,
                       izipcut=0.99, hvcut=0.99):
    """Generate set of four R-Z heatmaps showing the bad efficiency region
       for each SNOLAB detector type."""
    detCodes = {700:"iZIP7", 701:"iZIP7Si", 710:"HV100mm", 711:"HV100mmSi"}

    fig, ax = plt.subplots(2,2,figsize=(9,6), dpi=300)
    fig.set_tight_layout(True)
    fig.suptitle(title)
    
    for iplot in range(4):
        thePlot = ax.flatten()[iplot]
        dType = list(detCodes.keys())[iplot]
        mapEfficiency(files, dType, usefig=thePlot, izipcut=izipcut, hvcut=hvcut)
                
    if (figfile): plt.savefig(figfile)
    return

def stackEfficiencySets(filesets, title="SNOLAB Phonon Collection Efficiency", figfile=None):
    """Generate stacked histogram of phonon collection efficiency for SNOLAB
       detector types, assuming that 'filesets' is a dictionary in the form
       of { detType:files } pairs.  The numeric detType should be used for the
       keys, so that it can be passed to mapEfficiency() without conversion."""
    detCodes = {700:"iZIP7", 701:"iZIP7Si", 710:"HV100mm", 711:"HV100mmSi"}

    effByDet = {}
    for det, files in filesets.items():
        geom,events,hits = load(files)
        effByDet[detCodes[det]] = events["PhEff"][events["DetType"]==det]

    plt.hist(effByDet.values(),label=effByDet.keys(),bins=100,histtype='bar',stacked=True)
    plt.legend()
    plt.title(title)
    plt.xlabel("Phonon Efficiency")

    if figfile: plt.savefig(figfile)
    return


def mapEfficiencySets(filesets, title="SNOLAB SimFiducial Regions", figfile=None,
                       izipcut=0.99, hvcut=0.99):
    """Generate set of four R-Z heatmaps assuming that 'filesets' is a dictionary
       in the form of { detType:files } pairs.  The numeric detType should be used
       for the keys, so that it can be passed to mapEfficiency() without conversion."""
    fig, ax = plt.subplots(2,2,figsize=(9,6), dpi=300)
    fig.set_tight_layout(True)
    fig.suptitle(title)
    
    for iplot in range(4):
        thePlot = ax.flatten()[iplot]
        dType = list(filesets.keys())[iplot]
        files = filesets[dType]
        mapEfficiency(files, dType, usefig=thePlot, izipcut=izipcut, hvcut=hvcut)
                
    if (figfile): plt.savefig(figfile)
    return

def showSimFid(dettype, rzplot):
    """Draw SimFiducial boundaries for given detector in current plot.
       Boundary points determined by McKay Osborne for SNOLAB detectors,
       https://docs.google.com/presentation/d/1wNvbBrwvAZE4pgOmLAUJpYLo0MN_hKlQd5F4AZCioFA
       https://docs.google.com/presentation/d/1EHVi-MykxzPVz7w-hcqLK7hpJiTKG1yru4hrfmdoEUs/edit#slide=id.g2667a5c8d84_0_32
       using uniform electric field (not EPot files).

       NOTE: Coordinates are given in millimeters, not meters."""
    simfid = { 700: { "Rupper": 46.8, "Zupper": 16.4, "Rlower": 29.0, "Zlower": -16.4 },
               701: { "Rupper": 45.8, "Zupper": 16.4, "Rlower": 47.0, "Zlower": -16.4 },
               710: { "Rupper": 48.0, "Zupper": 16.4, "Rlower": 48.0, "Zlower": -16.4 },
               711: { "Rupper": 49.0, "Zupper": 16.4, "Rlower": 48.0, "Zlower": -16.4 },
             }
    if dettype not in simfid: return
    xpoints = [ 0., simfid[dettype]["Rupper"], simfid[dettype]["Rlower"], 0.]
    ypoints = [ simfid[dettype]["Zupper"], simfid[dettype]["Zupper"],
                simfid[dettype]["Zlower"],simfid[dettype]["Zlower"] ]
    rzplot.plot(xpoints, ypoints, "--", c="maroon")
    return