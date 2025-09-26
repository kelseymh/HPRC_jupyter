#!/usr/bin/env python
#
# Usage: import trajectoryPlotter
#
# Provides functions to access mcHitCounter or mcTrajectory data and draw
# full track trajectories.
#
# 20250925  Michael Kelsey

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from cats.cdataframe import CDataFrame


# Load trajectories for given file(s)
def getTracks(files):
    """Load mcTrajectory TTree from specified file or files."""
    if not files or len(files)==0: return None
        
    # Older datasets use mcHitCounter instead of mcTrajectory
    try:
        trajCDF = CDataFrame("G4SimDir/mcTrajectory", files)
    except:
        trajCDF = CDataFrame("G4SimDir/mcHitCounter", files)

    branches = ["EventNum","Track","Step","X1","Y1","Z1","X3","Y3","Z3","KE","Yield","PName","Charge"]
    steps = trajCDF.AsNumpy(branches)
    fixStepOverflow(steps)
    convertToMM(steps)

    steps["R1"] = np.sqrt(steps["X1"]**2+steps["Y1"]**2)
    steps["R3"] = np.sqrt(steps["X3"]**2+steps["Y3"]**2)
    
    return steps

def fixStepOverflow(steps):
    """Correct for integer overflow of Geant4 CurrentStep counter."""
    steps["Step"][(steps<0)] += 2**32
    return steps

def convertToMM(steps):
    """Convert from CDMS units of meters to mm, for better plotting."""
    for xyz in ["X1","Y1","Z1","X3","Y3","Z3"]: steps[xyz] *= 1e3
    return steps

# Get list of unique tracks in event
def listTraj(steps):
    """Extract table of (event,track) pairs from trajectory/step data."""
    evTrk = []

    events = np.unique(steps["EventNum"])
    for ev in events:
        isEv = steps["EventNum"]==ev
        evTrk.append((ev,trk) for trk in np.unique(steps["Track"][isEv]))

    return evTrk

def trajCut(steps, event, track):
    """Define a Numpy index mask to select the desired trajectory steps."""
    return (steps["EventNum"]==event) & (steps["Track"]==track)

# Build plottable trajectory for given track
def trajectoryXYZ(steps, event, track):
    """Construct three arrays of X, Y and Z coordinates for track."""
    theTrack = trajCut(event,track)
    step0 = steps["Step"]==(steps["Step"][theTrack].min())

    trajX = [ steps["X1"][theTrack & step0] ]
    trajY = [ steps["Y1"][theTrack & step0] ]
    trajZ = [ steps["Z1"][theTrack & step0] ]
    for i in steps["Step"][theTrack]:
        stepi = steps["Step"]==i
        trajX.append(steps["X3"][theTrack & stepi])
        trajY.append(steps["Y3"][theTrack & stepi])
        trajZ.append(steps["Z3"][theTrack & stepi])

    return trajX,trajY,trajZ

def trajectoryRZ(steps, event, track):
    """Construct two arrays of R,Z coordinates for track."""
    theTrack = trajCut(event,track)
    step0 = steps["Step"]==(steps["Step"][theTrack].min())

    trajR = [ steps["R1"][theTrack & step0] ]
    trajZ = [ steps["Z1"][theTrack & step0] ]
    for i in steps["Step"][theTrack]:
        stepi = steps["Step"]==i
        trajR.append(steps["R3"][theTrack & stepi])
        trajZ.append(steps["Z3"][theTrack & stepi])

    return trajR,trajZ


# Utilities for mapping track types to colors in each plot

def chgColors():
    """Return a colormap for orange e-, green phonon, cyan h+.  To use
       this, use 'c=steps["Charge"]+1, cmap=chgColors()', in plot."""
    return colors.ListedColormap(['o','g','c'])

def trackLegend():
    """Generate a custom legend for the three trajectory colors."""
    cmap = chgColors()
    types = [Line2D([0],[0], lw=4, color=cmap(0), label="Electron"),
             Line2D([0],[0], lw=4, color=cmap(1), label="Phonon"),
             Line2D([0],[0], lw=4, color=cmap(2), label="Hole")]
    plt.legend(handles=types)


# Draw requested track in current figure (will create new one if needed

def drawXYZ(steps, event, track):
    """Plot 3D trajectory in current figure."""
    x,y,z = trajectoryXYZ(steps, event, track)
    plt.scatter(x,y,z,c=steps["Charge"][+1,cmap=chgColors())

def drawXZ(steps, event, track):
    """Plot 2D trajectory projected onto X-Z plane."""
    x,y,z = trajectoryXYZ(steps, event, track)
    plt.scatter(x,z,c=steps["Charge"]+1,cmap=chgColors())

def drawYZ(steps, event, track):
    """Plot 2D trajectory projected onto Y-Z plane."""
    x,y,z = trajectoryXYZ(steps, event, track)
    plt.scatter(y,z,c=steps["Charge"]+1,cmap=chgColors())

def drawAllXZ(steps):
    """Plot all of the trajectories in the dataset in X-Z plane."""
    for evTrk in listTraj(steps):
        drawXZ(steps, *evTrk)

def drawAllYZ(steps):
    """Plot all of the trajectories in the dataset in Y-Z plane."""
    for evTrk in listTraj(steps):
        drawYZ(steps, *evTrk)
