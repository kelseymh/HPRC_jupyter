#!/usr/bin/env python
#
# Usage: charge_bounces.py
#
# Implemented in place of a Jupyter notebook due to Grace slowness.  Edit
# the dataset references below to select which DMC files to use.
#
# 20250925  Michael Kelsey

import trajectoryPlotter as tplot
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (5,3.75)
plt.rcParams['figure.autolayout'] = True

# Trajectories from testing G4CMP chargeBounces setting
datadir = "data/chargeBounces"
bounce0   = datadir+"/bounce0_51250925_0000.root"
bounce10  = datadir+"/bounce10_51250925_0000.root"
bounce100 = datadir+"/bounce100_51250925_0000.root"

# Make plot of desired trajectories

def trajectories(files,prefix="bounces",title="All tracks",zoom=None):
    """Make and save plot of all trajectories in given dataset; 'name'
       is used for output filename, 'title' is used for plotting.
       Optional zoom should be a list of lists: ((xmin,xmax),(ymin,ymax))."""
    print("Generating",prefix,"YZ trajectories",flush=True)
    
    # Zero bounce case, whole detector
    plt.figure()
    plt.title(title)
    plt.xlabel("Y [mm]")
    plt.ylabel("Z [mm]")

    traj = tplot.getTracks(files)
    tplot.drawAllYZ(traj)
    if zoom is not None:
        plt.xlim(zoom[0])
        plt.ylim(zoom[1])
        
    tplot.trackLegend()
    plt.savefig(datadir+"/"+prefix+"_YZ.png")


# Generate all three figures

tplot.debug = True
trajectories(bounce0,"bounce0","CDMSlite Zero bounces (immediate absorption)")
trajectories(bounce10, "bounce10", "CDMSlite 10 bounces before absorption")
trajectories(bounce100,"bounce100","CDMSlite 100 bounces before absorption")
