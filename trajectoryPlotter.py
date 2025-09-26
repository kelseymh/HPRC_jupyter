#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from imports import *
get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


# Phonon trajectory datasets
datadir = "data/phononTracks/"
walk0mm  = datadir+"phonon-1meV-0mm_51250806_0000.root"
walk34mm = datadir+"phonon-1meV-34mm_51250806_0000.root"
walkbulk = datadir+"phonon-1meV-bulk_51250807_0000.root"
ball0mm  = datadir+"phonon-0.1meV-0mm_51250806_0000.root"
ball34mm = datadir+"phonon-0.1meV-34mm_51250806_0000.root"
ballbulk = datadir+"phonon-0.1meV-bulk_51250807_0000.root"


# In[ ]:


# Load trajectories for given file
def getTracks(file):
    """Load mcTrajectory TTree from specified file or files."""
    if not file: return None
        
    branches = ["EventNum","Track","Step","X1","Y1","Z1","X3","Y3","Z3","KE","Yield","PName"]
    steps = CDataFrame("G4SimDir/mcTrajectory", file).AsNumpy(branches)
    steps["R1"] = np.sqrt(steps["X1"]**2+steps["Y1"]**2)
    steps["R3"] = np.sqrt(steps["X3"]**2+steps["Y3"]**2)
    
    return steps


# In[ ]:


# Get list of unique tracks in event
def listTraj(steps):
    """Extract table of (event,track) pairs from trajectory/step data."""
    evTrk = []

    events = np.unique(steps["EventNum"])
    for ev in events:
        isEv = steps["EventNum"]==ev
        evTrk.append((ev,trk) for trk in np.unique(steps["Track"][isEv]))

    return evTrk


# In[ ]:


# Build plottable trajectory for given track
def trajectoryXYZ(steps, event, track):
    """Construct three arrays of X, Y and Z coordinates for track."""
    theTrack = (steps["EventNum"]==event) & (steps["Track"]==track)
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
    theTrack = (steps["EventNum"]==event) & (steps["Track"]==track)
    step0 = steps["Step"]==(steps["Step"][theTrack].min())

    trajR = [ steps["R1"][theTrack & step0] ]
    trajZ = [ steps["Z1"][theTrack & step0] ]
    for i in steps["Step"][theTrack]:
        stepi = steps["Step"]==i
        trajR.append(steps["R3"][theTrack & stepi])
        trajZ.append(steps["Z3"][theTrack & stepi])

    return trajR,trajZ


# In[ ]:


bouncing = getTracks(ball0mm)
x,y,z = trajectoryXYZ(bouncing,0,0)
plt.scatter(x,y,z)


# In[ ]:




