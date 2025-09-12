#!/usr/bin/env python

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from cats.cdataframe import CDataFrame
import math
import os, os.path
import glob
import re

class chargeSteps:
    """Driver class to load data and make plots of trajectory
       statistics for charge track bounding.

       Constructor arguments:
       ddir:  Directory containing ROOT files at various voltages
       dtype: Detector type name (HV100mm, iZIP5, etc.)
       infix (optional): Matching string for filenames
       volts (optional): List of voltages to be read in (filename matches)

       See help(chargeSteps) for a complete list of member functions.
    """

    ### Constructor and supporting functions
    def __init__(self, ddir, dtype, infix='', volts=None):
        """chargeSteps constructor."""
        self.datadir = ddir
        self.det = dtype
        self.voltage = volts
        self.__loadDatasets(infix)
        return

    def __loadDatasets(self, infix=''):
        """Collect single-track ROOT files for selected detector type.
           Optional `infix` selects particular named files.  Optional
           `volts` is a list of requested voltage settings; if not given,
           all matching files are used, and the voltage list generated."""

        self.data = { "U": self.__useFiles("U",infix),
                      "E": self.__useFiles("E",infix)
                    }

        # Collect all of the keys from both lists
        if self.voltage is None:
            self.voltage = sorted(set([*self.data["U"]]+[*self.data["E"]]))
            
        self.vmin = min(self.voltage)
        self.vmax = max(self.voltage)
        return

    def __useFiles(self, vtype, infix):
        """Generate list of files and associated voltages for given voltage
           type, either 'U' (uniform) or 'E' (EPot files).  Globbing infix
           and list of voltages are passed from constructor."""

        flist = {}
        allf = glob.glob(f'{self.datadir}/{self.det}-*V{vtype}{infix}_5125*_0*.root')
            
        if len(allf)>0:
            if self.voltage is not None: 
                allv = self.voltage
            else:
                allv = sorted([int(re.search(r'-(\d+)V.',f).groups()[0]) for f in allf])

            flist = {v:glob.glob(f'{self.datadir}/{self.det}-{v}V{vtype}{infix}_5125*_0*.root')
                     for v in allv}

        # Remove empty entries from flist for clarity
        flist = {v:f for v,f in flist.items() if len(f)>0}
        return flist

    def __flatDataset(self, vtype):
        """Returns a flattened list of all the files for all voltages."""
        return sum(self.dataset[vtype].values()) if self.dataset[vtype] else None

    
    # Report configuration
    def __repr__(self):
        return f"<chargeSteps('{self.datadir}','{self.det}')>"

    def __str__(self):
        import json
        return f"""{self.__repr__()}
Voltage: {self.voltage}
Uniform: {json.dumps(self.data["U"],indent=4)}
EPot:    {json.dumps(self.data["E"],indent=4)}
"""

    # Apply restricted voltage range for plotting
    def limitVoltage(self, vmin=None,vmax=None,range=None):
        """Restrict plotting to smaller set of voltages.  User may provide
           either minimum and maximum values, or range=(min,max).  Returns
           subset of self.voltage list within the given range."""
        
        if vmin is None: vmin = self.vmin
        if vmax is None: vmax = self.vmax
        if range is None and (vmin is not None and vmax is not None):
            range = (vmin,vmax)

        npv = np.array(self.voltage)
        return npv[(vmin<=npv) & (npv<=vmax)]
        
    # Set up N-tuple with individual step information
    def getSteps(self, files):
        """Load mcHitCounter with individual steps, to get trajectory length."""
        if not files: return None

        branches = ["EventNum","Track","Step","StepLen","Charge","PName"]
        steps = CDataFrame("G4SimDir/mcHitCounter", files).AsNumpy(branches)
        return steps

    # Set up N-tuple with computed path length quantities
    def getHits(self, files):
        """Load mczip0 for DMC tracks, with track start/end information."""
        if not files: return None
        
        branches = ["EventNum","Track","Step","StepLen","X1","Y1","Z1","X3","Y3","Z3","Charge","PName"]
        hits = CDataFrame("G4SimDir/mczip0", files).AsNumpy(branches)
        hits["Flight"] = np.sqrt((hits["X3"]-hits["X1"])**2+(hits["Y3"]-hits["Y1"])**2+(hits["Z3"]-hits["Z1"])**2)
        hits["R1"] = np.sqrt(hits["X1"]**2+hits["Y1"]**2)
        hits["R3"] = np.sqrt(hits["X3"]**2+hits["Y3"]**2)
        return hits

    def addStepSum(self, hits, files):
        """Use step-counter data to get trajectory length, replace Hits["StepLen"]."""
        if not hits or not files: return hits
        
        # Get cumulative path length from step counter (not done properly in mczip)
        steps = getSteps(files)
        for ev in hits["EventNum"]:
            for trk in hits["Track"]:
            # Separate indices (cuts) for each N-tuple
                ihit = ((hits["EventNum"]==ev) & (hits["Track"]==trk))
                istep = ((steps["EventNum"]==ev) & (steps["Track"]==trk))
                hits["StepLen"][ihit] = steps["StepLen"][istep].sum()
            
        return hits

    def getTrajectory(self, files, event=None, track=None):
        """Load mcTrajectory or mcHitCounter for DMC charges, with all steps."""
        if not files: return None

        # Older datasets use mcHitCounter instead of mcTrajectory
        try:
            trajCDF = CDataFrame("G4SimDir/mcTrajectory", files)
        except:
            trajCDF = CDataFrame("G4SimDir/mcHitCounter", files)

        filter = "true"
        if event: filter = f"EventNum=={event}"
        if track:
            if len(filter)>0: filter += f" && "
            filter += f" Track=={track}"
        
        branches=["EventNum","Track","Step","StepLen","X1","Y1","Z1","X3","Y3","Z3","Charge","PName"]
        traj = trajCDF.Filter(filter).AsNumpy(branches)
        traj["R1"] = np.sqrt(traj["X1"]**2+traj["Y1"]**2)
        traj["R3"] = np.sqrt(traj["X3"]**2+traj["Y3"]**2)
        return traj
        
    def getGeometry(self, files=None):
        """Load basic geometry information from first file in set.
           NOTE: Assumes all files are for single detector."""
        if files is None:
            if   len(self.data["U"])>0: files = self.flattenList(self.data["U"].values())
            elif len(self.data["E"])>0: files = self.flattenList(self.data["E"].values())
            else: return None

        # Only load minimal information for now (gross dimensions)
        branches = ["DetType","TypeName","Material","Radius","Height","Voltage"]
        geom = CDataFrame(f"G4SettingsInfoDir/Geometry", files[0]).AsNumpy(branches)
        return geom

    def getJobInfo(self, files=None):
        """Load job information (events, time) from specified files."""
        if files is None: return None
            
        branches = ["Events","Requested","Elapsed","Threads","UserCPU"]
        jobinfo = CDataFrame("G4SettingsInfoDir/runtime", files).AsNumpy(branches)
        return jobinfo
            
    @staticmethod
    def flattenList(files):
        """Convert list-of-list from dataset to flat list for chaining."""
        flat = []
        if len(files)>0:
            for f in files: flat += f
        return flat

    ### Formatting strings for clean plots and filenames
    @staticmethod
    def TitleField(vtype):
        """Convert field type (U or E) into readable name."""
        if vtype=="U": return "Uniform Field"
        if vtype=="E": return "EPot Files"
        return ""

    @staticmethod
    def ShortField(vtype):
        """Convert field type (U or E) into readable name."""
        if vtype=="U": return "Uniform"
        if vtype=="E": return "EPot"
        return ""

    @staticmethod
    def ValUnits(stepcut):
        """Return a useful string for limiting steps in k or M."""
        if stepcut<1e6:  return f"{stepcut/1000:.0f}k"
        if stepcut<1e9:  return f"{stepcut/1e6:.0f}M"
        if stepcut<1e12: return f"{stepcut/1e9:.0f}G"
        if stepcut<1e15: return f"{stepcut/1e12:.0f}T"
        if stepcut<1e99: return f"{stepcut:.0f}"
        return ""
        
    def Title(self, vtype, volt=None):
        """Construct consistent plot title string for detector type and bias."""
        flong = self.TitleField(vtype)
        if volt is not None: flong = f'{volt}V {flong}'
        return ' '.join((self.det,flong))

    def Filename(self, name, vtype="", volt=None):
        """Construct consistent filenames with detector type and bias."""
        fname = f'{name}_{self.det}'
        if volt is not None: fname += f'_{volt}V'
        if vtype != "": fname += "-"+self.ShortField(vtype)
        return os.path.join(self.datadir,fname+".png")

    def DetectorRZ(self, file):
        """Draw detector boundary lines to get full axis range."""
        if not file: return

        geom = self.getGeometry(file)
        rmax = 1e3*geom["Radius"][0]
        zmax = 1e3*geom["Height"][0]/2.
        plt.plot((0,rmax,rmax,0),(-zmax,-zmax,zmax,zmax),'--k',linewidth=1)

        geom = None        # Ensure that allocated memory is freed
        return

    ### Generate trajectories in R-Z plane
    def PlotTrajectory(self, vtype, volt, event=None, track=None,
                       rlim=None, zlim=None):
        """Draw trajectory for individual track.  If event==None, all events
           are selected, same for track."""
        fileset = self.data[vtype]
        if not fileset or not fileset[volt]: return

        # Color trajectories based on charge
        colors = { 1: 'tab:blue', -1: 'tab:orange' }

        # Start with detector outline to get full coordinates
        self.DetectorRZ(fileset[volt])
        plt.xlabel("R [mm]")
        plt.ylabel("Z [mm]")
        
        traj = self.getTrajectory(fileset[volt], event, track)

        # Loop over unique event/track combinations
        evlist = np.unique(traj["EventNum"])
        for ev in evlist:
            evcut = traj["EventNum"]==ev
            trklist = np.unique(traj["Track"][evcut])
            for trk in trklist:
                trajcut = evcut & (traj["Track"]==trk)
                rpt = 1e3*traj["R3"][trajcut]
                zpt = 1e3*traj["Z3"][trajcut]
                chg = traj["Charge"][trajcut][0]   # Value, not array
                plt.plot(rpt,zpt,"-",color=colors[chg])

        # Zoom in on specific region if requested
        if rlim: plt.xlim(*rlim)
        if zlim: plt.ylim(*zlim)

        # Label and save final figure
        self.TrajectoryLegend()
        plt.title(self.Title(vtype,volt))

        ftitle="Trajectories"
        if rlim or zlim: ftitle += "-zoom"
        plt.savefig(self.Filename(ftitle,vtype,volt))
        return

    ### Add labels for track types
    @staticmethod
    def TrajectoryLegend():
        """Create legend for trajectory plots, colored by charge."""
        from matplotlib.lines import Line2D
        lines = [ Line2D([0],[0], color='tab:orange', lw=3, label='Electron'),
                  Line2D([0],[0], color='tab:blue', lw=3, label='Hole') ]
        plt.legend(handles=lines,loc='lower left')
        return


    ### Plot generators
    def TrajectoryStack(self, vtype):
        """Overlay track length distribution for all voltages."""
        fileset = self.data[vtype]
        if not fileset or len(fileset)==0: return
    
        for v in reversed(self.voltage):
            if v not in fileset: continue
            hits = self.getHits(fileset[v])["StepLen"]
            plt.hist(hits*1e3,label=f"{v}V",bins=50)
            hits = None        # Ensure that allocated memory is freed
        
        plt.title(self.Title(vtype))
        plt.xlabel('Trajectory length [mm]')
        plt.legend()
        plt.savefig(self.Filename("TrajLength",vtype))
        return

    def StepStack(self, vtype, log=False):
        """Overlay step-count distribution for all voltages."""
        fileset = self.data[vtype]
        if not fileset: return

        binning = []     # Will contain result of first (vmax) plot
    
        for v in reversed(self.voltage):
            if v not in fileset: continue
            hits = self.getHits(fileset[v])

            if len(binning)==0:
                _,binning,_ = plt.hist(hits["Step"],label=f"{v}V",bins=200,log=log)
            else:
                plt.hist(hits["Step"],label=f"{v}V",bins=binning,log=log)
            
            hits = None        # Ensure that allocated memory is freed
        
        plt.title(self.Title(vtype))
        plt.xlabel('Number of Steps per track')
        plt.legend()
        plt.savefig(self.Filename("Steps",vtype))
        return

    def PathLengths(self, vtype, v):
        """Scatter plot relating flight distance with path length."""
        fileset = self.data[vtype]
        if not fileset or not v in fileset: return

        hits = self.getHits(fileset[v])
        elec = hits["Charge"]<0
        hole = hits["Charge"]>0
    
        plt.scatter(hits["Flight"][elec]*1e3,hits["StepLen"][elec]*1e3,label='Electrons')
        plt.scatter(hits["Flight"][hole]*1e3,hits["StepLen"][hole]*1e3,label='Holes')
        plt.legend()
    
        plt.title(self.Title(vtype,v))
        plt.xlabel('Flight distance [mm]')
        plt.ylabel('Trajectory length [mm]')
        plt.savefig(self.Filename("PathLengths",vtype,v))

        hits = None        # Ensure that allocated memory is freed
        return

    def PathRatio(self, vtype, v):
        """Distribution of trajectory to flight ratio."""
        fileset = self.data[vtype]
        if not fileset or not v in fileset: return
        
        hits = self.getHits(fileset[v])
        plt.hist((hits["StepLen"]/hits["Flight"]),bins=50)
        plt.title(self.Title(vtype,v))
        plt.xlabel('Trajectory/Flight')
        plt.savefig(self.Filename("PathRatio",vtype,v))

        hits = None        # Ensure that allocated memory is freed
        return

    def PathVsSteps(self, vtype, v):
        """Relationship between trajectory ratio and number of steps."""
        fileset = self.data[vtype]
        if not fileset or not v in fileset: return
        
        hits = self.getHits(fileset[v])
        elec = hits["Charge"]<0
        hole = hits["Charge"]>0
        plt.scatter(hits["Step"][elec],(hits["StepLen"]/hits["Flight"])[elec],label='Electrons')
        plt.scatter(hits["Step"][hole],(hits["StepLen"]/hits["Flight"])[hole],label='Holes')
        plt.title(self.Title(vtype,v))
        plt.xlabel('Number of Steps per track')
        plt.ylabel('Trajectory/Flight')
        plt.legend()
        plt.savefig(self.Filename("Path-Steps",vtype,v))

        hits = None        # Ensure that allocated memory is freed
        return

    def StepsVsRZ(self, file, end="Final", stepcut=1e99):
        """Heatmap of steps per track vs. hit position for single file.
           Specify which end of track as "Initial" or "Final"."""
        if not file: return

        if end == "Initial": suffix="1"
        if end == "Final":   suffix="3"
        rbranch = f"R{suffix}"
        zbranch = f"Z{suffix}"

        # Load hits for requested file, define cuts
        hits = self.getHits(file)
        cut = hits["Step"]<stepcut
        elec = hits["Charge"]<0
        hole = hits["Charge"]>0

        # Show only electron hits (hole cut is there for switching)
        rhits = hits[rbranch][elec&cut]*1e3
        zhits = hits[zbranch][elec&cut]*1e3
        plt.scatter(rhits, zhits, c=hits["Step"][elec&cut], norm='log')
        plt.xlabel(f"{end} R [mm]")
        plt.ylabel(f"{end} Z [mm]")

        plt.colorbar(label='Number of Steps')

        hits = None        # Ensure that allocated memory is freed
        return
    
    def StepsVsFinalRZ(self, vtype, v, stepcut=1e99):
        """Final position (r,z) colored by number of steps."""
        fileset = self.data[vtype]
        if not fileset or not v in fileset: return

        self.DetectorRZ(fileset[v])
        self.StepsVsRZ(fileset[v],"Final",stepcut)
        plt.title(self.Title(vtype,v)+" Electrons")

        plt.savefig(self.Filename("StepsRZ",vtype,v))
        return
    
    def StepsVsVData(self, vtype, stepcut=1e9):
        """Return maximum steps as array vs. voltage, for use with fitting."""
        fileset = self.data[vtype]
        if not fileset: return None
            
        maxstep = []
        vlist = []
        for v in self.voltage:
            if v not in fileset: continue            
            hits = self.getHits(fileset[v])
            maxstep.append(hits["Step"][hits["Step"]<stepcut].max())
            vlist.append(v)
            hits = None
                
        return vlist,maxstep

    def PlotStepsVsVoltage(self, vtype, stepcut=1e99):
        """Scatter of maximum steps vs. voltage onto current plot.
           Optional vrange should be (vmin, vmax)."""
        vlist,maxStep = self.StepsVsVData(vtype, stepcut)
        if not maxStep or len(maxStep)==0: return

        legend = self.det + " " + self.TitleField(vtype)
        plt.scatter(vlist,maxStep,label=legend)

        ylbl = "Maximum Number of Steps"
        if stepcut<1e99: ylbl += f" (below {self.ValUnits(stepcut)})"
        plt.ylabel(ylbl)
        plt.xlabel("Bias Voltage [V]")

        maxStep = None
        return

    def PlotStepsVsField(self, vtype, stepcut=1e99):
        """Scatter of maximum steps, with voltage scaled by thickness."""
        vlist,maxStep = self.StepsVsVData(vtype, stepcut)
        if not maxStep or len(maxStep)==0: return

        thick = self.getGeometry()["Height"][0]
        efield = np.array(vlist)/thick
        
        legend = self.det + " " + self.TitleField(vtype)
        plt.scatter(efield,maxStep,label=legend)

        ylbl = "Maximum Number of Steps"
        if stepcut<1e99: ylbl += f" (below {self.ValUnits(stepcut)})"
        plt.ylabel(ylbl)
        plt.xlabel("Electric Field (bias/thickness) [V/m]")

        maxStep = None
        return
    
    def PlotVoltageCurve(self, curve, stepcut=1e99):
        """Draw polynomial as function of voltage on current plot"""
        if not curve: return

        xpts,ypts = curve.linspace(500,(self.vmin,self.vmax))
        plt.plot(xpts,ypts,label=f"{curve}")
        if stepcut<1e99: plt.ylim(0,stepcut)
        return
        
    def StepsVsVoltage(self, stepcut=1e99, curve=None, suffix=''):
        """Overlaid scatter plots of maximum steps vs. voltage.
           Optional curve should be an np.Polynomial."""
        self.PlotStepsVsVoltage("E", stepcut)
        self.PlotStepsVsVoltage("U", stepcut)
        self.PlotVoltageCurve(curve, stepcut)

        plt.title(f"{self.det} Charged Tracks")
        plt.legend()

        fname = "MaxSteps"
        if stepcut<1e9: fname += f"below{self.ValUnits(stepcut)}"
        elif suffix != '': fname += '-'+suffix
        fname += f"-Voltage"
        plt.savefig(self.Filename(fname))
        return

    def CPUvsVData(self, vtype):
        """Extract CPU time per event as Numpy array vs. voltage."""
        fileset = self.data[vtype]
        if not fileset: return None
            
        cpu = []
        vlist = []
        for v in self.voltage:
            if v not in fileset: continue            
            hits = self.getJobInfo(fileset[v])
            cpu.append(hits["UserCPU"]/hits["Events"])
            vlist.append(v)
            hits = None
                
        return vlist,cpu
        
    def PlotCPUvsVoltage(self, vtype):
        """Scatter plot of CPU time per event vs. voltage."""
        vlist,cpu = self.CPUvsVData(vtype)
        if not cpu or len(cpu)==0: return

        legend = self.det + " " + self.TitleField(vtype)
        plt.scatter(vlist,cpu,label=legend)
        plt.ylabel("CPU time per event [s]")
        plt.xlabel("Bias Voltage [V]")

        cpu = None
        return

    def CPUvsVoltage(self, suffix=''):
        """Overlay scatter plots of CPU/event vs. voltage."""
        self.PlotCPUvsVoltage("E")
        self.PlotCPUvsVoltage("U")

        plt.title(f"{self.det} Charged Tracks")
        plt.legend()

        fname = "CPUtime"
        if suffix != '': fname += '-'+suffix
        fname += f"-Voltage"
        plt.savefig(self.Filename(fname))
        return

    # TODO:  Implement fitting and bounding here as functions

### END OF chargeSteps CLASS
