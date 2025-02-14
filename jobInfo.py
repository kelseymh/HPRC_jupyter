"""Retrieve general information about DMC simulation from ROOT file.
   Requires: Numpy, ROOT, CATs from SuperCDMS

   Usage: import dmc
          myJob = dmc.jobInfo("<root-file>")
          print(myJob.versions())
          print(myJob.macro())

   If a list of files or a fileglob is provided, the functions will return
   list of arrays, with each entry corresponding to a single file.

   Michael Kelsey <kelsey@tamu.edu>, Texas A&M University 2025
"""

# 20250201  Inspired by functionality in Emre DuVarci's OF1X2.py fitter

import numpy as np
from cats.cdataframe import CDataFrame

### MAIN DRIVER CLASS -- ONE PER FILE ###

class jobInfo:
    """Load non-physics information from DMCintermediate (or any SuperSim) ROOT
       file.  Create one instance of this for each ROOT file or set of ROOT
       files to be processed.
    
       Constructor arguments:
       files = Name of ROOT file or list of names (to use as TChain)
       verbose = (optional) True reports function calls, progress
    """

    # Constructor
    def __init__(self, files, verbose=False):
        """Constructor: Stores filenames for loading data later."""
        self.verbose = verbose
        self.files = self._expand(files)
        return
        
    # Useful diagnostic functions
    def setVerbose(self, verbose=True):
        self.verbose = verbose
            
    def printVerbose(self, *args, **kwargs):
        """Calls print() only if verbose=True."""
        if self.verbose: print(*args, **kwargs)

    def __repr__(self):
        return f"<jobInfo('{self.files}','{self.verbose}')"

    # Information retrieval functions

    def versions(self):
        """Return dictionary of packages and versions, unwrapped from Numpy."""
        self.printVerbose(f"jobInfo.versions ({len(self.files)} files)")

        # Just load everything all at once
        varray = CDataFrame("G4SettingsInfoDir/Versions", self.files).AsNumpy()
        self.printVerbose(f"  Versions has {len(varray)} packages")

        # If user specifed single file, unwrap the single-element arrays
        if len(self.files)==1:
            for key,val in varray.items(): varray[key] = val[0]

        return varray

    def macro(self):
        """Return stored SuperSim job macro from each file."""
        self.printVerbose(f"jobInfo.macro ({len(self.files)} files)")

        # TODO: Should we build a dict using filenames instead?
        mset = [ self._fileMacro(f) for f in self.files ]
        return mset if len(mset)>1 else mset[0]

    def events(self):
        """Return list of event numbers in each file."""
        self.printVerbose(f"jobInfo.events ({len(self.files)} files)")

        # TODO: Should we build a dict using filenames instead?
        eset = [ self._fileEvents(f) for f in self.files ]
        return eset if len(eset)>1 else eset[0]

    # Internal support functions, called above

    def _fileMacro(self,file):
        """Load TMacro from specified single file, return as long string."""
        self.printVerbose("jobInfo._getMacro",file)

        import ROOT		# Avoid ROOT contaminating user's namespace
        f = ROOT.TFile.Open(file)
        macro = f.Get("G4SettingsInfoDir/SuperSim_Macro")
        f.Close()

        self.printVerbose("  macro has",len(macro.GetListOfLines()),"lines")

        # Unpack ROOT TMacro object into multiline string.
        return '\n'.join([str(l) for l in macro.GetListOfLines()])

    def _fileEvents(self,file):
        """Get list of event numbers from specified single file."""
        self.printVerbose("jobInfo._fileEvents",file)

        reader = CDataFrame("G4SimDir/mcevent", file)
        return reader.AsNumpy(["EventNum"])["EventNum"]
        
    def _expand(self,fileglobs):
        """Input may be a single string or a list of strings.  Any strings
           containing a glob wildcard expression will be expanded.  The
           result is returned as a list."""
        self.printVerbose("jobInfo._expand",fileglobs)

        import glob
        fileset = []
        if isinstance(fileglobs, str):
            fileset = glob.glob(fileglobs)
        else:
            for f in fileglobs: fileset += glob.glob(f)

        self.printVerbose("  after expansion: ",fileset)

        fileset = np.unique(fileset).tolist() 	# Strip out duplicates
        self.printVerbose("  unique names: ",fileset)

        return fileset

### END OF jobInfo CLASS ###
