#!/usr/bin/env python
# Usage: ./stuck_charges.py
#
# Exported from Jupyter notebook to generate all of the desired plots
# for G4CMP-358 validation.  Probably not worth adding in command-line
# argument processing.
#
# 20250911  Michael Kelsey

from datetime import datetime

print(datetime.now(),"stuck_charges.py starting...",flush=True)
from imports import *
from chargeSteps import chargeSteps
import numpy as np

print(datetime.now(),"Imports completed.",flush=True)
help(chargeSteps)

# SuperSim V16-00-00 + elog/2125 runs (pick one!)
##print(datetime.now(),"Using V16-00-00 samples",flush=True)
##data = chargeSteps('data/EPot_highV/V16-00-00/100eh', 'HV100mm')
##data = chargeSteps('data/EPot_highV/V16-00-00/100eh-HV100mmSi', 'HV100mmSi')
##data = chargeSteps('data/EPot_highV/V16-00-00/100eh-HVeV', 'HVeV')

# New TrackLimiter runs (G4CMP-358)
##print(datetime.now(),"Using G4CMP358 samples",flush=True)
##data = chargeSteps('data/EPot_highV/G4CMP358', 'HV100mm', '-max*k')
##data = chargeSteps('data/EPot_highV/G4CMP358', 'iZIP5', '-max*k',
##                   volts=[1,4,10,25,50,75,100,200,300,400,500,600])

# SuperSim develop runs (noLimits)
print(datetime.now(),"Using noLimits samples",flush=True)
data = chargeSteps('data/EPot_highV/noLimits', 'HV100mm')
##data = chargeSteps('data/EPot_highV/noLimits', 'iZIP5')
##data = chargeSteps('data/EPot_highV/noLimits', 'HVeV')

# New intervalley scattering code (G4CMP-404)
##print(datetime.now(),"Using G4CMP404 samples",flush=True)
##data = chargeSteps('data/EPot_highV/G4CMP404','HV100mm','-newIV')

# Fix to TimeStepper for stopping voltage (G4CMP-482)
##print(datetime.now(),"Using G4CMP482 samples",flush=True)
##data = chargeSteps('data/EPot_highV/G4CMP482','HV100mm','-stepperFix')
##data = chargeSteps('data/EPot_highV/G4CMP482','HV100mm','-max*k')

print(datetime.now(),"CPUvsVoltage...",flush=True)
data.CPUvsVoltage()

print(datetime.now(),"StepStack...",flush=True)
data.StepStack("E", log=True)
data.StepStack("U", log=True)

print(datetime.now(),"PathLengths...",flush=True)
data.PathLengths("E", 10)
data.PathLengths("E", data.vmax)
data.PathLengths("U", data.vmax)

#print(datetime.now(),"PathRatio...",flush=True)
#data.PathRatio("E", data.vmax)
#data.PathRatio("U", data.vmax)

#print(datetime.now(),"PathVsSteps...",flush=True)
#data.PathVsSteps("E", 200)
#data.PathVsSteps("E", data.vmax)
#data.PathVsSteps("U", data.vmax)

print(datetime.now(),"StepsVsFinalRZ('E')...",flush=True)
data.StepsVsFinalRZ("E", 10)
data.StepsVsFinalRZ("E", 100)
data.StepsVsFinalRZ("E", 200)
#data.StepsVsFinalRZ("E", 300)
#data.StepsVsFinalRZ("E", 400)
data.StepsVsFinalRZ("E", data.vmax)

print(datetime.now(),"StepsVsFinalRZ('U')...",flush=True)
data.StepsVsFinalRZ("U", 10)
#data.StepsVsFinalRZ("U", 100)
#data.StepsVsFinalRZ("U", 200)
#data.StepsVsFinalRZ("U", 500)
#data.StepsVsFinalRZ("U", 800)
data.StepsVsFinalRZ("U", data.vmax)

print(datetime.now(),"Fitting uniform field...",flush=True)
vlist,maxUhit = np.array(data.StepsVsVData("U"))
curve = np.polynomial.Polynomial((100,0,0.1,-0.001,2.5e-6))   # Initial guess
vrange = (min(vlist),max(vlist))
quad = curve.fit(vlist, maxUhit/1e3, 4, domain=[])   # Fit in thousands
print(datetime.now(),quad)
data.StepsVsVoltage(curve=quad*1e3,suffix="FitUniform")
maxUhit=None
curve=None
quad=None

print(datetime.now(),"Fitting EPot...",flush=True)
vlist,maxEhit = np.array(data.StepsVsVData("E"))
curve = np.polynomial.Polynomial((100,0,0.1,0.,0.))   # Initial guess
vrange = (min(vlist),max(vlist))
equad = curve.fit(vlist, maxEhit/1e3, 4, domain=[])   # Fit in thousands
print(datetime.now(),equad)
data.StepsVsVoltage(curve=equad*1e3,suffix="FitEPot")
maxEhit=None
curve=None
equad=None

print(datetime.now(),"StepsVsVoltage w/curves...",flush=True)
coeffK = (36.,12.,0.53,-2.15e-5,3.85e-6)       # Rounded, for kilosteps
fit = 1e3*np.polynomial.Polynomial(coeffK,domain=vrange,window=vrange)
data.StepsVsVoltage()

# Bounding coefficient attempt for iZIP5
upperK = (40.,15.,0.6,-2e-5,4e-6)
bound = 1e3*np.polynomial.Polynomial(upperK,vrange,window=vrange)
data.StepsVsVoltage(curve=bound,suffix="Bound")

# Make consolidated plots overlaying three detector types
thedir = 'data/EPot_highV/noLimits'
snolab = chargeSteps(thedir, 'HV100mm', volts=(1,4,10,25,50,75,100))
soudan = chargeSteps(thedir, 'iZIP5', volts=(1,4,10,25,50,75,100))
hvev   = chargeSteps(thedir, 'HVeV')

print(datetime.now(),"Consolidated StepsVsVoltage",flush=True)
##soudan.PlotStepsVsVoltage("U")
soudan.PlotStepsVsVoltage("E")
##snolab.PlotStepsVsVoltage("U")
snolab.PlotStepsVsVoltage("E")
plt.plot((0,100),(1e8,1e8),'--',label="100M Step Limit")
###plt.ylim(-1e7,2e8)   # Why does Python think it's on a logy scale?
###hvev.PlotStepsVsVoltage("U")
plt.title("Charged Tracks, Germanium Detectors, EPot Files")
plt.legend(loc="upper left")
plt.savefig(thedir+"/MaxSteps-Voltage_allDet.png")

print(datetime.now(),"Consolidated StepsVsField",flush=True)
snolab.PlotStepsVsField("U")
soudan.PlotStepsVsField("U")
###hvev.PlotStepsVsField("U")
plt.xlim(0,10000)
plt.title("Charged Tracks, Uniform Field Model")
plt.legend()
plt.savefig(thedir+"/MaxSteps-EField_allDet.png")
