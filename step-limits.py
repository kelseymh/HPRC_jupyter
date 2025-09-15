#!/usr/bin/env python
# Usage: ./step-limits.py
#
# Plot effect of step limits on ER events for G4CMP-358 validation.
#
# 20250914  Michael Kelsey

from datetime import datetime

print(datetime.now(),"step-limits.py starting...",flush=True)
from imports import *
from chargeSteps import chargeSteps
from dmc import phonon_eff
import numpy as np

print(datetime.now(),"Imports completed.",flush=True)

print(datetime.now(),"Using 10keV samples",flush=True)
data = chargeSteps('data/EPot_highV/10keVER', 'HV100mm')
##data = chargeSteps('data/EPot_highV/10keVER', 'iZIP5')

print(datetime.now(),"CPU time",flush=True)
data.CPUvsVoltage()

print(datetime.now(),"Collection Efficiency",flush=True)

effset = { "E": { "Volts": [], "Eff": [] },
           "U": { "Volts": [], "Eff": [] }
}

for vtype in ["E","U"]:
    dset = data.data[vtype]
    if len(dset)==0: continue

    _,ax = plt.subplots(figsize=(4.5,3))
    ax.title(f"{data.det} 10 keV ER with Charge Step Limit")
    ax.xlabel("Collection Efficiency (PhononE/Eexpected)")
    ax.ylabel("Events / 0.01%")

    for v in data.voltage:
        if v not in dset: continue

        print(datetime.now(),f" loading {v}V{vtype} files...",flush=True)
        geom,events,hits = phonon_eff.load(dset[v])
        lbl = f"{v}V"

        effhist,bins = np.histogram(events["PhEff"],bins=150,range=(0.9,1.05))
        ax.plot(effhist,bins,label=lbl,alpha=0.3)

        effset[vtype]["Volts"] += v

        imax = np.argmax(effhist)
        effset[vtype]["Eff"] += (bins[imax]+bins[imax+1])/2.

        geom=None	# Avoid memory leaks
        events=None
        hits=None

    ax.legend(loc="upper left")
    plt.savefig(data.Filename("Efficiency",vtype))
# End loop over data sets

_,ax = plt.subplots(figsize=(4.5,3))
ax.title("10 keV ER Efficiency with Charge Step Limit")
ax.xlabel("Voltage [V]")
ax.ylabel("Collection Efficiency (mode)")

for vtype in effset.keys():
    ax.scatter(effset[vtype]["Volts"],eff[vtype]["Eff"],label=data.title(vtype))

ax.legend()
plt.savefig(data.datadir+"/Efficiency-Voltage_allDet.png")
