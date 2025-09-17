#!/usr/bin/env python
# Usage: ./step-limits.py
#
# Plot effect of step limits on ER events for G4CMP-358 validation.
#
# 20250914  Michael Kelsey

from datetime import datetime

print(datetime.now(),"step-limits.py starting...",flush=True)
from chargeSteps import chargeSteps
print(datetime.now(),"chargeSteps done",flush=True)
from dmc import phonon_eff
print(datetime.now(),"phonon_eff done",flush=True)
import numpy as np
print(datetime.now(),"numpy done",flush=True)
import matplotlib.pyplot as plt
print(datetime.now(),"matplotlib done",flush=True)
###from imports import *
print(datetime.now(),"imports skipped",flush=True)

print(datetime.now(),"Imports completed.",flush=True)

print(datetime.now(),"Using 10keV samples",flush=True)
data = chargeSteps('data/EPot_highV/10keVER', 'HV100mm', infix='-max*k')
##data = chargeSteps('data/EPot_highV/10keVER', 'iZIP5', infix='-max*k')

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
    ax.set_title(f"{data.det} 10 keV ER with Charge Step Limit")
    ax.set_xlabel("Collection Efficiency (PhononE/Eexpected)")
    ax.set_ylabel("Events / 0.01%")

    for v in data.voltage:
        if v not in dset: continue

        print(datetime.now(),f" loading {v}V{vtype} files...",flush=True)
        geom,events,hits = phonon_eff.load(dset[v])
        lbl = f"{v}V"
        ax.hist(events["PhEff"],bins=150,range=(0.9,1.05),label=lbl,alpha=0.3)

        effhist,bins = np.histogram(events["PhEff"],bins=150,range=(0.9,1.05))
        imax = np.argmax(effhist)
        effset[vtype]["Eff"].append((bins[imax]+bins[imax+1])/2.)
        effset[vtype]["Volts"].append(v)

        geom=None	# Avoid memory leaks
        events=None
        hits=None

    ax.legend(loc="upper left")
    plt.savefig(data.Filename("Efficiency",vtype))
# End loop over data sets

_,ax = plt.subplots(figsize=(4.5,3))
ax.set_title("10 keV ER Efficiency with Charge Step Limit")
ax.set_xlabel("Voltage [V]")
ax.set_ylabel("Collection Efficiency (mode)")

for vtype in effset.keys():
    ax.scatter(effset[vtype]["Volts"],effset[vtype]["Eff"],
               label=data.Title(vtype))

ax.legend()
plt.savefig(data.datadir+f"/Efficiency-Voltage_{data.det}.png")
