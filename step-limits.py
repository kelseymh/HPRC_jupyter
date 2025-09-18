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
plt.rcParams['figure.figsize'] = (5,3.75)
plt.rcParams['figure.autolayout'] = True
print(datetime.now(),"matplotlib done",flush=True)
###from imports import *
print(datetime.now(),"imports skipped",flush=True)

print(datetime.now(),"Imports completed.",flush=True)

print(datetime.now(),"Using 10keV samples",flush=True)
hv100mm = chargeSteps('data/EPot_highV/10keVER', 'HV100mm', infix='-max*k')
izip5   = chargeSteps('data/EPot_highV/10keVER', 'iZIP5', infix='-max*k')

for data in (hv100mm,izip5):
    print(datetime.now(),f"{data.det} CPU time",flush=True)
    plt.clf()      # Discard any previously created figure
    plt.figure()
    data.CPUvsVoltage()

    print(datetime.now(),"Collection Efficiency",flush=True)

    effset = { "E": { "Volts": [], "Eff": [] },
               "U": { "Volts": [], "Eff": [] }
    }

    for vtype in ["E","U"]:
        dset = data.data[vtype]
        if len(dset)==0: continue

        plt.clf()      # Discard any previously created figure
        plt.figure()
        plt.title(data.Title(vtype)+"\nwith 50M Charge Step Limit")
        plt.xlabel("Collection Efficiency (PhononE/Eexpected)")
        plt.ylabel("Events / 0.1%")
        
        for v in data.voltage:
            if v not in dset: continue
            
            print(datetime.now(),f" loading {v}V{vtype} files...",flush=True)
            geom,events,hits = phonon_eff.load(dset[v])
            lbl = f"{v}V"
            plt.hist(events["PhEff"],bins=60,range=(0.95,1.01),
                     label=lbl,alpha=0.3)
            
            effhist,bins = np.histogram(events["PhEff"],bins=60,
                                        range=(0.95,1.01))
            imax = np.argmax(effhist)
            effset[vtype]["Eff"].append((bins[imax]+bins[imax+1])/2.)
            effset[vtype]["Volts"].append(v)
            
            geom=None	# Avoid memory leaks
            events=None
            hits=None
            
        plt.legend(loc="upper left")
        plt.savefig(data.Filename("Efficiency",vtype))
        # End loop over data sets

    plt.clf()      # Discard any previously created figure
    plt.figure()

    for vtype in effset.keys():
        plt.scatter(effset[vtype]["Volts"],effset[vtype]["Eff"],
                    label=data.Title(vtype))

    plt.title(f"{data.det} 10 keV ER Efficiency\nwith Charge Step Limit")
    plt.xlabel("Voltage [V]")
    plt.ylabel("Collection Efficiency (mode)")
    plt.ylim(0.94,1.02)
    plt.legend(loc="lower left")
    plt.savefig(data.datadir+f"/Efficiency-Voltage_{data.det}.png")
    # End loop over detector types

### TO DO:  Combine hv100mm and izip5 into single set of efficiency plots
