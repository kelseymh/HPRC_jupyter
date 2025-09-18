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
    data.CPUvsVoltage()

    print(datetime.now(),"Collection Efficiency",flush=True)
    data.OverlayEffVsVoltage("E","with 50M Charge Step Limit")
    data.OverlayEffVsVoltage("U","with 50M Charge Step Limit")
    data.EffVsVoltage("with 50M Charge Step Limit")
    # End loop over detector types

### TO DO:  Combine hv100mm and izip5 into single set of efficiency plots
