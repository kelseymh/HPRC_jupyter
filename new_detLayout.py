#!/usr/bin/env python3
#
# Creates detector channel figures suitable for slides.  See usage() for
# details.
#
# 20250115  Michael Kelsey -- Based on detectorLayout.py from Warren Perry

def usage():
    print("""
new_detLayout.py [-h] [-l] [-d det] [-z side] [-c|-s type] ROOTfile [image_file]

Arguments:
    ROOTfile	Name of ROOT file containing desired geometry to draw 
    image_file	Name of output (PNG or EPS), will use .png if no filetype

If image_file is omitted, will create <detType>.png file.

Options:
    -h		Display this usage information
    -l		Include channel labels on diagram
    -c, -s	Channel type to draw (1 or "TES", 2 or "FET")
    -d		Detector number from ROOT file to draw
    -v		Turn on verbose output (currently just reports arguments)
    -z		Detector face to draw (+1 or "top", -1 or "bottom")

If '-c' option is omitted, both phonon and charge channels will be drawn.
If '-d' option is omitted, will use detNum=0 (assumes a single-detector sim).
If '-z' option is omitted, both faces will be drawn in a 'two-frame' figure.

For interactive use in Jupyter, users should use the dmc.ChannelShape module.

Requires: Matplotlib, Numpy, SciPy, ROOT
""")

### CONFIGURATION ###

import matplotlib.pyplot as plt
import os, sys
import dmc

### MAIN PROGRAM ###

def main():
    settings = getargs()

    myDet = dmc.DetectorShape().Load(settings['file'],settings['detnum'])

    # Use detector subgroup if user requested specific channels
    drawDet = myDet
    if settings['sensor']=='TES': drawDet = myDet.TES
    if settings['sensor']=='FET': drawDet = myDet.FET

    # For two-sided detectors, if requested do double plot
    if settings['zside']==0 and drawDet.hasZside(+1) and drawDet.hasZside(-1):
        # NOTE: Extra inch in figsize() ensures both panels are 5x5 square
        fig,(ax1,ax2) = plt.subplots(1,2, figsize=(11,5), dpi=200)
        drawDetector(ax1, drawDet, +1, settings)
        drawDetector(ax2, drawDet, -1, settings)
    else:				# Figure out which side to draw
        zside = -1 if settings['zside']!=1 else +1

        if not drawDet.hasZside(zside):	# Nothing to draw, just quit
            print(myDet, f" has no {('+' if zside>0 else '-')}Z",
                  f"{settings['sensor']} channels", file=sys.stderr)
            sys.exit(0)

        ax = plt.figure(figsize=(5,5), dpi=200)
        drawDetector(ax, drawDet, zside, settings)

    # Finished with drawing, create output for user and exit
    if settings['image']=="": settings['image'] = myDet.detName+".png"
    plt.savefig(settings['image'])
    print(f"Detector diagram written to {settings['image']}")

### End of main()

def drawDetector(ax, det, zside, settings):
    """Draw single side of specified detector in active figure."""
    det.Draw(ax, zside, frame=True, labels=settings['labels'])

    # Construct title string from detector information being plotted
    title = f"{det.getDetName(det.detType)} {'+' if zside>0 else '-'}Z "

    types = list({ch.chantype for ch in det.values()})
    if len(types)==1: title += dmc.ChannelShape.getTypeName(types[0])
    title += " channels"

    ax.set_title(title)

def getargs():
    """Return arguments from the command line as a dictionary, for easier use.
       Output: settings = { 'file':	<name-of-ROOT-file>,
                            'image':    <name-of-output-image>,
                            'detnum':   <detector number> from -d,
                            'zside':    <+1|-1|0> from -z
                            'sensor':   <sensor type> from -c,
                            'labels':   <True|False> from -l,
                            'verbose':  <True|False> from -v,
                          }
    """
    import getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'c:d:s:z:lvh')
    except getopt.GetoptError as err:
        sys.exit(1)

    if args is None or len(args)==0:
        usage()
        sys.exit(2)

    settings = { 'file': args[0],
                 'image': args[1] if len(args)>1 else '',
                 'detnum': 0,
                 'zside': 'both',	# Will convert to 0 below
                 'sensor': '',
                 'labels': False,
                 'verbose': False,
               }

    for o,a in opts:	# Capture user input here, translate below
        if   o == '-c': settings['sensor']  = a
        elif o == '-s': settings['sensor']  = a
        elif o == '-d': settings['detnum']  = int(a)
        elif o == '-l': settings['labels']  = True
        elif o == '-v': settings['verbose'] = True
        elif o == '-z': settings['zside']   = a
        elif o == '-h':
            usage()
            exit()

    # Make side selection integer to pass into DetectorShape.Draw()
    if isinstance(settings['zside'], str):
        sides = { "top":1, "+":1, "bottom":-1, "-":-1, "both":0 }
        zarg = settings['zside'].casefold()
        if zarg in sides: settings['zside']=sides[zarg]

    # Make channel type consistent to do selection
    if isinstance(settings['sensor'], int):
        types = { 1:"TES", 2:"FET" }
        settings['sensor'] = types[settings['sensor']]
    else:
        settings['sensor'] = settings['sensor'].upper()

    # Report configuration if requested
    if settings['verbose']:
        print(f"{sys.argv[0]} settings:\n{settings}")

    return settings
### End of getargs() ###


### COMMAND-LINE CALL ###

if __name__ == "__main__":  main()
    
### END OF FILE ###
