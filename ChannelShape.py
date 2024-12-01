import numpy as np
from math import *

# This is a Python version of supersim/CDMSinfo/CDMSChannelShape
#
# NOTE:  The functionality of the classes in this package is not complete.
#        Look for "TODO" comments below to see what still needs to be
#        implemented.
#
# 20231019  Add additional optional flats to support full detector
# 20240615  New DetectorShape dictionary class for full detector
# 20240630  Add optional (default None) channel type to ChannelShape

class ChannelShape:
    """Implementation of CDMSChannelShape in Python.  Once created, users may
       call chan.contains((x,y)) or chan.contains((x,y,z)), which will return
       True/False if the given position is within the channel boundary or not.
    
       Other useful functions:
       chan.deltaphi(): Returns azimuthal width of channel in radians
       chan.phiBetween(phi): Returns True if phi is within channel range
       ChannelShape.phinorm(phi): Normalizes phi between -pi and pi.
       ChannelShape.cart2pol((x,y)): Returns (r,phi) corresponding to (x,y)
       ChannelShape.pol2cart((r,phi)): Returns (x,y) corresponding to (r,phi)

       NOTE: Angles must be specified in radians; use math.radians() to convert

       Data members are created with constructor arguments, in this order:
       .name = String identifier for channel (e.g., "PAS1")
       .chantype = Code or none for channel sensor type (TES=1, FET=2)
       .rmin = minimum radius of circular channel
       .rmax = maximum radius of circular channel
       .phimin = minimum (in CCW sense) phi angle of circular channel
       .phimax = maximum (in CCW sense) phi angle of circular channel
       .zside = +1 for top face, -1 for bottom face
       .zsurface = thickness to treat as "surface" for deciding .contains()
       .xRmin = minimum absolute X value for inner left/right flats
       .xRmax = maximum absolute X value for outer left/right flats
       .yRmin = minimum absolute Y value for inner top/bottom flats
       .yRmax = maximum absolute Y value for outer top/bottom flats
       .phiFlat (Rmax,phi) = additional flat(s) at fixed angle
       .phiminXCcut = for C-shaped square channels, X limit at phimin
       .phiminYCcut = for C-shaped square channels, Y limit at phimin
       .phimaxXCcut = for C-shaped square channels, X limit at phimax
       .phimaxYCcut = for C-shaped square channels, Y limit at phimax

       A set of flats are also defined to help with drawing and with boundary
       points:

       .flats = [ ChannelFlat(phi,R,rFlat), ... ]

       Since a given channel could have both inner and outer flats, the above
       cannot be a simple dictionary in phi; it could be a dictionary on
       (R,phi) coordinates, or we will provide a search function.
    """

    # This function may be called with either a single point or a set of points
    def contains(self, pos):
        """Return True if (x,y,z) position is within channel boundary"""
        assert isinstance(pos, (tuple, list, np.ndarray))
        
        if (isinstance(pos, tuple)):
            return self._containsPoint(pos)
        else:
            return self._containsNumpy(pos)
    
    # Angular width of channel
    def deltaphi(self):
        """Return angular width of channel, taking account of periodicity"""
        dphi = self.phimax-self.phimin;
        if (abs(dphi) < 1e-6): return 2.*pi    # [0,0] -> 2pi
        if (dphi > 0.): return dphi            # [-pi/2,pi/2] -> pi
        return (dphi+twopi)                    # [pi/2,0] -> 3pi/2

    # Add single flat at specified phi angle and radius
    def addFlat(self, Rflat, phi):
        """Slice off a chord at phi with center radius Rflat"""
        self.phiFlat.append((Rflat,phi))
        self._appendFlat(phi, self.rmax, Rflat)
        return self

    # Check if channel is of TES or FET type
    def isTES(self):
        return self.chantype==1

    def isFET(self):
        return self.chantype==2
        
    # Constructor
    def __init__(self, cname, ctype=None, rlo=0., rhi=0., flo=0., fhi=0., z=0.,
                 xin=0., xout=0., yin=0., yout=0., zmin=0.,
                 floXC=1000., floYC=1000., fhiXC=1000., fhiYC=1000.):
        """Constructor to initialize all 'public' data members"""
        self.name = cname
        self.chantype = ctype if ctype else self._chanType(cname)
        self.rmin = rlo
        self.rmax = rhi
        self.phimin = self.phinorm(flo)
        self.phimax = self.phinorm(fhi)
        self.zside = z
        self.zsurface = zmin
        self.xRmin = xin
        self.xRmax = xout
        self.yRmin = yin
        self.yRmax = yout
        self.phiFlat = []
        self.phiminXCcut = floXC
        self.phiminYCcut = floYC
        self.phimaxXCcut = fhiXC
        self.phimaxYCcut = fhiYC
        self.cosPhimin = cos(flo)     # Precompute trig values for convenience
        self.sinPhimin = sin(flo)
        self.cosPhimax = cos(fhi)
        self.sinPhimax = sin(fhi)

        # Populate flats for use in drawing and reflections
        self._fillFlats()
    ### End of __init__()

    # Guess channel type (TES or FET) based on the channel name string
    # TES channels are named "P$" or "P$S#", or are named "Ch#" (HVeV).
    # FET channels are named "Q$#", but "PTOP" and "PBTM" are also used.
    @staticmethod
    def _chanType(cname):
        TEStype,FETtype = 1,2		# Pythonic replacement for enums
        if (cname[0]=="Q"): return FETtype
        if (cname=="PTOP" or cname=="PBTM"): return FETtype
        if (cname[0]=="P" or cname[0]=="C"): return TEStype
        return None


    # Nicely formatted output: __str__ for print(), __repr__ for e.g., lists
    def __repr__(self):
        return (f"<ChannelShape('{self.name}',{self.chantype},{self.rmin}"+
                f",{self.rmax},{self.phimin},{self.phimax},{self.zside})>")

    def __str__(self):
        _typenames = { 1: "TES", 2: "FET", None: "NA" }	# Map sensor type code to name
        desc  = f"Channel {self.name} ({_typenames[self.chantype]}):"
        desc += f" r[{self.rmin}..{self.rmax} mm]"
        desc += f" phi[{degrees(self.phimin)}..{degrees(self.phimax)} deg]"
        desc += (" +z" if self.zside>0 else " -z")

        if (self.phiminXCcut<1000. or self.phiminYCcut<1000.):
            desc += f"\n  phiMin C-cut @ x {self.phiminXCcut}, y {self.phiminYCcut} mm"
        if (self.phimaxXCcut<1000. or self.phimaxYCcut<1000.):
            desc += f"\n  phiMax C-cut @ x {self.phimaxXCcut}, y {self.phimaxYCcut} mm"

        if (self.flats and len(self.flats)>0):
            desc += f"\n  {len(self.flats)} flats:"
            for flat in self.flats:
                if (self.flats.index(flat)>0): desc += "\n          "
                desc += f" {flat}"
        
        import textwrap
        return textwrap.dedent(desc)
        
    # Populate default X- and Y-flats as ChannelFlat objects
    def _fillFlats(self):
        self.flats = []
        if (self.xRmin > 0 and self.xRmin < self.rmin):
            self._appendFlat(0., self.rmin, self.xRmin)
            self._appendFlat(pi, self.rmin, self.xRmin)

        if (self.xRmax > 0 and self.xRmax < self.rmax):
            self._appendFlat(0., self.rmax, self.xRmax)
            self._appendFlat(pi, self.rmax, self.xRmax)

        if (self.yRmin > 0 and self.yRmin < self.rmin):
            self._appendFlat(pi/2, self.rmin, self.yRmin)
            self._appendFlat(-pi/2, self.rmin, self.yRmin)

        if (self.yRmax > 0 and self.yRmax < self.rmax):
            self._appendFlat(pi/2, self.rmax, self.yRmax)
            self._appendFlat(-pi/2, self.rmax, self.yRmax)
    
    # Attach separate ChannelFlat object to aid with boundary handling
    def _appendFlat(self, phi, radius, rflat):
        """Add new entry to list of flat cuts for use with drawing"""
        if not self.flats: self.flats = []
        self.flats.append(ChannelFlat(phi,radius,rflat))
        ### TODO: Sort flats to go from phi==0 at rmax around counterclockwise
        
    # Implementation of contains() for use with a single point (2D or 3D)
    def _containsPoint(self, pos):
        """Return True if position (x,y) or (x,y,z) is within boundary"""
        r,phi = self.cart2pol(pos)
 
        return ((len(pos)==2 or self._goodZ(pos[2])) and        # 2D position or near surface
                (r==0 or self._phiBetween(phi)) and             # Phi range for wedges
                ((r <= self.rmax) and self._outerFlat(pos)) and # Within outer rim of annulus
                ((r >= self.rmin) or self._innerFlat(pos)) and  # Outboard of inner rim of annulus
                self._extraFlat(pos) and                        # Within any diagonal flats
                self._goodXCcut(pos) and self._goodYCcut(pos)   # Check for C-shaped channels
               )

    # Implementation of contains() for use with Numpy (manual version of "vectorize()")
    def _containsNumpy(self, posarray):
        """Return Numpy array of True/False for each listed position"""
        return np.asarray([self._containsPoint(ipos) for ipos in posarray])
        
    # Utility functions used mostly by contains(); in C++, these would be "private:"
    def _phiBetween(self, phi):
        """Test if phi is within the channel's phi boundaries, allowing for
           negative angles"""
        return self.inPhi(phi, self.phimin, self.phimax)

    def _goodZ(self, z):
        """Ensure that Z coordinate is 'close enough' to surface"""
        return (z*self.zside > 0. and abs(z) >= self.zsurface);

    # Handle outer and inner flats on circular channels
    def _outerFlat(self, pos):
        """Check whether (x,y) position is within flats on outside of channel"""
        x,y = pos[:2]
        return ( (self.yRmax==0. or (y<self.yRmax and y>-self.yRmax)) and
                 (self.xRmax==0. or (x<self.xRmax and x>-self.xRmax)) )

    def _innerFlat(self, pos):
        """Check whether (x,y) position is outboard of flats on inner rim of channel"""
        x,y = pos[:2]
        return ( (self.yRmin>0. and (y>self.yRmin or y<-self.yRmin)) or      # Flats add to inside
                 (self.xRmin>0. and (x>self.xRmin or x<-self.xRmin)) )

    def _extraFlat(self, pos):
        """Check whether (x,y) position is within additional diagonal flats"""
        if len(self.phiFlat)==0: return True
        rpos,fpos, = self.cart2pol(pos)
        # Comparison uses rotated frame with flat normal to +X' axis
        for fl in self.phiFlat:
            x,y = self.pol2cart((rpos, fpos-fl[1]))    # Rotate point to X' frame
            if (x > fl[0]): return False
        return True
            
    # Handle flat boundaries for C-shaped channels
    def _goodXCcut(self, pos):
        """Check for a C-shaped channel with a vertical edge (cut at X value)"""
        x,y = pos[:2]
        # Only look at X coordinate if Y is in same half-plane as phi end
        return ( ( self.phiminXCcut>=1000. or 
                   ((y<0 or x<self.phiminXCcut) if (self.sinPhimin>=0.)
                        else (y>0 or x>self.phiminXCcut)) )
                and
                 ( self.phimaxXCcut>=1000. or 
                   ((y<0 or x>self.phimaxXCcut) if (self.sinPhimax>=0.)
                        else (y>0 or x<self.phimaxXCcut)) )
               )
        
    def _goodYCcut(self, pos):
        """Check for a C-shaped channel with a horizontal edge (cut at Y value)"""
        x,y = pos[:2]
        # Only look at Y coordinate if X is in same half-plane as phi end
        return ( ( self.phiminYCcut>=1000. or 
                   ((x<0 or y>self.phiminYCcut) if (self.cosPhimin>=0.)
                        else (x>0 or y<self.phiminYCcut)) )
                and
                 ( self.phimaxYCcut>=1000. or 
                   ((x<0 or y<self.phimaxYCcut) if (self.cosPhimax>=0.)
                        else (x>0 or y>self.phimaxYCcut)) )
               )

    # Guess channel type (TES or FET) based on the channel name string
    # TES channels are named "P$" or "P$S#", or are named "Ch#" (HVeV).
    # FET channels are named "Q$#", but "PTOP" and "PBTM" are also used.
    @staticmethod
    def _chanType(cname):
        TEStype,FETtype = 1,2		# Pythonic replacement for enums
        if (cname[0]=="Q"): return FETtype
        if (cname=="PTOP" or cname=="PBTM"): return FETtype
        if (cname[0]=="P" or cname[0]=="C"): return TEStype
        return None

    @staticmethod
    def phinorm(phi):
        """Normalize specified angle to [-pi,pi] range"""
        # C++ has: return ((phi>pi) ? phi-twopi : (phi<-pi) ? phi+twopi : phi);
        return (phi-2.*pi if phi>pi else phi+2.*pi if phi<-pi else phi)

    @staticmethod
    def inPhi(phi, phimin, phimax):
        """Return true if phi is contained in the specified phi range.  Handles
           situation around 180 degrees, where phimin<pi, while phimax>-pi
        """
        phi0 = ChannelShape.phinorm(phi)
        return ((phimin<=phi0 and phi0<=phimax) if (phimin<phimax)
                    else (phi0>=phimin or phi0<=phimax))
        
    # Conversions between 2D polar and Cartesian coordinates, for convenience
    @staticmethod
    def cart2pol(pos):
        """Returns r,phi coordinates of input Cartesian two-vector"""
        x,y = pos[:2]
        return (sqrt(x**2+y**2), atan2(y,x))

    @staticmethod
    def pol2cart(pos):
        """Returns x,y coordinates of input Polar two-vector"""
        r,phi = pos[:2]
        return (r*cos(phi), r*sin(phi))

### END OF ChannelShape ###


class ChannelFlat:
    """Utility class to store geometric information about flats: radius and
       angle at center of flat, angular range, and transverse length.

       Constructor takes cylindrical radius (channel or detector) and the
       central position of the flat (r,phi).  Data members for reference are:
       .phi    = phi position at midpoint of flat
       .rEdge  = radius of cylindrical edge (detector or individual channel)
       .rFlat  = radius at midpoint of flat
       .phiMin = minimum phi position where flat touches cylinder
       .phiMax = maximum phi position where flat touches cylinder
       .length = length of flat (useful for drawing detector or channel shape)

       NOTE: For flat at 180 degrees, .phiFMin will be close to +pi, .phiFMax
       will be close to -pi.
    """

    def __init__(self, phi, redge, rflat):
        self.rEdge = redge
        self.rFlat = rflat
        self.phi   = ChannelShape.phinorm(phi)

        theta = acos(rflat/redge)
        self.phiMin = ChannelShape.phinorm(self.phi-theta)
        self.phiMax = ChannelShape.phinorm(self.phi+theta)
        self.length = 2*sqrt(redge**2 - rflat**2)

    def rAtPhi(self, phi):
        """Returns radial position of point at given phi position"""
        # FIXME: If outside the flat range, should it return None or full radius?
        return (self.rFlat / cos(phi-self.phi)) if self.inPhi(phi) else None
            
    def inPhi(self, phi):
        """Returns true if given phi position is within range of flat"""
        return ChannelShape.inPhi(phi, self.phiMin, self.phiMax)

    # Nicely formatted output: __str__ for print(), __repr__ for e.g., lists
    def __repr__(self):
        return f"<ChannelFlat({self.phi},{self.rEdge},{self.rFlat})>"

    def __str__(self):
        desc = f"Flat @ phi {degrees(self.phi)}"
        desc += f" [{degrees(self.phiMin)}..{degrees(self.phiMax)}] deg"
        desc += f", r {self.rFlat} mm"
        return desc
        
    ### TODO:  Define all six comparison operators to allow automatic sorting
    #def __lt__(self):

### END OF ChannelFlat ###


class DetectorShape(dict):
    """Complete set of channels to define detector geometry for channel mapping.
       Includes shapes for each individual channel (phonon or charge), along
       with a "Det" channel shape corresponding to the detector outline as
       viewed from above.

       Constructor takes either the name or type code for the detector, and an
       optional predefined dictionary of channels.  Both the detector type and
       corresponding name will be stored as data members (.detType and .detName,
       respectively), for convenience.
       
       If a dictionary is not provided, user must call one of the following
       functions to build a set of channels for use elsewhere.  If the user
       assigns the channels individually, the overall "Det" channel should be
       done last, so that `findChannel()` treats it as a fallback.

       DetectorShape().Load(DMCfile, detnum=0)  Get channels from ROOT file

       DetectorShape().Add(ChannelShape)        Register ChannelShape object

       If you the second method to populate the geometry manually, you should
       call .FillGroups() when you're finished.  This method will construct
       a separate sets of TES and FET channels, which can be used below.

       .findChannel(pos) returns the name of the channel where position (3D) is
       located.  If none of the readout channels contain the position, then
       "Det" should be returned if the position is within the detector radius.

       .TES.findChannel() or .FET.findChannel() will limit the search to just
       those channels, and will not include the full "Det".
    """

    def __init__(self, det=None, chanset=None):
        """User may pass 'det' as either a detector number or name string,
           both will be stored.  If a dictionary is passed as 'chanset',
           it will also be stored."""

        if chanset is None:
            super(DetectorShape,self).__init__()
        else:
            super(DetectorShape,self).__init__(chanset)

        if isinstance(det, str):
            self.detName = det
            self.detType = self.getDetType(det)
        elif isinstance(det, (int,float)):
            self.detName = self.getDetName(det)
            self.detType = det

        if (len(self)>0): self.FillGroups()

    # NOTE: Can't use function aliasing here, since _loadXXX defined later
    def Load(self, dmcfile, detnum=0):
        """Usage: DetectorShape().Load(<root-file>, <detectorID>)
           Open specified DMC/SuperSim ROOT file and extract channel data for
           detector <detnum>."""
        self._loadTTree(dmcfile, detnum)     	# Switch to _loadRDF later
        self.FillGroups()			# Create .TES and .FET subsets
        return self
        
    def Add(self, chan):
        """Insert specified ChannelShape object into dictionary"""
        assert(isinstance(chan, ChannelShape))
        self[chan.name] = chan
        return self

    # Figure out which channel (if any) contains the specified location
    def findChannel(self, pos):
        """Return name of channel containing (x,y) or (x,y,z) position.  If
           position is outside of defined channel boundaries, but inside
           overall detector, should return 'Det' (assuming 'Det' exists and
           was registered last).
        """
        matches = []
        for ch in self:
            if self[ch].chantype and self[ch].contains(pos): matches.append(ch)
                
        if len(matches)==1: return matches[0]
        elif len(matches)>1: return matches
        else: return None
        
    # This function may be called with either a single point or a set of points
    def contains(self, pos):
        """Return True if the specified (x,y,z) position is within the overall
           detector boundary"""
        found = False
        if ("Det" in self):
            found = self["Det"].contains(pos)
        else:
            for chan in self:
                found |= chan.contains(pos)
        return found

    # Find point on detector edge nearest to specified point
    def getPointOnEdge(self, pos, direction=None):
        """Return (x,y) point on detector rim nearest to specified input point.
           If optional direction vector is provided, The position where
           'direction' intersects with the detector will be returned."""

        # If no direction specified, use radial vector
        if (direction is None):
            _,phi = ChannelShape.cart2pol(pos)
            return self.getPointAtPhi(phi)
        else:
            return self._getIntersectingPoint(pos, direction)
        
    # Find outward normal from detector at given point on edge
    def getOutwardNormal(self, pos):
        """Return (dx,dy) outward pointing normal at the specified position on
           the detector edge.   If position is not on detector edge, results
           are unreliable.  Generally, the output of `getPointOnEdge()` should
           be used as input."""
        r,phi = ChannelShape.cart2pol(pos)
    
        # Most of detector is cylindrical, so start with the radial vector
        norm = ChannelShape.pol2cart((1.,phi))

        # Check if position is in regions of flats
        for flat in self["Det"].flats:
            # NOTE: "Det" only has flats at rmax, so this will be unique
            if (flat.inPhi(phi)):
                # Normal to flat is central phi
                norm = ChannelShape.pol2cart((1.,flat.phi))

        return norm

    # Scan flats and args to determine radius of edge at given phi
    def getPointAtPhi(self, phi):
        """Find position at edge of detector nearest to phi position."""
        r = self["Det"].rmax		# Most of detecor is cylindrical

        for flat in self["Det"].flats:	# Scan flats for exceptional areas
            # NOTE: "Det" only has flats at rmax, so this will be unique
            if (flat.inPhi(phi)): r = flat.rAtPhi(phi)
           
        return ChannelShape.pol2cart((r,phi))

    # Find intersection of detector shape with line (point, direction)
    def _getIntersectingPoint(self, pos, pdir):
        """Find intersection of line from pos along pdir with detector.
           Assumes flats are small perturbations, finding the phi at
           intersection with a circle, then looking for flats at that phi.
        """
        rdet = self["Det"].rmax

        # Determinant method courtesy of Wolfram MathWorld
        # https://mathworld.wolfram.com/Circle-LineIntersection.html
        dx,dy = pdir					# dx=x2-x1, dy=y2-y1
        dr = sqrt(dx**2+dy**2)
        Det = pos[0]*(pos[1]+dy) - (pos[0]+dx)*pos[1]   # Det = x1y2 - x2y1
        disc = rdet**2 * dr**2 - Det**2

        if (disc == 0.):	# Single solution at tangent point
            edgePos = (Det*dy/dr**2, -Det*dx/dr**2)
        elif (disc > 0):	# Two solutions, pick the closest
            x1 = (Det*dy + copysign(1.,dy)*dx*sqrt(disc))/dr**2
            y1 = (-Det*dx + abs(dy)*sqrt(disc))/dr**2
            dist1 = sqrt((x1-pos[0])**2+(y1-pos[1])**2)
            
            x2 = (Det*dy - copysign(1.,dy)*dx*sqrt(disc))/dr**2
            y2 = (-Det*dx - abs(dy)*sqrt(disc))/dr**2
            dist2 = sqrt((x2-pos[0])**2+(y2-pos[1])**2)

            edgePos = (x1,y1) if (dist1<dist2) else (x2,y2)
        else:			# No solution, point of closest approach
            _,dirphi = ChannelShape.cart2pol(pdir)
            _,posphi = ChannelShape.cart2pol(pos)

            dirphi  = ChannelShape.phinorm(dirphi)
            antiphi = ChannelShape.phinorm(posphi+pi)

            # If 0 < dirphi < pi+posphi, 0 < perp < pi, else -pi < perp < 0
            dist = Det/dr           
            if (ChannelShape.inPhi(dirphi, 0., antiphi)):
                outperp = dirphi - pi/2.
            elif (ChannelShape.inPhi(dirphi, 0., antiphi)):
                outperp = dirphi + pi/2.
            else:
                print("_getIntersectingPoint failed completely!")
                print("pos ",pos," posphi ",posphi," dir",dir," dirphi",dirphi)
                print("antiphi ",antiphi," dist ",dist)

            edgePos = ChannelShape.pol2cart((dist+rdet, outperp))

        # Check if phi "interesection" is at a flat, adjust r accordingly
        _,phi = ChannelShape.cart2pol(edgePos)
        for flat in self["Det"].flats:
            # NOTE: "Det" only has flats at rmax, so this will be unique
            if (flat.inPhi(phi)):
                edgePos = ChannelShape.pol2cart((flat.rAtPhi(phi), phi))
        
        return edgePos

    # Internal functions to read channel data from ROOT file (see Load() above)
    def _loadTTree(self, dmcfile, detnum=0):
        """Open specified DMC/SuperSim ROOT file and extract channel data for
           detector <detnum>.
           NOTE:  This version is using plain pyROOT to access TFile and TTree
           interfaces directly.  This should be replaced."""
 
        geom = None
        import ROOT     # Avoids having ROOT intercept the command line
        tfile = ROOT.TFile(dmcfile)
        for det in tfile.Get("G4SettingsInfoDir/Geometry"):
            if (det.DetNum == detnum):
                geom = det
                break

        if geom is None:
            print(f"ERROR Loading Geometry from {dmcfile}.")
            return self

        # In _loadRDF(), the branches are '["Branch"][0]'
        self.detName = geom.DetName
        self.detType = geom.DetType
        
        tomm = 1e3    # Convert m to mm
        
        for ich in range(geom.Channels):
            chantype = geom.ChanType[ich] if hasattr(geom,'ChanType') else None

            # FIXME: Need to separate charge from phonon channels
            ch = ChannelShape(geom.ChanName[ich], chantype,
                              geom.ChanRMin[ich]*tomm,
                              geom.ChanRMax[ich]*tomm,
                              radians(geom.ChanPhiMin[ich]),
                              radians(geom.ChanPhiMax[ich]),
                              geom.ChanSide[ich],
                              geom.ChanxRMin[ich]*tomm,
                              geom.ChanxRMax[ich]*tomm,
                              geom.ChanyRMin[ich]*tomm,
                              geom.ChanyRMax[ich]*tomm)
            self[ch.name] = ch

        # Full detector is cylindrical with flats to align with X and Y
        det = ChannelShape("Det", None, 0., geom.Radius*tomm, -pi, pi, +1,
                           0., geom.Axis1Len*tomm/2.,
                           0., geom.Axis2Len*tomm/2.)
        
        for ifl in range(geom.ExtraFlats):
            det.addFlat(geom.FlatR[ifl]*tomm, geom.FlatPhi[ifl])

        self[det.name] = det
        return self
        
    def _loadRDF(self, dmcfile, detnum=0):
        """Open specified DMC/SuperSim ROOT file and extract channel data for
           detector <detnum>.
           NOTE:  This version uses CATs' CDataFrame for access to the data
           via a Numpy dictionary."""
        from cats.cdataframe import CDataFrame

        # Load entire Geometry TTree for one detector; most branches are
        # channel dimensions
        ### FIXME: CDataFrame doesn't convert arrays of strings!
        geom = CDataFrame("G4SettingsInfoDir/Geometry",dmcfile)\
                   .Filter(f"DetNum=={detnum}").AsNumpy()
        self.detName = geom["DetName"][0]
        self.detType = geom["DetType"][0]
        
        tomm = 1e3    # Convert m to mm
        
        for ich in range(geom["Channels"][0]):
            chantype = geom["ChanType"][0][ich] if "ChanType" in geom else None

            # FIXME: Need to separate charge from phonon channels
            ch = ChannelShape(geom["ChanName"][0][ich], chantype, 
                              geom["ChanRMin"][0][ich]*tomm,
                              geom["ChanRMax"][0][ich]*tomm,
                              radians(geom["ChanPhiMin"][0][ich]),
                              radians(geom["ChanPhiMax"][0][ich]),
                              geom["ChanSide"][0][ich],
                              geom["ChanxRMin"][0][ich]*tomm,
                              geom["ChanxRMax"][0][ich]*tomm,
                              geom["ChanyRMin"][0][ich]*tomm,
                              geom["ChanyRMax"][0][ich]*tomm)
            self[ch.name] = ch

        # Full detector is cylindrical with flats to align with X and Y
        det = ChannelShape("Det", None, 0., geom["Radius"][0]*tomm, -pi, pi, +1,
                           0., geom["Axis1Len"][0]*tomm/2.,
                           0., geom["Axis2Len"][0]*tomm/2.)
        
        for ifl in range(geom["ExtraFlats"][0]):
            det.addFlat(geom["FlatR"][0][ifl]*tomm, geom["FlatPhi"][0][ifl])

        self[det.name] = det
        return self

    # Scan through defined channels and create TES-only and FET-only subsets
    def FillGroups(self):
        """Scans the dictionary of all channels and creates DetectorShape data
           members for just TES or just FET channels, based on chan.chantype.
           If channel type wasn't set, these objects will be empty."""

        self.TES = DetectorShape(self.detType)
        self.FET = DetectorShape(self.detType)
        for chan in self.values():
            if (chan.isTES()): self.TES.Add(chan)
            if (chan.isFET()): self.FET.Add(chan)

        return self

    # Convert between detector type codes and name strings
    _dnames = {4:"oZIP", 10:"iZIP1", 12:"iZIP2", 11:"iZIP5", 21:"CDMSlite1",
               22:"CDMSlite2", 700:"iZIP7", 701:"iZIP7Si", 710:"HV100mm",
               711:"HV100mmSi"}
    _dtypes = None			# Will be filled by getDetType()
    
    @staticmethod
    def getDetName(dtype):
        return DetectorShape._dnames[dtype] if isinstance(dtype, (int,float)) else None

    @staticmethod
    def getDetType(dname):
        if DetectorShape._dtypes is None:
            DetectorShape._dtypes = {DetectorShape._dnames[t]:t for t in DetectorShape._dnames}

        return DetectorShape._dtypes[dname] if isinstance(dname, str) else None

### END OF DetectorShape ###
