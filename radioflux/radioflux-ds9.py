#!/usr/bin/env python

import pyregion
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys
from radioflux import *
import argparse
import re

parser = argparse.ArgumentParser(description='ds9 plugin to measure fluxes from FITS files.')
parser.add_argument('-s','--subtract', dest='bgsub', action='store_const', const=1,default=0,help='Subtract background')
args = parser.parse_args()

filename=sys.stdin.readline().rstrip()
# discard any ds9 qualifiers, since we can't use them
filename=re.sub(r'\[.*\]','',filename)
f_region=sys.stdin.readline().rstrip()
b_region=sys.stdin.readline().rstrip()

bgsub=args.bgsub

print '-------------------------------------------------------------'

print 'Filename is',filename
print 'FG region is <<'+f_region+'>>'
print 'BG region is <<'+b_region+'>>'

if f_region=="":
    print "FATAL ERROR: code can't work without at least one foreground region."
    sys.exit()

try:
    f=fits.open(filename)
except IOError as e:
    print 'FATAL ERROR: ',e
    sys.exit()

try:
    rm=radiomap(f,verbose=True)
except RadioError as e:
    print 'FATAL ERROR: ',e
    sys.exit()

#print 'Frequency is %g Hz' % rm.frq 

# now apply regions

if b_region:
    bg_ir=pyregion.parse(b_region).as_imagecoord(rm.headers[0])
    bg=applyregion(rm,bg_ir)
    print 'Pixels in background region',bg.pixels
    for i in range(rm.nchans):
        print '%8.4g Hz Background rms is %f Jy/beam' % (rm.frq[i],bg.rms[i])
        print '             Background mean is',bg.mean[i],'Jy/beam'
    noise=bg.rms
else:
    if bgsub:
        raise Error('Background subtraction requested but no bg region')
    noise=None

fg_ir=pyregion.parse(f_region).as_imagecoord(rm.headers[0])
if bgsub:
    fg=applyregion(rm,fg_ir,offsource=noise,background=bg.mean)
else:
    fg=applyregion(rm,fg_ir,offsource=noise)

print 'Pixels in foreground region',fg.pixels
for i in range(rm.nchans):
    freq=rm.frq[i]
    if noise:
        print '%8.4g Hz Region flux is %f +/- %f Jy' % (freq,fg.flux[i],fg.error[i])
    else:
        print '%8.4g Hz Region flux is %f Jy' % (freq,fg.flux[i])
