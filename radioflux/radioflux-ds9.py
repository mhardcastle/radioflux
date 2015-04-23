#!/usr/bin/env python

import pyregion
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys
from radioflux import *
import argparse

parser = argparse.ArgumentParser(description='ds9 plugin to measure fluxes from FITS files.')
parser.add_argument('-s','--subtract', dest='bgsub', action='store_const', const=1,default=0,help='Subtract background')
args = parser.parse_args()

filename=sys.stdin.readline().rstrip()
f_region=sys.stdin.readline().rstrip()
b_region=sys.stdin.readline().rstrip()

bgsub=args.bgsub

print '-------------------------------------------------------------'

print 'Filename is',filename
print 'FG region is',f_region
print 'BG region is',b_region

f=fits.open(filename)
rm=radiomap(f,verbose=1)

print 'Frequency is %g Hz' % rm.frq 

# now apply regions

if b_region:
    bg_ir=pyregion.parse(b_region).as_imagecoord(rm.prhd)
    bg=applyregion(rm,bg_ir)
    print 'Pixels in background region',bg.pixels
    print 'Background rms is',bg.rms,'Jy/beam'
    print 'Background mean is',bg.mean,'Jy/beam'
    noise=bg.rms
else:
    if bgsub:
        raise Error('Background subtraction requested but no bg region')
    noise=0

fg_ir=pyregion.parse(f_region).as_imagecoord(rm.prhd)
if bgsub:
    fg=applyregion(rm,fg_ir,offsource=noise,background=bg.mean)
else:
    fg=applyregion(rm,fg_ir,offsource=noise)

print 'Pixels passed by mask',fg.pixels

if noise:
    print 'Region flux is',fg.flux,'+/-',fg.error,'Jy'
else:
    print 'Region flux is',fg.flux,'Jy'
