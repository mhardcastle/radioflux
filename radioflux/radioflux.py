#!/usr/bin/env python

import pyregion
import scipy.stats
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys
import warnings

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return f[0].header,f[0].data

    w = wcs.WCS(f[0].header)
    wn=wcs.WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header,f[0].data[slice]

class RadioError(Exception):
    """Base class for exceptions in this module."""
    pass

class radiomap:
    """ Process a fits file as though it were a radio map, calculating beam areas etc """
    def __init__(self, fitsfile,**extras):
        # Catch warnings to avoid datfix errors
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gfactor=2.0*np.sqrt(2.0*np.log(2.0))
            verbose=('verbose' in extras)
            self.f=fitsfile[0]
            self.prhd=fitsfile[0].header
            self.units=self.prhd.get('BUNIT')
            if self.units is None:
                self.units=self.prhd.get('UNIT')
            if self.units!='JY/BEAM' and self.units!='Jy/beam':
                print 'Warning: units are',self.units,'but code expects JY/BEAM'
            self.bmaj=self.prhd.get('BMAJ')
            self.bmin=self.prhd.get('BMIN')
            if self.bmaj is None:
                # Try RESOL1 and RESOL2
                self.bmaj=self.prhd.get('RESOL1')
                self.bmin=self.prhd.get('RESOL2')
            if self.bmaj is None:
                if verbose:
                    print 'Can\'t find BMAJ in headers, checking history'
                try:
                    history=self.prhd['HISTORY']
                    for line in history:
                        if 'HISTORY' in line:
                            continue # stops it finding nested history
                        if 'CLEAN BMAJ' in line:
                            bits=line.split()
                            self.bmaj=float(bits[3])
                            self.bmin=float(bits[5])
                except KeyError:
                    pass
                                
            if self.bmaj is None:
                raise RadioError('No beam information found')

            # Various possibilities for the frequency. It's possible
            # that a bad (zero) value will be present, so keep
            # checking if one is found.

            self.frq=self.prhd.get('RESTFRQ')
            if self.frq is None or self.frq==0:
                self.frq=self.prhd.get('RESTFREQ')
            if self.frq is None or self.frq==0:
                self.frq=self.prhd.get('FREQ')
            if self.frq is None or self.frq==0:
                i=1
                while True:
                    keyword='CTYPE'+str(i)
                    ctype=self.prhd.get(keyword)
                    if ctype is None:
                        break
                    if ctype=='FREQ':
                        self.frq=self.prhd.get('CRVAL'+str(i))
                        break
                    i+=1

            if self.frq is None:
                print('Warning, can\'t get frequency -- set to zero')
                self.frq=0

            w=wcs.WCS(self.prhd)
            cd1=-w.wcs.cdelt[0]
            cd2=w.wcs.cdelt[1]
            if ((cd1-cd2)/cd1)>1.0001 and ((bmaj-bmin)/bmin)>1.0001:
                raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

            self.bmaj/=cd1
            self.bmin/=cd2
            if verbose:
                print 'beam is',self.bmaj,'by',self.bmin,'pixels'

            self.area=2.0*np.pi*(self.bmaj*self.bmin)/(gfactor*gfactor)
            if verbose:
                print 'beam area is',self.area,'pixels'

            self.fhead,self.d=flatten(fitsfile)

class applyregion:
    """ apply a region from pyregion to a radiomap """
    def __init__(self,rm,region,**extras):
        bgval=0;
        if 'background' in extras:
            bgval=extras['background']
        mask=region.get_mask(hdu=rm.f,shape=np.shape(rm.d))
        self.pixels=np.sum(mask)
        data=np.extract(mask,rm.d)
        data-=bgval
        self.rms=scipy.stats.nanstd(data)
        self.flux=data[np.logical_not(np.isnan(data))].sum()/rm.area
        self.mean=scipy.stats.nanmean(data)
        if 'offsource' in extras:
            self.error=extras['offsource']*np.sqrt(self.pixels/rm.area)

# Command-line running

def printflux(filename,rm,region,noise,bgsub,background=0,label=''):
    if bgsub:
        fg=applyregion(rm,region,offsource=noise,background=background)
    else:
        fg=applyregion(rm,region,offsource=noise)

    if noise:
        print filename,label,'%g' % rm.frq,fg.flux,fg.error
    else:
        print filename,label,'%g' % rm.frq,fg.flux

def flux_for_files(files,fgr,bgr=None,individual=False,bgsub=False,action=printflux):
    """Determine the flux in a region file for a set of files. This is the
    default action for the code called on the command line, but
    may be useful to other code as well.

    Keyword arguments:
    files -- list of files (mandatory)
    fdr -- foreground region name (mandatory)
    bgr -- background region name (optional)
    individual -- separate region into individual sub-regions
    bgsub -- subtract background
    action -- what to do once fluxes are measured: allows a user-defined action
              which must be a drop-in replacement for printflux
    """

    for filename in files:
        fitsfile=fits.open(filename)
        rm=radiomap(fitsfile)
        if bgr:
            bg_ir=pyregion.open(bgr).as_imagecoord(rm.fhead)
            bg=applyregion(rm,bg_ir)
            noise=bg.rms
            background=bg.mean
        else:
            if bgsub:
                raise RadioError('Background subtraction requested but no bg region')
            noise=0
            background=0

        fg_ir=pyregion.open(fgr).as_imagecoord(rm.fhead)

        if individual:
            for n,reg in enumerate(fg_ir):
                fg=pyregion.ShapeList([reg])
                action(filename,rm,fg,noise,bgsub,background,label=n+1)
        else:
            action(filename,rm,fg_ir,noise,bgsub,background)

        fitsfile.close()

        
if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Measure fluxes from FITS files.')
    parser.add_argument('files', metavar='FILE', nargs='+',
                        help='FITS files to process')
    parser.add_argument('-f','--foreground', dest='fgr', action='store',default='ds9.reg',help='Foreground region file to use')
    parser.add_argument('-b','--background', dest='bgr', action='store',default='',help='Background region file to use')
    parser.add_argument('-i','--individual', dest='indiv', action='store_true',default=False,help='Break composite region file into individual regions')
    parser.add_argument('-s','--subtract', dest='bgsub', action='store_true',default=False,help='Subtract background')

    args = parser.parse_args()

    flux_for_files(args.files,args.fgr,args.bgr,args.indiv,args.bgsub)
