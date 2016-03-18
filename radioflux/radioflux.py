#!/usr/bin/env python

import pyregion
import scipy.stats
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys
import warnings

def flatten(f,channel=0,freqaxis=0):
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

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)
        
# slice=(0,)*(naxis-2)+(np.s_[:],)*2
    return header,f[0].data[slice]

class RadioError(Exception):
    """Base class for exceptions in this module."""
    pass

class radiomap:
    """ Process a fits file as though it were a radio map, calculating beam areas etc """
    def __init__(self, fitsfile,verbose=False):
        # Catch warnings to avoid datfix errors
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gfactor=2.0*np.sqrt(2.0*np.log(2.0))
            self.f=fitsfile[0]
            self.prhd=fitsfile[0].header

            # Get units and resolution
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
                except KeyError:
                    history=None
                if history is not None:
                    for line in history:
                        if 'HISTORY' in line:
                            continue # stops it finding nested history
                        if 'CLEAN BMAJ' in line:
                            bits=line.split()
                            self.bmaj=float(bits[3])
                            self.bmin=float(bits[5])
                                
            if self.bmaj is None:
                raise RadioError('No beam information found')

            w=wcs.WCS(self.prhd)
            cd1=-w.wcs.cdelt[0]
            cd2=w.wcs.cdelt[1]
            if ((cd1-cd2)/cd1)>1.0001 and ((self.bmaj-self.bmin)/self.bmin)>1.0001:
                raise RadioError('Pixels are not square (%g, %g) and beam is elliptical' % (cd1, cd2))

            self.bmaj/=cd1
            self.bmin/=cd2
            if verbose:
                print 'beam is',self.bmaj,'by',self.bmin,'pixels'

            self.area=2.0*np.pi*(self.bmaj*self.bmin)/(gfactor*gfactor)
            if verbose:
                print 'beam area is',self.area,'pixels'

            # Now check what sort of a map we have
            naxis=len(fitsfile[0].data.shape)
            if verbose: print 'We have',naxis,'axes'
            self.cube=False
            if naxis<2 or naxis>4:
                raise RadioError('Too many or too few axes to proceed (%i)' % naxis)
            if naxis>2:
                # a cube, what sort?
                frequency=0
                self.cube=True
                freqaxis=-1
                stokesaxis=-1
                for i in range(3,naxis+1):
                    ctype=self.prhd.get('CTYPE%i' % i)
                    if 'FREQ' in ctype:
                        freqaxis=i
                    elif 'STOKES' in ctype:
                        stokesaxis=i
                    else:
                        raise RadioError('Unknown CTYPE %i = %s' % (i,ctype))
                if verbose:
                    print 'This is a cube with freq axis %i and Stokes axis %i' % (freqaxis, stokesaxis)
                if stokesaxis>0:
                    nstokes=self.prhd.get('NAXIS%i' % stokesaxis)
                    if nstokes>1:
                        raise RadioError('Multiple Stokes parameters present, not handled')
                if freqaxis>0:
                    nchans=self.prhd.get('NAXIS%i' % freqaxis)
                    if verbose:
                        print 'There are %i channel(s)' % nchans
                    self.nchans=nchans
            else:
                self.nchans=1
                    

            # Various possibilities for the frequency. It's possible
            # that a bad (zero) value will be present, so keep
            # checking if one is found.

            if not(self.cube) or freqaxis<0:
                # frequency, if present, must be in another keyword
                frequency=self.prhd.get('RESTFRQ')
                if frequency is None or frequency==0:
                    frequency=self.prhd.get('RESTFREQ')
                if frequency is None or frequency==0:
                    frequency=self.prhd.get('FREQ')
                if frequency is None or frequency==0:
                    # It seems some maps present with a FREQ ctype
                    # even if they don't have the appropriate axes!
                    # The mind boggles.
                    for i in range(5):
                        type_s=self.prhd.get('CTYPE%i' % i)
                        if type_s is not None and type_s[0:4]=='FREQ':
                            frequency=self.prhd.get('CRVAL%i' % i)
                self.frq=[frequency]
                # now if there _are_ extra headers, get rid of them so pyregion WCS can work
                for i in range(3,5):
                    self.prhd.remove('CTYPE%i' %i)
                    self.prhd.remove('CRVAL%i' %i)
                    self.prhd.remove('CDELT%i' %i)
                    self.prhd.remove('CRPIX%i' %i)
                    self.prhd.remove('CROTA%i' %i)
                self.headers=[self.prhd]
                self.d=[fitsfile[0].data]
            else:
                # if this is a cube, frequency/ies should be in freq header
                basefreq=self.prhd.get('CRVAL%i' % freqaxis)
                deltafreq=self.prhd.get('CDELT%i' % freqaxis)
                self.frq=[basefreq+deltafreq*i for i in range(nchans)]
                self.d=[]
                self.headers=[]
                for i in range(nchans):
                    header,data=flatten(fitsfile,freqaxis=freqaxis,channel=i)
                    self.d.append(data)
                    self.headers.append(header)
            for i,f in enumerate(self.frq):
                if f is None:
                    print('Warning, can\'t get frequency %i -- set to zero' % i)
                    self.frq[i]=0
            if verbose:
                print 'Frequencies are',self.frq,'Hz'

#            self.fhead,self.d=flatten(fitsfile)

class applyregion:
    """ apply a region from pyregion to a radiomap """
    def __init__(self,rm,region,background=None,offsource=None):
        self.rms=[]
        self.flux=[]
        self.mean=[]
        self.error=[]
        bgval=0
        if background is not None:
            bgval=background
        for i,d in enumerate(rm.d):
            mask=region.get_mask(hdu=rm.f,shape=np.shape(d))
            self.pixels=np.sum(mask)
            data=np.extract(mask,d)
            data-=bgval
            self.rms.append(scipy.stats.nanstd(data))
            self.flux.append(data[np.logical_not(np.isnan(data))].sum()/rm.area)
            self.mean.append(scipy.stats.nanmean(data))
            if offsource is not None:
                self.error.append(offsource[i]*np.sqrt(self.pixels/rm.area))

# Command-line running

def printflux(filename,rm,region,noise,bgsub,background=0,label=''):
    if bgsub:
        fg=applyregion(rm,region,offsource=noise,background=background)
    else:
        fg=applyregion(rm,region,offsource=noise)

    for i in range(rm.nchans):
        freq=rm.frq[i]
        if noise:
            print filename,label,'%8.4g %10.6f %10.6f' % (freq,fg.flux[i],fg.error[i])
        else:
            print filename,label,'%8.4g %10.6f' % (freq,fg.flux[i])

def flux_for_files(files,fgr,bgr=None,individual=False,bgsub=False,action=printflux,verbose=False):
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
        rm=radiomap(fitsfile,verbose=verbose)
        if bgr:
            bg_ir=pyregion.open(bgr).as_imagecoord(rm.headers[0])
            bg=applyregion(rm,bg_ir)
            noise=bg.rms
            background=bg.mean
        else:
            if bgsub:
                raise RadioError('Background subtraction requested but no bg region')
            noise=None
            background=None

        fg_ir=pyregion.open(fgr).as_imagecoord(rm.headers[0])

        if individual:
            for n,reg in enumerate(fg_ir):
                fg=pyregion.ShapeList([reg])
                r=action(filename,rm,fg,noise,bgsub,background,label=n+1)
        else:
            r=action(filename,rm,fg_ir,noise,bgsub,background)

        fitsfile.close()
        return r
        
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
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true',default=False,help='Be verbose')

    args = parser.parse_args()

    flux_for_files(args.files,args.fgr,args.bgr,args.indiv,args.bgsub,verbose=args.verbose)
