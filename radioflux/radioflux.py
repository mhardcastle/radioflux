#!/usr/bin/python

import pyregion
from astropy.io import fits
from astropy import wcs
import numpy as np
import sys
import warnings

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
            if self.units!='JY/BEAM':
                print 'Warning: units are',self.units,'but code expects JY/BEAM'
            self.bmaj=self.prhd.get('BMAJ')
            self.bmin=self.prhd.get('BMIN')
            if self.bmaj==None:
                if verbose:
                    print 'Can\'t find BMAJ in headers, checking history'
                history=self.prhd['HISTORY']
                for line in history:
                    if 'CLEAN BMAJ' in line:
                        bits=line.split()
                        self.bmaj=float(bits[3])
                        self.bmin=float(bits[5])
                                
            if self.bmaj==None:
                raise RadioError('No beam information found')

# Various possibilities for the frequency
            self.frq=self.prhd.get('RESTFRQ')
            if self.frq==None:
                self.frq=self.prhd.get('RESTFREQ')
                if self.frq==None:
                    i=1
                    while True:
                        keyword='CTYPE'+str(i)
                        try:
                            ctype=self.prhd.get(keyword)
                        except KeyError:
                            break
                        if ctype=='FREQ':
                            self.frq=self.prhd.get('CRVAL'+str(i))
                            break
                        i+=1
                    else:
                        print('Warning, can\'t get frequency -- set to zero')
                        self.frq=0

            w=wcs.WCS(self.prhd)
            cd1=-w.wcs.cdelt[0]
            cd2=w.wcs.cdelt[1]
            if ((cd1-cd2)/cd1)>1.0001:
                raise RadioError('Pixels are not square (%g, %g)' % (cd1, cd2))

            self.bmaj/=cd1
            self.bmin/=cd1
            if verbose:
                print 'beam is',self.bmaj,'by',self.bmin,'pixels'

            self.area=2.0*np.pi*(self.bmaj*self.bmin)/(gfactor*gfactor)
            if verbose:
                print 'beam area is',self.area,'pixels'

# now a bit of a hack to extract the image. Some AIPS-type images are data cubes with two null axes; other FITS images are just 2D images.

            axes=len(np.shape(self.f.data))
            if axes==4:
                self.d=self.f.data[0,0]
            elif axes==3:
                self.d=self.f.data[0]
            else:
                self.d=self.f.data

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
        self.rms=data.std()
        self.flux=data.sum()/rm.area
        self.mean=data.mean()
        if 'offsource' in extras:
            self.error=extras['offsource']*np.sqrt(self.pixels/rm.area)

# Command-line running

def printflux(rm,region,noise,bgsub,background=0,label=''):
    if bgsub:
        fg=applyregion(rm,region,offsource=noise,background=background)
    else:
        fg=applyregion(rm,region,offsource=noise)

    if noise:
        print filename,label,'%g' % rm.frq,fg.flux,fg.error
    else:
        print filename,label,'%g' % rm.frq,fg.flux


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Measure fluxes from FITS files.')
    parser.add_argument('files', metavar='FILE', nargs='+',
                        help='FITS files to process')
    parser.add_argument('-f','--foreground', dest='fgr', action='store',default='ds9.reg',help='Foreground region file to use')
    parser.add_argument('-b','--background', dest='bgr', action='store',default='',help='Background region file to use')
    parser.add_argument('-i','--individual', dest='indiv', action='store_const', const=1,default=0,help='Break composite region file into individual regions')
    parser.add_argument('-s','--subtract', dest='bgsub', action='store_const', const=1,default=0,help='Subtract background')

    args = parser.parse_args()

    c=0
    for filename in args.files:
        fitsfile=fits.open(filename)
        rm=radiomap(fitsfile)
        if args.bgr:
            bg_ir=pyregion.open(args.bgr).as_imagecoord(rm.prhd)
            bg=applyregion(rm,bg_ir)
            noise=bg.rms
            background=bg.mean
        else:
            if args.bgsub:
                raise RadioError('Background subtraction requested but no bg region')
            noise=0
            background=0

        fg_ir=pyregion.open(args.fgr).as_imagecoord(rm.prhd)

        if args.indiv:
            for n,reg in enumerate(fg_ir):
                fg=pyregion.ShapeList([reg])
                printflux(rm,fg,noise,args.bgsub,background,label=n+1)
        else:
            printflux(rm,fg_ir,noise,args.bgsub,background)

#        if args.bgsub:
#            fg=applyregion(rm,fg_ir,offsource=noise,background=bg.mean)
#        else:
#            fg=applyregion(rm,fg_ir,offsource=noise)
#        
#        if noise:
#            print filename,'%g' % rm.frq,fg.flux,fg.error
#        else:
#            print filename,'%g' % rm.frq,fg.flux
        fitsfile.close()

