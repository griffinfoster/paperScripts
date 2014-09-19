#!/usr/bin/env python
"""
Apply a beam correction to a Stokes FITS file, assumes the same pixel size and resolution
"""

import sys,os,shutil
import pyfits as pf
import numpy as np
import pywcs
import pyrap.tables as pt

import fitsUtil

stokesDict={'I':0,'Q':1,'U':2,'V':3}
cplxDict={'XX*':0,'XY*':1,'YX*':2,'YY*':3}

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -b BEAM_FITS STOKES_FITS')
    o.set_description(__doc__)
    o.add_option('-b','--beam', dest='beamFile', default=None,
        help='Interpolated Stokes beam FITS file, same number of pixels and resolution as the input FITS file')
    o.add_option('-i','--invert',dest='invert',action='store_true',
        help='Do not apply inverted beam, apply the beam directly')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.beamFile is None:
        print 'Error: Missing beam FITS'
        exit(1)
    print 'Loading',opts.beamFile
    beam,beamHdr,beamAxis=fitsUtil.readFITS(opts.beamFile,hdr=True,axis=True)

    for fid,fn in enumerate(args):
        print 'Processing',fn
        im,fHdr,fAxis=fitsUtil.readFITS(fn,hdr=True,axis=True)
        if im.shape==beam.shape:
            if opts.invert: outIm=im*beam
            else: outIm=im*(1./beam)
        else:
            print 'ERROR: beam shape different than image shape'
            exit(1)

        #write outIm to new FITS file
        path,srcFn=os.path.split(os.path.realpath(fn))
        dstFn=srcFn.split('.fits')[0]+'.beam.fits'
        print 'Creating new FITS file:',os.path.join(path,dstFn)
        shutil.copyfile(os.path.join(path,srcFn), os.path.join(path,dstFn))
        fh=pf.open(os.path.join(path,dstFn),mode='update')
        fh[0].data=outIm
        fh.flush()
        fh.close()

