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
    o.set_usage('%prog [options] -b BEAM_NPY STOKES_FITS')
    o.set_description(__doc__)
    o.add_option('-b','--beam', dest='beamFile', default=None,
        help='Interpolated Stokes beam NPY file, same number of pixels and resolution as the input FITS file')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.beamFile is None:
        print 'Error: Missing beam NPY'
        exit(1)
    print 'Loading',opts.beamFile
    invBeam=np.load(opts.beamFile)
    invStokesBeam=np.zeros((4,invBeam.shape[1],invBeam.shape[2],invBeam.shape[3]))

    for fid,fn in enumerate(args):
        print 'Processing',fn
        im,fHdr,fAxis=fitsUtil.readFITS(fn,hdr=True,axis=True)
        if im.shape==invStokesBeam.shape:
            #convert invBeam to Stokes
            invStokesBeam[0]=invBeam[0]+invBeam[5]
            invStokesBeam[1]=invBeam[0]-invBeam[5]
            invStokesBeam[2]=((invBeam[3]+1j*invBeam[4])+(invBeam[1]+1j*invBeam[2])).real
            invStokesBeam[3]=((invBeam[3]+1j*invBeam[4])-(invBeam[1]+1j*invBeam[2])).imag
        else:
            print 'ERROR: beam shape different than image shape'
            exit(1)

        #print im.shape
        #print invStokesBeam.shape
        #apply inverted beam to image cube
        beamCorrIm=invStokesBeam*im

        #write outIm to new FITS file
        path,srcFn=os.path.split(os.path.realpath(fn))
        dstFn=srcFn.split('.fits')[0]+'.beam.fits'
        print 'Creating new FITS file:',os.path.join(path,dstFn)
        shutil.copyfile(os.path.join(path,srcFn), os.path.join(path,dstFn))
        fh=pf.open(os.path.join(path,dstFn),mode='update')
        fh[0].data=beamCorrIm
        fh.flush()
        fh.close()

