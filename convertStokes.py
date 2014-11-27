#!/usr/bin/env python
"""
Convert WSCLEAN output FITS files to Stokes
"""

import sys,os,shutil
import pyfits as pf
import numpy as n
import pywcs
import pyrap.tables as pt

import fitsUtil

stokesDict={'I':0,'Q':1,'U':2,'V':3}
cplxDict={'XX*':0,'XY*':1,'YX*':2,'YY*':3}

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options]')
    o.set_description(__doc__)
    o.add_option('--xx',dest='xximg',default=None,
        help='XX FITS image, required if not using CASA image MS')
    o.add_option('--xyr',dest='xyrimg',default=None,
        help='XY real FITS image, required if not using CASA image MS')
    o.add_option('--xyi',dest='xyiimg',default=None,
        help='XY imaginary FITS image, required if not using CASA image MS')
    o.add_option('--yy',dest='yyimg',default=None,
        help='YY FITS image, required if not using CASA image MS')
    o.add_option('-o','--out',dest='outfn',default=None,
        help='Output FITS filename')
    opts, args = o.parse_args(sys.argv[1:])

    xximg,xxhdr=fitsUtil.readFITS(opts.xximg,hdr=True)
    xyrimg,xyrhdr=fitsUtil.readFITS(opts.xyrimg,hdr=True)
    xyiimg,xyihdr=fitsUtil.readFITS(opts.xyiimg,hdr=True)
    yyimg,yyhdr=fitsUtil.readFITS(opts.yyimg,hdr=True)

    im=n.zeros((4,xximg.shape[1],xximg.shape[2],xximg.shape[3]),dtype=n.complex64)
    im[cplxDict['XX*']]=xximg
    im[cplxDict['XY*']]=.5*(xyrimg+(1j*xyiimg))
    im[cplxDict['YX*']]=.5*(xyrimg+(1j*xyiimg)).conjugate()
    im[cplxDict['YY*']]=yyimg

    im=n.transpose(im,(2,3,0,1)) #put the image into a useful shape (x pixel, y pixel, stokes-vector, freq)

    #convert to Stokes Parameters
    outIm=n.zeros_like(im, dtype=n.float32)
    outIm[:,:,0]=.5*(im[:,:,3]+im[:,:,0]).real #I
    outIm[:,:,1]=.5*(im[:,:,3]-im[:,:,0]).real #Q
    outIm[:,:,2]=.5*(im[:,:,2]+im[:,:,1]).real #U
    outIm[:,:,3]=.5*(im[:,:,2]-im[:,:,1]).imag #V

    outIm=n.transpose(outIm,(2,3,0,1)) #put the ouput image back into the correct order

    #write outIm to new FITS file
    path,srcFn=os.path.split(os.path.realpath(opts.xximg))
    dstFn=opts.outfn
    print 'Creating new FITS file:',os.path.join(path,dstFn)
    shutil.copyfile(os.path.join(path,srcFn), os.path.join(path,dstFn))
    fh=pf.open(os.path.join(path,dstFn),mode='update')
    hdr=fh[0].header
    hdr.update('CRVAL4',1)
    fh[0].data=outIm
    fh.flush()
    fh.close()

