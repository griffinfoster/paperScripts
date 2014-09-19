#!/usr/bin/env python
"""
Apply a D-Projection correction to a CASA image or set of FITS images
and output full Stokes image
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
    o.add_option('-d','--dproj', dest='dprojFile', default=None,
        help='Numpy de-projection file to apply, must be of the same (l,m) dimensions as the FITS images')
    o.add_option('-f','--freqs',dest='freqs',default=None,
        help='Comma separated list of frequency IDs to apply de-projection to, default: all')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.dprojFile is None:
        print 'Error: Missing de-projection file'
        exit(1)
    print 'Loading',opts.dprojFile
    dproj=n.load(opts.dprojFile)

    nfiles=len(args)
    if opts.freqs is None: freqs=None
    else: freqs=map(int,opts.freqs.split(','))

    xximg,xxhdr=fitsUtil.readFITS(opts.xximg,hdr=True)
    xyrimg,xyrhdr=fitsUtil.readFITS(opts.xyrimg,hdr=True)
    xyiimg,xyihdr=fitsUtil.readFITS(opts.xyiimg,hdr=True)
    yyimg,yyhdr=fitsUtil.readFITS(opts.yyimg,hdr=True)

    im=n.zeros((4,xximg.shape[1],xximg.shape[2],xximg.shape[3]),dtype=n.complex64)
    im[cplxDict['XX*']]=xximg
    im[cplxDict['XY*']]=.5*(xyrimg+(1j*xyiimg))
    im[cplxDict['YX*']]=.5*(xyrimg+(1j*xyiimg)).conjugate()
    im[cplxDict['YY*']]=yyimg

    im=n.transpose(im,(2,3,0,1)) #put the image into a useful shape (DEC pixel, RA pixel, stokes-vector, freq)
    cplxIm=n.zeros_like(im)
    if freqs is None: freqIDs=range(im.shape[3])
    else: freqIDs=freqs
    #this is a good ol' hack to handle the multidimensional matrix multiplication
    #for i in range(4): cplxIm[:,:,i]=n.sum(dproj[:,:,i,:,n.newaxis]*im,axis=2)

    #apply inverted d-projection matrix and convert to Stokes Parameters
    for f in freqIDs:
        #for i in range(4): cplxIm[:,:,i,f]=n.sum(dproj[:,:,i,:]*im[:,:,:,f],axis=2)
        for i in range(4): cplxIm[:,:,i,f]=n.sum(dproj[:,:,i,:]*im[:,:,:,f],axis=2)

    outIm=n.zeros_like(im, dtype=n.float32)
    outIm[:,:,0]=cplxIm[:,:,0].real #I
    outIm[:,:,1]=cplxIm[:,:,1].real #Q
    outIm[:,:,2]=cplxIm[:,:,2].real #U
    outIm[:,:,3]=cplxIm[:,:,3].imag #V

    outIm=n.transpose(outIm,(2,3,0,1)) #put back into the FITS shape

    #write outIm to new FITS file
    path,srcFn=os.path.split(os.path.realpath(opts.xximg))
    dstFn=opts.outfn
    print 'Creating new FITS file:',os.path.join(path,dstFn)
    shutil.copyfile(os.path.join(path,srcFn), os.path.join(path,dstFn))
    fh=pf.open(os.path.join(path,dstFn),mode='update')
    hdr=fh[0].header
    hdr.update('CRVAL4',1) #make first Stokes ID Stokes I
    hdr.update('CRVAL4',4) #fix header to 4 Stokes values
    fh[0].data=outIm
    fh.flush()
    fh.close()

