#!/usr/bin/env python
"""
Convert a NPY file to a FITS file using a FITS file as a template
"""

import sys,os,shutil
import pyfits as pf
import numpy as np

import fitsUtil

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] NPY_FILE')
    o.set_description(__doc__)
    o.add_option('-t','--template',dest='template',default=None,
        help='Template fits file')
    o.add_option('-o','--out',dest='outfn',default=None,
        help='Output FITS filename, default: inpuy NPY filename with .fits extension')
    opts, args = o.parse_args(sys.argv[1:])

    npyFn=args[0]
    imCube=np.load(npyFn)
    print 'NPY shape:', imCube.shape

    imFits,imHdr=fitsUtil.readFITS(opts.template,hdr=True)
    print 'FITS shape:', imFits.shape

    #write imCube to FITS file
    path,srcFn=os.path.split(os.path.realpath(npyFn))
    if opts.outfn is None:
        dstFn=srcFn.split('.npy')[0]+'.fits'
    else: dstFn=opts.outfn
    print 'Creating new FITS file:',os.path.join(path,dstFn)
    tpath,tempFn=os.path.split(os.path.realpath(opts.template))
    shutil.copyfile(os.path.join(tpath,tempFn), os.path.join(path,dstFn))
    fh=pf.open(os.path.join(path,dstFn),mode='update')
    #hdr=fh[0].header
    #hdr.update('CRVAL4',1)
    fh[0].data=imCube
    fh.flush()
    fh.close()

