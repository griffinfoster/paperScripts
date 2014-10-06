#!/usr/bin/env python
"""
Compute the IXR from a NPY array, assumes a known order of (complex[xx,xyr,xyi,yxr,yxi,yy], freq, ra, dec)
"""

import os,sys
import numpy as np

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] NPY_FILE')
    o.set_description(__doc__)
    o.add_option('-o','--out',dest='outfn',default=None,
        help='Output FITS filename, default: inpuy NPY filename with .fits extension')
    opts, args = o.parse_args(sys.argv[1:])

    npyFn=args[0]
    imCube=np.load(npyFn)
    print 'NPY shape:', imCube.shape

    #compute IXR
    kappaCube=np.zeros((1,imCube.shape[1],imCube.shape[2],imCube.shape[3]))
    #very slow set of loops
    for ll in range(imCube.shape[2]):
        for mm in range(imCube.shape[3]):
            for freq in range(imCube.shape[1]):
                U,s,V=np.linalg.svd(np.matrix([
                    [imCube[0,freq,ll,mm],imCube[1,freq,ll,mm]+1j*imCube[2,freq,ll,mm]],
                    [imCube[3,freq,ll,mm]+1j*imCube[4,freq,ll,mm],imCube[5,freq,ll,mm]]]))
                kappaCube[0,freq,ll,mm]=s[0]/s[1]
    ixrJ=((kappaCube-1)/(kappaCube+1))**2.
    ixrJdB=10.*np.log10(ixrJ)
    print ixrJdB.shape

    #write ixrJdB to FITS file
    path,srcFn=os.path.split(os.path.realpath(npyFn))
    if opts.outfn is None:
        dstFn=srcFn.split('.npy')[0]+'.ixr.npy'
    else: dstFn=opts.outfn
    print 'Creating new NPY file:',os.path.join(path,dstFn)
    np.save(os.path.join(path,dstFn),ixrJdB)

