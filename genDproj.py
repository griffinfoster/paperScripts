#!/usr/bin/env python
"""
From a FITS image generate a dipole projection correction factor for each pixel based on FITS header for resolution and extent
"""

import sys,os,shutil
import pyfits as pf
import numpy as n
import pywcs

import fitsUtil

def d2r(a):
    """Degrees to Radians"""
    return n.pi*a/180.

def r2d(a):
    """Radians to Degrees"""
    return 180.*a/n.pi

def d2ra(a):
    """Degrees to RA 24 Hours"""
    return 24.*a/360.

def r2ra(a):
    """Radians to RA 24 Hours"""
    return 24.*a/(2.*n.pi)

def deprojectMatrix(lat,dec,ha):
    """Compute the inverted d-projection matrix"""
    dproj=n.matrix([[n.cos(lat)*n.cos(dec)+n.sin(lat)*n.sin(dec)*n.cos(ha), -1.*n.sin(lat)*n.sin(ha)],[n.sin(dec)*n.sin(ha),n.cos(ha)]], dtype=n.complex)    #dipole projection matrix
    stokesDproj=n.kron(dproj,dproj) #to convert a Jones to Mueller matrix compute the auto-Kronecker product
    #invDproj=n.linalg.inv(stokesDproj) #invert [(PP outer PP*)]
    invDproj=n.linalg.pinv(stokesDproj) #invert [(PP outer PP*)] using the pseudo-inverse
    invMat=invStokes*invDproj
    return invMat

Stokes=.5*n.matrix([[1.,1.,0.,0.],[0.,0.,1.,1.j],[0.,0.,1.,-1.j],[1.,-1.,0.,0.]])   #note 1/2
invStokes=n.matrix([[1.,0.,0.,1.],[1.,0.,0.,-1.],[0.,1.,1.,0.],[0.,-1.j,1.j,0.]])

#assume dec phase direction doesn't change from zenith
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('-o','--ofn',dest='ofn',default='dproj_matrix',
        help='Output filename of Numpy array, default: dproj_matrix')
    opts, args = o.parse_args(sys.argv[1:])

    nfiles=len(args)
    fn=args[0]

    im,hdr,axisInfo=fitsUtil.readFITS(fn,hdr=True,axis=True) #load template FITS
    print axisInfo
    wcs=hdr['wcs']
    #wcs.wcs.print_contents()

    #create an output data cube
    outShape=[]
    axisDict={}
    for aid,aname in enumerate(axisInfo):
        if aname.startswith('RA'):
            outShape.append(im.shape[aid])
            axisDict['RA']=aid
        if aname.startswith('DEC'):
            outShape.append(im.shape[aid])
            axisDict['DEC']=aid
        if aname.startswith('FREQ'):
            outShape.append(1)
            axisDict['FREQ']=aid
        if aname.startswith('STOKES'):
            outShape.append(4)
            axisDict['STOKES']=aid
    
    #compute the sky coordinates for each pixel
    mms=n.repeat(n.arange(im.shape[axisDict['RA']]),im.shape[axisDict['DEC']]).reshape((im.shape[axisDict['RA']],im.shape[axisDict['DEC']]))
    lls=mms.T
    ras,decs=wcs.wcs_pix2sky(n.array(mms.flatten(), n.float_),n.array(lls.flatten(), n.float_),1)  #compute RA/DEC 
    ras=ras.reshape(im.shape[axisDict['RA']],im.shape[axisDict['DEC']])
    decs=decs.reshape(im.shape[axisDict['RA']],im.shape[axisDict['DEC']])
    cra,cdec=ras[im.shape[axisDict['RA']]/2+1,im.shape[axisDict['DEC']]/2+1],decs[im.shape[axisDict['RA']]/2+1,im.shape[axisDict['DEC']]/2+1] #centre ra,dec
    print 'Using centre position (RA: %f, DEC: %f)'%(cra,cdec)
    lst=d2r(cra)    #set central RA to LST such that the centre pixel has HA=0
    lat=d2r(cdec)   #set the array latitude to zenith declination

    #operate on each pixel, this is gonna be slow
    #the complex image and FITS image have the same shape, except the complex image has an additional dimension prefix
    print 'Computing per pixel dipole de-projection, this will take a minute...'
    dproj=n.zeros((im.shape[axisDict['RA']],im.shape[axisDict['DEC']],4,4),dtype=complex)
    #This is the loop that slows everything down, we need to speed this up via parallelization
    for ll,mm in zip(lls.flatten(),mms.flatten()):
        if mm%100==0 and ll==0: print mm
        if True in n.isnan([ras[ll,mm],decs[ll,mm]]): continue
        else:
            dec=d2r(decs[ll,mm])
            ha=lst-d2r(ras[ll,mm])
            dproj[mm,ll]=n.asarray(deprojectMatrix(lat,dec,ha))

    #save Numpy array to file
    ofn=opts.ofn+'.npy'
    n.save(open(ofn,'wb'),dproj)

