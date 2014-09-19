#!/usr/bin/env python
"""
Using a FEKO beam model and a template FITS file, generate an interpolated Stokes beam
"""

import sys,os,shutil
import pyfits as pf
import numpy as np
import pywcs
from scipy import interpolate

import fitsUtil
import fmt

def d2r(a):
    """Degrees to Radians"""
    return np.pi*a/180.

def r2d(a):
    """Radians to Degrees"""
    return 180.*a/np.pi

def d2ra(a):
    """Degrees to RA 24 Hours"""
    return 24.*a/360.

def r2ra(a):
    """Radians to RA 24 Hours"""
    return 24.*a/(2.*np.pi)

def genInterpStokesBeam(fekoX,fekoY,fitsfn,nstokes=-1,invert=False):
    """Generate an interpolated Stokes beam from FEKO models for an image cube
    fekoX/fekoY: feko beam model class for X and Y polarizations
    fitsfn: FITS files
    nstokes: override number of stokes parameters in the FITS header"""

    im,hdr,axisInfo=fitsUtil.readFITS(fn,hdr=True,axis=True) #load template FITS
    wcs=hdr['wcs']
    #wcs.wcs.print_contents()
    #print hdr
    #print axisInfo
    #print im.shape

    #for fid,f in enumerate(fekoX.fields):
    #    print 'Field ID %i: %f MHz'%(fid,f.freq*10e-7)

    axisDict={}
    for aid,aname in enumerate(axisInfo):
        if aname.startswith('RA'): axisDict['RA']=aid
        if aname.startswith('DEC'): axisDict['DEC']=aid
        if aname.startswith('FREQ'): axisDict['FREQ']=aid
        if aname.startswith('STOKES'): axisDict['STOKES']=aid
    #print axisDict

    freqs=np.arange(hdr['nchans'])*hdr['bw']+hdr['centreFreq']
    #print freqs
    if nstokes<0: stokes=np.arange(hdr['nstokes'])+hdr['stokes0']
    else: stokes=np.arange(nstokes,dtype='int')+int(hdr['stokes0'])

    #(l,m) grid for spatial interpolation
    phi=fekoX.fields[0].phi*np.pi/180.
    theta=fekoX.fields[0].theta*np.pi/180.
    #print phi.shape
    #Generate a WCS structure with the same resolution as the template FITS file, but centered at the North Pole
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crval = [0.,90.]
    wcs.wcs.crpix = [hdr['raPix'],hdr['decPix']]
    wcs.wcs.cdelt = [hdr['dra'],hdr['ddec']]
    wcs.wcs.ctype = ["RA---SIN", "DEC--SIN"]
    #compute the sky coordinates for each pixel
    mms=np.repeat(np.arange(im.shape[axisDict['RA']]),im.shape[axisDict['DEC']]).reshape((im.shape[axisDict['RA']],im.shape[axisDict['DEC']]))
    lls=mms.T
    imPhi,imTheta=wcs.wcs_pix2sky(np.array(mms.flatten(), np.float_),np.array(lls.flatten(), np.float_),1)  #compute phi(RA)/theta(DEC)
    imPhi=imPhi.reshape(im.shape[axisDict['RA']],im.shape[axisDict['DEC']])*np.pi/180.
    #flip the theta values
    imTheta-=90.
    imTheta*=-1.
    imTheta=imTheta.reshape(im.shape[axisDict['RA']],im.shape[axisDict['DEC']])*np.pi/180.

    interpBeamIm=np.zeros_like(im)
    for sid in stokes:
        for fid,freq in enumerate(freqs):
            print 'Interpolating Beam: Stokes %s Frequency %f MHz'%(stokesID[sid-1],freq*10e-7)
            sbStart=freq-hdr['bw']/2.
            sbStop=freq+hdr['bw']/2.
            #Select out FEKO freqs to use
            fekoIds=[]
            for fekoid,fekof in enumerate(fekoX.fields):
                if fekof.freq <= sbStop and fekof.freq >= sbStart: fekoIds.append(fekoid)
            #Frequency edge cases
            fekoStart=min(fekoIds)-1
            fekoStop=max(fekoIds)+1
            if fekoStart < 0:
                print 'WARNING: beam model does not go below frequency of image frequency, using lowest beam model frequency for images below %f MHz'%(fekoX.fields[min(fekoIds)].freq*10e-7)
            else: fekoIds.append(fekoStart)
            if fekoStop > len(fekoX.fields):
                print 'WARNING: beam model does not go above frequency of image frequency, using highest beam model frequency for images above %f MHz'%(fekoX.fields[max(fekoIds)].freq*10e-7)
            else: fekoIds.append(fekoStop)
            fekoIds=sorted(fekoIds)
            #print fekoIds

            #generate an interpolated beam for Stokes sid at frequency freq
            #compute FEKO Stokes beams
            fekoFreqs=[]
            fekoStokes=[]
            for fekoid in fekoIds:
                fekoFreqs.append(fekoX.fields[fekoid].freq)
                fieldX=fekoX.fields[fekoid]
                fieldY=fekoY.fields[fekoid]
                if stokesID[sid-1]=='I':
                    gxx=fieldX.etheta*np.conj(fieldX.etheta)+fieldX.ephi*np.conj(fieldX.ephi)
                    gyy=fieldY.etheta*np.conj(fieldY.etheta)+fieldY.ephi*np.conj(fieldY.ephi)
                    sbeam=(gxx+gyy).real
                if stokesID[sid-1]=='Q':
                    gxx=fieldX.etheta*np.conj(fieldX.etheta)+fieldX.ephi*np.conj(fieldX.ephi)
                    gyy=fieldY.etheta*np.conj(fieldY.etheta)+fieldY.ephi*np.conj(fieldY.ephi)
                    sbeam=(gxx-gyy).real
                if stokesID[sid-1]=='U':
                    gxy=fieldX.etheta*np.conj(fieldY.etheta)+fieldX.ephi*np.conj(fieldY.ephi)
                    gyx=fieldY.etheta*np.conj(fieldX.etheta)+fieldY.ephi*np.conj(fieldX.ephi)
                    sbeam=(gxy+gyx).real
                if stokesID[sid-1]=='V':
                    gxy=fieldX.etheta*np.conj(fieldY.etheta)+fieldX.ephi*np.conj(fieldY.ephi)
                    gyx=fieldY.etheta*np.conj(fieldX.etheta)+fieldY.ephi*np.conj(fieldX.ephi)
                    sbeam=(gxy-gyx).imag
                if invert: fekoStokes.append(1./sbeam)
                else: fekoStokes.append(sbeam)
            fekoStokes=np.array(fekoStokes)
            #print np.max(fekoStokes), np.min(fekoStokes)
            #print fekoStokes.shape
            #print fekoFreqs

            #freq interpolate: weighted average of freq planes, based on absolute distance from centre freq
            print 'Frequency Interpolation'
            interpSamples=10 #number of interpolated values to compute across the channel band
            interpFreqs=np.linspace(sbStart,sbStop,interpSamples,endpoint=True)
            #print interpFreqs
            freqInterpBeam=np.zeros((fekoStokes.shape[1]))
            #this for loop can probably be sped up with np.apply_along_axis()
            for bid in range(fekoStokes.shape[1]):
                freqInterpBeam[bid]=np.sum(np.interp(interpFreqs, fekoFreqs, fekoStokes[:,bid], left=fekoStokes[0,bid], right=fekoStokes[-1,bid]))/interpSamples
            #print freqInterpBeam.shape

            #spatial interpolate: triangular averaging
            print 'Spatial Interpolation'
            interpFunc=interpolate.interp2d(phi, theta, freqInterpBeam, kind='linear',bounds_error=False, fill_value=np.nan)
            #TODO:this is a very slow loop
            for ll in range(im.shape[axisDict['DEC']]):
                for mm in range(im.shape[axisDict['RA']]):
                    if np.isnan(imPhi[ll,mm]) or np.isnan(imTheta[ll,mm]): interpBeamIm[sid-1,fid,mm,ll]=np.nan
                    else: interpBeamIm[sid-1,fid,mm,ll]=interpFunc(imPhi[ll,mm],imTheta[ll,mm])

    return interpBeamIm

Stokes=.5*np.matrix([[1.,1.,0.,0.],[0.,0.,1.,1.j],[0.,0.,1.,-1.j],[1.,-1.,0.,0.]])   #note 1/2
invStokes=np.matrix([[1.,0.,0.,1.],[1.,0.,0.,-1.],[0.,1.,1.,0.],[0.,-1.j,1.j,0.]])
stokesID=['I','Q','U','V']

#assume dec phase direction doesn't change from zenith
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] --xpol=FEKO_FILE --ypol=FEKO_FILE FITS_IMAGE')
    o.set_description(__doc__)
    o.add_option('--xpol',dest='xpol',default=None,
        help='FEKO model for the X linear polarization')
    o.add_option('--ypol',dest='ypol',default=None,
        help='FEKO model for the Y linear polarization')
    o.add_option('-o','--ofn',dest='ofn',default=None,
        help='Output filename of FITS file, default: beam_[FREQ]_[RES].fits')
    o.add_option('-i','--invert',dest='invert',action='store_true',
        help='Interpolate the inverted beam, use if output is being used to correct image')
    o.add_option('-c','--clip',dest='clip',default=None,
        help='Clip the beam below this level, before interpolation, if using the invert option, this is the maximum value to clip above')
    opts, args = o.parse_args(sys.argv[1:])

    fn=args[0]

    #Read in FEKO format beam models
    fekoX=fmt.FEKO(opts.xpol)
    fekoY=fmt.FEKO(opts.ypol)

    beam=genInterpStokesBeam(fekoX,fekoY,fn,nstokes=4,invert=opts.invert)
    path,srcFn=os.path.split(os.path.realpath(fn))
    if opts.ofn is None:
        centreFreq=fitsUtil.getKeyword(fn,'CRVAL3')
        res=fitsUtil.getKeyword(fn,'CDELT1')
        dstFn='beam_f%.2f_r%i.fits'%(centreFreq*10e-7,int(np.abs(res*3600.)))
    else: dstFn=opts.ofn
    print 'Creating new FITS file:',os.path.join(path,dstFn)
    shutil.copyfile(os.path.join(path,srcFn), os.path.join(path,dstFn))
    fh=pf.open(os.path.join(path,dstFn),mode='update')
    fh[0].data=beam
    #Change RA/DEC header items to point at North Pole
    hdr=fh[0].header
    hdr.update('CRVAL1',0.) #RA
    hdr.update('CRVAL2',90.) #DEC
    fh.flush()
    fh.close()

