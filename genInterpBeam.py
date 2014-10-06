#!/usr/bin/env python
"""
Using a FEKO beam model and a template FITS file, generate an interpolated Stokes beam (and possibly invert) and write to a numpy file
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

def genInterpBeam(fekoX,fekoY,fitsfn,invert=False):
    """Generate an interpolated complex beam from FEKO models for an image cube
    fekoX/fekoY: feko beam model class for X and Y polarizations
    fitsfn: FITS files
    """

    im,hdr,axisInfo=fitsUtil.readFITS(fitsfn,hdr=True,axis=True) #load template FITS
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

    #(l,m) grid for spatial interpolation
    phi=fekoX.fields[0].phi*np.pi/180.
    theta=fekoX.fields[0].theta*np.pi/180.
    #print phi.shape
    #print np.max(theta), np.max(phi)
    
    #drop beam points near zenith as this can create artefacts
    deltaTheta=(np.pi/2.)/100. #delta of theta away from pi/2 to ignore
    thetaIdx=np.argwhere(theta<(np.pi/2.-deltaTheta))
    theta0=theta[thetaIdx].flatten()
    phi0=phi[thetaIdx].flatten()

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

    interpBeamIm=np.zeros((len(polID),im.shape[axisDict['FREQ']],im.shape[axisDict['RA']],im.shape[axisDict['DEC']]))
    #print interpBeamIm.shape
    for fid,freq in enumerate(freqs):
        freqInterpBeams={}
        for pid,plabel in enumerate(polID):
            print 'Interpolating Beam: %s Frequency %f MHz'%(plabel,freq*10e-7)
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

            #generate an interpolated beam for pid at frequency freq
            #compute FEKO complex beams
            fekoFreqs=[]
            fekoCmplxBeams=[]
            for fekoid in fekoIds:
                fekoFreqs.append(fekoX.fields[fekoid].freq)
                fieldX=fekoX.fields[fekoid]
                fieldY=fekoY.fields[fekoid]
                if plabel=='XX':
                    gxx=fieldX.etheta*np.conj(fieldX.etheta)+fieldX.ephi*np.conj(fieldX.ephi)
                    fekoCmplxBeams.append(gxx.real)
                elif plabel=='YY':
                    gyy=fieldY.etheta*np.conj(fieldY.etheta)+fieldY.ephi*np.conj(fieldY.ephi)
                    fekoCmplxBeams.append(gyy.real)
                elif plabel=='XYr':
                    gxy=fieldX.etheta*np.conj(fieldY.etheta)+fieldX.ephi*np.conj(fieldY.ephi)
                    fekoCmplxBeams.append(gxy.real)
                elif plabel=='XYi':
                    gxy=fieldX.etheta*np.conj(fieldY.etheta)+fieldX.ephi*np.conj(fieldY.ephi)
                    fekoCmplxBeams.append(gxy.imag)
                elif plabel=='YXr':
                    gyx=fieldY.etheta*np.conj(fieldX.etheta)+fieldY.ephi*np.conj(fieldX.ephi)
                    fekoCmplxBeams.append(gyx.real)
                elif plabel=='YXi':
                    gyx=fieldY.etheta*np.conj(fieldX.etheta)+fieldY.ephi*np.conj(fieldX.ephi)
                    fekoCmplxBeams.append(gyx.imag)
            fekoCmplxBeams=np.array(fekoCmplxBeams)
            #print np.max(fekoCmplxBeams), np.min(fekoCmplxBeams)
            #print fekoCmplxBeams.shape
            #print fekoFreqs

            #freq interpolate: weighted average of freq planes, based on absolute distance from centre freq
            print 'Frequency Interpolation'
            interpSamples=10 #number of interpolated values to compute across the channel band
            interpFreqs=np.linspace(sbStart,sbStop,interpSamples,endpoint=True)
            #print interpFreqs
            freqInterpBeam=np.zeros((fekoCmplxBeams.shape[1]))
            #this for loop can probably be sped up with np.apply_along_axis()
            for bid in range(fekoCmplxBeams.shape[1]):
                freqInterpBeam[bid]=np.sum(np.interp(interpFreqs, fekoFreqs, fekoCmplxBeams[:,bid], left=fekoCmplxBeams[0,bid], right=fekoCmplxBeams[-1,bid]))/interpSamples
            #print freqInterpBeam.shape

            freqInterpBeam0=freqInterpBeam[thetaIdx].flatten()
            #print theta.shape,theta0.shape
            freqInterpBeams[plabel]=freqInterpBeam0

        for pid,plabel in enumerate(polID):
            print 'Spatial Interpolation',plabel
            #interpFunc=interpolate.interp2d(phi, theta, freqInterpBeam, kind='linear',bounds_error=False, fill_value=np.nan)
            #interpFunc=interpolate.interp2d(phi0, theta0, freqInterpBeam0, kind='linear',bounds_error=False, fill_value=np.nan)
            interpFunc=interpolate.interp2d(phi0, theta0, freqInterpBeams[plabel], kind='linear',bounds_error=False, fill_value=np.nan)
            #TODO:this is a very slow loop
            for ll in range(im.shape[axisDict['DEC']]):
                for mm in range(im.shape[axisDict['RA']]):
                    if np.isnan(imPhi[ll,mm]) or np.isnan(imTheta[ll,mm]): interpBeamIm[pid,fid,mm,ll]=np.nan
                    else: interpBeamIm[pid,fid,mm,ll]=interpFunc(imPhi[ll,mm],imTheta[ll,mm])

        if invert:
            print 'Inverting frequency-interpolated beam'
            cplxBeamIm=np.array([ [interpBeamIm[0,fid],interpBeamIm[1,fid]+1j*interpBeamIm[2,fid]], [interpBeamIm[3,fid]+1j*interpBeamIm[4,fid],interpBeamIm[5,fid]]])
            invCplxBeam=np.zeros_like(cplxBeamIm)
            #slow loop
            for ll in range(im.shape[axisDict['DEC']]):
                for mm in range(im.shape[axisDict['RA']]):
                    if np.isnan(cplxBeamIm[:,:,ll,mm]).any(): invCplxBeam[:,:,ll,mm]=np.nan
                    else: invCplxBeam[:,:,ll,mm]=np.linalg.pinv(np.matrix(cplxBeamIm[:,:,ll,mm])) #pseudo-inverse
                    #else: invCplxBeam[:,:,ll,mm]=np.linalg.inv(np.matrix(cplxBeamIm[:,:,ll,mm]))
            interpBeamIm[0,fid]=invCplxBeam[0,0].real
            interpBeamIm[1,fid]=invCplxBeam[0,1].real
            interpBeamIm[2,fid]=invCplxBeam[0,1].imag
            interpBeamIm[3,fid]=invCplxBeam[1,0].real
            interpBeamIm[4,fid]=invCplxBeam[1,0].imag
            interpBeamIm[5,fid]=invCplxBeam[1,1].real

    return interpBeamIm

Stokes=.5*np.matrix([[1.,1.,0.,0.],[0.,0.,1.,1.j],[0.,0.,1.,-1.j],[1.,-1.,0.,0.]])   #note 1/2
invStokes=np.matrix([[1.,0.,0.,1.],[1.,0.,0.,-1.],[0.,1.,1.,0.],[0.,-1.j,1.j,0.]])
stokesID=['I','Q','U','V']
polID=['XX','XYr','XYi','YXr','YXi','YY']


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
        help='Output filename of numpy file, default: beam_[FREQ]_[RES].npy')
    o.add_option('-i','--invert',dest='invert',action='store_true',
        help='Interpolate the inverted beam, use if output is being used to correct image')
    o.add_option('-S','--stokes',dest='stokes',action='store_true',
        help='Output the Stokes beams instead of the complex gain beams')
    opts, args = o.parse_args(sys.argv[1:])

    fn=args[0]

    #Read in FEKO format beam models
    fekoX=fmt.FEKO(opts.xpol)
    fekoY=fmt.FEKO(opts.ypol)

    beam=genInterpBeam(fekoX,fekoY,fn,invert=opts.invert)

    #convert output to Stokes if option set
    stokesStr=''
    if opts.stokes:
        print 'Generating Stokes beams'
        stokesBeam=np.zeros((4,beam.shape[1],beam.shape[2],beam.shape[3]))
        cplxBeam=np.zeros((4,beam.shape[1],beam.shape[2],beam.shape[3]),dtype='complex')
        cplxBeam[0]=beam[0]
        cplxBeam[1]=beam[1]+1j*beam[2]
        cplxBeam[2]=beam[3]+1j*beam[4]
        cplxBeam[3]=beam[5]
        stokesBeam[0]=(cplxBeam[0]+cplxBeam[3]).real
        stokesBeam[1]=(cplxBeam[0]-cplxBeam[3]).real
        stokesBeam[2]=(cplxBeam[1]+cplxBeam[2]).real
        stokesBeam[3]=(cplxBeam[1]-cplxBeam[2]).imag
        stokesStr='.stokes'
        beam=stokesBeam

    path,srcFn=os.path.split(os.path.realpath(fn))
    if opts.ofn is None:
        centreFreq=fitsUtil.getKeyword(fn,'CRVAL3')
        res=fitsUtil.getKeyword(fn,'CDELT1')
        if opts.invert: invertStr='.invert'
        else: invertStr=''
        dstFn='beam_f%.2f_r%i%s%s.npy'%(centreFreq*10e-7,int(np.abs(res*3600.)),invertStr,stokesStr)
    else: dstFn=opts.ofn
    print 'Creating numpy file:',os.path.join(path,dstFn)
    np.save(os.path.join(path,dstFn),beam)

