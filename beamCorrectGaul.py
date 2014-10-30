#!/usr/bin/env python

import os,sys
import numpy as np
import pyfits as pf
import pywcs
import fitsUtil

stokesDict={'I':0,'Q':1,'U':2,'V':3}
cplxDict={'XX*':0,'XY*':1,'YX*':2,'YY*':3}

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] -b BEAM_NPY -f FITS_IMAGE GAUL_FILE')
    o.set_description(__doc__)
    o.add_option('-b','--beam', dest='beam', default=None,
        help='Beam as a numpy array, assumes the beam is complex and has the same resolution and number of pixels as in the FITS image')
    o.add_option('-f','--fits', dest='fits', default=None,
        help='FITS image the sky model is derived from, required to get phase centre and resolution information')
    o.add_option('-i','--invert', dest='invert', action='store_true',
        help='Use if input complex beam has already been inverted')
    opts, args = o.parse_args(sys.argv[1:])

    #check for beam
    if opts.beam is None:
        print 'Error: Missing beam NPY'
        exit(1)
    print 'Loading',opts.beam
    beam=np.load(opts.beam)
    stokesBeam=np.zeros((1,beam.shape[1],beam.shape[2],beam.shape[3]))

    #check for fits image
    if opts.fits is None:
        print 'Error: Missing FITS image'
        exit(1)
    print 'Loading',opts.fits

    #check beam and fits of the same shape
    im,fHdr,fAxis=fitsUtil.readFITS(opts.fits,hdr=True,axis=True)
    if im.shape==stokesBeam.shape:
        stokesBeam[0]=beam[0]+beam[5] #Stokes I beam
    else:
        print 'ERROR: beam shape different than %s image shape'%opts.fits
        exit()

    #GAUL file columns
    raIdx=4
    decIdx=6
    tfIdx=8
    etfIdx=9
    pfIdx=10
    epfIdx=11

    gaulFile=args[0]
    fh=open(gaulFile,'r')
    gaulOutStr=''
    for line in fh.readlines():
        if line.startswith('# Gaus_id'): #header string
            hdrStr=line.split('\n')[0]
            hdrStr+=' beam_gain int_Total_flux, int_E_Total_flux int_Peak_flux int_E_Peak_flux\n' 
            gaulOutStr+=hdrStr
        elif line=='\n' or line.startswith('#'): gaulOutStr+=line
        else:
            wsline=" ".join(line.split())
            lineArr=wsline.split(' ')
            #get RA,DEC, Total Flux, Total FLux Error, Peak Flux, Peak FLux Error
            #print lineArr[raIdx],lineArr[decIdx],lineArr[tfIdx],lineArr[etfIdx],lineArr[pfIdx],lineArr[epfIdx]
            pxPair=fHdr['wcs'].wcs_sky2pix(np.array(lineArr[raIdx],np.float_),np.array(lineArr[decIdx],np.float_),1)
            raPix=int(pxPair[0][0])
            decPix=int(pxPair[1][0])
            gain=stokesBeam[0,0,raPix,decPix]
            outLine=line.split('\n')[0]
            if not opts.invert: gain=1./gain
            outLine+=' %f %f %f %f %f\n'%(gain,gain*float(lineArr[tfIdx]),gain*float(lineArr[etfIdx]),gain*float(lineArr[pfIdx]),gain*float(lineArr[epfIdx]))
            gaulOutStr+=outLine
    ogaulFile=gaulFile.split('/')[-1].split('.gaul')[0]+'.beamCorrect.gaul'
    ofh=open(ogaulFile,'w')
    ofh.write(gaulOutStr)
    ofh.close()
    fh.close()

