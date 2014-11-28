#!/usr/bin/env python
"""
Scale a FITS image to an absolute flux scale based on closest near by sources
"""

import sys,os,shutil
import pyfits as pf
import numpy as np
import pywcs

import fitsUtil

def computePos(srcDict):
    """Return the RA,Dec in radians"""
    raElms=map(float,srcDict['RA'].split(':'))
    srcRa=(15.*(raElms[0]+(raElms[1]/60.)+(raElms[2]/3600.)))*np.pi/180.
    decElms=map(float,srcDict['Dec'].split(':'))
    srcDec=(np.abs(decElms[0])+(decElms[1]/60.)+(decElms[2]/3600.))*np.pi/180.
    if decElms[0]<0.: srcDec*=-1.
    return srcRa,srcDec

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] STOKES_FITS')
    o.set_description(__doc__)
    o.add_option('-r','--radius',dest='radius',type='float',default=45.0,
        help='Radius, in degrees, to set maximum distance from zenith for calibration source selection default: 45')
    o.add_option('-s','--scale',dest='scale',type='float',default=None,
        help='Instead of looking nearby sources to set a flux scale, scale the data by this value')
    opts, args = o.parse_args(sys.argv[1:])

    #calibation source dictionary, positions in J2000, flux and spi taken from Culgoora, except Pic A which was from Jacobs et al. PAPER catalgoue 
    calSrcs={
            'pic_a':  {'RA':'05:19:46.3', 'Dec':'-45:49:31', 'Freq':160e6, 'S':452.0,  'SPI':-0.76}, #Pictor A / Cul 0518-458
            'tau_a':  {'RA':'05:31:30.4', 'Dec':'+21:58:55', 'Freq':160e6, 'S':1256.0, 'SPI':-0.39}, #Taurus A / Cul 0531+219
            'hydra_a':{'RA':'09:15:42.4', 'Dec':'-11:53:13', 'Freq':160e6, 'S':243.0,  'SPI':-1.12}, #Hydra A / Cul 0915-118
            '3c273':  {'RA':'12:26:32.1', 'Dec':'+02:19:14', 'Freq':160e6, 'S':102.0,  'SPI':-0.61}, #3C273 / Cul 1226+023
            'virgo_a':{'RA':'12:28:17.6', 'Dec':'+12:39:47', 'Freq':160e6, 'S':566.0,  'SPI':-1.54}, #Virgo A / Cul 1228+126
            'cent_a': {'RA':'13:23:34.3', 'Dec':'-42:44:15', 'Freq':160e6, 'S':1104.0, 'SPI':-0.52}, #Centaurus A / Cul 1322-427
            'herc_a': {'RA':'16:48:39.3', 'Dec':'+05:04:17', 'Freq':160e6, 'S':378.0,  'SPI':-1.15}, #Hercules A / Cul 1648+050
            '3c353':  {'RA':'17:17:54.1', 'Dec':'-00:55:40', 'Freq':160e6, 'S':276.0,  'SPI':-0.49}, #3C353 / Cul 1717-009
            '3c444':  {'RA':'22:11:41.7', 'Dec':'-17:16:48', 'Freq':160e6, 'S':79.2,   'SPI':-0.88}, #3C444 / Cul 2211-172
            'fornax_a':   {'RA':'03:22:36.5', 'Dec':'-37:14:21', 'Freq':408e6, 'S':259,   'SPI':-0.51}, #Fornax A / NGC 1316
            'B1610-607':  {'RA':'16:15:01.6', 'Dec':'-60:55:18', 'Freq':1430e6,'S':43.,   'SPI':-1.36}, #PKS B1610-607
            'B1932-464':  {'RA':'19:35:56.5', 'Dec':'-46:20:40', 'Freq':80e6,  'S':147.0, 'SPI':-0.8}, #PKS B1932-464
            'B2356-611':  {'RA':'23:58:59.3', 'Dec':'-60:54:56', 'Freq':408e6, 'S':57.,   'SPI':-0.76}, #PKS B2356-611
            }

    haLimit=opts.radius*np.pi/180.
    for fid,fn in enumerate(args):
        print 'Processing',fn
        im,fHdr,fAxis=fitsUtil.readFITS(fn,hdr=True,axis=True)
        #get PSF scales
        bmin=fitsUtil.getKeyword(fn,'BMIN')
        bmaj=fitsUtil.getKeyword(fn,'BMAJ')
        bpa=fitsUtil.getKeyword(fn,'BPA')
        #get image center frequency
        freq=fitsUtil.getKeyword(fn,'CRVAL3')
        zenRa=fHdr['ra']*np.pi/180.
        zenDec=fHdr['dec']*np.pi/180.

        #Gaussian PSF function
        def gaussFunc(x,y):
            rotv=np.dot(np.matix([[np.cos(bpa*np.ip/180.),-1.*np.sin(bpa*np.ip/180.)],[np.sin(bpa*np.ip/180.),np.cos(bpa*np.ip/180.)]]),np.array([x,y]))
            return 1. * np.exp(-1.*( ((rotv[0]**2.)/(2.*bmaj**2.)) + ((rotv[y]**2.)/(2.*bmin**2.)) ))

        #Find calibration sources to use
        usableSrcs=[]
        for key,val in calSrcs.iteritems():
            srcRa,srcDec=computePos(val)
            if srcRa>np.pi: HA=(zenRa-(srcRa-(2.*np.pi)))
            else: HA=(zenRa-srcRa)

            #source is usable if it is with an HA limit distance
            if np.abs(HA) < haLimit: usableSrcs.append(key)
            
        #if no sources, print error and exit
        if len(usableSrcs)<1 and opts.scale is None:
            print 'WARNING: no available calibration source, skipping'
            continue

        if opts.scale is None:
            sFactors=[]
            distances=[]
            for key in usableSrcs:
                srcRa,srcDec=computePos(calSrcs[key])
                circDist = np.arccos(np.sin(zenDec)*np.sin(srcDec)+np.cos(zenDec)*np.cos(srcDec)*np.cos(np.abs(zenRa-srcRa))) #great-circle distance
                distances.append(circDist)
                #print 'circle distance:',circDist
                pxPair=fHdr['wcs'].wcs_sky2pix(np.array(srcRa*180./np.pi),np.array(srcDec*180./np.pi),1)
                #search a region around the source center
                pixPerPSF=(bmaj+bmin)/(np.abs(fHdr['ddec'])+np.abs(fHdr['dra']))
                xmin=int(pxPair[0][0]-5*pixPerPSF/2.)
                xmax=int(pxPair[0][0]+5*pixPerPSF/2.)
                ymin=int(pxPair[1][0]-5*pixPerPSF/2.)
                ymax=int(pxPair[1][0]+5*pixPerPSF/2.)
                subIm=im[0,0,ymin:ymax,xmin:xmax]
                #TODO: compute the PSF weighted average instead of the peak
                srcFlux=np.max(subIm)
                #print srcFlux
                expectFlux=calSrcs[key]['S']*((freq/calSrcs[key]['Freq'])**calSrcs[key]['SPI'])
                #expectFlux*=(((bmaj+bmin)*pixPerPSF)*np.sqrt(2.*np.pi))#multiply by the inverse of normalization, to approximately remove the effect of the PSF
                sFactors.append(expectFlux/srcFlux)
            #sFactor=np.average(np.array(sFactors))
            #select out the source closest to Zenith
            sfid=np.argmin(np.array(distances))
            sFactor=sFactors[sfid]
        else: sFactor=opts.scale #manual scale factor
        print 'SRCS:',usableSrcs
        print 'DISTANCES:',distances
        print 'SCALES:',sFactors
        print 'SCALE:',sFactor
    
        #write outIm to new FITS file
        path,srcFn=os.path.split(os.path.realpath(fn))
        dstFn=srcFn.split('.fits')[0]+'.abs.fits'
        print 'Creating new FITS file:',os.path.join(path,dstFn)
        shutil.copyfile(os.path.join(path,srcFn), os.path.join(path,dstFn))
        fh=pf.open(os.path.join(path,dstFn),mode='update')
        fh[0].data*=sFactor
        fh.flush()
        fh.close()
            
