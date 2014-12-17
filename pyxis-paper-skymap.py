# pyxis-paper-skymap.py: a Pyxis recipe to generate PAPER PSA32 Skymaps

# All recipes start with this line, which initializes Pyxis
import Pyxis

# import some Pyxides modules.
# Note that Pyxides is implicitly added to the include path
# when you import Pyxis above
import ms,imager,std,lsm,stefcal

# import some other Python modules that we make use of below
import pyfits
import pyrap.tables
import os.path
import glob
import shutil

#SCRIPTS_DIR='/home/foster/PAPER/pipelines/paper32imaging/scripts'
SCRIPTS_DIR='/home/foster/sw/paperScripts'

# ELWOOD specific settings
#mqt.MULTITHREAD = 4
#local imagers
imager.LWIMAGER_PATH = "/home/foster/meqtrees/git/lwimager/build/synthesis/lwimager"
WSCLEAN_PATH="/home/foster/sw/wsclean-code/wsclean/build/wsclean"

#Hack to force decimal numbers to use '.' instead of ','
x.sh('export LC_NUMERIC=en_GB.utf8')

## default MS -- use this if nothing specified on command line
## commented out in this case, as we always want to specify stuff on the command line,
## or else go with the default list below
#v.MS = 'zen.2455819.55853.uvcRREM.MS'
## default list of measurement sets for e.g. runall()
#v.MS_List = glob.glob('/home/foster/PAPER/pipelines/paper32imaging/data/zen.*.uvcRREM.MS')
v.MS_List = glob.glob('/home/foster/PAPER/pipelines/paper32imaging/data2/zen.2455819.*.uvcRREM.MS')
#v.MS_List = glob.glob('/home/foster/PAPER/pipelines/paper32imaging/data2/zen.2455819.4[0-5]*.uvcRREM.MS')

ms.DDID = 0
ms.FIELD = 0

## destination directory for plots, images, etc.
v.DESTDIR_Template = '${OUTDIR>/}phase-nobeam-sun-output-plots${-stage<STAGE}'
#v.DESTDIR_Template = '${OUTDIR>/}phase-nobeam-natural-sun-output-plots${-stage<STAGE}'
#v.DESTDIR_Template = '${OUTDIR>/}phase-sun-output-plots${-stage<STAGE}'
#v.DESTDIR_Template = '${OUTDIR>/}phase-sun-natural-output-plots${-stage<STAGE}'
## base filename for these files
v.OUTFILE_Template = '${DESTDIR>/}${MS:BASE}${_s<STEP}${_<LABEL}'

lsm.PYBDSM_POLARIZED = False 

## when running stefcal be verbose
stefcal.STEFCAL_TDLOPTS = "stefcal_verbose=3"

## 0 means we make a single image from all channels. Use 1 to
## make a per-channel cube.
#imager.IMAGE_CHANNELIZE = 0 
### lwimager options for this field 
##imager.npix	    = 2048
##imager.cellsize	= "202arcsec"
##imager.mode	    = "channel"
##imager.stokes	  = "IQUV"
##imager.weight	  = "natural"
##imager.wprojplanes = 128
#
##imager.npix	    = 512
##imager.cellsize	= "400arcsec"
##imager.mode	    = "channel"
##imager.stokes	  = "IQUV"
##imager.weight	  = "natural"
##imager.wprojplanes = 64
#
##imager.npix	    = 1024
##imager.cellsize	= "280arcsec" #4.67 arcmin res
##imager.mode	    = "channel"
##imager.stokes	  = "IQUV"
##imager.weight	  = "natural"
##imager.wprojplanes = 64
#
#imager.npix	    = 1024
#imager.cellsize	= "330arcsec" #5.46 arcmin res
#imager.mode	    = "channel"
#imager.stokes	  = "IQUV"
#imager.weight	  = "natural"
#imager.wprojplanes = 128
#
## cleaning options
#imager.niter	    = 300
#imager.gain	      = .1
#imager.threshold = 0

## default logfile is based on MS name
v.LOG_Template = "${OUTDIR>/}log-${MS?pyxis:BASE}.log"

#Channel/Band selection
chanTotal=201 #total number of channels in MS
bw=100. #total bandwidth in MHz
sdf=bw/chanTotal #channel resolution
startChan=20 #First channel to use in pipeline, inclusive
endChan=181 #last channel to use in pipeline, exclusive
subband=10. #bandwidth (in MHz) of the final image slices
nchan=endChan-startChan #number of channels to use in pipeline
clippedBand=nchan*sdf #bandwidth to use in pipeline
nsubbands=int(clippedBand/subband) #number of subbands to image
#CHANRANGE:channel range, as first,last[,step], or list of such tuples per DDID, or None for all
ms.CHANRANGE=startChan,endChan
# Default set of interferometers to use, see MeqTree documentation for notation
#ms.IFRS = "S83,-5*"

###################################################################################################
#swift dirty images
def swift_image_pipeline(goto_step=0.):

    if not os.path.isdir(v.DESTDIR):
        os.makedirs(v.DESTDIR)

    ##copy data to corrected data column to run clean
    #if goto_step <= 0:
    #    v.STEP=1
    #    info("########## copying data to corrected data column step %i"%v.STEP)
    #    pper("MS",copy_data_to_corrected_data)

    #run a basic clean
    if goto_step <= 1:
        v.STEP=1
        info("########## making initial clean image step %i"%v.STEP)
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 250,
            'threshold': 1.,
            'pol': 'I',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': 1}
        for msFile in sorted(v.MS_List):
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)


###################################################################################################
#selfcal/imaging/source-finding pipeline function
def selfcal_loop_pipeline(goto_step=0.):

    if not os.path.isdir(v.DESTDIR):
        os.makedirs(v.DESTDIR)

    #copy data to corrected data column to run clean
    if goto_step <= 0.:
        v.STEP=1
        info("########## adding columns %i"%v.STEP)
        for msFile in sorted(v.MS_List):pyrap.tables.addImagingColumns(msFile)

    if goto_step <= 0.5:
        v.STEP=1
        info("########## copying data to corrected data column step %i"%v.STEP)
        pper("MS",copy_data_to_corrected_data)

    #run a basic clean to perform source finding on
    if goto_step <= 1:
        v.STEP=1
        info("########## making initial clean image step %i"%v.STEP)
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 250,
            'threshold': 1.,
            'pol': 'I',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': 1}
        for msFile in sorted(v.MS_List):
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #derive an initial sky model from the clean'd image
    if goto_step <= 1.5:
        v.STEP=1
        fitsfiles=glob.glob(v.DESTDIR+'/*_s1-image.fits')
        for ff in sorted(fitsfiles):
            olsm=ff.split('-image.fits')[0]+'.lsm'
            lsm.pybdsm_search(image=ff,output=olsm,thresh_pix=20,thresh_isl=15)
            lsm.tigger_convert(olsm)

    #run stefcal using the initial sky model
    if goto_step <= 2:
        v.STEP=2
        lsm0files=glob.glob(v.DESTDIR+'/*_s1.lsm.html')
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 500,
            'threshold': 1.,
            'pol': 'I',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': 1}
        for msFile in sorted(v.MS_List):
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            lsm0=filter(lambda x: x.split('/')[-1].startswith(msBase),lsm0files)[0]
            v.LSM=lsm0
            v.MS=msFile
            stefcal.STEFCAL_STEP_INCR = 0 # we set the step counter explicitly
            info("########## initial calibration of %s using sky model %s step %i"%(msFile,lsm0,v.STEP))
            stefcal.stefcal(output='CORR_DATA',dirty=False,restore=False,section='calico-stefcal')
            info("########## making clean image step %i"%v.STEP)
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #derive an second sky model from the clean'd image
    if goto_step <= 2.5:
        v.STEP=2
        fitsfiles=glob.glob(v.DESTDIR+'/*_s2-image.fits')
        for ff in sorted(fitsfiles):
            olsm=ff.split('-image.fits')[0]+'.lsm'
            lsm.pybdsm_search(image=ff,output=olsm,thresh_pix=15,thresh_isl=10)
            lsm.tigger_convert(olsm)

    #run stefcal using the second sky model
    if goto_step <= 3:
        v.STEP=3
        lsm1files=glob.glob(v.DESTDIR+'/*_s2.lsm.html')
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 1000,
            'threshold': 1.,
            'pol': 'I',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': 1}
        for msFile in sorted(v.MS_List):
            x.sh('export LC_NUMERIC=en_GB.utf8')
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            lsm1=filter(lambda x: x.split('/')[-1].startswith(msBase),lsm1files)[0]
            v.LSM=lsm1
            v.MS=msFile
            stefcal.STEFCAL_STEP_INCR = 0 # we set the step counter explicitly
            info("########## second calibration of %s using sky model %s step %i"%(msFile,lsm1,v.STEP))
            stefcal.stefcal(output='CORR_DATA',dirty=False,restore=False,section='calico-stefcal')
            info("########## making clean image step %i"%v.STEP)
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #derive an sky model from the clean'd image
    if goto_step <= 3.5:
        v.STEP=3
        fitsfiles=glob.glob(v.DESTDIR+'/*_s3-image.fits')
        for ff in sorted(fitsfiles):
            olsm=ff.split('-image.fits')[0]+'.lsm'
            lsm.pybdsm_search(image=ff,output=olsm,thresh_pix=12,thresh_isl=5)
            lsm.tigger_convert(olsm)

    #run stefcal using the third sky model and produce final complex, multi-frequency images
    if goto_step <= 4:
        v.STEP=4
        lsm2files=glob.glob(v.DESTDIR+'/*_s3.lsm.html')
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 1000,
            'threshold': 1.,
            'pol': 'xx,xy,yx,yy',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': nsubbands,
            'joinpolarizations': None}
        for msFile in sorted(v.MS_List):
            x.sh('export LC_NUMERIC=en_GB.utf8')
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            lsm2=filter(lambda x: x.split('/')[-1].startswith(msBase),lsm2files)[0]
            v.LSM=lsm2
            v.MS=msFile
            stefcal.STEFCAL_STEP_INCR = 0 # we set the step counter explicitly
            info("########## second calibration of %s using sky model %s step %i"%(msFile,lsm2,v.STEP))
            stefcal.stefcal(output='CORR_DATA',dirty=False,restore=False,section='calico-stefcal')
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #clean up
    if goto_step <= 5:
        gaulfiles=glob.glob('*.gaul')
        for gf in gaulfiles: shutil.move(gf,v.DESTDIR)

###################################################################################################
# Image domain corrections pipeline function
def image_domain_corrections_pipeline(goto_step=0):

    useStep=3   #use this step ID for selecting FITS files
    #generate the dproj matrix based on an image template
    if goto_step <= 1.:
        v.STEP=3
        fitsfiles=glob.glob(v.DESTDIR+'/*s%i-MFS-XX-image.fits'%(useStep))
        hdulist=pyfits.open(fitsfiles[0])
        hdr=hdulist[0].header
        dprojStr='paper32_%iarcsec_%ipx.dproj'%(int(abs(hdr['CDELT1']*3600.)),hdr['NAXIS1'])
        info('########## Using %s.npy for d-projection correction'%dprojStr)
        hdulist.close()
        if not os.path.isfile(dprojStr+'.npy'):
            info("########## generating D-Projection matrix from template, outputing %s step %i"%(dprojStr,v.STEP))
            generate_dproj_matrix(fitsfiles[0],dprojFn=dprojStr)

    #apply dproj to compleximage
    if goto_step <= 1.5:
        v.STEP=3
        info("########## applying D-Projection correction step %i"%v.STEP)
        #fitsfiles=glob.glob(v.DESTDIR+'/*s%i.restored.fits'%(useStep))
        fitsXX=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XX-image.fits'%(useStep)))
        fitsYY=sorted(glob.glob(v.DESTDIR+'/*s%i-*-YY-image.fits'%(useStep)))
        fitsXYr=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XY-image.fits'%(useStep)))
        fitsXYi=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XYi-image.fits'%(useStep)))
        for fxx,fyy,fxyr,fxyi in zip(fitsXX,fitsYY,fitsXYr,fitsXYi):
            ofn=fxx.split('/')[-1]
            ofn=ofn.split('XX')[0]+'stokes.fits'
            info("################ D-Projection: %s"%(ofn))
            x.sh('%s/applyDproj.py -d %s --xx=%s --yy=%s --xyr=%s --xyi=%s -o %s'%(SCRIPTS_DIR,dprojStr+'.npy',fxx,fyy,fxyr,fxyi,ofn))

    #generate beam model based on image template
    #FEKO beam model location
    fekoXpol='/home/foster/PAPER/pipelines/paper32imaging/data/PAPER_FF_X.ffe'
    fekoYpol='/home/foster/PAPER/pipelines/paper32imaging/data/PAPER_FF_Y.ffe'
    if goto_step <= 2.:
        v.STEP=4
        globTemplate=['/*-MFS-stokes.fits']
        for sbid in range(nsubbands): globTemplate.append('/*-%04i-stokes.fits'%sbid)
        for gtStr in globTemplate:
            fitsfiles=sorted(glob.glob(v.DESTDIR+gtStr))
            hdulist=pyfits.open(fitsfiles[0])
            hdr=hdulist[0].header
            beamStr='paper32_%.2fMHz_bw%04.1fMHz_%iarcsec_%ipx.beam.npy'%(hdr['CRVAL3']*10e-7,hdr['CDELT3']*10e-7,int(abs(hdr['CDELT1']*3600.)),int(hdr['NAXIS1']))
            invBeamStr='paper32_%.2fMHz_bw%04.1fMHz_%iarcsec_%ipx.invert.beam.npy'%(hdr['CRVAL3']*10e-7,hdr['CDELT3']*10e-7,int(abs(hdr['CDELT1']*3600.)),int(hdr['NAXIS1']))
            hdulist.close()
            info('########## Using %s for beam template'%invBeamStr)
            if not os.path.isfile(v.DESTDIR+'/'+invBeamStr):
                info("########## generating interpolated beam from template FITS, outputing %s step %i"%(v.DESTDIR+'/'+beamStr,v.STEP))
                x.sh('%s/genInterpBeam.py -i --xpol=%s --ypol=%s --ofn=%s %s'%(SCRIPTS_DIR,fekoXpol,fekoYpol,invBeamStr,fitsfiles[0]))
            info('########## Using %s for beam template'%beamStr)
            if not os.path.isfile(v.DESTDIR+'/'+beamStr):
                info("########## generating interpolated beam from template FITS, outputing %s step %i"%(v.DESTDIR+'/'+beamStr,v.STEP))
                x.sh('%s/genInterpBeam.py --xpol=%s --ypol=%s --ofn=%s %s'%(SCRIPTS_DIR,fekoXpol,fekoYpol,beamStr,fitsfiles[0]))

    #apply inverted beam to Stokes image
    if goto_step <= 2.5:
        v.STEP=4
        info("########## applying beam correction step %i"%v.STEP)
        fitsfiles=sorted(glob.glob(v.DESTDIR+'/*-stokes.fits'))
        for ff in fitsfiles:
            hdulist=pyfits.open(ff)
            hdr=hdulist[0].header
            beamStr='paper32_%.2fMHz_bw%04.1fMHz_%iarcsec_%ipx.invert.beam.npy'%(hdr['CRVAL3']*10e-7,hdr['CDELT3']*10e-7,int(abs(hdr['CDELT1']*3600.)),int(hdr['NAXIS1']))
            hdulist.close()
            x.sh('%s/applyInterpBeam.py -b %s %s'%(SCRIPTS_DIR,v.DESTDIR+'/'+beamStr,ff))

    
    #absolute calibration
    if goto_step <= 3.:
        v.STEP=5
        info("########## Correcting for Absolute calibration %i"%v.STEP)
        fitsfiles=sorted(glob.glob(v.DESTDIR+'/*-stokes.beam.fits'))
        for ff in fitsfiles:
            x.sh('%s/setAbsFlux.py %s'%(SCRIPTS_DIR,ff))

    #make fits freqeuncy image cube
    if goto_step <= 4.:
        v.STEP=6
        info("########## generating multi-frequency FITS files step %i"%v.STEP)
        mfitsfiles=sorted(glob.glob(v.DESTDIR+'/*0000-stokes.beam.abs.fits'))
        for mff in mfitsfiles:
            mfreqList=mff.split('-0000-')
            mfreqFn=mfreqList[0]+'.mfreq.'+mfreqList[1]
            x.sh('cp %s %s'%(mff,mfreqFn)) #copy MFS file to .mfreq.fits
            hdulist=pyfits.open(mfreqFn,mode='update') #open .mfreq.fits
            im=hdulist[0].data
            mfreqIm=np.zeros((im.shape[0],nsubbands,im.shape[2],im.shape[3]))
            for sbid in range(nsubbands):
                ff=mff.split('-0000-')[0]+'-%04i-'%sbid+'stokes.beam.abs.fits'
                subhdulist=pyfits.open(ff)
                imdata=subhdulist[0].data
                mfreqIm[:,sbid]=subhdulist[0].data[:,0]
                subhdulist.close()
            hdr=hdulist[0].header
            hdulist[0].data=mfreqIm
            hdulist.flush()
            hdulist.close()

    #TODO: cleanup intermediate files

###################################################################################################
# HEALPIX skymap pipeline function
def map_pipeline(goto_step=0):
    #postFixStr='-stokes'
    #postFixStr='-stokes.beam'
    postFixStr='-stokes.abs'
    #postFixStr='-stokes.beam.abs'

    #convert FITS to healpix
    if goto_step <= 1.:
        v.STEP=7
        info("########## generating healpix maps step %i"%v.STEP)
        fitsfiles=sorted(glob.glob(v.DESTDIR+'/*%s.fits'%postFixStr))
        for ff in fitsfiles:
            info("################ %s"%ff)
            for sid in ['I','Q','U','V']:
                hpxfn=ff.split('.fits')[0]+'.s%s.hpx'%sid
                #nside: 512 -> 7 arcmin res
                #nside: 1024 -> 3.5 arcmin res
                #x.sh('%s/mk_map_mod.py %s -n --min_alt=30 -S %s --nside=512 -m %s'%(SCRIPTS_DIR,ff,sid,hpxfn))
                x.sh('%s/mk_map_mod.py %s -n --min_alt=15 -S %s --nside=512 -m %s'%(SCRIPTS_DIR,ff,sid,hpxfn))

    #summed healpix map
    if goto_step <= 2.:
        v.STEP=8
        for sid in ['I','Q','U','V']:
            outHpxFn=v.DESTDIR+'/skymap_paper32_cMFS_s%s.hpx'%(sid)
            x.sh('%s/save_map.py -s %s %s'%(SCRIPTS_DIR,outHpxFn,v.DESTDIR+'/*MFS%s.s%s.hpx'%(postFixStr,sid)))
            for sbid in range(nsubbands):
                outHpxFn=v.DESTDIR+'/skymap_paper32_c%04i_s%s.hpx'%(sbid,sid)
                x.sh('%s/save_map.py -s %s %s'%(SCRIPTS_DIR,outHpxFn,v.DESTDIR+'/*-%04i%s.s%s.hpx'%(sbid,postFixStr,sid)))

    #TODO: cleanup intermediate files

###################################################################################################
#Generate images from CORRECTED DATA
def gen_images_from_corrected_data(goto_step=1.):
    useStep=1
    #generate multi-frequency, complex images
    if goto_step <= 1.:
        v.STEP=1
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 10,
            'threshold': 1.,
            'pol': 'xx,xy,yx,yy',
            'weight': 'natural',
            #'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': nsubbands,
            'joinpolarizations': None}
        for msFile in sorted(v.MS_List):
            x.sh('export LC_NUMERIC=en_GB.utf8')
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #convert complex images to Stokes
    if goto_step <= 2.:
        v.STEP=2
        info("########## convertinf to Stokes step %i"%v.STEP)
        fitsXX=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XX-image.fits'%(useStep)))
        fitsYY=sorted(glob.glob(v.DESTDIR+'/*s%i-*-YY-image.fits'%(useStep)))
        fitsXYr=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XY-image.fits'%(useStep)))
        fitsXYi=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XYi-image.fits'%(useStep)))
        for fxx,fyy,fxyr,fxyi in zip(fitsXX,fitsYY,fitsXYr,fitsXYi):
            ofn=fxx.split('/')[-1]
            ofn=ofn.split('XX')[0]+'stokes.fits'
            #x.sh('%s/applyDproj.py -d %s --xx=%s --yy=%s --xyr=%s --xyi=%s -o %s'%(SCRIPTS_DIR,dprojStr+'.npy',fxx,fyy,fxyr,fxyi,ofn))
            x.sh('/home/foster/sw/paperScripts/convertStokes.py --xx=%s --yy=%s --xyr=%s --xyi=%s -o %s'%(fxx,fyy,fxyr,fxyi,ofn))

###################################################################################################
#Self-cal from phase-only corrected gain, and sun removal
def phase_precal_selfcal_pipeline(goto_step=3.):
    useStep=1
    skipCalMSlist=[ 'zen.2455819.25927.uvcRREM',
                    'zen.2455819.27319.uvcRREM',
                    'zen.2455819.32190.uvcRREM',
                    'zen.2455819.34278.uvcRREM',
                    'zen.2455819.34974.uvcRREM',
                    'zen.2455819.38454.uvcRREM',
                    'zen.2455819.40542.uvcRREM',
                    'zen.2455819.41238.uvcRREM',
                    'zen.2455819.41934.uvcRREM',
                    'zen.2455819.42630.uvcRREM',
                    'zen.2455819.73947.uvcRREM',
                    'zen.2455819.74643.uvcRREM',
                    'zen.2455819.75339.uvcRREM',
                    'zen.2455819.76035.uvcRREM',
                    'zen.2455819.78819.uvcRREM',
                    'zen.2455819.79515.uvcRREM',
                    'zen.2455819.80211.uvcRREM',
                    'zen.2455819.80906.uvcRREM'] #MS files to skip calibration on
    skipCleanMSlist=['zen.2455819.79515.uvcRREM'] #problematic MS files when imaging multi-freq

    if not os.path.isdir(v.DESTDIR):
        os.makedirs(v.DESTDIR)

    if goto_step <= 0.:
        v.STEP=1
        info("########## adding columns %i"%v.STEP)
        for msFile in sorted(v.MS_List):pyrap.tables.addImagingColumns(msFile)

    #copy data to corrected data column to run clean
    if goto_step <= 0.5:
        v.STEP=1
        info("########## copying data to corrected data column step %i"%v.STEP)
        pper("MS",copy_data_to_corrected_data)

    #run a basic clean to perform source finding on
    if goto_step <= 1:
        v.STEP=1
        info("########## making initial clean image step %i"%v.STEP)
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 250,
            'threshold': 100.,
            'pol': 'I',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': 1}
        for msFile in sorted(v.MS_List):
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #derive a shallow depth sky model from the clean'd image
    if goto_step <= 1.5:
        v.STEP=1
        fitsfiles=glob.glob(v.DESTDIR+'/*_s1-image.fits')
        for ff in sorted(fitsfiles):
            olsm=ff.split('-image.fits')[0]+'.lsm'
            lsm.pybdsm_search(image=ff,output=olsm,thresh_pix=20,thresh_isl=15)
            lsm.tigger_convert(olsm)

    #run stefcal using the sky model, produce a Stokes I image for source finding
    if goto_step <= 2:
        v.STEP=2
        lsm0files=glob.glob(v.DESTDIR+'/*_s1.lsm.html')
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            'niter': 250,
            'threshold': 100.,
            'pol': 'I',
            #'weight': 'natural',
            'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': 1}
        for msFile in sorted(v.MS_List):
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            if msBase in skipCalMSlist:
                info('Skipping this MS file')
            else:
                x.sh('export LC_NUMERIC=en_GB.utf8')
                lsm0=filter(lambda x: x.split('/')[-1].startswith(msBase),lsm0files)[0]
                v.LSM=lsm0
                v.MS=msFile
                stefcal.STEFCAL_STEP_INCR = 0 # we set the step counter explicitly
                info("########## initial calibration of %s using sky model %s step %i"%(msFile,lsm0,v.STEP))
                stefcal.stefcal(output='CORR_DATA',dirty=False,restore=False,section='calico-stefcal')
            info("########## making clean image step %i"%v.STEP)
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #derive a deeper depth sky model from the clean'd image
    if goto_step <= 2.5:
        v.STEP=2
        fitsfiles=glob.glob(v.DESTDIR+'/*_s2-image.fits')
        for ff in sorted(fitsfiles):
            olsm=ff.split('-image.fits')[0]+'.lsm'
            lsm.pybdsm_search(image=ff,output=olsm,thresh_pix=10,thresh_isl=7)
            lsm.tigger_convert(olsm)

    #generate multi-frequency, complex images for image domain corrections
    if goto_step <= 3:
        v.STEP=3
        wscleanDict={
            'name':'wsclean',
            'datacolumn':'CORRECTED_DATA',
            'size': '1024 1024',
            #'scale': 330./3600.,
            'scale': 0.11,
            #'niter': 1000,
            'niter': 1,
            'threshold': 100.,
            'pol': 'xx,xy,yx,yy',
            'weight': 'natural',
            #'weight': 'briggs 0.0',
            'nwlayers':32,
            'channelrange': '%i %i'%(startChan,endChan),
            'channelsout': nsubbands,
            'joinpolarizations': None}
        for msFile in sorted(v.MS_List):
            x.sh('export LC_NUMERIC=en_GB.utf8')
            msBase=msFile.split('/')[-1].split('.MS')[0] 
            if msBase in skipCleanMSlist: wscleanDict['niter']=1
            else: wscleanDict['niter']=1
            #else: wscleanDict['niter']=1000
            info("########## making clean image step %i"%v.STEP)
            wscleanDict['name']=msBase+'_s%i'%v.STEP
            run_wsclean(msFile,wscleanDict)
        fitsfiles=glob.glob('*.fits')
        for ff in fitsfiles: shutil.move(ff,v.DESTDIR)

    #TODO: cleanup intermediate files

###################################################################################################
#Utility Functions

def reset_ms ():
    """reset_ms: make a clean start. 
    If FULL_RESET is set (or MS does not exist): does a complete refresh. Removes the measurement set given by the MS variable, 
    If FULL_RESET is not set, then simply clears out flagsets raised during calibration."""
    if is_true("$FULL_RESET"):
        pyrap.tables.addImagingColumns(MS)
        # add a bitflag column (used by MeqTrees)
        ms.addbitflagcol();
        # use the flag-ms.py tool to copy the FLAG column to the "legacy" bitflag
        ms.flagms("-Y +L -f legacy -c")
    else:
        info("$MS exists, reusing (but clearing calibration-related flags). Use FULL_RESET=1 to unpack fresh MS from tarball");
        ms.flagms("-r threshold,fmthreshold,aoflagger,stefcal");

def copy_corrected_data_to_data():
    """
    Copy the contents of the CORRECTED_DATA column to the DATA column
    """
    ms.copycol(msname="$MS",fromcol="CORRECTED_DATA",tocol="DATA")

def copy_corrected_data_to_data_all():
    info("########## copying data to corrected data column step %i"%v.STEP)
    pper("MS",copy_corrected_data_to_data)

def copy_data_to_corrected_data():
    """
    Copy the contents of the DATA column to the CORRECTED_DATA column
    """
    ms.copycol(msname="$MS",fromcol="DATA",tocol="CORRECTED_DATA")

####################################################################################################
##Imaging Functions

def run_wsclean(msFile, wscleanDict):
    wscleanStr='%s '%WSCLEAN_PATH
    for key, value in wscleanDict.iteritems():
        if value is None:
            wscleanStr+='-%s '%(key)
        else:
            wscleanStr+='-%s %s '%(key,str(value))
    wscleanStr+='%s'%(msFile)
    info("making deconvolved image using WSCLEAN: %s"%wscleanStr)
    x.sh(wscleanStr)

###################################################################################################
#D-Projection Function

def generate_dproj_matrix(fitsFile,dprojFn='paper32_rsa_dproj'):
    x.sh('%s/genDproj.py %s -o %s'%(SCRIPTS_DIR,fitsFile,dprojFn))


