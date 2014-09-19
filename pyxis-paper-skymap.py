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

SCRIPTS_DIR='/home/foster/PAPER/pipelines/paper32imaging/scripts'

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
v.MS_List = glob.glob('/home/foster/PAPER/pipelines/paper32imaging/data/zen.*.uvcRREM.MS')

ms.DDID = 0
ms.FIELD = 0

## destination directory for plots, images, etc.
v.DESTDIR_Template = '${OUTDIR>/}output-plots${-stage<STAGE}'
## base filename for these files
v.OUTFILE_Template = '${DESTDIR>/}${MS:BASE}${_s<STEP}${_<LABEL}'

lsm.PYBDSM_POLARIZED = True 

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
#selfcal/imaging/source-finding pipeline function

def cal_pipeline(goto_step=0.):
    #TODO:rfi flagging

    if not os.path.isdir(v.DESTDIR):
        os.makedirs(v.DESTDIR)

    #copy data to corrected data column to run clean
    if goto_step <= 0:
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
            'scale': 330./3600.,
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
            'scale': 330./3600.,
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
            'scale': 330./3600.,
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
            'scale': 330./3600.,
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
# HEALPIX skymap pipeline function

def map_pipeline(goto_step=2):

    useStep=4   #use this step ID for selecting FITS files
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
            beamStr='paper32_%.2fMHz_bw%04.1fMHz_%iarcsec_%ipx.fits'%(hdr['CRVAL3']*10e-7,hdr['CDELT3']*10e-7,int(abs(hdr['CDELT1']*3600.)),int(hdr['NAXIS1']))
            hdulist.close()
            info('########## Using %s for beam template'%beamStr)
            if not os.path.isfile(v.DESTDIR+'/'+beamStr):
                info("########## generating interpolated beam from template FITS, outputing %s step %i"%(v.DESTDIR+'/'+beamStr,v.STEP))
                x.sh('%s/genInterpBeam.py -i --xpol=%s --ypol=%s --ofn=%s %s'%(SCRIPTS_DIR,fekoXpol,fekoYpol,beamStr,fitsfiles[0]))

    #apply inverted beam to Stokes image
    if goto_step <= 2.5:
        v.STEP=4
        info("########## applying beam correction step %i"%v.STEP)
        ##fitsfiles=glob.glob(v.DESTDIR+'/*s%i.restored.fits'%(useStep))
        #fitsXX=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XX-image.fits'%(useStep)))
        #fitsYY=sorted(glob.glob(v.DESTDIR+'/*s%i-*-YY-image.fits'%(useStep)))
        #fitsXYr=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XY-image.fits'%(useStep)))
        #fitsXYi=sorted(glob.glob(v.DESTDIR+'/*s%i-*-XYi-image.fits'%(useStep)))
        #for fxx,fyy,fxyr,fxyi in zip(fitsXX,fitsYY,fitsXYr,fitsXYi):
        #    ofn=fxx.split('/')[-1]
        #    ofn=ofn.split('XX')[0]+'stokes.fits'
        #    info("################ D-Projection: %s"%(ofn))
        #    x.sh('%s/applyDproj.py -d %s --xx=%s --yy=%s --xyr=%s --xyi=%s -o %s'%(SCRIPTS_DIR,dprojStr+'.npy',fxx,fyy,fxyr,fxyi,ofn))

    #apply inverted beam
    #make multifreq fits file
    #run pybdsm on fits file to produce final sky model
    #convert image to healpix
    #sum up maps

#    #TODO:fits files are per freq, per complex polarization, should I run this in the skymap pipeline?
#    #derive an final sky model from the clean'd image
#    if goto_step <= 4.5:
#        v.STEP=4
#        fitsfiles=glob.glob(v.DESTDIR+'/*_s4-image.fits')
#        for ff in sorted(fitsfiles):
#            olsm=ff.split('-image.fits')[0]+'.lsm'
#            lsm.pybdsm_search(image=ff,output=olsm,thresh_pix=12,thresh_isl=5)
#            lsm.tigger_convert(olsm)
#
#    #healpix skymap
#    if goto_step <= 2.:
#        v.STEP=4
#        info("########## generating healpix maps step %i"%v.STEP)
#        #fitsfiles=glob.glob(v.DESTDIR+'/*.dirty.dproj.fits')
#        fitsfiles=glob.glob(v.DESTDIR+'/*s%s.restored.dproj.fits'%(useStep))
#        for ff in fitsfiles:
#            info("################ %s"%ff)
#            for sid in ['I','Q','U','V']:
#                hpxfn=ff.split('.fits')[0]+'.s%s.hpx'%sid
#                #nside: 512 -> 7 arcmin res
#                #nside: 1024 -> 3.5 arcmin res
#                x.sh('/home/foster/PAPER/scripts/mk_map_mod.py %s -n --min_alt=30 -S %s --nside=512 -m %s'%(ff,sid,hpxfn))
#
#    #summed healpix map
#    if goto_step <= 3.:
#        v.STEP=5
#        for sid in ['I','Q','U','V']:
#            #x.sh('/home/foster/PAPER/scripts/save_map.py -s %s %s'%(v.DESTDIR+'/paper_32_combined_%s.hpx',v.DESTDIR+'/*.dirty.dproj.s%s.hpx')%(sid,sid))
#            x.sh('/home/foster/PAPER/scripts/save_map.py -s %s %s'%(v.DESTDIR+'/paper_32_combined_%s_s%i.hpx'%(sid,useStep),v.DESTDIR+'/*s%i.restored.dproj.s%s.hpx')%(useStep,sid))

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


