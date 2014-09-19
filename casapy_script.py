import sys,os
import cmath
from casa import table as tb
import numpy as np
from ephem import *
#import ephem
import fmt
import math
import matplotlib
import matplotlib.pyplot as p
#
# process a snapshot of PAPER data
#

def make_primary_beams(image_name,lst,stokes_choice,beam_filename):
#
# define the frequencies of the simulated beams
# 
 freq_beams = np.zeros(11)
 freq_beams[0] = 100e6
 for j in range(1,11):
  freq_beams[j] = freq_beams[j - 1] + 10e6
#
# read the input cube
#
 tb.open(image_name)
 q_cube = tb.getcol('map')
 tb.close()
#  ra = np.ndarray(shape=(image.shape[0],image.shape[1]),dtype=float)
#  dec = np.ndarray(shape=(image.shape[0],image.shape[1]),dtype=float)
 ia.open(image_name)
 summary = ia.summary()		# read the image summary
 cube_shape = summary['shape']		# read the image shape: RA, DEC, Stokes, Freq
 ra = np.ndarray(shape=(cube_shape[0],cube_shape[1]),dtype=float)
 dec = np.ndarray(shape=(cube_shape[0],cube_shape[1]),dtype=float)
 nchan = cube_shape[3]			# number of channels in the cube
 start_freq = summary['refval'][3]	# start frequencies in Hz
 df = summary['incr'][3]		# frequency increment in Hz
#
# ra and dec will contain the RA and DEC corresponding to the pixel values in the image
#
#
 for j in range(0,cube_shape[0]):
  for k in range(0,cube_shape[1]):
   a=ia.toworld([j,k,0,0])	# a dictionary is returned with the world coordinates of that pixel
   b=a['numeric']		# the array is extracted from the dictionary a 
   ra[j,k] = b[0]		# save the RA for pixel j,k
   dec[j,k] = b[1]		# save the DEC for pixel j,k
#   print ra[j,k] * 12/np.pi,dec[j,k] * 180/np.pi,j,k
#
 print 'RA and DEC calculated'
 ia.close()
#
#
#
# read the beams
#  
 fekoX=fmt.FEKO('/home/gianni/PAPER/beams/fitBeam/data/PAPER_FF_X.ffe')
 fekoY=fmt.FEKO('/home/gianni/PAPER/beams/fitBeam/data/PAPER_FF_Y.ffe')   
 feko_xpol=fekoX.fields[0]
 feko_ypol=fekoY.fields[0]
 phi=feko_xpol.phi*np.pi/180.					# phi is the azimuth
 theta=feko_xpol.theta*np.pi/180.				# theta is the zenith angle
 theta = np.pi/2 - theta					# pyephem wants the elevation rather than the zenith angle
 ra_beams = np.ndarray(shape=(phi.shape[0]),dtype=float)	# array of RA corresponding to (phi,theta)
 dec_beams = np.ndarray(shape=(phi.shape[0]),dtype=float)	# array of DEC corresponding to (phi,theta)
#
# compute complex beams
#
 gxx=feko_xpol.etheta*np.conj(feko_xpol.etheta)+feko_xpol.ephi*np.conj(feko_xpol.ephi)
#
# define an array that will contain all the simulated beams
#
 beams = np.ndarray(shape=(gxx.shape[0],4,freq_beams.shape[0]),dtype=float)
#
# read all the beam models and save them in the beams array
#
 for j in range(0,freq_beams.shape[0]):
  feko_xpol = fekoX.fields[j]		# these two lines give an error, I need to check with Griffin how to fix it
  feko_ypol = fekoY.fields[j]
#
  gxx = feko_xpol.etheta*np.conj(feko_xpol.etheta)+feko_xpol.ephi*np.conj(feko_xpol.ephi)
  gyy = feko_ypol.etheta*np.conj(feko_ypol.etheta)+feko_ypol.ephi*np.conj(feko_ypol.ephi)
  gxy = feko_xpol.etheta*np.conj(feko_ypol.etheta)+feko_xpol.ephi*np.conj(feko_ypol.ephi)
  gyx = feko_ypol.etheta*np.conj(feko_xpol.etheta)+feko_ypol.ephi*np.conj(feko_xpol.ephi)
#
# make the stokes beams
#
  beams[:,0,j] = (gxx+gyy).real
#  beams[:,0,j] = beams[:,0,j] / np.max(beams[:,0,j])	# normalize the beams to be 1 at zenith
  beams[:,1,j] = (gxx-gyy).real
  beams[:,2,j] = (gxy+gyx).real
  beams[:,3,j] = (gxy-gyx).imag
#
 norm_beam = np.max(beams[:,0,5]) 			# beam peak at 150 MHz
 beams[:,0,:] = beams[:,0,:] / norm_beam		# normalize the beams to be 1 at zenith at 150 MHz
#
 print norm_beam,np.max(beams[:,0,:])
 print 'Beams read'
#
# bring the beam to RA,DEC coordinates
#
# Create an observer 
#
 paper = Observer()
#
# Set the observer at the Karoo
#
 paper.lat, paper.long, paper.elevation = '-30:43:17', '21:25:40.08', 0.0
# j0 = ephem.julian_date(0)			# invoked this way if    import ephem
 j0 = julian_date(0)				# invoked this way if    from ephem import *
# paper.date = float(lst)
 paper.date = float(lst) - j0 + 5./60/24	# I think I need this. At http://stackoverflow.com/questions/8962426/convert-topocentric-coordinates-azimuth-elevation-to-equatorial-coordinates
  						# they seem to suggest that this is needed in order to get the right date and I seem to get the right 
						# RA if I include this. The factor 5./60/24 needs to be added as the lst reported in the filename refers to
						# the beginning of the integration which is ~10 min long
 for j in range(0,ra_beams.shape[0]):
  a = paper.radec_of(phi[j],theta[j])
  ra_beams[j] = a[0]				# RA is in radians
  dec_beams[j] = a[1]				# DEC is in radians
#
#
# now interpolate the beams in frequency
#
 interp_beam = np.ndarray(shape=(beams.shape[0],beams.shape[1],nchan),dtype=float)
 cube_freq = start_freq
 for chan in range(0,nchan):
  a = np.max(q_cube[:,:,0,chan,0])
  b = np.min(q_cube[:,:,0,chan,0])
  if (a != 0 and b != 0):			# if the image is not empty, then proceed
   freq_dist = np.abs(cube_freq - freq_beams)
   freq_dist_s = np.sort(freq_dist)
   w = np.where(freq_dist == freq_dist_s[0])
   if freq_dist_s[0] == 0:
#
# if the beam is simulated at the exact frequency channel, then do not interpolate
#
    for j in range(0,4):
     interp_beam[:,j,chan] = beams[:,j,w[0][0]]
#
# if they are not, perform a weighted average of the two closest beams in frequency. The weights are the inverse of the frequency distance squared
#
   else:
    w1 = np.where(freq_dist == freq_dist_s[1])
    for j in range(0,4):
     interp_beam[:,j,chan] = (beams[:,j,w[0][0]] * freq_dist_s[0]**(-2) + beams[:,j,w1[0][0]] * freq_dist_s[1]**(-2)) / (freq_dist_s[0]**(-2) + freq_dist_s[1]**(-2))
#
  cube_freq = cube_freq + df
#
 print 'Beams interpolated in frequency'
#
# now interpolate the beam at the observed RA,DEC
#
 interp_beam_maps_q = np.ndarray(shape=(ra.shape[0],ra.shape[1],1,nchan,1),dtype=float) 
#
 for j in range(0,ra.shape[0]):
  for k in range(0,ra.shape[1]):
#
# interpolating amongst the three closest points
#
#   x = np.cos(ra[j,k])*np.cos(ra_beams)*np.cos(dec[j,k])*np.cos(dec_beams)
#   y = np.sin(ra[j,k])*np.sin(ra_beams)*np.cos(dec[j,k])*np.cos(dec_beams)
#   z = np.sin(dec[j,k])*np.sin(dec_beams)
#   dist = np.sqrt(x**2 + y**2 + z**2)
#
   dist = np.sqrt((ra[j,k] - ra_beams)**2 + (dec[j,k] - dec_beams)**2)
   dist_s = np.sort(dist)
   w0 = np.where(dist == dist_s[0])
   w1 = np.where(dist == dist_s[1])
   w2 = np.where(dist == dist_s[2])
#
   interp_beam_maps_q[j,k,0,:,0] = interp_beam[w0[0][0],stokes_choice,:] / dist_s[0] + interp_beam[w1[0][0],stokes_choice,:] / dist_s[1] + interp_beam[w2[0][0],stokes_choice,:] / dist_s[2]
   interp_beam_maps_q[j,k,0,:,0] = interp_beam_maps_q[j,k,0,:,0] / (dist_s[0]**(-1) + dist_s[1]**(-1) + dist_s[2]**(-1))
#
# nearest neighbour interpolation
#
#   dist = np.sqrt((ra[j,k] - ra_beams)**2 + (dec[j,k] - dec_beams)**2)
#   dist_s = np.sort(dist)
#   w0 = np.where(dist == dist_s[0])
#   interp_beam_maps_q[j,k,0,:,0] = interp_beam[w0[0][0],stokes_choice,:]
#
 print 'Beams interpolated in angle'
# 
# store the beams into an image
#
# beam_filename = image_name.strip('.image')
 cmd = 'rm -rf ' + beam_filename + '.beams'
 os.system(cmd)
 cmd = 'cp -r ' + image_name + ' ' + beam_filename + '.beams'
 os.system(cmd)
#
 tb.open(beam_filename + '.beams',nomodify=False)
 tb.putcol('map',interp_beam_maps_q)
 tb.close()
 tb.done()
 ia.open(beam_filename + '.beams')
 ia.tofits(beam_filename + '_beams.fits',overwrite=true)
 ia.close()

 




















 


#if __name__ == '__main__':
# from optparse import OptionParser
# o = OptionParser()
# o.set_usage('%prog [options] input_MS_file out_MS_file')
# o.set_description(__doc__)
# opts, args = o.parse_args(sys.argv[1:])
#
filename = 'zen.2455819.50285.uvcRREM.MS'
filename = 'zen.2455819.50981.uvcRREM.MS'
filename = 'zen.2455819.51677.uvcRREM.MS'
filename = 'zen.2455819.52373.uvcRREM.MS'
filename = 'zen.2455819.53069.uvcRREM.MS'
filename = 'zen.2455819.53765.uvcRREM.MS'
#
filenames = ['zen.2455819.50285.uvcRREM.MS','zen.2455819.50981.uvcRREM.MS','zen.2455819.51677.uvcRREM.MS','zen.2455819.52373.uvcRREM.MS','zen.2455819.53069.uvcRREM.MS','zen.2455819.53765.uvcRREM.MS','zen.2455819.54461.uvcRREM.MS','zen.2455819.55157.uvcRREM.MS','zen.2455819.55853.uvcRREM.MS','zen.2455819.56548.uvcRREM.MS','zen.2455819.57244.uvcRREM.MS','zen.2455819.57940.uvcRREM.MS','zen.2455819.58636.uvcRREM.MS','zen.2455819.59332.uvcRREM.MS','zen.2455819.60028.uvcRREM.MS','zen.2455819.60724.uvcRREM.MS','zen.2455819.61420.uvcRREM.MS','zen.2455819.62116.uvcRREM.MS','zen.2455819.64204.uvcRREM.MS','zen.2455819.64900.uvcRREM.MS','zen.2455819.65596.uvcRREM.MS','zen.2455819.66292.uvcRREM.MS','zen.2455819.66988.uvcRREM.MS','zen.2455819.67684.uvcRREM.MS','zen.2455819.68380.uvcRREM.MS','zen.2455819.69075.uvcRREM.MS','zen.2455819.69771.uvcRREM.MS']
#
# initialize the frequencies at which the beams are computed from the FEKO simulations
#
#freq_beams = np.zeros(11)
#freq_beams[0] = 100e6
#for i in range(1,11):
# freq_beams[i] = freq_beams[i-1] + 10e6
#
#
first_pass = 0
bandpass_normalize = 0
recal_image = 0
selfcal = 0
pol_cube = 0
rotate = 0
cube = 0
beam_correct = 1
rmsynthesis = 0
peel = 0
#
#for i in range(0,len(filenames)):
#for i in range(0,1):
for i in range(11,12):			# this frame has Fornax A transiting
#
 filename = filenames[i]
 a = filename.strip('zen.')
 lst = a.strip('.uvcRREM.MS')
#
# make_primary_beams('cube_I_' + lst + '.image',lst,0,'beams_StokesI_10MHz')
# make_primary_beams('cube_I_' + lst + '.image',lst,1,'beams_StokesQ_10MHz')
# make_primary_beams('cube_I_' + lst + '.image',lst,2,'beams_StokesU_10MHz')
# make_primary_beams('cube_I_' + lst + '.image',lst,3,'beams_StokesV_10MHz')
#
#
#
 if first_pass == 1:
#
# specific PAPER options
#
#cmd = 'rm -rf '+filename
#os.system(cmd)
#cmd = 'python /home/gianni/Python_code/swap_column.py /home/gianni/PAPER/psa32/data/Original/'+filename+' '+filename 
#os.system(cmd)
#ft(vis=filename,model='First_model/2455819.50285.model',usescratch=True)
#
  cmd = 'rm -rf ' + filename
  os.system(cmd)
  cmd = 'cp -r /home/gianni/PAPER/psa32/data/Original/'+ filename + ' ' + filename
  os.system(cmd)
  tb.open(filename,nomodify=False)
  data = tb.getcol('DATA')
  tb.putcol('CORRECTED_DATA',data)				# the CORRECTED_DATA column is initialized to DATA
  tb.close()
  tb.done()
#
  tflagdata(vis=filename,mode='manual',spw='0:167~202',action='apply',datacolumn='DATA')
  flagdata(vis=filename,mode='clip',clipminmax=[0,300],correlation='ABS_XY,YX',action='apply',datacolumn='DATA')
  flagdata(vis=filename,mode='clip',clipminmax=[0,1000],correlation='ABS_XX,YY',action='apply',datacolumn='DATA')
#
  out_fits  = lst+'.fits'
  out_fits_res  = lst+'_residual.fits'
  out_fits_model  = lst+'_model.fits'
#
  cmd = 'rm -rf '+ lst + '.image ' + lst + '.model ' + lst + '.residual ' + lst + '.flux' + lst + '.psf'
  os.system(cmd)
#
# make an image
#
  myphasecenter='J2000 20h00m00s -30d43m17.5s'
#clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=0000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=64,gain=0.2,stokes='I',phasecenter=myphasecenter)
  clean(vis=filename,imagename=lst,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
  ia.open(lst+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
  ia.open(lst+'.residual')
  ia.tofits(out_fits_res,overwrite=true)
  ia.close()
  ia.open(lst+'.model')
  ia.tofits(out_fits_model,overwrite=true)
  ia.close()
  ia.open(lst+'.psf')
  ia.tofits('psf.fits',overwrite=true)
  ia.close()
#
  bandpass_table = 'b_'+lst+'.cal'
  bandpass(vis=filename,solint='inf',combine='scan',caltable=bandpass_table)
#
#
#
#
#
 if bandpass_normalize == 1:
  bandpass_table = 'b_'+lst+'.cal'
#
  tb.open(bandpass_table,nomodify=False)
  gain=tb.getcol('CPARAM')
  gain_norm=gain
# 
# for each antenna and for each channel of the bandpass I divide out by the modulo of the complex number
#  
  for i in range(0,gain.shape[2]):
   for j in range(0,gain.shape[1]):
    a = gain[0,j,i]
#
# if the real part of the antenna gain is set to 1 it means that that antenna a/o channel is flag, so don't bother looking at it
#   
    if a.real != 1:
     gain_norm[0,j,i] = gain_norm[0,j,i]/abs(gain_norm[0,j,i]) 
     gain_norm[1,j,i] = gain_norm[1,j,i]/abs(gain_norm[1,j,i])
#     print abs(gain_norm[i,j,0]),abs(gain[i,j,0]),i,j
#
#
# put back the normalized bandpass
#    
  tb.putcol('CPARAM',gain_norm)
  tb.close()
  tb.done()
#
  applycal(vis=filename,gaintable=bandpass_table)
#
  tb.open(filename,nomodify=False)
  corrected = tb.getcol('CORRECTED_DATA')
  tb.putcol('DATA',corrected)					# now the CORRECTED_DATA column contains the recalibrated visibilities
  tb.close()
  tb.done()
#
#
#
#
 if recal_image == 1:
  out_image = lst+'_recal'
  out_fits  = out_image+'.fits'
  out_fits_res  = out_image+'_residual.fits'
  out_fits_model  = out_image+'_model.fits'
#
  cmd = 'rm -rf '+ out_image + '.image ' + out_image + '.model ' + out_image + '.residual ' + out_image + '.flux' + out_image + '.psf'
  os.system(cmd)
#
# make an image
#
  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
#  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[1024,1024],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
  ia.open(out_image+'.residual')
  ia.tofits(out_fits_res,overwrite=true)
  ia.close()
  ia.open(out_image+'.model')
  ia.tofits(out_fits_model,overwrite=true)
  ia.close()
  ia.open(out_image+'.psf')
  ia.tofits('psf.fits',overwrite=true)
  ia.close()
#
  cmd = 'python /home/gianni/Python_code/flagging_residuals_GB.py ' + filename + ' --showfit -n 3e4 -d --nomodel'
#  os.system(cmd)
#
#
#python ~/Python_code/flagging_residuals_GB.py zen.2455819.50285.uvcRREM.MS -n 3e4 -d --nomodel
#
#
#
#
 if selfcal == 1:
  ft(vis=filename,model=lst+'_recal.model',usescratch=True)
  bandpass_table = 'b1_'+lst+'.cal'
#  bandpass(vis=filename,solint='inf',combine='scan',caltable=bandpass_table)
#  applycal(vis=filename,gaintable=bandpass_table)
  gaincal_table = 'g_'+lst+'.cal'
  gaincal(vis=filename,gaintype='G',caltable=gaincal_table,calmode='p')
  applycal(vis=filename,gaintable=gaincal_table)
#
#
#
#
#
 if pol_cube == 1:
  stokes_par = ['I','Q','U','V']
  cellsize = '12arcmin'
  npix = 512
  cellsize = '6arcmin'
  npix = 512
#
  for j in range (0,len(stokes_par)):  
   out_image = 'pol_cube_' + lst + '_' + stokes_par[j]
   out_fits  = out_image+'.fits'
   out_fits_res  = out_image+'_residual.fits'
   out_fits_model  = out_image+'_model.fits'
#
# make a polarization cube
#
#
   cmd = 'rm -rf '+out_image+'*'
   os.system(cmd)
#   clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='3Jy',niter=0000,mode='channel',nchan=203,weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes=stokes_par[j],uvtaper=True,outertaper=['30arcmin'])
   if stokes_par[j] == 'I':
    clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='50Jy',niter=100,mode='channel',nchan=203,weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes=stokes_par[j],uvtaper=True,outertaper=['30arcmin'])
   else: clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='20Jy',niter=0000,mode='channel',nchan=203,weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes=stokes_par[j],uvtaper=True,outertaper=['30arcmin'])
#
   ia.open(out_image+'.image')
   ia.tofits(out_fits,overwrite=true)
   ia.close()
   ia.open(out_image+'.residual')
   ia.tofits(out_fits_res,overwrite=true)
   ia.close()
   ia.open(out_image+'.model')
   ia.tofits(out_fits_model,overwrite=true)
   ia.close()
#
  out_image = 'pol_cube_' + lst
  out_fits  = out_image+'.fits'
  out_fits_res  = out_image+'_residual.fits'
  out_fits_model  = out_image+'_model.fits'
#
# make a polarization cube
#
#
  cmd = 'rm -rf ' + out_image + '.model '+ out_image + '.image '+ out_image + '.residual '
  os.system(cmd)
  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='10Jy',niter=10000,mode='mfs',weighting='uniform',cell=[cellsize,cellsize],imsize=[2*npix,2*npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='IQUV',uvtaper=True,outertaper=['30arcmin'])
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
  ia.open(out_image+'.residual')
  ia.tofits(out_fits_res,overwrite=true)
  ia.close()
  ia.open(out_image+'.model')
  ia.tofits(out_fits_model,overwrite=true)
  ia.close()
#
#
#
#
#
 if rotate == 1:
  cellsize = '12arcmin'
  npix = 512
#
  out_image = lst+'_rotated'
  out_fits  = out_image+'.fits'
  out_fits_res  = out_image+'_residual.fits'
  out_fits_model  = out_image+'_model.fits'
#
  cmd = 'rm -rf '+ out_image + '.image ' + out_image + '.model ' + out_image + '.residual ' + out_image + '.flux' + out_image + '.psf'
  os.system(cmd)
#
# make an image
#
  myphasecenter='J2000 3h00m00s -30d43m17.5s'
  myphasecenter='J2000 1h15m00s -30d43m17.5s'
  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=128,gain=0.2,stokes='I',phasecenter=myphasecenter)
  #clean(vis=filename,imagename=lst,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
  ia.open(out_image + '.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
  ia.open(out_image + '.residual')
  ia.tofits(out_fits_res,overwrite=true)
  ia.close()
  ia.open(out_image +'.model')
  ia.tofits(out_fits_model,overwrite=true)
  ia.close()
  ia.open(out_image + '.psf')
  ia.tofits('psf.fits',overwrite=true)
  ia.close()
#
#
#
#
#
 if cube == 1:
  cellsize = '3arcmin'
  npix = 1024
#  cellsize = '6arcmin'
#  npix = 512
#  cellsize = '24arcmin'
#  npix = 128
#
  chan1 = 0
  chan2 = 9
  midchan = 5
  dchan = 10
  tot_chan = 10
#
  chan1 = 20
  chan2 = 29
  midchan = 10
  dchan = 20
  tot_chan = 8
#
  out_image = 'cube_I_' + lst
  cmd = 'rm -rf ' + out_image + '*'
  os.system(cmd)
  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='10Jy',niter=000,mode='channel',nchan=tot_chan,start=1,width=dchan,weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2)
  tb.open(out_image + '.image',nomodify=False)
  cube = tb.getcol('map')
  tb.close()
#
  k = 0
  for nchan in range(0,tot_chan):
#  for nchan in range(0,0):
   out_image1 = 'temp_' + lst + '_' + str(midchan)
#   out_fits  = out_image+'.fits'
#   out_fits_res  = out_image+'_residual.fits'
#   out_fits_model  = out_image+'_model.fits'
#
# make a cube
#
#
   cmd = 'rm -rf ' + out_image1 + '*'
   os.system(cmd)
   clean(vis=filename,imagename=out_image1,imagermode='csclean',psfmode='clark',threshold='10Jy',niter=10000,mode='mfs',weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,spw='0:' + str(chan1) + '~' + str(chan2),uvtaper=True,outertaper=['30arcmin'])
   tb.open(out_image1 + '.image',nomodify=False)
   temp = tb.getcol('map')
   tb.close()
   cube[:,:,0,k,0] = temp[:,:,0,0,0]
   cmd = 'rm -rf ' + out_image1 + '*'
   os.system(cmd)
#
   chan1 = chan1 + dchan
   chan2 = chan2 + dchan
   midchan = midchan + dchan
   k += 1
   
#
  tb.open(out_image + '.image',nomodify=False)
  tb.putcol('map',cube)
  tb.close()
  out_fits  = out_image + '.fits'
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
#   ia.open(out_image+'.residual')
#   ia.tofits(out_fits_res,overwrite=true)
#   ia.close()
#   ia.open(out_image+'.model')
#   ia.tofits(out_fits_model,overwrite=true)
#   ia.close()
#
#
#  make_primary_beams(out_image + '.image',lst,0,'beams_StokesI_10MHz')

#
  out_image = 'cube_' + lst
  out_fits  = out_image+'.fits'
  out_fits_res  = out_image+'_residual.fits'
  out_fits_model  = out_image+'_model.fits'
#
  cmd = 'rm -rf '+out_image+'*'
#  os.system(cmd)
#
#  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='IQUV')
##  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='3Jy',niter=10000,mode='mfs',weighting='natural',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='IQUV',uvtaper=True,outertaper=['60arcmin'])
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
  ia.open(out_image+'.residual')
  ia.tofits(out_fits_res,overwrite=true)
  ia.close()
  ia.open(out_image+'.model')
  ia.tofits(out_fits_model,overwrite=true)
  ia.close()
#
#
#
#
#
 if beam_correct == 1:
  cellsize = '3arcmin'
  npix = 1024
#
  inp_image = 'cube_I_' + lst
  out_image = inp_image + '_beam_corrected'
  cmd = 'rm -rf ' + out_image + '.image' 
  os.system(cmd)
  cmd = 'cp -r ' + inp_image + '.image ' + out_image + '.image'
  os.system(cmd)
  tb.open(inp_image + '.image')
  inp_cube = tb.getcol('map')
  tb.close()
#
# read the beams. You need to know upfront that you image was made by binning ~10 MHz (i.e. the image and the beams images have the same channel numbers)
#
  tb.open('beams_StokesI_10MHz.beams')
  beams = tb.getcol('map')
  tb.close()
#  
  for j in range(0,inp_cube.shape[3]):
   inp_cube[:,:,0,j,0] = inp_cube[:,:,0,j,0] / beams[:,:,0,j,0]
#
  tb.open(out_image + '.image',nomodify=False)
  tb.putcol('map',inp_cube)
  tb.close()
  ia.open(out_image + '.image')
  ia.tofits(out_image + '.fits',overwrite=true)
  ia.close()
#
  out_image = 'cube_I_' + lst + '_freq_averaged_beam_corrected'
  cmd = 'rm -rf ' + out_image + '*'
  os.system(cmd)
  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='10Jy',niter=0000,mode='mfs',weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2)
  tb.open(out_image + '.image',nomodify=False)
  image = tb.getcol('map')
  image = image - image
  print np.max(image),np.min(image)
  for j in range(0,inp_cube.shape[3]):
   image[:,:,0,0,0] = image[:,:,0,0,0] + inp_cube[:,:,0,j,0]
  image = image / inp_cube.shape[3]
  print np.max(image),np.min(image)
  tb.putcol('map',image)
  tb.close()
  ia.open(out_image + '.image')
  ia.tofits(out_image + '.fits',overwrite=true)
  ia.close()
#
  out_image = 'cube_I_' + lst + '_freq_averaged'
  cmd = 'rm -rf ' + out_image + '*'
  os.system(cmd)
  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='5Jy',niter=10000,mode='mfs',weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2)
  ia.open(out_image + '.image')
  ia.tofits(out_image + '.fits',overwrite=true)
  ia.close()
#
#
#
#
#
 if rmsynthesis == 1:
#
# RM synthesis parameters
#
  rm_min = -100		# minimum RM in rad/m^2
  rm_step = 1		# RM step in rad/m^2
  rm_frames = 201	# total number of RM frames in the output cube
  rm_step = 5		# RM step in rad/m^2
  rm_frames = 41	# total number of RM frames in the output cube
  c_speed = 299792458.0	# speed of light in meters
#
  cellsize = '6arcmin'
  npix = 512
#
# read the Stokes Q cube
#
  out_image = 'pol_cube_' + lst + '_Q'
  tb.open(out_image+'.image')
  image = tb.getcol('map')
  tb.close()
#
# read the frequency info of the cube
#
  ia.open(out_image+'.image')
  summary = ia.summary()		# read the image summary
  cube_shape = summary['shape']		# read the image shape: RA, DEC, Stokes, Freq
  nchan = cube_shape[3]			# number of channels in the cube
  start_freq = summary['refval'][3]	# start frequencies in Hz
  df = summary['incr'][3]		# frequency increment in Hz
#
#
# prepare the frequency axis
#
#
  cube_freq = np.zeros(nchan)
  flag = np.zeros(nchan)
  freq = start_freq
#
# simply convert the frequencies into lambda^2
#  
  for n in range(0,nchan):
   cube_freq[n] = freq
   a = np.max(image[:,:,0,n,0])
   b = np.min(image[:,:,0,n,0])
   if (a != 0 and b != 0):				# if the image is not empty, then proceed
    flag[n] = 1
#   print cube_freq/1e6,flag 
   freq = freq + df
#
  w = np.where(flag == 1) 
  cube_lambda = np.ndarray(shape=w[0][0],dtype=float)
  cube_lambda = (c_speed/cube_freq[w[0]])**2
#
  cube_lambda = cube_lambda - np.mean(cube_lambda) 	# derotate the polarization angles to the average lambda_2
#
# make a plot of the RMSF
#
  rm_rmsf = np.zeros(101)
  rmsf = np.zeros(101)
  rmsfq = np.zeros(101)
  rmsfu = np.zeros(101)
  rm1 = -50		# rad/m^2
  for j in range(0,rmsf.shape[0]):
   rm_rmsf[j] = rm1
   arg = -2.0 * rm_rmsf[j] * cube_lambda
   rmsf[j] = np.sqrt((np.sum(np.cos(arg)) / cube_lambda.shape[0])**2 + (np.sum(np.sin(arg)) / cube_lambda.shape[0])**2)
   rmsfq[j] = np.sum(np.cos(arg)) / cube_lambda.shape[0]
   rmsfu[j] = np.sum(np.cos(arg)) / cube_lambda.shape[0]
   rm1 += 1
#
#   p.plot(rm_rmsf,rmsf)
#
# make the Stokes Q cube with only the frames that are meaningful 
#
  q_cube = np.ndarray(shape=(cube_shape[0],cube_shape[1],1,cube_lambda.shape[0],1),dtype=float)
  q_cube[:,:,0,:,0] = image[:,:,0,w[0],0]
#
# read the Stokes U cube
#
  out_image = 'pol_cube_' + lst + '_U'
  tb.open(out_image+'.image')
  image = tb.getcol('map')
  tb.close()
#
# make the Stokes U cube with only the frames that are meaningful 
#
  u_cube = np.ndarray(shape=(cube_shape[0],cube_shape[1],1,cube_lambda.shape[0],1),dtype=float)
  u_cube[:,:,0,:,0] = image[:,:,0,w[0],0]
#
# is there a faster and smarter way than the following to make a fake cube?
#
#  out_image = 'P_RM_' + lst
#  cmd = 'rm -rf '+ out_image + '.model'
#  os.system(cmd)
#  clean(vis=filename,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='20Jy',niter=0000,mode='channel',nchan=rm_frames,weighting='uniform',cell=[cellsize,cellsize],imsize=[npix,npix],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
#  cmd = 'cp -r ' + out_image +'.image Q_RM_' + lst +'.image'
#  os.system(cmd)
#  cmd = 'cp -r ' + out_image +'.image U_RM_' + lst +'.image'
#  os.system(cmd)
#  cmd = 'rm -rf ' + out_image +'.flux'
#  os.system(cmd)
#  cmd = 'rm -rf ' + out_image +'.model'
#  os.system(cmd)
#  cmd = 'rm -rf ' + out_image +'.residual'
#  os.system(cmd)
#  cmd = 'rm -rf ' + out_image +'.psf'
#  os.system(cmd)
#
# read the fake cubes
#
  out_image = 'P_RM_' + lst
  tb.open(out_image+'.image')
  p_rm = tb.getcol('map')
  tb.close()
  out_image = 'Q_RM_' + lst
  tb.open(out_image+'.image')
  q_rm = tb.getcol('map')
  tb.close()
  out_image = 'U_RM_' + lst
  tb.open(out_image+'.image')
  u_rm = tb.getcol('map')
  tb.close()
#  print p_rm.shape
#  print q_rm.shape
#  print u_rm.shape
#
#
# now do the RM synthesis
#
#
  p_in = q_cube + u_cube*1j
  rm = rm_min
  for n in range(0,rm_frames):
   arg = -2.0 * rm * cube_lambda * 1j
#   print arg.shape[0]
#
   for j in range(0,q_cube.shape[0]):
    for k in range(0,q_cube.shape[1]):
     a = np.sum(p_in[j,k,0,:,0] * np.exp(arg)) / cube_lambda.shape[0]
     p_rm[j,k,0,n,0] = np.abs(a)
     q_rm[j,k,0,n,0] = a.real
     u_rm[j,k,0,n,0] = a.imag
#
   print rm
   rm = rm + rm_step 
#
#
# write the RM synthesis cubes to disk
#
  out_image = 'P_RM_' + lst
  out_fits  = out_image+'.fits'
  tb.open(out_image+'.image',nomodify=False)
  tb.putcol('map',p_rm)
  tb.close()
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
#
  out_image = 'Q_RM_' + lst
  out_fits  = out_image+'.fits'
  tb.open(out_image+'.image',nomodify=False)
  tb.putcol('map',q_rm)
  tb.close()
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
#
  out_image = 'U_RM_' + lst
  out_fits  = out_image+'.fits'
  tb.open(out_image+'.image',nomodify=False)
  tb.putcol('map',u_rm)
  tb.close()
  ia.open(out_image+'.image')
  ia.tofits(out_fits,overwrite=true)
  ia.close()
















#
#
#
#
#
if peel == 1:
 cl1_name = 'PKS2356-61_'+lst+'.cl'
 cl2_name = 'PicA_'+lst+'.cl'
#  
 file_in = 'test.MS'
 os.system('rm -rf test.MS')
 os.system('cp -r zen.2455819.50285.uvcRREM.MS test.MS')
 cmd = 'rm -rf ' + cl1_name + ' ' + cl2_name
 os.system(cmd)
#
 cl.addcomponent(flux=26, fluxunit='Jy', polarization='Stokes',dir='J2000 5h19m35.4s -45d53m39s', shape='point')
#you can add more components if you wish by calling addcomponent repeatedly with different params
#save it to disk
 cl.rename(cl2_name)
 cl.close()
#
 cl.addcomponent(flux=48, fluxunit='Jy', polarization='Stokes',dir='J2000 23h58m26.7s -60d52m28s', shape='point')
#you can add more components if you wish by calling addcomponent repeatedly with different params
#save it to disk
 cl.rename(cl1_name)
 cl.close()
#
 ft(vis=file_in,complist=cl2_name,usescratch=True)
 uvsub(vis=file_in)
 ft(vis=file_in,complist=cl1_name,usescratch=True)
 uvsub(vis=file_in)
# uvsub(vis=file_in,reverse=True)
#
 out_image = 'test_image'
 out_fits  = out_image+'.fits'
 out_fits_res  = out_image+'_residual.fits'
 out_fits_model  = out_image+'_model.fits'
#
# make an image
#
#
 cmd = 'rm -rf '+out_image+'*'
 os.system(cmd)
 clean(vis=file_in,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='3Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
 ia.open(out_image+'.image')
 ia.tofits(out_fits,overwrite=true)
 ia.close()
 ia.open(out_image+'.residual')
 ia.tofits(out_fits_res,overwrite=true)
 ia.close()
 ia.open(out_image+'.model')
 ia.tofits(out_fits_model,overwrite=true)
 ia.close()
#
# start peeling PKS
#
 ft(vis=file_in,model=out_image+'.model',usescratch=True)
 uvsub(vis=file_in)						# subtract the whole sky model
 ft(vis=file_in,complist=cl1_name,usescratch=True)		# FT the PKS model and add it to the MODEL column
 uvsub(vis=file_in,reverse=True)				# add PKS to the residual visibilities 
#
 tb.open(file_in,nomodify=False)
 data = tb.getcol('DATA')
 corrected = tb.getcol('CORRECTED_DATA')
 tb.putcol('DATA',corrected)					# now the CORRECTED_DATA column contains the residual visibilities + PKS
 tb.close()
 tb.done()
#
# bandpass(vis=file_in,solint='inf',combine='scan',caltable='bandpass.cal')
 caltable = 'g_' + 'PKS2356-61_' + lst + '_phase.cal'
 gaincal(vis=file_in,gaintype='G',combine='scan',caltable=caltable,calmode='p')
 applycal(vis=file_in,gaintable=caltable)
#
#
 out_image = 'test1_image'
 out_fits  = out_image+'.fits'
 out_fits_res  = out_image+'_residual.fits'
 out_fits_model  = out_image+'_model.fits'
#
# make an image
#
#
 cmd = 'rm -rf '+out_image+'*'
 os.system(cmd)
 clean(vis=file_in,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='3Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
 ia.open(out_image+'.image')
 ia.tofits(out_fits,overwrite=true)
 ia.close()
 ia.open(out_image+'.residual')
 ia.tofits(out_fits_res,overwrite=true)
 ia.close()
 ia.open(out_image+'.model')
 ia.tofits(out_fits_model,overwrite=true)
 ia.close()
#
# restore the original columns
#
 tb.open(file_in,nomodify=False)
 tb.putcol('DATA',data)
 tb.putcol('CORRECTED_DATA',corrected)
 tb.close()
 tb.done()
#
# start peeling PicA
#
 ft(vis=file_in,model=out_image+'.model',usescratch=True)
 uvsub(vis=file_in)						# subtract the whole sky model
 ft(vis=file_in,complist=cl2_name,usescratch=True)		# FT the PicA model and add it to the MODEL column
 uvsub(vis=file_in,reverse=True)				# add PicA to the residual visibilities 
#
 tb.open(file_in,nomodify=False)
 data = tb.getcol('DATA')
 corrected = tb.getcol('CORRECTED_DATA')
 tb.putcol('DATA',corrected)					# now the CORRECTED_DATA column contains the residual visibilities + PicA
 tb.close()
 tb.done()
#
 caltable = 'g_' + 'PicA_' + lst + '_phase.cal'
 gaincal(vis=file_in,gaintype='G',combine='scan',caltable=caltable,calmode='p')
 applycal(vis=file_in,gaintable=caltable)
#
#
 out_image = 'test2_image'
 out_fits  = out_image+'.fits'
 out_fits_res  = out_image+'_residual.fits'
 out_fits_model  = out_image+'_model.fits'
#
# make an image
#
#
 cmd = 'rm -rf '+out_image+'*'
 os.system(cmd)
 clean(vis=file_in,imagename=out_image,imagermode='csclean',psfmode='clark',threshold='3Jy',niter=10000,mode='mfs',weighting='uniform',cell=['3arcmin','3arcmin'],imsize=[2048,2048],gridmode='widefield',wprojplanes=32,gain=0.2,stokes='I')
 ia.open(out_image+'.image')
 ia.tofits(out_fits,overwrite=true)
 ia.close()
 ia.open(out_image+'.residual')
 ia.tofits(out_fits_res,overwrite=true)
 ia.close()
 ia.open(out_image+'.model')
 ia.tofits(out_fits_model,overwrite=true)
 ia.close()
#
# restore the original columns
#
 tb.open(file_in,nomodify=False)
 tb.putcol('DATA',data)
 tb.putcol('CORRECTED_DATA',corrected)
 tb.close()
 tb.done()
