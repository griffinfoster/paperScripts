#!/usr/bin/env python
"""
Plot a HEALPIX FITS map
"""

import sys,os
import numpy as n

import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as p
from mpl_toolkits.basemap import Basemap

import healpy as h

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] FITS_MAP')
    o.set_description(__doc__)
    o.add_option('-m','--mode',dest='mode',default='H',
        help='Scale mode: H(histogram), I(linear), L(log), S(sqrt); default: H')
    o.add_option('-p','--proj',dest='proj',default='CG',
        help='Projection of map, E(cliptic), G(alactic), C(equatorial), default: CG')
    o.add_option('-t','--title',dest='title',default='',
        help='Map title, default: Sky Map')
    o.add_option('-s','--savefig',dest='savefig',default=None,
        help='Save map to file')
    o.add_option('--max',dest='max',default=None,type='float',
        help='Max flux value, default: None')
    o.add_option('--min',dest='min',default=None,type='float',
        help='Min flux value, default: None')
    o.add_option('--dec_min',dest='dec_min',default=None,type='float',
        help='Min declination to plot, in degrees, default: None')
    o.add_option('--dec_max',dest='dec_max',default=None,type='float',
        help='Max declination to plot, in degrees, default: None')
    o.add_option('-w','--weight',dest='weight_map',action='store_true',
        help='Plot the sample weighting as a map instead of the data')
    o.add_option('--cmap',dest='cmap',default='jet',
        help='Set a matplotlib colormap(jet,hot,gist_heat,gist_earth,PRGn,RdYlBu,spectral), default:jet')
    o.add_option('--cbar',dest='cbar',default=True,action='store_false',
        help='Disable colorbar')
    o.add_option('--fill',dest='fill',default=False,action='store_true',
        help='FIll UNSEEN/NaN values with zeros')
    o.add_option('--rot',dest='rot',default=None,type='str',
        help='Rotate the map such that the center is at (lon,lat)')
    opts, args = o.parse_args(sys.argv[1:])

    cbar=opts.cbar
    unit='Flux (Jy/Beam)'
    if opts.weight_map: unit='Samples'
    if opts.rot:
        lon,lat=opts.rot.split(',')
        rot=(float(lon),float(lat),0.)
    else: rot=None

    xsize=5000

    m=None
    w=None
    for fn in args:
        print 'Opening:',fn
        if m is None: m,w,hdr=h.read_map(fn,field=(0,1),h=True)
        else:
            m0,w0,hdr=h.read_map(fn,field=(0,1),h=True)
            m+=m0
            w+=w0
    #set map projection
    coord_proj=opts.proj.upper()

    m/=w #divide by the pixel weights
    
    print 'Map :: min=%f :: max=%f'%(n.nanmin(m), n.nanmax(m))

    ma_m=n.ma.masked_array(m,n.isnan(m))

    if opts.weight_map: m=w

    if opts.mode.lower()=='s':m=n.sqrt(m-n.min(n.nan_to_num(m)))

    #mask out declination regions
    nside=h.pixelfunc.get_nside(m)
    if opts.dec_min is None: dec_min=180.
    else: dec_min=90.-opts.dec_min
    if opts.dec_max is None: dec_max=0.
    else: dec_max=90.-opts.dec_max
    theta_min=(dec_min/180.)*n.pi
    theta_max=(dec_max/180.)*n.pi
    ring_min=h.pixelfunc.ang2pix(nside,theta_min,0.)
    ring_max=h.pixelfunc.ang2pix(nside,theta_max,0.)
    m[ring_min:]=h.pixelfunc.UNSEEN
    m[:ring_max]=h.pixelfunc.UNSEEN

    #replace nan values with zero
    if opts.fill:
        minVal=opts.max
        m+=minVal
        m=n.nan_to_num(m)
        m-=minVal

    fig=p.figure(1,figsize=(12,8))

    if opts.mode.lower()=='h':
        h.mollview(m,coord=coord_proj,norm='hist',max=opts.max,min=opts.min,cmap=p.get_cmap(opts.cmap),cbar=cbar,fig=1,unit=unit,title=opts.title,xsize=xsize,rot=rot)
    elif opts.mode.lower()=='i':
        h.mollview(m,coord=coord_proj,max=opts.max,min=opts.min,cmap=p.get_cmap(opts.cmap),cbar=cbar,fig=1,unit=unit,title=opts.title,xsize=xsize,rot=rot)
    elif opts.mode.lower()=='l':
        h.mollview(m,coord=coord_proj,norm='log',cmap=p.get_cmap(opts.cmap),cbar=cbar,fig=1,unit=unit,title=opts.title,xsize=xsize,rot=rot)
    elif opts.mode.lower()=='s':
        h.mollview(m,coord=coord_proj,max=opts.max,min=opts.min,cmap=p.get_cmap(opts.cmap),cbar=cbar,fig=1,unit=unit,title=opts.title,xsize=xsize)

    #set grid lines
    h.graticule(dpar=5,dmer=5)

    if not(opts.savefig is None):
        p.savefig(opts.savefig)
    else: p.show()
    
