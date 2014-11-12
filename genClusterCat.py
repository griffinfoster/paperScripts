#!/usr/bin/env python
"""Generate a clustered catalogue from pyBDSM Gaul files"""

import numpy as np
import cPickle as pkl
import os,sys
import ephem

def combineSrcs(srcDict0,srcDict1):
    """Combine two source dictionaries and return a single dictionary
    assumes both dictionaries have the same keys"""
    combineSrcDict={}
    for key in srcDict0.iterkeys():
        if key.startswith('MJD'): combineSrcDict[key]=srcDict0[key]
        elif key.startswith('S_Code'): combineSrcDict[key]=srcDict0[key]
        elif key in ['Total_flux','E_Total_flux','Peak_flux','E_Peak_flux']:
            combineSrcDict[key]=srcDict0[key]+srcDict1[key]
        else:
            combineSrcDict[key]=np.mean([srcDict0[key],srcDict1[key]])
    return combineSrcDict

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options]')
    o.set_description(__doc__)
    o.add_option('-o','--out',dest='outfn',default=None,
        help='Output catalogue filename')
    o.add_option('--fmt',dest='fmt',default='pickle',
        help='Output format: pickle, ascii default:pickle')
    o.add_option('-r','--radius',dest='radius',default=20.,type='float',
        help='Clustering radius, in arcmin default: 20')
    opts, args = o.parse_args(sys.argv[1:])

    #get header labels from first file, assumes all files have the same format
    fh=open(args[0],'r')
    for lines in fh.readlines():
        if lines.startswith('# Gaus_id'): #Start of header line
            hdrItems=lines.split('\n')[0].split(' ')
            hdrItems=hdrItems[1:]
    fh.close()
    #New header items added: MJD, Alt, Az
    hdrItems.extend(['MJD', 'Alt', 'Az'])
    for hid,hdrVal in enumerate(hdrItems):
        if hdrVal.startswith('RA'): raIdx=hid
        if hdrVal.startswith('DEC'): decIdx=hid

    #HARDCODE: assumes a PAPER-32 SA observer
    obs=ephem.Observer()
    obs.lat='-30:43:17.5'
    obs.long='21:25:41.9'

    mjdIdx=len(hdrItems)-3 #There might be a prettier way to get this number...

    nfiles=len(args)
    for fid,fn in enumerate(args):
        print 'Loading %s (%i of %i)'%(fn,fid+1,nfiles)
        gaulArr=np.loadtxt(fn, dtype='str')
        mjd=float(fn.split('/')[-1].split('_')[0]) #HARDCODE: this assumes a file format [MJD]_[text]
        obs.date=ephem.date(mjd - 2415020.) #Convert Julian date to ephem date, measured from noon, Dec. 31, 1899.
        #update the source list to include mjd, alt, az
        #updateGaulArr=np.chararray((gaulArr.shape[0],gaulArr.shape[1]+3))
        updateGaulArr=np.zeros((gaulArr.shape[0],gaulArr.shape[1]+3),dtype='S13')
        for srcIdx in range(gaulArr.shape[0]):
            #make a source to compute the Alt,Az
            ra=float(gaulArr[srcIdx][raIdx])*np.pi/180.
            dec=float(gaulArr[srcIdx][decIdx])*np.pi/180.
            src=ephem.FixedBody()
            src._ra=ra
            src._dec=dec
            src.compute(obs)
            updateGaulArr[srcIdx,:-3]=gaulArr[srcIdx]
            #insert MJD, Alt,Az into array
            updateGaulArr[srcIdx,-3:]=[str(mjd),str(float(src.alt)),str(float(src.az))] #a bit of format jiu jistu
        if fid==0: paperSrcs=updateGaulArr
        else: paperSrcs=np.concatenate((paperSrcs,updateGaulArr))

    clusterDict={}
    clusterIdx=0
    for sid in range(paperSrcs.shape[0]):
        #Look for the source in all the clusters created so far
        srcInCluster=False
        for ckey,cval in clusterDict.iteritems(): #for every cluster
            if sid in cval['idx']: #see if the source is in the cluster
                srcInCluster=True
                continue #exit for loop
        if srcInCluster: continue #exit for loop
        else:
            #create a new cluster
            clusterDict[clusterIdx]={'idx':[sid]}
            srciDict={}
            for lbl,attr in zip(hdrItems,paperSrcs[sid]):
                try:
                    srciDict[lbl]=float(attr)
                except ValueError:
                    srciDict[lbl]=attr
            clusterDict[clusterIdx]['src%i'%sid]=srciDict

            #find all other sources that belong in that cluster and add them to it
            for jidx in range(sid+1,paperSrcs.shape[0]):
                dist=np.sqrt(((float(paperSrcs[sid][raIdx])-float(paperSrcs[jidx][raIdx]))**2.)+((float(paperSrcs[sid][decIdx])-float(paperSrcs[jidx][decIdx]))**2.))*60.
                if (dist < opts.radius):
                    srcjDict={}
                    for lbl,attr in zip(hdrItems,paperSrcs[jidx]):
                        try:
                            srcjDict[lbl]=float(attr)
                        except ValueError:
                            srcjDict[lbl]=attr
                    if srciDict['MJD']==srcjDict['MJD']: #if sources have the same MJD, then they should be treated as a single source
                        print 'Combining sources into one for same Cluster, same MJD (ids: %i, %i)'%(sid,jidx)
                        bb=combineSrcs(clusterDict[clusterIdx]['src%i'%sid],srcjDict)
                        clusterDict[clusterIdx]['src%i'%sid]=combineSrcs(clusterDict[clusterIdx]['src%i'%sid],srcjDict)
                    else:
                        clusterDict[clusterIdx]['idx'].append(jidx)
                        clusterDict[clusterIdx]['src%i'%jidx]=srcjDict

            clusterIdx+=1
    print 'Found %i clusters for %i sources'%(clusterIdx,paperSrcs.shape[0])

    #compute the average RA,DEC
    for ckey,cval in clusterDict.iteritems():
        #print ckey
        #print clusterDict[ckey].keys()
        raList=[]
        decList=[]
        for src in clusterDict[ckey].iterkeys():
            if src.startswith('src'):
                raList.append(clusterDict[ckey][src]['RA'])
                decList.append(clusterDict[ckey][src]['DEC'])
        clusterDict[ckey]['ara']=np.mean(np.array(raList)) #average RA
        clusterDict[ckey]['adec']=np.mean(np.array(decList)) #average Dec

    #write to output file
    if opts.outfn is None: outfn='clusterCat.pkl'
    else: outfn=opts.outfn
    print 'Writing catalogue dictionary to %s'%outfn
    pkl.dump(clusterDict,open(outfn,'wb'))

