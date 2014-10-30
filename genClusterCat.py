#!/usr/bin/env python
"""Generate a clustered catalogue from pyBDSM Gaul files"""

import numpy as np
import cPickle as pkl
import os,sys

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
    for hid,hdrVal in enumerate(hdrItems):
        if hdrVal.startswith('RA'): raIdx=hid
        if hdrVal.startswith('DEC'): decIdx=hid

    nfiles=len(args)
    for fid,fn in enumerate(args):
        print 'Loading %s (%i of %i)'%(fn,fid+1,nfiles)
        gaulArr=np.loadtxt(fn, dtype='str')
        if fid==0: paperSrcs=gaulArr
        else: paperSrcs=np.concatenate((paperSrcs,gaulArr))

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
            srcDict={}
            for lbl,attr in zip(hdrItems,paperSrcs[sid]):
                try:
                    srcDict[lbl]=float(attr)
                except ValueError:
                    srcDict[lbl]=attr
            clusterDict[clusterIdx]['src%i'%sid]=srcDict

            #find all other sources that belong in that cluster and add them to it
            for jidx in range(sid+1,paperSrcs.shape[0]):
                dist=np.sqrt(((float(paperSrcs[sid][raIdx])-float(paperSrcs[jidx][raIdx]))**2.)+((float(paperSrcs[sid][decIdx])-float(paperSrcs[jidx][decIdx]))**2.))*60.
                if (dist < opts.radius):
                    clusterDict[clusterIdx]['idx'].append(jidx)
                    for lbl,attr in zip(hdrItems,paperSrcs[jidx]):
                        try:
                            srcDict[lbl]=float(attr)
                        except ValueError:
                            srcDict[lbl]=attr
                    clusterDict[clusterIdx]['src%i'%jidx]=srcDict

            clusterIdx+=1
    print 'Found %i clusters for %i sources'%(clusterIdx,paperSrcs.shape[0])

    #compute the average RA,DEC
    for ckey,cval in clusterDict.iteritems():
        print ckey
        print clusterDict[ckey].keys()
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

