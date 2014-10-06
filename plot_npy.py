#!/usr/bin/env python
"""
Plot a NPY file
"""

import sys,os
import numpy as np
import pylab as p

if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] NPY_FILE')
    o.set_description(__doc__)
    o.add_option('-x','--xaxis',dest='xaxis',default=0,type='int',
        help='Axis to plot along X axis, default:0')
    o.add_option('-y','--yaxis',dest='yaxis',default=1,type='int',
        help='Axis to plot along Y axis, default:1')
    opts, args = o.parse_args(sys.argv[1:])

    fn=args[0]
    imCube=np.load(fn)
    print 'NPY shape:', imCube.shape
    print 'X Axis: %i\t Y Axis: %i\t'%(opts.xaxis,opts.yaxis)

    #p.imshow()
    #p.show()

