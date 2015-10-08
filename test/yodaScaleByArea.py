#!/usr/bin/env python

import sys, os
if len(sys.argv) < 3:
    print "Usage: %s INPUT.yoda OUTPUT.yoda" % sys.argv[0]
    print "\nScale all histograms by unit area"
    sys.exit()

from math import *
import yoda
src = yoda.readYODA(sys.argv[1])

out = {}
for k in src.viewkeys():
    h = src[k]
    out[k] = yoda.Scatter2D(h.path, h.title)
    area = h.integral(True)
    for b in h.bins:
        x = (b.xEdges[0]+b.xEdges[1])/2
        xErr = (b.xEdges[1]-b.xEdges[0])/2
        y, yErr = b.height, b.heightErr
        out[k].addPoint(x, y/area, xErr, yErr/area)

yoda.writeYODA(out, sys.argv[2])
