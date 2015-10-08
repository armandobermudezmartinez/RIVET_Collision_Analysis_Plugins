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
    out[k] = src[k]
    if out[k].integral() == 0: continue
    out[k].normalize()

yoda.writeYODA(out, sys.argv[2])
