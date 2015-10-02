#!/usr/bin/env python

from math import *
import yoda
src = yoda.readYODA("MC.yoda")

hists = {}

## Get histogram pairs
for k in src.viewkeys():
    module, hName = k.split('/')[-2:]
    if 'CMS_2015_I1370682' not in module: continue

    if hName not in hists: hists[hName] = []

    hists[hName].append(src[k])

## Get correction factors
corrFactors = {}
for k in hists.keys():
    if len(hists[k]) != 2: continue
    h1, h2 = hists[k]
    # put particle level front, parton level latter.
    if 'internal' in h1.path: h2, h1 = h1, h2

    cfs = []
    for i, b in enumerate(h1.bins):
        #x = (b.xEdges[0]+b.xEdges[1])/2
        #xErr = (b.xEdges[1]-b.xEdges[0])/2
        y1, y1Err = h1.bins[i].height, h1.bins[i].heightErr
        y2, y2Err = h2.bins[i].height, h2.bins[i].heightErr

        r, rErr = 1, 1000
        if y1 > 0: r = y2/y1 ## We are calculating particle -> parton in full PS
        if y1 > 0 and y2 > 0: rErr = r*sqrt((y1Err/y1)**2+(y2Err/y2)**2)

        cfs.append([r, rErr])
    corrFactors[h1.path] = cfs[:]

## Print out correction functions in C++ code format
for path in corrFactors.keys():
    hName = path.split('/')[-1]

    print '//', hName
    print '{',
    for x in corrFactors[path]: print ("%f," % x[0]),
    print '};'
