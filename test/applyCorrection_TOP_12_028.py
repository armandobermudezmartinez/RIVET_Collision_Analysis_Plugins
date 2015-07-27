#!/usr/bin/env python

#from ROOT import *
import yoda
src = yoda.readYODA("out.yoda")
out = {}

corrFactors = {
    'h14_diffXSecTopSemiLepPartontopPt' : [1.161660, 1.136480, 1.020996, 0.895649, 0.772136, 0.685911, 0.559711, 0.566430],
    'h16_diffXSecTopSemiLepPartontopY' : [2.101211, 1.099831, 0.937698, 0.883005, 0.868135, 0.882153, 0.878180, 0.941096, 1.095958, 2.056497],
    'h20_diffXSecTopSemiLepPartonttbarPt' : [1.602612, 0.913407, 0.816876, 0.849766, 0.889415, 0.857082],
    'h21_diffXSecTopSemiLepPartonttbarY' : [2.461665, 1.147150, 0.908031, 0.848166, 0.814687, 0.803214, 0.824948, 0.947269, 1.122359, 2.428979],
    'h22_diffXSecTopSemiLepPartonttbarMass' : [1.498358, 1.362128, 1.024490, 0.819021, 0.646227, 0.475925, 0.372441],

    'h23_diffXSecTopDiLepPartontopPt' : [0.933825, 1.069645, 1.051336, 0.919932, 0.774565],
    'h25_diffXSecTopDiLepPartontopY' : [1.682022, 1.002849, 0.925246, 0.924734, 0.880097, 0.901330, 1.042041, 1.733911],
    'h29_diffXSecTopDiLepPartonttbarPt' : [1.129278, 0.908123, 0.933110, 0.963850],
    'h30_diffXSecTopDiLepPartonttbarY' : [2.401265, 1.140515, 0.937143, 0.889803, 0.833903, 0.946386, 1.179555, 2.445021],
    'h31_diffXSecTopDiLepPartonttbarMass' : [0.803342, 1.136017, 1.206834, 1.037619, 1.081579, 0.741247],
}

for k in src.viewkeys():

    hName = k.split('/')[-1]
    if hName not in corrFactors:
        out[k] = src[k]
    else:
        f = corrFactors[hName]
        #out[k] = src[k].clone()
        out[k] = yoda.Scatter2D(src[k].path, src[k].title)
        #out[k].reset()
        for i, b in enumerate(src[k].bins):
            x = (b.xEdges[0]+b.xEdges[1])/2
            xErr = (b.xEdges[1]-b.xEdges[0])/2
            #out[k].addPoint(x, b.area*f[i], xErr, b.areaErr*f[i])
            out[k].addPoint(x, b.height*f[i], xErr, b.heightErr*f[i])
            #out[k].fillBin(i, f[i]*b.sumW)

yoda.writeYODA(out, "corr.yoda")
