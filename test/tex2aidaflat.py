#!/usr/bin/env python
#TEX input converter for TOP-12-028 paper

BASEDIR='papers/TOP-12-028/arxiv-0/'
INPUT=BASEDIR+'supplemental_material.tex'

files = []
for l in open(INPUT).readlines():
    l = l.strip()
    if not l.startswith('\\input{figures'): continue
    files.append(l[7:-1])

for i, f in enumerate(files):
    name = f.split('/')[-1].split('.')[0]
    lines = open(BASEDIR+f).readlines()
    nbins = len(lines)

    print "# BEGIN HISTOGRAM /REF/CMS_TOP_12_028/h%2d_%s" % (i, name)
    print "AidaPath=/REF/CMS_TOP_12_028/h%2d_%s" % (i, name)
    print "Title=%s" % name
    print "## Area: 1"
    print "## Num bins: %d" % nbins
    print "## xlow   xhigh     val     err"
    for l in lines:
        l = l.strip().replace('$', '')
        item = [x.strip() for x in l[:-2].split('&')]
        xlow, xhigh = [float(x) for x in item[0][1:-1].split(',')]
        #xval = item[1]
        yval = item[-4].replace('\\cdot', '*').replace('{', '').replace('}', '')
        yval = yval.replace('10^', '1E')
        yval = eval(yval)
        stat = float(item[-3])
        syst = float(item[-2])
        total = float(item[-1])

        err = yval*total/100

        print "%4.10e %4.10e %4.10e %4.10e" % (xlow, xhigh, yval, err)
    print "# END HISTOGRAM"
    print

