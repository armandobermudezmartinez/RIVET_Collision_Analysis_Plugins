#! /usr/bin/env python

import yoda
from math import sqrt
from collections import OrderedDict

data = yoda.read('CMS_2017_I1608166_hepdata.yoda')

tabs = OrderedDict()
tabs['01_pi+'] = data['/REF/CMS_2017_I1608166/d01-x01-y01']
tabs['02_k+']  = data['/REF/CMS_2017_I1608166/d01-x01-y02']
tabs['03_p']   = data['/REF/CMS_2017_I1608166/d01-x01-y03']
tabs['04_pi-'] = data['/REF/CMS_2017_I1608166/d02-x01-y01']
tabs['05_k-']  = data['/REF/CMS_2017_I1608166/d02-x01-y02']
tabs['06_p~']  = data['/REF/CMS_2017_I1608166/d02-x01-y03']

def fixXErrs(scatter):
    for p in scatter.points():
        if p.xErrAvg() == 0.:
            p.setXErrs(0.025)

for t in data:
    fixXErrs(data[t])


def addScatters(scatter1, scatter2, path = '', title = ''):
    scatter_out = yoda.Scatter2D(path, title)
    for p in scatter1.points():
        for q in scatter2.points():
            if p.x() == q.x():
                pqy_err2 = 0.
                for unc in ['stat', 'sys']:
                    pqy_err2 += p.yErrsFromSource(unc)[0]**2
                    pqy_err2 += q.yErrsFromSource(unc)[0]**2
                scatter_out.addPoint(p.x(), p.y() + q.y(), p.xErrs()[0], sqrt(pqy_err2))
                break
    return scatter_out

def divideScatters(scatter1, scatter2, path = '', title = ''):
    scatter_out = yoda.Scatter2D(path, title)
    for p in scatter1.points():
        for q in scatter2.points():
            if p.x() == q.x():
                pq = p.y() / q.y()
                if b'stat' in scatter1.variations():
                    py_err2 = 0.
                    qy_err2 = 0.
                    for unc in ['stat', 'sys']:
                        py_err2 += p.yErrsFromSource(unc)[0]**2
                        qy_err2 += q.yErrsFromSource(unc)[0]**2
                else:
                    py_err2 = p.yErrAvg()**2
                    qy_err2 = q.yErrAvg()**2
                py_err2 /= p.y()**2
                qy_err2 /= q.y()**2
                
                pq_err = sqrt(py_err2 + qy_err2) * pq
                scatter_out.addPoint(p.x(), pq, p.xErrs()[0], pq_err)
                break
    return scatter_out

tabs['12_pi-/pi+'] = divideScatters(tabs['04_pi-'], tabs['01_pi+'], '/REF/CMS_2017_I1608166/d100-x01-y01', 'pi-/pi+')
tabs['13_k-/k+']   = divideScatters(tabs['05_k-'],  tabs['02_k+'] , '/REF/CMS_2017_I1608166/d100-x01-y02', 'k-/k+')
tabs['14_p~/p']    = divideScatters(tabs['06_p~'],  tabs['03_p']  , '/REF/CMS_2017_I1608166/d100-x01-y03', 'p~/p')

sum_pions   = addScatters(tabs['04_pi-'], tabs['01_pi+'])
sum_kaons   = addScatters(tabs['05_k-'],  tabs['02_k+'])
sum_protons = addScatters(tabs['06_p~'],  tabs['03_p'])

tabs['15_k/pi'] = divideScatters(sum_kaons,   sum_pions, '/REF/CMS_2017_I1608166/d101-x01-y01', 'k/pi')
tabs['16_p/pi'] = divideScatters(sum_protons, sum_pions, '/REF/CMS_2017_I1608166/d101-x01-y02', 'p/pi')


yoda.write(tabs, 'CMS_2017_I1608166.yoda')