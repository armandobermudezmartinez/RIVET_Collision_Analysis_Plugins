#! /usr/bin/env python

import yoda
from collections import OrderedDict
from math import sqrt

tabs = OrderedDict()

count = 0
for obs in ['ph_energy', 'ph_angle']:
    count += 1
    tabs[obs] = yoda.Scatter2D('/REF/DELPHI_1994_I375963/d0%i-x01-y01' % count, obs)
    
    data = []
    with open('DELPHI_1994_I375963_%s.csv' % obs) as f:
        for line in f:
            if line.startswith('#'):
                continue
            dataline = []
            for s in line.split(','):
                if s == '':
                    dataline.append(-999.)
                else:
                    dataline.append(float(s))
            data.append(dataline)
            
    for i in range(len(data)-1):
        x  = (data[i+1][6] + data[i][6]) / 2.
        ex = (data[i+1][6] - data[i][6]) / 2.
        y = data[i][1] / (2*ex)
        eff_err = 0.08 * data[i][1]
        eyplus  = sqrt((data[i][3] - data[i][1])**2 + eff_err**2) / (2*ex)
        eyminus = sqrt((data[i][1] - data[i][5])**2 + eff_err**2) / (2*ex)
        tabs[obs].addPoint(x, y, (ex, ex), (eyminus, eyplus))


yoda.write(tabs, 'DELPHI_1994_I375963.yoda')