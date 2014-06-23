from __future__ import division
import sys
import numpy as np
import matplotlib
if "show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab

datafile = 'data/shape-p/box-plot--p-0.30-0.30-0.00-0.00-15.00'

def ignoreme(value):
    return 0.0

def readbox(name):
    data = np.loadtxt(name, converters = {0: ignoreme, 1: ignoreme})
    data = data[:,2:]
    good = np.zeros((7, len(data[0,:])))
    nsections = len(data[:,0])//7
    for i in range(7):
        # add up over all sections!
        good[i,:] = np.sum(data[i*nsections:(i+1)*nsections,:], axis=0)
    return good[:,0:]

boxfull = readbox(datafile + '-full_array.dat')
boxexact = readbox(datafile + '-exact.dat')

proteins = ['nATP', 'nADP', 'nE', 'ND', 'NDE', 'NflD', 'NflE']


def print_analysis(which):
    pylab.figure()
    pylab.title(proteins[which])
    pylab.plot(boxfull[which,:])
    pylab.plot(boxexact[which,:])

    print 'working on', proteins[which]
    print '=================='
    #print 'full mean', boxfull[which,:].mean()
    #print 'full stdev', boxfull[which,:].std()
    num_uncorrelated = len(boxfull[which,:])/200.0
    error = boxfull[which,:].std()/np.sqrt(num_uncorrelated)
    #print 'full error assuming uncorrelated', error
    #print 'exact mean', boxexact[which,:].mean()
    #print ''
    print 'num_uncorrelated', num_uncorrelated, 'fractional uncertainty', error/boxfull[which,:].mean()
    print 'off by', (boxfull[which,:].mean() - boxexact[which,:].mean())/error, 'standard errors'

for i in range(5):
    print_analysis(i)

pylab.show()
