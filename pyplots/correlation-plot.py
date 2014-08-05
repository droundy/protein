from __future__ import division
import sys
import os.path
import numpy as np
import matplotlib
if "show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab

data_file = './data/shape-%s/plots/corr/correlation-%s-%s-%s-%s-%s-%s.dat'% \
    (sys.argv[1],sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])

data = np.loadtxt(data_file)

time_array = data[0]
correlation_array = data[1]
correlation_array_short = data[2]


pylab.figure()
pylab.subplot(111)
pylab.title('Full data')
pylab.plot(time_array,correlation_array,color='g')
pylab.xlim(0,time_array[-1]/2)
pylab.ylim(np.amin(correlation_array),np.amax(correlation_array)/8)
#pylab.xlabel('time (sec)')
#pylab.ylabel('Auto Correlation')

pylab.subplots_adjust(hspace=0.4)

pylab.subplot(111)
pylab.title('Shortened data')
pylab.plot(time_array,correlation_array_short,color='r' )
#pylab.xlim(0,time_array[-1]/2)
#pylab.ylim(np.amin(correlation_array),np.amax(correlation_array))
#pylab.xlabel('time (sec)')
pylab.ylabel('Auto Correlation')


# if plot_exact:
#     time_array =  np.arange(0,.5*len(boxexact[0,:])*dt,dt)
#     correlation_array =  np.zeros(len(time_array))
#     for i in range(len(time_array)):
#         correlation_array[i] = correlation(boxexact,3,time_array[i])
#     pylab.subplots_adjust(hspace=0.4)
#     #pylab.figure()
#     pylab.subplot(212)
#     pylab.title('Continuous')
#     pylab.plot(time_array,correlation_array)
#     pylab.xlim(time_array[-1])
#     pylab.xlabel('time (sec)')
#     #pylab.ylabel('Auto Correlation')

printout_file = 'data/shape-%s/plots/correlation-box-%s-%s-%s-%s-%s-%s.pdf'% \
    (sys.argv[1],sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
pylab.savefig(printout_file)
