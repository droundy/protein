from __future__ import division
import sys
import os.path
import numpy as np
import matplotlib
if "show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab
import re

job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                    sys.argv[5],sys.argv[6],sys.argv[7])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
print job_string

data_file = '.' + job_string + 'correlation.dat'
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


printout_file = '.' + job_string + 'plots/correlation.pdf'
pylab.savefig(printout_file)
