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

start = float(sys.argv[8])
end  = float(sys.argv[9])

difD = 2.5
dx = .05
time_step = .1*dx*dx/difD;#sec
print_denominator = 1000;
dt = time_step*print_denominator

job_string_full = "/data/shape-%s/%s-%s-%s-%s-%s-full_array/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                    sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string_full = p.sub('_',job_string_full)

job_string_exact = "/data/shape-%s/%s-%s-%s-%s-%s-exact/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                    sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string_exact = p.sub('_',job_string_exact)

data_file_full = '.' + job_string_full + 'correlation.dat'
data_file_exact = '.' + job_string_exact + 'correlation.dat'
data_full = np.loadtxt(data_file_full)
data_exact = np.loadtxt(data_file_exact)
print 'loading data from ',data_file_full
print 'loading data from ',data_file_exact


print 'start time = ',start,', end_time = ',end
if (len(data_full[0,:])*dt < end or len(data_exact[0,:])*dt < end):
    print 'we dont have enough data for this end time'
    print 'len(data[0,:])*dt = ',min(len(data_full[0,:])*dt,len(data_exact[0,:])*dt)
    exit(1)
data_start_index = int(start/dt)
data_stop_index = int(end/dt)

print 'this is the time covered by the entire data set that youre using'
data_full = data_full[:,data_start_index:data_stop_index]
data_exact = data_exact[:,data_start_index:data_stop_index]

time_array = data_full[0]

correlation_arrays = [0]*(len(data_full)-1)
for i in range(len(correlation_arrays)):
    correlation_arrays[i] = data_full[i+1]

colors = ['m','g','r','y','b','c']
labels = ['Long Array Left','Long Array Mid','Long Array Right','Short Array Left','Short Array Mid','Short Array Right']
pylab.figure()
for i in range(len(correlation_arrays)):
    if (i != 0) and (i != 3):
        pylab.subplot(211)
        pylab.title('Stochastic Data')
        pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
        pylab.xlim(0,time_array[-1])
        pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
        pylab.ylabel('Auto Correlation')

correlation_arrays = [0]*(len(data_exact)-1)
for i in range(len(correlation_arrays)):
    correlation_arrays[i] = data_exact[i+1]

colors = ['m','g','r','y','b','c']
labels = ['Long Array Left','Long Array Mid','Long Array Right','Short Array Left','Short Array Mid','Short Array Right']
for i in range(len(correlation_arrays)):
    if (i != 0) and (i != 3):
        pylab.subplot(212)
        pylab.title('Determinisitc Data')
        pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
        pylab.xlim(0,time_array[-1])
        pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
        pylab.ylabel('Auto Correlation')

plt.xlabel('time (sec)')
plt.legend()

printout_file = '.' + job_string_full + 'plots/correlation.pdf'
print 'printing to ',printout_file
pylab.savefig(printout_file)

printout_file = '.' + job_string_exact + 'plots/correlation.pdf'
print 'printing to ',printout_file
pylab.savefig(printout_file)





#pylab.ylabel('Auto Correlation')
#pylab.subplots_adjust(hspace=0.4)
#pylab.subplot(111)
#pylab.title('Shortened data_full')
#pylab.plot(time_array,correlation_array_short,color='r' )
#pylab.xlim(0,time_array[-1]/2)
#pylab.ylim(np.amin(correlation_array),np.amax(correlation_array))
#pylab.xlabel('time (sec)')
