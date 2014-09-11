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

if (sys.argv[10] != "auto") and (sys.argv[10] != "rl"):
    print 'The last command line argument, after the time, should be either \"auto\" or \"rl\"'
    print 'for auto-correlation or right-left correlation'
    exit(1)
corr_type = sys.argv[10]

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

data_file_full = ''
data_file_exact = ''
if corr_type == "auto":
    data_file_full = '.' + job_string_full + 'auto-correlation.dat'
    data_file_exact = '.' + job_string_exact + 'auto-correlation.dat'
elif corr_type == "rl":
    data_file_full = '.' + job_string_full + 'correlation-right-left.dat'
    data_file_exact = '.' + job_string_exact + 'correlation-right-left.dat'

data_full = np.loadtxt(data_file_full)
data_exact = np.loadtxt(data_file_exact)

print 'loading data from ',data_file_full
print 'loading data from ',data_file_exact

print data_exact.shape
print data_full.shape
print data_exact
print data_full


print 'start time = ',start,', end_time = ',end
if (len(data_full[0,:])*dt < end or len(data_exact[0,:])*dt < end):
    print 'we dont have enough data for this end time'
    print 'len(data[0,:])*dt = ',min(len(data_full[0,:])*dt,len(data_exact[0,:])*dt)
    exit(1)
data_start_index = int(start/dt)
data_stop_index = int(end/dt)

data_full = data_full[:,data_start_index:data_stop_index]
data_exact = data_exact[:,data_start_index:data_stop_index]

time_array = data_full[0]

correlation_arrays = [0]*(len(data_full)-1)
for i in range(len(correlation_arrays)):
    correlation_arrays[i] = data_full[i+1]

correlation_arrays = correlation_arrays/np.amax(correlation_arrays)

colors = ['m','g','r','y','b','c']
labels = ['Long Array Left','Long Array Mid','Long Array Right','Short Array Left','Short Array Mid','Short Array Right']

if (corr_type == "auto"):
    pylab.figure()
    for i in range(len(correlation_arrays)):
        if (i == 2) or (i == 5):
            pylab.subplot(211)
            pylab.title('Stochastic Data, Auto-correlation, '+sys.argv[1]+' '+sys.argv[5])
            pylab.title(sys.argv[1])
            pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
            pylab.xlim(0,time_array[-1])
            pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
            pylab.ylabel('Auto Correlation')

    correlation_arrays = [0]*(len(data_exact)-1)
    for i in range(len(correlation_arrays)):
        correlation_arrays[i] = data_exact[i+1]

    correlation_arrays = correlation_arrays/np.amax(correlation_arrays)
    for i in range(len(correlation_arrays)):
        if (i == 2) or (i == 5):
            pylab.subplot(212)
            pylab.title('Determinisitc Data')
            pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
            pylab.xlim(0,time_array[-1])
            pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
            pylab.ylabel('Auto Correlation')

    plt.xlabel('time (sec)')
    plt.legend()

colors = ['b','g']
labels = ['Long Array','Short Array']
if (corr_type == "rl"):
    pylab.figure()
    for i in range(len(correlation_arrays)):
        pylab.subplot(211)
        pylab.title('Stochastic Data, Right-left-correlation, '+sys.argv[1]+' '+sys.argv[5])
        pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
        pylab.xlim(0,time_array[-1])
        pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
        pylab.ylabel('Left-Right Correlation')

    correlation_arrays = [0]*(len(data_exact)-1)
    for i in range(len(correlation_arrays)):
        correlation_arrays[i] = data_exact[i+1]

    correlation_arrays = correlation_arrays/np.amax(correlation_arrays)
    for i in range(len(correlation_arrays)):
        pylab.subplot(212)
        pylab.title('Deterministic Data')
        pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
        pylab.xlim(0,time_array[-1])
        pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
        pylab.ylabel('Left-Right Correlation')

    plt.xlabel('time (sec)')
    plt.legend()


printout_file = '.' + job_string_full + 'plots/correlation-' + corr_type + '.pdf'
print 'printing to ',printout_file
pylab.savefig(printout_file)

printout_file = '.' + job_string_exact + 'plots/correlation-' + corr_type + '.pdf'
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
