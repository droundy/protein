from __future__ import division
import sys, os
import os.path
import numpy as np
import matplotlib
if "show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab
import re
import math as m

start = float(sys.argv[8])
end  = float(sys.argv[9])

difD = 2.5
dx = .05
time_step = .1*dx*dx/difD;#sec
print_denominator = 1000;
dt = time_step*print_denominator

# def notify_reading_file(f):
#     #print 'reading', f
#     #assert os.system('git add -f %s' % f) == 0, "Couldn't git add %s!" % f
#     pass


job_string = "data/shape-%s/%s-%s-%s-%s-%s-full_array/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                   sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
data_file = job_string + 'fast-correlation-right-left.dat'
print 'loading file ',data_file
data_full = np.loadtxt(data_file)

job_string = "data/shape-%s/%s-%s-%s-%s-%s-exact/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                   sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
data_file = job_string + 'fast-correlation-right-left.dat'
print 'loading file ',data_file
data_exact = np.loadtxt(data_file)

# plt.figure()
# plt.plot(data_exact[0],data_exact[1])
# plt.show()
# exit(0)

print 'start time = ',start,', end_time = ',end
if (len(data_full[0,:])*dt < end or len(data_exact[0,:])*dt < end):
    print 'we dont have enough data for this end time'
    print 'len(data[0,:])*dt = ',min(len(data_full[0,:])*dt,len(data_exact[0,:])*dt)
    exit(1)

data_full = data_full[:,int(start/dt):int(end/dt)]
data_exact = data_exact[:,int(start/dt):int(end/dt)]
time_array = data_full[0]

#######################################
#Following is for the decay modelling:
max_Cs = []
max_times = []
num_calc_periods = int(sys.argv[10])

for i in range(3,len(time_array[3:-3])):
    if (data_full[1][i] >= data_full[1][i-1]) and (data_full[1][i] > data_full[1][i+1]):
        if (data_full[1][i] >= data_full[1][i-2]) and (data_full[1][i] >= data_full[1][i+2]):
            if (data_full[1][i] > data_full[1][i-3]) and (data_full[1][i] > data_full[1][i+3]):
                max_Cs.append(data_full[1][i])
                max_times.append(time_array[i])
                if len(max_Cs) > num_calc_periods:
                    break

if len(max_Cs) <= num_calc_periods:
    print '\nYou need more data there havent been %d maxima in the correlatation yet\n'%(num_calc_periods)
    exit(0)

period = (max_times[num_calc_periods] - max_times[0]) / num_calc_periods
w = 2*m.pi/period

print 'max_times ',max_times
print 'max_Cs ',max_Cs
rate = m.log(max_Cs[0]/max_Cs[num_calc_periods])/(max_times[num_calc_periods]-max_times[0])

print 'Using the equation C = -cos(wt)*exp(rate*t), the period is = ',period,'and the rate is = ',rate

decay = np.zeros_like(time_array) #same size as short correlation array
Amplitude = max_Cs[0]/m.exp(-rate*max_times[0])#fix to the first peak
for i in range(len(decay)):
    decay[i] = -Amplitude*m.cos(w*time_array[i])*m.exp(-rate*time_array[i])

####################################

normalization = np.amax(np.absolute(data_exact[1]))
data_full[1] = data_full[1]/normalization
data_exact[1] = data_exact[1]/normalization
decay = decay/normalization

# data_full[2] = data_full[2]/np.amax(np.absolute(data_full[2]))
# data_exact[1] = data_exact[1]/np.amax(np.absolute(data_exact[1]))
# data_exact[2] = data_exact[2]/np.amax(np.absolute(data_exact[2]))

plot_final_version = True

colors = ['g:','b-','c--']
labels = ['','Long Array','Short Array']
size_scale = 1.6
pylab.figure(figsize=(4*size_scale,3*size_scale))

num_plots = len(data_full)
if plot_final_version:
    num_plots = 2

for i in range(1,num_plots):
    print 'period ',period
    pylab.plot(time_array,data_full[i],colors[i],label='stochastic')
    pylab.xlim(0,time_array[-1])
    pylab.ylim(np.amin(data_full[i]),np.amax(data_full[i]))

if not plot_final_version:
    pylab.plot(time_array,decay,colors[0],label='decay')

pylab.title(r'%s %s $\tau_c = %.3g$ sec $T = %.3g$ sec (ratio %.2g)'
            % (sys.argv[1], sys.argv[5], 1/rate, period, 1/(rate*period)))
if sys.argv[1] == 'randst' and sys.argv[5] == '95.00':
    #pylab.title(r'mannik A, $\tau_c = %.3g$ sec $T = %.3g$ sec (ratio %.2g)'
    #            % (1/rate, period, 1/(rate*period)))
    pylab.title(r'Temporal correlation of mannik A')
elif sys.argv[1] == 'stad' and sys.argv[4] == '1.32':
    #pylab.title(r'stadium A, $\tau_c = %.3g$ sec $T = %.3g$ sec (ratio %.2g)'
    #            % (1/rate, period, 1/(rate*period)))
    pylab.title(r'Temporal correlation of stadium A')
elif sys.argv[1] == 'p':
    #pylab.title(r'pill-shaped cell, $\tau_c = %.3g$ sec $T = %.3g$ sec (ratio %.2g)'
    #            % (1/rate, period, 1/(rate*period)))
    pylab.title(r'Temporal correlation of wild-type pill shape')



    #data_full = data_full/np.amax(data_full)
#pylab.subplot(212)
colors = ['k:','r:','m--.']
for i in range(1,num_plots):
    pylab.plot(time_array,data_exact[i],colors[i],label='deterministic')
    pylab.xlim(0,time_array[-1])
    pylab.ylim(np.amin(data_exact[i]),np.amax(data_exact[i]))

pylab.ylabel('$C(t)$')
plt.xlabel(r'$\tau$ (sec)')
plt.legend()

for string in ["full_array","exact"]:
    job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                       sys.argv[5],sys.argv[6],string)
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)
    if not os.path.exists(job_string + 'plots'):
        print "making directory "+job_string + 'plots'+" because doesnt exist"
        os.makedirs(job_string + 'plots')
    printout_file = job_string + 'plots/correlation.pdf'
    print 'printing to ',printout_file
    pylab.savefig(printout_file)

for i in sys.argv:
    print i
if "show" in sys.argv:
    plt.show()
