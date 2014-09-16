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
import math as m

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

#########################################

num_of_data_files = 2
def create_data_array(sim_type):
    data_array = []
    init_weight_array = []
    weight_array = []
    for i in range(num_of_data_files):
        job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                                sys.argv[5],sys.argv[6],sim_type)
        p = re.compile('[.]')
        job_string = p.sub('_',job_string)
        data_file = ''
        if i == num_of_data_files-1:
            if corr_type == "auto":
                data_file = job_string + 'fast-auto-correlation.dat'
            elif corr_type == "rl":
                data_file = job_string + 'fast-correlation-right-left.dat'
        else:
            if corr_type == "auto":
                data_file = '../new-protein-%d/'%(i) + job_string + 'fast-auto-correlation.dat'
            elif corr_type == "rl":
                data_file = '../new-protein-%d/'%(i) + job_string + 'fast-correlation-right-left.dat'
        new_data = np.loadtxt(data_file,skiprows=1)
        time_of_total_data = np.loadtxt(data_file)[0][1] - np.loadtxt(data_file)[0][0]
        init_weight_array.append(-new_data[0] + time_of_total_data) # the new_data[0] is an array of the times
        data_array.append(new_data)

    print data_file

    shortest_data_len = 10000000
    for d in data_array:
        if d.shape[1] < shortest_data_len:
            shortest_data_len = d.shape[1]
    print 'shortest data length ',shortest_data_len

    weight_array = []
    for i in range(num_of_data_files):
        weight_array.append(np.zeros(shortest_data_len))

    for i in range(len(weight_array[0])):
        weight_total = 0
        for data_set in range(len(weight_array)):
            weight_total += init_weight_array[data_set][i]
        for data_set in range(len(weight_array)):
            weight_array[data_set][i] = init_weight_array[data_set][i]/weight_total

    data = np.zeros_like(data_array[0][:,:shortest_data_len])
    data[0] = np.copy(data_array[0][0,:shortest_data_len])

    for data_set in range(0,len(data_array)):
        for j in range(1,len(data_array[0])):
            for i in range(shortest_data_len):
                data[j][i] += data_array[data_set][j][i]*weight_array[data_set][i]
    return data

##############################################3

data_full = create_data_array("full_array")
data_exact = create_data_array("exact")

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

#######################################

#Following is for the decay modelling:
max_Cs = []
max_times = []
num_calc_periods = 5
for i in range(3,len(time_array[3:-3])):
    if (correlation_arrays[0][i] >= correlation_arrays[0][i-1]) and (correlation_arrays[0][i] > correlation_arrays[0][i+1]):
        if (correlation_arrays[0][i] >= correlation_arrays[0][i-2]) and (correlation_arrays[0][i] >= correlation_arrays[0][i+2]):
            if (correlation_arrays[0][i] > correlation_arrays[0][i-3]) and (correlation_arrays[0][i] > correlation_arrays[0][i+3]):
                max_Cs.append(correlation_arrays[0][i])
                max_times.append(time_array[i])
                if len(max_Cs) > num_calc_periods:
                    break



if len(max_Cs) < num_calc_periods:
    print '\nYou need more data there havent been %d maxima in the correlatation yet\n'%(num_calc_periods)
    exit(0)

period = (max_times[num_calc_periods] - max_times[0]) / num_calc_periods
w = 2*m.pi/period
rate = m.log(max_Cs[0]/max_Cs[num_calc_periods])/(max_times[num_calc_periods]-max_times[0])

print 'Using the equation C = -cos(wt)*exp(rate*t), the period is = ',period,'and the rate is = ',rate

decay = np.zeros_like(time_array) #same size as short correlation array
Amplitude = max_Cs[0]/m.exp(-rate*max_times[0])#fix to the first peak
for i in range(len(decay)):
    decay[i] = -Amplitude*m.cos(w*time_array[i])*m.exp(-rate*time_array[i])

normalization = np.amax(np.absolute(correlation_arrays))
correlation_arrays = correlation_arrays/normalization
decay = decay/normalization
####################################



# colors = ['m','g','r','y','b','c']
# labels = ['Long Array Left','Long Array Mid','Long Array Right','Short Array Left','Short Array Mid','Short Array Right']

# if (corr_type == "auto"):
#     pylab.figure()
#     for i in range(len(correlation_arrays)):
#         if (i == 2) or (i == 5):
#             pylab.subplot(211)
#             pylab.title('Stochastic Data, Auto-correlation, '+sys.argv[1]+' '+sys.argv[5]+' rate = '+str(rate)+' period = '+str(period))
#             pylab.title(sys.argv[1])
#             pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
#             pylab.xlim(0,time_array[-1])
#             pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
#             pylab.ylabel('Auto Correlation')

#     correlation_arrays = [0]*(len(data_exact)-1)
#     for i in range(len(correlation_arrays)):
#         correlation_arrays[i] = data_exact[i+1]

#     #correlation_arrays = correlation_arrays/np.amax(correlation_arrays)
#     for i in range(len(correlation_arrays)):
#         if (i == 2) or (i == 5):
#             pylab.subplot(212)
#             pylab.title('Determinisitc Data')
#             pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
#             pylab.xlim(0,time_array[-1])
#             pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
#             pylab.ylabel('Auto Correlation')

#     plt.xlabel('time (sec)')
#     plt.legend()

colors = ['b','g','r']
labels = ['Long Array','Short Array']
if (corr_type == "rl"):
    pylab.figure()
    for i in range(len(correlation_arrays)):
        pylab.subplot(211)
        pylab.title('Stochastic Data, Right left correlatation, '+sys.argv[1]+' '+sys.argv[5]+' rate = '+str(rate)+' period = '+str(period))
        pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
        pylab.xlim(0,time_array[-1])
        pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
        pylab.ylabel('Left-Right Correlation')
    pylab.plot(time_array,decay,color=colors[2],label=labels[i])

    correlation_arrays = [0]*(len(data_exact)-1)
    for i in range(len(correlation_arrays)):
        correlation_arrays[i] = data_exact[i+1]

    #correlation_arrays = correlation_arrays/np.amax(correlation_arrays)
    for i in range(len(correlation_arrays)):
        pylab.subplot(212)
        pylab.title('Deterministic Data')
        pylab.plot(time_array,correlation_arrays[i],color=colors[i],label=labels[i])
        pylab.xlim(0,time_array[-1])
        pylab.ylim(np.amin(correlation_arrays),np.amax(correlation_arrays))
        pylab.ylabel('Left-Right Correlation')

    plt.xlabel('time (sec)')
    plt.legend()

for string in ["full_array","exact"]:
    job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                            sys.argv[5],sys.argv[6],string)
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)

    printout_file = job_string + 'plots/correlation-' + corr_type + '.pdf'
    print 'printing to ',printout_file
    pylab.savefig(printout_file)

for i in sys.argv:
    print i
if "show" in sys.argv:
    plt.show()




#pylab.ylabel('Auto Correlation')
#pylab.subplots_adjust(hspace=0.4)
#pylab.subplot(111)
#pylab.title('Shortened data_full')
#pylab.plot(time_array,correlation_array_short,color='r' )
#pylab.xlim(0,time_array[-1]/2)
#pylab.ylim(np.amin(correlation_array),np.amax(correlation_array))
#pylab.xlabel('time (sec)')
