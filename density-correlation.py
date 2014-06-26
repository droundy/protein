from __future__ import division
import sys
import numpy as np
import matplotlib
if "show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab

datafile = 'data/shape-%s/box-plot--%s-%s-%s-%s-%s-%s'% \
    (sys.argv[1],sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
datafile = 'data/shape-p/box-plot--p-3.00-0.50-0.00-0.00-15.00'

start = float(sys.argv[7])
end  = float(sys.argv[8])

def ignoreme(value):
    return 0.0

def readbox(name):
    data = np.loadtxt(name, converters = {0: ignoreme, 1: ignoreme})
    data = data[:,2:]
    summed_data = np.zeros((7, len(data[0,:])))
    nsections = len(data[:,0])//7
    for i in range(7):
        # add up over all sections!
        summed_data[i,:] = np.sum(data[i*nsections:(i+1)*nsections,:], axis=0)
    shortened_data = np.copy(summed_data[:,int(len(data[0,:])*start/10.0):int(len(data[0,:])*end/10.0)])
    return shortened_data[:,:]

boxfull = readbox(datafile + '-full_array.dat')
boxexact = readbox(datafile + '-exact.dat')

#proteins = ['nATP', 'nADP', 'nE', 'ND', 'NDE', 'NflD', 'NflE']


def mean_density(data,protein):
    total = 0
    for i in data[protein,:]:
        total += i
    return total/len(data[protein,:])


dt = .1
def correlation(data,protein,time):
    integral = 0
    total_time = dt*float(len(data[protein,:])-int(time/dt))
    mean = mean_density(data,protein)
    for i in range(len(data[protein,:])-int(time/dt)):
        integral += (data[protein,i] - mean) \
            * (data[protein,i+int(time/dt)] - mean) * dt
    return integral/(mean*mean*total_time)

difD = 2.5
dx = .05
time_step = .1*dx*dx/difD;#sec
print_denominator = 1000;
#  if (i%print_denominator==0)
dt = time_step*print_denominator

time_array =  np.arange(0,.5*len(boxfull[0,:])*dt,dt)
correlation_array =  np.zeros(len(time_array))
for i in range(len(time_array)):
    correlation_array[i] = correlation(boxfull,3,time_array[i])

#Running on the first half of the data to compare between the two:
boxfull_short = np.copy(boxfull[:,:len(boxfull[0,:])/2])

time_array_short =  np.arange(0,.5*len(boxfull_short[0,:])*dt,dt)
correlation_array_short =  np.zeros(len(time_array_short))
for i in range(len(time_array_short)):
    correlation_array_short[i] = correlation(boxfull_short,3,time_array_short[i])




pylab.figure()
pylab.subplot(311)
pylab.title('Full data')
pylab.plot(time_array,correlation_array)
#pylab.xlabel('time (sec)')
#pylab.ylabel('Auto Correlation')

pylab.subplots_adjust(hspace=0.4)

pylab.subplot(312)
pylab.title('Shortened data')
pylab.plot(time_array_short,correlation_array_short)
#pylab.xlabel('time (sec)')
pylab.ylabel('Auto Correlation')


time_array =  np.arange(0,.5*len(boxexact[0,:])*dt,dt)
correlation_array =  np.zeros(len(time_array))
for i in range(len(time_array)):
    correlation_array[i] = correlation(boxexact,3,time_array[i])

pylab.subplots_adjust(hspace=0.4)
#pylab.figure()
pylab.subplot(313)
pylab.title('Continuous')
pylab.plot(time_array,correlation_array)
pylab.xlabel('time (sec)')
#pylab.ylabel('Auto Correlation')
printout_file = 'data/shape-%s/plots/correlation-box-%s-%s-%s-%s-%s-%s.pdf'% \
    (sys.argv[1],sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
pylab.savefig(printout_file)

