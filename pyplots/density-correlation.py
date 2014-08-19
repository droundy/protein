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

def ignoreme(value):
    return 0.0

def readbox(name):
    data = np.loadtxt(name, converters = {0: ignoreme, 1: ignoreme})
    data = data[:,2:]
    #data = data[:,:1000]
    summed_data = np.zeros((7, len(data[0,:])))
    nsections = len(data[:,0])//7
    for i in range(7):
        # add up over all sections!
        summed_data[i,:] = np.sum(data[i*nsections:(i+1)*nsections,:], axis=0)
    shortened_data = np.copy(summed_data[:,int(len(data[0,:])*start/10.0):int(len(data[0,:])*end/10.0)])
    return shortened_data[:,:]

plot_exact = False
job_string = "/data/shape-%s/%s-%s-%s-%s-%s-exact/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                       sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
print '.' + job_string + 'box-plot.dat'
if os.path.isfile('.' + job_string + 'box-plot.dat'):
    boxexact = readbox('.' + job_string + 'box-plot.dat')
    plot_exact = True
else:
    print '.' + job_string + 'box-plot.dat  doesnt exist so we wont print it'

job_string = "/data/shape-%s/%s-%s-%s-%s-%s-full_array/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                            sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
boxfull = readbox('.' + job_string + 'box-plot.dat')


#proteins = ['nATP', 'nADP', 'nE', 'ND', 'NDE', 'NflD', 'NflE']


def mean_density(data,protein):
    total = 0
    for i in data[protein,:]:
        total += i
    return float(total)/float(len(data[protein,:]))


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

printout_file = '.' + job_string + 'correlation.dat'

print int(len(correlation_array)/2)
print len(correlation_array_short)
print len(time_array)
print len(time_array_short)

p_file = open(printout_file,'w')
for t in time_array_short:
    t = float(t)
    p_file.write('%g '%t)
p_file.write('\n')
correlation_array = correlation_array[:len(correlation_array_short)]
for c in correlation_array:
    p_file.write('%g '%c)
p_file.write('\n')
for c in correlation_array_short:
    p_file.write('%g '%c)
p_file.close()


