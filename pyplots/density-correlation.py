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

end  = float(sys.argv[9])

def ignoreme(value):
    return 0.0

protein_name = "NflD"

difD = 2.5
dx = .05
time_step = .1*dx*dx/difD;#sec
print_denominator = 1000;
dt = time_step*print_denominator


#following version of readbox selects out NflD
def readbox(name):
    data = np.loadtxt(name, converters = {0: ignoreme, 1: ignoreme})
    data_stop_index = int(end/dt)
    print 'start time = 0, end_time = ',end
    print 'this is the time covered by the entire data set that youre using'
    data = data[:,2:data_stop_index+2]
    shortened_data = np.zeros((3, len(data[0,:])))
    first_column_data = np.genfromtxt(name, dtype='str',usecols=(0,))
    row_num = 0
    while first_column_data[row_num] != protein_name:
        row_num += 1
    shortened_data = np.copy(data[row_num:row_num+3,:])
    return shortened_data[:,:]


job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                       sys.argv[5],sys.argv[6],sys.argv[7])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
print '.' + job_string + 'box-plot.dat'
if not os.path.isfile('.' + job_string + 'box-plot.dat'):
    print '.' + job_string + 'box-plot.dat  doesnt exist so we wont print it'

boxfull = readbox('.' + job_string + 'box-plot.dat')

def mean_density(data,sec):
    total = 0
    for i in data[sec,:]:
        total += i
    return float(total)/float(len(data[sec,:]))

def correlation(data,sec,time):
    integral = 0
    total_time = dt*float(len(data[sec,:])-int(time/dt))
    mean = mean_density(data,sec)
    for i in range(len(data[sec,:])-int(time/dt)):
        integral += (data[sec,i] - mean) \
            * (data[sec,i+int(time/dt)] - mean) * dt
    return integral/(mean*mean*total_time)


time_array =  np.arange(0,.5*len(boxfull[0,:])*dt,dt)
correlation_array =  np.zeros((3,len(time_array)))
for i in range(len(time_array)):
    if (i%50 == 0):
        print 'time = ',i*dt,' we are ',int(100.0*2.0*float(i)/float(len(boxfull[0,:]))),' percent done with long array'
    for sec in range(len(correlation_array[:,0])):
            correlation_array[sec,i] = correlation(boxfull,sec,time_array[i])
    #print correlation_array[i]

#Running on the first half of the data to compare between the two:
boxfull_short = np.copy(boxfull[:,:len(boxfull[0,:])/2])

time_array_short =  np.arange(0,.5*len(boxfull_short[0,:])*dt,dt)
correlation_array_short =  np.zeros((3,len(time_array_short)))
for i in range(len(time_array_short)):
    if i%50 == 0:
        print 'time = ',i*dt,' we are ',int(100.0*4.0*float(i)/float(len(boxfull[0,:]))),' percent done with short array'
    for sec in range(len(correlation_array_short[:,0])):
            correlation_array_short[sec,i] = correlation(boxfull_short,sec,time_array_short[i])

printout_file = '.' + job_string + 'correlation.dat'

p_file = open(printout_file,'w')
for t in time_array_short:
    t = float(t)
    p_file.write('%g '%t)
p_file.write('\n')
correlation_array = correlation_array[:len(correlation_array_short)]
for i in range(len(correlation_array[:,0])):
    for c in correlation_array[i,:]:
        p_file.write('%g '%c)
    p_file.write('\n')
for i in range(len(correlation_array[:,0])):
    for c in correlation_array_short[i,:]:
        p_file.write('%g '%c)
    p_file.write('\n')
p_file.close()
