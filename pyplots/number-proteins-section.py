from __future__ import division
import sys
import os.path
import numpy as np
import matplotlib
#if "show" not in sys.argv:
#    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab
import re


def ignoreme(value):
    return 0.0

protein_name = "NflD"

difD = 2.5
dx = .05
time_step = .1*dx*dx/difD;#sec
print_denominator = 1000;
dt = time_step*print_denominator

constant_start_time = 10

def readbox(name):
    data = np.loadtxt(data_file, converters = {0:ignoreme, 1:ignoreme})
    end = (len(data[0])*dt)-10 #a bit less than full data set by default
    data_stop_index = int(end/dt)
    data_start_index = int(constant_start_time/dt)
    print 'start time = ',constant_start_time,', end_time = ',end
    print 'this is the time covered by the entire data set that youre using'
    data = data[:,(data_start_index+2):(data_stop_index+2)]
    shortened_data = np.zeros((3, len(data[0,:])))
    first_column_data = np.genfromtxt(name, dtype='str',usecols=(0,))
    row_num = 0
    while first_column_data[row_num] != protein_name:
        row_num += 1
    shortened_data = np.copy(data[row_num:row_num+3,:])
    return shortened_data[:,:]






job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                   sys.argv[5],sys.argv[6],sys.argv[7])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
data_file = '../new-protein-2/' + job_string + 'box-plot.dat'
if not os.path.isfile(data_file):
        print '\n',data_file, ' doesnt exist so we cant load it'
        exit(0)
print data_file
data = readbox(data_file)
print 'time if data ',dt*len(data[0,:])


plt.figure()
plt.plot(np.arange(len(data[0,:])),data[0,:],color='r',label='left')
plt.plot(np.arange(len(data[0,:])),data[1,:],color='b',label='middle')
plt.plot(np.arange(len(data[0,:])),data[2,:],color='g',label='right')
plt.legend()
plt.show()
#            exit(0)


# python pyplots/number-proteins-section.py randst 0.25 18.50 18.50 95.00 15.00 full_array 0.00
# python pyplots/number-proteins-section.py randst 0.25 18.50 18.50 95.00 15.00 exact 0.00

# python pyplots/number-proteins-section.py randst 0.25 18.60 28.60 94.00 15.00 full_array 0.00
# python pyplots/number-proteins-section.py randst 0.25 18.60 28.60 94.00 15.00 exact 0.00

# pyplots/number-proteins-section.py stad 0.25 2.35 1.32 0.00 15.00 full_array 0.00
# python pyplots/number-proteins-section.py stad 0.25 2.35 1.32 0.00 15.00 exact 0.00

# pyplots/number-proteins-section.py stad 0.25 2.92 1.18 0.00 15.00 full_array 0.00
# pyplots/number-proteins-section.py stad 0.25 2.92 1.18 0.00 15.00 exact 0.00

# pyplots/number-proteins-section.py p 3.00 0.50 0.00 0.00 15.00 full_array 0.00
# pyplots/number-proteins-section.py p 3.00 0.50 0.00 0.00 15.00 exact 0.00
