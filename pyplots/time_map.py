from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys
import math
import re

# WIP

def gaussian_smear(data,wavelength):
    new = np.zeros_like(data)
    N_A = 1.3
    sigma = .21*wavelength/N_A #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    print "new way sigma ",sigma
    dis = int(3*sigma/0.05) #for now
    for x in range(new.shape[0]):
        for y in range(new.shape[1]):
            for i in np.arange(-dis,dis,1):
                for j in np.arange(-dis,dis,1):
                    if (x+i >= 0 and x+i < new.shape[0]-1 and y+j >= 0 and y+j < new.shape[1]-1):
                        new[x+i,y+j] += data[x,y]*math.exp( -(i*i+j*j)*.05*.05/2.0/sigma/sigma )
    return new


proteinList = ['NflD']

for protein in proteinList:
    job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (load.f_shape,load.f_param1,load.f_param2,
                                                       load.f_param3,load.f_param4,load.f_param5,load.sim_type)
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)
    data_filename = job_string + protein + '/' + 'time-map.dat'
    data_file = np.loadtxt(data_filename)

    print "starting time_map.p for ",data_filename

    data = gaussian_smear(data_file,.509)
    timemax = np.max(data)

    plt.figure()
    plt.pcolormesh(data,vmin=0)
    plt.axes().set_aspect('equal')
    plt.xlim(0,data.shape[1])
    plt.ylim(0,data.shape[0])
    plt.xlabel("Z grid position")
    plt.ylabel("Y grid position")
    plt.title("Time averaged view of %s"%(protein))
    plt.colorbar()
    plt.savefig(load.print_string("time-map",protein))
