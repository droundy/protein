from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import file_loader as load
import sys
import math

# WIP

def gaussian_smear(data,wavelength):
    new = np.zeros_like(data)
    n_sin_theta = 1.5
    sigma = wavelength/2.0/n_sin_theta #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    dis = int(3*sigma/0.05) #for now
    for x in range(new.shape[0]):
        for y in range(new.shape[1]):
            for i in np.arange(-dis,dis,1):
                for j in np.arange(-dis,dis,1):
                    if (x+i >= 0 and x+i < new.shape[0]-1 and y+j >= 0 and y+j < new.shape[1]-1):
                        new[x+i,y+j] += data[x,y]*math.exp( -(i*i+j*j)*.05*.05/sigma/sigma )
    return new


proteinList = ["NflD"]

for protein in proteinList:
    data_file = np.loadtxt("./data/shape-%s/time-map-%s-%s-%s-%s-%s-%s-%s-%s.dat"%(load.f_shape,protein,load.f_shape,load.f_param1,load.f_param2,load.f_param3,load.f_param4,load.f_param5,load.sim_type))

    data = gaussian_smear(data_file,.6)
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
