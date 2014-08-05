from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load
import Image
import math

f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
dens_factor = sys.argv[6]
sim_type = sys.argv[7]
f_param6 = sys.argv[8]
f_param7 = sys.argv[9]

dump_time_step = 0.5
skip_times = 2 # only plot every skip_times of the snapshots (which are seperated by dump_time_steps (.5) simulation seconds.)
start_time = float(f_param6)
end_time = float(f_param7)


def gaussian_smear(data,wavelength):
    new = np.zeros_like(data)
    new = new[0:len(new)/skip_times]
    n_sin_theta = 1.5
    sigma = wavelength/2.0/n_sin_theta #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    dis = int(3*sigma/0.05) #for now
    for num in range(new.shape[0]):
        print "num ",num," of ",new.shape[0]-1
        for x in range(new.shape[1]):
            print x," of ",new.shape[1]-1
            for y in range(new.shape[2]):
                for i in np.arange(-dis,dis,1):
                    for j in np.arange(-dis,dis,1):
                        if (x+i >= 0 and x+i < new.shape[1]-1 and y+j >= 0 and y+j < new.shape[2]-1):
                            new[num,x+i,y+j] += data[num*skip_times,x,y]*math.exp( -(i*i+j*j)*.05*.05/sigma/sigma )
    return new


#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(protein):
    return max([v.max() for v in protein])

proteinList = ["nADP","nATP","ND","NDE","nE",]
#proteinList = ["nATP"]
proteinLabels = ["MinD:ADP (cyto)","MinD:ATP (cyto)","MinD:ATP (mem)","MinD:MinE:ATP (mem)","MinD:E (cyto)",]

proteins = [0]*len(proteinList)

for i in range(len(proteinList)):
    proteins[i] = load.data(proteinList[i],sim_type,start_time,end_time)

plt.figure(figsize=(9,3.5))
numtimes = int(end_time/dump_time_step)- int(start_time/dump_time_step)
numproteins = len(proteins)
dZ = proteins[0].datashape[1]*1.05
dY = proteins[0].datashape[0]*2.00/skip_times

#axes = plt.subplot(1, 1, 1, axisbg='black')

for i in range(len(proteins)):
    print proteinList[i]
    nlevels = 20
    Z, Y = np.meshgrid(np.arange(0,proteins[i].dataset[0].shape[1],1),
                       np.arange(0,proteins[i].dataset[0].shape[0],1))
#generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    smeared_data = gaussian_smear(proteins[i].dataset,.6) #this is in microns green light at 500nm,
    # if (proteinList[i]=="ND" or proteinList[i]=="NDE"):
    #     maxval = timemax(proteins[i].dataset)/4
    # else:
    #     maxval = timemax(proteins[i].dataset)
    maxval = timemax(smeared_data)
    mylevels = np.linspace(0,(1+1.0/nlevels)*maxval,nlevels)
    for k in range(len(smeared_data)):
        page = smeared_data[k]
        page[page>maxval] = maxval
        plt.contourf(Y+k*dY, proteins[0].datashape[1]-Z+i*dZ, page, cmap=plt.cm.hot_r, levels=mylevels)
plt.axes().get_yaxis().set_ticks([(i+0.5)*dZ for i in range(len(proteinList))])
plt.axes().get_yaxis().set_ticklabels(proteinLabels)
plt.axes().get_xaxis().set_ticks([(0.5+k)*dY for k in range(len(smeared_data))[::int(2.5*skip_times)]])
plt.axes().get_xaxis().set_ticklabels(range(len(smeared_data))[::int(2.5*skip_times)])
plt.axes().set_aspect('equal')
plt.axes().get_yaxis().set_ticks_position('none')
plt.axes().get_xaxis().set_ticks_position('none')
plt.xlabel('time (s)')

# I got the colorbar bit here from
# http://matplotlib.org/users/tight_layout_guide.html, which explains
# the difficulties with tight_layout() and colorbar positioning.

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="1%")
cbar = plt.colorbar(cax=cax,ticks=[0,maxval])
cbar.ax.set_yticklabels(['0', 'max'])

#plt.colorbar(ticks=[])
plt.tight_layout()
plt.savefig(load.print_string("single-image-plot",""))
