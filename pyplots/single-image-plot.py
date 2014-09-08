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
import re

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

job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (load.f_shape,load.f_param1,load.f_param2,
                                                   load.f_param3,load.f_param4,load.f_param5,sim_type)
p = re.compile('[.]')
job_string = p.sub('_',job_string)

proteinList = ["nADP","nATP","ND","NDE","nE",]
proteinLabels = ["MinD:ADP (cyto)","MinD:ATP (cyto)","MinD:ATP (mem)","MinD:MinE:ATP (mem)","MinD:E (cyto)",]

plt.figure(figsize=(9,3.5))
times  = np.arange(float(start_time),float(end_time),1.0)
print "Starting single image plot."
print times


#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(data):
    return max([v.max() for v in data])

plt.figure(figsize=(9,3.5))
numtimes = int(end_time/dump_time_step)- int(start_time/dump_time_step)
numproteins = len(proteinList)

#axes = plt.subplot(1, 1, 1, axisbg='black')

for i in range(len(proteinList)):
    print proteinList[i]
    nlevels = 20
#generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    smeared_data = [0]*len(times)
    for j in range(len(times)):
        image_data_file = job_string +str(proteinList[i])+'/images/single-'+str(times[j])+'.dat'
        print image_data_file
        smeared_data[j] = np.loadtxt(image_data_file) #this is in microns green light at 500nm,
    dZ = smeared_data[0].shape[1]*1.05
    dY = smeared_data[0].shape[0]*2.00/skip_times
    Z, Y = np.meshgrid(np.arange(0,smeared_data[0].shape[1],1),
                       np.arange(0,smeared_data[0].shape[0],1))
    maxval = timemax(smeared_data)
    mylevels = np.linspace(0,(1+1.0/nlevels)*maxval,nlevels)
    for k in range(len(smeared_data)):
        #print i, ' ',k
        page = smeared_data[k]
        page[page>maxval] = maxval
        cdata = np.array([[0  ,1,1,1],
                          [.01 ,1,1,1],
                          [.25,0.8,.8,1],
                          [.5 ,0,.8,.8],
                          [.7 ,1,1,0],
                          [.85 ,1,0,0],
                          [1  ,0,0,0]])
        cdict = {'red':   [], 'green': [], 'blue':  []}
        for xi in range(cdata.shape[0]):
            cdict['red']   += [(cdata[xi, 0], cdata[xi, 1], cdata[xi, 1])]
            cdict['green'] += [(cdata[xi, 0], cdata[xi, 2], cdata[xi, 2])]
            cdict['blue']  += [(cdata[xi, 0], cdata[xi, 3], cdata[xi, 3])]
        cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)
        plt.contourf(Y+k*dY, smeared_data[0].shape[1]-Z+i*dZ, page, cmap=cmap, levels=mylevels)
        #plt.contourf(Y+k*dY, smeared_data[0].shape[1]-Z+i*dZ, page, cmap=plt.cm.hot_r, levels=mylevels)

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
#plt.savefig(load.print_string("real-gauss-single-image-plot",""))
