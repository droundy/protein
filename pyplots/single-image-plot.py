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
import mycolormap

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
        #image_data_file = job_string +str(proteinList[i])+'/images/real-gauss-single-'+str(times[j])+'.dat'
        temp = np.loadtxt(image_data_file) #this is in microns green light at 500nm,
        smeared_data[j] = np.zeros_like(temp[::4,::4])
        smeared_data[j][:,:] = temp[::4,::4] # keep every fourth data point
    dZ = smeared_data[0].shape[1]*1.05
    dY = smeared_data[0].shape[0]*2.00/skip_times
    Z, Y = np.meshgrid(np.arange(0,smeared_data[0].shape[1],1),
                       np.arange(0,smeared_data[0].shape[0],1))
    maxval = timemax(smeared_data)
    print 'maxval for', proteinList[i], 'is', maxval
    mylevels = np.linspace(0,(1+1.0/nlevels)*maxval,nlevels)
    mylevels = mycolormap.movie_levels
    for k in range(len(smeared_data)):
        #print i, ' ',k
        page = smeared_data[k]
        page[page>maxval] = maxval

        plt.contourf(Y+k*dY, smeared_data[0].shape[1]-Z+i*dZ, page, cmap=mycolormap.cmap, levels=mylevels) # FIXME levels=mylevels
        #plt.contourf(Y+k*dY, smeared_data[0].shape[1]-Z+i*dZ, page, cmap=plt.cm.hot_r, levels=mylevels)

plt.axes().get_yaxis().set_ticks([(i+0.5)*dZ for i in range(len(proteinList))])
plt.axes().get_yaxis().set_ticklabels(proteinLabels)
plt.axes().get_xaxis().set_ticks([(0.5+k)*dY for k in range(len(smeared_data))[::int(2.5*skip_times)]])
plt.axes().get_xaxis().set_ticklabels(range(len(smeared_data))[::int(2.5*skip_times)])
plt.axes().set_aspect('equal')
plt.axes().get_yaxis().set_ticks_position('none')
plt.axes().get_xaxis().set_ticks_position('none')
if sim_type == 'full_array':
    plt.axes().set_title('Stochastic Model',x=.2)
elif sim_type == 'exact':
    plt.axes().set_title('Deterministic Model',x=.2)
else:
    print '\nSomething went wrong while writing the title of the plots\n\n'
    exit(1)
plt.xlabel('time (s)')

# I got the colorbar bit here from
# http://matplotlib.org/users/tight_layout_guide.html, which explains
# the difficulties with tight_layout() and colorbar positioning.

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "3%", pad="1%")
# cbar = plt.colorbar(cax=cax,ticks=[0,maxval])
# cbar.ax.set_yticklabels(['0', 'max'])
cbar = plt.colorbar(cax=cax)
cbar.set_label('molecules/$\mu m^2$')

#plt.colorbar(ticks=[])
plt.tight_layout()
print 'saving as', load.print_string("single-image-plot","")
plt.savefig(load.print_string("single-image-plot",""), transparent=True)
plt.savefig(load.print_string("single-image-plot","")[:-4]+'.eps', transparent=True)
plt.savefig(load.print_string("single-image-plot","")[:-4]+'.png', transparent=True, dpi=500)
#plt.savefig(load.print_string("real-gauss-single-image-plot",""))
