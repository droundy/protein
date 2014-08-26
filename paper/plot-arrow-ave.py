from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
import matplotlib.image as mpimg
import os
import sys
import time
import imp
#load = imp.load_source('load', 'pyplots/file_loader.py')
import math
import matplotlib.gridspec as gridspec
import re

#create data objects (see file_loader.py)





dx =0.05
dump_time_step = 0.5
protein_name = "NflD"
# start_time = float(f_param6)
# input_end_time = float(f_param7)
#end_time = int(input_end_time - input_end_time%10)

arg_set = ["0.25-18.50-18.50-95.00-15.00-full_array",
           "0.25-18.60-28.60-94.00-15.00-full_array",
           "0.25-18.50-18.50-95.00-15.00-exact",
           "0.25-18.60-28.60-94.00-15.00-exact"]

bound_times = ["500","1500","500","540","500","520","500","520"]

arrow_files = []
contour_values = []
for i in range(len(arg_set)):
    job_string = "data/shape-randst/%s/" % (arg_set[i])
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)
    arrow_files += [job_string +'ave-time/ave-time-arrow-'+str(int(bound_times[i*2]))+'-'+protein_name+'.dat']
    contour_values += [ job_string +'ave-time/contour-values-' + protein_name +'-'+ str(int(bound_times[i*2]))+'-' \
                            +str(bound_times[i*2+1])+'.dat']
    print i
    print arrow_files[i]
    print contour_values[i]



plt.figure(1)

#plt.clf()
print "len",len(arg_set)

# half_width = float(c_data.shape[1]/max_width)/2.0
# print half_width
#gs.update(left=.5-half_width,right=.5+half_width)
#plt.title("%f"%half_width)

hts = []
for i in range(len(contour_values)):
    c_data = np.loadtxt(contour_values[i])
    hts += [float(c_data.shape[1])]

if len(contour_values) > 1:
    hts += [hts[1]]

gs = gridspec.GridSpec(len(arg_set)+2,2)
plt.subplots_adjust(left=.1,right=.9,top=1.0,bottom=-.9)#left=.3,right=.7)
#gs.tight_layout(w_pad=.2)

for arg_num in range(len(arg_set)+2):
#    plt.subplot2grid((len(arg_set)+1,1),(arg_num,0),aspect='equal')
    #gs_in = gridspec.GridSpec(arg_num,0)
    #gs_in.update(left=.3+arg_num*.1,right=.8)
    if arg_num == len(arg_set):
        plt.subplot(gs[0,0],aspect='equal')
        m1=mpimg.imread('mannik-1.png')
        #dx = 7.0/380
        print "shape, ",m1.shape[1]
        y1 = -dx*np.arange(0, m1.shape[0])
        x1 = dx*np.arange(0, m1.shape[1])
        x1, y1 = np.meshgrid(x1, y1)
        xcenter = x1[m1[:,:,0] > 0].mean()
        y1off = 1.7#1.7 and there was only one yoff, not two
        plt.contourf(x1-xcenter, y1-y1off, m1[:,:,1], levels=[0, 0.9], colors=('r', 'w'))
        plt.contourf(x1-xcenter, y1-y1off, m1[:,:,0], levels=[0, 0.9], colors=('k', 'w'))
    elif arg_num == len(arg_set)+1:
        plt.subplot(gs[0,1],aspect='equal')
        m1=mpimg.imread('mannik-2.png')
        dx = 7.0/380
        print m1.shape[1]
        y1 = -dx*np.arange(0, m1.shape[0])
        x1 = dx*np.arange(0, m1.shape[1])
        x1, y1 = np.meshgrid(x1, y1)
        xcenter = x1[m1[:,:,0] > 0].mean()
        y1off = 1.7#1.7 and there was only one yoff, not two
        plt.contourf(x1-xcenter, y1-y1off, m1[:,:,1], levels=[0, 0.9], colors=('r', 'w'))
        plt.contourf(x1-xcenter, y1-y1off, m1[:,:,0], levels=[0, 0.9], colors=('k', 'w'))
    else:
        plt.subplot(gs[int(arg_num/2)+1,arg_num%2],aspect='equal')
        print int(arg_num/2)+1,"",arg_num%2
        c_data = np.loadtxt(contour_values[arg_num])
        a_data = np.loadtxt(arrow_files[arg_num])
        start_time = float(bound_times[arg_num*2])
        end_time = float(bound_times[arg_num*2+1])

        a_data = a_data[int((start_time-float(bound_times[arg_num*2]))/dump_time_step):int((end_time-start_time)/dump_time_step)]

        last_time = a_data[0,0]
        for i in range(1,len(a_data[:,0])):
            if (a_data[i,0] <= last_time):
                print "The times in the arrow file are not in chronological order!  Something's wrong! Exiting."
                print i
                print last_time
                print a_data[i,0]
                exit(0)

        time_max = np.max(c_data)
        arrow_cutoff = 3.0*(np.max(a_data[:,1]))/5.0

        high_maximas = np.zeros(0)
        times = np.zeros(0)
        x_vals = np.zeros(0)
        y_vals = np.zeros(0)
        last_x = 0
        last_y = 0
        index = 0
        # print np.max(a_data[:,1])
        # print "arrow cutoff ",arrow_cutoff
        for i in range(len(a_data[:,1])):
            if (a_data[i,1] > arrow_cutoff):
                if ( (a_data[i,2]*dx != last_x or a_data[i,3]*dx != last_y) \
                         and math.sqrt((a_data[i,2]*dx - last_x)**2 + (a_data[i,3]*dx-last_y)**2) > .2):
                    x_vals = np.append(x_vals,a_data[i,2]*dx)
                    y_vals = np.append(y_vals,a_data[i,3]*dx)
                    times = np.append(times,a_data[i,0])
                    last_x = a_data[i,2]*dx
                    last_y = a_data[i,3]*dx

        print "arg ",arg_num


    #cbar = plt.colorbar(CS)



        Ny = len(c_data[:,0])
        Nz = len(c_data[0,:])

        Z, Y = np.meshgrid(np.arange(0,(c_data.shape[1]-.9)*dx,dx),np.arange(0,(c_data.shape[0]-.9)*dx,dx))
    #plt.contourf(originx+cell_x,originy+cell_y,c_data, linewidths=2,levels=[.99])
        plt.contourf(Z, Y, c_data, cmap=plt.cm.jet,origin='lower',levels=np.arange(0,time_max+1.0,1))


        for i in range(len(x_vals)-1):
            plt.annotate('%g'%i,xy=(y_vals[i+1],x_vals[i+1]),xytext=(y_vals[i],x_vals[i]),
                         fontsize=11,
                         arrowprops=dict(color='red',shrink=0.01, width=.3, headwidth=5.))
            plt.clim(0,time_max)
    plt.axis('off')
    #plt.savefig('./paper/plot-ave-%d.pdf'%arg_num)



#ax.get_xaxis().set_visible(True)
#plt.xlim((0,dx*c_data.shape[1]))
#plt.ylim((0,dx*c_data.shape[0]))
plt.subplots_adjust(hspace=0.02)
# left  = 0.125  # the left side of the subplots of the figure
# right = 0.9    # the right side of the subplots of the figure
# bottom = 0.1   # the bottom of the subplots of the figure
# top = 0.9      # the top of the subplots of the figure
# wspace = 0.2   # the amount of width reserved for blank space between subplots
# hspace = 0.2   # the amount of height reserved for white space between subplots
#plt.tight_layout()
plt.xlabel("Z grid position")
plt.ylabel("Y grid position")
#plt.title("Local temporal maxima, global spatial maxima view of MinD")

plt.savefig('./paper/plot-ave.pdf')

plt.show()

