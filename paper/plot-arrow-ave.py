from __future__ import division
import numpy as np
import matplotlib
#matplotlib.use('Agg')
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


arrow_file = open('arrow-file.txt','w')


# arg_set = ["randst/0.25-18.50-18.50-95.00-15.00-full_array",
#            "randst/0.25-18.60-28.60-94.00-15.00-full_array",
#            "stad/0.25-3.00-1.18-0.00-15.00-full_array",
#            "randst/0.25-18.50-18.50-95.00-15.00-exact",
#            "randst/0.25-18.60-28.60-94.00-15.00-exact",
#            "stad/0.25-3.00-1.18-0.00-15.00-exact"]
arg_set = ["randst/0.25-18.50-18.50-95.00-15.00-full_array",
           "randst/0.25-18.60-28.60-94.00-15.00-full_array",
           "randst/0.25-18.50-18.50-95.00-15.00-exact",
           "randst/0.25-18.60-28.60-94.00-15.00-exact"]
rot_theta = [math.pi/2.0,math.pi/3.0,math.pi/4.0,math.pi/5.0,7.0*math.pi/6.0,21.0*math.pi/9.0]
rot_theta = [0.0*math.pi/8.0,0,0,0,0,0]

bound_times = [500,1000,500,980,500,1000,500,1000]
arrow_cutoffs = [2.4/5.0,2.5/5.0,2.2/5.0,2.9/5.0]
arrow_cutoffs = [180]*2 + [57]*2
arrow_goals = [20, 16, 16, 14]

arrow_files = []
contour_values = []
for i in range(len(arg_set)):
    job_string = "data/shape-%s/" % (arg_set[i])
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)
    arrow_files += [job_string +'ave-time/ave-time-arrow-'+str(int(bound_times[i*2]))+'-'+protein_name+'.dat']
    contour_values += [ job_string +'ave-time/contour-values-' + protein_name +'-'+ str(int(bound_times[i*2]))+'-' \
                            +str(bound_times[i*2+1])+'.dat']



plt.figure(1, figsize=(8,5))

#plt.clf()


for i in range(len(contour_values)):
    c_data = np.loadtxt(contour_values[i])

gs = gridspec.GridSpec(len(arg_set)+2,2)
plt.subplots_adjust(left=.1,right=.9,top=1.0,bottom=-.9)#left=.3,right=.7)

for arg_num in range(len(arg_set)+2):
    mannik_micron = 5.0/1012
    if arg_num == len(arg_set):
        gs.update(wspace=0)
        plt.subplot(gs[0,0],aspect='equal')
        m1=mpimg.imread('mannik-1.png')
        print "shape, ",m1.shape[1]
        y1 = -np.arange(0, m1.shape[0])*mannik_micron
        x1 = np.arange(0, m1.shape[1])*mannik_micron
        x1, y1 = np.meshgrid(x1, y1)
        xrot = x1*math.cos(rot_theta[arg_num]) - y1*math.sin(rot_theta[arg_num])
        yrot = x1*math.sin(rot_theta[arg_num]) + y1*math.cos(rot_theta[arg_num])
        xcenter = x1[m1[:,:,0] > 0].mean()
        y1off = 1.7#1.7 and there was only one yoff, not two
        plt.contourf(xrot-xcenter, yrot-y1off, m1[:,:,1], levels=[0, 0.9], colors=('r', 'w'))
        plt.contourf(xrot-xcenter, yrot-y1off, m1[:,:,0], levels=[0, 0.9], colors=('k', 'w'))
    elif arg_num == len(arg_set)+1:
        plt.subplot(gs[0,1],aspect='equal')
        m1=mpimg.imread('mannik-2.png')
        print m1.shape[1]
        y1 = -np.arange(0, m1.shape[0])*mannik_micron
        x1 = np.arange(0, m1.shape[1])*mannik_micron
        x1, y1 = np.meshgrid(x1, y1)
        xrot = x1*math.cos(rot_theta[arg_num]) - y1*math.sin(rot_theta[arg_num])
        yrot = x1*math.sin(rot_theta[arg_num]) + y1*math.cos(rot_theta[arg_num])
        xcenter = xrot[m1[:,:,0] > 0].mean()
        y1off = 1.7#1.7 and there was only one yoff, not two
        plt.contourf(xrot-xcenter, yrot-y1off, m1[:,:,1], levels=[0, 0.9], colors=('r', 'w'))
        plt.contourf(xrot-xcenter, yrot-y1off, m1[:,:,0], levels=[0, 0.9], colors=('k', 'w'))
    else:
        ax = plt.subplot(gs[int(arg_num/2)+1,arg_num%2],aspect='equal')
        if arg_num == len(arg_set)-2:
            ax.add_artist(AnchoredSizeBar(
                    ax.transData,
                    5, # length of the bar in the data reference
                    "5$\mu$", # label of the bar
                    # bbox_to_anchor=(0.,0.,1.,1.),
                    loc=8, # 'best', # location (lower right)
                    borderpad=-1.8, sep=3,
                    frameon=False
                    ))
        print int(arg_num/2)+1,"",arg_num%2
        first_c_data = np.loadtxt(contour_values[arg_num])
        start_time = float(bound_times[arg_num*2])
        end_time = float(bound_times[arg_num*2+1])

        y_buf = 10
        x_buf = 10
        if arg_num == len(arg_set)-1 or arg_num == len(arg_set)-3:
            y_buf = 8
            x_buf = 8
        c_data = np.zeros((first_c_data.shape[0]+x_buf*2,first_c_data.shape[1]+y_buf*2))
        for i in range(first_c_data.shape[0]):
            for j in range(first_c_data.shape[1]):
                c_data[i+x_buf][j+y_buf] = first_c_data[i][j]

        a_data = np.loadtxt(arrow_files[arg_num])
        a_data = np.loadtxt(arrow_files[arg_num])
        a_data[:,2] += x_buf
        a_data[:,3] += y_buf

        a_data = a_data[int((start_time-float(bound_times[arg_num*2]))/dump_time_step):int((end_time-start_time)/dump_time_step)]
        temp_maxes = np.zeros_like(a_data)
        j = 0
        for i in range(1,len(a_data)-1):
            if a_data[i-1,1] < a_data[i,1] < a_data[i+1,1]:
                temp_maxes[j] = a_data[i]
                j+=1
        temp_maxes = temp_maxes[:j,:]
        print temp_maxes.shape
        print a_data.shape
        a_data = temp_maxes

        last_time = a_data[0,0]
        for i in range(1,len(a_data[:,0])):
            if (a_data[i,0] <= last_time):
                print "The times in the arrow file are not in chronological order!  Something's wrong! Exiting."
                print i
                print last_time
                print a_data[i,0]
                exit(0)

        time_max = np.max(c_data)
        # arrow_cutoff = arrow_cutoffs[arg_num]*(np.max(a_data[:,1]))
        # print 'silly max is', (np.max(a_data[:,1]))
        # arrow_cutoff = arrow_cutoffs[arg_num]

        arrow_cutoff = 0
        arrow_goal = arrow_goals[arg_num]
        x_vals = [0]*(arrow_goal+1)
        while len(x_vals) > arrow_goal:
            arrow_cutoff += 0.1
            #print 'trying arrow_cutoff', arrow_cutoff, 'with num arrows', len(x_vals)
            large_a_data = np.zeros((len(a_data[a_data[:,1]>arrow_cutoff,:]),4))
            large_a_data[:,:] = a_data[a_data[:,1]>arrow_cutoff,:]

            clean_data = np.zeros_like(large_a_data)
            j = 0
            tooclose = 1.0
            i = 0
            while i < len(large_a_data[:,1]):
                n = 1
                while i+n < len(large_a_data[:,1]) - 1 and \
                      dx*math.sqrt((large_a_data[i+n-1,2] - large_a_data[i+n,2])**2 + (large_a_data[i+n-1,3]-large_a_data[i+n,3])**2) < tooclose:
                    n+=1
                maxn = np.argmax(large_a_data[i:i+n,1])
                clean_data[j,:] = large_a_data[i+maxn,:]
                i+=n
                j+=1
            clean_data = clean_data[1:j-1,:]

            high_maximas = np.zeros(0)
            times = np.zeros(0)
            x_vals = np.zeros(0)
            y_vals = np.zeros(0)
            for i in range(len(clean_data[:,1])-1):
                x_vals = np.append(x_vals,clean_data[i,2]*dx)
                y_vals = np.append(y_vals,clean_data[i,3]*dx)
                times = np.append(times,clean_data[i,0])
        print 'finished with arrow_cutoff', arrow_cutoff, 'and num arrows', len(x_vals)
        print 'len x_vals',len(x_vals)

        # x_cen = len(c_data[:,0])*dx/2.0
        # y_cen = len(c_data[0,:])*dx/2.0
        # x_vals = x_vals[0:3]
        # y_vals = y_vals[0:3]
        # x_vals = x_vals - x_cen
        # y_vals = y_vals - y_cen
        # print x_vals
        # rot_thet = -9.0*math.pi/8.0
        # rot = 9.0*math.pi/8.0
        # x_vals, y_vals = x_cen + (x_vals-x_cen)*math.cos(-rot_thet) + (y_vals-y_cen)*math.sin(-rot_thet), \
        #     y_cen - (x_vals-x_cen)*math.sin(-rot_thet) + (y_vals-y_cen)*math.cos(-rot_thet)
        # print x_vals

        arrow_file.write('starting the arrows for %s'%(arg_set[arg_num]))
        for i in range(len(x_vals)):
            arrow_file.write('t %g, x %g, y %g\n'%(times[i],x_vals[i],y_vals[i]))

        Z, Y = np.meshgrid(np.arange(0,(c_data.shape[1]-.9)*dx,dx),np.arange(0,(c_data.shape[0]-.9)*dx,dx))
        #Z, Y = np.meshgrid(np.arange(-c_data.shape[1]*dx/2.0,(c_data.shape[1])*dx/2.0,dx),np.arange(-c_data.shape[0]*dx/2.0,(c_data.shape[0])*dx/2.0,dx))
        print Z
        print "len(z) ",len(Z)
        print "len(z[0]) ",len(Z[0])
        Zrot = Z*math.cos(rot_theta[arg_num]) - Y*math.sin(rot_theta[arg_num])
        Yrot = Z*math.sin(rot_theta[arg_num]) + Y*math.cos(rot_theta[arg_num])
        print "len(roty) ",len(Yrot)
        print "len(roty[0]) ",len(Yrot[0])
        cdata = np.array([[0  ,1,1,1],
                          [.01 ,1,1,1],
                          [.25,0.8,.8,1],
                          [.5 ,0,.8,.8],
                          [.7 ,1,1,0],
                          [.85 ,1,0,0],
                          [1  ,0,0,0]])
        cdict = {'red':   [],
                 'green': [],
                 'blue':  []}
        for i in range(cdata.shape[0]):
            print 'color', i
            cdict['red']   += [(cdata[i, 0], cdata[i, 1], cdata[i, 1])]
            cdict['green'] += [(cdata[i, 0], cdata[i, 2], cdata[i, 2])]
            cdict['blue']  += [(cdata[i, 0], cdata[i, 3], cdata[i, 3])]
        cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)

        plt.contourf(Zrot, Yrot, c_data, cmap=cmap, origin='lower',levels=np.arange(0,time_max+1.0,1))
        #plt.colorbar()
        #plt.bar(.2,6,.2)
        # for i in Zrot:
        #     print i
        # for i in Yrot:
        #     print i
        print x_vals
        print y_vals
        for i in range(len(x_vals)-1):
            plt.annotate('',xy=(y_vals[i+1],x_vals[i+1]),xytext=(y_vals[i],x_vals[i]),
                         fontsize=11,
                         arrowprops=dict(color='red',shrink=0.01, width=.3, headwidth=5.))
            plt.clim(0,time_max)
    plt.axis('off')
    #plt.savefig('./paper/plot-ave-%d.pdf'%arg_num)

arrow_file.close()
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

#plt.show()

