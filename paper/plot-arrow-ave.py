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
arg_set = ["randst/0.25-18.50-18.50-95.00-15.00-exact",
           "randst/0.25-18.50-18.50-95.00-15.00-full_array",
           "randst/0.25-18.60-28.60-94.00-15.00-exact",
           "randst/0.25-18.60-28.60-94.00-15.00-full_array",
           "stad/0.25-2.35-1.32-0.00-15.00-full_array",
           "stad/0.25-2.35-1.32-0.00-15.00-exact",
           "stad/0.25-2.92-1.18-0.00-15.00-full_array",
           "stad/0.25-2.92-1.18-0.00-15.00-exact",
           "p/3.00-0.50-0.00-0.00-15.00-exact",
           "p/3.00-0.50-0.00-0.00-15.00-full_array"]

#rot_theta = [math.pi/2.0,math.pi/3.0,math.pi/4.0,math.pi/5.0,7.0*math.pi/6.0,21.0*math.pi/9.0]

bound_times = [500,1000,500,1000,500,1000,500,980,500,1000,500,780,500,1000,500,780,500,700,500,700]
arrow_goals = [20, 16, 16, 14, 19, 18, 18, 16, 5, 5]

col_0x = 1.0
col_1x = 8.0
col_2x = 15.0
col_3x = 21.0
col_4x = 26.5

row_0y = 1.0
row_1y = 6.0

x_position_m1 = col_0x
y_position_m1 = 13.0
x_position_m2 = col_1x
y_position_m2 = 13.0
x_position_m3 = col_4x
y_position_m3 = 13.0

X_position = [col_0x,col_0x,col_1x,col_1x,col_2x,col_2x,col_3x,col_3x,col_4x,col_4x]
Y_position = [row_0y,row_1y,row_0y,row_1y,row_0y,row_1y,row_0y,row_1y,row_0y,row_1y]

theta_0 = 0#8.0*math.pi/8.0
theta_1 = 0#6.78*math.pi/8.0
theta_2 = 0#1.1235813*math.pi/8.0
theta_3 = 0#1.1235813*math.pi/8.0
theta_4 = 0

rot_theta = [theta_0,theta_0,theta_1,theta_1,theta_2,theta_2,theta_3,theta_3,theta_4,theta_4]
mannik_rot_theta = [theta_0,theta_1,theta_4]

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
plt.axes().set_aspect('equal', 'datalim')

data = [0]*(len(arg_set)+2)
print len(data)
for arg_num in range(len(arg_set)):
    data[arg_num] = np.loadtxt(contour_values[arg_num])


mannik_micron = 5.0/1012
m1=mpimg.imread('mannik-1.png')
m1x = -np.arange(-m1.shape[0]/2, m1.shape[0]/2)*mannik_micron
m1y = np.arange(-m1.shape[1]/2, m1.shape[1]/2)*mannik_micron
M1x, M1y = np.meshgrid(m1y, m1x)
M1xrot = M1x*math.cos(mannik_rot_theta[0]) - M1y*math.sin(mannik_rot_theta[0]) + x_position_m1
M1yrot = M1x*math.sin(mannik_rot_theta[0]) + M1y*math.cos(mannik_rot_theta[0]) + y_position_m1
plt.contourf(M1xrot, M1yrot, m1[:,:,0], levels=[0, 0.9], colors=('r', 'w'))
plt.contourf(M1xrot, M1yrot, m1[:,:,1], levels=[0, 0.9], colors=('k', 'w'))


m2=mpimg.imread('mannik-2.png')
m2x = -np.arange(-m2.shape[0]/2, m2.shape[0]/2)*mannik_micron
m2y = np.arange(-m2.shape[1]/2, m2.shape[1]/2)*mannik_micron
M2x, M2y = np.meshgrid(m2y, m2x)
M2xrot = M2x*math.cos(mannik_rot_theta[1]) - M2y*math.sin(mannik_rot_theta[1]) + x_position_m2
M2yrot = M2x*math.sin(mannik_rot_theta[1]) + M2y*math.cos(mannik_rot_theta[1]) + y_position_m2
plt.contourf(M2xrot, M2yrot, m2[:,:,0], levels=[0, 0.9], colors=('r', 'w'))
plt.contourf(M2xrot, M2yrot, m2[:,:,1], levels=[0, 0.9], colors=('k', 'w'))

# mannik_micron3 = 4.0/400
# m3=mpimg.imread('mannik-3.png')
# m3x = -np.arange(-m3.shape[0]/2, m3.shape[0]/2)*mannik_micron3
# m3y = np.arange(-m3.shape[1]/2, m3.shape[1]/2)*mannik_micron3
# M3x, M3y = np.meshgrid(m3y, m3x)
# M3xrot = M3x*math.cos(mannik_rot_theta[2]) - M3y*math.sin(mannik_rot_theta[2]) + x_position_m3
# M3yrot = M3x*math.sin(mannik_rot_theta[2]) + M3y*math.cos(mannik_rot_theta[2]) + y_position_m3
# plt.contourf(M3xrot, M3yrot, m3[:,:,0], levels=[0, 0.9], colors=('r', 'w'))
# plt.contourf(M3xrot, M3yrot, m3[:,:,1], levels=[0, 0.9], colors=('k', 'w'))


arrow_file = open('arrow_printout.txt', 'w')

#rotate and plot sim data
for arg_num in range(len(arg_set)):
    x = np.arange(-data[arg_num].shape[0]/2, data[arg_num].shape[0]/2)*dx
    y = np.arange(-data[arg_num].shape[1]/2, data[arg_num].shape[1]/2)*dx
    X,Y = np.meshgrid(y,x)
    Xrot = X*math.cos(rot_theta[arg_num]) - Y*math.sin(rot_theta[arg_num]) + X_position[arg_num]
    Yrot = X*math.sin(rot_theta[arg_num]) + Y*math.cos(rot_theta[arg_num]) + Y_position[arg_num]
    cdata = np.array([[0  ,1,1,1],
                      [.01 ,1,1,1],
                      [.25,0.8,.8,1],
                      [.5 ,0,.8,.8],
                      [.7 ,1,1,0],
                      [.85 ,1,0,0],
                      [1  ,0,0,0]])
    cdict = {'red':   [], 'green': [], 'blue':  []}
    for i in range(cdata.shape[0]):
        cdict['red']   += [(cdata[i, 0], cdata[i, 1], cdata[i, 1])]
        cdict['green'] += [(cdata[i, 0], cdata[i, 2], cdata[i, 2])]
        cdict['blue']  += [(cdata[i, 0], cdata[i, 3], cdata[i, 3])]
    cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)
    time_max = np.max(data[arg_num])

    plt.contourf(Xrot, Yrot, data[arg_num], cmap=cmap, origin='lower',levels=np.arange(0,time_max+1.0,1))

    start_time = float(bound_times[arg_num*2])
    end_time = float(bound_times[arg_num*2+1])

    a_data = np.loadtxt(arrow_files[arg_num])

    a_data = a_data[int((start_time-float(bound_times[arg_num*2]))/dump_time_step):int((end_time-start_time)/dump_time_step)]
    temp_maxes = np.zeros_like(a_data)
    j = 0
    for i in range(1,len(a_data)-1):
        if a_data[i-1,1] < a_data[i,1] < a_data[i+1,1]:
            temp_maxes[j] = a_data[i]
            j+=1
    temp_maxes = temp_maxes[:j,:]
    a_data = temp_maxes

    last_time = a_data[0,0]
    for i in range(1,len(a_data[:,0])):
        if (a_data[i,0] <= last_time):
            print "The times in the arrow file are not in chronological order!  Something's wrong! Exiting."
            print i
            print last_time
            print a_data[i,0]
            exit(0)

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
        amount = np.zeros(0)
        x_vals = np.zeros(0)
        y_vals = np.zeros(0)
        for i in range(len(clean_data[:,1])-1):
            x_vals = np.append(x_vals,clean_data[i,3]*dx)
            y_vals = np.append(y_vals,clean_data[i,2]*dx)
            amount = np.append(amount,clean_data[i,1]*dx)
            times = np.append(times,clean_data[i,0])
            #print 'finished with arrow_cutoff', arrow_cutoff, 'and num arrows', len(x_vals)

    arrow_file.write(arg_set[arg_num])
    for i in range(len(times)):
        arrow_file.write(' %g %g %g %g\n'%(times[i],amount[i],x_vals[i],y_vals[i]))


    x_vals += y[0]
    y_vals += x[0]
    x_vals,y_vals = x_vals*math.cos(rot_theta[arg_num]) - y_vals*math.sin(rot_theta[arg_num]), \
        +x_vals*math.sin(rot_theta[arg_num]) + y_vals*math.cos(rot_theta[arg_num])

    x_vals += X_position[arg_num]
    y_vals += Y_position[arg_num]

    #plt.colorbar()
    #plt.bar(.2,6,.2)

    for i in range(len(x_vals)-1):
        plt.annotate('',xy=(x_vals[i+1],y_vals[i+1]),xytext=(x_vals[i],y_vals[i]),
                     fontsize=11,
                     arrowprops=dict(color='red',shrink=0.01, width=.3, headwidth=5.))
        plt.clim(0,time_max)
arrow_file.close()
plt.axis('off')
plt.savefig('./paper/plot-ave.pdf')

#plt.show()


# arrow_file.close()
# #ax.get_xaxis().set_visible(True)
# #plt.xlim((0,dx*c_data.shape[1]))
# #plt.ylim((0,dx*c_data.shape[0]))
# plt.subplots_adjust(hspace=0.02)
# # left  = 0.125  # the left side of the subplots of the figure
# # right = 0.9    # the right side of the subplots of the figure
# # bottom = 0.1   # the bottom of the subplots of the figure
# # top = 0.9      # the top of the subplots of the figure
# # wspace = 0.2   # the amount of width reserved for blank space between subplots
# # hspace = 0.2   # the amount of height reserved for white space between subplots
# #plt.tight_layout()
# plt.xlabel("Z grid position")
# plt.ylabel("Y grid position")
# #plt.title("Local temporal maxima, global spatial maxima view of MinD")

