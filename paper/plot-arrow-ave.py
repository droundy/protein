# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
import sys
import time
import imp
#load = imp.load_source('load', 'pyplots/file_loader.py')
import math
import matplotlib.gridspec as gridspec
import re

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
           "stad/0.25-2.35-1.32-0.00-15.00-exact",
           "stad/0.25-2.35-1.32-0.00-15.00-full_array",
           "stad/0.25-2.92-1.18-0.00-15.00-exact",
           "stad/0.25-2.92-1.18-0.00-15.00-full_array",
           "p/3.00-0.50-0.00-0.00-15.00-exact",
           "p/3.00-0.50-0.00-0.00-15.00-full_array"]

bound_times = [500,1000,500,1000,500,1000,500,980,500,780,500,1000,500,780,500,1000,500,700,500,700]
bound_times = [500,850,500,850,500,850,500,850,500,850,500,850,500,850,500,850,500,840,500,840]
arrow_goals = [20, 16, 16, 14, 19, 18, 18, 16, 11, 11]
arrow_goals = [13, 11, 12, 12, 9, 16, 9, 15, 18, 18]

left_annotate_x = -9.0
bottom_annotate_y = -1.7

row_0y = 1.0
row_1y = 6.0
row_my = 10.5

col_0x = -5.0
col_1x = 1.0
col_2x = 8.0
col_3x = 14.5
col_4x = col_3x + 8.5

pill_bar_x = col_4x + 3
color_bar_x = col_4x - 5

x_position_m1 = col_0x
y_position_m1 = row_my
x_position_m2 = col_2x
y_position_m2 = row_my

viewport_xmin = left_annotate_x - 6
viewport_xmax = pill_bar_x + 2
viewport_ymin = bottom_annotate_y
viewport_ymax = row_my + 2

def in_viewport_x(x):
    return (x - viewport_xmin)/(viewport_xmax - viewport_xmin)
def in_viewport_y(y):
    return (y - viewport_ymin)/(viewport_ymax - viewport_ymin)

X_position = [col_0x,col_0x,col_2x,col_2x,col_1x,col_1x,col_3x,col_3x,col_4x,col_4x]
Y_position = [row_0y,row_1y,row_0y,row_1y,row_0y,row_1y,row_0y,row_1y,row_0y,row_1y]

theta_0 = 0#8.0*math.pi/8.0
theta_1 = 0#6.78*math.pi/8.0
theta_2 = 0#1.1235813*math.pi/8.0
theta_3 = 0#1.1235813*math.pi/8.0
theta_4 = 0

arrow_files = []
contour_values = []

for i in range(len(arg_set)):
    job_string = "data/shape-%s/" % (arg_set[i])
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)
    arrow_files += [job_string +'ave-time/ave-time-arrow-'+str(int(bound_times[i*2]))+'-'+protein_name+'.dat']
    contour_values += [ job_string +'ave-time/contour-values-' + protein_name +'-'+ str(int(bound_times[i*2]))+'-' \
                            +str(bound_times[i*2+1])+'.dat']


plt.figure(1, figsize=(8,3.4))

ax = plt.axes([0,0,1,1]) # for nice figure
#ax = plt.axes([0.1,0.1,.9,.9]) # for debugging positions
ax.set_aspect('equal')
ax.set_axis_off()

data = [0]*(len(arg_set)+2)
print len(data)
for arg_num in range(len(arg_set)):
    data[arg_num] = np.loadtxt(contour_values[arg_num])


arrow_file = open('arrow_printout.txt', 'w')

#rotate and plot sim data
for arg_num in range(len(arg_set)):
    x = np.arange(0., data[arg_num].shape[0])*dx
    y = np.arange(0., data[arg_num].shape[1])*dx
    X,Y = np.meshgrid(y,x)

    # Find the quadrupole moment and diagonlize to get a rotation matrix
    cellshape = np.zeros_like(X)
    cellshape[:,:] = data[arg_num]>0
    meanX = np.sum(X*cellshape)/np.sum(cellshape)
    meanY = np.sum(Y*cellshape)/np.sum(cellshape)
    X -= meanX
    Y -= meanY
    Qxx = np.sum(X*X*cellshape)/np.sum(cellshape)
    Qyy = np.sum(Y*Y*cellshape)/np.sum(cellshape)
    Qxy = np.sum(X*Y*cellshape)/np.sum(cellshape)
    D = np.sqrt((Qxx-Qyy)**2 + 4*Qxy**2)
    eigvec1 = np.matrix([[2*Qxy],[Qyy-Qxx+D]])
    eigvec1 /= np.sqrt(np.transpose(eigvec1)*eigvec1)
    eigvec2 = -np.matrix([[2*Qxy],[Qyy-Qxx-D]])
    eigvec2 /= np.sqrt(np.transpose(eigvec2)*eigvec2)
    R = np.matrix([[eigvec1[0,0], eigvec2[0,0]],
                   [eigvec1[1,0], eigvec2[1,0]]])
    R = np.matrix([[ 0, 1],
                   [-1, 0]])*R
    print 'R', arg_set[arg_num], '\n', R
    Q = np.matrix([[Qxx, Qxy],
                   [Qxy, Qyy]])
    print 'Q\n', Q
    print 'diagonalized Q\n', np.transpose(R)*Q*R

    # Now let's rotate our coordinate system.
    Xrot = R[0,0]*X + R[1,0]*Y + X_position[arg_num]
    Yrot = R[0,1]*X + R[1,1]*Y + Y_position[arg_num]

    cdata = np.array([[0  ,1,1,1],
                      [.02 ,1,1,1],
                      [.25,0.2,.2,1],
                      [.35 ,0,.8,0],
                      [.5 ,1,1,0],
                      [.7 ,1,0,0],
                      [1  ,0,0,0]])
    cdict = {'red':   [], 'green': [], 'blue':  []}
    for i in range(cdata.shape[0]):
        cdict['red']   += [(cdata[i, 0], cdata[i, 1], cdata[i, 1])]
        cdict['green'] += [(cdata[i, 0], cdata[i, 2], cdata[i, 2])]
        cdict['blue']  += [(cdata[i, 0], cdata[i, 3], cdata[i, 3])]
    cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)
    time_max = np.max(data[arg_num])

    levels = np.arange(0,time_max+1.0,1)
    levels = np.arange(0, 41.0, 1)
    if 'p/' in arg_set[arg_num]:
        levels = np.arange(0, 81.0, 1)
    cs = ax.contourf(Xrot, Yrot, data[arg_num], cmap=cmap, origin='lower',levels=levels)

    if arg_num in [1,8]:
        if arg_num == 1:
            axColor = plt.axes([in_viewport_x(color_bar_x), in_viewport_y(row_0y),
                                0.01, in_viewport_y(row_1y-row_0y)])
        else:
            axColor = plt.axes([in_viewport_x(pill_bar_x), in_viewport_y(row_0y),
                                0.01, in_viewport_y(row_1y-row_0y)])
        plt.colorbar(cs, cax=axColor)
    ax.bar(.2,6,.2, .2)

    ax.add_artist(matplotlib.patches.Rectangle((col_1x-2.5,row_my), 5, 0.1, facecolor='b'))
    ax.text(col_1x, row_my - 0.3, r'$5\mu$',
            verticalalignment='top', horizontalalignment='center')

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


    x_vals -= meanX
    y_vals -= meanY
    x_vals,y_vals = R[0,0]*x_vals + R[1,0]*y_vals + X_position[arg_num], R[0,1]*x_vals + R[1,1]*y_vals + Y_position[arg_num]

    for i in range(len(x_vals)-1):
        ax.annotate('',xy=(x_vals[i+1],y_vals[i+1]),xytext=(x_vals[i],y_vals[i]),
                     fontsize=11,
                     arrowprops=dict(arrowstyle='->', shrinkA=0, shrinkB=0, linewidth=.5)) # ,color='black',shrink=0.01, width=.3, headwidth=5.
arrow_file.close()


mannik_micron = 5.0/1012
m1=mpimg.imread('mannik-1.png')
m1x = -np.arange(-m1.shape[0]/2, m1.shape[0]/2)*mannik_micron
m1y = np.arange(-m1.shape[1]/2, m1.shape[1]/2)*mannik_micron
M1x, M1y = np.meshgrid(m1y, m1x)
R1 = np.matrix([[ 0.74914988, -0.66240052],
                [ 0.66240052,  0.74914988]])
M1xrot = R1[0,0]*M1x + R1[1,0]*M1y + x_position_m1
M1yrot = R1[0,1]*M1x + R1[1,1]*M1y + y_position_m1

ax.contourf(M1xrot, M1yrot, m1[:,:,1], levels=[0, 0.9], colors=('r', 'w'))
ax.contourf(M1xrot, M1yrot, m1[:,:,0], levels=[0, 0.9], colors=('k', 'w'))


m2=mpimg.imread('mannik-2.png')
m2x = -np.arange(-m2.shape[0]/2, m2.shape[0]/2)*mannik_micron
m2y = np.arange(-m2.shape[1]/2, m2.shape[1]/2)*mannik_micron
M2x, M2y = np.meshgrid(m2y, m2x)
R2 = np.matrix([[ 0.99256182, -0.12174166],
                [ 0.12174166,  0.99256182]])

M2xrot = R2[0,0]*M2x + R2[1,0]*M2y + x_position_m2
M2yrot = R2[0,1]*M2x + R2[1,1]*M2y + y_position_m2

ax.contourf(M2xrot, M2yrot, m2[:,:,1], levels=[0, 0.9], colors=('r', 'w'))
ax.contourf(M2xrot, M2yrot, m2[:,:,0], levels=[0, 0.9], colors=('k', 'w'))


ax.text(left_annotate_x, row_0y, 'deterministic',
        verticalalignment='center', horizontalalignment='right')
ax.text(left_annotate_x, row_1y, 'stochastic',
        verticalalignment='center', horizontalalignment='right')
ax.text(left_annotate_x, row_my, 'experiment',
        verticalalignment='center', horizontalalignment='right')

ax.text(col_0x, bottom_annotate_y, u'shape A',
        verticalalignment='center', horizontalalignment='center')
ax.text(col_1x, bottom_annotate_y, 'stadium A',
        verticalalignment='center', horizontalalignment='center')
ax.text(col_2x, bottom_annotate_y, u'shape B',
        verticalalignment='center', horizontalalignment='center')
ax.text(col_3x, bottom_annotate_y, 'stadium B',
        verticalalignment='center', horizontalalignment='center')
ax.text(col_4x, bottom_annotate_y, 'natural pill',
        verticalalignment='center', horizontalalignment='center')

ax.set_xlim(viewport_xmin,viewport_xmax)
ax.set_ylim(viewport_ymin,viewport_ymax)

plt.savefig('./paper/plot-ave.pdf', facecolor='white')

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

