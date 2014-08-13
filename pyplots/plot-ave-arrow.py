from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load
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

#create data objects (see file_loader.py)

dx =0.05
dump_time_step = 0.5
protein_name = "NflD"
start_time = float(f_param6)
input_end_time = float(f_param7)
end_time = int(input_end_time - input_end_time%10)

arrow_file = './data/shape-'+f_shape+'/plots/ave-time/maxima-arrow-'+str(int(start_time))+'-'+protein_name+ \
    '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+dens_factor+'-'+sim_type+'.dat'

contour_values = './data/shape-'+f_shape+'/plots/ave-time/contour-values-'+str(int(start_time))+'-' \
    +str(end_time)+'-'+protein_name+'-'+f_shape+'-'+f_param1+'-'+f_param2+'-' \
    +f_param3+'-'+f_param4+'-'+dens_factor+'-'+sim_type+'.dat'

print arrow_file
print contour_values

c_data = np.loadtxt(contour_values)
print os.path.isfile(arrow_file)
a_data = np.loadtxt(arrow_file)

print a_data
print "this ",a_data[:,0]
print len(a_data[:,0])

print "hello",str((end_time-start_time)*dump_time_step)
a_data = a_data[0:int((end_time-start_time)/dump_time_step)]

print a_data[:,0]
print len(a_data[:,0])

last_time = a_data[0,0]
for i in range(1,len(a_data[:,0])):
    if (a_data[i,0] <= last_time):
        print "The times in the arrow file are not in chronological order!  Something's wrong! Exiting."
        print i
        print last_time
        print a_data[i,0]
        exit(0)


time_max = np.max(c_data)
arrow_cutoff = 3.2*(np.max(a_data[:,1]))/5.0

high_maximas = np.zeros(0)
times = np.zeros(0)
x_vals = np.zeros(0)
y_vals = np.zeros(0)
last_x = 0
last_y = 0
index = 0
print np.max(a_data[:,1])
print "arrow cutoff ",arrow_cutoff
for i in range(len(a_data[:,1])):
    if (a_data[i,1] > arrow_cutoff):
        if ( (a_data[i,2]*dx != last_x or a_data[i,3]*dx != last_y) \
                 and math.sqrt((a_data[i,2]*dx - last_x)**2 + (a_data[i,3]*dx-last_y)**2) > .2):
            x_vals = np.append(x_vals,a_data[i,2]*dx)
            y_vals = np.append(y_vals,a_data[i,3]*dx)
            times = np.append(times,a_data[i,0])
            last_x = a_data[i,2]*dx
            last_y = a_data[i,3]*dx
print x_vals
print y_vals

#don't know why but the stad 5.00 nees to switch this second value to c_data.shape[0]-.9 not c_data.shape[0]
Z, Y = np.meshgrid(np.arange(0,(c_data.shape[1]-.9)*dx,dx),np.arange(0,(c_data.shape[0]-.9)*dx,dx))

print '*******************************************************************'
print c_data.shape
print '*******************************************************************'
zwidth = Z.max() - Z.min()
ywidth = Y.max() - Y.min()
figwidth = 6
barwidth = 0.2*6
plt.figure(figsize=(figwidth,(figwidth - barwidth)*ywidth/zwidth)) # leave room for the colorbar!
plt.clf()

print "start"
print len(Z)
print len(Z[0])
print c_data.shape[1]
print c_data.shape[0]
print len(c_data)
print len(c_data[0])

plt.axes().set_aspect('equal', 'datalim')
cdata = np.array([[0  ,1,1,1],
                  [.1 ,1,1,1],
                  [.25,0.8,.8,1],
                  [.5 ,0,.8,.8],
                  [.7 ,1,1,0],
                  [.9 ,1,0,0],
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
CS = plt.contourf(Z, Y, c_data, cmap=cmap,origin='lower',levels=np.arange(0,time_max+1.0,1))
cbar = plt.colorbar(CS)

for i in range(len(x_vals)-1):
    plt.annotate('%g'%i,xy=(y_vals[i+1],x_vals[i+1]),xytext=(y_vals[i],x_vals[i]),
                 fontsize=9,
                 arrowprops=dict(color='red',shrink=0.01, width=.3, headwidth=5.))
plt.clim(0,time_max)
#plt.axis('off')
save_file_name = './data/shape-'+f_shape+'/plots/plot-time-averaged-arrow-'+str(int(start_time))+'-' \
                +str(end_time)+'-'+protein_name+'-'+f_shape+'-'+str(int(100*float(f_param1)))+'-'+str(int(100*float(f_param2))) \
                +'-'+str(int(100*float(f_param3)))+'-'+str(int(100*float(f_param4)))+'-'+str(int(100*float(dens_factor)))+'-'+sim_type+'.pdf'
print save_file_name
plt.savefig(save_file_name)
