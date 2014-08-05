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

#create data objects (see file_loader.py)



  # \includegraphics[width=\columnwidth]{../data/shape-randst/plots/arrow-1500-95-stoch}
  # \includegraphics[width=\columnwidth]{../data/shape-randst/plots/arrow-1850-95-stoch}
  # \includegraphics[width=\columnwidth]{../data/shape-randst/plots/arrow-2500-95-stoch}



dx =0.05
dump_time_step = 0.5
protein_name = "NflD"
# start_time = float(f_param6)
# input_end_time = float(f_param7)
#end_time = int(input_end_time - input_end_time%10)

arg_set = ["randst-0.25-15.00-15.00-95.00-15.00",
           "randst-0.25-18.50-18.50-95.00-15.00",
           "randst-0.25-25.00-25.00-95.00-15.00"]


bound_times = ["300","1500","300","1100","300","920"]

arrow_files = []
contour_values = []
for i in range(len(arg_set)):
    arrow_files += ['./data/shape-randst/plots/ave-time/maxima-arrowNflD-'+arg_set[i]+'-full_array.dat']
    contour_values += ['./data/shape-randst/plots/ave-time/contour-values-'+bound_times[i*2]+'-' \
        +bound_times[i*2+1]+'-NflD-'+arg_set[i]+'-full_array.dat']


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
hts += [hts[1]]

gs = gridspec.GridSpec(len(arg_set)+1,1,
                       height_ratios=hts,
                       )
for arg_num in range(len(arg_set)+1):
#    plt.subplot2grid((len(arg_set)+1,1),(arg_num,0),aspect='equal')
    plt.subplot(gs[arg_num,0],aspect='equal')
    #gs_in = gridspec.GridSpec(arg_num,0)
    #gs_in.update(left=.3+arg_num*.1,right=.8)
    if arg_num == len(arg_set):
            m1=mpimg.imread('mannik-1.png')
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
    #plt.set_aspect('equal')
        # rect = [.2, .2+arg_num*.3, .6, .3]
        # plt.axes(rect)
        c_data = np.loadtxt(contour_values[arg_num])
        a_data = np.loadtxt(arrow_files[arg_num])
        start_time = float(bound_times[arg_num*2])
        end_time = float(bound_times[arg_num*2+1])




        # plt.subplots_adjust(hspace=0.1,left=(.5-half_width),right=(.5+half_width))
        # plt.subplots_adjust(hspace=0.1,left=.4,right=.5)
        # left  = 0.125  # the left side of the subplots of the figure
        # right = 0.9    # the right side of the subplots of the figure

        a_data = a_data[int((start_time-300)/dump_time_step):int((end_time-start_time)/dump_time_step)]

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







#plt.title("Time-averaged density from "+str(start_time)+" to "+str(input_end_time))
#plt.xlabel("Z position")
#plt.ylabel("Y position")



# for i in range(len(arg_set)):
#     cell_membrane_initial = np.loadtxt("./data/shape-randst/membrane_files/sections-%s.dat" % arg_set[i])

#     #http://stackoverflow.com/questions/8486294/how-to-add-an-extra-column-to-an-numpy-array
#     #this is algorithm want
#     Ny = len(cell_membrane_initial[:,0])
#     Nz = len(cell_membrane_initial[0,:])

#     temp_array = np.zeros((Ny+1,Nz+1))
#     temp_array[:-1,:-1] = cell_membrane_initial
#     cell_membrane = temp_array

#     row_of_zeros = np.zeros(len(cell_membrane[:,0]))
#     np.vstack((cell_membrane[:,0],row_of_zeros))
#     np.vstack((row_of_zeros,cell_membrane[:,0]))

#     row_of_zeros = np.zeros(len(cell_membrane[0,:]))
#     np.vstack((cell_membrane[0,:],row_of_zeros))
#     np.vstack((row_of_zeros,cell_membrane[0,:]))

#     arg_list = list(arg_set[i])
#     for j in range(len(arg_list)-1):
#         if arg_list[j] == "9" and arg_list[j+1] == "4":
#             np.transpose(cell_membrane)
#             print "There is a 94 in the argument so were transposing!"

#     print './data/shape-randst/arrow-plot-NflD-%s.dat' % arg_set[i]

#     tails = np.loadtxt('./data/shape-randst/arrow-plot-NflD-%s.dat' % arg_set[i])
#     if len(tails)==0:
#         print "The tails in this file has no tails!!"
#         #exit(1)
#     for k in range(len(tails[:,0])-1,0,-1):
#         if tails[k,3] <0:
#             tails = np.delete(tails,k,0)
#         if tails[k,2] < tails[k-1,2]+3:
#             if tails[k,3] > tails[k-1,3]:
#                 tails = np.delete(tails,k-1,0)
#             else:
#                 tails = np.delete(tails,k,0)
#     tails = tails*dx
#     tails = tails[9:,:] # cut out the first few arrows
#     ax = plt.subplot(1,1,1)
#     cell_membrane[cell_membrane>0] = 1
#     cell_x = np.linspace(0,dx*len(cell_membrane[0,:]),len(cell_membrane[0,:]))
#     cell_y = np.linspace(0,dx*len(cell_membrane[:,0]),len(cell_membrane[:,0]))
#     cell_x, cell_y = np.meshgrid(cell_x, cell_y)
#     centerx = cell_x[cell_membrane == 1].mean()
#     centery = cell_y[cell_membrane == 1].mean()
#     originx = (1 - (i//2))*8 - centerx -8
#     originy = (i%2)*4.5 + (i%2)**2*0.6 - centery + 1
#     ax.set_aspect('equal', 'datalim')
#     font=FontProperties()
#     font.set_family('serif')
#     for j in range(len(tails)-1):
#         radial_length = np.sqrt((tails[j][0]-dx*Ny/2.0)**2+(tails[j][1]-dx*Nz/2.0)**2)
#         dir_z = (tails[j][0]-dx*Nz/2.0)/radial_length
#         dir_y = (tails[j][1]-dx*Ny/2.0)/radial_length
#         number_zpos = tails[j][0]+.18*dir_z
#         number_ypos = tails[j][1]+.18*dir_y
#         ax.annotate('',xy=(originx+tails[j+1][1],originy+tails[j+1][0]),xytext=(originx+tails[j][1],originy+tails[j][0]),
#                     fontsize=7,
#                     arrowprops=dict(color='red',shrink=0.01, width=.3, headwidth=5.))
#   #      ax.text(number_zpos,number_ypos,"%d"%(j+1),fontsize=30,fontproperties=font,color='red')
#     plt.contour(originx+cell_x,originy+cell_y,cell_membrane, linewidths=2,levels=[.99])
#     #plt.ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
#     # ax.get_yaxis().set_visible(True)
#     ax.add_artist(AnchoredSizeBar(
#             ax.transData,
#             2.13, # length of the bar in the data reference
#             "2.13$\mu$", # label of the bar
#             bbox_to_anchor=(44,40),
#             loc=6, # 'best', # location (lower right)
#             pad=-.4, borderpad=0.25, sep=3,
#             frameon=False
#             ))
#     ax.add_artist(AnchoredSizeBar(
#             ax.transData,
#             2.13, # length of the bar in the data reference
#             "2.13$\mu$", # label of the bar
#             bbox_to_anchor=(220,53),
#             loc=6, # 'best', # location (lower right)
#             pad=-.4, borderpad=0.25, sep=3,
#             frameon=False
#             ))
# # ax.get_xaxis().set_visible(True)
#     #plt.xlim((0,dx*cell_membrane.shape[1]))
#     #plt.ylim((0,dx*cell_membrane.shape[0]))
#     #plt.xlabel("Z grid position")
#     #plt.ylabel("Y grid position")
#     #plt.title("Local temporal maxima, global spatial maxima view of %s"%(protein))
# plt.savefig("data/shape-randst/plots/paper-arrow-plot.pdf")


# plt.show()
