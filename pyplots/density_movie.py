from __future__ import division
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import time
import file_loader as load


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

protein_name = "NflD"

dump_time_step = 0.5
input_start_time = float(f_param6)
input_end_time = float(f_param7)
total_number_of_files = 0

if (input_start_time > input_end_time):
    print "Your start time is greater than your end time"
    exit(1)

for i in range(0,20000):
    fname = "./data/shape-%s/%s-%s-%s-%s-%s-%s-%s-%03d-%s.dat"%(f_shape,protein_name,f_shape,f_param1,
                                                                f_param2,f_param3,f_param4,dens_factor,i,sim_type)
    total_number_of_files += 1
    if (not os.path.isfile(fname)):
        break
if (input_end_time >= total_number_of_files*dump_time_step):
    print "This end_time is too great, there are only enough files to support a end_time less than ", total_number_of_files*dump_time_step
    exit(1)


time_left = input_end_time - input_start_time
video_limit = 700 #time of each gif created
video_list = []#list of video_limit long movies
video_number = 0

while (time_left > 0):
    next_end_time = input_start_time + (video_number+1)*video_limit
    if next_end_time > input_end_time:
        next_end_time = input_end_time
    print "next_end_time = ",next_end_time
    video_list = video_list + [load.data(protein=protein_name, sim_type=sim_type,start_time = input_start_time + video_number*video_limit,
                                         end_time = next_end_time)]
    time_left -= video_limit
    video_number += 1

print "video number = ",video_number
print "time left = ",time_left
print "size video_list = ",len(video_list)

#  NflE = load.data(protein="NflE")
# nATP = load.data(protein="nATP")
# nE = load.data(protein="nE")
# nADP = load.data(protein="nADP")
# ND = load.data(protein="ND")
# NDE = load.data(protein="NDE")

#computes the maximum of a two dimensional array - REPLACE with method
def maxnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = max(page[i])
    maxval = max(Z)
    return maxval

#computes the minimum of a two dimensional array - REPLACE with method
def minnum(page):
    Z = [0. for i in range(page.shape[0])]
    for i in range(page.shape[0]):
        Z[i] = min(page[i])
        if Z[i] == 0:
            Z[i] = 50000000 #eh
    minval = min(Z)
    return minval

#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = maxnum(protein.dataset[i+1])
    maxval = max(Z)
    return maxval

#computes the global minimum over a set of two dimensional arrays (stored as files)
def timemin(protein):
    Z = [0. for i in range(protein.tsteps)]
    for i in (range(protein.tsteps-1)):
        Z[i+1] = minnum(protein.dataset[i+1])
    minval = min(Z)
    return minval

#function to carry out the animation generation
def contourplt(protein,video_number):
    #get extrema values for color bar (global extrema in time)
    maxval = timemax(protein)/2.0
    minval = timemin(protein)
    plt.figure(1)

    #shell command to clean up any previous .png's, just in case (perhaps a process was cancelled midway)
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+load.debug_str+load.hires_str+load.slice_str+"tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".png")

    Z, Y = np.meshgrid(np.arange(0,(len(protein.dataset[1][0])-.9)*0.05,0.05), np.arange(0,len(protein.dataset[1])*0.05,0.05))
    #generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    for k in range(len(protein.dataset)): #fig.dpi method
        page = protein.dataset[k]
        print './data/shape-'+f_shape+'/plots/'+load.debug_str+load.hires_str+load.slice_str+'tmp_'+str(k)+'-'+str(protein.protein)+ \
                        '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+dens_factor+'-'+sim_type+'.png'
        plt.clf()
        plt.axes().set_aspect('equal', 'datalim')
        CS = plt.contourf(Z, Y, page, cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval+10,1))
        cbar = plt.colorbar(CS)
        plt.clim(minval,maxval)
        plt.title(str(protein.protein)+" volume density at time: "+str(k*dump_time_step)+" sec")
        #plt.clabel(CS, fontsize=9, inline=1)
        plt.xlabel("Z position")
        plt.ylabel("Y position")
        plt.savefig('./data/shape-'+f_shape+'/plots/'+load.debug_str+load.hires_str+load.slice_str+'tmp_%06d'%k+'-'+str(protein.protein)+ \
                        '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+dens_factor+'-'+sim_type+'.png',dpi=50)

    #shell command to convert all of the recently generated .png's to a single .gif using convert utility
    os.system("convert -delay 10 ./data/shape-"+f_shape+"/plots/"+load.debug_str+load.hires_str+load.slice_str+"tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".png ./data/shape-"+f_shape \
    +"/plots/"+load.debug_str+load.hires_str+load.slice_str+"density_movie-"+str(protein.protein)+"-"+str(video_number)+"-"+f_shape \
    +"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".gif")

    # shell command to clean up the .png's
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+load.debug_str+load.hires_str+load.slice_str+"tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".png")
    return 0

# contourplt(NDE)
for j in range(video_number):
    contourplt(video_list[j],j)

# contourplt(NflE)
# contourplt(nATP)
# contourplt(nE)
# contourplt(nADP)
# contourplt(ND)
