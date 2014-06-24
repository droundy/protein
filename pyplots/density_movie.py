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


#computes the global maximum over a set of two dimensional arrays (stored as files)
def timemax(data):
    spatial_maxes = np.zeros(data.shape[0])
    for i in range(len(spatial_maxes)):
        spatial_maxes[i] = np.max(data[i])
        print spatial_maxes[i]
    return np.max(spatial_maxes)

#computes the global minimum over a set of two dimensional arrays (stored as files)
def timemin(data):
    spatial_mins = np.zeros(data.shape[0])
    for i in range(len(spatial_mins)):
        spatial_mins[i] = np.min(data[i])
        print spatial_mins[i]
    return np.min(spatial_mins)


def gaussian_smear(data,wavelength):
    new = np.zeros_like(data)
    print 'shape '+ str(new.shape)
    print 'len [0,0] ' +str(len(new[0,0]))
    print 'shape [0] '+ str(new[0].shape)
    print 'shape [0,0] '+ str(new[0,0].shape)
    n_sin_theta = 1.5
    sigma = 1.5*wavelength/2.0/n_sin_theta #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    dis = int(2*sigma/0.05) #for now
    for num in range(new.shape[0]):
        for x in range(new.shape[1]):
            for y in range(new.shape[2]):
                for i in np.arange(-dis,dis,1):
                    for j in np.arange(-dis,dis,1):
                        if (x+i >= 0 and x+i < new.shape[1]-1 and y+j >= 0 and y+j < new.shape[2]-1):
                            new[num,x+i,y+j] += data[num,x,y]*math.exp( -(abs(i)+abs(j))*.05/sigma/sigma )
    return new

#function to carry out the animation generation
def contourplt(protein,video_number):
    #get extrema values for color bar (global extrema in time)
    ###############The max vals here are really messy still need to do it properly for smeared data
    smeared_data = gaussian_smear(protein.dataset,.5) #green light at 500nm
    maxval = timemax(smeared_data)
    minval = timemin(smeared_data)

    plt.figure(1)

    #shell command to clean up any previous .png's, just in case (perhaps a process was cancelled midway)
    os.system("rm -f ./data/shape-"+f_shape+"/plots/"+load.debug_str+load.hires_str+load.slice_str+"tmp_*" \
    +"-"+str(protein.protein)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
    +"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".png")

    Z, Y = np.meshgrid(np.arange(0,(len(protein.dataset[1][0])-.9)*0.05,0.05), np.arange(0,len(protein.dataset[1])*0.05,0.05))
    #generate a sequence of .png's for each file (printed time step). these will be used to create a gif.
    for k in range(len(protein.dataset)): #fig.dpi method

        #page = protein.dataset[k]
        print './data/shape-'+f_shape+'/plots/'+load.debug_str+load.hires_str+load.slice_str+'tmp_'+str(k)+'-'+str(protein.protein)+ \
                        '-'+f_shape+'-'+f_param1+'-'+f_param2+'-'+f_param3+'-'+f_param4+'-'+dens_factor+'-'+sim_type+'.png'
        plt.clf()
        plt.axes().set_aspect('equal', 'datalim')
        CS = plt.contourf(Z, Y, smeared_data[k], cmap=plt.cm.jet,origin='lower',levels=np.arange(minval,maxval+10,1))
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
