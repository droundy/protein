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

#create data objects (see file_loader.py)

protein_name = "NflD"

dump_time_step = 0.5
input_start_time = float(f_param6)
input_end_time = float(f_param7)
total_number_of_files = 0

if (input_start_time > input_end_time):
    print "Your start time is greater than your end time"
    exit(1)


job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (load.f_shape,f_param1,f_param2,
                                                    f_param3,f_param4,dens_factor,sim_type)

p = re.compile('[.]')
job_string = p.sub('_',job_string)
dir_name = job_string + 'ave-time'

if not os.path.exists(dir_name):
    print "making directory "+dir_name+" because doesnt exist"
    os.makedirs(dir_name)

for i in range(0,20000):
    fname = job_string + protein_name + '/movie-frame-%05d.dat'%(i)

    total_number_of_files += 1
    if (not os.path.isfile(fname)):
        print fname + ' is not a file'
        break
if (input_end_time >= total_number_of_files*dump_time_step):
    print "For ", fname
    print input_end_time
    print total_number_of_files*dump_time_step
    print "This end_time is too great, there are only enough files to support a end_time less than ", total_number_of_files*dump_time_step
    exit(1)


def gaussian_smear(data,wavelength,protein):
    new = np.zeros_like(data[0])
    unsmeared = np.zeros_like(data[0])
    N_A = 1.3
    dx = 0.05 # microns per grid spacing
    sigma = .21*wavelength/N_A #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    print "sigma ",sigma
    dis = int(3*sigma/dx) #for now
    arrow_file = job_string +'ave-time/ave-time-arrow-'+str(int(input_start_time))+'-'+str(protein)+'.dat'
    print arrow_file
    p_file = open(arrow_file,'w')
    last_max_x = 0
    last_max_y = 0
    for num in range(data.shape[0]):
        max_data = np.zeros_like(data[0])
        norm = 0
        for i in np.arange(-dis,dis,1):
            xstart = 0
            xstop = len(new[:,0])
            if i > 0:
                xstart += i
            else:
                xstop += i
            for j in np.arange(-dis,dis,1):
                ystart = 0
                ystop = len(new[0,:])
                if j > 0:
                    ystart += j
                else:
                    ystop += j
                weight = math.exp(-(i*i+j*j)*dx*dx/2.0/sigma/sigma )
                max_data[xstart:xstop,ystart:ystop] += data[num,xstart-i:xstop-i,ystart-j:ystop-j] * weight
                norm += weight
        max_data /= norm # normalize!
        new += max_data
        unsmeared += data[num,:,:]

        maxima = 0
        max_x = 0
        max_y = 0
        for x in range(max_data.shape[0]):
            for y in range(max_data.shape[1]):
                if max_data[x,y] > maxima:
                    maxima = max_data[x,y]
                    max_x = x
                    max_y = y
        p_file.write('%g %g %g %g\n'%(input_start_time+num*dump_time_step,maxima,max_x,max_y))
        if (num%20 == 0 and num > 1):
            contour_values = job_string +'ave-time/contour-values-' + str(protein) +'-'+ str(int(input_start_time))+'-' \
                +str(int(input_start_time+num*dump_time_step))+'.dat'
            unsmeared_values = job_string +'ave-time/unsmeared-values-' + str(protein) +'-'+ str(int(input_start_time))+'-' \
                +str(int(input_start_time+num*dump_time_step))+'.dat'
            print contour_values
            c_file = open(contour_values,'w')
            unsmeared_file = open(unsmeared_values,'w')
            for x in range(new.shape[0]):
                for y in range(new.shape[1]):
                    c_file.write("%g "%(new[x,y]/num/dx**2))
                    unsmeared_file.write("%g "%(unsmeared[x,y]/num/dx**2))
                c_file.write('\n')
                unsmeared_file.write('\n')
            c_file.close()
            unsmeared_file.close()
    p_file.close()
    return new/data.shape[0]


data = load.data(protein=protein_name, sim_type=sim_type,start_time = input_start_time, end_time = input_end_time)
smeared_data = gaussian_smear(data.dataset,.509,data.protein) #this is in microns green light at 500nm,
