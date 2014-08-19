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
start_time = float(f_param6)
end_time = float(f_param7)

job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (load.f_shape,load.f_param1,load.f_param2,
                                                   load.f_param3,load.f_param4,load.f_param5,sim_type)
p = re.compile('[.]')
job_string = p.sub('_',job_string)

proteinList = ["nADP","nATP","ND","NDE","nE",]
proteins = [0]*len(proteinList)

for i in range(len(proteins)):
    directory = '.' + job_string + str(proteinList[i]) + '/images'
    print directory
    if (not os.path.isdir(directory)):
        os.makedirs(directory)

# n_sin_theta = 1.5
# wavelength = .6
# sigma = wavelength/2.0/n_sin_theta #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
# dis = int(3.5*sigma/0.05) #for now
# print "wavelength",wavelength
# print "sigma",sigma
# print "dis",dis

# for i in np.arange(-dis,dis,1):
#     for j in np.arange(-dis,dis,1):
#      #   if 10.0*math.exp( -(i*i+j*j)*.05*.05/2.0/sigma/sigma ) > 7.0:
#         print "dis ",math.sqrt((i*i+j*j)),"val ",10.0*math.exp( -(i*i+j*j)*.05*.05/2.0/sigma/sigma )
# print "dis",dis
# exit(0)

def gaussian_smear(data,wavelength):
    new = np.zeros_like(data)
    print "here ",new.shape[0]
    n_sin_theta = 10.5
    #sigma below is tha Abbe resolution.  n is the diffraction of the medium, sin_theta the aperture angle
    #So far I've wiki'd the guassian form of the Airy disk for this
    sigma = wavelength/2.0/n_sin_theta #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    dis = int(3.5*sigma/0.05) #for now
    for num in range(new.shape[0]):
        print "num ",num," of ",new.shape[0]-1
        for x in range(new.shape[1]):
            #print x," of ",new.shape[1]-1
            for y in range(new.shape[2]):
                for i in np.arange(-dis,dis,1):
                    for j in np.arange(-dis,dis,1):
                        if (x+i >= 0 and x+i < new.shape[1]-1 and y+j >= 0 and y+j < new.shape[2]-1
                            and math.sqrt(i*i+j*j) <= dis ):
                            new[num,x+i,y+j] += data[num,x,y]*math.exp( -(i*i+j*j)*.05*.05/2.0/sigma/sigma )
    return new


#computes the global maximum over a set of two dimensional arrays (stored as files)

for i in range(len(proteinList)):
    proteins[i] = load.data(proteinList[i],sim_type,start_time,end_time)

times  = np.arange(float(start_time),float(end_time),dump_time_step)

for i in range(len(proteins)):
    print proteinList[i]
    smeared_data = gaussian_smear(proteins[i].dataset,.6) #this is in microns green light at 500nm,
    for j in range(len(times)):
        image_data_file = '.' + job_string +str(proteinList[i])+'/images/single-'+str(times[j])+'.dat'
        p_file = open(image_data_file,'w')
        p_file.close()
        p_file = open(image_data_file,'a')
        for y in range(smeared_data.shape[1]):
            for x in range(smeared_data.shape[2]):
                p_file.write('%g '%smeared_data[j,y,x])
            p_file.write('\n')
        p_file.close()

