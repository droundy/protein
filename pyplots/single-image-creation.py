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
    print os.path.isdir(directory)
    if (not os.path.isdir(directory)):
        os.makedirs(directory)



def gaussian_smear(data,wavelength):
    new = np.zeros_like(data)
    N_A = 1.3
    dx = 0.05 # microns per grid spacing
    sigma = .21*wavelength/N_A #n_sin_theta can reach 1.4 to 1.6 in modern optics according to wikipedia
    dis = int(3*sigma/dx) #for now
    for num in range(new.shape[0]):
        print 'num ',num
        norm = 0
        for i in np.arange(-dis,dis,1):
            xstart = 0
            xstop = len(new[num,:,0])
            if i > 0:
                xstart += i
            else:
                xstop += i
            for j in np.arange(-dis,dis,1):
                ystart = 0
                ystop = len(new[num,0,:])
                if j > 0:
                    ystart += j
                else:
                    ystop += j
                weight = math.exp(-(i*i+j*j)*dx*dx/2.0/sigma/sigma )
                new[num,xstart:xstop,ystart:ystop] += data[num,xstart-i:xstop-i,ystart-j:ystop-j] * weight
                norm += weight
        new[num,:,:] /= norm # normalize!
    return new

#computes the global maximum over a set of two dimensional arrays (stored as files)

for i in range(len(proteinList)):
    proteins[i] = load.data(proteinList[i],sim_type,start_time,end_time)

times  = np.arange(float(start_time),float(end_time),dump_time_step)

# proteinList = ["nADP","nATP","ND","NDE","nE",]
dx = 0.05
unit_conversions = [dx,dx,1.0/dx**2,1.0/dx**2,dx]

for i in range(len(proteins)):
    print proteinList[i]
    smeared_data = gaussian_smear(proteins[i].dataset,.650) #this is in microns. green light at 509nm,red at 650nm
    for j in range(len(times)):
        image_data_file = '.' + job_string +str(proteinList[i])+'/images/single-'+str(times[j])+'.dat'
        #image_data_file = '.' + job_string +str(proteinList[i])+'/images/real-gauss-single-'+str(times[j])+'.dat'
        p_file = open(image_data_file,'w')
        p_file.close()
        p_file = open(image_data_file,'a')
        for y in range(smeared_data.shape[1]):
            for x in range(smeared_data.shape[2]):
                p_file.write('%g '%(smeared_data[j,y,x]*unit_conversions[i]))
            p_file.write('\n')
        p_file.close()

