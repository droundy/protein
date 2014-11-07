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
print "hello"

#create data objects (see file_loader.py)

protein_name = "NflD"

dump_time_step = 0.5
start_time = int(float(f_param6))
end_time = int(float(f_param7))
total_number_of_files = 0

if (start_time > end_time):
    print "Your start time is greater than your end time"
    exit(1)

min_file_num = 100000
max_file_num = 0

job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (load.f_shape,f_param1,f_param2,
                                                    f_param3,f_param4,dens_factor,sim_type)
p = re.compile('[.]')
job_string = p.sub('_',job_string)

for i in range(0,20000):
    fname = '.' + job_string + protein_name + '/images/frame-%05d.png'%(i)
    total_number_of_files += 1
    if (os.path.isfile(fname)):
        if i < min_file_num:
            min_file_num = i
        if i > max_file_num:
            max_file_num = i
print "min file number available = ",min_file_num
print "max file number available = ",max_file_num

if (end_time >= total_number_of_files*dump_time_step):
    print "files ", total_number_of_files
    print "For ", fname
    print "The end_time is too great or start_time too early..."
    print "there are only enough files to support a total time of less than ", total_number_of_files*dump_time_step
    exit(1)

num_in_a_video = 2000
number_of_videos = 1 + int( (end_time-start_time)/float(dump_time_step)/float(num_in_a_video) )

print num_in_a_video
print number_of_videos

for video_num in range(number_of_videos):
    start_file = int( start_time/dump_time_step + video_num*num_in_a_video )
    end_file = int( start_time/dump_time_step + (video_num + 1)*num_in_a_video )
#    print '\n\n HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhere\n\n'
    if (video_num + 1 == number_of_videos):
        end_file = int(end_time/float(dump_time_step))
#        print '\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHhello? end_file = %d ',end_file
    file_names = ""
    for i in range(start_file,end_file):
        file_names += '-i .' + job_string + protein_name + '/images/frame-%05d.png '%(i)
    cell_shape = ''
    if load.f_shape == 'p':
        cell_shape = 'pill-'
    elif load.f_shape == 'randst':
        if f_param4 == '95.00':
            cell_shape = 'shape-A-'
        if f_param4 == '94.00':
            cell_shape = 'shape-B-'
    elif load.f_shape == 'stad':
        if f_param2 == '2.35':
            cell_shape = 'stadium-A-'
        if f_param2 == '2.92':
            cell_shape = 'stadium-B-'
    print cell_shape
    os.system('avconv -r 8 '+file_names+ './movies/density-movie-total-MinD-'+cell_shape\
                  +str(int(start_file*dump_time_step))+ '-' +str(int(end_file*dump_time_step))+'.avi')
    # os.system('convert -delay 20 '+file_names+ '.' + job_string + 'plots/density-movie-total-MinD-'+\
    #               +str(int(start_file*dump_time_step))+ '-' +str(int(end_file*dump_time_step))+'.mp4')
    # os.system('ffmpeg -i '+file_names+ '-r 24 ./movies/density-movie-total-MinD-'+cell_shape\
    #               +str(int(start_file*dump_time_step))+ '-' +str(int(end_file*dump_time_step))+'.avi')

    # convert -delay 30 '+file_names+ './movies/density-movie-total-MinD-'+cell_shape\
    #               +str(int(start_file*dump_time_step))+ '-' +str(int(end_file*dump_time_step))+'.avi')


# ./data/shape-"+f_shape+"/plots/"+load.debug_str+load.hires_str \
#                   +load.slice_str+"density_movie-"+str(protein_name)+"-"+str(int(start_file*dump_time_step))+"-" \
#                   +str(int(end_file*dump_time_step))+"-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3 \
#                   +"-"+f_param4+"-"+dens_factor+"-"+sim_type+".mp4")

# contourplt(NDE)
# for j in range(video_number):
#     print "huh"
#     print f_shape
#     print len(f_shape)
#     contourplt(video_list[j],j)

# contourplt(NflE)
# contourplt(nATP)
# contourplt(nE)
# contourplt(nADP)
# contourplt(ND)
