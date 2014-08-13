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
for i in range(0,20000):
    fname = "./data/shape-"+f_shape+"/plots/images/"+load.debug_str+load.hires_str+load.slice_str+"tp_"+str(i) \
        +"-"+str(protein_name)+"-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".png"
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

num_in_a_video = 400
number_of_videos = 1 + int( (end_time-start_time)/float(dump_time_step)/float(num_in_a_video) )

print num_in_a_video
print number_of_videos

for video_num in range(number_of_videos):
    start_file = int( start_time/dump_time_step + video_num*num_in_a_video )
    end_file = int( start_time/dump_time_step + (video_num + 1)*num_in_a_video )
    if (video_num + 1 == number_of_videos):
        end_file = int(end_time/float(dump_time_step))
    file_names = ""
    for i in range(start_file,end_file):
        file_names += "./data/shape-"+f_shape+"/plots/images/"+load.debug_str+load.hires_str \
            +load.slice_str+"tp_%d"%i +"-"+str(protein_name)+"-"+f_shape+"-"+f_param1+"-"+f_param2 \
            +"-"+f_param3+"-"+f_param4+"-"+dens_factor+"-"+sim_type+".png "
    os.system("convert -delay 10 "+file_names+ "./data/shape-"+f_shape+"/plots/"+load.debug_str+load.hires_str \
                  +load.slice_str+"density_movie-"+str(protein_name)+"-"+str(int(start_file*dump_time_step))+"-" \
                  +str(int(end_file*dump_time_step))+"-"+f_shape+"-"+f_param1+"-"+f_param2+"-"+f_param3 \
                  +"-"+f_param4+"-"+dens_factor+"-"+sim_type+".mp4")

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