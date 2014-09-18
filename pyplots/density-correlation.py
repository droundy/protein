from __future__ import division
import sys
import os.path
import numpy as np
import matplotlib
if "show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab
import re


print 'For this script, the ninth argument should be the end time of the data set'
print 'but you can write use no ninth argument if want to just use the full data set'
constant_start_time = 100.0 #this won't change, will always cut off first 100 seconds of data
end = 0

def ignoreme(value):
    return 0.0

protein_name = "NflD"

difD = 2.5
dx = .05
time_step = .1*dx*dx/difD;#sec
print_denominator = 1000;
dt = time_step*print_denominator

longest_data_len = 0


def readbox(name):
    #notify_reading_file(data_file)
    data = np.loadtxt(data_file, converters = {0:ignoreme, 1:ignoreme})
    global end
    end = (len(data[0])*dt)-10 #a bit less than full data set by default
    if len(sys.argv) > 9 and end > float(sys.argv[9]):
        end  = float(sys.argv[9])
    data_stop_index = int(end/dt)
    data_start_index = int(constant_start_time/dt)
    print 'start time = ',constant_start_time,', end_time = ',end
    print 'this is the time covered by the entire data set that youre using'
    data = data[:,(data_start_index+2):(data_stop_index+2)]
    shortened_data = np.zeros((3, len(data[0,:])))
    print name
    first_column_data = np.genfromtxt(name, dtype='str',usecols=(0,))
    row_num = 0
    while first_column_data[row_num] != protein_name:
        #print first_column_data[row_num]
        row_num += 1
    shortened_data = np.copy(data[row_num:row_num+3,:])
    return shortened_data[:,:]



def mean_density(data,sec):
    total = 0
    for i in data[sec,:]:
        total += i
    return float(total)/float(len(data[sec,:]))



data_files = [10,0,1,2,3,4,5]
if sys.argv[1] == "stad" and sys.argv[3] == "2.92":
    data_files = [0,1,2,3,4,5]
taus = np.arange(0,1000,dt)

total_times_full = np.arange(0,1000,dt)
total_times_full_short = np.arange(0,1000,dt)
total_times_exact = np.arange(0,1000,dt)
total_times_exact_short = np.arange(0,1000,dt)

corrs_full  = np.zeros_like(taus)
corrs_full_short  = np.zeros_like(taus)
corrs_exact = np.zeros_like(taus)
corrs_exact_short = np.zeros_like(taus)

total_calcs = len(data_files)*4
num_calcs_so_far = 0

for sim_type in ["full_array","exact"]:
    for data_length in ["short","long"]:
        for data_file_num in data_files:
            print '\nwe are %.0f%% done!' % (num_calcs_so_far*100.0/total_calcs),\
                'were on data_file #',data_file_num,' ',sim_type,' ',data_length
            num_calcs_so_far += 1
            job_string = "data/shape-%s/%s-%s-%s-%s-%s-%s/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                               sys.argv[5],sys.argv[6],sim_type)
            p = re.compile('[.]')
            job_string = p.sub('_',job_string)
            data_file = ''
            if data_file_num == 10:
                data_file = job_string + 'box-plot.dat'
            else:
                data_file = '../new-protein-%d/'%(data_file_num) + job_string + 'box-plot.dat'
            if not os.path.isfile(data_file):
                print '\n',data_file, ' doesnt exist so we cant load it'
                exit(0)
            data = readbox(data_file)
            if len(data[0,:]) > longest_data_len:
                longest_data_len = len(data[0,:])
            if data_length == "short":
                data = data[:,:len(data[0,:])/2]
            for i in range( min( len(data[0,:])-1, len(taus)-1 ) ):
                integral = 0.0
                total_time = dt*(float(len(data[0,:])-i))
                mean_left = mean_density(data,0)
                mean_right = mean_density(data,2)
                integral = np.dot( data[0,:(len(data[0,:])-i)]-mean_left, data[2,i:]-mean_right)
                if sim_type == "full_array":
                    if data_length == "long":
                        corrs_full[i] += integral
                        total_times_full += total_time
                    elif data_length == "short":
                        corrs_full_short[i] += integral
                        total_times_full_short += total_time
                elif sim_type == "exact":
                    if data_length == "long":
                        corrs_exact[i] += integral
                        total_times_exact += total_time
                    elif data_length == "short":
                        corrs_exact_short[i] += integral
                        total_times_exact_short += total_time

corrs_exact /= total_times_exact
corrs_full /= total_times_full
corrs_exact_short /= total_times_exact_short
corrs_full_short /= total_times_full_short

print 'longest_data_len',longest_data_len
##for the full:
job_string = "data/shape-%s/%s-%s-%s-%s-%s-full_array/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                   sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
printout_file_name = job_string + 'fast-correlation-right-left.dat'
print 'printing out to ',printout_file_name
printout_file = open(printout_file_name,'w')

#following: first line of data file will tell how much data is being used by full box data array
for tau in taus:
    printout_file.write('%g '%tau)
printout_file.write('\n')

for c in corrs_full[:len(corrs_full_short)]:
    printout_file.write('%g '%c)
printout_file.write('\n')

for c in corrs_full_short:
    printout_file.write('%g '%c)

printout_file.close()

plt.figure()
plt.plot(taus,corrs_full)
plt.xlim(0,longest_data_len*dt)
plt.show()


##for the exact:
job_string = "data/shape-%s/%s-%s-%s-%s-%s-exact/" % (sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],
                                                   sys.argv[5],sys.argv[6])
p = re.compile('[.]')
job_string = p.sub('_',job_string)
printout_file_name = job_string + 'fast-correlation-right-left.dat'
print 'printing out to ',printout_file_name
printout_file = open(printout_file_name,'w')

#following: first line of data file will tell how much data is being used by full box data array
for tau in taus:
    printout_file.write('%g '%tau)
printout_file.write('\n')

for c in corrs_exact[:len(corrs_full_short)]:
    printout_file.write('%g '%c)
printout_file.write('\n')

for c in corrs_exact_short:
    printout_file.write('%g '%c)

printout_file.close()
