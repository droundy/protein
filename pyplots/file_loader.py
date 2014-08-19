import numpy as np
import sys
import glob
import os
import re

#get the shape an physical parameters of the cell (last argument is density lopsidedness)
f_shape = sys.argv[1]
f_param1 = sys.argv[2]
f_param2 = sys.argv[3]
f_param3 = sys.argv[4]
f_param4 = sys.argv[5]
f_param5 = sys.argv[6]
sim_type = sys.argv[7]

git_add_files = False # set this to True to automagically git add dat
                      # files that we need to do the plots

dx=.05

#for the specific file reading, these are the protein strings
proteinList = ["D_nATP","D_nADP","D_E_NDE","E_nE","D_ND","NflD","NflE"]

#initialize as empty lists incase arguments not present in sys.argv
hires_str=""
slice_str=""
debug_str=""

#change them to print the correct file names if arguments present
if "-hires" in sys.argv:
    dx=.025
    hires_str="hires-"
if "-slice" in sys.argv:
    slice_str="slice-"
if "-debug" in sys.argv:
    debug_str="debug-"

expected_arg_number = 7
for arg in sys.argv:
    if "-" in arg:
        expected_arg_number+=1

if len(sys.argv)>expected_arg_number:
    f_param6 = round(float(sys.argv[8]))
    f_param7 = round(float(sys.argv[9]))

filename_tuple = (f_shape,hires_str,slice_str,f_shape,f_param1,f_param2,f_param3,f_param4,f_param5)



#data object, to create in another plot: proteinname=data(protein="proteinname")
class data(object):
    def __init__(self,protein,sim_type,start_time,end_time):
        self.protein = protein
        self.filenames = get_filenames(protein,sim_type,start_time,end_time)
        self.tsteps = len(self.filenames)
        self.dataset = np.array([np.loadtxt(file) for file in self.filenames])
        # for file in self.filenames:
        #     print 'loading', file
        self.datashape = self.dataset[0].shape
        self.axes = [[i*dx for i in range(self.datashape[1])],[j*dx for j in range(self.datashape[0])]]

#loads the files for creating the data object.
def get_filenames(protein,sim_type,start_time,end_time):
    dump_time_step = 0.5 # seconds
    dat_filenames = []
    if end_time == 0:
        end_time = 10000
    for f_num in np.arange(start_time,end_time,dump_time_step):
        file_num = round(f_num/dump_time_step)
        job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (f_shape,f_param1,f_param2,
                                                            f_param3,f_param4,f_param5,sim_type)
        p = re.compile('[.]')
        job_string = p.sub('_',job_string)
        fname = '.' + job_string + protein + '/movie-frame-%05d.dat'%(file_num)
        print fname
        if os.path.isfile(fname):
            dat_filenames.append(fname)
    if git_add_files:
        for filename in dat_filenames:
            os.system('git add -f %s' % filename)
    if (dat_filenames == []):
        print "File loading error: filename list is empty."
        exit(1)
    else:
        return dat_filenames

#function for easier plot name printing. probably should be renamed itself.
def print_string(plot_name,protein):
    job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (f_shape,sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sim_type)
    p = re.compile('[.]')
    job_string = p.sub('_',job_string)
    filename = ''
    if (protein == ''):
        filename = '.' + job_string + 'plots/' + plot_name + '.pdf'
    else:
        filename = '.' + job_string + 'plots/' + protein + '-' + plot_name + '.pdf'
    print "printing to ",filename
    return filename
