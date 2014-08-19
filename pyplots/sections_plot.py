from __future__ import division
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import pylab
import file_loader as load
import re

job_string = "/data/shape-%s/%s-%s-%s-%s-%s-%s/" % (load.f_shape,load.f_param1,load.f_param2,
                                                    load.f_param3,load.f_param4,load.f_param5,load.sim_type)
p = re.compile('[.]')
job_string = p.sub('_',job_string)
data = np.loadtxt('.' + job_string + '/sections.dat')

plt.figure()
plt.pcolormesh(data)
plt.savefig('.' + job_string + '/plots/sections.pdf')
