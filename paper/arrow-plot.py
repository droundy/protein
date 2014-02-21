from __future__ import division
import matplotlib
import sys
if "show" not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

#membrane.dat or arrow.dat printed transposed, gotta align them:

dx=0.05

#The following two sets are the potentially relevent Mannik shape
#data.  For the 95 shape, the 11.00 11.00 data is real Mannik size and
#for the 94.00 the 13.00 20.00 data is real Mannik size

arg_set = ["randst-0.25-11.00-11.00-95.00-15.00","randst-0.25-20.00-20.00-95.00-15.00",
           "randst-0.25-21.00-21.00-95.00-15.00","randst-0.25-28.00-28.00-95.00-15.00",
           "randst-0.25-13.00-20.00-94.00-15.00","randst-0.25-18.00-25.00-94.00-15.00",
           "randst-0.25-23.00-33.00-94.00-15.00","randst-0.25-30.00-40.00-94.00-15.00"]

plt.figure()

for i in range(len(arg_set)):
    cell_membrane_initial = np.loadtxt("./data/shape-randst/membrane_files/sections-%s.dat" % arg_set[i])

    #http://stackoverflow.com/questions/8486294/how-to-add-an-extra-column-to-an-numpy-array
    #this is algorithm want
    Ny = len(cell_membrane_initial[:,0])
    Nz = len(cell_membrane_initial[0,:])

    temp_array = np.zeros((Ny+1,Nz+1))
    temp_array[:-1,:-1] = cell_membrane_initial
    cell_membrane = temp_array

    row_of_zeros = np.zeros(len(cell_membrane[:,0]))
    np.vstack((cell_membrane[:,0],row_of_zeros))
    np.vstack((row_of_zeros,cell_membrane[:,0]))

    row_of_zeros = np.zeros(len(cell_membrane[0,:]))
    np.vstack((cell_membrane[0,:],row_of_zeros))
    np.vstack((row_of_zeros,cell_membrane[0,:]))

    arg_list = list(arg_set[i])
    for j in range(len(arg_list)-1):
        if arg_list[j] == "9" and arg_list[j+1] == "4":
            np.transpose(cell_membrane)
            print "There is a 94 in the argument so were transposing!"

    print './data/shape-randst/arrow-plot-NflD-%s.dat' % arg_set[i]

    tails = np.loadtxt('./data/shape-randst/arrow-plot-NflD-%s.dat' % arg_set[i])
    if len(tails)==0:
        print "The tails in this file has no tails!!"
        exit(1)
    print tails
    for k in range(len(tails[:,0])-1,0,-1):
        if tails[k,3] <0:
            tails = np.delete(tails,k,0)
            print k
        else:
            print tails[k,3]
        if tails[k,2] < tails[k-1,2]+3:
            if tails[k,3] > tails[k-1,3]:
                tails = np.delete(tails,k-1,0)
            else:
                tails = np.delete(tails,k,0)
    print tails
    tails = tails*dx
    tails = tails[22:,:] # cut out just the first arrow
    ax = plt.subplot(2,4,i)
    ax.set_aspect('equal', 'datalim')
    font=FontProperties()
    font.set_family('serif')
    for j in range(len(tails)-1):
        radial_length = np.sqrt((tails[j][0]-dx*Ny/2.0)**2+(tails[j][1]-dx*Nz/2.0)**2)
        dir_z = (tails[j][0]-dx*Nz/2.0)/radial_length
        dir_y = (tails[j][1]-dx*Ny/2.0)/radial_length
        number_zpos = tails[j][0]+.18*dir_z
        number_ypos = tails[j][1]+.18*dir_y
        ax.annotate('',xy=(tails[j+1][1],tails[j+1][0]),xytext=(tails[j][1],tails[j][0]),
                    fontsize=7,
                    arrowprops=dict(color='red',shrink=0.01, width=.3, headwidth=5.))
        #ax.text(number_zpos,number_ypos,"%d"%(j+1),fontsize=30,fontproperties=font,color='red')
    cell_membrane[cell_membrane>0] = 1
    cell_x = np.linspace(0,dx*len(cell_membrane[0,:]),len(cell_membrane[0,:]))
    cell_y = np.linspace(0,dx*len(cell_membrane[:,0]),len(cell_membrane[:,0]))
    plt.contour(cell_x,cell_y,cell_membrane, linewidths=2,levels=[.99])
    #plt.ax.quiver(X,Y,U,V,scale_units='xy',angles='xy',scale=1)
    ax.get_yaxis().set_visible(True)
    ax.get_xaxis().set_visible(True)
    #plt.xlim((0,dx*cell_membrane.shape[1]))
    #plt.ylim((0,dx*cell_membrane.shape[0]))
    plt.xlabel("Z grid position")
    plt.ylabel("Y grid position")
    #plt.title("Local temporal maxima, global spatial maxima view of %s"%(protein))
plt.savefig("data/shape-randst/plots/paper-arrow-plot.pdf")


plt.show()



# x = np.arange(0,10,.1)
# y = np.array(x*x)
# y2 = np.array(x*x*x)
# y3 = np.array(x*x*x*x)

# f, axarr = plt.subplots(3,2)
# axarr[0,0].plot(x,)
# axarr[0,1].plot(x,y2)
# axarr[1,0].plot(x,y3)
# axarr[1,1].scatter(x,y)
# axarr[2,0].scatter(x,y)
# axarr[2,1].scatter(x,y)


# plt.show()
