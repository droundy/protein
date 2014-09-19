import numpy as np
import matplotlib

cdata = np.array([[0  ,1,1,1],
                  [.02 ,1,1,1],
                  [.05,0.7,.3,1],
                  [.15,0.2,.2,1],
                  [.3 ,0,.8,0],
                  [.4 ,1,1,0],
                  [.6 ,1,0,0],
                  [1  ,0,0,0]])
cdict = {'red':   [], 'green': [], 'blue':  []}
for i in range(cdata.shape[0]):
    cdict['red']   += [(cdata[i, 0], cdata[i, 1], cdata[i, 1])]
    cdict['green'] += [(cdata[i, 0], cdata[i, 2], cdata[i, 2])]
    cdict['blue']  += [(cdata[i, 0], cdata[i, 3], cdata[i, 3])]
cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)

pancake_levels = np.arange(0, 961, 20)
pill_levels = np.arange(0, 2001, 50)
movie_levels = np.arange(0, 1801, 10)
