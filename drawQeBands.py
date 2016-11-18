import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

# cd /Users/cjjia/Documents/Research/CTHFAM/wien2k/PythonScripts
#dir = '/Users/cjjia/Documents/Research/Thermoelectric/CoO_5units/'
#seedname = 'CoO_5units'


def drawQeBands(ax, dir, seedname, xstart, xend, ystart, yend, Ef, steps, labels):
    f = open(dir + seedname + '.qebands.dat', "r")
    line = f.readline()
    nbnd = int(line.split()[2][:-1])
    nks = int(line.split()[4])
    print 'nbnd = ', nbnd
    print 'nks = ', nks
    bands = np.zeros((nks, nbnd))
    klist = []
    for i in range(nks):
        line = f.readline()
        #print i, line
        klist.append(map(float, line.split()[:]))
        x = [i for j in range(nbnd)]
        E = []
        jend = nbnd / 10 + 1
        if nbnd % 10 == 0:
            jend = nbnd / 10
        for j in range(jend):
            line = f.readline()
            Es = line.split()[:]
            if(len(Es[0]) == 8):
                E = E + map(float, line.split()[:])
            else:
                for k in range(len(line)/8):
                    #print 'print here ', line[8*k: 8*(k+1)]
                    E.append(float(line[8*k: 8*(k+1)]))
        Ep = map(lambda x: x - Ef, E)
        #ax.plot(x,Ep,'ob', markersize = 3.0)
        bands[i] = np.array(Ep)
    x = np.arange(nks)
    for i in range(12):
        x[nks-12+i] = x[nks-12+i] + 0.5 * i
        
    for i in range(nbnd):
        if i == 480:
            print bands[:-3,i]
            print min(bands[:-3,i])
        ax.plot(x, bands[:,i], 'k')
    ax.set_ylim((ystart,yend))
    ax.set_xlim((xstart,xend))
    for i in range(len(steps)):
        xstart = xstart + steps[i]
        ax.plot([xstart, xstart], [ystart, yend], 'k-')
    ticks = []
    t = 0
    ticks.append(t)
    for i in range(len(steps)):
        t = t + steps[i]
        ticks.append(t)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    #plt.show()
    #return klist
    x = np.arange(nks)
    y = np.arange(nks)
    y[:] = 0
    ax.plot(x,y,'--k')
    return klist, bands
    

# *****************  Example  ****************************
'''
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

mydir = '/Users/cjjia/Documents/Research/ExcitonicInsulator/QE/TiS2/'
Ef = 7.7659
seedname = 'TiS2'
steps = [10, 10, 10, 10, 40, 10, 10]
labels = ['G', 'M', 'K', 'G', 'A', 'L', 'H', 'A']
klist, bands = drawqebands(ax, mydir, seedname, 0, 100, -6, 2, Ef, steps, labels)
'''