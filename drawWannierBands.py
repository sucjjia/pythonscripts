import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

# cd /Users/cjjia/Documents/Research/CTHFAM/wien2k/PythonScripts

#dir = '/Users/cjjia/Documents/Research/Thermoelectric/CoO_5units/'
#seedname = 'CoO_5units'




def dist(dx, dy, n):
    a = np.array([2.0, 0.0])
    b = np.array([-1.0, 1.732])
    r = dx * a + dy * b
    tol = 0.1
    return np.sqrt(r[0] ** 2 + r[1] ** 2) <= 2 * n + tol
    
def drawWannierBands(ax, dir, seedname, xstart, xend, ystart, yend, Ef, klist):
    pi = 3.1415926
    dim = 5000
    dats = np.zeros((5,5, dim, 3))
    tdim = 0
    f = open(dir + seedname + '_hr.dat', "r")
    f.readline()
    f.readline()
    line = f.readline()
    #print int(line.split()[0])/15
    for i in range(int(line.split()[0])/15 + 1):
        f.readline()
    while True:
        line = f.readline()
        if line:
            nx, ny, nz, oa, ob, val = line.split()[:-1]
            if dist(int(nx),int(ny),3) and abs(float(val)) > 0.01 \
                and int(oa) in [4] and int(ob) in [1]:          
                print line[:-1]
                oa = int(oa)
                ob = int(ob)
                dats[oa-1,ob-1, tdim, 0] = float(nx)
                dats[oa-1,ob-1, tdim, 1] = float(ny)
                dats[oa-1,ob-1, tdim, 2] = float(val)
                tdim = tdim + 1
        else:
            break
    print 'tdim = ,', tdim
    f.close()
    for ai in np.arange(xstart, xend + 1):
        kx, ky, kz = klist[ai]
        #print kx, ky
        tv = np.zeros((5,5))
        tv = np.array(tv, dtype=complex)
        for b1 in range(5):
            for b2 in range(5):
                for i in range(tdim): 
                    ix = dats[b1,b2, i, 0]
                    iy = dats[b1,b2, i, 1]
                    val = dats[b1,b2, i, 2]
                    tv[b1,b2] = tv[b1,b2] + val * np.exp((ix*kx + iy*ky)*2*pi*1j)
        #ax = fig.add_subplot(111)
        eigs = LA.eig(tv)
        kxlist = np.zeros(5)
        kxlist[:] = ai
        ax.scatter(kxlist[:3], eigs[0].real[:3] - Ef, s = 20.0, facecolors='none', 
                   edgecolors='r')
        #ax.plot(kxlist[:3], eigs[0].real[:3] - Ef, 'o', color = 'b', mfc = 'none') 
        #print ai, eigs[0].real[:3] 
    ax.set_ylim((ystart,yend))
    ax.set_xlim((xstart,100))
    #plt.show()                


# Start the main program...
