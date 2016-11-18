import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

# cd /Users/cjjia/Documents/Research/CTHFAM/wien2k/PythonScripts

#dir = '/Users/cjjia/Documents/Research/Thermoelectric/CoO_5units/'
#seedname = 'CoO_5units'




def drawqebands(ax, dir, seedname, xstart, xend, ystart, yend, Ef, steps, labels):
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
    
def drawscf(ax, dir, seedname, xstart, xend, ystart, yend, Ef, xs, nlines, color):
    f = open(dir + seedname + '.scf.out.energy', "r")
    for i in range(len(xs)):
        line = f.readline()
        line = f.readline()
        Elist = []
        for k in range(nlines):
            line = f.readline()
            Elist = Elist + (map(float, line.split()[:]))
        a = xs[i]
        print a
        x = [a for j in range(len(Elist))]
        line = f.readline()
        Ep = map(lambda x: x - Ef, Elist) 
        #ax.scatter(x, Ep, edgecolors = 'k', facecolors='r', s = 20.0)   
        ax.plot(x, Ep,'o', markersize = 4.5, color = color)
    ax.set_ylim((ystart,yend))
    ax.set_xlim((xstart,xend))

def dist(dx, dy, n):
    a = np.array([2.0, 0.0])
    b = np.array([-1.0, 1.732])
    r = dx * a + dy * b
    tol = 0.1
    return np.sqrt(r[0] ** 2 + r[1] ** 2) <= 2 * n + tol
    
def drawwannierbands(ax, dir, seedname, xstart, xend, ystart, yend, Ef, klist):
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

# $$$ $$$$$

# Start the main program...

def effectiveMass(ax, bands, n, step, index, direction):
    nbnd = len(bands[0,:])
    nks = int(bands[:,0])
    if direction == 'r' or 'right':
        x = np.arange(n, n+step+1)
    elif direction == 'l' or 'left':
        x = np.arange(n, n-step-1, -1)
    for i in range(nbnd):
        subbands = [bands[j,i] for j in x]
        ax.plot(x, subbands[:], 'k')
    y = [bands[j,n] for j in x]
    
    w1 = np.exp(-(x - x[0])**2/(sigma**2))
    p = np.polyfit(x, y, 2, w=w1)
    y1p = p[2] + p[1] * x + p[0] * x**2
    ax.plot(x, y1p, '--')
 




fig = plt.figure(figsize=(16,8))

ax = fig.add_subplot(171)

dir = '/Users/cjjia/Documents/Research/Diamondoid/QE/CuSada_nc_oncv_Nov2015/'
Ef = 3.575
seedname = 'CuSada_cori_484'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926*0.8,3.1415926*0.8+0.01,3.1415926*0.8/60)
y1 = bands[:,479]
y2 = bands[:,480]
yp1 = np.zeros(121)
yp2 = np.zeros(121)
yp1[:61] = y1[61::-1]
yp2[:61] = y2[61::-1]
yp1[60:] = y1[:]
yp2[60:] = y2[:]
#yp1[:17] = y1[::-1]
#yp1[17:33] = y1[1:]
#yp2[:17] = y2[::-1]
#yp2[17:33] = y2[1:]
ax.plot(x,yp1)
ax.plot(x,yp2)
ax.set_ylim((-1.0, 2.0))
ax.set_xlim((0, 3.1416))
ticks = [0, 3.1416]
labels = ['G', 'Y']
ax.set_xticks(ticks)
ax.set_xticklabels(labels)


ax.set_ylabel('Energy (eV)')


sigma = 0.4

fac = 1.05459**2 / 1.602 / (6.8)**2 / 0.091095 / 2

w1 = np.exp(-(x - x[60])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
ax.set_ylim((-1.0, 2.0))
#ax.set_xlim((-3.1416/2, 3.1416/2))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 0.2, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")


w2 = np.exp(-(x - 3.1416/2)**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_ylim((-1.0, 2.0))
#ax.set_xlim((-3.1416/2, 3.1416/2))
print p[0],p[1]/(2*p[0])
ax.text(1.57, 1.3, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('CuSada: G-Y')


ax = fig.add_subplot(172)

dir = '/Users/cjjia/Documents/Research/Diamondoid/QE/4DiCu_nc_oncv_Nov2015/'
Ef = 2.160
seedname = '4DiCu_484_2'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926,2*3.1415926+0.01,3.1415926/40)
y1 = bands[:,399]
y2 = bands[:,400]
yp1 = np.zeros(121)
yp2 = np.zeros(121)
yp1[:41] = y1[40::-1]
yp2[:41] = y2[40::-1]
yp1[40:81] = y1[:41]
yp2[40:81] = y2[:41]
yp1[80:] = y1[40::-1]
yp2[80:] = y2[40::-1]
ax.plot(x,yp1)
ax.plot(x,yp2)
ax.set_xlim((0, 3.1416))


w1 = np.exp(-(x - x[80])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
print p[0], p[1]/(2*p[0])
ax.text(1.57, 0.5, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")


w2 = np.exp(-(x - x[40])**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_xlim((0, 3.1416))
ax.set_ylim((-1,4))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 2.0, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('4DiCu: G-Y')
ax.set_xticks(ticks)
ax.set_xticklabels(labels)

#ax = fig.add_subplot(122)


sigma = 1

ax = fig.add_subplot(173)

dir = '/Users/cjjia/Documents/Research/Diamondoid/XS_QE/AgSAda/HaoYan/'
Ef = 4.0794
seedname = 'AgSAda'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926,1*3.1415926+0.01,3.1415926/16)
y1 = bands[:,639]
y2 = bands[:,640]
yp1 = np.zeros(33)
yp2 = np.zeros(33)
yp1[16:33] = y1[9:26]
yp2[16:33] = y2[9:26]
yp1[:17] = y1[16::-1]
yp2[:17] = y2[16::-1]
ax.plot(x,yp1)
ax.plot(x,yp2)
ax.set_xlim((0, 3.1416))


w1 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
p[0] = p[0] * 1.1
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
print p[0], p[1]/(2*p[0])
ax.text(1.57, -1.0, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")


w2 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
p[0] = p[0] * 1.2
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_xlim((0, 3.1416))
ax.set_ylim((-2.5,2.5))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 1.0, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('AgSAda: G-Y')
ax.set_xticks(ticks)
ax.set_xticklabels(labels)



sigma = 1

ax = fig.add_subplot(174)

dir = '/Users/cjjia/Documents/Research/Diamondoid/XS_QE/Ag-S-mDia/'
Ef = 1.6691
seedname = 'AgSmDia.0'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926,1*3.1415926+0.01,3.1415926/16)
y1 = bands[:,199]
y2 = bands[:,200]
yp1 = np.zeros(33)
yp2 = np.zeros(33)
yp1[16:33] = y1[20:37]
yp2[16:33] = y2[20:37]
yp1[:17] = y1[16::-1]
yp2[:17] = y2[16::-1]
ax.plot(x,yp1)
ax.plot(x,yp2)
ax.set_xlim((0, 3.1416))


w1 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
p[0] = p[0] * 1.3
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
print p[0], p[1]/(2*p[0])
ax.text(1.57, -1.0, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")


w2 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
p[0] = p[0] * 1.4
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_xlim((0, 3.1416))
ax.set_ylim((-2.5,2.5))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 1.0, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('AgSmDia: G-X')
labels = ['G', 'X']
ax.set_xticks(ticks)
ax.set_xticklabels(labels)


ax = fig.add_subplot(175)

dir = '/Users/cjjia/Documents/Research/Diamondoid/XS_QE/AgS-M1/HaoYan/'
Ef = 2.1277
seedname = 'AgSM1'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926,1*3.1415926+0.01,3.1415926/12)
y1 = bands[:,591]
y2 = bands[:,592]
yp1 = np.zeros(25)
yp2 = np.zeros(25)
yp1[12:25] = y1[24:37]
yp2[12:25] = y2[24:37]
yp1[:13] = y1[12::-1]
yp2[:13] = y2[12::-1]
ax.plot(x,yp1)
ax.plot(x,yp2)
ax.set_xlim((0, 3.1416))


w1 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
p[0] = p[0] * 1.0
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
print p[0], p[1]/(2*p[0])
ax.text(1.57, 0.2, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")


w2 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
p[0] = p[0] * 1.0
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_xlim((0, 3.1416))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 2.5, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('AgSM1: Y-A')
labels = ['Y', 'A']
ax.set_xticks(ticks)
ax.set_xticklabels(labels)
ax.set_ylim((-1,4))


ax = fig.add_subplot(176)

dir = '/Users/cjjia/Documents/Research/Diamondoid/XS_QE/AgS-M9/HaoYan/'
Ef = 4.2871
seedname = 'AgSM9'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926,3.1415926+0.01,3.1415926/16)
y1 = bands[:,147]
y2 = bands[:,148]
yp1 = np.zeros(33)
yp2 = np.zeros(33)
yp1[16:33] = y1[9:26]
yp2[16:33] = y2[9:26]
yp1[:17] = y1[16::-1]
yp2[:17] = y2[16::-1]
ax.plot(x,yp1)
ax.plot(x,yp2)
ax.set_ylim((-2.0, 2.0))
ax.set_xlim((0, 3.1416))
ticks = [0, 3.1416]
labels = ['G', 'Y']
ax.set_xticks(ticks)
ax.set_xticklabels(labels)


sigma = 0.4

w1 = np.exp(-(x - x[28])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
#ax.set_xlim((-3.1416/2, 3.1416/2))
print p[0], p[1]/(2*p[0])
ax.text(1.57, -1.2, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

sigma = 1

w2 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
p[0] = p[0] * 1.3
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_xlim((0, 3.1416))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 1.2, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('AgSM9: G-Y')
ax.set_ylim((-2.0, 2.5))



ax = fig.add_subplot(177)

dir = '/Users/cjjia/Documents/Research/Diamondoid/XS_QE/Hg-S-M1/Hg-S-M1/'
Ef = 2.1332
seedname = 'HgSM1'
bands = effectiveMass(dir, seedname, Ef)
x = np.arange(-3.1415926,3.1415926+0.01,3.1415926/16)
y1 = bands[:,119]
y2 = bands[:,120]
yp1 = np.zeros(33)
yp2 = np.zeros(33)
yp1[16:33] = y1[0:17]
yp2[16:33] = y2[0:17]
yp1[:17] = y1[16::-1]
yp2[:17] = y2[16::-1]
ax.plot(x,yp1)
ax.plot(x,yp2)
#ax.set_ylim((-2.0, 2.0))
ax.set_xlim((0, 3.1416))
ticks = [0, 3.1416]
labels = ['X', 'G']
ax.set_xticks(ticks)
ax.set_xticklabels(labels)


sigma = 1

w1 = np.exp(-(x - x[16])**2/(sigma**2))
p = np.polyfit(x, yp1, 2, w=w1)
y1p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y1p, '--')
#ax.set_xlim((-3.1416/2, 3.1416/2))
print p[0], p[1]/(2*p[0])
ax.text(1.57, -0.5, 'm_h = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

sigma = 0.6

w2 = np.exp(-(x - x[23])**2/(sigma**2))
p = np.polyfit(x, yp2, 2, w=w2)
p[0] = p[0] * 1.05
y2p = p[2] + p[1] * x + p[0] * x**2
ax.plot(x, y2p, '--')
ax.set_xlim((0, 3.1416))
print p[0], p[1]/(2*p[0])
ax.text(1.57, 2, 'm_e = ' + str(round(fac/p[0],2)), ha="center", va="bottom", size="large")

ax.set_title('HgSM1: X-G')
ax.set_ylim((-2.0, 3))





plt.show()