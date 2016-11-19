import numpy as np
import matplotlib.pyplot as plt
    


def plotKpdos(ax, mydir, myfile, orbitals, natoms, Ef, Estart, Eend):
    
    f = open(mydir + myfile + '.scf', 'r')
    line = f.readline()
    ntypes = len(orbitals)
    types = []
    #atoms = ['C'] * natoms
    while(line.split()[0] != 'ATOMIC_SPECIES'):
        line = f.readline()
    for i in range(ntypes):
        line = f.readline()
        types.append(line.split()[0]) 
    print 'ntypes ', ntypes   
    listofAtomseachType = np.zeros((ntypes, natoms)).astype(int)
    lengtheachType = np.zeros(ntypes).astype(int)
    while(len(line.split()) > 0 and line.split()[0] != 'ATOMIC_POSITIONS'):        
        line = f.readline()
    for i in range(natoms):
        line = f.readline()
        print 'line.split()[0] = ', line.split()[0], i
        typei = types.index(line.split()[0])
        #print 'typei = ', typei
        #print 'lengtheachType[typei] = ', lengtheachType[typei]
        listofAtomseachType[typei,lengtheachType[typei]] = i
        lengtheachType[typei] += 1
    
    #print atoms
    data = []
    for i in range(ntypes):
        for j in range(lengtheachType[i]):
            for k in range(len(orbitals[i])):
                filename = mydir + myfile + '.pdos_atm#' + str(listofAtomseachType[i,j]+1) + '(' + types[i] + \
                       ')_wfc#' + str(k+1) + '(' + orbitals[i][k] + ')'
                if j == 0 and k == 0:
                    data = np.loadtxt(filename) 
                else:
                    mydata = np.loadtxt(filename)
                    data[:,1] = data[:,1] + mydata[:,1]
        ax.plot(data[:,1], data[:,0] - Ef, label = types[i])
    plt.legend() 
    plt.ylim(Estart, Eend)   


orbitals = []
orbitals.append(['s', 'p']) #C
orbitals.append(['s', 'p']) #S
orbitals.append(['d', 's']) #Ag
orbitals.append(['s', 'p']) #B
orbitals.append(['s'])      #H

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
mydir = '/Users/cjjia/Documents/Research/Diamondoid/XS_QE/AgS-M9/kpdos/'
myfile = 'AgSM9'
Ef = 3.5141
plotKpdos(ax, mydir, myfile, orbitals, 100, Ef, -4, 4)
plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    