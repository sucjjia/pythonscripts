import numpy as np
import matplotlib.pyplot as plt



mydir = '/Users/cjjia/Documents/Research/Glassy_Carbon/Progress_MD_sructure/40GPAvcrelax/'
myfile = 'gc'
natoms = 250
        
def RelaxToXtl(mydir, myfile, natoms):
    v = np.zeros((3,3))
    a0 = 0.0
    pos = []
    
    f = open(mydir + myfile + '.rx.out', 'r')
    line = f.readline()
    while(line):
        line = f.readline()
        if line[:15] == 'CELL_PARAMETERS':
            a0 = float(line.split()[2][:8])
            for i in range(3):
                line = f.readline()
                v[i][0:3] = np.array(map(float, line.split()[0:3])) * a0
        if line[:16] == 'ATOMIC_POSITIONS':
            for i in range(natoms):
                line = f.readline()
                pos.append(line)
    f.close()
    
    f = open(mydir + myfile + '.xtl', 'w')
    f.write('TITLE \n')
    f.write('CELL\n')
    a1 = np.sqrt(v[0][0] **2 + v[0][1] **2 + v[0][2] **2)  
    a2 = np.sqrt(v[1][0] **2 + v[1][1] **2 + v[1][2] **2)  
    a3 = np.sqrt(v[2][0] **2 + v[2][1] **2 + v[2][2] **2)  
    ang1 = np.arccos(np.dot(v[1][:], v[2][:]) / (a2*a3)) / np.pi * 180
    ang2 = np.arccos(np.dot(v[0][:], v[2][:]) / (a1*a3)) / np.pi * 180
    ang3 = np.arccos(np.dot(v[0][:], v[1][:]) / (a1*a2)) / np.pi * 180
    f.write(str(a1) + ' ')
    f.write(str(a2) + ' ')
    f.write(str(a3) + ' ')
    f.write(str(ang1) + ' ')
    f.write(str(ang2) + ' ')
    f.write(str(ang3) + '\n')
    f.write('SYMMETRY NUMBER 1\n')
    f.write('SYMMETRY LABEL  P1\n')
    f.write('ATOMS\n')
    f.write('NAME         X           Y           Z\n')
    for i in range(natoms):
        f.write(pos[i])
    f.write('EOF\n')
    f.close()