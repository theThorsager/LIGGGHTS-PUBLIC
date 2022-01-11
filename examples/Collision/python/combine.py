import csv
import sys
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline


from os import listdir
from os.path import isfile, join

def get_octant(path):
    with open(path, newline='') as csvfile:
        linereader = csv.reader(csvfile, delimiter=',')
    
        parameters = next(linereader)
        a,b,c,n1,n2 = [float(s) for s in parameters]
    
        parameters = next(linereader)
        xN,yN,zN = [int(s) for s in parameters]

        xSide = []
        ySide = []
        zSide = []
        
        for i in range(yN):
            l = [int(s) for s in next(linereader)]
            xSide.append(l)
        for i in range(xN):
            l = [int(s) for s in next(linereader)]
            ySide.append(l)
        for i in range(xN):
            l = [int(s) for s in next(linereader)]
            zSide.append(l)

        return xSide, ySide, zSide, xN, yN, zN

def add_matrix(A, B):
    for i in range(len(B)):
        for j in range(len(B[0])):
            A[i][j] = A[i][j] + B[i][j]



if len(sys.argv) < 2:
    print("Please provide a path to the .csv files")
    sys.exit()
mypath = sys.argv[1]

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

xSide = np.zeros([100, 100], dtype=int)
ySide = np.zeros([100, 100], dtype=int)
zSide = np.zeros([100, 100], dtype=int)

for file in onlyfiles:
    xSidet, ySidet, zSidet, xN, yN, zN = get_octant(mypath + file)
    add_matrix(xSide, xSidet)
    add_matrix(ySide, ySidet)
    add_matrix(zSide, zSidet)

with open(mypath + onlyfiles[0], newline='') as csvfile:
    linereader = csv.reader(csvfile, delimiter=',')

    parameters = next(linereader)
    a,b,c,n1,n2 = [float(s) for s in parameters]

f = open("combinedfile.csv", "w")

f.write("{}, {}, {}, {}, {}\n".format(a,b,c,n1,n2));
f.write("{}, {}, {}\n".format(xN, yN, zN));

for y in range(0,yN):
    f.write(str(xSide[y][0]))
    for z in range(1,zN):
        f.write(","+str(xSide[y][z]))
    f.write("\n")

for x in range(0,xN):
    f.write(str(ySide[x][0]))
    for z in range(1,zN):
        f.write(","+str(ySide[x][z]))
    f.write("\n")

for x in range(0,xN):
    f.write(str(zSide[x][0]))
    for y in range(1,yN):
        f.write(","+str(zSide[x][y]))
    f.write("\n")

f.close()









