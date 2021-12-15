import csv
import sys
import argparse
import pathlib
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline


from os import listdir
from os.path import isfile, join


parser = argparse.ArgumentParser(description='Convert a directory of exctencive collision data files into a file with the number of collision per face.')
parser.add_argument('-f', '--file', dest='file', type=pathlib.Path, required=True,
                    help='file directory path, it will use all .csv files in the directory.')
parser.add_argument('-o', '--output', dest='output', type=pathlib.Path,
                    default="octo_counts.csv",
                    help='file output')
parser.add_argument('-ns', dest='ns', type=int, nargs=3, metavar='n',  help='Optional number of buckets in each direction')
parser.add_argument('-N',  dest='N', default=300,  type=int, help='Approximate total number of buckets if no exact ones are provided.')
parser.add_argument('-s', '--shape', dest='shape', type=float, nargs=5, metavar='s', help='Optional overwrite for shape and blockiness values of the shape. Will be read from the files elsewise.') 

args = parser.parse_args()

mypath = args.file

pos = []

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and ".csv" in f]

for file in onlyfiles:
    with open(str(mypath) +"/"+ file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in spamreader:
            if not "#" in row[0]:
                pos.append([row[2], row[3], row[4]])
            else:
                shapeX = float(row[1])
                shapeY = float(row[2])
                shapeZ = float(row[3])
                block1 = float(row[4])
                block2 = float(row[5])

if args.shape:
    shapeX = args.shape[0]
    shapeY = args.shape[1]
    shapeZ = args.shape[2]
    block1 = args.shape[3]
    block2 = args.shape[4]

if args.ns:
    nx = args.ns[0]
    ny = args.ns[1]
    nz = args.ns[2]
else:
    temp = shapeY + shapeZ + shapeY*shapeZ/shapeX
    nx = int(np.sqrt(args.N * shapeX / temp))
    ny = int(nx * shapeY / shapeX)
    nz = int(nx * shapeZ / shapeX)

x_oct = np.zeros((ny,nz), dtype=int)
y_oct = np.zeros((nx,nz), dtype=int)
z_oct = np.zeros((nx,ny), dtype=int)

pos = np.array(pos, dtype=float)

# unit cube projection
pos = np.abs(pos)
shape = np.array([1/shapeX, 1/shapeY, 1/shapeZ])
pos = pos * shape

pmax = 1 / np.amax(pos, 1)
pos = (pos.T * pmax).T

# unit cube indexing 
for vec in pos:
    if vec[0] >= vec[1] and vec[0] >= vec[2]:
        y = int(min(ny-1, np.floor(vec[1]*ny)));
        z = int(min(nz-1, np.floor(vec[2]*nz)));
        x_oct[y][z]+=1
    elif vec[1] >= vec[2]:
        x = int(min(nx-1, np.floor(vec[0]*nx)));
        z = int(min(nz-1, np.floor(vec[2]*nz)));
        y_oct[x][z]+=1
    else: 
        x = int(min(nx-1, np.floor(vec[0]*nx)));
        y = int(min(ny-1, np.floor(vec[1]*ny)));
        z_oct[x][y]+=1

f = open(args.output, "w")

f.write("{}, {}, {}, {}, {}\n".format(shapeX, shapeY, shapeZ, block1, block2));
f.write("{}, {}, {}\n".format(nx, ny, nz));

for y in range(0,ny):
    f.write(str(x_oct[y][0]))
    for z in range(1,nz):
        f.write(","+str(x_oct[y][z]))
    f.write("\n")

for x in range(0,nx):
    f.write(str(y_oct[x][0]))
    for z in range(1,nz):
        f.write(","+str(y_oct[x][z]))
    f.write("\n")

for x in range(0,nx):
    f.write(str(z_oct[x][0]))
    for y in range(1,ny):
        f.write(","+str(z_oct[x][y]))
    f.write("\n")

f.close()


















