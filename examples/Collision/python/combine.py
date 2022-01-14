import csv
import sys
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pathlib
import glob

from os import listdir
from os.path import isfile, join


def get_args():
    parser = argparse.ArgumentParser(description='Combine multiple octant maps into one map')
    parser.add_argument('-f', '--files', dest='files', nargs='+', metavar="path", type=str,
                        help='Path of files to combine. Can use * as wildcard for folder/file names. List of paths is also accepted')
    parser.add_argument('-d', '--dir', dest='dir',
                        help='Path to directory of files to combine')
    parser.add_argument('-o', '--out', dest='out', default='combinedfile.csv',
                        help='Path to outputfile')

    args = parser.parse_args()

    if args.files and args.dir:
        print('--file and --dir are mutually exclusive')
        sys.exit(1)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not args.files and not args.dir:
        print('Input is needed')
        sys.exit(1)

    return args

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


args = get_args()

if args.dir:
    onlyfiles = [args.dir + f for f in listdir(args.dir) if isfile(join(args.dir, f))]
elif args.files:
    # Needs specials handling as user can pass wildcards in two ways: -f example_dir/*.csv or -f 'example_dir/*.csv'
    # Linux automatically expands unquoted wildcards as example_dir/*.csv while windows does not
    # Any unexpanded wildcards need to be expanded
    expanded = [glob.glob(f) for f in args.files]
    onlyfiles = [f for f_list in expanded for f in f_list if isfile(f)]

xSide = np.zeros([100, 100], dtype=int)
ySide = np.zeros([100, 100], dtype=int)
zSide = np.zeros([100, 100], dtype=int)

for file in onlyfiles:
    xSidet, ySidet, zSidet, xN, yN, zN = get_octant(file)
    add_matrix(xSide, xSidet)
    add_matrix(ySide, ySidet)
    add_matrix(zSide, zSidet)

with open(onlyfiles[0], newline='') as csvfile:
    linereader = csv.reader(csvfile, delimiter=',')

    parameters = next(linereader)
    a,b,c,n1,n2 = [float(s) for s in parameters]

with open(args.out,'w', newline='') as f:
    f.write("{}, {}, {}, {}, {}\n".format(a,b,c,n1,n2))
    f.write("{}, {}, {}\n".format(xN, yN, zN))

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










