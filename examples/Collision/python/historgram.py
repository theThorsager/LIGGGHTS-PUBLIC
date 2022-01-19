import csv
import sys
import argparse
import pathlib
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

from os import listdir
from os.path import isfile, join

parser = argparse.ArgumentParser(description='Creates a log plotted histogram for the relative velocities in both the normal and tengent direction.')
parser.add_argument('-f', '--file', dest='file', type=pathlib.Path, required=True,
                    help='file directory path, it will use all .csv files in the directory.')
parser.add_argument('-n', dest='n', type=int, default=40, help='Number of buckets to use for the histogram')

args = parser.parse_args()

mypath = args.file

vel = []
tvel = []

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and ".csv" in f]

for file in onlyfiles:
    with open(str(mypath) +"/"+ file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in spamreader:
            if not "#" in row[0]:
                vel.append(row[0])
                tvel.append(row[1])

vel = np.array(vel, dtype=float)
tvel = np.array(tvel, dtype=float)

nbins = 40
# Histogram
fig = plt.figure()
plt.hist(vel, density=False, bins=nbins)  # density=False would make counts
plt.ylabel('# of collisions')
plt.yscale('log')
plt.xlabel('normal velocity')

plt.show(block=True)
# Histogram
fig = plt.figure()
plt.hist(tvel, density=False, bins=nbins)  # density=False would make counts
plt.ylabel('# of collisions')
plt.yscale('log')
plt.xlabel('tangential velocity')

plt.show(block=True)


