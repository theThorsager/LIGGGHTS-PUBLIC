import csv
import sys
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline


from os import listdir
from os.path import isfile, join

if len(sys.argv) < 2:
    print("Please provide a path to the .csv files")
    sys.exit()
mypath = sys.argv[1]

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

vel = []
tvel = []
posX = []
posY = []
posZ = []

for file in onlyfiles:
    with open(mypath + file, newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in spamreader:
            #print(', '.join(row))
            vel.append(row[0])
            tvel.append(row[1])
            posX.append(row[2])
            posY.append(row[3])
            posZ.append(row[4])

vel = np.array(vel, dtype=float)
tvel = np.array(tvel, dtype=float)
posX = np.array(posX, dtype=float)
posY = np.array(posY, dtype=float)
posZ = np.array(posZ, dtype=float)

# Histogram
fig = plt.figure()
plt.hist(vel, density=False, bins=30)  # density=False would make counts
plt.ylabel('# of collisions')
plt.yscale('log')
plt.xlabel('normal velocity')

plt.show(block=True)
# Histogram
fig = plt.figure()
plt.hist(tvel, density=False, bins=30)  # density=False would make counts
plt.ylabel('# of collisions')
plt.yscale('log')
plt.xlabel('tangential velocity')

plt.show(block=True)

print("Number of Collisions: " + str(len(vel)))
# Copy to Symmetries
# Mirror X
posX = np.concatenate((posX, [-i for i in posX]));
posY = np.concatenate((posY, posY))
posZ = np.concatenate((posZ, posZ))
vel =  np.concatenate((vel , vel))
# Mirror Y
posX = np.concatenate((posX, posX))
posY = np.concatenate((posY, [-i for i in posY]))
posZ = np.concatenate((posZ, posZ))
vel =  np.concatenate((vel , vel))
# Mirror Z
posX = np.concatenate((posX, posX))
posY = np.concatenate((posY, posY))
posZ = np.concatenate((posZ, [-i for i in posZ]))
vel =  np.concatenate((vel , vel))

# Point Cloud
fig = plt.figure()
ax = plt.axes(projection='3d')
col = vel / vel.max() * 255
every = 1000
ax.scatter3D(posX[::every], posY[::every], posZ[::every], c=vel[::every], cmap='autumn', depthshade=0);
size = 0.01
dumbx = [size, size, size, size, -size, -size, -size, -size]
dumby = [size, size, -size, -size, -size, -size, size, size]
dumbz = [size, -size, size, -size, size, -size, size, -size]
dcol =  [1, 1, 1, 1, 1, 1, 1, 1]
ax.scatter3D(dumbx, dumby, dumbz, c=dcol)
plt.xlim(-size, size)
plt.ylim(-size, size)
#plt.zlim(-size, size)

plt.show(block=True)


