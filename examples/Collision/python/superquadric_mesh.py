import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math as m
import argparse
import pathlib
import sys
import os
import csv


def get_args():
    parser = argparse.ArgumentParser(description='Create cube-spherical mesh')
    parser.add_argument('-f', '--file', dest='file', type=pathlib.Path,
                        help='file path')

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return args


def fexp(x, p):
    return np.sign(x)*(abs(x)**p)


def plot_superquadric(fig=None,ax=None,show=True,n=16):

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(projection='3d')

    # Account for start and end being the same in a circle
    n_samples = n + 1

    eta = np.linspace(-np.pi/2, np.pi/2, n_samples, endpoint=True)
    omega = np.linspace(-np.pi, np.pi, n_samples, endpoint=True)
    eta, omega = np.meshgrid(eta, omega)

    a = 1#0.005
    b = 1#0.005
    c = 1#0.00225

    #a = a/2
    #b = b/2
    #c = c/2

    e1 = 0.5 # theta roundness
    e2 = 0.16667 # phi roundness
    n1 = 2/e1 # Blockiness parameter
    n2 = 2/e2 # Blockiness parameter

    x = a * fexp(np.cos(eta), e1) * fexp(np.cos(omega), e2)
    y = b * fexp(np.cos(eta), e1) * fexp(np.sin(omega), e2)
    z = c * fexp(np.sin(eta), e1)

    ax.plot_surface(x, y, z, color='b',zorder=1)
    #ax.plot_wireframe(x, y, z, color='b')
    if show:
        plt.show()

    print(f(x,y,z,a,b,c,n1,n2))

    return fig,ax


def g(x):
    return np.tan(np.pi/4*x)


def plot_grid(fig=None,ax=None,show=True):

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(projection='3d')

    a = 1#0.005
    b = 1#0.005
    c = 1#0.00225

    #a = a/2
    #b = b/2
    #c = c/2

    e1 = 0.5 # theta roundness
    e2 = 0.16667 # phi roundness
    n1 = 2/e1 # Blockiness parameter
    n2 = 2/e2 # Blockiness 
    #e1 = 2/n1
    #e2 = 2/n2

    unit = 1
    gridsize = 10

    u = np.linspace(0, unit, gridsize)
    v = np.linspace(0, unit, gridsize)
    u,v = np.meshgrid(u,v)

    x = u
    y = v
    z = np.ones(np.shape(u))
    #ax.plot_wireframe(x, y, z, color='r')
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    ax.plot_wireframe(x, y, z, color='m')

    x = u
    y = np.ones(np.shape(u))
    z = v
    #ax.plot_wireframe(x, y, z, color='r')
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    ax.plot_wireframe(x, y, z, color='m')

    x = np.ones(np.shape(u))
    y = u
    z = v
    #ax.plot_wireframe(x, y, z, color='r')
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    ax.plot_wireframe(x, y, z, color='m')


    # Shift grid out a bit so it is not covered by superquadric
    x = 1.05 * x
    y = 1.05 * y
    z = 1.05 * z

    #ax.plot_wireframe(x, y, z, color='r')

    if show:
        plt.show()

    return fig,ax


# The superquadric equation
def f(x,y,z,a,b,c,n1,n2):
    f = np.power(np.power(np.abs(x/a), n2) + np.power(np.abs(y/b), n2),n1/n2) + np.power(np.abs(z/c), n1) - 1
    return f

def get_beta(x,y,z,a,b,c,n1,n2):
    r = f(x,y,z,a,b,c,n1,n2)
    beta = np.power(r+1,-1/n1)
    return beta


def get_projected_coord(x,y,z,a,b,c,n1,n2):
    beta = get_beta(x,y,z,a,b,c,n1,n2)
    x = beta * x
    y = beta * y
    z = beta * z
    return x,y,z


def get_octant(path):
    with open(path, newline='') as csvfile:
        linereader = csv.reader(csvfile, delimiter=',')
        parmeters = next(linereader)
        xN,yN,zN = [int(s) for s in parmeters]
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


if __name__ == "__main__":

    args = get_args()
    get_octant(args.file)

    fig,ax = plot_superquadric(show=False,n=8)
    plot_grid(fig,ax)

    a = 1#0.005
    b = 1#0.005
    c = 1#0.00225

    a = a/2
    b = b/2
    c = c/2

    thetaroundness = 2
    phiroundness = 6
    n1 = 2/thetaroundness
    n2 = 2/phiroundness

    x,y,z = get_projected_coord(1,1,1,a,b,c,n1,n2)
    val = f(x,y,z,a,b,c,n1,n2)

    print(x,y,z)
    print(val)

    x,y,z = get_projected_coord(0,0,0.5,a,b,c,n1,n2)
    val = f(x,y,z,a,b,c,n1,n2)

    print(x,y,z)
    print(val)

    x,y,z = get_projected_coord(0.01,0.01,0.01,a,b,c,n1,n2)
    val = f(x,y,z,a,b,c,n1,n2)

    print(x,y,z)
    print(val)



'''
  different approach but just as valid
  n = 50;
  etamax = np.pi/2;
  etamin = -np.pi/2;
  wmax = np.pi;
  wmin = -np.pi;
  deta = (etamax-etamin)/n;
  dw = (wmax-wmin)/n;
  [i,j] = meshgrid(1:n+1,1:n+1)
  eta = etamin + (i-1) * deta;
  w   = wmin + (j-1) * dw;
  x = a(1) .* sign(cos(eta)) .* abs(cos(eta)).^epsilon(1) .* sign(cos(w)) .* abs(cos(w)).^epsilon(1);
  y = a(2) .* sign(cos(eta)) .* abs(cos(eta)).^epsilon(2) .* sign(sin(w)) .* abs(sin(w)).^epsilon(2);
  z = a(3) .* sign(sin(eta)) .* abs(sin(eta)).^epsilon(3);
'''
