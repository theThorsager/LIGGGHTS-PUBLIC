import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math as m
import argparse
import pathlib
from mayavi import mlab
import sys
import os
import csv
from typing import NamedTuple


class Superquadric(NamedTuple):
    a: float
    b: float
    c: float
    e1: float
    e2: float
    n1: float
    n2: float

class Spherecube(NamedTuple):
    xN: float
    yN: float
    zN: float


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


def plot_superquadric_mayavi(superquadric,n=16):

    # Account for start and end being the same in a circle
    n_samples = n + 1

    eta = np.linspace(-np.pi/2, np.pi/2, n_samples, endpoint=True)
    omega = np.linspace(-np.pi, np.pi, n_samples, endpoint=True)
    eta, omega = np.meshgrid(eta, omega)

    a = superquadric.a
    b = superquadric.b
    c = superquadric.c

    e1 = superquadric.e1
    e2 = superquadric.e2
    n1 = superquadric.n1
    n2 = superquadric.n2

    x = a * fexp(np.cos(eta), e1) * fexp(np.cos(omega), e2)
    y = b * fexp(np.cos(eta), e1) * fexp(np.sin(omega), e2)
    z = c * fexp(np.sin(eta), e1)


    fig = mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

    ax_scale = [1.0, 1.0, 1.0] # Need to change later
    ax_ranges = [-2, 2, -2, 2, -2, 2]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)

    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    #superquad_surf.actor.actor.scale = ax_scale
    #mlab.outline(superquad_surf, color=(.7, .7, .7), extent=ax_extent)
    #mlab.axes(superquad_surf, color=(.7, .7, .7), extent=ax_extent, ranges=ax_ranges, xlabel='x', ylabel='y', zlabel='z')

    mlab.view(.0, -5.0, 4)
    mlab.show()

    #return fig


def prep_scalar_value(side):
    # Each square consists of 2 triangles
    scalar = np.transpose(side)
    scalar = np.array(scalar).flatten()
    scalar = np.append(scalar,scalar)

    return scalar

def plot_grid_mayavi(superquadric, spherecube, xSide, ySide, zSide):

    a = superquadric.a
    b = superquadric.b
    c = superquadric.c

    e1 = superquadric.e1
    e2 = superquadric.e2
    n1 = superquadric.n1
    n2 = superquadric.n2

    unit = 1
    xN = spherecube.xN + 1 # For fence post problem
    yN = spherecube.yN + 1 # For fence post problem
    zN = spherecube.zN + 1 # For fence post problem

    fig = mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

    ax_scale = [1.0, 1.0, 1.0] # Need to change later
    ax_ranges = [-2, 2, -2, 2, -2, 2]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)

    # Create meshgrid for side X,-X
    u = np.linspace(-unit, unit, yN)
    v = np.linspace(-unit, unit, zN)
    u,v = np.meshgrid(u,v)

    # Positive x side
    x = np.ones(np.shape(u))
    y = u
    z = v
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    w = prep_scalar_value(xSide)
    #superquad_surf = mlab.mesh(x, y, z, representation='wireframe', colormap='Oranges')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf.mlab_source.dataset.cell_data.scalars = w
    superquad_surf.mlab_source.dataset.cell_data.scalars.name = 'X Cell data'
    superquad_surf.mlab_source.dataset.cell_data.update()
    superquad_surf.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(superquad_surf,cell_scalars='X Cell data') 
    surf = mlab.pipeline.surface(mesh2,colormap='Blues')

    # Negative x side
    x = -x
    #superquad_surf = mlab.mesh(x, y, z, representation='wireframe', colormap='Blues')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf.mlab_source.dataset.cell_data.scalars = w
    superquad_surf.mlab_source.dataset.cell_data.scalars.name = '-X Cell data'
    superquad_surf.mlab_source.dataset.cell_data.update()
    superquad_surf.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(superquad_surf,cell_scalars='-X Cell data') 
    surf = mlab.pipeline.surface(mesh2,colormap='Blues')

    # Create meshgrid for side X,-X
    u = np.linspace(-unit, unit, xN)
    v = np.linspace(-unit, unit, zN)
    u,v = np.meshgrid(u,v)

    # Positive y side
    x = u
    y = np.ones(np.shape(u))
    z = v
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    w = prep_scalar_value(ySide)
    #superquad_surf = mlab.mesh(x, y, z, representation='wireframe', colormap='Oranges')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf.mlab_source.dataset.cell_data.scalars = w
    superquad_surf.mlab_source.dataset.cell_data.scalars.name = 'Y Cell data'
    superquad_surf.mlab_source.dataset.cell_data.update()
    superquad_surf.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(superquad_surf,cell_scalars='Y Cell data') 
    surf = mlab.pipeline.surface(mesh2,colormap='Blues')

    # Negative y side
    y = -y
    #superquad_surf = mlab.mesh(x, y, z, representation='wireframe', colormap='Blues')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf.mlab_source.dataset.cell_data.scalars = w
    superquad_surf.mlab_source.dataset.cell_data.scalars.name = '-Y Cell data'
    superquad_surf.mlab_source.dataset.cell_data.update()
    superquad_surf.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(superquad_surf,cell_scalars='-Y Cell data') 
    surf = mlab.pipeline.surface(mesh2,colormap='Blues')


    # Create meshgrid for side X,-X
    u = np.linspace(-unit, unit, xN)
    v = np.linspace(-unit, unit, yN)
    u,v = np.meshgrid(u,v)

    # Positive z side
    x = u
    y = v
    z = np.ones(np.shape(u))
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    w = prep_scalar_value(zSide)
    #superquad_surf = mlab.mesh(x, y, z, representation='wireframe', colormap='Oranges')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf.mlab_source.dataset.cell_data.scalars = w
    superquad_surf.mlab_source.dataset.cell_data.scalars.name = 'Z Cell data'
    superquad_surf.mlab_source.dataset.cell_data.update()
    superquad_surf.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(superquad_surf,cell_scalars='Z Cell data') 
    surf = mlab.pipeline.surface(mesh2,colormap='Blues')

    # Negative z side
    z = -z
    #superquad_surf = mlab.mesh(x, y, z, representation='wireframe', colormap='Blues')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf.mlab_source.dataset.cell_data.scalars = w
    superquad_surf.mlab_source.dataset.cell_data.scalars.name = '-Z Cell data'
    superquad_surf.mlab_source.dataset.cell_data.update()
    superquad_surf.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(superquad_surf,cell_scalars='-Z Cell data') 
    surf = mlab.pipeline.surface(mesh2,colormap='Blues')

    #mlab.colorbar(title='Phase', orientation='vertical', nb_labels=3)
    #mlab.colorbar()

    #mlab.view(.0, -5.0, 4)
    mlab.show()


def plot_superquadric_matplotlib(fig=None,ax=None,show=True,n=16):

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

    return fig,ax


def g(x):
    return np.tan(np.pi/4*x)


def plot_grid_matplotlib(fig=None,ax=None,show=True):

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
        
        # Since values are octant symmetric, divide them by 8
        for i in range(yN):
            l = [float(s)/8 for s in next(linereader)]
            xSide.append(l)
        for i in range(xN):
            l = [float(s)/8 for s in next(linereader)]
            ySide.append(l)
        for i in range(xN):
            l = [float(s)/8 for s in next(linereader)]
            zSide.append(l)

        # Expand dimension from octant to full cube
        sc = Spherecube(2*xN,2*yN,2*zN)

        # Expand sides from octant to full  
        xSide = [*xSide[::-1], *xSide]
        xSide = [[*i[::-1],*i] for i in xSide]

        ySide = [*ySide[::-1], *ySide]
        ySide = [[*i[::-1],*i] for i in ySide]

        zSide = [*zSide[::-1], *zSide]
        zSide = [[*i[::-1],*i] for i in zSide]

        return sc, xSide, ySide, zSide


if __name__ == "__main__":

    args = get_args()
    sc, xSide, ySide, zSide = get_octant(args.file)

    a = 1#0.005
    b = 1#0.005
    c = 1#0.00225

    #scale = 0.5
    #a = scale*a
    #b = scale*b
    #c = scale*c

    e1 = 0.5 # theta roundness
    e2 = 0.16667 # phi roundness
    n1 = 2/e1 # Blockiness parameter
    n2 = 2/e2 # Blockiness 

    superquad = Superquadric(a=a,b=b,c=c,e1=e1,e2=e2,n1=n1,n2=n2)

    #plot_superquadric_mayavi(superquad,n=8)
    plot_grid_mayavi(superquad, sc, xSide, ySide, zSide)


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
