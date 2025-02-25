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
from enum import Enum

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

class Area(Enum):
    none = 'none'
    norm = 'norm'
    si = 'SI' # SI units

    def __str__(self):
        return self.value


def get_args():
    parser = argparse.ArgumentParser(description='Create cube-spherical mesh')
    parser.add_argument('-f', '--file', dest='file', type=pathlib.Path,
                        required=True, help='file path')
    # Frequency per particle per time per square-area
    # If parameter is left out then default value is used instead
    parser.add_argument('--particle-num', dest='part_num', type=int, default=1,
                        help='Number of particles observed in simulation')
    parser.add_argument('--time', dest='time', type=float, default=1,
                        help='Duration of observation in simulation, measured in seconds')
    parser.add_argument('--area', dest='area', type=Area, choices=list(Area), default=Area.none,
                        help='How surface area is accounted for in frequency calculations')

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

    #superquad_surf = mlab.mesh(x, y, z, colormap='Blues')
    superquad_surf = mlab.mesh(x, y, z, colormap='Blues',representation='wireframe')

    mlab.view(.0, -5.0, 4)
    mlab.show()

    #return fig


def prep_scalar_value(side):
    # Each square consists of 2 triangles
    #scalar = np.transpose(side)
    scalar = np.array(side).flatten()
    scalar = np.append(scalar,scalar)

    return scalar


def triangle_area(p1,p2,p3):
    a = p1-p2
    b = p3-p2
    area = np.cross(a,b)
    return np.linalg.norm(area)/2


def meshface_area(xMesh,yMesh,zMesh):
    N,M = np.shape(xMesh)
    w = []
    for i in range(N-1):
        wTemp = []
        for j in range(M-1):
            p1 = np.array([xMesh[i][j],yMesh[i][j],zMesh[i][j]])
            p2 = np.array([xMesh[i+1][j],yMesh[i+1][j],zMesh[i+1][j]])
            p3 = np.array([xMesh[i][j+1],yMesh[i][j+1],zMesh[i][j+1]])
            p4 = np.array([xMesh[i+1][j+1],yMesh[i+1][j+1],zMesh[i+1][j+1]])
            wTemp.append(triangle_area(p1,p2,p3)+triangle_area(p4,p2,p3))
        w.append(wTemp)
    return np.array(w)


def create_meshgrid(u, uN, v, vN):
    u = np.linspace(-u, u, uN)
    v = np.linspace(-v, v, vN)
    u,v = np.meshgrid(u,v)
    return u,v

def plot_cell_data(mesh, scalars, name='none'):
    mesh.mlab_source.dataset.cell_data.scalars = scalars
    mesh.mlab_source.dataset.cell_data.scalars.name = name
    mesh.mlab_source.dataset.cell_data.update()
    mesh.actor.mapper.scalar_mode = 'use_cell_data'
    mesh2 = mlab.pipeline.set_active_attribute(mesh,cell_scalars=name)
    return mesh2


def plot_grid_mayavi(superquadric, spherecube, xSide, ySide, zSide, arg_dict={}):

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

    fig = mlab.figure('Collision frequency', fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

    ax_scale = [1.0, 1.0, 1.0] # Need to change later
    ax_ranges = [-2, 2, -2, 2, -2, 2]
    ax_extent = ax_ranges * np.repeat(ax_scale, 2)
    vmin = 0
    vmax = np.max([np.max(np.array(xSide).flatten()),np.max(np.array(ySide).flatten()),np.max(np.array(zSide).flatten())])


    # Create meshgrid for side X,-X
    u = np.linspace(-b, b, yN)
    v = np.linspace(-c, c, zN)
    u,v = np.meshgrid(u,v)

    # Positive x side
    x = np.ones(np.shape(u))*a
    y = u
    z = v
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    x_pos = mlab.mesh(x, y, z, colormap='Blues')
    # Negative x side
    x = -x
    x_neg = mlab.mesh(x, y, z, colormap='Blues')
    x_area = meshface_area(x,y,z)

    # Create meshgrid for side Y,-Y
    u = np.linspace(-a, a, xN)
    v = np.linspace(-c, c, zN)
    u,v = np.meshgrid(u,v)

    # Positive y side
    x = u
    y = np.ones(np.shape(u))*b
    z = v
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    y_pos = mlab.mesh(x, y, z, colormap='Blues')
    # Negative x side
    y = -y
    y_neg = mlab.mesh(x, y, z, colormap='Blues')
    y_area = meshface_area(x,y,z)

    # Create meshgrid for side Z,-Z
    u = np.linspace(-a, a, xN)
    v = np.linspace(-b, b, yN)
    u,v = np.meshgrid(u,v)

    # Positive z side
    x = u
    y = v
    z = np.ones(np.shape(u))*c
    x,y,z = get_projected_coord(x,y,z,a,b,c,n1,n2)
    z_pos = mlab.mesh(x, y, z, colormap='Blues')
    # Negative x side
    z = -z
    z_neg = mlab.mesh(x, y, z, colormap='Blues')
    z_area = meshface_area(x,y,z)


    #area_normalisition = False
    area_handling = arg_dict.get('area')
    unit = None
    if area_handling == Area.norm:
        # Normalised in comparison to largest surface face
        max_area = np.max([np.max(x_area),np.max(y_area),np.max(z_area)])
        x_area = x_area/max_area
        y_area = y_area/max_area
        z_area = z_area/max_area

        x_w = prep_scalar_value(xSide/x_area)
        y_w = prep_scalar_value(ySide/y_area)
        z_w = prep_scalar_value(zSide/z_area)

        vmax = np.max([np.max(x_w),np.max(y_w),np.max(z_w)])

    elif area_handling == Area.si: 
        max_area = np.max([np.max(x_area),np.max(y_area),np.max(z_area)])
        min_area = np.min([np.min(x_area),np.min(y_area),np.min(z_area)])
        # Select appropriate prefix for surface area: mm^2, cm^2, m^2
        
        unit_coversion = 1
        if min_area < 0.000001: # 0.000001 m^2 == 1 mm^2
            unit_coversion = 1000000
            unit = '$mm^2$'
        elif min_area < 0.0001: # 0.0001 m^2 == 1 cm^2
            unit_coversion = 10000
            unit = '$cm^2$'
        else: # Use m^2
            unit_coversion = 1
            unit = '$m^2$'

        x_area = x_area*unit_coversion
        y_area = y_area*unit_coversion
        z_area = z_area*unit_coversion

        x_w = prep_scalar_value(xSide/x_area)
        y_w = prep_scalar_value(ySide/y_area)
        z_w = prep_scalar_value(zSide/z_area)

        vmax = np.max([np.max(x_w),np.max(y_w),np.max(z_w)])

    else:
        x_w = prep_scalar_value(xSide)
        y_w = prep_scalar_value(ySide)
        z_w = prep_scalar_value(zSide)


    # Calculate frequency per particle per time
    part_num = arg_dict.get('part_num',1)
    time = arg_dict.get('time',1)
    temp = part_num*time

    x_w = x_w/temp
    y_w = y_w/temp
    z_w = z_w/temp
    vmax = np.max([np.max(x_w),np.max(y_w),np.max(z_w)])


    # Plot heatmapped sides
    mesh2 = plot_cell_data(x_pos, x_w, name='X Cell data')
    surf = mlab.pipeline.surface(mesh2, vmin=vmin, vmax=vmax)
    mesh2 = plot_cell_data(x_neg, x_w, name='-X Cell data')
    surf = mlab.pipeline.surface(mesh2, vmin=vmin, vmax=vmax)
    
    mesh2 = plot_cell_data(y_pos, y_w, name='Y Cell data')
    surf = mlab.pipeline.surface(mesh2, vmin=vmin, vmax=vmax)
    mesh2 = plot_cell_data(y_neg, y_w, name='-Y Cell data')
    surf = mlab.pipeline.surface(mesh2, vmin=vmin, vmax=vmax)
    
    mesh2 = plot_cell_data(z_pos, z_w, name='Z Cell data')
    surf = mlab.pipeline.surface(mesh2, vmin=vmin, vmax=vmax)
    mesh2 = plot_cell_data(z_neg, z_w, name='-Z Cell data')
    surf = mlab.pipeline.surface(mesh2, vmin=vmin, vmax=vmax)

    if unit:
        title = 'Collisions \nper {} \nper second\n'.format(unit)
    else:
        title = 'Collisions \nper second\n'

    #mlab.colorbar(title='Phase', orientation='vertical', nb_labels=3)
    #mlab.colorbar(object=surf, title='Frequency', orientation='vertical', nb_labels=5)
    mlab.colorbar(object=surf, title=title, orientation='vertical')
    #mlab.colorbar()

    mlab.show()


def plot_superquadric_matplotlib(superquad, fig=None,ax=None,show=True,n=16):

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(projection='3d')

    # Account for start and end being the same in a circle
    n_samples = n + 1

    eta = np.linspace(-np.pi/2, np.pi/2, n_samples, endpoint=True)
    omega = np.linspace(-np.pi, np.pi, n_samples, endpoint=True)
    eta, omega = np.meshgrid(eta, omega)

    a = superquad.a
    b = superquad.b
    c = superquad.c

    e1 = superquad.e1
    e2 = superquad.e2
    n1 = superquad.n1
    n2 = superquad.n2

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


def plot_grid_matplotlib(superquad, fig=None,ax=None,show=True):

    if not fig:
        fig = plt.figure()
    if not ax:
        ax = fig.add_subplot(projection='3d')

    a = superquad.a
    b = superquad.b
    c = superquad.c

    e1 = superquad.e1
    e2 = superquad.e2
    n1 = superquad.n1
    n2 = superquad.n2

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


def get_projected_coords(x,y,z,superquad):
    return get_projected_coord(x,y,z,superquad.a,superquad.b,superquad.c,superquad.n1,superquad.n2)


def get_octant(path):
    with open(path, newline='') as csvfile:
        linereader = csv.reader(csvfile, delimiter=',')

        parameters = next(linereader)
        a,b,c,n1,n2 = [float(s) for s in parameters]

        superquad = Superquadric(a=a,b=b,c=c,e1=2/n1,e2=2/n2,n1=n1,n2=n2)

        parameters = next(linereader)
        xN,yN,zN = [int(s) for s in parameters]

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
        spherecube = Spherecube(2*xN,2*yN,2*zN)

        # Expand sides from octant to full  
        xSide = [*xSide[::-1], *xSide]
        xSide = [[*i[::-1],*i] for i in xSide]
        xSide = np.transpose(xSide)

        ySide = [*ySide[::-1], *ySide]
        ySide = [[*i[::-1],*i] for i in ySide]
        ySide = np.transpose(ySide)


        zSide = [*zSide[::-1], *zSide]
        zSide = [[*i[::-1],*i] for i in zSide]
        zSide = np.transpose(zSide)


        return spherecube, superquad, xSide, ySide, zSide


if __name__ == "__main__":

    args = get_args()
    sc, sq, xSide, ySide, zSide = get_octant(args.file)

    plot_grid_mayavi(sq, sc, xSide, ySide, zSide, arg_dict=vars(args))


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
