/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */



#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "atom.h"
#include "pair.h"
#include "update.h"
#include "error.h"
#include "fix_collision_tracker.h"
#include "pair_gran.h"
#include "contact_models.h"
#include "force.h"

#include "neigh_list.h"

#include "pointers.h"
#include "lammps.h"
#include "contact_interface.h"
#include "property_registry.h"
#include "settings.h"
#include "contact_model_constants.h"
#include "contact_model_base.h"
#include "surface_model_base.h"
#include "normal_model_base.h"
#include "tangential_model_base.h"
#include "rolling_model_base.h"
#include "cohesion_model_base.h"
#include "style_surface_model.h"
#include "style_normal_model.h"
#include "style_tangential_model.h"
#include "style_rolling_model.h"
#include "style_cohesion_model.h"

#include "fix.h"
#include "fix_insert.h"
#include "fix_insert_stream.h"
#include "fix_insert_stream_predefined.h"
#include "fix_wall_gran.h"
#include "fix_contact_history_mesh.h"

#include "math_extra_liggghts_nonspherical.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LIGGGHTS::ContactModels;

/* ---------------------------------------------------------------------- */

FixCollisionTracker::FixCollisionTracker(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix collision/tracker command, not enough arguments");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 1));
  particles_were_in_contact_offset = pair_gran->get_history_offset("particles_were_in_contact", "0");
  contact_point_offset = pair_gran->get_history_offset("cpx", "0");
  pre_particles_were_in_contact_offset = pair_gran->get_history_offset("pre_particles_were_in_contact", "0"); // = pair_gran->add_history_value("pre_particles_were_in_contact", "0");

  // local array setup
  size_local_rows = 0;
  size_local_cols = 2;

  int memsize = 2 * sizeof(double*);
  array_local = (double**)memory->smalloc(memsize, "fix/collision/tracker");

  local_freq = 1;
  local_flag = 1;
  nevery = atoi(arg[3]);
  global_freq = 1; // ?
  array_offset = 0;
  
  // Cube projection parameters and 2d arrays to 0;
  cube_projection = 0;
  x_nsplit = 0;
  y_nsplit = 0;
  z_nsplit = 0;
  x_octsurface = NULL;
  y_octsurface = NULL;
  z_octsurface = NULL;
  
  MPI_Comm_rank(world,&me); 
  fp = NULL;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      if (1){ //me == 0) {
        int n = strlen(arg[iarg+1]) + 1;
        filename = new char[n];
        strcpy(filename,arg[iarg+1]);
        writetofile = 1;
      }
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"cpoctant") == 0)
    {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix cpoctant command");

      cube_projection = 1;

      // x,y,z sides
      x_nsplit = atoi(arg[iarg+1]);
      if (x_nsplit < 1) error->all(FLERR,"x axis split < 1");

      y_nsplit = atoi(arg[iarg+2]);
      if (y_nsplit < 1) error->all(FLERR,"y axis split < 1");
      
      z_nsplit = atoi(arg[iarg+3]);
      if (z_nsplit < 1) error->all(FLERR,"z axis split < 1");

      memory->create(x_octsurface, y_nsplit, z_nsplit, "xprojection");
      memory->create(y_octsurface, x_nsplit, z_nsplit, "yprojection");
      memory->create(z_octsurface, x_nsplit, y_nsplit, "zprojection");

      // Zero out newly allocated array
      // Using a trick based on how memory->creates allocates memory
      memset(x_octsurface[0], 0, sizeof(**x_octsurface)*y_nsplit*z_nsplit);
      memset(y_octsurface[0], 0, sizeof(**y_octsurface)*x_nsplit*z_nsplit);
      memset(z_octsurface[0], 0, sizeof(**z_octsurface)*x_nsplit*y_nsplit);

      iarg += 4;
    }
  }

}

/* ---------------------------------------------------------------------- */

FixCollisionTracker::~FixCollisionTracker()
{
  memory->sfree(array_local);
  memory->destroy(x_octsurface);
  memory->destroy(y_octsurface);
  memory->destroy(z_octsurface);
}

/* ---------------------------------------------------------------------- */

int FixCollisionTracker::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::end_of_step()
{
  if (!writetofile)
    return;
  
  openfile();
  if (!fp)
  {
    // log error?
    return;
  }

  if(cube_projection)
  {
    fprintf(fp,"%d, %d, %d\n", x_nsplit, y_nsplit, z_nsplit);
    print_cube_projection(fp);
  }
  else
  {
    double * rel = rel_vels.data();
    double * col = lcol.data();
    int n = rel_vels.size() / 2;

    for (int i = 0; i < n; ++i)
      fprintf(fp,"%f %f %f %f %f\n", rel[2*i], rel[2*i+1], col[3*i], col[3*i+1], col[3*i+2]);
  }
  fflush(fp);
  fclose(fp);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::openfile()
{
  // if one file per timestep, replace '*' with current timestep  
  char *filecurrent = filename;
  
  char *filestar = filecurrent;
  filecurrent = new char[strlen(filestar) + 16 + 4];
  char *ptr = strchr(filestar,'*');
  if (ptr) {
    *ptr = '\0';
    sprintf(filecurrent,"%s" BIGINT_FORMAT "_%i" "%s",
            filestar,update->ntimestep, me, ptr+1);
    *ptr = '*';
      
    fp = fopen(filecurrent,"w");
  } else {
    sprintf(filecurrent,"%i_%s", me, filestar);
    fp = fopen(filename, "a");
  }

  if (fp == NULL) error->one(FLERR,"Cannot open dump file");

  delete [] filecurrent;
}

void FixCollisionTracker::print_cube_projection(FILE *fp)
{
  // Print x side
  for(int y = 0; y < y_nsplit; y++)
  {
    fprintf(fp,"%d", x_octsurface[y][0]);
    for(int z = 1; z < z_nsplit; z++)
    {
      fprintf(fp,",%d", x_octsurface[y][z]);
    }
    fprintf(fp,"\n");
  }

  // Print y side
  for(int x = 0; x < x_nsplit; x++)
  {
    fprintf(fp,"%d", y_octsurface[x][0]);
    for(int z = 1; z < z_nsplit; z++)
    {
      fprintf(fp,",%d", y_octsurface[x][z]);
    }
    fprintf(fp,"\n");
  }

  // Print z side
  for(int x = 0; x < x_nsplit; x++)
  {
    fprintf(fp,"%d", z_octsurface[x][0]);
    for(int y = 1; y < y_nsplit; y++)
    {
      fprintf(fp,",%d", z_octsurface[x][y]);
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::post_force(int vflag)
{
  if ((update->ntimestep - 1) % nevery == 0)
  {
    array_offset = 0;
    rel_vels.clear();
    lcol.clear();
  }
  // Get local atom information
  int *mask = atom->mask;

  PairGran *pg = pair_gran;
  double ** first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;
  int ** firstneigh = pg->list->firstneigh;

  int inum = pg->list->inum;
  int * ilist = pg->list->ilist;
  const int dnum = pg->dnum();
  int * numneigh = pg->list->numneigh;
  
  SurfacesIntersectData sidata;

  // particle - particle
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    sidata.i = i;
  
    double * const all_contact_hist = first_contact_hist ? first_contact_hist[i] : NULL;

    int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      sidata.j = j;

      sidata.contact_history = all_contact_hist ? &all_contact_hist[dnum*jj] : NULL;
     
      resolve_contact_status(sidata);
    }
  }

  // particle - wall (mesh is included) 
  //...

  // This is as weird as it looks, but LIGGGHTS want their array this way
  array_local[0] = rel_vels.data() + array_offset;
  size_local_rows = (rel_vels.size() - array_offset) / 2;
  array_offset += size_local_rows * 2;
}

/* ---------------------------------------------------------------------- */

// Should only be called once
inline bool FixCollisionTracker::check_collision(SurfacesIntersectData& sidata)
{
  double *const particles_were_in_contact = &sidata.contact_history[particles_were_in_contact_offset];
  double *const pre_particles_were_in_contact = &sidata.contact_history[pre_particles_were_in_contact_offset];

  bool pre_intersect = *pre_particles_were_in_contact == SURFACES_INTERSECT;
  bool intersect = *particles_were_in_contact == SURFACES_INTERSECT;
 
  *pre_particles_were_in_contact = *particles_were_in_contact; 
  
  return intersect && !pre_intersect;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::resolve_contact_status(SurfacesIntersectData& sidata)
{
  if (!check_collision(sidata))
    return;
    
  bool jislocal = sidata.j < atom->nlocal;

  compute_normal(sidata);
  double vels[2]; 
  compute_relative_velocity(sidata, vels, vels+1);

  double point_of_contact[6];
  compute_local_contact(sidata, point_of_contact, point_of_contact+3);

  rel_vels.push_back(vels[0]);
  rel_vels.push_back(vels[1]);
  for (int i = 0; i < 3; ++i)
    lcol.push_back(point_of_contact[i]);

  if (jislocal)
  {
    rel_vels.push_back(vels[0]);
    rel_vels.push_back(vels[1]); 
    for (int i = 3; i < 6; ++i)
      lcol.push_back(point_of_contact[i]);
  }

  // Projecting contact point to octant of unit cube
  if(!cube_projection)
    return;

  double result[3];
  unit_cube_oct_projection(sidata.i, point_of_contact, result);
  unit_cube_oct_indexing(result);
  if (jislocal)
  {
    unit_cube_oct_projection(sidata.j, point_of_contact+3, result);
    unit_cube_oct_indexing(result);
  }
} 

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::compute_relative_velocity(SurfacesIntersectData& sidata, double* normal, double* tangent)
{
  double *const prev_step_point = &sidata.contact_history[contact_point_offset];
  int iPart = sidata.i;
  int jPart = sidata.j;

  // v1 = v + cross(w,(p-x))
  double a_i[3]; 
  double b_i[3];
  double a_j[3]; 
  double b_j[3];
  vectorSubtract3D(prev_step_point, atom->x[iPart], a_i);
  vectorSubtract3D(prev_step_point, atom->x[jPart], a_j);
  vectorCross3D(atom->omega[iPart], a_i, b_i);
  vectorCross3D(atom->omega[jPart], a_j, b_j);
  vectorAdd3D(atom->v[iPart], b_i, a_i);
  vectorAdd3D(atom->v[jPart], b_j, a_j);

  // relv = dot((v1-v2),sidata.en)
  double rel_v[3];
  vectorSubtract3D(a_i, a_j, rel_v);
  double res = vectorDot3D(rel_v, sidata.en);
  
  *normal = res > 0 ? res : -res;
  *tangent = sqrt(vectorDot3D(rel_v, rel_v) - res * res);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::compute_normal(SurfacesIntersectData& sidata)
{
  double *const prev_step_point = &sidata.contact_history[contact_point_offset];
  
  int iPart = sidata.i;
  int jPart = sidata.j;

  Superquadric particle_i;
  Superquadric particle_j;
  particle_i.set(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
  particle_j.set(atom->x[jPart], atom->quaternion[jPart], atom->shape[jPart], atom->blockiness[jPart]);

  particle_i.shape_function_gradient_global(prev_step_point, particle_i.gradient);
  particle_j.shape_function_gradient_global(prev_step_point, particle_j.gradient);

  vectorSubtract3D(particle_i.gradient,particle_j.gradient, sidata.en);
  vectorNormalize3D(sidata.en);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::compute_local_contact(SurfacesIntersectData& sidata, double *iResult, double *jResult)
{
  double *const prev_step_point = &sidata.contact_history[contact_point_offset];
  int iPart = sidata.i;
  int jPart = sidata.j;

  double iLocal[3];
  double jLocal[3];

  vectorSubtract3D(prev_step_point, atom->x[iPart], iLocal);
  vectorSubtract3D(prev_step_point, atom->x[jPart], jLocal);

  MathExtraLiggghtsNonspherical::rotate_global2local(atom->quaternion[iPart], iLocal, iResult);
  MathExtraLiggghtsNonspherical::rotate_global2local(atom->quaternion[jPart], jLocal, jResult);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::unit_cube_oct_projection(int iPart, double *contact, double *result)
{
  double* shape = atom->shape[iPart];

  vectorCopy3D(contact, result);
  vectorAbs3D(result);

  // Normalise superquadric particle shape
  result[0] /= shape[0];
  result[1] /= shape[1];
  result[2] /= shape[2];
  
  double maxVal = vectorMax3D(result);
  vectorScalarDiv3D(result, maxVal);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::unit_cube_oct_indexing(double *cube_proj)
{
  if (cube_proj[0] >= cube_proj[1] && cube_proj[0] >= cube_proj[2]) // Project to x side
  {
    //calculating indexes
    int y = std::min(y_nsplit-1, (int)floor(y_nsplit*cube_proj[1]));
    int z = std::min(z_nsplit-1, (int)floor(z_nsplit*cube_proj[2]));
    x_octsurface[y][z]++; //(y,z)
  }
  else if(cube_proj[1] >= cube_proj[2]) // Project to y side
  {
    //calculating indexes
    int x = std::min(x_nsplit-1, (int)floor(x_nsplit*cube_proj[0]));
    int z = std::min(z_nsplit-1, (int)floor(z_nsplit*cube_proj[2]));
    y_octsurface[x][z]++; //(x,z)
  }
  else // Project to z side
  {
    //calculating indexes
    int x = std::min(x_nsplit-1, (int)floor(x_nsplit*cube_proj[0]));
    int y = std::min(y_nsplit-1, (int)floor(y_nsplit*cube_proj[1]));
    z_octsurface[x][y]++; //(x,y)
  } 
}

/* ---------------------------------------------------------------------- */

double* FixCollisionTracker::get_triangle_contact_history(TriMesh *mesh, FixContactHistoryMesh *fix_contact, int iPart, int iTri)
{
  // get contact history of particle iPart and triangle idTri
  // NOTE: depends on naming in fix_wall_gran!

  std::string fix_nneighs_name("n_neighs_mesh_");
  fix_nneighs_name += mesh->mesh_id();
  FixPropertyAtom* fix_nneighs = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name.c_str(),"property/atom","scalar",0,0,this->style));

  int idTri = mesh->id(iTri);
  const int nneighs = fix_nneighs->get_vector_atom_int(iPart);
  for (int j = 0; j < nneighs; ++j)
  {
    if (fix_contact->partner(iPart, j) == idTri)
    {
      return fix_contact->contacthistory(iPart, j);
    }
  }
  return NULL;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::print_atom_pair_info(int i, int j)
{
  double **x = atom->x;
  double **v = atom->v;
  printf("Atom i[%d]: x[%f,%f,%f]; Atom j[%d]: x[%f,%f,%f]\n", i,x[i][0],x[i][1],x[i][2],j,x[j][0],x[j][1],x[j][2]);
  printf("Atom i[%d]: v[%f,%f,%f]; Atom j[%d]: v[%f,%f,%f]\n", i,v[i][0],v[i][1],v[i][2],j,v[j][0],v[j][1],v[j][2]);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::print_atom_info(int i)
{
  double **x = atom->x;
  double **v = atom->v;
  printf("Atom i[%d]: x[%f,%f,%f]; v[%f,%f,%f]\n", i,x[i][0],x[i][1],x[i][2],v[i][0],v[i][1],v[i][2]);
}

/* ---------------------------------------------------------------------- */
