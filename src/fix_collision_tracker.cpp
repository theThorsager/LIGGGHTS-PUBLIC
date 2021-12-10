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

#define DELTA 1000

/* ---------------------------------------------------------------------- */

FixCollisionTracker::FixCollisionTracker(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix collision/tracker command, not enough arguments");

  // int iarg = 6;

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 1));
  particles_were_in_contact_offset = pair_gran->get_history_offset("particles_were_in_contact", "0");
  contact_point_offset = pair_gran->get_history_offset("cpx", "0");
  pre_particles_were_in_contact_offset = pair_gran->get_history_offset("pre_particles_were_in_contact", "0"); // = pair_gran->add_history_value("pre_particles_were_in_contact", "0");

  
  ncollisions = 0;

  size_local_rows = 0;
  size_local_cols = 2;

  int memsize = 2 * sizeof(double*);
  array_local = (double**)memory->smalloc(memsize, "fix/collision/tracker");
//  memory->create(array_local, memsize, "fix/collision/tracker");
//  vector_local_size = DELTA;

  local_freq = 1;
  local_flag = 1;
  nevery = atoi(arg[3]);
  global_freq = 1; // ?
  array_offset = 0;
  
  // We will add options later
  time_step_counter = 0;

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
    } else
    {
        error->all(FLERR,"Illegal fix print command");
    }
  }

}

FixCollisionTracker::~FixCollisionTracker()
{
  memory->sfree(array_local);
}
/* ---------------------------------------------------------------------- */

int FixCollisionTracker::setmask()
{
  int mask = 0;
//  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::init()
{
  // Don't know why this did not work in the constructor
 // pre_particles_were_in_contact_offset = pair_gran->add_history_value("pre_particles_were_in_contact", "0");
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::end_of_step()
{
  if (writetofile){ //me == 0) {
    openfile();
    if (fp)
    {
      double * rel = rel_vels.data();
      double * col = lcol.data();
      int n = rel_vels.size() / 2;

      for (int i = 0; i < n; ++i)
        fprintf(fp,"%f %f %f %f %f\n", rel[2*i], rel[2*i+1], col[3*i], col[3*i+1], col[3*i+2]);
   
      fflush(fp);
      fclose(fp);
    }
  }
//  printf("removing data at: %i\n", update->ntimestep);
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

  // delete string with timestep replaced

  delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */
/*
void FixCollisionTracker::pre_force(int vflag)
{
  PairGran *pg = pair_gran;
  double ** first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;
  int ** firstneigh = pg->list->firstneigh;

  int inum = pg->list->inum;
  int * ilist = pg->list->ilist;
  const int dnum = pg->dnum();
  int * numneigh = pg->list->numneigh;
  
  SurfacesIntersectData sidata;
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
  
    double * const all_contact_hist = first_contact_hist ? first_contact_hist[i] : NULL;

    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      sidata.contact_history = all_contact_hist ? &all_contact_hist[dnum*jj] : NULL;
      double *const prething = &sidata.contact_history[particles_were_in_contact_offset];
      prev_intersections.push_back(prething);
    } 
  }
}
*/
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
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *tag = atom->tag;

  // pair_gran_base.h in compute_force()
  // Get granular pair information and neighboor information
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
     
      // if (collision)
      //        compute impact velocity and add to bin
      //        compute impact angle and add to bin
      //
      print_contact_status(sidata);//, iter++);  
      print_atom_pair_info(i,j);
    }
  }

  // particle - wall (mesh is included)
  int n_wall_fixes = modify->n_fixes_style("wall/gran");

  for (int ifix = 0; ifix < n_wall_fixes; ++ifix)
  {
    FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));

    if (fwg->is_mesh_wall())
    {
      int n_FixMesh = fwg->n_meshes();

      for (int iMesh = 0; iMesh < n_FixMesh; iMesh++)
      {
        TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
        int nTriAll = mesh->sizeLocal() + mesh->sizeGhost(); // Need to check how to handle ghost particles
        FixContactHistoryMesh *fix_contact = fwg->mesh_list()[iMesh]->contactHistory();
        if (!fix_contact) continue;

        // get neighborList and numNeigh
        FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();
        if (!meshNeighlist) continue;

        // loop owned trinagles //(and ghost triangles)
        for (int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for (int iCont = 0; iCont < numneigh; iCont++) {

            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if (iPart >= nlocal) continue;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
          
            double *contact_history = get_triangle_contact_history(mesh, fix_contact, iPart, iTri);
            if (contact_history)
            {
            //  printf("triangle: %d; particle: %d; \n", iTri, iPart);
            //  print_atom_info(iPart);
            }
          
          }
        }
      }
    }
    else // Is a primitive wall
    {
      /*
      double **c_history = get_primitive_wall_contact_history(fwg);
      if (c_history)
      {
        // loop neighbor list
        int *neighborList;
        int nNeigh = fwg->primitiveWall()->getNeighbors(neighborList);

        for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++)
        {
          int iPart = *neighborList;
          // do not need to handle ghost particles
          if (iPart >= nlocal) continue;
          if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;

          double *contact_history = c_history[iPart];

          if (contact_history)
          {
            // Do stuff
          }
        }
      }
      */
    }
  }

  // This is as weird as it looks, but LIGGGHTS want their array this way
  array_local[0] = rel_vels.data() + array_offset;
  size_local_rows = (rel_vels.size() - array_offset) / 2;
  array_offset += size_local_rows * 2;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::print_atom_pair_info(int i, int j)
{
  double **x = atom->x;
  double **v = atom->v;
  //double **f = atom->f;
 // printf("Atom i[%d]: x[%f,%f,%f]; Atom j[%d]: x[%f,%f,%f]\n", i,x[i][0],x[i][1],x[i][2],j,x[j][0],x[j][1],x[j][2]);
 // printf("Atom i[%d]: v[%f,%f,%f]; Atom j[%d]: v[%f,%f,%f]\n", i,v[i][0],v[i][1],v[i][2],j,v[j][0],v[j][1],v[j][2]);

}

void FixCollisionTracker::print_atom_info(int i)
{
  double **x = atom->x;
  double **v = atom->v;
  //double **f = atom->f;
  printf("Atom i[%d]: x[%f,%f,%f]; v[%f,%f,%f]\n", i,x[i][0],x[i][1],x[i][2],v[i][0],v[i][1],v[i][2]);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::print_contact_status(SurfacesIntersectData& sidata)//, int iter) //, IContactHistorySetup* hsetup)
{
  double *const particles_were_in_contact = &sidata.contact_history[particles_were_in_contact_offset];
  double *const pre_particles_were_in_contact = &sidata.contact_history[pre_particles_were_in_contact_offset];
  double *const prev_step_point = &sidata.contact_history[contact_point_offset];
//  fprintf(screen, "The offset is: %i\n", particles_were_in_contact_offset);
//  fprintf(screen , "if statment %i\n", *particles_were_in_contact == SURFACES_INTERSECT);
  
  /*
  if (*particles_were_in_contact == SURFACES_INTERSECT) {
    fprintf(screen , "Particles %i and %i are in contact\n", sidata.i, sidata.j);
  } else {
    fprintf(screen , "Particles %i and %i are not in contact\n", sidata.i, sidata.j);
  }
  */
  bool pre_intersect = *pre_particles_were_in_contact == SURFACES_INTERSECT;
  bool intersect = *particles_were_in_contact == SURFACES_INTERSECT;

//  printf("Particles were intersecting: %i, are: %i\n", pre_intersect, intersect); 
  if (intersect && !pre_intersect)
  {
    ++ncollisions;
    compute_normal(sidata);
    compute_relative_velocity(sidata);

//    rel_vels.push_back(rel_v);
//    rel_vels.push_back(rel_v);

    double point_of_contact[6];
    compute_local_contact(sidata, point_of_contact, point_of_contact+3);
    int nind = 6;
    if (sidata.j < atom->nlocal)
      nind = 3;

    for (int i = 0; i < nind; ++i)
      lcol.push_back(point_of_contact[i]);
    /*
    vector_local[size_local_rows++] = rel_v;
    if (size_local_rows == vector_local_size) 
    {
      vector_local_size += DELTA;
      memory->grow(vector_local, vector_local_size, "fix/collision/tracker");
    }
    */
   // printf("The relative velocity of the particles at the time of impact was: %f \n", rel_v);
  }
  
//  bool intersect = checkSurfaceIntersect(sidata);
//  printf("Particles(%d,%d) intersect %d\n", sidata.i, sidata.j, intersect);
//  printf("Intersection point: %f,%f,%f\n", prev_step_point[0], prev_step_point[1], prev_step_point[2]);  
//  printf("Total number of collisions: %i\n", ncollisions); 

  *pre_particles_were_in_contact = *particles_were_in_contact; 
} 

/* ---------------------------------------------------------------------- */

double FixCollisionTracker::compute_scalar()
{
  return InternalValue;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::compute_relative_velocity(SurfacesIntersectData& sidata)
{
  double *const prev_step_point = &sidata.contact_history[contact_point_offset];
  int iPart = sidata.i;
  int jPart = sidata.j;

  // v1 = v + cross(w,(p-x))
  double a_i[3]; 
  double b_i[3];
  vectorSubtract3D(prev_step_point, atom->x[iPart], a_i);
  vectorCross3D(atom->omega[iPart], a_i, b_i);
  vectorAdd3D(atom->v[iPart], b_i, a_i);
  // ...
  double a_j[3]; 
  double b_j[3];
  vectorSubtract3D(prev_step_point, atom->x[jPart], a_j);
  vectorCross3D(atom->omega[jPart], a_j, b_j);
  vectorAdd3D(atom->v[jPart], b_j, a_j);

  // relv = dot((v1-v2),sidata.en)
  double rel_v[3];
  vectorSubtract3D(a_i, a_j, rel_v);
  double res = vectorDot3D(rel_v, sidata.en);
  double rel_vel = res > 0 ? res : -res;

  double tan_res = sqrt(vectorDot3D(rel_v, rel_v) - res * res);

  rel_vels.push_back(rel_vel);
  rel_vels.push_back(tan_res);

  if (jPart < atom->nlocal) {
    rel_vels.push_back(rel_vel);
    rel_vels.push_back(tan_res);
  }

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
