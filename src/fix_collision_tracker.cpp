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
  pre_particles_were_in_contact_offset = 1; // = pair_gran->add_history_value("pre_particles_were_in_contact", "0");

  
  ncollisions = 0;

  size_local_rows = 0;
  size_local_cols = 0;

//  memory->create(vector_local, DELTA, "fix/collision/tracker");
//  vector_local_size = DELTA;

  local_freq = 1;
  local_flag = 1;
  nevery = atoi(arg[3]);

  // We will add options later
  time_step_counter = 0;

  MPI_Comm_rank(world,&me); 

  fp = NULL;
  //char *title = NULL;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix print command");
      if (me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == NULL) {
          char str[512];
          sprintf(str,"Cannot open fix print file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
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
  memory->destroy(vector_local);
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

void FixCollisionTracker::init()
{
  // Don't know why this did not work in the constructor
  pre_particles_were_in_contact_offset = pair_gran->add_history_value("pre_particles_were_in_contact", "0");
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::end_of_step()
{

  // print to file
  // make a copy of string to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

//  modify->clearstep_compute();

//  strcpy(copy,string);
//  input->substitute(copy,work,maxcopy,maxwork,0);

//  modify->addstep_compute(update->ntimestep + nevery);

  if (me == 0 && fp) {

    std::vector<double>::iterator rel_it = rel_vels.begin();
    std::vector<double>::iterator lcol_it = lcol.begin();

    while (rel_it < rel_vels.end())
      fprintf(fp,"%f, %f, %f, %f\n", rel_it++, lcol_it++, lcol_it++, lcol_it++);
   
    fflush(fp);
    rel_vels.clear();
    lcol.clear();
  }

}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::post_force(int inumber)
{
  //time_step_counter++;
  //printf("timestep %d\n", time_step_counter);
  
  // pair_gran_base.h in compute_force()

  PairGran *pg = pair_gran;
  double ** first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;
  int ** firstneigh = pg->list->firstneigh;

  int inum = pg->list->inum;
  int * ilist = pg->list->ilist;
  const int dnum = pg->dnum();
  int * numneigh = pg->list->numneigh;
  double **v = atom->v;
  double **f = atom->f;
  
  SurfacesIntersectData sidata;
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
      print_contact_status(sidata);  
      print_atom_pair_info(i,j);
    }
  }

/*
  double **velocity = atom->v;
  int n = atom->nlocal;

  double speed[n];

  InternalValue = 0;
  fprintf(screen , "Current velocities are: ");
  for (int i = 0; i < n; ++i)
  {
    speed[i] = sqrt(velocity[i][0]*velocity[i][0] + 
                    velocity[i][1]*velocity[i][1] + 
                    velocity[i][2]*velocity[i][2]); 
   
    InternalValue += speed[i]; 
    fprintf(screen , "%f, ", speed[i]);
  } 
  fprintf(screen, "\n");

  // fprintf(screen , "current iteration value is %f\n", InternalValue);
*/
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::print_atom_pair_info(int i, int j)
{
  double **x = atom->x;
  double **v = atom->v;
  //double **f = atom->f;
  printf("Atom i[%d]: x[%f,%f,%f]; Atom j[%d]: x[%f,%f,%f]\n", i,x[i][0],x[i][1],x[i][2],j,x[j][0],x[j][1],x[j][2]);
  printf("Atom i[%d]: v[%f,%f,%f]; Atom j[%d]: v[%f,%f,%f]\n", i,v[i][0],v[i][1],v[i][2],j,v[j][0],v[j][1],v[j][2]);

}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::print_contact_status(SurfacesIntersectData& sidata) //, IContactHistorySetup* hsetup)
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
  bool intersect = *particles_were_in_contact == SURFACES_INTERSECT;
  bool pre_intersect = *pre_particles_were_in_contact == SURFACES_INTERSECT;
 
  if (intersect && !pre_intersect)
  {
    ++ncollisions;
    compute_normal(sidata);
    double rel_v = compute_relative_velocity(sidata);

    rel_vels.push_back(rel_v);
    /*
    vector_local[size_local_rows++] = rel_v;
    if (size_local_rows == vector_local_size) 
    {
      vector_local_size += DELTA;
      memory->grow(vector_local, vector_local_size, "fix/collision/tracker");
    }
    */
    printf("The relative velocity of the particles at the time of impact was: %f \n", rel_v);
  }
  
//  bool intersect = checkSurfaceIntersect(sidata);
  printf("Particles(%d,%d) intersect %d\n", sidata.i, sidata.j, intersect);
  printf("Intersection point: %f,%f,%f\n", prev_step_point[0], prev_step_point[1], prev_step_point[2]);  
  printf("Total number of collisions: %i\n", ncollisions);  


  *pre_particles_were_in_contact = *particles_were_in_contact;
} 

/* ---------------------------------------------------------------------- */

double FixCollisionTracker::compute_scalar()
{
  return InternalValue;
}

/* ---------------------------------------------------------------------- */

double FixCollisionTracker::compute_relative_velocity(SurfacesIntersectData& sidata)
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
  return vectorDot3D(rel_v, sidata.en);
}

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


bool FixCollisionTracker::checkSurfaceIntersect(SurfacesIntersectData & sidata)
    {
      sidata.is_non_spherical = true;
      bool particles_in_contact = false;
      //double *const prev_step_point = &sidata.contact_history[contact_point_offset]; //contact points
      //double *const inequality_start = &sidata.contact_history[inequality_start_offset];
      //double *const particles_were_in_contact = &sidata.contact_history[particles_were_in_contact_offset];

      const int iPart = sidata.i;
      const int jPart = sidata.j;

      #ifdef LIGGGHTS_DEBUG
        if(std::isnan(vectorMag3D(atom->x[iPart])))
          error->one(FLERR,"atom->x[iPart] is NaN!");
        if(std::isnan(vectorMag4D(atom->quaternion[iPart])))
          error->one(FLERR,"atom->quaternion[iPart] is NaN!");
        if(std::isnan(vectorMag3D(atom->x[jPart])))
          error->one(FLERR,"atom->x[jPart] is NaN!");
        if(std::isnan(vectorMag4D(atom->quaternion[jPart])))
          error->one(FLERR,"atom->quaternion[jPart] is NaN!");
      #endif

      particle_i.set(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
      particle_j.set(atom->x[jPart], atom->quaternion[jPart], atom->shape[jPart], atom->blockiness[jPart]);

      /*
      unsigned int int_inequality_start = MathExtraLiggghtsNonspherical::round_int(*inequality_start);
      bool obb_intersect = false;
      if(*particles_were_in_contact == SURFACES_INTERSECT)
        obb_intersect = true; //particles had overlap on the previous time step, skipping OBB intersection check
      else
        obb_intersect = MathExtraLiggghtsNonspherical::obb_intersect(&particle_i, &particle_j, int_inequality_start);
      */
      bool obb_intersect = false;
      obb_intersect = MathExtraLiggghtsNonspherical::obb_intersect(&particle_i, &particle_j);
      if(obb_intersect) {//OBB intersect particles in possible contact

        double fi, fj;
        const double ri = cbrt(particle_i.shape[0]*particle_i.shape[1]*particle_i.shape[2]);
        const double rj = cbrt(particle_j.shape[0]*particle_j.shape[1]*particle_j.shape[2]);
        double ratio = ri / (ri + rj);

        /*
        if(*particles_were_in_contact == SURFACES_FAR)
          MathExtraLiggghtsNonspherical::calc_contact_point_if_no_previous_point_avaialable(sidata, &particle_i, &particle_j, sidata.contact_point, fi, fj, this->error);
        else
          MathExtraLiggghtsNonspherical::calc_contact_point_using_prev_step(sidata, &particle_i, &particle_j, ratio, update->dt, prev_step_point, sidata.contact_point, fi, fj, this->error);
        */
        MathExtraLiggghtsNonspherical::calc_contact_point_if_no_previous_point_avaialable(sidata, &particle_i, &particle_j, sidata.contact_point, fi, fj, this->error);
        //LAMMPS_NS::vectorCopy3D(sidata.contact_point, prev_step_point); //store contact point in contact history for the next DEM time step

        #ifdef LIGGGHTS_DEBUG
          if(std::isnan(vectorMag3D(sidata.contact_point)))
            error->one(FLERR,"sidata.contact_point is NaN!");
        #endif
        
        particles_in_contact = std::max(fi, fj) < 0.0;

        if(particles_in_contact) {
          //*particles_were_in_contact = SURFACES_INTERSECT;
        } 
        else
        {
          //*particles_were_in_contact = SURFACES_CLOSE;
        }
       } 
       else
       {
         //*particles_were_in_contact = SURFACES_FAR;
       }
      //*inequality_start = static_cast<double>(int_inequality_start);
      return particles_in_contact;
    }
