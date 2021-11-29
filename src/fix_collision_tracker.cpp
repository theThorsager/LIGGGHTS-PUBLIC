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


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LIGGGHTS::ContactModels;

/* ---------------------------------------------------------------------- */

FixCollisionTracker::FixCollisionTracker(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix collision/tracker command, not enough arguments");

  // int iarg = 6;

  nevery = atoi(arg[3]);

  InternalValue = 0.0;
  scalar_flag = 1;
  // global_freq = nevery;
  // extscalar = 1;
}

/* ---------------------------------------------------------------------- */

int FixCollisionTracker::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */
/*
void FixCollisionTracker::init()
{
}
*/
/* ---------------------------------------------------------------------- */

void FixCollisionTracker::end_of_step()
{

  // SurfacesIntersectData
  // sidata = ...
  // pair_gran_base.h
  // surface_model_superquadric.h

  // sidata;
  // Pair_Gran gran = static_cast<Pair_Gran>force->pair;
  //
  // pair_gran_base.h in compute_force()
  PairGran* pair_gran = static_cast<PairGran*>(force->pair_match("gran", 1));


  // No... Granular gran = new Granular(lmp, pair_gran, 4832579328);

  PairGran *pg = pair_gran;
  SurfacesIntersectData sidata;
  double ** first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;

  int ** firstneigh = pg->list->firstneigh;

  int inum = pg->list->inum;
  int * ilist = pg->list->ilist;
  const int dnum = pg->dnum();
  int * numneigh = pg->list->numneigh;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];

  
    double * const all_contact_hist = first_contact_hist ? first_contact_hist[i] : NULL;

    sidata.i = i;

    int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      sidata.j = j;


      sidata.contact_history = all_contact_hist ? &all_contact_hist[dnum*jj] : NULL;
  
      print_contact_status(sidata, pair_gran);  
    }
  }
//  int particles_were_in_contact_offset = hsetup->add_history_value("particles_were_in_contact","0");
  // enum {SURFACES_FAR, SURFACES_CLOSE, SURFACES_INTERSECT};
//  double *const particles_were_in_contact = &sidata.contact_history[particles_were_in_contact_offset];
  // if (*particles_were_in_contact == SURFACES_INTERSECT)



  // forevery particle
  // check collision with every other particle
  //
  // It saves the previous intersection point somewhere (not necessarelly a intersection tho)
  // Hopefully has a bool if there was a proper intersection



//  int particles_were_in_contact_offset = gran->add_history_value("particles_were_in_contact","0");
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

void FixCollisionTracker::print_contact_status(SurfacesIntersectData& sidata, IContactHistorySetup* hsetup)
{
  int particles_were_in_contact_offset = hsetup->add_history_value("particles_were_in_contact","0");
  double *const particles_were_in_contact = &sidata.contact_history[particles_were_in_contact_offset];

  fprintf(screen , "Particles are: %i\n", *particles_were_in_contact);
  fprintf(screen , "ENUM are: %i\n", SURFACES_INTERSECT);
  fprintf(screen , "if statment %i\n", *particles_were_in_contact == SURFACES_INTERSECT);


  /*
  if (*particles_were_in_contact == SURFACES_INTERSECT) {
    fprintf(screen , "Particles %i and %i are in contact\n", sidata.i, sidata.j);
  } else {
    fprintf(screen , "Particles %i and %i are not in contact\n", sidata.i, sidata.j);
  }
  */
} 

/* ---------------------------------------------------------------------- */

double FixCollisionTracker::compute_scalar()
{
  return InternalValue;
}
