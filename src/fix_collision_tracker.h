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

#ifdef FIX_CLASS

FixStyle(collision/tracker,FixCollisionTracker)

#else

#ifndef LMP_FIX_COLLISION_TRACKER_H
#define LMP_FIX_COLLISION_TRACKER_H

#include "fix.h"
#include "contact_interface.h"
#include "pair_gran.h"
#include "tri_mesh.h"
#include "fix_contact_history_mesh.h"
#include <vector>

using namespace LIGGGHTS;
using namespace ContactModels;

namespace LAMMPS_NS {

class FixCollisionTracker : public Fix {
 public:
  FixCollisionTracker(class LAMMPS *, int, char **);
  ~FixCollisionTracker();
  int setmask();
  void post_force(int);
  void end_of_step();

  void compute_normal(SurfacesIntersectData &);
  void compute_relative_velocity(SurfacesIntersectData &, double*, double*);
 private:

  enum {SURFACES_FAR, SURFACES_CLOSE, SURFACES_INTERSECT};
  int pre_particles_were_in_contact_offset;
  int particles_were_in_contact_offset;
  int contact_point_offset;
  class PairGran *pair_gran;

  std::vector<double> rel_vels;
  std::vector<double> lcol;
  double* rel_all; 
  double* col_all;
  int array_offset = 0;

  bool cube_projection;
  int x_nsplit;
  int y_nsplit;
  int z_nsplit;
  int** x_octsurface;
  int** y_octsurface;
  int** z_octsurface;
  int** x_octsurface_all;
  int** y_octsurface_all;
  int** z_octsurface_all;

  char* filename;  
  int me;
  FILE *fp;
  int writetofile = 0;

  void openfile();
  void resolve_contact_status(SurfacesIntersectData &); 
  bool check_collision(SurfacesIntersectData&);
  void print_atom_pair_info(int i, int j);
  void print_atom_info(int i);

  void compute_local_contact(SurfacesIntersectData& sidata, double *iResult, double *jResult);

  double* get_triangle_contact_history(TriMesh *mesh, FixContactHistoryMesh *fix_contact, int iPart, int iTri);

  void unit_cube_oct_projection(int iPart, double *contact, double *result);
  void unit_cube_oct_indexing(double *cube_projection);
  void print_cube_projection(FILE *fp);
};

}

#endif
#endif
