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
#include <string>
#include <unordered_map>
#include <vector>
#include "atom.h"
#include "pair.h"
#include "update.h"
#include "error.h"
#include "fix_collision_tracker.h"
#include "pair_gran.h"
#include "force.h"
#include "group.h"

#include "neigh_list.h"

#include "contact_interface.h"

#include <unordered_map>

#include "fix.h"
#include "fix_wall_gran.h"
#include "fix_contact_history_mesh.h"
#include "tri_mesh.h"
#include "primitive_wall.h"

#include "superquadric.h"
#include "math_extra_liggghts_nonspherical.h"

#include <limits>

#if !defined(_WINDOWS) && !defined(__MINGW32__)
#include <sys/stat.h>
#endif

using namespace LIGGGHTS;
using namespace ContactModels;
using namespace FixConst;

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
  writeraw = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"raw") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix collision tracker: raw command");
      if (me == 0) {
        int n = strlen(arg[iarg+1]) + 1;
        rawname = new char[n];
        strcpy(rawname,arg[iarg+1]);
      }
      SetGroupShapeBlockiness();
      writeraw = 1;
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"cpoctant") == 0)
    {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix cpoctant command");

      cube_projection = 1;

      if (iarg+3 > narg || atoi(arg[iarg+2]) == 0) // single number
      {
        int N = atoi(arg[iarg+1]);
        if (N < 1) error->all(FLERR, "number of buckets N is less than 1");
      
        SetGroupShapeBlockiness();
        //double* shape = atom->shape[0];
        double factor = shape[1] + shape[2] + shape[1]*shape[2]/shape[0];
        x_nsplit = (int) sqrt(N * shape[0] / factor);
        x_nsplit = x_nsplit < 1 ? 1 : x_nsplit;
        y_nsplit = (int) (x_nsplit * shape[1] / shape[0]);
        y_nsplit = y_nsplit < 1 ? 1 : y_nsplit;
        z_nsplit = (int) (x_nsplit * shape[2] / shape[0]);
        z_nsplit = z_nsplit < 1 ? 1 : z_nsplit;
        iarg += 2;
      }
      else if (iarg + 3 < narg)
      {
        // x,y,z sides
        x_nsplit = atoi(arg[iarg+1]);
        if (x_nsplit < 1) error->all(FLERR,"x axis split < 1");

        y_nsplit = atoi(arg[iarg+2]);
        if (y_nsplit < 1) error->all(FLERR,"y axis split < 1");
        
        z_nsplit = atoi(arg[iarg+3]);
        if (z_nsplit < 1) error->all(FLERR,"z axis split < 1");

        iarg += 4;
      }
      else
        error->all(FLERR,"Illegal fix cpoctant command");

      // check for more arguments
      while (iarg < narg && strcmp(arg[iarg], "range") == 0)
      {
        if (iarg+2<narg && strcmp(arg[iarg+1], "full")==0)
        {
          int n = strlen(arg[iarg+2]) + 1;
          octfilenames.push_back(new char[n]);
          strcpy((octfilenames.back()),arg[iarg+2]);
          rangefrom.push_back(-1.);
          rangeto.push_back(std::numeric_limits<double>::max());
          iarg += 3;
        }
        else if (iarg+3<narg)
        {
          int n = strlen(arg[iarg+3]) + 1;
          octfilenames.push_back(new char[n]);
          strcpy((octfilenames.back()),arg[iarg+3]);
          rangefrom.push_back(std::stod(arg[iarg+1]));
          rangeto.push_back(std::stod(arg[iarg+2]));
          iarg += 4;
        }
        else
          iarg++;
      }
      if (octfilenames.size() == 0)
      {
        // create a default file
        int n = strlen("cpoctant_full*.csv") + 1;
        octfilenames.push_back(new char[n]);
        strcpy((octfilenames.back()),"cpoctant_full*.csv");
        rangefrom.push_back(-1.);
        rangeto.push_back(std::numeric_limits<double>::max());
      }
      nfiles = octfilenames.size();

      memory->create(x_octsurface, nfiles, y_nsplit, z_nsplit, "xprojection");
      memory->create(y_octsurface, nfiles, x_nsplit, z_nsplit, "yprojection");
      memory->create(z_octsurface, nfiles, x_nsplit, y_nsplit, "zprojection");

      // Zero out newly allocated array
      // Using a trick based on how memory->creates allocates memory
      memset(x_octsurface[0][0], 0, sizeof(***x_octsurface)*y_nsplit*z_nsplit*nfiles);
      memset(y_octsurface[0][0], 0, sizeof(***y_octsurface)*x_nsplit*z_nsplit*nfiles);
      memset(z_octsurface[0][0], 0, sizeof(***z_octsurface)*x_nsplit*y_nsplit*nfiles);

      if (me == 0)
      {
        memory->create(x_octsurface_all, y_nsplit, z_nsplit, "xprojection_all");
        memory->create(y_octsurface_all, x_nsplit, z_nsplit, "yprojection_all");
        memory->create(z_octsurface_all, x_nsplit, y_nsplit, "zprojection_all");
      }
      else
      {
        x_octsurface_all = y_octsurface_all = z_octsurface_all = x_octsurface[0];  // so that they can acess a first index for the MPI buiss
      }
    }
    else if(strcmp(arg[iarg],"type") == 0)
    {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix collision tracker: type command");
      
      if(strcmp(arg[iarg+1],"wall") == 0) {
        store_particle = 0; 
        store_wall = 1; 
      }
      if(strcmp(arg[iarg+1],"particle") == 0) {
        store_particle = 1; 
        store_wall = 0; 
      }
      if(strcmp(arg[iarg+1],"permesh") == 0) {
        permesh = 1;
      }
      else
        error->all(FLERR,"Illegal fix collision tracker: type command");
      
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "othergroup") == 0)
    {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix collision tracker: othergroup command");

      int igroup = group->find(arg[iarg+1]);
      if (igroup == -1) error->all(FLERR,"Could not find fixcollisiontrackers othergroup fix group ID");
      othergroupbit = group->bitmask[igroup];
      useothergroup = 1;

      iarg += 2;
    }
    else
      iarg+=1;
  }

  if (me == 0)
  {
    rel_all = memory->create(rel_all, 2, "Fix_Collision_Tracker_velocities");
    col_all = memory->create(col_all, 3, "Fix_Collision_Tracker_collision");
  }


  // check whether the folder is accessible, not available on windows
#if !defined(_WINDOWS) && !defined(__MINGW32__)
  if (me == 0)
  {
    if (writeraw)
      create_folder(rawname);
    if (cube_projection)
      for (int i = 0; i < nfiles; ++i)
        create_folder(octfilenames[i]);
  }
#endif

  n_meshes = 0;
  int n_wall_fixes = modify->n_fixes_style("wall/gran");
  for (int ifix = 0; ifix < n_wall_fixes; ++ifix)
  {
      FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
      if (fwg->is_mesh_wall())
      {
        int n_FixMesh = fwg->n_meshes();
        n_meshes += n_FixMesh;
      }
  }

  // perform initial allocation of atom-based array
  // register with Atom class
  mesh_contact = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  int nlocal = atom->nlocal;

  // setting no last mesh contact
  for (int i = 0; i < nlocal; i++)
  {
    for(int k = 0; k < n_meshes; k++)
    {
      mesh_contact[i][k] = -1;
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::create_folder(std::string fname)
{
#if !defined(_WINDOWS) && !defined(__MINGW32__)
  std::size_t last_slash = fname.rfind("/");
  // check if we use directories at all
  if (last_slash != std::string::npos)
  {
      std::size_t next_slash = fname.find("/", 1);
      while (next_slash != std::string::npos)
      {
          std::string curdir = fname.substr(0, next_slash);
          struct stat statbuf;
          const bool exists = (stat(curdir.c_str(), &statbuf) != -1) && S_ISDIR(statbuf.st_mode);
          if (!exists)
              mkdir(curdir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP);
          next_slash = fname.find("/", next_slash+1);
      }
  }
#endif
}

/* ---------------------------------------------------------------------- */

FixCollisionTracker::~FixCollisionTracker()
{
  memory->sfree(array_local);
  memory->destroy(x_octsurface);
  memory->destroy(y_octsurface);
  memory->destroy(z_octsurface);
  if (me == 0)
  {
    memory->destroy(rel_all);
    memory->destroy(col_all);
    if (cube_projection)
    {
      memory->destroy(x_octsurface_all);
      memory->destroy(y_octsurface_all);
      memory->destroy(z_octsurface_all);
    }
  }

  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);

  // delete locally stored array
  memory->destroy(mesh_contact);

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
  //if (cube_projection || writeraw)
  //  SetGroupShapeBlockiness();  // Does it need to be run each time we print?

  if(cube_projection)
  {
    for (int i = 0; i < nfiles; ++i)
    {
      openfile(octfilenames[i]);
      print_cube_projection(fp, i);
 
      if (me == 0)
      {
        fflush(fp);
        fclose(fp);
      }
    }
    // Zero out array
    // Using a trick based on how memory->creates allocates memory
    memset(x_octsurface[0][0], 0, sizeof(***x_octsurface)*y_nsplit*z_nsplit*nfiles);
    memset(y_octsurface[0][0], 0, sizeof(***y_octsurface)*x_nsplit*z_nsplit*nfiles);
    memset(z_octsurface[0][0], 0, sizeof(***z_octsurface)*x_nsplit*y_nsplit*nfiles);
  }
  if (writeraw)
  {
    openfile(rawname);

    double * rel = rel_vels.data();
    double * col = lcol.data();
    int n = rel_vels.size() / 2;

    int n_all, size;
    MPI_Reduce(&n, &n_all, 1, MPI_INT, MPI_SUM, 0, world);
    MPI_Comm_size(world, &size);
    if (me == 0)
    {
      rel_all = memory->grow(rel_all, n_all*2, "Fix_Collision_Tracker_velocities");
      col_all = memory->grow(col_all, n_all*3, "Fix_Collision_Tracker_collision");
    }
    int ns[size];
    MPI_Gather(&n, 1, MPI_INT, ns, 1, MPI_INT, 0, world);
    int vns[size], cns[size], dvns[size], dcns[size];
    vns[0] = 2*ns[0];
    cns[0] = 3*ns[0];
    dvns[0] = dcns[0] = 0;
    for (int i = 1; i < size; ++i)
    {
      vns[i] = 2*ns[i];
      cns[i] = 3*ns[i];
      dvns[i] = dvns[i-1]+2*ns[i-1];
      dcns[i] = dcns[i-1]+3*ns[i-1];
    }
    MPI_Gatherv(rel, n*2, MPI_DOUBLE, rel_all, vns, dvns, MPI_DOUBLE, 0, world);
    MPI_Gatherv(col, n*3, MPI_DOUBLE, col_all, cns, dcns, MPI_DOUBLE, 0, world);

    if (me == 0){
      fprintf(fp, "# %f %f %f %f %f\n", shape[0], shape[1], shape[2], blockiness[0], blockiness[1]);

      for (int i = 0; i < n_all; ++i)
        fprintf(fp,"%f %f %f %f %f\n", rel_all[2*i], rel_all[2*i+1], col_all[3*i], col_all[3*i+1], col_all[3*i+2]);
      
      fflush(fp);
      fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::openfile(char* filename)
{
  if (me != 0)
    return;

  // if one file per timestep, replace '*' with current timestep  
  char *filecurrent = filename;
  
  char *filestar = filecurrent;
  filecurrent = new char[strlen(filestar) + 16 + 4];
  char *ptr = strchr(filestar,'*');
  if (ptr) {
    *ptr = '\0';
    sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
            filestar,update->ntimestep, ptr+1);
    *ptr = '*';
      
    fp = fopen(filecurrent,"w");
  } else {
    sprintf(filecurrent, BIGINT_FORMAT "_%s", update->ntimestep, filestar);
    fp = fopen(filecurrent, "w");
  }

  if (fp == NULL) error->one(FLERR,"Cannot open collision tracking file.");

  delete [] filecurrent;
}

void FixCollisionTracker::print_cube_projection(FILE *fp, int i)
{
  MPI_Reduce(x_octsurface[i][0], x_octsurface_all[0], y_nsplit*z_nsplit, MPI_INT, MPI_SUM, 0, world); 
  MPI_Reduce(y_octsurface[i][0], y_octsurface_all[0], x_nsplit*z_nsplit, MPI_INT, MPI_SUM, 0, world); 
  MPI_Reduce(z_octsurface[i][0], z_octsurface_all[0], y_nsplit*x_nsplit, MPI_INT, MPI_SUM, 0, world); 

  if (me != 0)
    return;

  fprintf(fp,"%f, %f, %f, %f, %f\n", shape[0], shape[1], shape[2], blockiness[0], blockiness[1]);
  fprintf(fp,"%d, %d, %d\n", x_nsplit, y_nsplit, z_nsplit);
  // Print x side
  for(int y = 0; y < y_nsplit; y++)
  {
    fprintf(fp,"%d", x_octsurface_all[y][0]);
    for(int z = 1; z < z_nsplit; z++)
    {
      fprintf(fp,",%d", x_octsurface_all[y][z]);
    }
    fprintf(fp,"\n");
  }

  // Print y side
  for(int x = 0; x < x_nsplit; x++)
  {
    fprintf(fp,"%d", y_octsurface_all[x][0]);
    for(int z = 1; z < z_nsplit; z++)
    {
      fprintf(fp,",%d", y_octsurface_all[x][z]);
    }
    fprintf(fp,"\n");
  }

  // Print z side
  for(int x = 0; x < x_nsplit; x++)
  {
    fprintf(fp,"%d", z_octsurface_all[x][0]);
    for(int y = 1; y < y_nsplit; y++)
    {
      fprintf(fp,",%d", z_octsurface_all[x][y]);
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::post_force(int vflag)
{
  // clear the buffers every time if it only needs to be avavibel for other fixes, or after every time we write out
  if (!writeraw || (update->ntimestep - 1) % nevery == 0)
  {
    array_offset = 0;
    rel_vels.clear();
    lcol.clear();
  }
  // Get local atom information
  int *mask = atom->mask;

  if (store_particle)
  {
    PairGran *pg = pair_gran;
    double ** first_contact_hist = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;
    int ** firstneigh = pg->list->firstneigh;

    int inum = pg->list->inum;
    int * ilist = pg->list->ilist;
    const int dnum = pg->dnum();
    int * numneigh = pg->list->numneigh;
    
    SurfacesIntersectData sidata;

    // particle - particle
    for (int ii = 0; ii < inum; ii++)
    {
      const int i = ilist[ii];
      sidata.i = i;
    
      double * const all_contact_hist = first_contact_hist ? first_contact_hist[i] : NULL;

      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        if (!(mask[j] & groupbit) && !(mask[i] & groupbit))
          continue;

        sidata.j = j;

        sidata.contact_history = all_contact_hist ? &all_contact_hist[dnum*jj] : NULL;
      
        if (check_collision(sidata))
          resolve_contact_status(sidata);
      }
    }
  }

  if (store_wall)
  {
    // particle - wall (only mesh and mesh walls, primitive walls are not included) 
    // Based on public fork https://github.com/ParticulateFlow/LIGGGHTS-PFM (2021-12-26)
    int mesh_number = 0;
    int nlocal = atom->nlocal;
    int n_wall_fixes = modify->n_fixes_style("wall/gran");
    //printf("wall_coll n_wall_fixes:%d\n",n_wall_fixes);
    for (int ifix = 0; ifix < n_wall_fixes; ++ifix)
    {
      FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));

      if (fwg->is_mesh_wall())
      {
        int n_FixMesh = fwg->n_meshes();
        //printf("wall_coll n_FixMesh:%d\n",n_FixMesh);
        for (int iMesh = 0; iMesh < n_FixMesh; iMesh++)
        {
          TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
          FixContactHistoryMesh *fix_contact = fwg->mesh_list()[iMesh]->contactHistory();
          if (!fix_contact) continue; // Skip if null

          // get neighborList and numNeigh
          FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();
          if (!meshNeighlist) continue; // Skip if null

          // Do things for moving wall here
          double *** vMesh = NULL;
          MultiVectorContainer<double,3,3> *vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
          if(vMeshC)
            vMesh = vMeshC->begin();

          for (int iPart = 0; iPart < nlocal; ++iPart)
          {
            if (!(mask[iPart] & groupbit)) continue;
            
            const int nneighs = fix_contact->nneighs(iPart);
            // Disregard any collisions if there has already been contact with this mesh at any triangle
            /*
            bool prev_contact = 0;
            if (permesh)
            {
              for (int j = 0; j < nneighs; ++j)
              {
                const int idTri = fix_contact->partner(iPart, j);
                if (idTri != -1)
                {
                  double * contact_history = fix_contact->contacthistory(iPart, j);
                  printf("idTri: %d \n", idTri);
                  if (contact_history && contact_history[pre_particles_were_in_contact_offset] == 1)
                  {
                    prev_contact = 1;
                    break;
                  }
                }
              }
            }
            */
            bool prev_contact = 0;
            bool current_contact = 0;
            //printf("mesh_contact[%d][%d] = %d; timestep: %d\n",iPart,mesh_number,mesh_contact[iPart][mesh_number],update->ntimestep);
            if (mesh_contact[iPart][mesh_number] == 1)
            {
              // Previous contact was had
              prev_contact = 1;
            }

            for (int j = 0; j < nneighs; ++j)
            {
              const int idTri = fix_contact->partner(iPart, j);
              //printf("idTri: %d\n",idTri);
              if (idTri != -1)
              {
                double * contact_history = fix_contact->contacthistory(iPart, j);

                if (contact_history && contact_history[pre_particles_were_in_contact_offset] == 0)
                {
                  int iTri = mesh->map(idTri, 0); // can it be something else than 0? the implication is that one idTri can have more than one iTri
                  
                  Superquadric particle(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
                  double delta[3], contact_point[3], bary[3];
                  // Recalculating contact_point. It would be preferable to access already calculated values
                  mesh->resolveTriSuperquadricContact(iTri, delta, contact_point, particle, bary);

                  // Negative Baryocentric coordinates are outside of the triangle and not actual contact
                  if(bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0 )
                  {
                    if (!prev_contact)
                    {
                      resolve_mesh_contact_status(vMesh, iPart, iTri, bary, contact_point);
                    }
                    current_contact = 1;

                    //printf("iPart: %d, iTri: %d; (%f,%f,%f); (%f,%f,%f)\n", iPart, iTri, contact_point[0],contact_point[1],contact_point[2],bary[0],bary[1],bary[2]);

                    // contact_history is not used by wall - particle calculations, but it is zeroed out between contacts
                    // We take advantage of that by setting pre_particles_were_in_contact_offset to 1
                    contact_history[pre_particles_were_in_contact_offset] = 1;
                    // Check against multiple collisions happening when particle slides across mesh wall is missing
                    // Baryocentric coordinates may be used for that check
                  }
                }
                else if (contact_history && contact_history[pre_particles_were_in_contact_offset] == 1)
                {
                  current_contact = 1;
                }
              }
            }

            if (current_contact)
            {
              mesh_contact[iPart][mesh_number] = 1;
            }
            else
            {
              mesh_contact[iPart][mesh_number] = 0;
            }
          }
          // changing to a different mesh
          mesh_number++;
        }
      }
      else
      {
        if (fwg->dnum() > 0)
        {
          char *hist_name = new char[strlen(fwg->id)+1+10];
          strcpy(hist_name,"history_");
          strcat(hist_name,fwg->id);
          FixPropertyAtom *fix_history_primitive =
                  static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name,"property/atom","vector",fwg->dnum(),0,fwg->style));
          delete []hist_name;
          double ** wall_history = fix_history_primitive->array_atom;

          int* neighbourList;
          int nNeigh = fwg->primitiveWall()->getNeighbors(neighbourList);

          for (int iCont = 0; iCont < nNeigh; ++iCont, ++neighbourList)
          {
            int iPart = *neighbourList;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;


            double * contact_history = wall_history[iPart];
            if (contact_history && contact_history[pre_particles_were_in_contact_offset] == 0)
            {
              double delta[3];
              double deltan = fwg->primitiveWall()->resolveContact(atom->x[iPart], vectorMax3D(atom->shape[iPart]), delta);

              if (deltan <= 0)
              {
                double sphere_contact_point[3];
                vectorAdd3D(atom->x[iPart], delta, sphere_contact_point);
                double closestPoint[3], point_of_lowest_potential[3];
                Superquadric particle(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
                bool intersectflag = particle.plane_intersection(delta, sphere_contact_point, closestPoint, point_of_lowest_potential);
                
                if (intersectflag)
                {
                  resolve_primitive_contact_status(iPart, closestPoint);
            
                  // Less than sure this isn't used...  
                  contact_history[pre_particles_were_in_contact_offset] = 1;
                  // But it is set to zero when the contact is lifted
                }
              }
            }
          }
        }
      }
    }
  }
  // This is as weird as it looks, but LIGGGHTS want their array this way
  array_local[0] = rel_vels.data() + array_offset;
  size_local_rows = (rel_vels.size() - array_offset) / 2;
  array_offset += size_local_rows * 2;
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::SetGroupShapeBlockiness()
{
  // All particles in the same group are assumed to have the same shape
  // This is for printing reasons, does not affect # of collisions
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int i;
  for (i = 0; i < nlocal; ++i)
    if (mask[i] & groupbit)
      break;
  
  bool found = 1;
  if (i == nlocal) // no instance of group on this processor
  {
    i = 0;
    found = 0;
  }

  double shp[3] = {-10.0,-10.0,-10.0};
  double blo[2] = {0.0,0.0};
  if (found)
  {
    vectorCopyN(atom->shape[i], shp, 3);
    vectorCopyN(atom->blockiness[i], blo, 2);
  } 

  int size = 0;
  MPI_Comm_size(world, &size);

  double shapes[size*3];
  double blocks[size*2];
  
  MPI_Gather(shp, 3, MPI_DOUBLE, shapes, 3, MPI_DOUBLE, 0, world);
  MPI_Gather(blo, 2, MPI_DOUBLE, blocks, 2, MPI_DOUBLE, 0, world);
  if (me == 0)
  {
    int j;
    for (j = 0; j < size; ++j)
    {
      if (shapes[j*3] > -1)
      {
        vectorCopyN(&shapes[j*3], shape, 3);
        vectorCopyN(&blocks[j*2], blockiness, 2);
        break;
      }
    }
    if (j == size) // set some values to avoid issues if non are found
    {
      vectorCopyN(shapes, shape, 3);
      vectorCopyN(blocks, blockiness, 2);
    }
  }
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

void FixCollisionTracker::store_data(int iPart, double velnormal, double veltangent, double* contact_point)
{
  rel_vels.push_back(velnormal);
  rel_vels.push_back(veltangent);
  lcol.push_back(contact_point[0]);
  lcol.push_back(contact_point[1]);
  lcol.push_back(contact_point[2]);

  if (cube_projection)
  {           
    double result[3];
    unit_cube_oct_projection(iPart, contact_point, result);
    for (int i = 0; i < nfiles; ++i)
    {
      if (rangefrom[i] <= velnormal && rangeto[i] >= velnormal)
        unit_cube_oct_indexing(result, i);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::resolve_primitive_contact_status(int iPart, double* contact_point)
{
  double normal[3];
  compute_normal_wall(iPart, contact_point, normal);

  // v1 = v + cross(w,(p-x))  speed particle
  double a_i[3];
  double b_i[3];
  vectorSubtract3D(contact_point, atom->x[iPart], a_i);
  vectorCross3D(atom->omega[iPart], a_i, b_i);
  vectorAdd3D(atom->v[iPart], b_i, a_i);


  double res = vectorDot3D(a_i, normal);

  double velnormal = res > 0 ? res : -res;
  double veltangent = sqrt(vectorDot3D(a_i, a_i) - res * res);

  double iLocal[3], iResult[3];
  vectorSubtract3D(contact_point, atom->x[iPart], iLocal);
  MathExtraLiggghtsNonspherical::rotate_global2local(atom->quaternion[iPart], iLocal, iResult);
  store_data(iPart, velnormal, veltangent, iResult);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::resolve_contact_status(SurfacesIntersectData& sidata)
{
  bool jislocal = sidata.j < atom->nlocal;

  compute_normal(sidata);
  double vels[2]; 
  compute_relative_velocity(sidata, vels, vels+1);

  double point_of_contact[6];
  compute_local_contact(sidata, point_of_contact, point_of_contact+3);

  // Get local atom information
  int *mask = atom->mask;

  if ((mask[sidata.i] & groupbit) && (!useothergroup || mask[sidata.j] & othergroupbit))
    store_data(sidata.i, vels[0], vels[1], point_of_contact);
  if (jislocal && (mask[sidata.j] & groupbit) && (!useothergroup || mask[sidata.i] & othergroupbit))
    store_data(sidata.j, vels[0], vels[1], point_of_contact+3);
} 

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::resolve_mesh_contact_status(double ***vMesh, int iPart, int iTri, double* bary, double* contact_point)
{
  double normal[3];
  compute_normal_wall(iPart, contact_point, normal);

  // v1 = v + cross(w,(p-x))  speed particle
  double a_i[3];
  double b_i[3];
  vectorSubtract3D(contact_point, atom->x[iPart], a_i);
  vectorCross3D(atom->omega[iPart], a_i, b_i);
  vectorAdd3D(atom->v[iPart], b_i, a_i);
  // speed wall
  double v_wall[3] = {0.,0.,0.};
  if(vMesh)
  {
      for(int i = 0; i < 3; i++)
          v_wall[i] = (bary[0]*vMesh[iTri][0][i] +
                      bary[1]*vMesh[iTri][1][i] +
                      bary[2]*vMesh[iTri][2][i] );
  }

  double rel_v[3];
  vectorSubtract3D(a_i, v_wall, rel_v);
  double res = vectorDot3D(rel_v, normal);

  double velnormal = res > 0 ? res : -res;
  double veltangent = sqrt(vectorDot3D(rel_v, rel_v) - res * res);

  double iLocal[3], iResult[3];
  vectorSubtract3D(contact_point, atom->x[iPart], iLocal);
  MathExtraLiggghtsNonspherical::rotate_global2local(atom->quaternion[iPart], iLocal, iResult);
  store_data(iPart, velnormal, veltangent, iResult);
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

void FixCollisionTracker::compute_normal_wall(const int iPart, const double* contact_point, double* normal)
{
  Superquadric particle_i(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);

  particle_i.shape_function_gradient_global(contact_point, particle_i.gradient);

  vectorCopy3D(particle_i.gradient, normal);
  vectorNormalize3D(normal);
}

/* ---------------------------------------------------------------------- */

void FixCollisionTracker::compute_normal(SurfacesIntersectData& sidata)
{
  double *const prev_step_point = &sidata.contact_history[contact_point_offset];
  
  int iPart = sidata.i;
  int jPart = sidata.j;

  Superquadric particle_i(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
  Superquadric particle_j(atom->x[jPart], atom->quaternion[jPart], atom->shape[jPart], atom->blockiness[jPart]);

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

void FixCollisionTracker::unit_cube_oct_indexing(double *cube_proj, int i)
{
  if (cube_proj[0] >= cube_proj[1] && cube_proj[0] >= cube_proj[2]) // Project to x side
  {
    //calculating indexes
    int y = std::min(y_nsplit-1, (int)floor(y_nsplit*cube_proj[1]));
    int z = std::min(z_nsplit-1, (int)floor(z_nsplit*cube_proj[2]));
    x_octsurface[i][y][z]++; //(y,z)
  }
  else if(cube_proj[1] >= cube_proj[2]) // Project to y side
  {
    //calculating indexes
    int x = std::min(x_nsplit-1, (int)floor(x_nsplit*cube_proj[0]));
    int z = std::min(z_nsplit-1, (int)floor(z_nsplit*cube_proj[2]));
    y_octsurface[i][x][z]++; //(x,z)
  }
  else // Project to z side
  {
    //calculating indexes
    int x = std::min(x_nsplit-1, (int)floor(x_nsplit*cube_proj[0]));
    int y = std::min(y_nsplit-1, (int)floor(y_nsplit*cube_proj[1]));
    z_octsurface[i][x][y]++; //(x,y)
  } 
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


/* ----------------------------------------------------------------------
   Handling of local atom-based array. Based on other fixes
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixCollisionTracker::memory_usage()
{
  double bytes = atom->nmax*n_meshes * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixCollisionTracker::grow_arrays(int nmax)
{
  memory->grow(mesh_contact,nmax,n_meshes,"collision/tracker:mesh_contact");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixCollisionTracker::copy_arrays(int i, int j, int delflag)
{
  for(int k = 0; k < n_meshes; k++)
  {
    mesh_contact[j][k] = mesh_contact[i][k];
  }
  //mesh_contact[j][0] = mesh_contact[i][0];
  //mesh_contact[j][1] = mesh_contact[i][1];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixCollisionTracker::pack_exchange(int i, double *buf)
{
  for(int k = 0; k < n_meshes; k++)
  {
    buf[k] = mesh_contact[i][k];
  }
  return n_meshes;

  //buf[0] = mesh_contact[i][0];
  //buf[1] = mesh_contact[i][1];
  //return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixCollisionTracker::unpack_exchange(int nlocal, double *buf)
{
  for(int k = 0; k < n_meshes; k++)
  {
    mesh_contact[nlocal][k] = buf[k];
  }
  return n_meshes;

  //mesh_contact[nlocal][0] = buf[0];
  //mesh_contact[nlocal][1] = buf[1];
  //return 2;
}