c LAMMPS 2001 - Molecular Dynamics Simulator
c Sandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
c Steve Plimpton, sjplimp@sandia.gov
c
c Copyright (1998-2001) Sandia Corporation.  Under the terms of Contract
c DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
c certain rights in this software.  This software is distributed under 
c the GNU General Public License.
c
c See the README file in the top-level LAMMPS directory.

c -------------------------------------------------------------------------
c rRESPA partial force routines

      subroutine force_stretch(vflag)
      use global
      use mpi
      implicit none

c argument variables

      integer vflag

c local variables

      integer ierror
      double precision time1,time2

      call zero_force_virial

      time1 = mpi_wtime()

      if (nbonds.gt.0) then
        if (bondstyle == 1) then
          call bond_harmonic(vflag)
        else if (bondstyle == 2) then
          call bond_fene_standard(vflag)
        else if (bondstyle == 3) then
          call bond_fene_shift(vflag)
        else if (bondstyle == 4) then
          call bond_nonlinear(vflag)
        else
          call bond_class2(vflag)
        endif
      endif

      time2 = mpi_wtime()
      time_bond = time_bond + time2-time1

      if (newton_bond == 1) then
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
        call reverse_comm
        time2 = mpi_wtime()
        time_fcomm = time_fcomm + time2-time1
      endif

      call copy_force(nlocal,f_stretch)
      call copy_virial(nlocal,vir_stretch)

      return
      end

      
c -------------------------------------------------------------------------

      subroutine force_intra(vflag)
      use global
      use mpi
      implicit none

c argument variables

      integer vflag

c local variables

      integer ierror
      double precision time1,time2,time3,time4

      call zero_force_virial

      time1 = mpi_wtime()

      if (nangles.gt.0) then
        if (anglestyle == 1) then
          call angle_harmonic(vflag)
        else if (anglestyle == 2) then
          call angle_class2(vflag)
        else if (anglestyle == 3) then
          call angle_charmm(vflag)
        else if (anglestyle == 4) then
          call angle_cosine(vflag)
        endif
      endif
      
      time2 = mpi_wtime()
      
      if (ndihedrals.gt.0) then
        if (dihedstyle == 1) then
          call dihedral_harmonic(vflag)
        else if (dihedstyle == 2) then
          call dihedral_class2(vflag)
        else if (dihedstyle == 3) then
          call dihedral_multiharmonic(vflag)
        else if (dihedstyle == 4) then
          call dihedral_charmm(vflag)
        endif
      endif

      time3 = mpi_wtime()

      if (nimpropers.gt.0) then
        if (improstyle == 1) then
          call improper_harmonic(vflag)
        else if (improstyle == 2) then
          call improper_cvff(vflag)
        else if (improstyle == 3) then
          call improper_class2(vflag)
          call angleangle_class2(vflag)
        endif
      endif

      time4 = mpi_wtime()

      time_angle = time_angle + time2-time1
      time_dihedral = time_dihedral + time3-time2
      time_improper = time_improper + time4-time3

      if (newton_bond == 1) then
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
        call reverse_comm
        time2 = mpi_wtime()
        time_fcomm = time_fcomm + time2-time1
      endif

      call copy_force(nlocal,f_intra)
      call copy_virial(nlocal,vir_intra)

      return
      end


c -------------------------------------------------------------------------

      subroutine force_short(vflag)
      use global
      use mpi
      implicit none

c argument variables

      integer vflag

c local variables

      integer ierror
      double precision time1,time2

      call zero_force_virial

#ifdef SYNC
      call mpi_barrier(mpi_comm_world,ierror)
#endif
      time1 = mpi_wtime()

      if (nonstyle == 0) then
        if (coulstyle == 1) then
          call lj_cut_coul_cut(vflag)
        else if (coulstyle == 2) then
          call lj_cut_coul_smooth(vflag)
        else if (coulstyle == 3 .or. coulstyle == 4) then
          call lj_cut_coul_long(vflag)
        else if (coulstyle == 6) then
          call lj_cut_coul_debye(vflag)
        endif
      else if (nonstyle == 1) then
        if (coulstyle == 0) then
          call lj_cut(vflag)
        else if (coulstyle == 1) then
          call lj_cut_coul_cut(vflag)
        else if (coulstyle == 2) then
          call lj_cut_coul_smooth(vflag)
        else if (coulstyle == 3 .or. coulstyle == 4) then
          call lj_cut_coul_long(vflag)
        else if (coulstyle == 6) then
          call lj_cut_coul_debye(vflag)
        endif
      else if (nonstyle == 2) then
        if (coulstyle == 0) then
          call lj_smooth(vflag)
        else if (coulstyle == 1) then
          call lj_smooth_coul_cut(vflag)
        else if (coulstyle == 2) then
          call lj_smooth_coul_smooth(vflag)
        else if (coulstyle == 3 .or. coulstyle == 4) then
          call lj_smooth_coul_long(vflag)
        endif
      else if (nonstyle == 3) then
        call lj_shift(vflag)
      else if (nonstyle == 4) then
        call soft(vflag)
      else if (nonstyle == 5) then
        if (coulstyle == 0) then
          call ljclass2_cut(vflag)
        else if (coulstyle == 1) then
          call ljclass2_cut_coul_cut(vflag)
        else if (coulstyle == 3 .or. coulstyle == 4) then
          call ljclass2_cut_coul_long(vflag)
        endif
      else if (nonstyle == 6) then
        if (coulstyle == 3 .or. coulstyle == 4) then
          call lj_charmm_coul_long(vflag)
        else if (coulstyle == 5) then
          call lj_charmm_coul_charmm(vflag)
        endif
      endif

      time2 = mpi_wtime()
      time_nonbond = time_nonbond + time2-time1

      if (newton_nonbond == 1) then
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
        call reverse_comm
        time2 = mpi_wtime()
        time_fcomm = time_fcomm + time2-time1
      endif

      call copy_force(nlocal,f_short)
      call copy_virial(nlocal,vir_short)

      return
      end


c -------------------------------------------------------------------------

      subroutine force_long(vflag)
      use global
      use mpi
      implicit none

c argument variables

      integer vflag

c local variables

      integer ierror
      double precision time1,time2

      call zero_force_virial

      if (coulstyle == 3 .or. coulstyle == 4) then
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
        if (coulstyle == 3) then
          call ewald(vflag)
        else
          call pppm(vflag)
        endif
        time2 = mpi_wtime()
        time_long = time_long + time2-time1
      endif

      call copy_force(nlocal,f_long)

      return
      end 
