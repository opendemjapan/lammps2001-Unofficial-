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
c timestep integrator

      subroutine integrate
      use global
      use mpi
      implicit none

c local variables

      integer i,vflag,ierror,nflag
      real*8 p_omega_dt(3),time1,time2,time3,time4,time5
      real*8 target

      call mpi_barrier(mpi_comm_world,ierror)
      time_loop = mpi_wtime()

      do itime = 1,nsteps

        ntimestep = ntimestep + 1

c should virial be computed this timestep
c yes, if NPT or NPH or if printing thermodynamics
c otherwise no
c vflag = 2 means compute virial via ghost forces in store_virial
c only do this if newton is fully on so ghost forces are accurate

        vflag = 0
        if (ensemble >= 3) then
          vflag = 1
        else if (nthermo_next == ntimestep) then
          vflag = 1
        endif
        if (vflag == 1 .and. newton == 3) vflag = 2

c 1st half of velocity Verlet update

        if (ensemble == 1) then
          call nve_v(f,dthalf)
        else if (ensemble == 2) then
          call nvt_eta(dthalf,1)
          call nvt_v(f,dthalf,1)
        else if (ensemble == 3) then
          call nph_eta
          call npt_omegadot(dthalf,p_omega_dt)
          call npt_v(f,dthalf,1)
          call npt_x(dthalf,p_omega_dt,1)
        else if (ensemble == 4) then
          call nvt_eta(dthalf,1)
          call npt_omegadot(dthalf,p_omega_dt)
          call npt_v(f,dthalf,1)
          call npt_x(dthalf,p_omega_dt,1)
        endif

        call nve_x

c nflag = 0/1 = no/yes neighbor list rebuild on this step
c force rebuild on steps restart files are written
c   this (1) insures atoms in file will have PBC applied and
c        (2) enables exact restarts since atoms will migrate to new procs
c            before being written to file (same as when file is read)
c else check for neighbor list rebuild if neighdelay steps have occurred

        nflag = 0
        neighago = neighago + 1

        if (ntimestep == nrestart_next) then
          nflag = 1
        else if (neighago >= neighdelay .and.
     $         mod(neighago,neighfreq) == 0) then
          if (neightrigger == 0) then
            nflag = 1
          else
            call check_neighbor(nflag)
          endif
        endif

c ghost coord communication if no neighbor list rebuild or
c neighbor list rebuild w/ communication:
c   (1) remap atoms into periodic box
c   (2) shrink-wrap box around atoms if any dims are nonperiodic
c   (3) reset comm and neighbor-bin patterns if box size has changed
c   (4) move atoms to new procs via exchange routine
c   (5) acquire ghost atoms via borders routine
c   (6) build all new neighbor lists

#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
          
        if (nflag == 0) then
          call communicate
          time2 = mpi_wtime()
          time_comm = time_comm + time2-time1
        else
          call pbc
          if (nonperiodic) call setup_box
          if (nonperiodic .or. pressstyle > 0 .or. volstyle > 0) then
            call setup_comm
            if (neighstyle == 1) call setup_neigh
          endif
          call exchange
          time2 = mpi_wtime()
          call borders
          time3 = mpi_wtime()
#ifdef BONDSWAP
          if (numdiag > 0) call diagnostic
#endif
          call neighbor
          time4 = mpi_wtime()
          if (nmolecular > 0) then
            if (nbonds > 0) call neighbor_bond
            if (nangles > 0) call neighbor_angle
            if (ndihedrals > 0) call neighbor_dihedral
            if (nimpropers > 0) call neighbor_improper
            if (nshake == 1) call neighbor_shake
            time5 = mpi_wtime()
            time_neigh2 = time_neigh2 + time5-time4
          endif

          time_exch = time_exch + time2-time1
          time_comm = time_comm + time3-time2
          time_neigh1 = time_neigh1 + time4-time3
        endif

c zero forces and virial

        call zero_force_virial

c bonded forces

        if (nmolecular == 1) then

#ifdef SYNC
          call mpi_barrier(mpi_comm_world,ierror)
#endif
          time1 = mpi_wtime()

          if (nbonds > 0) then
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

          time3 = mpi_wtime()

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

          time4 = mpi_wtime()

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

          time5 = mpi_wtime()

          time_bond = time_bond + time2-time1
          time_angle = time_angle + time3-time2
          time_dihedral = time_dihedral + time4-time3
          time_improper = time_improper + time5-time4

        endif

c make copy of bond forces/virial so can compute nonbond virial via ghost atoms
c zero force, since must accumulate f afresh in nonbond/long-range routines

        if (vflag == 2) then
          call copy_force(nlocal+nghost,fhold)
          call copy_virial(nlocal+nghost,virialhold)
          call zero_force_virial
        endif

c nonbond forces

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

c compute nonbond virial via ghost forces
c also recombines previous copies of forces

        if (vflag == 2) call store_virial

c long-range forces

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

c reverse force communication

        if (newton >= 1) then
#ifdef SYNC
          call mpi_barrier(mpi_comm_world,ierror)
#endif
          time1 = mpi_wtime()
          call reverse_comm
          time2 = mpi_wtime()
          time_fcomm = time_fcomm + time2-time1
        endif

c global temperature control via Langevin drag

        if (tempstyle == 3) call langevin(f,dthalf)

c constraints, including SHAKE

        if (nfixes > 0) call fix_apply(1,f,dthalf)
        if (nshake == 1) then
          if (ntimestep == nshake_next) call shake_bond_stats
#ifdef SYNC
          call mpi_barrier(mpi_comm_world,ierror)
#endif
          time1 = mpi_wtime()
          call shake(1,vflag)
          time2 = mpi_wtime()
          time_shake = time_shake + time2-time1
        endif

c 2nd half of velocity Verlet update

        if (ensemble == 1) then
          call nve_v(f,dthalf)
        else if (ensemble == 2) then
          call nvt_v(f,dthalf,2)
          call temperature
          call nvt_eta(dthalf,2)
        else if (ensemble == 3) then
          call npt_v(f,dthalf,2)
          call npt_x(dthalf,p_omega_dt,2)
          call temperature
          call nph_eta
          call pressure
          call npt_omegadot(dthalf,p_omega_dt)
        else if (ensemble == 4) then
          call npt_v(f,dthalf,2)
          call npt_x(dthalf,p_omega_dt,2)
          call temperature
          call nvt_eta(dthalf,2)
          call pressure
          call npt_omegadot(dthalf,p_omega_dt)
        endif

c global temperature control via velocity rescale and replace

        if (tempstyle == 1) then
          if (mod(itime,t_every) == 0) then
            call temperature
            target = t_start + float(itime)/nsteps*(t_stop-t_start)
            if (abs(t_current-target) > t_window) then
              target = t_current - t_fraction*(t_current-target)
              call temp_rescale(target)
            endif
          endif
        else if (tempstyle == 2) then
          if (mod(itime,t_every) == 0) then
            target = t_start + float(itime)/nsteps*(t_stop-t_start)
            call temp_replace(target)
          endif
        endif

c global volume expand/contract

        if (volstyle == 1) call vol_rescale

c outputs

        if (ntimestep == noutput_next) then
          if (ntimestep == nthermo_next) call thermo(0)
          if (ntimestep == ndumpatom_next) call dump_atom(0)
          if (ntimestep == ndumpvel_next) call dump_vel
          if (ntimestep == ndumpforce_next) call dump_force(0)
          if (ntimestep == nrestart_next) call write_restart(0)
          if (ntimestep == ndiag_next) call diagnostic
          noutput_next = nthermo_next
          noutput_next = min(noutput_next,ndumpatom_next)
          noutput_next = min(noutput_next,ndumpvel_next)
          noutput_next = min(noutput_next,ndumpforce_next)
          noutput_next = min(noutput_next,nrestart_next)
          noutput_next = min(noutput_next,ndiag_next)
        endif

      enddo

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      time_loop = time1-time_loop

      return
      end
