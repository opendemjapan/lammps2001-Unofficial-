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

c Contributing author: Mark Stevens (Sandia)

c -------------------------------------------------------------------------
c 4-level timestep integrator using rRESPA
c outer loop (long) is conceptually equivalent to regular timestep
c  in standard integrator
c 3 inner-nested timestep loops (short, intra, stretch)
c  conceptually replace nve_x and reneighboring in standard integrator
c communicate routine is only called once - in inner loop which is the
c  only place positions are ever updated
c reverse communication occurs whenever forces are computed, within each
c  of the 4 force routines
c neighbor check only needed at nonbond short-range level
c  prior to force_short computation since bonded forces are correct 
c  no matter how far atoms move
c virial flag is only turned on for last inner loop iterations
c NVT integrator is coupled to nonbond short-range forces

      subroutine integrate_respa
      use global
      use mpi
      implicit none

c local variables

      integer vflag,v1flag,ierror,jshort,jintra,jstretch
      integer nflag
      real*8 p_omega_dt,time1,time2,time3,time4,time5,target
      dimension p_omega_dt(3)

      call mpi_barrier(mpi_comm_world,ierror)
      time_loop = mpi_wtime()

c long-range loop
c 1st half of long-range update
c inner iteration on short,intra,bond
c long-range forces
c 2nd half of long-range update

      do itime = 1,nsteps

        ntimestep = ntimestep + 1

c should virial be computed this timestep

        vflag = 0
        if (ensemble >= 3) then
          vflag = 1
        else if (nthermo_next == ntimestep) then
          vflag = 1
        endif

c 1st half of update

        if (ensemble == 2) then
          call nvt_eta(dthalf_long,1)
          call nvt_v(f_long,dthalf_long,1)
        else
          call nve_v(f_long,dthalf_long)
        endif

c nonbond loop
c 1st half of nonbond update
c inner iteration on intra,bond
c neighboring
c nonbond forces
c 2nd half of nonbond update

        do jshort = 1,nshort 

          if (ensemble <= 2) then
            call nve_v(f_short,dthalf_short)
          else if (ensemble == 3) then
            call nph_eta
            call npt_omegadot(dthalf_short,p_omega_dt)
            call npt_v(f_short,dthalf_short,1)
            call npt_x(dthalf_short,p_omega_dt,1)
          else if (ensemble == 4) then
            call nvt_eta(dthalf_short,1)
            call npt_omegadot(dthalf_short,p_omega_dt)
            call npt_v(f_short,dthalf_short,1)
            call npt_x(dthalf_short,p_omega_dt,1)
          endif

c 3/4-body loop
c 1st half of 3/4-body update
c inner iteration on bond
c 3/4-body forces
c 2nd half of 3/4-body update

          do jintra = 1,nintra

            call nve_v(f_intra,dthalf_intra)

c bond loop
c 1st half of bond update
c communicate new positions
c bond forces
c 2nd half of bond update

            do jstretch = 1,nstretch

              call nve_v(f_stretch,dthalf)
              call nve_x

#ifdef SYNC
              call mpi_barrier(mpi_comm_world,ierror)
#endif
              time1 = mpi_wtime()
              call communicate
              time2 = mpi_wtime()
              time_comm = time_comm + time2-time1

              v1flag = vflag
              if (jstretch < nstretch) v1flag = 0
              call force_stretch(v1flag)

              if (nfixes_respa > 0) call fix_apply(0,f_stretch,dthalf)
              call nve_v(f_stretch,dthalf)

            enddo

            v1flag = vflag
            if (jintra < nintra) v1flag = 0
            call force_intra(v1flag)

            if (nfixes_respa > 0) call fix_apply(0,f_intra,dthalf_intra)
            call nve_v(f_intra,dthalf_intra)

          enddo

c force neighbor list rebuild on steps restart files are written
c   this (1) insures atoms in file will have PBC applied and
c        (2) enables exact restarts since atoms will migrate to new procs
c            before being written to file (same as when file is read)
c else check for neighbor list rebuild if neighdelay steps have occurred
          
          nflag = 0
          neighago = neighago + 1

          if (ntimestep == nrestart_next .and. jshort == nshort) then
            nflag = 1
          else if (neighago >= neighdelay .and.
     $           mod(neighago,neighfreq) == 0) then
            if (neightrigger == 0) then
              nflag = 1
            else
              call check_neighbor(nflag)
            endif
          endif

c neighbor list rebuild w/ communication
c (1) remap atoms into periodic box
c (2) shrink-wrap box around atoms if any dims are nonperiodic
c (3) reset comm and neighbor-bin patterns if box size has changed
c (4) move atoms to new procs via exchange routine
c (5) acquire of ghost atoms via borders routine
c (6) build all new neighbor lists
          
          if (nflag == 1) then
#ifdef SYNC
            call mpi_barrier(mpi_comm_world,ierror)
#endif
            time1 = mpi_wtime()
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
            call neighbor
            time4 = mpi_wtime()
            if (nbonds > 0) call neighbor_bond
            if (nangles > 0) call neighbor_angle
            if (ndihedrals > 0) call neighbor_dihedral
            if (nimpropers > 0) call neighbor_improper
            time5 = mpi_wtime()
            time_exch = time_exch + time2-time1
            time_comm = time_comm + time3-time2
            time_neigh1 = time_neigh1 + time4-time3
            time_neigh2 = time_neigh2 + time5-time4
          endif

          call force_short(vflag)

c global temperature control via Langevin drag
c applied on short timestep

          if (tempstyle == 3) call langevin(f_short,dthalf_short)
          if (nfixes > 0) call fix_apply(1,f_short,dthalf_short)

          if (ensemble <= 2) then
            call nve_v(f_short,dthalf_short)
          else if (ensemble == 3) then
            call npt_v(f_short,dthalf_short,2)
            call npt_x(dthalf_short,p_omega_dt,2)
            call temperature
            call nph_eta
            call pressure
            call npt_omegadot(dthalf_short,p_omega_dt)
          else if (ensemble == 4) then
            call npt_v(f_short,dthalf_short,2)
            call npt_x(dthalf_short,p_omega_dt,2)
            call temperature
            call nvt_eta(dthalf_short,2)
            call pressure
            call npt_omegadot(dthalf_short,p_omega_dt)
          endif

        enddo

        call force_long(vflag)

        if (nfixes_respa > 0) call fix_apply(0,f_long,dthalf_long)

c 2nd half of update

        if (ensemble == 2) then
          call nvt_v(f_long,dthalf_long,2)
          call temperature
          call nvt_eta(dthalf_long,2)
        else
          call nve_v(f_long,dthalf_long)
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
