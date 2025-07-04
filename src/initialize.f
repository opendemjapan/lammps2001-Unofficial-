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
c initialize run parameters

      subroutine initialize
      use global
      use mpi
      implicit none

c local variables

      integer ierror,i,j

c identify self

      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,node,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)

c print out version and determine starting time
      
      iversion = 2001
      if (node == 0) then
        time_total = mpi_wtime()
        write (6,*) 'LAMMPS 2001 (Nov 2001)'
      endif
      call mpi_bcast(time_total,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

c open log file

      if (node == 0) then
        open (unit=1,file='log.lammps',status='unknown')
        close (1,status='delete')
        open (unit=1,file='log.lammps',status='new')
        write (1,*) 'LAMMPS 2001 (Nov 2001)'
      endif

c constants - initial default is conventional units
c efactor is now in agreement with CHARMM

      boltz = 0.001987191D0
      dtfactor = 48.88821D0
      pfactor = 68589.796D0
      efactor = 332.0636D0
      two_1_3 = 2.0D0**(1.0D0/3.0D0)

c values that must be set before read_data or read_restart commands

      units = 0

      extra_own = 1.5
      if (nprocs == 1) extra_own = 1.0

      extra_ghost = 1.5
      extra_neigh = 1.5
      extra_buf = 1.5

      idimension = 3
      pgrid(1) = 0
      pgrid(2) = 0
      pgrid(3) = 0
      perflagx = 0
      perflagy = 0
      perflagz = 0
      nonperiodic = .FALSE.
      newton = 3
      newton_bond = 1
      newton_nonbond = 1
      trueflag = 0
      cutmax = 0.0D0
      restart_version = iversion

      skin = 2.0D0
      neighstyle = 1
      neighfreq = 1
      neighdelay = 10
      neightrigger = 1

c integrator settings

      dt = 1.0D0
      nrespa = 0

c 1-2, 1-3, 1-4 neighbor weightings
c CHARMM default

      amberflag = 0
      special(1) = 0.0D0
      special(2) = 0.0D0
      special(3) = 0.0D0

c output settings

      nthermo = 0
      thermostyle = 0

      ndumpatom = 0
      dumpatomfileflag = 0
      ndumpvel = 0
      dumpvelfileflag = 0
      ndumpforce = 0
      dumpforcefileflag = 0

      nrestart = 0

      call diag_init

c ensemble settings

      tempstyle = 0
      pressstyle = 0
      volstyle = 0

c force field settings
c cutlj in conventional units

      nonstyle = 1
      cutlj = 10.0D0
      offsetflag = 0
      mixflag = 0
      mixstyle = 1

      coulstyle = 1
      cutcoul = 10.0D0
      meshflag = 0
      orderflag = 5
      dielectric = 1.0D0

      bondstyle = 1
      anglestyle = 1
      dihedstyle = 1
      improstyle = 1

c no slab correction for Ewald/PPPM

      slabflag = 0
      slab_volfactor = 1.0

c no group yet chosen for velocity creation

      rotationflag = 0
      creategroup = 0

c no fixes

      nfixes = 0
      nfixes_respa = 0
      nshake = 0
      nshakebonds = 0

c initial Ewald array size

      kmax = 0
      
c optimization settings

      optstyle = 1
      optflag = 1

c no comm & neighbor arrays allocated

      maxswap = 0
      maxbin = 0

c no commands have been read-in yet

      firstflag = 0
      readflag = 0

      return
      end
