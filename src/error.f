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
c parallel error handler routine called from all procs
c close all files before exit

      subroutine err(str)
      use global
      use mpi
      implicit none

c argument variables

      character*(*) str

c local variables 

      integer idiag,ierror
      real*8 time

c close all files

      if (dumpatomfileflag.gt.0) call dumpatomclose
      if (dumpvelfileflag.gt.0) call dumpvelclose
      if (dumpforcefileflag.gt.0) call dumpforceclose
      do idiag = 1,numdiag
        if (diagfileflag(idiag).gt.0) call diag_close(idiag)
      enddo
      if (node.eq.0) close (1)

c final timer

      call mpi_barrier(mpi_comm_world,ierror)
      time = mpi_wtime()
      time_total = time-time_total

      if (str(1:8).eq.'All Done') then
        if (node.eq.0) write (6,*) 'Total time =',time_total
      else
        if (node.eq.0) write (6,*) 'ERROR: ',str
      endif

c close down MPI

      call mpi_finalize(ierror)

      call exit(0)

      return
      end


c -------------------------------------------------------------------------
c serial error handler routine called from one proc

      subroutine abort(str)
      use global
      use mpi
      implicit none

c argument variables

      character*(*) str

c local variables 

      integer ierror

      write (6,*) 'Abort: ',str

c mpi_abort should cause all procs to halt

      call mpi_abort(mpi_comm_world,1,ierror)
      call exit(0)

      return
      end
