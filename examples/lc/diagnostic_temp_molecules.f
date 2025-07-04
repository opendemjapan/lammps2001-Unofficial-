c LAMMPS 2001 - Molecular Dynamics Simulator, Copyright 1998, 1999, 2000
c Authored by Steve Plimpton
c   (505) 845-7873, sjplimp@cs.sandia.gov
c   Dept 9221, MS 1111, Sandia National Labs, Albuquerque, NM 87185-1111
c See the README file for more information

c -------------------------------------------------------------------------
c user-specified diagnostic routines
c user needs to
c   (1) modify diag_init routine in 2 places (see ***)
c   (2) modify diagnostic routine in 1 place (see ***)
c   (3) add own routines for computing diagnostics (see samples at ***)
c       can put them in this file or in separate file(s)
c       if in separate files, have to add files to Makefile

c ----------------------------------------------------------------------
c initializization of total # of diagnostic routines and their names

      subroutine diag_init
      use global
      implicit none

c local variables

      integer idiag

c *** start of user-modified section
c set numdiag = # of diagnostic routines defined
c numdiag must be <= maxdiag setting in param.h

      numdiag = 1

c *** end of user-modified section

c error check and initialization

      if (numdiag.gt.maxdiag)
     $     call error('Too many diagnostic routines - boost maxdiag')

      do idiag = 1,numdiag
        ndiag(idiag) = 0
        diagfileflag(idiag) = 0
      enddo

c *** start of user-modified section
c assign name of each diagnostic routine from 1 to numdiag
c these names are used in diagnostic commands in input script
c each name must be a lowercase word (under_scores OK), 16 characters max

      diagnames(1) = 'temp_molecules'

c *** end of user-modified section

      return
      end


c ----------------------------------------------------------------------
c driver that calls user diagnostic routines

      subroutine diagnostic
      use global
      implicit none

c local variables

      integer idiag

c determine which user routines should be called this timestep
c set diagcall = 0 if not called, 1 if called
c only call if this diagnostic is currently enabled
c only call if diagnext = current timestep
c allow routine to be called even if previously called on this timestep,
c   this allows user to do setup at beginning of a new run,
c   diagnostic routine itself can choose to exit on this condition
c update diagnext to the next time this routine should be called

      do idiag = 1,numdiag
        diagcall(idiag) = 0
        if (ndiag(idiag).gt.0) then
          if (ntimestep.eq.diagnext(idiag)) then
            diagcall(idiag) = 1
            diagnext(idiag) = min(ntimestep+ndiag(idiag),ntime_last)
          endif
        endif
      enddo

c *** start of user-modified section
c call each user diagnostic routine from 1 to numdiag
c argument of called routine is same as index of diagcall
c routine names can be arbitrary
c routines must be provided by user

      if (diagcall(1).eq.1) call diag_temp_molecules(1)

c *** end of user-modified section

c ndiag_next = earliest timestep any diagnostic routine should next be called
c only check enabled diagnostics

      ndiag_next = ntime_last + 1
      do idiag = 1,numdiag
        if (ndiag(idiag).gt.0)
     $       ndiag_next = min(ndiag_next,diagnext(idiag))
      enddo

      return
      end


c ----------------------------------------------------------------------
c *** start of user-modified section

c sample user routines, called from "diagnostic" routine up above
c one routine for each diagnostic from 1 to numdiag
c idiag argument can be used to index diagnostic params and file info
c   diagnparams(idiag) = # of user-specified parameters for this diagnostic
c   diagparam(n,idiag) = 1 to n parameters for this diagnostic
c   diagfileflag(idiag) = 0,1 if diagnostic file is closed/open
c also, itime variable tells when in the run this routine is called
c   itime = 0, before run (from start routine)
c   itime < nsteps, during run (from integrate routine)
c   itime = nsteps, last call of run (from integrate routine)

c 1st example routine

      subroutine diag_temp_molecules(idiag)
      use mpi
      use global
      implicit none

c argument variables

      integer idiag

c local variables

      integer ifile,i,ierror
      real*8 t_molecules,t_other,tmp
      integer ncount,itmp

c check if this routine was called previously on same timestep
c could have occurred in a previous run
c if so, probably want to just exit

      if (ntimestep.eq.diagprev(idiag)) return

c ifile = handle for file opened when the "diagnostic" command was issued

      ifile = 20 + idiag

c set diagprev 

      diagprev(idiag) = ntimestep

c itime = 0 on 1st call to diagnostic in each run
c can do setup here if wish, error check on params, open own files, etc

      if (diagnparams(idiag) /= 2)
     $     call error('Incorrect # of params for diag_temp_molecules')

c -------------------------------------------
c user code to compute the desired diagnostic
c this is just an example of simple code

c write file header if this is 1st call to this diagnostic

      continue

c do a computation using params, write result to screen and log file

      lo_molecule = nint(diagparam(1,1))
      hi_molecule = nint(diagparam(2,1))

      t_molecules = 0.0
      t_other = 0.0
      ncount = 0
      do i = 1,nlocal
        if (molecule(i) >= lo_molecule .and.
     $       molecule(i) <= hi_molecule) then
          t_molecules = t_molecules +
     $         (v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)) *
     $         mass(type(i))
          ncount = ncount + 1
        else
          t_other = t_tother +
     $         (v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)) *
     $         mass(type(i))
        endif
      enddo

      call mpi_allreduce(ncount,itmp,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      ncount = itmp
      call mpi_allreduce(t_molecules,tmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      t_molecules = tmp / (boltz*idimension*ncount)
      call mpi_allreduce(t_other,tmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      t_other = tmp / (boltz*idimension*(natoms-ncount))

 900  format (a,i10,f12.4,f12.4)
      if (node == 0) then
        write (6,900) 'Temp of requested & other molecules =',
     $       ntimestep,t_molecules,t_other
        write (1,900) 'Temp of requested & other molecules =',
     $       ntimestep,t_molecules,t_other
      endif

c end of user code for the diagnostic
c -------------------------------------------

c itime = nsteps on final call to diagnostic in each run
c can do clean-up here if wish
c note that this call occurs within integrate,
c   so you may wish to have also computed the diagnostic

      continue

      return
      end

c *** end of user-modified section

c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c user should not modify routines from here to end-of-file

c ----------------------------------------------------------------------
c open diagnostic file

      subroutine diag_open(idiag)
      use global
      implicit none

c argument variables

      integer idiag

c local variables

      integer ifile

      ifile = 20 + idiag

      if (node.eq.0) then
        open (ifile,file=diagfile(idiag),status='unknown')
        close (ifile,status='delete')
        open (ifile,file=diagfile(idiag),status='new')
      endif

c no diagnostics have been written into new file
      
      diagprev(idiag) = -1

      return
      end


c ----------------------------------------------------------------------
c close diagnostic file

      subroutine diag_close(idiag)
      use global
      implicit none

c argument variables

      integer idiag

c local variables

      integer ifile

      ifile = 20 + idiag
      if (node.eq.0) close (ifile)

      return
      end


c ----------------------------------------------------------------------
c lookup the name of a diagnostic and return an index
c error return of 0 means name not found

      subroutine diag_lookup(name,idiag)
      use global
      implicit none

c argument variables

      character*(*) name
      integer idiag

      do idiag = 1,numdiag
        if (name.eq.diagnames(idiag)) return
      enddo

      idiag = 0

      return
      end
