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
c dump atom positions to a file
c iflag = 0/1 = called from integrator/minimizer

      subroutine dump_atom(iflag)
      use global
      use mpi
      implicit none

c argument variables

      integer iflag

c local variables

      integer iperatom,ntotal,nbuf,most
      integer ierror,j,i,nlocal_tmp,inode,irequest,itmp
      integer istatus(mpi_status_size)
      real*8 time1,time2
      real*8, allocatable :: buf(:)
      character*32 format_str

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()

c do not write this timestep if previously written, just update counter
c unless called from minimizer, then do dump regardless

      if (iflag == 0 .and. ntimestep.eq.ndumpatom_prev) then
        ndumpatom_next = min(ntimestep+ndumpatom,ntime_last)
        return
      endif

c format string depends on total # of atoms

      if (natoms < 1000) then
        format_str = "(i3,i3,3f9.5,3i5)"
      else if (natoms < 1000000) then
        format_str = "(i6,i3,3f9.5,3i5)"
      else
        format_str = "(i9,i3,3f9.5,3i5)"
      endif

c compute current ntotal = total number of atoms
c if atoms have been lost, could be different than natoms

      call mpi_allreduce(nlocal,ntotal,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

c write timestep header

      if (node.eq.0) then
        write (11,*) 'ITEM: TIMESTEP'
        write (11,*) ntimestep
        write (11,*) 'ITEM: NUMBER OF ATOMS'
        write (11,*) ntotal
        write (11,*) 'ITEM: BOX BOUNDS'
        write (11,*) box(1,1),box(2,1)
        write (11,*) box(1,2),box(2,2)
        if (idimension.eq.3) write (11,*) box(1,3),box(2,3)
        if (node.eq.0) write (11,*) 'ITEM: ATOMS'
      endif

c allocate a temporary buffer for the snapshot info
c big enough for most atoms on any proc
c iperatom = number of datums per atom, depends on trueflag

      iperatom = 5
      if (mod(trueflag,2).eq.1) iperatom = 8

      call mpi_allreduce(nlocal,most,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      nbuf = iperatom*most
      allocate(buf(nbuf))

c buffer up my atoms
c add trueflags to buffer if requested
c scale atom coords by box size so that they print as 0.0 -> 1.0
c IMPORTANT NOTE: an atom can be slightly outside the box 
c                 since PBC enforcement is only occasional

      if (mod(trueflag,2) == 0) then
        j = 0
        do i = 1,nlocal
          buf(j+1) = tag(i)
          buf(j+2) = type(i)
          buf(j+3) = (x(1,i)-box(1,1))/xprd
          buf(j+4) = (x(2,i)-box(1,2))/yprd
          buf(j+5) = (x(3,i)-box(1,3))/zprd
          j = j + iperatom
        enddo
      else
        j = 0
        do i = 1,nlocal
          buf(j+1) = tag(i)
          buf(j+2) = type(i)
          buf(j+3) = (x(1,i)-box(1,1))/xprd
          buf(j+4) = (x(2,i)-box(1,2))/yprd
          buf(j+5) = (x(3,i)-box(1,3))/zprd
          buf(j+6) = mod(true(i),1000) - 500
          buf(j+7) = mod(true(i)/1000,1000) - 500
          buf(j+8) = true(i)/1000000 - 500
          j = j + iperatom
        enddo
      endif
      nlocal_tmp = nlocal

c node 0 pings each node, receives their buffer, writes to file
c  all other nodes wait for ping, send buffer to node 0

      if (node == 0) then

        do inode = 0,nprocs-1
          if (inode.ne.0) then
            call mpi_irecv(buf,nbuf,mpi_double_precision,
     $           inode,0,mpi_comm_world,irequest,ierror)
            call mpi_send(itmp,0,mpi_integer,inode,0,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_double_precision,
     $           nlocal_tmp,ierror)
            nlocal_tmp = nlocal_tmp/iperatom
          endif

          j = 0
          if (iperatom.eq.5) then
            do i = 1,nlocal_tmp
              write (11,format_str) nint(buf(j+1)),nint(buf(j+2)),
     $             buf(j+3),buf(j+4),buf(j+5)
              j = j + iperatom
            enddo
          else
            do i = 1,nlocal_tmp
              write (11,format_str) nint(buf(j+1)),nint(buf(j+2)),
     $             buf(j+3),buf(j+4),buf(j+5),
     $             nint(buf(j+6)),nint(buf(j+7)),nint(buf(j+8))
              j = j + iperatom
            enddo
          endif
        enddo

      else

        call mpi_recv(itmp,0,mpi_integer,0,0,
     $       mpi_comm_world,istatus,ierror)
        call mpi_rsend(buf,iperatom*nlocal_tmp,mpi_double_precision,
     $       0,0,mpi_comm_world,ierror)

      endif

c free temporary buffer

      deallocate(buf)

c update timestep counters

      ndumpatom_next = min(ntimestep+ndumpatom,ntime_last)
      ndumpatom_prev = ntimestep

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1

      return
      end


c -------------------------------------------------------------------------
c open a atom dump file

      subroutine dumpatomopen
      use global
      implicit none

      if (node.eq.0) then
        open (11,file=dumpatomfile,status='unknown')
        close (11,status='delete')
        open (11,file=dumpatomfile,status='new')
      endif

c no atom dumps have been written into new file
      
      ndumpatom_prev = -1

      return
      end


c -------------------------------------------------------------------------
c close a atom dump file

      subroutine dumpatomclose
      use global
      implicit none

      if (node.eq.0) close (11)

      return
      end
