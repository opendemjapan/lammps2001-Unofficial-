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
c dump atom velocities to a file

      subroutine dump_vel
      use global
      use mpi
      implicit none

c local variables

      integer, parameter :: iperatom = 5
      integer ntotal
      integer ierror,j,i,nlocal_tmp,inode,irequest,itmp,most,nbuf
      integer istatus(mpi_status_size)
      real*8 time1,time2,rinverse
      real*8, allocatable :: buf(:)
      character*32 format_str

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()

c do not write this timestep if previously written, just update counter

      if (ntimestep.eq.ndumpvel_prev) then
        ndumpvel_next = min(ntimestep+ndumpvel,ntime_last)
        return
      endif

c format string depends on total # of atoms

      if (natoms < 1000) then
        format_str = "(i3,i3,3f10.5)"
      else if (natoms < 1000000) then
        format_str = "(i6,i3,3f10.5)"
      else
        format_str = "(i9,i3,3f10.5)"
      endif

c compute current ntotal = total number of atoms
c if atoms have been lost, could be different than natoms

      call mpi_allreduce(nlocal,ntotal,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

c write timestep header

      if (node.eq.0) then
        write (12,*) 'ITEM: TIMESTEP'
        write (12,*) ntimestep
        write (12,*) 'ITEM: NUMBER OF ATOMS'
        write (12,*) ntotal
        write (12,*) 'ITEM: BOX BOUNDS'
        write (12,*) box(1,1),box(2,1)
        write (12,*) box(1,2),box(2,2)
        write (12,*) box(1,3),box(2,3)
        write (12,*) 'ITEM: VELOCITIES'
      endif

c allocate a temporary buffer for the snapshot info
c big enough for most atoms on any proc

      call mpi_allreduce(nlocal,most,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      nbuf = iperatom*most
      allocate(buf(nbuf))
      
c buffer up my atom velocities (scaled by dtfactor)

      rinverse = 1.0/dtfactor
      j = 0
      do i = 1,nlocal
        buf(j+1) = tag(i)
        buf(j+2) = type(i)
        buf(j+3) = v(1,i)*rinverse
        buf(j+4) = v(2,i)*rinverse
        buf(j+5) = v(3,i)*rinverse
        j = j + iperatom
      enddo
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
          do i = 1,nlocal_tmp
            write (12,format_str) nint(buf(j+1)),nint(buf(j+2)),
     $           buf(j+3),buf(j+4),buf(j+5)
            j = j + iperatom
          enddo
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

      ndumpvel_next = min(ntimestep+ndumpvel,ntime_last)
      ndumpvel_prev = ntimestep

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1

      return
      end


c -------------------------------------------------------------------------
c open a velocity dump file

      subroutine dumpvelopen
      use global
      implicit none

      if (node.eq.0) then
        open (12,file=dumpvelfile,status='unknown')
        close (12,status='delete')
        open (12,file=dumpvelfile,status='new')
      endif

c no velocity dumps have been written into new file
      
      ndumpvel_prev = -1

      return
      end


c -------------------------------------------------------------------------
c close a velocity dump file

      subroutine dumpvelclose
      use global
      implicit none

      if (node.eq.0) close (12)

      return
      end
