cLAMMPS 2001 - Molecular Dynamics Simulator
cSandia National Laboratories, www.cs.sandia.gov/~sjplimp/lammps.html
cSteve Plimpton, sjplimp@sandia.gov
c
cCopyright (1998-2001) Sandia Corporation.  Under the terms of Contract
cDE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
ccertain rights in this software.  This software is distributed under 
cthe GNU General Public License.
c
cSee the README file in the top-level LAMMPS directory.

c -------------------------------------------------------------------------
c write a restart file in native binary format
c iflag = 0/1 = called from integrator/minimizer

      subroutine write_restart(iflag)
      use global
      use mpi
      implicit none

c argument variables

      integer iflag

c local variables

      integer, parameter :: iperatom = 16
      integer ierror,nbyte,ibyte,i,j,jj,icnt,inode,ilen,n,nbuf,icount
      integer irequest,itmp,nlocal_tmp,k,numbond_tmp
      integer numangle_tmp,numdihed_tmp,numimpro_tmp
      real*8 time1,time2
      real*8, allocatable :: buf(:)
      integer istatus(mpi_status_size)
      character*80 filename,suffix
      
      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()

c create restart filename and open it
c style 1 = append timestep to restart_out
c style 2 = toggle between restart_out1 and restart_out2

      if (iflag == 1) then
        call writeint(ntimestep,suffix,ilen)
        filename = trim(restart_out)//'.'//suffix(1:ilen)
        filename = trim(filename)//'.min'
      else if (restartstyle == 1) then
        call writeint(ntimestep,suffix,ilen)
        filename = trim(restart_out)//'.'//suffix(1:ilen)
      else
        if (restartlast == 2) then
          filename = trim(restart_out1)
          restartlast = 1
        else
          filename = trim(restart_out2)
          restartlast = 2
        endif
      endif

      if (node == 0) open (unit=2,file=filename,
     $     form='unformatted',status='replace')

c write header
c include all parameters that must be set before "read restart" command
c   so can check for match when restart is done
c include all values that insure "perfect" restart capability

      if (node.eq.0) then
        write (2) iversion
        write (2) ntimestep,nprocs,pgrid(1),pgrid(2),pgrid(3)
        write (2) idimension,perflagx,perflagy,perflagz
        write (2) units,newton
        write (2) natoms,nbonds,nangles,ndihedrals,nimpropers
        write (2) ntypes
        if (nbonds.gt.0) write (2) nbondtypes
        if (nangles.gt.0) write (2) nangletypes
        if (ndihedrals.gt.0) write (2) ndihedtypes
        if (nimpropers.gt.0) write (2) nimprotypes
        write (2) maxbondper,maxangleper,maxdihedper,maximproper
        write (2) box
        write (2) eta,eta_dot
        write (2) omega(1),omega(2),omega(3),
     $       omega_dot(1),omega_dot(2),omega_dot(3)
      endif

c write all coefficients

      if (node.eq.0) then
        
        write (2) (mass(i),i=1,ntypes)

        write (2) nonstyle
        if (nonstyle.eq.1) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.2) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.3) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.4) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.5) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.6) then
          write (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff14_1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff14_2(i,j),j=i,ntypes),i=1,ntypes)
        endif

        if (nbonds.gt.0) then
          write (2) bondstyle
          if (bondstyle.eq.1) then
            write (2) ((bondcoeff(i,j),i=1,2),j=1,nbondtypes)
          else if (bondstyle.eq.2) then
            write (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          else if (bondstyle.eq.3) then
            write (2) ((bondcoeff(i,j),i=1,5),j=1,nbondtypes)
          else if (bondstyle.eq.4) then
            write (2) ((bondcoeff(i,j),i=1,3),j=1,nbondtypes)
          else if (bondstyle.eq.5) then
            write (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          endif
        endif

        if (nangles.gt.0) then
          write (2) anglestyle
          if (anglestyle.eq.1) then
            write (2) ((anglecoeff(i,j),i=1,2),j=1,nangletypes)
          else if (anglestyle.eq.2) then
            write (2) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
            write (2) ((bondbondcoeff(i,j),i=1,3),j=1,nangletypes)
            write (2) ((bondanglecoeff(i,j),i=1,4),j=1,nangletypes)
          else if (anglestyle.eq.3) then
            write (2) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
          else if (anglestyle.eq.4) then
            write (2) ((anglecoeff(i,j),i=1,1),j=1,nangletypes)
          endif
        endif

        if (ndihedrals.gt.0) then
          write (2) dihedstyle
          if (dihedstyle.eq.1) then
            write (2) ((dihedcoeff(i,j),i=1,3),j=1,ndihedtypes)
          else if (dihedstyle.eq.2) then
            write (2) ((dihedcoeff(i,j),i=1,6),j=1,ndihedtypes)
            write (2) ((midbondtorsioncoeff(i,j),i=1,4),
     $           j=1,ndihedtypes)
            write (2) ((endbondtorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            write (2) ((angletorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            write (2) ((angleangletorsioncoeff(i,j),i=1,3),
     $           j=1,ndihedtypes)
            write (2) ((bondbond13coeff(i,j),i=1,3),
     $           j=1,ndihedtypes)
          else if (dihedstyle.eq.3) then
            write (2) ((dihedcoeff(i,j),i=1,5),j=1,ndihedtypes)
          else if (dihedstyle.eq.4) then
            write (2) ((dihedcoeff(i,j),i=1,4),j=1,ndihedtypes)
          endif
        endif

        if (nimpropers.gt.0) then
          write (2) improstyle
          if (improstyle.eq.1) then
            write (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
          else if (improstyle.eq.2) then
            write (2) ((improcoeff(i,j),i=1,3),j=1,nimprotypes)
          else if (improstyle.eq.3) then
            write (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
            write (2) ((angleanglecoeff(i,j),i=1,6),j=1,nimprotypes)
          endif
        endif

      endif

c allocate a temporary buffer for the restart info
c big enough for max atoms on any proc + 1 leading loc for # of atoms

      n = 1
      do i = 1,nlocal
        n = n + iperatom + 3*numbond(i) + 4*numangle(i) +
     $       5*numdihed(i) + 5*numimpro(i)
      enddo

      call mpi_allreduce(n,nbuf,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      allocate(buf(nbuf))

c buffer up my atoms
c buf(1) = # of atoms

      j = 1
      buf(j) = nlocal
      do i = 1,nlocal
        buf(j+1) = x(1,i)
        buf(j+2) = x(2,i)
        buf(j+3) = x(3,i)
        buf(j+4) = tag(i)
        buf(j+5) = molecule(i)
        buf(j+6) = type(i)
        buf(j+7) = q(i)
        buf(j+8) = v(1,i)
        buf(j+9) = v(2,i)
        buf(j+10) = v(3,i)
        buf(j+11) = true(i)
        buf(j+12) = fix(i)
        buf(j+13) = numbond(i)
        buf(j+14) = numangle(i)
        buf(j+15) = numdihed(i)
        buf(j+16) = numimpro(i)
        j = j + iperatom
        do jj = 1,numbond(i)
          buf(j+1) = bondatom1(jj,i)
          buf(j+2) = bondatom2(jj,i)
          buf(j+3) = bondtype(jj,i)
          j = j + 3
        enddo
        do jj = 1,numangle(i)
          buf(j+1) = angleatom1(jj,i)
          buf(j+2) = angleatom2(jj,i)
          buf(j+3) = angleatom3(jj,i)
          buf(j+4) = angletype(jj,i)
          j = j + 4
        enddo
        do jj = 1,numdihed(i)
          buf(j+1) = dihedatom1(jj,i)
          buf(j+2) = dihedatom2(jj,i)
          buf(j+3) = dihedatom3(jj,i)
          buf(j+4) = dihedatom4(jj,i)
          buf(j+5) = dihedtype(jj,i)
          j = j + 5
        enddo
        do jj = 1,numimpro(i)
          buf(j+1) = improatom1(jj,i)
          buf(j+2) = improatom2(jj,i)
          buf(j+3) = improatom3(jj,i)
          buf(j+4) = improatom4(jj,i)
          buf(j+5) = improtype(jj,i)
          j = j + 5
        enddo
      enddo
      icount = j

c node 0 pings each node, receives their buffer, writes to file atom by atom
c  all other nodes wait for ping, send buffer to node 0

      if (node == 0) then

        do inode = 0,nprocs-1
          if (inode /= 0) then
            call mpi_irecv(buf,nbuf,mpi_double_precision,
     $           inode,0,mpi_comm_world,irequest,ierror)
            call mpi_send(itmp,0,mpi_integer,inode,0,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
          endif

          j = 1
          nlocal_tmp = nint(buf(j))
          do i = 1,nlocal_tmp
            write (2) (buf(j+k),k=1,iperatom)
            j = j + iperatom
            numbond_tmp = nint(buf(j-3))
            numangle_tmp = nint(buf(j-2))
            numdihed_tmp = nint(buf(j-1))
            numimpro_tmp = nint(buf(j))
            do jj = 1,numbond_tmp
              write (2) buf(j+1),buf(j+2),buf(j+3)
              j = j + 3
            enddo
            do jj = 1,numangle_tmp
              write (2) buf(j+1),buf(j+2),buf(j+3),buf(j+4)
              j = j + 4
            enddo
            do jj = 1,numdihed_tmp
              write (2) buf(j+1),buf(j+2),buf(j+3),
     $             buf(j+4),buf(j+5)
              j = j + 5
            enddo
            do jj = 1,numimpro_tmp
              write (2) buf(j+1),buf(j+2),buf(j+3),
     $             buf(j+4),buf(j+5)
              j = j + 5
            enddo
          enddo

        enddo

      else

        call mpi_recv(itmp,0,mpi_integer,0,0,
     $       mpi_comm_world,istatus,ierror)
        call mpi_rsend(buf,icount,mpi_double_precision,
     $       0,0,mpi_comm_world,ierror)

      endif

c free temporary buffer

      deallocate(buf)

c close restart file

      if (node == 0) close (2)

c update timestep counter

      nrestart_next = min(ntimestep+nrestart,ntime_last)

      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      time_io = time_io + time2-time1
      
      return
      end
