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

c -----------------------------------------------------------------------
c interprocessor communication via 6-way stencil

c -------------------------------------------------------------------------
c swap slabs of atoms with other processors in all 3 directions
c done every timestep
c send/receive message always to another proc - never to self
      
      subroutine communicate
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,m,ii,ierror
      integer iflag,jflag,kflag,isize_send,isize_recv
      integer irequest,istatus(mpi_status_size)

      do k = 1,nswap

c sending to another proc, pack up message buffer with current atom coords

        if (spart(k) /= node) then

          j = 0

c not sending across periodic boundary

          if (commflagall(k) == 0) then
            do ii = nsfirst(k),nslast(k)
              i = slist(ii)
              buf1(j+1) = x(1,i)
              buf1(j+2) = x(2,i)
              buf1(j+3) = x(3,i)
              j = j + 3
            enddo

          else

c sending across periodic boundary, add/subtract box length

            iflag = commflag(1,k)
            jflag = commflag(2,k)
            kflag = commflag(3,k)
            do ii = nsfirst(k),nslast(k)
              i = slist(ii)
              buf1(j+1) = x(1,i) + iflag*xprd
              buf1(j+2) = x(2,i) + jflag*yprd
              buf1(j+3) = x(3,i) + kflag*zprd
              j = j + 3
            enddo

          endif

c send/receive message

          isize_send = nslast(k) - nsfirst(k) + 1
          isize_recv = nrlast(k) - nrfirst(k) + 1

#ifdef SENDRECV
          call mpi_sendrecv(buf1,3*isize_send,
     $         mpi_double_precision,spart(k),0,
     $         x(1,nrfirst(k)),3*isize_recv,
     $         mpi_double_precision,rpart(k),0,
     $         mpi_comm_world,istatus,ierror)
#endif
#ifndef SENDRECV
          call mpi_irecv(x(1,nrfirst(k)),3*isize_recv,
     $         mpi_double_precision,rpart(k),0,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf1,3*isize_send,
     $         mpi_double_precision,spart(k),0,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
#endif

        else

c only swapping with self, avoid buffering, do direct copy

          m = nrfirst(k)
          iflag = commflag(1,k)
          jflag = commflag(2,k)
          kflag = commflag(3,k)
          do ii = nsfirst(k),nslast(k)
            i = slist(ii)
            x(1,m) = x(1,i) + iflag*xprd
            x(2,m) = x(2,i) + jflag*yprd
            x(3,m) = x(3,i) + kflag*zprd
            m = m + 1
          enddo

        endif
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c reverse swap slabs of forces with other processors in all 3 directions
c sum to existing forces
c done every timestep
c send/receive message always to another proc - never to self
      
      subroutine reverse_comm
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,m,ii,ierror
      integer iflag,jflag,kflag,isize_send,isize_recv
      integer irequest,istatus(mpi_status_size)

      do k = nswap,1,-1

c sending to another proc, pack up message buffer with current forces

        if (spart(k) /= node) then

c send/receive message

          isize_send = nrlast(k) - nrfirst(k) + 1
          isize_recv = nslast(k) - nsfirst(k) + 1

#ifdef SENDRECV
          call mpi_sendrecv(f(1,nrfirst(k)),3*isize_send,
     $         mpi_double_precision,rpart(k),0,
     $         buf1,3*isize_recv,
     $         mpi_double_precision,spart(k),0,
     $         mpi_comm_world,istatus,ierror)
#endif
#ifndef SENDRECV
          call mpi_irecv(buf1,3*isize_recv,
     $         mpi_double_precision,spart(k),0,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(f(1,nrfirst(k)),3*isize_send,
     $         mpi_double_precision,rpart(k),0,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
#endif

c add received forces to existing forces

          j = 0
          do ii = nsfirst(k),nslast(k)
            i = slist(ii)
            f(1,i) = f(1,i) + buf1(j+1)
            f(2,i) = f(2,i) + buf1(j+2)
            f(3,i) = f(3,i) + buf1(j+3)
            j = j + 3
          enddo

        else

c only swapping with self, avoid buffering, do direct summation

          m = nrfirst(k)
          do ii = nsfirst(k),nslast(k)
            i = slist(ii)
            f(1,i) = f(1,i) + f(1,m)
            f(2,i) = f(2,i) + f(2,m)
            f(3,i) = f(3,i) + f(3,m)
            m = m + 1
          enddo

        endif
      enddo

      return
      end


c -------------------------------------------------------------------------
c send out atoms that have left my box, receive ones entering my box
c  since last reneighboring
c done in all 3 directions
c all info associated with an atom migrates with the atom
c called right before each reneighboring
      
      subroutine exchange
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,m,n,jj,nexch,irequest,ierror,icnt,itmp
      real*8 blo,bhi,coord
      integer istatus(mpi_status_size)
      integer packatom

c clear all atoms from global list before exchanged atoms migrate

      do i = 1,nlocal+nghost
        localptr(tag(i)) = 0
      enddo

c loop over dimensions

      do k = 1,3
        
        blo = border(1,k)
        bhi = border(2,k)
        
c fill buffer with atoms leaving my box, update local list
c use copyatom to replace any atom put in buffer with one from end of list

        nexch = 0
        j = 1
        i = 1
        do while (i <= nlocal)
          if (x(k,i) < blo .or. x(k,i) >= bhi) then
            nexch = nexch + 1
            if (j <= maxbuf) j = j + packatom(i,buf1(j))
            call copyatom(nlocal,i)
            nlocal = nlocal - 1
          else
            i = i + 1
          endif
        enddo
        j = j - 1
        
        max_exch = max(max_exch,nexch)
        max_buf = max(max_buf,2*j)
        
        if (j > maxbuf/2) then
          write (6,*) 'Too many atoms sent -- ',
     $         'boost extra_buf:',node,ntimestep,nexch,maxbuf
          call abort('Buffer overflow')
        endif

c send/recv atoms to neighbors
c if no neighbors, just set icnt = 0, this effectively deletes atoms
c   that have left box
c if same neighbor on both sides (pgrid = 2), just send atoms once
c if 2 neighbors send atoms both directions

        if (pgrid(k) == 1) then
          icnt = 0
        else

#ifdef SENDRECV
          call mpi_sendrecv(buf1,j,mpi_double_precision,mpart(1,k),0,
     $         buf2,maxbuf,mpi_double_precision,mpart(2,k),0,
     $         mpi_comm_world,istatus,ierror)
#endif
#ifndef SENDRECV
          call mpi_irecv(buf2,maxbuf,
     $         mpi_double_precision,mpart(2,k),0,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf1,j,mpi_double_precision,mpart(1,k),0,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
#endif

          call mpi_get_count(istatus,mpi_double_precision,
     $         icnt,ierror)
        endif

        if (pgrid(k) > 2) then

#ifdef SENDRECV
          call mpi_sendrecv(buf1,j,mpi_double_precision,mpart(2,k),0,
     $         buf2(icnt+1),maxbuf-icnt,mpi_double_precision,
     $         mpart(1,k),0,mpi_comm_world,istatus,ierror)
#endif
#ifndef SENDRECV             
          call mpi_irecv(buf2(icnt+1),maxbuf-icnt,
     $         mpi_double_precision,mpart(1,k),0,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf1,j,mpi_double_precision,mpart(2,k),0,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
#endif

          call mpi_get_count(istatus,mpi_double_precision,
     $         itmp,ierror)
          icnt = icnt + itmp
        endif
        
c check incoming atoms to see if they are in my box
c if they are, add to my list

        j = 1
        do while (j <= icnt)
          coord = buf2(j+k)
          if (coord >= blo .and. coord < bhi) then
            nlocal = nlocal + 1
            if (nlocal <= maxown) call unpackatom(nlocal,buf2(j))
          endif
          j = j + nint(buf2(j))
        enddo
        
        if (nlocal > maxown) then
          write (6,*) 'Too many atoms recvd -- ',
     $         'boost extra_own:',node,ntimestep,nlocal,maxown
          call abort('Atom overflow')
        endif
        
      enddo

      return
      end


c -------------------------------------------------------------------------
c make lists of nearby atoms to send to neighboring nodes at every timestep
c one list made for every swap that will be made
c as list is made, actually do swaps
c this does equivalent of a communicate (so don't need to explicitly
c  call communicate routine on reneighboring timestep)
c tag, type, charge, 3 atom positions are what is swapped
c  since they are needed for force computation
c this routine is called before every reneighboring
      
      subroutine borders
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,m,n,kk,iflag,jflag,kflag
      integer nrecv,maxrecv,npnt,nsend,nstatus,ierror
      real*8 blo,bhi
      integer irequest(4),istatus(mpi_status_size,4)

      nswap = 0
      nghost = 0
      npnt = 0

c nstatus = # of messages to be exchanged, one less if no charge in system

      nstatus = 4
      if (ncharge == 0) nstatus = 3
      
      do k = 1,3
        do kk = 0,2*need(k)-1
          
c buffer up all atoms (own + ghost) that are inside slab boundaries
c store atoms in slist for communicate routine to use in future timesteps
          
          nswap = nswap + 1
          nsfirst(nswap) = npnt + 1
          nrfirst(nswap) = nlocal + nghost + 1
          blo = slablo(nswap)
          bhi = slabhi(nswap)
          iflag = commflag(1,nswap)
          jflag = commflag(2,nswap)
          kflag = commflag(3,nswap)

          nsend = 0
          j = 0
          do i = 1,nlocal+nghost
            if (x(k,i) >= blo .and. x(k,i) < bhi) then
              npnt = npnt + 1
              if (npnt <= maxghost) slist(npnt) = i
              nsend = nsend + 1
              if (j <= maxbuf) then
                buf1(j+1) = x(1,i) + iflag*xprd
                buf1(j+2) = x(2,i) + jflag*yprd
                buf1(j+3) = x(3,i) + kflag*zprd
                ibuf1(nsend) = tag(i)
                ibuf2(nsend) = type(i)
                if (ncharge == 1) buf2(nsend) = q(i)
              endif
              j = j + 3
            endif
          enddo

          max_bord = max(max_bord,nsend)
          max_buf = max(max_buf,j)
          
          if (npnt > maxghost) then
            write (6,*) 'Too many ghosts sent -- ',
     $           'boost extra_ghost:',node,ntimestep,npnt,maxghost
            call abort('Ghost atom overflow')
          endif
          
          if (j > maxbuf) then
            write (6,*) 'Too many border atoms -- ',
     $           'boost extra_buf:',node,ntimestep,nsend,maxbuf
            call abort('Buffer overflow')
          endif

c swap atoms, put incoming ghosts at end of my position array
c if swapping with self, simply copy, do not send messages
c maxrecv = max ghosts I can recv
c check on it by handshaking first since MPI may otherwise crash

          maxrecv = maxghost - nghost

          if (spart(nswap) /= node) then

            call mpi_send(nsend,1,mpi_integer,spart(nswap),0,
     $           mpi_comm_world,ierror)
            call mpi_recv(nrecv,1,mpi_integer,rpart(nswap),0,
     $           mpi_comm_world,istatus,ierror)

            if (nrecv > maxrecv) then
              write (6,*) 'Too many ghosts recvd -- ',
     $             'boost extra_ghost:',
     $             node,ntimestep,nghost+nrecv,maxghost
              call abort('Ghost atom overflow')
            endif

#ifdef SENDRECV
            call mpi_sendrecv(buf1,3*nsend,
     $           mpi_double_precision,spart(nswap),0,
     $           x(1,nrfirst(nswap)),3*nrecv,
     $           mpi_double_precision,rpart(nswap),0,
     $           mpi_comm_world,istatus,ierror)
            call mpi_sendrecv(ibuf1,nsend,
     $           mpi_integer,spart(nswap),0,
     $           tag(nrfirst(nswap)),nrecv,
     $           mpi_integer,rpart(nswap),0,
     $           mpi_comm_world,istatus,ierror)
            call mpi_sendrecv(ibuf2,nsend,
     $           mpi_integer,spart(nswap),0,
     $           type(nrfirst(nswap)),nrecv,
     $           mpi_integer,rpart(nswap),0,
     $           mpi_comm_world,istatus,ierror)
            if (ncharge == 1)
     $           call mpi_sendrecv(buf2,nsend,
     $           mpi_double_precision,spart(nswap),0,
     $           q(nrfirst(nswap)),nrecv,
     $           mpi_double_precision,rpart(nswap),0,
     $           mpi_comm_world,istatus,ierror)
#endif
#ifndef SENDRECV
            call mpi_irecv(x(1,nrfirst(nswap)),3*nrecv,
     $           mpi_double_precision,rpart(nswap),0,
     $           mpi_comm_world,irequest(1),ierror)
            call mpi_irecv(tag(nrfirst(nswap)),nrecv,
     $           mpi_integer,rpart(nswap),0,
     $           mpi_comm_world,irequest(2),ierror)
            call mpi_irecv(type(nrfirst(nswap)),nrecv,
     $           mpi_integer,rpart(nswap),0,
     $           mpi_comm_world,irequest(3),ierror)
            if (ncharge == 1)
     $           call mpi_irecv(q(nrfirst(nswap)),nrecv,
     $           mpi_double_precision,rpart(nswap),0,
     $           mpi_comm_world,irequest(4),ierror)
            call mpi_send(buf1,3*nsend,
     $           mpi_double_precision,spart(nswap),0,
     $           mpi_comm_world,ierror)
            call mpi_send(ibuf1,nsend,
     $           mpi_integer,spart(nswap),0,
     $           mpi_comm_world,ierror)
            call mpi_send(ibuf2,nsend,
     $           mpi_integer,spart(nswap),0,
     $           mpi_comm_world,ierror)
            if (ncharge == 1)
     $           call mpi_send(buf2,nsend,
     $           mpi_double_precision,spart(nswap),0,
     $           mpi_comm_world,ierror)
            call mpi_waitall(nstatus,irequest,istatus,ierror)
#endif

          else

            nrecv = nsend
            if (nrecv.gt.maxrecv) then
              write (6,*) 'Too many ghosts recvd -- ',
     $             'boost extra_ghost:',k,nghost+nrecv,maxghost
              call abort('Ghost atom overflow')
            endif

            m = 0
            n = 0
            do i = nlocal+nghost+1,nlocal+nghost+nrecv
              x(1,i) = buf1(m+1)
              x(2,i) = buf1(m+2)
              x(3,i) = buf1(m+3)
              m = m + 3
              n = n + 1
              tag(i) = ibuf1(n)
              type(i) = ibuf2(n)
            enddo
            if (ncharge == 1) then
              n = 0
              do i = nlocal+nghost+1,nlocal+nghost+nrecv
                n = n + 1
                q(i) = buf2(n)
              enddo
            endif
          endif

          nghost = nghost + nrecv
          nslast(nswap) = npnt
          nrlast(nswap) = nlocal + nghost

        enddo
      enddo
      
c set global ptrs to all owned and ghost atoms
c loop in reverse order so that nearby images take precedence over far images
c  and owned atoms take precedence over images
c this enables valid lookups of bond topology atoms

      do i = nlocal+nghost,1,-1
        localptr(tag(i)) = i
      enddo

      return
      end


c -----------------------------------------------------------------------
c move all info about an atom from location i to location j in local arrays
c used to keep atom list densely packed from 1 to nlocal

      subroutine copyatom(i,j)
      use global
      implicit none

c argument variables

      integer i,j

c local variables

      integer m

c variables active in all cases

      x(1,j) = x(1,i)
      x(2,j) = x(2,i)
      x(3,j) = x(3,i)
      v(1,j) = v(1,i)
      v(2,j) = v(2,i)
      v(3,j) = v(3,i)
      tag(j) = tag(i)
      molecule(j) = molecule(i)
      type(j) = type(i)
      true(j) = true(i)

c variables only active sometimes

      q(j) = q(i)
      fix(j) = fix(i)

      if (nrespa /= 0) then
        f_stretch(1,j) = f_stretch(1,i)
        f_stretch(2,j) = f_stretch(2,i)
        f_stretch(3,j) = f_stretch(3,i)
        f_intra(1,j) = f_intra(1,i)
        f_intra(2,j) = f_intra(2,i)
        f_intra(3,j) = f_intra(3,i)
      endif

c topology variables only active if there are bonds
c else set them to zero

      if (nmolecular == 1) then

        if (nbonds > 0) then
          numbond(j) = numbond(i)
          do m = 1,numbond(j)
            bondatom1(m,j) = bondatom1(m,i)
            bondatom2(m,j) = bondatom2(m,i)
            bondtype(m,j) = bondtype(m,i)
          enddo
        endif

        if (nangles > 0) then
          numangle(j) = numangle(i)
          do m = 1,numangle(j)
            angleatom1(m,j) = angleatom1(m,i)
            angleatom2(m,j) = angleatom2(m,i)
            angleatom3(m,j) = angleatom3(m,i)
            angletype(m,j) = angletype(m,i)
          enddo
        endif

        if (ndihedrals > 0) then
          numdihed(j) = numdihed(i)
          do m = 1,numdihed(j)
            dihedatom1(m,j) = dihedatom1(m,i)
            dihedatom2(m,j) = dihedatom2(m,i)
            dihedatom3(m,j) = dihedatom3(m,i)
            dihedatom4(m,j) = dihedatom4(m,i)
            dihedtype(m,j) = dihedtype(m,i)
          enddo
        endif

        if (nimpropers > 0) then
          numimpro(j) = numimpro(i)
          do m = 1,numimpro(j)
            improatom1(m,j) = improatom1(m,i)
            improatom2(m,j) = improatom2(m,i)
            improatom3(m,j) = improatom3(m,i)
            improatom4(m,j) = improatom4(m,i)
            improtype(m,j) = improtype(m,i)
          enddo
        endif

        num1bond(j) = num1bond(i)
        num2bond(j) = num2bond(i)
        num3bond(j) = num3bond(i)
        do m = 1,num3bond(j)
          specbond(m,j) = specbond(m,i)
        enddo

        if (nshake == 1) then
          shakegroup(j) = shakegroup(i)
          do m = 1,shakegroup(j)
            shakepartner(m,j) = shakepartner(m,i)
          enddo
          do m = 1,shakegroup(j)-1
            shakebondtype(m,j) = shakebondtype(m,i)
          enddo
        endif

      else

        numbond(j) = 0
        numangle(j) = 0
        numdihed(j) = 0
        numimpro(j) = 0
        num1bond(j) = 0
        num2bond(j) = 0
        num3bond(j) = 0

      endif

      return
      end


c -----------------------------------------------------------------------
c pack all info about atom i into buf
c used to keep atom list packed from 1 to nlocal
c store total # of atom datums in 1st location of buf
c also return the total value

      integer function packatom(i,buf)
      use global
      implicit none

c argument variables

      integer i
      real*8 buf(*)

c local variables

      integer j,m

c variables active in all cases

      buf(2) = x(1,i)
      buf(3) = x(2,i)
      buf(4) = x(3,i)
      buf(5) = v(1,i)
      buf(6) = v(2,i)
      buf(7) = v(3,i)
      buf(8) = tag(i)
      buf(9) = molecule(i)
      buf(10) = type(i)
      buf(11) = true(i)
      j = 11

c variables only active sometimes

      if (ncharge == 1) then
        j = j + 1
        buf(j) = q(i)
      endif
      if (nfixes > 0) then
        j = j + 1
        buf(j) = fix(i)
      endif

      if (nrespa == 1) then
        buf(j+1) = f_stretch(1,i)
        buf(j+2) = f_stretch(2,i)
        buf(j+3) = f_stretch(3,i)
        buf(j+4) = f_intra(1,i)
        buf(j+5) = f_intra(2,i)
        buf(j+6) = f_intra(3,i)
        j = j + 6
      endif

c variables only active if there are bonds

      if (nmolecular == 1) then

        if (nbonds > 0) then
          j = j + 1
          buf(j) = numbond(i)
          do m = 1,numbond(i)
            buf(j+1) = bondatom1(m,i)
            buf(j+2) = bondatom2(m,i)
            buf(j+3) = bondtype(m,i)
            j = j + 3
          enddo
        endif

        if (nangles > 0) then
          j = j + 1
          buf(j) = numangle(i)
          do m = 1,numangle(i)
            buf(j+1) = angleatom1(m,i)
            buf(j+2) = angleatom2(m,i)
            buf(j+3) = angleatom3(m,i)
            buf(j+4) = angletype(m,i)
            j = j + 4
          enddo
        endif

        if (ndihedrals > 0) then
          j = j + 1
          buf(j) = numdihed(i)
          do m = 1,numdihed(i)
            buf(j+1) = dihedatom1(m,i)
            buf(j+2) = dihedatom2(m,i)
            buf(j+3) = dihedatom3(m,i)
            buf(j+4) = dihedatom4(m,i)
            buf(j+5) = dihedtype(m,i)
            j = j + 5
          enddo
        endif

        if (nimpropers > 0) then
          j = j + 1
          buf(j) = numimpro(i)
          do m = 1,numimpro(i)
            buf(j+1) = improatom1(m,i)
            buf(j+2) = improatom2(m,i)
            buf(j+3) = improatom3(m,i)
            buf(j+4) = improatom4(m,i)
            buf(j+5) = improtype(m,i)
            j = j + 5
          enddo
        endif

        buf(j+1) = num1bond(i)
        buf(j+2) = num2bond(i)
        buf(j+3) = num3bond(i)
        j = j + 3
        do m = 1,num3bond(i)
          j = j + 1
          buf(j) = specbond(m,i)
        enddo

        if (nshake == 1) then
          j = j + 1
          buf(j) = shakegroup(i)
          do m = 1,shakegroup(i)
            j = j + 1
            buf(j) = shakepartner(m,i)
          enddo
          do m = 1,shakegroup(i)-1
            j = j + 1
            buf(j) = shakebondtype(m,i)
          enddo
        endif

      endif

      buf(1) = j
      packatom = j

      return
      end


c -----------------------------------------------------------------------
c unpack info from buf for atom i

      subroutine unpackatom(i,buf)
      use global
      implicit none

c argument variables

      integer i
      real*8 buf(*)

c local variables

      integer j,m

c variables active in all cases

      x(1,i) = buf(2)
      x(2,i) = buf(3)
      x(3,i) = buf(4)
      v(1,i) = buf(5)
      v(2,i) = buf(6)
      v(3,i) = buf(7)
      tag(i) = nint(buf(8))
      molecule(i) = nint(buf(9))
      type(i) = nint(buf(10))
      true(i) = nint(buf(11))
      j = 11

c variables only active sometimes

      if (ncharge == 1) then
        j = j + 1
        q(i) = buf(j)
      else
        q(i) = 0.0
      endif

      if (nfixes > 0) then
        j = j + 1
        fix(i) = nint(buf(j))
      else
        fix(i) = 0
      endif

      if (nrespa == 1) then
        f_stretch(1,i) = buf(j+1)
        f_stretch(2,i) = buf(j+2)
        f_stretch(3,i) = buf(j+3)
        f_intra(1,i) = buf(j+4)
        f_intra(2,i) = buf(j+5)
        f_intra(3,i) = buf(j+6)
        j = j + 6
      endif

c topology variables only active if there are bonds
c else set them to zero

      if (nmolecular == 1) then

        if (nbonds > 0) then
          j = j + 1
          numbond(i) = nint(buf(j))
          do m = 1,numbond(i)
            bondatom1(m,i) = nint(buf(j+1))
            bondatom2(m,i) = nint(buf(j+2))
            bondtype(m,i) = nint(buf(j+3))
            j = j + 3
          enddo
        endif

        if (nangles > 0) then
          j = j + 1
          numangle(i) = nint(buf(j))
          do m = 1,numangle(i)
            angleatom1(m,i) = nint(buf(j+1))
            angleatom2(m,i) = nint(buf(j+2))
            angleatom3(m,i) = nint(buf(j+3))
            angletype(m,i) = nint(buf(j+4))
            j = j + 4
          enddo
        endif

        if (ndihedrals > 0) then
          j = j + 1
          numdihed(i) = nint(buf(j))
          do m = 1,numdihed(i)
            dihedatom1(m,i) = nint(buf(j+1))
            dihedatom2(m,i) = nint(buf(j+2))
            dihedatom3(m,i) = nint(buf(j+3))
            dihedatom4(m,i) = nint(buf(j+4))
            dihedtype(m,i) = nint(buf(j+5))
            j = j + 5
          enddo
        endif

        if (nimpropers > 0) then
          j = j + 1
          numimpro(i) = nint(buf(j))
          do m = 1,numimpro(i)
            improatom1(m,i) = nint(buf(j+1))
            improatom2(m,i) = nint(buf(j+2))
            improatom3(m,i) = nint(buf(j+3))
            improatom4(m,i) = nint(buf(j+4))
            improtype(m,i) = nint(buf(j+5))
            j = j + 5
          enddo
        endif

        num1bond(i) = nint(buf(j+1))
        num2bond(i) = nint(buf(j+2))
        num3bond(i) = nint(buf(j+3))
        j = j + 3
        do m = 1,num3bond(i)
          j = j + 1
          specbond(m,i) = nint(buf(j))
        enddo

        if (nshake == 1) then
          j = j + 1
          shakegroup(i) = nint(buf(j))
          do m = 1,shakegroup(i)
            j = j + 1
            shakepartner(m,i) = nint(buf(j))
          enddo
          do m = 1,shakegroup(i)-1
            j = j + 1
            shakebondtype(m,i) = nint(buf(j))
          enddo
        endif

      else

        numbond(i) = 0
        numangle(i) = 0
        numdihed(i) = 0
        numimpro(i) = 0
        num1bond(i) = 0
        num2bond(i) = 0
        num3bond(i) = 0

      endif

      return
      end
