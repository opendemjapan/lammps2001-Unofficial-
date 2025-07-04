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
c create 1st,2nd,3rd neighbor lists for each atom

      subroutine special_create
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,n,ierror,itag,iflag,ioriginal,itmp,nbuf,nbufall
      integer inum,num1,num2,num3,maxnum1,maxnum2,maxnum3,maxnum
      integer ipnt,inew,inext,iprev,mess_size,mess_tag
      integer loop,irequest,unique
      integer istatus(mpi_status_size)
      integer, allocatable :: onetwo(:,:),onethree(:,:),onefour(:,:)
      integer, allocatable :: buf(:),bufcopy(:)

      call mpi_barrier(mpi_comm_world,ierror)
      if (node == 0)
     $     write (6,*) 'Finding 1-2, 1-3, 1-4 neighbors ...'

c processor IDs for cyclic ring passing

      inext = node + 1
      iprev = node - 1
      if (inext == nprocs) inext = 0
      if (iprev < 0) iprev = nprocs - 1

c ---------------------------------------------------------
c compute num1bond = # of 1-2 neighbors of each of my atoms
c ---------------------------------------------------------

c count bond partners stored by atom itself

      do i = 1,nlocal
        num1bond(i) = numbond(i)
      enddo

c if not newton_bond, then done
c else this only counted 1/2 of all bonds, so count other half

      if (newton_bond == 1) then

c nbuf = largest buffer needed to hold bond partners from any proc

        nbuf = 0
        do i = 1,nlocal
          nbuf = nbuf + numbond(i)
        enddo

        call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)

        allocate(buf(nbufall))
        allocate(bufcopy(nbufall))

c fill buffer with global tags of bond partners of my atoms

        ipnt = 0
        do i = 1,nlocal
          do j = 1,numbond(i)
            buf(ipnt+1) = bondatom2(j,i)
            ipnt = ipnt + 1
          enddo
        enddo

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan tags for atoms I own
c when find one, increment num1bond count for that atom
c do not send/recv buffer to self

        mess_size = ipnt
        mess_tag = 1

        do loop = 1,nprocs
          do i = 1,mess_size
            inew = localptr(buf(i))
            if (inew /= 0 .and. inew <= nlocal)
     $           num1bond(inew) = num1bond(inew) + 1
          enddo
          if (node /= inext) then
            call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $           mpi_comm_world,irequest,ierror)
            call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
            do i = 1,mess_size
              buf(i) = bufcopy(i)
            enddo
          endif
        enddo

c free buffers

        deallocate(buf,bufcopy)

      endif

c -------------------------------------------------
c create list of 1-2 neighbors for each of my atoms
c -------------------------------------------------

c maxnum1 = most 1-2 neighbors of any atom

      itmp = 0
      do i = 1,nlocal
        itmp = max(itmp,num1bond(i))
      enddo

      call mpi_allreduce(itmp,maxnum1,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (node == 0) then
        write (6,*) 'Max possible # of 1-2 neighbors =',maxnum1
        write (1,*) 'Max possible # of 1-2 neighbors =',maxnum1
      endif

c onetwo = list of 1-2 neighbors of each atom

      allocate(onetwo(maxnum1,nlocal))

c zero out num1bond so can use it as accumulating counter

      do i = 1,nlocal
        num1bond(i) = 0
      enddo

c add bond partners stored by atom itself to onetwo list
c if newton_bond is off, bond partner can be either bondatom1 or bondatom2

      do i = 1,nlocal
        do j = 1,numbond(i)
          num1bond(i) = num1bond(i) + 1
          if (tag(i) == bondatom1(j,i))
     $         onetwo(num1bond(i),i) = bondatom2(j,i)
          if (tag(i) == bondatom2(j,i))
     $         onetwo(num1bond(i),i) = bondatom1(j,i)
        enddo
      enddo

c if not newton_bond, then done
c else this only stored 1/2 of all bonds, so store other half

      if (newton_bond == 1) then

c nbuf = largest buffer needed to hold pairs of bond partners from any proc

        nbuf = 0
        do i = 1,nlocal
          nbuf = nbuf + 2*numbond(i)
        enddo

        call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)

        allocate(buf(nbufall))
        allocate(bufcopy(nbufall))

c fill buffer with global tags of both bond partners for all my atoms

        ipnt = 0
        do i = 1,nlocal
          do j = 1,numbond(i)
            buf(ipnt+1) = bondatom1(j,i)
            buf(ipnt+2) = bondatom2(j,i)
            ipnt = ipnt + 2
          enddo
        enddo

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan 2nd-atom tags for atoms I own
c when find one, add 1st-atom tag to onetwo list for 2nd-atom
c do not send/recv buffer to self

        mess_size = ipnt
        mess_tag = 2

        do loop = 1,nprocs
          do i = 2,mess_size,2
            inew = localptr(buf(i))
            if (inew /= 0 .and. inew <= nlocal) then
              num1bond(inew) = num1bond(inew) + 1
              onetwo(num1bond(inew),inew) = buf(i-1)
            endif
          enddo
          if (node /= inext) then
            call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $           mpi_comm_world,irequest,ierror)
            call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
            do i = 1,mess_size
              buf(i) = bufcopy(i)
            enddo
          endif
        enddo

c free buffers

        deallocate(buf,bufcopy)

      endif

c ---------------------------------------------------------
c compute num2bond = # of 1-3 neighbors of each of my atoms
c ---------------------------------------------------------

c nbuf = largest buffer needed to hold info from any proc
c info for each atom = 2 scalars + list of 1-2 neighbors

      nbuf = 0
      do i = 1,nlocal
        nbuf = nbuf + 2 + num1bond(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      allocate(buf(nbufall))
      allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = # of 1-2 neighbors
c  (2) = counter for 1-3 neighbors, initialized to 0
c  (3:N) = list of 1-2 neighbors

      ipnt = 0
      do i = 1,nlocal
        buf(ipnt+1) = num1bond(i)
        buf(ipnt+2) = 0
        ipnt = ipnt + 2
        do j = 1,num1bond(i)
          buf(ipnt+1) = onetwo(j,i)
          ipnt = ipnt + 1
        enddo
      enddo

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan list of 1-2 neighbor tags for atoms I own
c when find one, increment 1-3 count in buf(i+2) by
c   # of 1-2 neighbors of my atom, subtracting one since my list 
c   will contain original atom
c do not send/recv buffer to self

      mess_size = ipnt
      mess_tag = 3

      do loop = 1,nprocs
        i = 0
        do while (i < mess_size)
          num1 = buf(i+1)
          num2 = buf(i+2)
          do j = 1,num1
            inew = localptr(buf(i+2+j))
            if (inew /= 0 .and. inew <= nlocal)
     $           num2 = num2 + num1bond(inew) - 1
          enddo
          buf(i+2) = num2
          i = i + 2 + num1
        enddo
        if (node /= inext) then
          call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
          do i = 1,mess_size
            buf(i) = bufcopy(i)
          enddo
        endif
      enddo

c extract num2bond from buffer that has cycled back to me
c num2bond = # of 1-3 neighbors of each of my atoms

      ipnt = 0
      do i = 1,nlocal
        num2bond(i) = buf(ipnt+2)
        ipnt = ipnt + 2 + num1bond(i)
      enddo

c free buffers

      deallocate(buf,bufcopy)

c -------------------------------------------------
c create list of 1-3 neighbors for each of my atoms
c -------------------------------------------------

c nbuf = largest buffer needed to hold info from any proc
c info for each atom = 4 scalars + list of 1-2 neighs + list of 1-3 neighs

      nbuf = 0
      do i = 1,nlocal
        nbuf = nbuf + 4 + num1bond(i) + num2bond(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      allocate(buf(nbufall))
      allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = global tag of original atom
c  (2) = # of 1-2 neighbors
c  (3) = # of 1-3 neighbors
c  (4) = counter for 1-3 neighbors, initialized to 0
c  (5:N) = list of 1-2 neighbors
c  followed by space for list of 1-3 neighbors

      ipnt = 0
      do i = 1,nlocal
        buf(ipnt+1) = tag(i)
        buf(ipnt+2) = num1bond(i)
        buf(ipnt+3) = num2bond(i)
        buf(ipnt+4) = 0
        ipnt = ipnt + 4
        do j = 1,num1bond(i)
          buf(ipnt+1) = onetwo(j,i)
          ipnt = ipnt + 1
        enddo
        ipnt = ipnt + num2bond(i)
      enddo

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan list of 1-2 neighbor tags for atoms I own
c when find one, add its neighbors to 1-3 list,
c   incrementing the count in buf(i+4), 
c   but excluding the atom whose tag = ioriginal,
c   may include other duplicates but they will be culled later
c do not send/recv buffer to self

      mess_size = ipnt
      mess_tag = 4

      do loop = 1,nprocs
        i = 0
        do while (i < mess_size)
          ioriginal = buf(i+1)
          num1 = buf(i+2)
          num2 = buf(i+3)
          inum = buf(i+4)
          do j = 1,num1
            inew = localptr(buf(i+4+j))
            if (inew /= 0 .and. inew <= nlocal) then
              do k = 1,num1bond(inew)
                if (onetwo(k,inew) /= ioriginal) then
                  buf(i+4+num1+inum+1) = onetwo(k,inew)
                  inum = inum + 1
                endif
              enddo
            endif
          enddo
          buf(i+4) = inum
          i = i + 4 + num1 + num2
        enddo
        if (node /= inext) then
          call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
          do i = 1,mess_size
            buf(i) = bufcopy(i)
          enddo
        endif
      enddo

c maxnum2 = most 1-3 neighbors of any atom

      itmp = 0
      do i = 1,nlocal
        itmp = max(itmp,num2bond(i))
      enddo

      call mpi_allreduce(itmp,maxnum2,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (node == 0) then
        write (6,*) 'Max possible # of 1-3 neighbors =',maxnum2
        write (1,*) 'Max possible # of 1-3 neighbors =',maxnum2
      endif

c onethree = list of 1-3 neighbors of each atom

      allocate(onethree(maxnum2,nlocal))

c fill onethree with buffer values that have been returned to me
c sanity check: accumlated count buf(j+4) should = num2bond for each atom

      j = 0
      do i = 1,nlocal
        if (buf(j+4) /= num2bond(i)) then
          write (6,*) 'Num2bond inconsistent:',
     $         node,i,num2bond(i),buf(j+4)
          call abort('Num2bond error')
        endif
        j = j + 4 + num1bond(i)
        do k = 1,num2bond(i)
          j = j + 1
          onethree(k,i) = buf(j)
        enddo
      enddo

c free buffers

      deallocate(buf,bufcopy)

c ---------------------------------------------------------
c compute num3bond = # of 1-4 neighbors of each of my atoms
c ---------------------------------------------------------

c nbuf = largest buffer needed to hold info from any proc
c info for each atom = 2 scalars + list of 1-3 neighbors

      nbuf = 0
      do i = 1,nlocal
        nbuf = nbuf + 2 + num2bond(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      allocate(buf(nbufall))
      allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = # of 1-3 neighbors
c  (2) = counter for 1-4 neighbors, initialized to 0
c  (3:N) = list of 1-3 neighbors

      ipnt = 0
      do i = 1,nlocal
        buf(ipnt+1) = num2bond(i)
        buf(ipnt+2) = 0
        ipnt = ipnt + 2
        do j = 1,num2bond(i)
          buf(ipnt+1) = onethree(j,i)
          ipnt = ipnt + 1
        enddo
      enddo

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan list of 1-3 neighbor tags for atoms I own
c when find one, increment 1-4 count in buf(i+2) by 
c   # of 1-2 neighbors of my atom,
c   may include duplicates and original atom but they will be culled later
c do not send/recv buffer to self

      mess_size = ipnt
      mess_tag = 5

      do loop = 1,nprocs
        i = 0
        do while (i < mess_size)
          num2 = buf(i+1)
          num3 = buf(i+2)
          do j = 1,num2
            inew = localptr(buf(i+2+j))
            if (inew /= 0 .and. inew <= nlocal)
     $           num3 = num3 + num1bond(inew)
          enddo
          buf(i+2) = num3
          i = i + 2 + num2
        enddo
        if (node /= inext) then
          call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
          do i = 1,mess_size
            buf(i) = bufcopy(i)
          enddo
        endif
      enddo

c extract num3bond from buffer that has cycled back to me
c num3bond = # of 1-4 neighbors of each of my atoms

      ipnt = 0
      do i = 1,nlocal
        num3bond(i) = buf(ipnt+2)
        ipnt = ipnt + 2 + num2bond(i)
      enddo

c free buffers

      deallocate(buf,bufcopy)

c -------------------------------------------------
c create list of 1-4 neighbors for each of my atoms
c -------------------------------------------------

c nbuf = largest buffer needed to hold info from any proc
c info for each atom = 4 scalars + list of 1-2 neighs + 
c    list of 1-3 neighs + list of 1-4 neighs

      nbuf = 0
      do i = 1,nlocal
        nbuf = nbuf + 4 + num1bond(i) + num2bond(i) + num3bond(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      allocate(buf(nbufall))
      allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = # of 1-2 neighbors
c  (2) = # of 1-3 neighbors
c  (3) = # of 1-4 neighbors
c  (4) = counter for 1-4 neighbors, initialized to 0
c  (5:N) = list of 1-2 neighbors
c  (N+1,:) = list of 1-3 neighbors
c  followed by space for list of 1-4 neighbors

      ipnt = 0
      do i = 1,nlocal
        buf(ipnt+1) = num1bond(i)
        buf(ipnt+2) = num2bond(i)
        buf(ipnt+3) = num3bond(i)
        buf(ipnt+4) = 0
        ipnt = ipnt + 4
        do j = 1,num1bond(i)
          buf(ipnt+1) = onetwo(j,i)
          ipnt = ipnt + 1
        enddo
        do j = 1,num2bond(i)
          buf(ipnt+1) = onethree(j,i)
          ipnt = ipnt + 1
        enddo
        ipnt = ipnt + num3bond(i)
      enddo

c cycle buffer around to all nodes in handshaked fashion
c when receive buffer, scan list of 1-3 neighbor tags for atoms I own
c when find one, add its neighbors to 1-4 list,
c   incrementing the count in buf(i+4)
c do not send/recv buffer to self

      mess_size = ipnt
      mess_tag = 6

      do loop = 1,nprocs
        i = 0
        do while (i < mess_size)
          num1 = buf(i+1)
          num2 = buf(i+2)
          num3 = buf(i+3)
          inum = buf(i+4)
          do j = 1,num2
            inew = localptr(buf(i+4+num1+j))
            if (inew /= 0 .and. inew <= nlocal) then
              do k = 1,num1bond(inew)
                buf(i+4+num1+num2+inum+1) = onetwo(k,inew)
                inum = inum + 1
              enddo
            endif
          enddo
          buf(i+4) = inum
          i = i + 4 + num1 + num2 + num3
        enddo
        if (node /= inext) then
          call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
          do i = 1,mess_size
            buf(i) = bufcopy(i)
          enddo
        endif
      enddo

c maxnum3 = most 1-4 neighbors of any atom

      itmp = 0
      do i = 1,nlocal
        itmp = max(itmp,num3bond(i))
      enddo

      call mpi_allreduce(itmp,maxnum3,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (node == 0) then
        write (6,*) 'Max possible # of 1-4 neighbors =',maxnum3
        write (1,*) 'Max possible # of 1-4 neighbors =',maxnum3
      endif

c onefour = list of 1-4 neighbors of each atom

      allocate(onefour(maxnum3,nlocal))

c fill onefour with buffer values that have been returned to me
c sanity check: accumlated count buf(j+4) should = num3bond for each atom

      j = 0
      do i = 1,nlocal
        if (buf(j+4) /= num3bond(i)) then
          write (6,*) 'Num3bond inconsistent:',
     $         node,i,num3bond(i),buf(j+4)
          call abort('Num3bond error')
        endif
        j = j + 4 + num1bond(i) + num2bond(i)
        do k = 1,num3bond(i)
          j = j + 1
          onefour(k,i) = buf(j)
        enddo
      enddo

c free buffers

      deallocate(buf,bufcopy)

c --------------------------
c allocate and fill specbond
c --------------------------

c max of my atoms = max of sum of 1-2, 1-3, 1-4 neighbor lists

      maxspec = 0
      do i = 1,nlocal
        maxspec = max(itmp,num1bond(i)+num2bond(i)+num3bond(i))
      enddo

c wipe-out localptr so it can be used as scratch space

      do i = 1,natoms
        localptr(i) = 0
      enddo

c precompute how many will be left after culling duplicates for each atom
c exclude original atom explicitly
c unique = # of unique special neighbors of this atom
c use localptr as scratch space to check for duplicates

      maxspec = 0

      do i = 1,nlocal

        unique = 0
        localptr(tag(i)) = 1

        do j = 1,num1bond(i)
          itag = onetwo(j,i)
          if (localptr(itag) == 0) then
            unique = unique + 1
            localptr(itag) = 1
          endif
        enddo
        do j = 1,num2bond(i)
          itag = onethree(j,i)
          if (localptr(itag) == 0) then
            unique = unique + 1
            localptr(itag) = 1
          endif
        enddo
        do j = 1,num3bond(i)
          itag = onefour(j,i)
          if (localptr(itag) == 0) then
            unique = unique + 1
            localptr(itag) = 1
          endif
        enddo

        maxspec = max(maxspec,unique)

        localptr(tag(i)) = 0
        do j = 1,num1bond(i)
          localptr(onetwo(j,i)) = 0
        enddo
        do j = 1,num2bond(i)
          localptr(onethree(j,i)) = 0
        enddo
        do j = 1,num3bond(i)
          localptr(onefour(j,i)) = 0
        enddo
        
      enddo

c maxspec = max of maxspec across all procs

      itmp = maxspec
      call mpi_allreduce(itmp,maxspec,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (node == 0) then
        write (6,*)
     $       'Max allocation for special neighbors =',maxspec
        write (1,*)
     $       'Max allocation for special neighbors =',maxspec
      endif

      allocate(specbond(maxspec,maxown))

c fill up specbond with 3 kinds of neighbors for each atom
c use localptr as scratch space to insure no duplicates appear in list
c also insure original atom tag is not stored in neighbor list
c this handles the duplicate entries produced by rings of bonded atoms
c adjust num1bond/num2bond/num3bond values to reflect 
c   removed duplicate entries and to be cummulative counters

      do i = 1,nlocal

        unique = 0
        localptr(tag(i)) = 1

        do j = 1,num1bond(i)
          itag = onetwo(j,i)
          if (localptr(itag) == 0) then
            unique = unique + 1
            specbond(unique,i) = itag
            localptr(itag) = 1
          endif
        enddo
        num1bond(i) = unique

        do j = 1,num2bond(i)
          itag = onethree(j,i)
          if (localptr(itag) == 0) then
            unique = unique + 1
            specbond(unique,i) = itag
            localptr(itag) = 1
          endif
        enddo
        num2bond(i) = unique

        do j = 1,num3bond(i)
          itag = onefour(j,i)
          if (localptr(itag) == 0) then
            unique = unique + 1
            specbond(unique,i) = itag
            localptr(itag) = 1
          endif
        enddo
        num3bond(i) = unique

        localptr(tag(i)) = 0
        do j = 1,num3bond(i)
          localptr(specbond(j,i)) = 0
        enddo

      enddo

c re-create localptr

      do i = 1,nlocal
        localptr(tag(i)) = i
      enddo

c free up local arrays used just in this routine

      deallocate(onetwo,onethree,onefour)

c print out final stats

      itmp = 0
      do i = 1,nlocal
        itmp = max(itmp,num3bond(i))
      enddo
      call mpi_allreduce(itmp,maxnum,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (node == 0) then
        write (6,*) 'Max actual special neighbors =',maxnum
        write (1,*) 'Max actual special neighbors =',maxnum
      endif
      
      return
      end
