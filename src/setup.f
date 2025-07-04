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
c shrink-wrap global simulation box boundaries around atoms
c only called if any boundaries are non-periodic

      subroutine setup_box
      use global
      use mpi
      implicit none

c local variables

      real*8, parameter :: small = 1.0D-4
      real*8, parameter :: big = 1.0D20
      real*8 extent(2,3),extent_all(2,3)
      integer i,ierror
      
      extent(1,1) = big
      extent(2,1) = -big
      extent(1,2) = big
      extent(2,2) = -big
      extent(1,3) = big
      extent(2,3) = -big

      do i = 1,nlocal
        if (x(1,i) < extent(1,1)) extent(1,1) = x(1,i)
        if (x(1,i) > extent(2,1)) extent(2,1) = x(1,i)
        if (x(2,i) < extent(1,2)) extent(1,2) = x(2,i)
        if (x(2,i) > extent(2,2)) extent(2,2) = x(2,i)
        if (x(3,i) < extent(1,3)) extent(1,3) = x(3,i)
        if (x(3,i) > extent(2,3)) extent(2,3) = x(3,i)
      enddo

c global min/max of each extent across all procs (min and max together)

      extent(1,1) = -extent(1,1)
      extent(1,2) = -extent(1,2)
      extent(1,3) = -extent(1,3)

      call mpi_allreduce(extent,extent_all,6,mpi_double_precision,
     $     mpi_max,mpi_comm_world,ierror)

      extent_all(1,1) = -extent_all(1,1)
      extent_all(1,2) = -extent_all(1,2)
      extent_all(1,3) = -extent_all(1,3)

c reset box bounds in each non-periodic dimension
c  adding small -> prevents exchange routine from trying to move an atom
c   exactly on a sub-domain boundary to a non-existent processor
c if using 2-d Ewald/PPPM, do not adjust z box bounds, even if non-periodic

      if (perflagx == 1) then
        box(1,1) = extent_all(1,1) - small
        box(2,1) = extent_all(2,1) + small
      endif
      if (perflagy == 1) then
        box(1,2) = extent_all(1,2) - small
        box(2,2) = extent_all(2,2) + small
      endif
      if (perflagz == 1 .and. slabflag == 0) then
        box(1,3) = extent_all(1,3) - small
        box(2,3) = extent_all(2,3) + small
      endif

c re-define simulation box

      xprd = box(2,1) - box(1,1)
      yprd = box(2,2) - box(1,2)
      zprd = box(2,3) - box(1,3)

      xprd_half = 0.5 * xprd
      yprd_half = 0.5 * yprd
      zprd_half = 0.5 * zprd

      return
      end


c ----------------------------------------------------------------------
c setup spatial-decomposition communication patterns

      subroutine setup_comm
      use global
      use mpi
      implicit none

c local variables

      integer k,kk,nbox
      real*8 blo,bhi,prd(3)
      integer boundary(3)

c border = my local box bounds

      border(1,1) = float(me(1))/pgrid(1) * xprd + box(1,1)
      border(2,1) = float(me(1)+1)/pgrid(1) * xprd + box(1,1)
      border(1,2) = float(me(2))/pgrid(2) * yprd + box(1,2)
      border(2,2) = float(me(2)+1)/pgrid(2) * yprd + box(1,2)
      border(1,3) = float(me(3))/pgrid(3) * zprd + box(1,3)
      border(2,3) = float(me(3)+1)/pgrid(3) * zprd + box(1,3)

c need = # of boxes I need info from in each of 6 directions

      need(1) = cutneigh / (xprd/pgrid(1)) + 1
      need(2) = cutneigh / (yprd/pgrid(2)) + 1
      need(3) = cutneigh / (zprd/pgrid(3)) + 1

c for 2d, don't need to communicate in z

      if (idimension == 2) need(3) = 0

c if non-periodic in a dim, do not communicate further than pgrid-1 away
c this allows very large cutoffs in non-periodic system

      if (perflagx == 1) need(1) = min(need(1),pgrid(1)-1)
      if (perflagy == 1) need(2) = min(need(2),pgrid(2)-1)
      if (perflagz == 1) need(3) = min(need(3),pgrid(3)-1)

c alloc comm memory only if nswap is larger than previous ones
c deallocate first if necessary

      nswap = 2 * (need(1)+need(2)+need(3))

      if (nswap > maxswap) then

        if (allocated(slablo)) then
          deallocate(slablo,slabhi,spart,rpart)
          deallocate(nsfirst,nslast,nrfirst,nrlast)
          deallocate(commflag,commflagall)
        endif

        allocate(slablo(nswap))
        allocate(slabhi(nswap))
        allocate(spart(nswap))
        allocate(rpart(nswap))
        allocate(nsfirst(nswap))
        allocate(nslast(nswap))
        allocate(nrfirst(nswap))
        allocate(nrlast(nswap))
        allocate(commflag(3,nswap))
        allocate(commflagall(nswap))

        maxswap = nswap

      endif

c setup 4 parameters for each exchange: (spart,rpart,slablo,slabhi)
c  spart(nswap) = proc to send to at each swap
c  rpart(nswap) = proc to recv from at each swap
c  slablo/slabhi(nswap) = slab boundaries (in correct dimension) of atoms
c   to send at each swap
c 1st part of if-statement is sending to the west (south,down)
c 2nd part of if-statement is sending to the east (north,up)
c nbox = box (0 - pgrid(k)-1) where atoms to be sent are coming from

      prd(1) = xprd
      prd(2) = yprd
      prd(3) = zprd
      nswap = 0
      
      do k = 1,3
        do kk = 0,2*need(k)-1
          nswap = nswap + 1
          if (mod(kk,2).eq.0) then
            spart(nswap) = mpart(1,k)
            rpart(nswap) = mpart(2,k)
            nbox = me(k) + kk/2
            blo = box(1,k) + float(nbox)/pgrid(k) * prd(k)
            bhi = border(1,k) + cutneigh
            bhi = min(bhi,box(1,k)+float(nbox+1)/pgrid(k)*prd(k))
          else
            spart(nswap) = mpart(2,k)
            rpart(nswap) = mpart(1,k)
            nbox = me(k) - kk/2
            bhi = box(1,k) + float(nbox+1)/pgrid(k) * prd(k)
            blo = border(2,k) - cutneigh
            blo = max(blo,box(1,k)+float(nbox)/pgrid(k)*prd(k))
          endif
          slablo(nswap) = blo
          slabhi(nswap) = bhi
        enddo
      enddo

c set commflag if atoms are being exchanged across a box boundary
c commflag(idim,nswap) =  0 -> not across a boundary
c                      =  1 -> add box-length to position when sending
c                      = -1 -> subtract box-length from position when sending

      nswap = 0
      do k = 1,3
        do kk = 0,2*need(k)-1
          nswap = nswap + 1
          commflagall(nswap) = 0
          commflag(1,nswap) = 0
          commflag(2,nswap) = 0
          commflag(3,nswap) = 0
          if (mod(kk,2).eq.0.and.me(k).eq.0) then
            commflagall(nswap) = 1
            commflag(k,nswap) = 1
          else if (mod(kk,2).eq.1.and.me(k).eq.pgrid(k)-1) then
            commflagall(nswap) = 1
            commflag(k,nswap) = -1
          endif
        enddo
      enddo

c set slablo = slabhi for those swaps that are across non-periodic boundaries
c this effectively prevents any atoms from being swapped
c only for procs who own sub-domains at non-periodic ends of global box

      boundary(1) = perflagx
      boundary(2) = perflagy
      boundary(3) = perflagz

      nswap = 0
      do k = 1,3
        do kk = 0,2*need(k)-1
          nswap = nswap + 1
          if (mod(kk,2).eq.0) then
            if (boundary(k) /= 0 .and. me(k) == 0)
     $           slabhi(nswap) = slablo(nswap)
          else
            if (boundary(k) /= 0 .and. me(k)+1 == pgrid(k))
     $           slabhi(nswap) = slablo(nswap)
          endif
        enddo
      enddo

      return
      end


c ----------------------------------------------------------------------
c setup neighbor binning parameters
c bin numbering is global: 0 = 0.0 to binsize
c                          1 = binsize to 2*binsize
c                          nbin-1 = prd-binsize to binsize
c                          nbin = prd to prd+binsize
c                          -1 = -binsize to 0.0
c coord = lowest and highest values of ghost atom coords I will have
c         add in "small" for round-off safety
c mbinlo = lowest global bin any of my ghost atoms could fall into
c mbinhi = highest global bin any of my ghost atoms could fall into
c mbin = number of bins I need in a dimension
c stencil() = bin offsets in 1-d sense for stencil of surrounding bins

      subroutine setup_neigh
      use global
      implicit none

c local variables

      real*8, parameter :: small = 1.0D-6
      integer i,j,k,nextx,nexty,nextz,mbinxhi,mbinyhi,mbinzhi
      real*8 coord,cutsq
      real*8, external :: bindist

c bins must evenly divide into box size, 
c   becoming larger than cutneigh if necessary
c binsize = 1/2 of cutoff is near optimal

      nbinx = 2.0 * xprd / cutneigh
      nbiny = 2.0 * yprd / cutneigh
      if (idimension == 2) then
        nbinz = 1
      else
        nbinz = 2.0 * zprd / cutneigh
      endif

      if (nbinx.eq.0) nbinx = 1
      if (nbiny.eq.0) nbiny = 1
      if (nbinz.eq.0) nbinz = 1
      
      binsizex = xprd/nbinx
      binsizey = yprd/nbiny
      binsizez = zprd/nbinz

      bininvx = 1.0 / binsizex
      bininvy = 1.0 / binsizey
      bininvz = 1.0 / binsizez

      coord = border(1,1) - cutneigh - small*xprd
      mbinxlo = int((coord-box(1,1))*bininvx)
      if (coord.lt.box(1,1)) mbinxlo = mbinxlo - 1
      coord = border(2,1) + cutneigh + small*xprd
      mbinxhi = int((coord-box(1,1))*bininvx)

      coord = border(1,2) - cutneigh - small*yprd
      mbinylo = int((coord-box(1,2))*bininvy)
      if (coord.lt.box(1,2)) mbinylo = mbinylo - 1
      coord = border(2,2) + cutneigh + small*yprd
      mbinyhi = int((coord-box(1,2))*bininvy)

      coord = border(1,3) - cutneigh - small*zprd
      mbinzlo = int((coord-box(1,3))*bininvz)
      if (coord.lt.box(1,3)) mbinzlo = mbinzlo - 1
      coord = border(2,3) + cutneigh + small*zprd
      mbinzhi = int((coord-box(1,3))*bininvz)

c extend bins by 1 in each direction to insure stencil coverage

      mbinxlo = mbinxlo - 1
      mbinxhi = mbinxhi + 1
      mbinx = mbinxhi - mbinxlo + 1

      mbinylo = mbinylo - 1
      mbinyhi = mbinyhi + 1
      mbiny = mbinyhi - mbinylo + 1

      mbinzlo = mbinzlo - 1
      mbinzhi = mbinzhi + 1
      mbinz = mbinzhi - mbinzlo + 1

c compute bin stencil of all bins whose closest corner to central bin
c   is within neighbor cutoff
c for non-Newton (newton_nonbond = 0),
c   stencil is all surrounding bins including self
c for Newton (newton_nonbond = 1),
c   stencil is bins to the "upper right" of central bin, does NOT include self
c for non-Newton stencil is all surrounding bins including self
c differing i,j,k loops take account of dimension = 2/3 and Newton status
c next(xyz) = how far the stencil could possibly extend

      nextx = cutneigh*bininvx
      if (nextx*binsizex.lt.cutneigh) nextx = nextx + 1
      nexty = cutneigh*bininvy
      if (nexty*binsizey.lt.cutneigh) nexty = nexty + 1
      nextz = cutneigh*bininvz
      if (nextz*binsizez.lt.cutneigh) nextz = nextz + 1

      cutsq = cutneigh*cutneigh

      if (idimension == 2) then
        if (newton_nonbond == 0) then

          nstencil = 0
          do j = -nexty,nexty
            do i = -nextx,nextx
              if (bindist(i,j,0) < cutsq) then
                nstencil = nstencil + 1
                if (nstencil <= 1000) stencil(nstencil) = j*mbinx + i
              endif
            enddo
          enddo

        else

          nstencil = 0
          do j = 0,nexty
            do i = -nextx,nextx
              if (j > 0 .or. (j == 0 .and. i > 0)) then
                if (bindist(i,j,0) < cutsq) then
                  nstencil = nstencil + 1
                  if (nstencil <= 1000) stencil(nstencil) = j*mbinx + i
                endif
              endif
            enddo
          enddo

        endif
      else
        if (newton_nonbond == 0) then

          nstencil = 0
          do k = -nextz,nextz
            do j = -nexty,nexty
              do i = -nextx,nextx
                if (bindist(i,j,k) < cutsq) then
                  nstencil = nstencil + 1
                  if (nstencil <= 1000) stencil(nstencil) =
     $                 k*mbiny*mbinx + j*mbinx + i
                endif
              enddo
            enddo
          enddo

        else

          nstencil = 0
          do k = 0,nextz
            do j = -nexty,nexty
              do i = -nextx,nextx
                if (k > 0 .or. j > 0 .or. (j == 0 .and. i > 0)) then
                  if (bindist(i,j,k) < cutsq) then
                    nstencil = nstencil + 1
                    if (nstencil <= 1000) stencil(nstencil) =
     $                   k*mbiny*mbinx + j*mbinx + i
                  endif
                endif
              enddo
            enddo
          enddo

        endif
      endif

      if (nstencil > 1000) then
        if (node == 0) write (6,*) '# of Stencils =',nstencil
        call err('Too many stencils - use N^2 neighboring')
      endif

c allocate bins only if array is larger than previous one
c deallocate first if necessary

      if (mbinx*mbiny*mbinz > maxbin) then
        if (allocated(binpnt)) deallocate(binpnt)
        allocate(binpnt(mbinx*mbiny*mbinz))
        maxbin = mbinx*mbiny*mbinz
      endif

      return
      end


c ----------------------------------------------------------------------
c compute closest distance between central bin (0,0,0) and bin (i,j,k)

      real*8 function bindist(i,j,k)
      use global
      implicit none

c argument variables

      integer i,j,k

c local variables

      real*8 delx,dely,delz

      if (i.gt.0) then
        delx = (i-1)*binsizex
      else if (i.eq.0) then
        delx = 0.0
      else
        delx = (i+1)*binsizex
      endif
      if (j.gt.0) then
        dely = (j-1)*binsizey
      else if (j.eq.0) then
        dely = 0.0
      else
        dely = (j+1)*binsizey
      endif
      if (k.gt.0) then
        delz = (k-1)*binsizez
      else if (k.eq.0) then
        delz = 0.0
      else
        delz = (k+1)*binsizez
      endif

      bindist = delx*delx + dely*dely + delz*delz

      return
      end
