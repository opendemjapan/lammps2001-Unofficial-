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
c check each timestep to see if neighbor lists should be constructed
c  return nflag=0 for no, nflag=1 for yes
c  for neightrigger=1, test max distance moved against skin distance

      subroutine check_neighbor(nflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nflag

c local variables

      integer i,iflag,ierror
      real*8 delx,dely,delz,rsq

      iflag = 0

      do i = 1,nlocal
        delx = x(1,i) - xhold(1,i)
        dely = x(2,i) - xhold(2,i)
        delz = x(3,i) - xhold(3,i)
        rsq = delx*delx + dely*dely + delz*delz
        if (rsq > triggersq) iflag = 1
      enddo

      call mpi_allreduce(iflag,nflag,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (nflag == 1 .and. neighago == max(neighfreq,neighdelay))
     $     ndanger = ndanger + 1

      return
      end


c -------------------------------------------------------------------------
c driver for neighbor-list construction routines
c  store current positions (neightrigger=1)
c  call appropriate routine (neighstyle=0,1)
c  store new maximums

      subroutine neighbor
      use global
      implicit none

c local variables

      integer i

      neighago = 0
      numneigh = numneigh + 1

c store current position to check for next reneighboring

      if (neightrigger == 1) then
        do i = 1,nlocal
          xhold(1,i) = x(1,i)
          xhold(2,i) = x(2,i)
          xhold(3,i) = x(3,i)
        enddo
      endif
      
c call appropriate neighbor routine

      if (neighstyle == 0) then
        call neighbor_nsq
      else
        call binsort
        if (newton_nonbond == 1) then
          call neighbor_bin_newton
        else
          call neighbor_bin_nonewton
        endif
      endif
      
c store max amounts

      max_nlocal = max(max_nlocal,nlocal)
      max_nghost = max(max_nghost,nghost)
      if (nswap > 0) max_slist = max(max_slist,nslast(nswap))
      max_neigh = max(max_neigh,nnlast(nlocal))

      return
      end
      
      
c -------------------------------------------------------------------------
c N^2 / 2 search for neighbor pairs
c pair added to list if atoms i and j are both owned and i < j
c no-Newton: pair added if j is ghost (also stored by proc owning j)
c Newton: if j is ghost only me or other proc adds pair,
c         decision based on mod(itag+jtag,2) tests
      
      subroutine neighbor_nsq
      use global
      implicit none

c local variables

      integer npnt,i,itag,itype,num1tmp,num2tmp,num3tmp
      integer j,jj,jtag,jtype,jpnt,iflag,iwhich,ierror
      real*8 xtmp,ytmp,ztmp,delx,dely,delz,rsq

      npnt = 0

      do i = 1,nlocal

        nnfirst(i) = npnt + 1
        itag = tag(i)
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        num1tmp = num1bond(i)
        num2tmp = num2bond(i)
        num3tmp = num3bond(i)

c loop over remaining atoms, owned and ghost

        do j = i+1,nlocal+nghost

          if (newton_nonbond == 1 .and. j > nlocal) then
            jtag = tag(j)
            if (itag > jtag) then
              if (mod(itag+jtag,2) == 0) cycle
            else if (itag < jtag) then
              if (mod(itag+jtag,2) == 1) cycle
            else
              delx = xtmp - x(1,j)
              dely = ytmp - x(2,j)
              delz = ztmp - x(3,j)
              if (x(3,j) < ztmp) then
                cycle
              else if (x(3,j) == ztmp .and. x(2,j) < ytmp) then
                cycle
              else if (x(3,j) == ztmp .and. x(2,j) == ytmp .and.
     $               x(1,j) < xtmp) then
                cycle
              endif
            endif
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq <= cutneighsq(itype,jtype)) then
            jtag = tag(j)

            iflag = 0
            do jpnt = 1,num3tmp
              if (specbond(jpnt,i) == jtag) then
                iflag = 1
                if (jpnt <= num1tmp) then
                  iwhich = 1
                else if (jpnt <= num2tmp) then
                  iwhich = 2
                else
                  iwhich = 3
                endif
                exit
              endif
            enddo

            if (iflag == 0) then
              npnt = npnt + 1
              nlist(npnt) = j
            else if (special(iwhich) == 0.0) then
              continue
            else if (special(iwhich) == 1.0) then
              npnt = npnt + 1
              nlist(npnt) = j
            else
              npnt = npnt + 1
              nlist(npnt) = iwhich*maxatomp1 + j-1
            endif

          endif

        enddo

        nnlast(i) = npnt

        if (npnt > maxneigh) then
          write (6,*) 'Too many neighbors -- ',
     $         'boost extra_neigh:',i,nlocal,npnt,maxneigh
          call abort('Neighbor overflow')
        endif

      enddo
      
      return
      end


c -------------------------------------------------------------------------
c binned neighbor list construction with partial Newton's 3rd law
c each owned atom i checks own bin and surrounding bins in non-Newton stencil
c pair stored once if i,j are both owned and i < j
c pair stored by me if j is ghost (also stored by node owning j)
      
      subroutine neighbor_bin_nonewton
      use global
      implicit none

c local variables

      integer i,ix,iy,iz,ib,npnt,itype,ibnew
      integer ixx,iyy,izz,num1tmp,num2tmp,num3tmp,k,j
      integer jtag,jtype,jpnt,iflag,iwhich,ierror
      integer, external :: coord2bin
      real*8 xtmp,ytmp,ztmp,delx,dely,delz,rsq
      
      npnt = 0

      do i = 1,nlocal

        nnfirst(i) = npnt + 1

        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        ib = coord2bin(xtmp,ytmp,ztmp)
        num1tmp = num1bond(i)
        num2tmp = num2bond(i)
        num3tmp = num3bond(i)

c loop over all atoms in surrounding bins in stencil including self
c only store pair if i < j
c stores own/own pairs only once
c stores all own/ghost pairs on both procs

        do k = 1,nstencil
          ibnew = ib + stencil(k)
          j = binpnt(ibnew)
          do while (j /= 0)
            if (j <= i) then
              j = bin(j)
              cycle
            endif

            jtype = type(j)
            delx = xtmp - x(1,j)
            dely = ytmp - x(2,j)
            delz = ztmp - x(3,j)
            rsq = delx*delx + dely*dely + delz*delz

            if (rsq <= cutneighsq(itype,jtype)) then

              jtag = tag(j)
              iflag = 0
              do jpnt = 1,num3tmp
                if (specbond(jpnt,i) == jtag) then
                  iflag = 1
                  if (jpnt <= num1tmp) then
                    iwhich = 1
                  else if (jpnt <= num2tmp) then
                    iwhich = 2
                  else
                    iwhich = 3
                  endif
                  exit
                endif
              enddo
              
              if (iflag == 0) then
                npnt = npnt + 1
                nlist(npnt) = j
              else if (special(iwhich) == 0.0) then
                continue
              else if (special(iwhich) == 1.0) then
                npnt = npnt + 1
                nlist(npnt) = j
              else
                npnt = npnt + 1
                nlist(npnt) = iwhich*maxatomp1 + j-1
              endif

            endif
            j = bin(j)
          enddo
        enddo
        
        nnlast(i) = npnt
        
        if (npnt > maxneigh) then
          write (6,*) 'Too many neighbors -- ',
     $         'boost extra_neigh:',i,nlocal,npnt,maxneigh
          call abort('Neighbor overflow')
        endif
        
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c binned neighbor list construction with full Newton's 3rd law
c every pair stored exactly once by some processor
c each owned atom i checks its own bin and other bins in Newton stencil
      
      subroutine neighbor_bin_newton
      use global
      implicit none

c local variables

      integer i,ix,iy,iz,ib,npnt,itype,ibnew
      integer ixx,iyy,izz,num1tmp,num2tmp,num3tmp,k,j
      integer jtag,jtype,jpnt,iflag,iwhich,ierror
      integer, external :: coord2bin
      real*8 xtmp,ytmp,ztmp,delx,dely,delz,rsq
      
      npnt = 0

      do i = 1,nlocal

        nnfirst(i) = npnt + 1

        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        ib = coord2bin(xtmp,ytmp,ztmp)
        num1tmp = num1bond(i)
        num2tmp = num2bond(i)
        num3tmp = num3bond(i)

c loop over rest of atoms in i's bin, ghosts are at end of linked list
c if j is owned atom, store it, since j is beyond i in linked list
c if j is ghost, only store if j coords are "above and to the right" of i

        j = bin(i)
        do while (j /= 0)
          if (j > nlocal) then
            if (x(3,j) < ztmp) then
              j = bin(j)
              cycle
            else if (x(3,j) == ztmp .and. x(2,j) < ytmp) then
              j = bin(j)
              cycle
            else if (x(3,j) == ztmp .and. x(2,j) == ytmp .and.
     $             x(1,j) < xtmp) then
              j = bin(j)
              cycle
            endif
          endif

          jtype = type(j)
          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz

          if (rsq <= cutneighsq(itype,jtype)) then

            jtag = tag(j)
            iflag = 0
            do jpnt = 1,num3tmp
              if (specbond(jpnt,i) == jtag) then
                iflag = 1
                if (jpnt <= num1tmp) then
                  iwhich = 1
                else if (jpnt <= num2tmp) then
                  iwhich = 2
                else
                  iwhich = 3
                endif
                exit
              endif
            enddo
            
            if (iflag == 0) then
              npnt = npnt + 1
              nlist(npnt) = j
            else if (special(iwhich) == 0.0) then
              continue
            else if (special(iwhich) == 1.0) then
              npnt = npnt + 1
              nlist(npnt) = j
            else
              npnt = npnt + 1
              nlist(npnt) = iwhich*maxatomp1 + j-1
            endif
            
          endif
          j = bin(j)
        enddo

c loop over all atoms in other bins in stencil, store every pair

        do k = 1,nstencil
          ibnew = ib + stencil(k)
          j = binpnt(ibnew)
          do while (j /= 0)
            jtype = type(j)
            delx = xtmp - x(1,j)
            dely = ytmp - x(2,j)
            delz = ztmp - x(3,j)
            rsq = delx*delx + dely*dely + delz*delz

            if (rsq <= cutneighsq(itype,jtype)) then

              jtag = tag(j)
              iflag = 0
              do jpnt = 1,num3tmp
                if (specbond(jpnt,i) == jtag) then
                  iflag = 1
                  if (jpnt <= num1tmp) then
                    iwhich = 1
                  else if (jpnt <= num2tmp) then
                    iwhich = 2
                  else
                    iwhich = 3
                  endif
                  exit
                endif
              enddo
              
              if (iflag == 0) then
                npnt = npnt + 1
                nlist(npnt) = j
              else if (special(iwhich) == 0.0) then
                continue
              else if (special(iwhich) == 1.0) then
                npnt = npnt + 1
                nlist(npnt) = j
              else
                npnt = npnt + 1
                nlist(npnt) = iwhich*maxatomp1 + j-1
              endif
              
            endif
            j = bin(j)
          enddo
        enddo
        
        nnlast(i) = npnt
        
        if (npnt > maxneigh) then
          write (6,*) 'Too many neighbors -- ',
     $         'boost extra_neigh:',i,nlocal,npnt,maxneigh
          call abort('Neighbor overflow')
        endif

      enddo
      
      return
      end


c -----------------------------------------------------------------------
c bin own and ghost atoms

      subroutine binsort
      use global
      implicit none

      integer i,ib
      integer, external :: coord2bin

      do i = 1,mbinx*mbiny*mbinz
        binpnt(i) = 0
      enddo
      
c bin ghost atoms 1st, so will be at end of linked list

      do i = nlocal+1,nlocal+nghost
        ib = coord2bin(x(1,i),x(2,i),x(3,i))
        bin(i) = binpnt(ib)
        binpnt(ib) = i
      enddo

c bin own atoms

      do i = 1,nlocal
        ib = coord2bin(x(1,i),x(2,i),x(3,i))
        bin(i) = binpnt(ib)
        binpnt(ib) = i
      enddo

      return
      end

      
c -------------------------------------------------------------------------
c construction of local bond list
c  for each bond store 2 local atoms in bond and ptr to bond parameters 
c  j12 = 0 check is to see if atom in bond is not in own or other list
c  if so is an error
c do not add bond to neighbor list if SHAKE will constrain it

      subroutine neighbor_bond
      use global
      implicit none

c local variables

      integer i,k,j1,j2,ierror

      nbondlocal = 0
      do i = 1,nlocal
        do k = 1,numbond(i)
          if (nshake == 1 .and.
     $         shakeablebond(bondtype(k,i)) == 1) cycle
          j1 = localptr(bondatom1(k,i))
          j2 = localptr(bondatom2(k,i))
          if (j1 == 0 .or. j2 == 0) then
            write (6,*) 'Bond partner missing:',node,
     $           bondatom1(k,i),bondatom2(k,i),numbond(i),ntimestep
            call abort('Bond error')
          endif
          if (newton_bond == 1 .or. (i <= j1 .and. i <= j2)) then
            nbondlocal = nbondlocal + 1
            bondlist(1,nbondlocal) = j1
            bondlist(2,nbondlocal) = j2
            bondlist(3,nbondlocal) = bondtype(k,i)
          endif
        enddo
      enddo

      max_bond = max(max_bond,nbondlocal)

      return
      end


c -------------------------------------------------------------------------
c construction of local angle list
c  for each angle store 3 local atoms in angle and ptr to angle parameters 
c  j123 = 0 check is to see if atom in angle is not in own or other list
c  if so is an error
c do not add angle to neighbor list if SHAKE will constrain it

      subroutine neighbor_angle
      use global
      implicit none

c local variables

      integer i,k,j1,j2,j3,ierror

      nanglelocal = 0
      do i = 1,nlocal
        do k = 1,numangle(i)
          if (nshake == 1 .and. shakeableangle == angletype(k,i)) cycle
          j1 = localptr(angleatom1(k,i))
          j2 = localptr(angleatom2(k,i))
          j3 = localptr(angleatom3(k,i))
          if (j1 == 0 .or. j2 == 0 .or. j3 == 0) then
            write (6,*) 'Angle partner missing:',node,
     $           angleatom1(k,i),angleatom2(k,i),angleatom3(k,i),
     $           numangle(i),ntimestep
            call abort('Angle error')
          endif
          if (newton_bond == 1 .or. 
     $         (i <= j1 .and. i <= j2 .and. i <= j3)) then
            nanglelocal = nanglelocal + 1
            anglelist(1,nanglelocal) = j1
            anglelist(2,nanglelocal) = j2
            anglelist(3,nanglelocal) = j3
            anglelist(4,nanglelocal) = angletype(k,i)
          endif
        enddo
      enddo

      max_angle = max(max_angle,nanglelocal)

      return
      end


c -------------------------------------------------------------------------
c construction of local dihedral list
c  for each dihed store 4 local atoms in dihed and ptr to dihed parameters 
c  j1234 = 0 check is to see if atom in dihed is not in own or other list
c  if so is an error

      subroutine neighbor_dihedral
      use global
      implicit none

c local variables

      integer i,k,j1,j2,j3,j4,ierror

      ndihedlocal = 0
      do i = 1,nlocal
        do k = 1,numdihed(i)
          j1 = localptr(dihedatom1(k,i))
          j2 = localptr(dihedatom2(k,i))
          j3 = localptr(dihedatom3(k,i))
          j4 = localptr(dihedatom4(k,i))
          if (j1 == 0 .or. j2 == 0 .or. j3 == 0 .or. j4 == 0) then
            write (6,*) 'Dihedral partner missing:',node,
     $           dihedatom1(k,i),dihedatom2(k,i),dihedatom3(k,i),
     $           dihedatom4(k,i),numdihed(i),ntimestep
            call abort('Dihedral error')
          endif
          if (newton_bond == 1 .or. 
     $         (i <= j1 .and. i <= j2 .and.
     $         i <= j3 .and. i <= j4)) then
            ndihedlocal = ndihedlocal + 1
            dihedlist(1,ndihedlocal) = j1
            dihedlist(2,ndihedlocal) = j2
            dihedlist(3,ndihedlocal) = j3
            dihedlist(4,ndihedlocal) = j4
            dihedlist(5,ndihedlocal) = dihedtype(k,i)
          endif
        enddo
      enddo

      max_dihed = max(max_dihed,ndihedlocal)

      return
      end


c -------------------------------------------------------------------------
c construction of local improper list
c  for each impro store 4 local atoms in impro and ptr to impro parameters 
c  j1234 = 0 check is to see if atom in impro is not in own or other list
c  if so is an error

      subroutine neighbor_improper
      use global
      implicit none

c local variables

      integer i,k,j1,j2,j3,j4,ierror

      nimprolocal = 0
      do i = 1,nlocal
        do k = 1,numimpro(i)
          j1 = localptr(improatom1(k,i))
          j2 = localptr(improatom2(k,i))
          j3 = localptr(improatom3(k,i))
          j4 = localptr(improatom4(k,i))
          if (j1 == 0 .or. j2 == 0 .or. j3 == 0 .or. j4 == 0) then
            write (6,*) 'Improper partner missing:',node,
     $           improatom1(k,i),improatom2(k,i),improatom3(k,i),
     $           improatom4(k,i),numimpro(i),ntimestep
            call abort('Improper error')
          endif
          if (newton_bond == 1 .or. 
     $         (i <= j1 .and. i <= j2 .and.
     $         i <= j3 .and. i <= j4)) then
            nimprolocal = nimprolocal + 1
            improlist(1,nimprolocal) = j1
            improlist(2,nimprolocal) = j2
            improlist(3,nimprolocal) = j3
            improlist(4,nimprolocal) = j4
            improlist(5,nimprolocal) = improtype(k,i)
          endif
        enddo
      enddo

      max_impro = max(max_impro,nimprolocal)

      return
      end


c -------------------------------------------------------------------------
c construction of local SHAKE list
c store entry in local list exactly once no matter how many atoms
c  in SHAKE group are on this proc
c stored by multiple procs if atoms in SHAKE group are on different procs
c for clusters of size 3, shakesize can be set to +3 or -3
c  +3 means normal bond-only SHAKE, -3 means bond/angle SHAKE

      subroutine neighbor_shake
      use global
      implicit none

c local variables

      integer i,j1,j2,j3,j4

      nshakelocal = 0
      do i = 1,nlocal
        if (shakegroup(i) > 0) then
          if (shakegroup(i) == 2) then
            j1 = localptr(shakepartner(1,i))
            j2 = localptr(shakepartner(2,i))
            if (j1 == 0 .or. j2 == 0) then
              write (6,*) 'SHAKE partner missing:',node,
     $             shakepartner(1,i),shakepartner(2,i),ntimestep
              call abort('SHAKE error')
            endif
            if (i <= j1 .and. i <= j2) then
              nshakelocal = nshakelocal + 1
              shakesize(nshakelocal) = 2
              shakeatom(1,nshakelocal) = j1
              shakeatom(2,nshakelocal) = j2
              shakebondlen(1,nshakelocal) =
     $             bondcoeff(shakewhichbondcoeff,shakebondtype(1,i))
            endif
          else if (shakegroup(i) == 3) then
            j1 = localptr(shakepartner(1,i))
            j2 = localptr(shakepartner(2,i))
            j3 = localptr(shakepartner(3,i))
            if (j1 == 0 .or. j2 == 0 .or. j3 == 0) then
              write (6,*) 'SHAKE partner missing:',node,
     $             shakepartner(1,i),shakepartner(2,i),
     $             shakepartner(3,i),
     $             ntimestep
              call abort('SHAKE error')
            endif
            if (i <= j1 .and. i <= j2 .and. i <= j3) then
              nshakelocal = nshakelocal + 1
              shakesize(nshakelocal) = 3
              shakeatom(1,nshakelocal) = j1
              shakeatom(2,nshakelocal) = j2
              shakeatom(3,nshakelocal) = j3
              shakebondlen(1,nshakelocal) =
     $             bondcoeff(shakewhichbondcoeff,shakebondtype(1,i))
              shakebondlen(2,nshakelocal) =
     $             bondcoeff(shakewhichbondcoeff,shakebondtype(2,i))
              if (shakeableangle > 0) then
                if (shakebondtype(1,i) == shakeableanglebond .and.
     $               shakebondtype(2,i) == shakeableanglebond)
     $               shakesize(nshakelocal) = -3
              endif
            endif
          else
            j1 = localptr(shakepartner(1,i))
            j2 = localptr(shakepartner(2,i))
            j3 = localptr(shakepartner(3,i))
            j4 = localptr(shakepartner(4,i))
            if (j1 == 0 .or. j2 == 0 .or. j3 == 0 .or. j4 == 0) then
              write (6,*) 'SHAKE partner missing:',node,tag(i),
     $             shakepartner(1,i),shakepartner(2,i),
     $             shakepartner(3,i),shakepartner(4,i),
     $             ntimestep
              call abort('SHAKE error')
            endif
            if (i <= j1 .and. i <= j2 .and. i <= j3 .and. i <= j4) then
              nshakelocal = nshakelocal + 1
              shakesize(nshakelocal) = 4
              shakeatom(1,nshakelocal) = j1
              shakeatom(2,nshakelocal) = j2
              shakeatom(3,nshakelocal) = j3
              shakeatom(4,nshakelocal) = j4
              shakebondlen(1,nshakelocal) =
     $             bondcoeff(shakewhichbondcoeff,shakebondtype(1,i))
              shakebondlen(2,nshakelocal) =
     $             bondcoeff(shakewhichbondcoeff,shakebondtype(2,i))
              shakebondlen(3,nshakelocal) =
     $             bondcoeff(shakewhichbondcoeff,shakebondtype(3,i))
            endif
          endif
        endif
      enddo

      return
      end
