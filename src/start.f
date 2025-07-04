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
c beginning of timestep

      subroutine start
      use global
      use mpi
      implicit none

c local variables

      integer i,iflag,jflag,ierror,itmp,iwhich,j,idiag,vflag
      real*8 cut,ratio,r6inv,t,tsq
      integer, allocatable :: scratch1(:),scratch2(:)

      if (node == 0) write (6,*) 'Setting up run ...'

      if (readflag == 0)
     $     call err('Data file must be read before run')

c other errors/warnings

      if ((neighfreq < neighdelay .and. 
     $     mod(neighdelay,neighfreq) /= 0) .or.
     $     (neighdelay > 0 .and. neighdelay < neighfreq)) then
        if (node == 0) write (6,*) 'WARNING: Mismatch between ',
     $       ' neighbor frequency and delay factors'
      endif

      if (coulstyle == 3 .or. coulstyle == 4) then
        if (idimension == 2)
     $     call err('Ewald/PPPM requires 3-d systems')
        if (nonperiodic .and. slabflag == 0)
     $     call err('3-d Ewald/PPPM requires all periodic BC')
      endif

      if (slabflag == 1 .and. (coulstyle /= 3 .and. coulstyle /= 4))
     $     call err('Slab setting requires Ewald or PPPM')

      if (slabflag == 1 .and.
     $     (perflagx == 1 .or. perflagy == 1 .or. perflagz == 0))
     $     call err(
     $     '2-d slab Ewald/PPPM requires periodic xy, non-periodic z')

      if (slabflag == 1 .and. abs(qsum) > 0.01) call err(
     $     'Cannot do 2-d slab Ewald/PPPM with non-neutral system')

      if (pressstyle > 0) then
        if (p_freq(1) /= 0.0 .and. perflagx == 1)
     $       call err('Pressure control requires periodic BC')
        if (p_freq(2) /= 0.0 .and. perflagy == 1)
     $       call err('Pressure control requires periodic BC')
        if (p_freq(3) /= 0.0 .and. perflagz == 1)
     $       call err('Pressure control requires periodic BC')
      endif

      if (volstyle > 0) then
        if (voldimx == 1 .and. perflagx == 1)
     $       call err('Volume control requires periodic BC')
        if (voldimy == 1 .and. perflagy == 1)
     $       call err('Volume control requires periodic BC')
        if (voldimz == 1 .and. perflagz == 1)
     $       call err('Volume control requires periodic BC')
      endif

      if (volstyle > 0 .and. pressstyle > 0) then
        if (voldimx == 1 .and. p_freq(1) /= 0.0) call err(
     $       'Volume control incompatible with pressure control')
        if (voldimy == 1 .and. p_freq(2) /= 0.0) call err(
     $       'Volume control incompatible with pressure control')
        if (voldimz == 1 .and. p_freq(3) /= 0.0) call err(
     $       'Volume control incompatible with pressure control')
      endif

      if ((nonstyle == 3 .or. nonstyle == 4) .and. coulstyle /= 0)
     $     call err('Mismatch between nonbond and coulomb styles')
      if (nonstyle == 5 .and. coulstyle == 2)
     $     call err('Mismatch between nonbond and coulomb styles')
      if (nonstyle == 6 .and.
     $     (coulstyle /= 3 .and. coulstyle /= 4 .and. coulstyle /= 5))
     $     call err('Mismatch between nonbond and coulomb styles')
      if (coulstyle == 5 .and. nonstyle /= 6)
     $     call err('Mismatch between nonbond and coulomb styles')
      if (coulstyle == 6 .and. (nonstyle /= 0 .and. nonstyle /= 1))
     $     call err('Mismatch between nonbond and coulomb styles')

      if (nshake == 1 .and. nrespa == 1)
     $     call err('Cannot combine SHAKE with rRESPA')

      if (dihedstyle == 4 .and. special(3) /= 0) call err(
     $     '1-4 special bond must be zero for CHARMM dihedral')

c set molecular and charge flags

      nmolecular = 0
      if (nbonds > 0 .or. nangles > 0 .or.
     $     ndihedrals > 0 .or. nimpropers > 0) nmolecular = 1

      iflag = 0
      do i = 1,nlocal
        if (q(i) /= 0.0) iflag = 1
      enddo
      call mpi_allreduce(iflag,ncharge,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      
c check for using Coulomb style with no charges

      if (coulstyle > 0) then
        if (ncharge == 0) then
          write (6,*) 'WARNING: Coulomb style set ',
     $         'but all charges = 0.0'
          write (6,*) '  Setting coulomb style to none ...'
          coulstyle = 0
          cutcoul = 0.0
        endif
      endif

c check for using rRESPA or SHAKE with no bonds

      if (nrespa == 1 .and. nmolecular == 0)
     $     call err('Using rREPSA with no molecular bonds')

      if (nshake == 1 .and. nmolecular == 0)
     $     call err('Using SHAKE with no molecular bonds')

c check for skipped constraint number

      if (nfixes > 0) then
        do i = 1,nfixes
          if (fixstyle(i) == 0) call err('Missing fix')
        enddo
      endif

c check no atom has been assigned to missing constraint
c each constraint must have at least one atom assigned to it
c   except SHAKE, which is fixstyle = 9

      if (nfixes > 0) then
        allocate(scratch1(nfixes))
        allocate(scratch2(nfixes))
        do i = 1,nfixes
          scratch1(i) = 0
        enddo
      endif

      iflag = 0
      do i = 1,nlocal
        itmp = fix(i)
        do while (itmp > 0)
          iwhich = mod(itmp,maxfix)
          if (iwhich > nfixes) then
            iflag = 1
          else
            scratch1(iwhich) = 1
          endif
          itmp = itmp/maxfix
        enddo
      enddo
        
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag > 0) call err('Atom assigned to missing fix')

      if (nfixes > 0) then
        call mpi_allreduce(scratch1,scratch2,nfixes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nfixes
          if (scratch2(i) == 0 .and. fixstyle(i) /= 9)
     $         call err('No atom assigned to a fix')
        enddo
        deallocate(scratch1,scratch2)
      endif

c check that no atom is assigned to same constraint more than once

      if (nfixes > 0) then

        allocate(scratch1(nfixes))

        iflag = 0
        do i = 1,nlocal
          do j = 1,nfixes
            scratch1(j) = 0
          enddo
          itmp = fix(i)
          do while (itmp > 0)
            iwhich = mod(itmp,maxfix)
            scratch1(iwhich) = scratch1(iwhich) + 1
            itmp = itmp/maxfix
          enddo
          do j = 1,nfixes
            if (scratch1(j) > 1) iflag = 1
          enddo
        enddo

        deallocate(scratch1)
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0) call err(
     $       'An atom is assigned more than once to same fix')
        
      endif

c initialize constraint pointers and counters

      if (nfixes > 0) call fix_setup

c print out constraint statistics, except for SHAKE

      if (node == 0) then
        do i = 1,nfixes
          if (fixstyle(i) /= 9) then
            write (6,*) fixcount(i),' atoms assigned to fix',i
            write (1,*) fixcount(i),' atoms assigned to fix',i
          endif
        enddo
      endif

c do not allow thermostatting fixes in conjunction with global temp control

      if (tempstyle > 0) then
        iflag = 0
        do i = 1,nfixes
          if (fixstyle(i) == 4 .or. fixstyle(i) == 5 .or. 
     $         fixstyle(i) == 6) iflag = 1
        enddo
        if (iflag > 0) call err(
     $       'Temperature control incompatible with temperature fix')
      endif

c check consistency of atom types and nonbond type pairs

      allocate(scratch1(ntypes))
      allocate(scratch2(ntypes))

      do i = 1,ntypes
        scratch1(i) = 0
      enddo

      iflag = 0
      do i = 1,nlocal
        if (type(i) <= 0 .or. type(i) > ntypes) then
          iflag = 1
        else
          scratch1(type(i)) = 1
        endif
      enddo
      
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag > 0) call err('Illegal atom type')

      call mpi_allreduce(scratch1,scratch2,ntypes,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      do i = 1,ntypes
        if (scratch2(i) == 0 .and. node == 0)
     $       write (6,*) 'WARNING: Non-existent atom type',i
      enddo

      deallocate(scratch1,scratch2)

      iflag = 0
      do i = 1,ntypes
        do j = i,ntypes
          if (nontypeflag(i,j) /= nonstyle) iflag = 1
        enddo
      enddo
      
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag > 0 .and. nonstyle > 0)
     $     call err('Undefined nonbond type pair')

c check consistency of bond types

      if (nbonds > 0) then

        allocate(scratch1(nbondtypes))
        allocate(scratch2(nbondtypes))

        do i = 1,nbondtypes
          scratch1(i) = 0
        enddo

        iflag = 0
        do i = 1,nlocal
          do j = 1,numbond(i)
            if (bondtype(j,i) <= 0 .or.
     $           bondtype(j,i) > nbondtypes) then
              iflag = 1
            else
              scratch1(bondtype(j,i)) = 1
            endif
          enddo
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0) call err('Illegal bond type')

        call mpi_allreduce(scratch1,scratch2,nbondtypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nbondtypes
          if (scratch2(i) == 0 .and. node == 0)
     $         write (6,*) 'WARNING: Non-existent bond type',i
        enddo

        deallocate(scratch1,scratch2)

        iflag = 0
        do i = 1,nbondtypes
          if (bondtypeflag(i) /= bondstyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0 .and. bondstyle > 0)
     $       call err('Undefined bond type')

      endif

c check consistency of angle types

      if (nangles > 0) then

        allocate(scratch1(nangletypes))
        allocate(scratch2(nangletypes))

        do i = 1,nangletypes
          scratch1(i) = 0
        enddo

        iflag = 0
        do i = 1,nlocal
          do j = 1,numangle(i)
            if (angletype(j,i) <= 0 .or. 
     $           angletype(j,i) > nangletypes) then
              iflag = 1
            else
              scratch1(angletype(j,i)) = 1
            endif
          enddo
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0) call err('Illegal angle type')

        call mpi_allreduce(scratch1,scratch2,nangletypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nangletypes
          if (scratch2(i) == 0 .and. node == 0)
     $         write (6,*) 'WARNING: Non-existent angle type',i
        enddo

        deallocate(scratch1,scratch2)

        iflag = 0
        do i = 1,nangletypes
          if (angletypeflag(i) /= anglestyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0 .and. anglestyle > 0)
     $       call err('Undefined angle type')

      endif

c check consistency of dihedral types

      if (ndihedrals > 0) then

        allocate(scratch1(ndihedtypes))
        allocate(scratch2(ndihedtypes))

        do i = 1,ndihedtypes
          scratch1(i) = 0
        enddo

        iflag = 0
        do i = 1,nlocal
          do j = 1,numdihed(i)
            if (dihedtype(j,i) <= 0 .or. 
     $           dihedtype(j,i) > ndihedtypes) then
              iflag = 1
            else
              scratch1(dihedtype(j,i)) = 1
            endif
          enddo
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0) call err('Illegal dihed type')

        call mpi_allreduce(scratch1,scratch2,ndihedtypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,ndihedtypes
          if (scratch2(i) == 0 .and. node == 0)
     $         write (6,*) 'WARNING: Non-existent dihedral type',i
        enddo

        deallocate(scratch1,scratch2)

        iflag = 0
        do i = 1,ndihedtypes
          if (dihedtypeflag(i) /= dihedstyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0 .and. dihedstyle > 0)
     $       call err('Undefined dihedral type')

      endif

c check consistency of improper types

      if (nimpropers > 0) then

        allocate(scratch1(nimprotypes))
        allocate(scratch2(nimprotypes))

        do i = 1,nimprotypes
          scratch1(i) = 0
        enddo

        iflag = 0
        do i = 1,nlocal
          do j = 1,numimpro(i)
            if (improtype(j,i) <= 0 .or. 
     $           improtype(j,i) > nimprotypes) then
              iflag = 1
            else
              scratch1(improtype(j,i)) = 1
            endif
          enddo
        enddo
        
        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0) call err('Illegal impro type')

        call mpi_allreduce(scratch1,scratch2,nimprotypes,mpi_integer,
     $       mpi_sum,mpi_comm_world,ierror)
        do i = 1,nimprotypes
          if (scratch2(i) == 0 .and. node == 0)
     $         write (6,*) 'WARNING: Non-existent improper type',i
        enddo

        deallocate(scratch1,scratch2)

        iflag = 0
        do i = 1,nimprotypes
          if (improtypeflag(i) /= improstyle) iflag = 1
        enddo

        call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        if (jflag > 0 .and. improstyle > 0)
     $       call err('Undefined improper type')

      endif

c allocate arrays for nonbond interactions

      call nonbond_memory

c pre-compute coeffs and cutoff distances for all nonbond/Coulombic styles
c  assumes noncoeffs are set for all type pairs with i <= j

      cutcoulsq = cutcoul*cutcoul
      if (coulstyle == 2) cutcoulintsq = cutcoulint*cutcoulint
      if (coulstyle == 5) cutcoulintsq = cutcoulint*cutcoulint
      triggersq = (skin/2.0)*(skin/2.0)
      cutforce = 0.0

      do i = 1,ntypes
        do j = i,ntypes

          if (nonstyle == 0) then
            cutljsq(i,j) = 0.0
            cut = cutcoul
          else if (nonstyle == 1) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            cutljsq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            if (offsetflag == 0 .or. noncoeff3(i,j) == 0.0) then
              offset(i,j) = 0.0
            else
              ratio = noncoeff2(i,j)/noncoeff3(i,j)
              offset(i,j) = 4.0*noncoeff1(i,j)*(ratio**12 - ratio**6)
            endif
            cut = max(cutcoul,noncoeff3(i,j))
          else if (nonstyle == 2) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            cutljinner(i,j) = noncoeff3(i,j)
            cutljinnersq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            cutljsq(i,j) = noncoeff4(i,j)*noncoeff4(i,j)
            if (noncoeff3(i,j) /= 0.0) then
              r6inv = 1.0/noncoeff3(i,j)**6
              t = noncoeff4(i,j) - noncoeff3(i,j)
              tsq = t*t
              ratio = noncoeff2(i,j)/noncoeff3(i,j)
              ljsw0(i,j) = 4.0*noncoeff1(i,j)*(ratio**12 - ratio**6)
              ljsw1(i,j) = r6inv*(lj1(i,j)*r6inv-lj2(i,j)) /
     $             noncoeff3(i,j)
              ljsw2(i,j) = -r6inv *
     $             (13.0*lj1(i,j)*r6inv - 7.0*lj2(i,j)) /
     $             cutljinnersq(i,j)
              ljsw3(i,j) = -(3.0/tsq) *
     $             (ljsw1(i,j) + 2.0/3.0*ljsw2(i,j)*t)
              ljsw4(i,j) = -1.0/(3.0*tsq) *
     $             (ljsw2(i,j) + 2.0*ljsw3(i,j)*t)
              offset(i,j) = ljsw0(i,j) - ljsw1(i,j)*t -
     $             ljsw2(i,j)*tsq/2.0 - ljsw3(i,j)*tsq*t/3.0 -
     $             ljsw4(i,j)*tsq*tsq/4.0
            else
              ljsw0(i,j) = 0.0
              ljsw1(i,j) = 0.0
              ljsw2(i,j) = 0.0
              ljsw3(i,j) = 0.0
              ljsw4(i,j) = 0.0
              offset(i,j) = 0.0
            endif
            cut = max(cutcoul,noncoeff4(i,j))
          else if (nonstyle == 3) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj5(i,j) = noncoeff3(i,j)
            cutljsq(i,j) = (noncoeff4(i,j)+noncoeff3(i,j))*
     $           (noncoeff4(i,j)+noncoeff3(i,j))
            if (offsetflag == 0 .or. noncoeff4(i,j) == 0.0) then
              offset(i,j) = 0.0
            else
              ratio = noncoeff2(i,j)/noncoeff4(i,j)
              offset(i,j) = 4.0*noncoeff1(i,j)*(ratio**12 - ratio**6)
            endif
            cut = max(cutcoul,noncoeff4(i,j)+noncoeff3(i,j))
          else if (nonstyle == 4) then
            lj1(i,j) = noncoeff1(i,j)
            lj2(i,j) = noncoeff2(i,j)
            lj3(i,j) = noncoeff3(i,j)
            cutljsq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            cut = noncoeff3(i,j)
          else if (nonstyle == 5) then
            lj1(i,j) = 18.0*noncoeff1(i,j)*noncoeff2(i,j)**9
            lj2(i,j) = 18.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 2.0*noncoeff1(i,j)*noncoeff2(i,j)**9
            lj4(i,j) = 3.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            cutljsq(i,j) = noncoeff3(i,j)*noncoeff3(i,j)
            if (offsetflag == 0 .or. noncoeff3(i,j) == 0.0) then
              offset(i,j) = 0.0
            else
              ratio = noncoeff2(i,j)/noncoeff3(i,j)
              offset(i,j) = noncoeff1(i,j) *
     $             (2.0*ratio**9 - 3.0*ratio**6)
            endif
            cut = max(cutcoul,noncoeff3(i,j))
          else if (nonstyle == 6) then
            lj1(i,j) = 48.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj2(i,j) = 24.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj3(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**12
            lj4(i,j) = 4.0*noncoeff1(i,j)*noncoeff2(i,j)**6
            lj14_1(i,j) = 48.0*noncoeff14_1(i,j)*noncoeff14_2(i,j)**12
            lj14_2(i,j) = 24.0*noncoeff14_1(i,j)*noncoeff14_2(i,j)**6
            lj14_3(i,j) = 4.0*noncoeff14_1(i,j)*noncoeff14_2(i,j)**12
            lj14_4(i,j) = 4.0*noncoeff14_1(i,j)*noncoeff14_2(i,j)**6
            cutljsq(i,j) = cutlj*cutlj
            cut = max(cutcoul,cutlj)
          endif

          cutforcesq(i,j) = cut*cut
          cutneighsq(i,j) = (cut+skin)*(cut+skin)
          cutforce = max(cutforce,cut)

        enddo
      enddo

      cutforce_sq = cutforce*cutforce
      cutneigh = cutforce + skin

c check if atom allocation used a sufficiently large cutoff value 
c  to estimate # of ghosts

      if (cutneigh > cutmemory .and. node == 0) then
        write (6,*) 'WARNING: Actual cutoff > memory-allocation cutoff'
        write (6,*) '         May run out of ghost-atom memory'
        write (6,*) '         Should use "maximum cutoff" command'
      endif

c zero out all potentials where user-specified cutoff = 0.0
c  setting values to -1.0 means no interaction will ever be computed
c   in force or neighbor routines

      do i = 1,ntypes
        do j = i,ntypes
          if (cutljsq(i,j) == 0.0) cutljsq(i,j) = -1.0
          if (nonstyle == 2) then
            if (cutljinnersq(i,j) == 0.0) cutljinnersq(i,j) = -1.0
          endif
          if (cutforcesq(i,j) == 0.0) then
            cutforcesq(i,j) = -1.0
            cutneighsq(i,j) = -1.0
          endif
        enddo
      enddo

      if (cutforce == 0.0) cutneigh = 0.0

c symmetrize nonbond j-i interaction to be equivalent to i-j
c  don't need to do this for noncoeff arrays since not used in force/energy

      do i = 1,ntypes
        do j = 1,i-1

          if (nonstyle == 0) then
            cutljsq(i,j) = cutljsq(j,i)
          else if (nonstyle == 1) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle == 2) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            cutljinner(i,j) = cutljinner(j,i)
            cutljinnersq(i,j) = cutljinnersq(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            ljsw0(i,j) = ljsw0(j,i)
            ljsw1(i,j) = ljsw1(j,i)
            ljsw2(i,j) = ljsw2(j,i)
            ljsw3(i,j) = ljsw3(j,i)
            ljsw4(i,j) = ljsw4(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle == 3) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            lj5(i,j) = lj5(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle == 4) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            cutljsq(i,j) = cutljsq(j,i)
          else if (nonstyle == 5) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            cutljsq(i,j) = cutljsq(j,i)
            offset(i,j) = offset(j,i)
          else if (nonstyle == 6) then
            lj1(i,j) = lj1(j,i)
            lj2(i,j) = lj2(j,i)
            lj3(i,j) = lj3(j,i)
            lj4(i,j) = lj4(j,i)
            lj14_1(i,j) = lj14_1(j,i)
            lj14_2(i,j) = lj14_2(j,i)
            lj14_3(i,j) = lj14_3(j,i)
            lj14_4(i,j) = lj14_4(j,i)
            cutljsq(i,j) = cutljsq(j,i)
          endif

          cutforcesq(i,j) = cutforcesq(j,i)
          cutneighsq(i,j) = cutneighsq(j,i)

        enddo
      enddo

c check bond lengths against max neighbor cutoff

      if (nbonds > 0) then

        if (bondstyle == 1) then
          do i = 1,nbondtypes
            if (cutneigh < bondcoeff(2,i)) then
              if (node == 0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < average length of bond type',i
            endif
          enddo
        else if (bondstyle == 2) then
          do i = 1,nbondtypes
            if (cutneigh < bondcoeff(2,i)) then
              if (node == 0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < maximum length of bond type',i
            endif
          enddo
        else if (bondstyle == 3) then
          do i = 1,nbondtypes
            if (cutneigh < bondcoeff(2,i)+bondcoeff(5,i)) then
              if (node == 0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < maximum length of bond type',i
            endif
          enddo
        else if (bondstyle == 4) then
          do i = 1,nbondtypes
            if (cutneigh < bondcoeff(2,i)+bondcoeff(3,i)) then
              if (node == 0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < maximum length of bond type',i
            endif
          enddo
        else if (bondstyle == 5) then
          do i = 1,nbondtypes
            if (cutneigh < bondcoeff(1,i)) then
              if (node == 0) write (6,*) 'WARNING: ',
     $             'Neighbor cutoff < average length of bond type',i
            endif
          enddo
        endif

      endif

c include dielectric constant in Coulomb prefactor

      coulpre = efactor/dielectric

c change special coeffs for Ewald/PPPM

      if (coulstyle == 3 .or. coulstyle == 4) then
        special(1) = special(1) + 2.0
        special(2) = special(2) + 2.0
        special(3) = special(3) + 2.0
      endif

c zero out velocity flags for subsequent runs

      do i = 1,nlocal
        velflag(i) = 0
      enddo

c timestep sizes

      dthalf = dt/2.0
      dthalf_intra = dthalf * nstretch
      dthalf_short = dthalf * nintra*nstretch
      dthalf_long = dthalf * nshort*nintra*nstretch

c ensemble type: (1) NVE, (2) NVT, (3) NPH, (4) NPT

      if (tempstyle /= 4 .and. pressstyle == 0) ensemble = 1
      if (tempstyle == 4 .and. pressstyle == 0) ensemble = 2
      if (tempstyle /= 4 .and. pressstyle > 0) ensemble = 3
      if (tempstyle == 4 .and. pressstyle > 0) ensemble = 4

c re-initialize Nose/Hoover eta and omega
c zero any component not being integrated on this run

      if (tempstyle /= 4) then
        eta = 0.0
        eta_dot = 0.0
      endif

      if (pressstyle == 0) then
        omega(1) = 0.0
        omega(2) = 0.0
        omega(3) = 0.0
        omega_dot(1) = 0.0
        omega_dot(2) = 0.0
        omega_dot(3) = 0.0
      endif

      if (pressstyle > 0) then
        if (p_freq(1) == 0.0) omega(1) = 0.0
        if (p_freq(1) == 0.0) omega_dot(1) = 0.0
        if (p_freq(2) == 0.0) omega(2) = 0.0
        if (p_freq(2) == 0.0) omega_dot(2) = 0.0
        if (p_freq(3) == 0.0) omega(3) = 0.0
        if (p_freq(3) == 0.0) omega_dot(3) = 0.0
      endif

c store initial box size if volume control enabled

      if (volstyle > 0) then
        volstart_xlo = box(1,1)
        volstart_xhi = box(2,1)
        volstart_ylo = box(1,2)
        volstart_yhi = box(2,2)
        volstart_zlo = box(1,3)
        volstart_zhi = box(2,3)
      endif

c initialize time counters
c  used by soft potential, fixes, and diagnostics

      itime = 0
      ntime_last = ntimestep + nsteps

c initialize counters before 1st communication, neighboring, I/O

      max_nlocal = 0
      max_nghost = 0
      max_neigh = 0
      max_slist = 0
      max_exch = 0
      max_buf = 0

      max_bond = 0
      max_angle = 0
      max_dihed = 0
      max_impro = 0

      numneigh = 0
      time_io = 0.0

c initialize all energies
c  in case force routines never called so thermo routine has good values

      e_long = 0.0
      e_vdwl = 0.0
      e_coul = 0.0

      e_bond = 0.0
      e_angle = 0.0
      e_dihedral = 0.0
      e_improper = 0.0

      e_bondbond = 0.0
      e_bondangle = 0.0
      e_midbondtorsion = 0.0
      e_endbondtorsion = 0.0
      e_angletorsion = 0.0
      e_angleangletorsion = 0.0
      e_bondbond13 = 0.0
      e_angleangle = 0.0

c allocate fhold storage if needed

      if (nrespa == 0 .and. newton == 3) allocate(fhold(3,maxatom))

c allocate rRESPA force arrays

      if (nrespa == 1) then
        allocate(f_stretch(3,maxown))
        allocate(f_intra(3,maxown))
        allocate(f_short(3,maxown))
        allocate(f_long(3,maxown))
      endif

c create special lists of 1-2, 1-3, 1-4 neighbors and SHAKE neighbors
      
      call special_create
      if (nshake == 1) then
        call shake_allocate
        call shake_create
      endif

c do PBC remapping since this may be continuation of a previous run
c if system is nonperiodic (in any dim) shrink-wrap global box around atoms
c setup communication pattern for current box
c setup neighbor bins
c acquire nearby atoms

      call pbc
      if (nonperiodic) call setup_box
      call setup_comm
      if (neighstyle == 1) call setup_neigh
      call exchange
      call borders

c create initial neighbor lists

      call neighbor

      if (nmolecular > 0) then
        if (nbonds > 0) call neighbor_bond
        if (nangles > 0) call neighbor_angle
        if (ndihedrals > 0) call neighbor_dihedral
        if (nimpropers > 0) call neighbor_improper
        if (nshake == 1) call neighbor_shake
      endif

c should virial be computed inside force routines or by store_virial
c use store_virial only if newton is fully on so ghost forces are accurate
c do not use store_virial if RESPA is enabled

      vflag = 1
      if (newton == 3 .and. nrespa == 0) vflag = 2

c first force calls to initialize velocity Verlet integrator

      call zero_force_virial

      if (nbonds > 0) then
        if (bondstyle == 1) then
          call bond_harmonic(vflag)
        else if (bondstyle == 2) then
          call bond_fene_standard(vflag)
        else if (bondstyle == 3) then
          call bond_fene_shift(vflag)
        else if (bondstyle == 4) then
          call bond_nonlinear(vflag)
        else if (bondstyle == 5) then
          call bond_class2(vflag)
        endif
      endif

      if (nrespa == 1) then
        if (newton_bond == 1) call reverse_comm
        call copy_force(nlocal,f_stretch)
        call copy_virial(nlocal,vir_stretch)
        if (nfixes_respa > 0) call fix_apply(0,f_stretch,dthalf)
        call zero_force_virial
      endif

      if (nangles > 0) then
        if (anglestyle == 1) then
          call angle_harmonic(vflag)
        else if (anglestyle == 2) then
          call angle_class2(vflag)
        else if (anglestyle == 3) then
          call angle_charmm(vflag)
        else if (anglestyle == 4) then
          call angle_cosine(vflag)
        endif
      endif

      if (ndihedrals > 0) then
        if (dihedstyle == 1) then
          call dihedral_harmonic(vflag)
        else if (dihedstyle == 2) then
          call dihedral_class2(vflag)
        else if (dihedstyle == 3) then
          call dihedral_multiharmonic(vflag)
        else if (dihedstyle == 4) then
          call dihedral_charmm(vflag)
        endif
      endif

      if (nimpropers > 0) then
        if (improstyle == 1) then
          call improper_harmonic(vflag)
        else if (improstyle == 2) then
          call improper_cvff(vflag)
        else if (improstyle == 3) then
          call improper_class2(vflag)
          call angleangle_class2(vflag)
        endif
      endif
      
      if (nrespa == 1) then
        if (newton_bond == 1) call reverse_comm
        call copy_force(nlocal,f_intra)
        call copy_virial(nlocal,vir_intra)
        if (nfixes_respa > 0) call fix_apply(0,f_intra,dthalf_intra)
        call zero_force_virial
      endif

      if (vflag == 2) then
        call copy_force(nlocal+nghost,fhold)
        call copy_virial(nlocal+nghost,virialhold)
        call zero_force_virial
      endif

      if (nonstyle == 0) then
        if (coulstyle == 0) then
          continue
        else if (coulstyle == 1) then
          call lj_cut_coul_cut(vflag)
        else if (coulstyle == 2) then
          call lj_cut_coul_smooth(vflag)
        else if (coulstyle == 3) then
          call ewald_coeff
          call lj_cut_coul_long(vflag)
        else if (coulstyle == 4) then
          call pppm_coeff(0)
          call lj_cut_coul_long(vflag)
        else if (coulstyle == 6) then
          call lj_cut_coul_debye(vflag)
        endif
      else if (nonstyle == 1) then
        if (coulstyle == 0) then
          call lj_cut(vflag)
        else if (coulstyle == 1) then
          call lj_cut_coul_cut(vflag)
        else if (coulstyle == 2) then
          call lj_cut_coul_smooth(vflag)
        else if (coulstyle == 3) then
          call ewald_coeff
          call lj_cut_coul_long(vflag)
        else if (coulstyle == 4) then
          call pppm_coeff(0)
          call lj_cut_coul_long(vflag)
        else if (coulstyle == 6) then
          call lj_cut_coul_debye(vflag)
        endif
      else if (nonstyle == 2) then
        if (coulstyle == 0) then
          call lj_smooth(vflag)
        else if (coulstyle == 1) then
          call lj_smooth_coul_cut(vflag)
        else if (coulstyle == 2) then
          call lj_smooth_coul_smooth(vflag)
        else if (coulstyle == 3) then
          call ewald_coeff
          call lj_smooth_coul_long(vflag)
        else if (coulstyle == 4) then
          call pppm_coeff(0)
          call lj_smooth_coul_long(vflag)
        endif
      else if (nonstyle == 3) then
        call lj_shift(vflag)
      else if (nonstyle == 4) then
        call soft(vflag)
      else if (nonstyle == 5) then
        if (coulstyle == 0) then
          call ljclass2_cut(vflag)
        else if (coulstyle == 1) then
          call ljclass2_cut_coul_cut(vflag)
        else if (coulstyle == 3) then
          call ewald_coeff
          call ljclass2_cut_coul_long(vflag)
        else if (coulstyle == 4) then
          call pppm_coeff(0)
          call ljclass2_cut_coul_long(vflag)
        endif
      else if (nonstyle == 6) then
        if (coulstyle == 3) then
          call ewald_coeff
          call lj_charmm_coul_long(vflag)
        else if (coulstyle == 4) then
          call pppm_coeff(0)
          call lj_charmm_coul_long(vflag)
        else if (coulstyle == 5) then
          call lj_charmm_coul_charmm(vflag)
        endif
      endif

      if (nrespa == 1) then
        if (newton_nonbond == 1) call reverse_comm
        call copy_force(nlocal,f_short)
        call copy_virial(nlocal,vir_short)
        if (tempstyle == 3) call langevin(f_short,dthalf_short)
        if (nfixes > 0) call fix_apply(1,f_short,dthalf_short)
        call zero_force_virial
      endif

      if (vflag == 2) call store_virial

      if (coulstyle == 3) then
        call ewald(vflag)
      else if (coulstyle == 4) then
        time_rho = 0.0
        time_poiss = 0.0
        time_field = 0.0
        call pppm(vflag)
      endif

      if (nrespa == 0) then
        if (newton >= 1) call reverse_comm
      else if (nrespa == 1) then
        call copy_force(nlocal,f_long)
        if (nfixes_respa > 0) call fix_apply(0,f_long,dthalf_long)
      endif

      if (nrespa == 0 .and. tempstyle == 3) call langevin(f,dthalf)
      if (nrespa == 0 .and. nfixes > 0) call fix_apply(1,f,dthalf)

      if (nshake == 1) then
        if (nshakestats > 0) call shake_bond_stats
        if (nshakestats == 0) nshake_next = ntime_last + 1
        call shake(0,vflag)
      endif

c all permanent arrays have now been allocated
c print out dynamic memory usage

      call memory_usage

c print run header

      if (node == 0) then
        write (6,*)
        write (1,*)
        if (thermostyle == 1) then
          write (6,*) 'Step Temp E_nbond E_bond E_long ',
     $         'E_total Pressure Volume'
          write (1,*) 'Step Temp E_nbond E_bond E_long ',
     $         'E_total Pressure Volume'
        endif
      endif

c initial thermodynamics

      call thermo(-1)
      if (nthermo == 0) nthermo_next = ntime_last

c initial dumps and restart file

      if (ndumpatom > 0) call dump_atom(0)
      if (ndumpatom == 0) ndumpatom_next = ntime_last + 1
      if (ndumpvel > 0) call dump_vel
      if (ndumpvel == 0) ndumpvel_next = ntime_last + 1
      if (ndumpforce > 0) call dump_force(0)
      if (ndumpforce == 0) ndumpforce_next = ntime_last + 1

      if (nrestart > 0)
     $     nrestart_next = min(ntimestep+nrestart,ntime_last)
      if (nrestart == 0) nrestart_next = ntime_last + 1

c initial diagnostics
c set all diagnext to this timestep to trigger calling of user routines

      do idiag = 1,numdiag
        diagnext(idiag) = ntimestep
      enddo
      if (numdiag > 0) call diagnostic
      if (numdiag == 0) ndiag_next = ntime_last + 1

      noutput_next = nthermo_next
      noutput_next = min(noutput_next,ndumpatom_next)
      noutput_next = min(noutput_next,ndumpvel_next)
      noutput_next = min(noutput_next,ndumpforce_next)
      noutput_next = min(noutput_next,nrestart_next)
      noutput_next = min(noutput_next,ndiag_next)
      
c zero out timers and counters

      time_nonbond = 0.0
      time_long = 0.0
      time_rho = 0.0
      time_poiss = 0.0
      time_field = 0.0

      time_bond = 0.0
      time_angle = 0.0
      time_dihedral = 0.0
      time_improper = 0.0

      time_neigh1 = 0.0
      time_neigh2 = 0.0
      time_exch = 0.0

      time_comm = 0.0
      time_fcomm = 0.0

      time_io = 0.0
      time_shake = 0.0

      numneigh = 0
      ndanger = 0

      return
      end
