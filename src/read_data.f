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
c read a data file

      subroutine read_data
      use global
      use mpi
      implicit none

c local variables

      logical match
      integer ierror,i,j,m
      integer massflag,natomflag,nbondflag,nangleflag
      integer ndihedflag,nimproflag,noncoeffflag,nbondcoeffflag
      integer nanglecoeffflag,ndihedcoeffflag,nimprocoeffflag
      integer nbondbondflag,nbondangleflag,nmidbondtorsionflag
      integer nendbondtorsionflag,nangletorsionflag
      integer nangleangletorsionflag,nbondbond13flag
      integer nangleangleflag
      real*8 small,rtmp
      character*80 str
      parameter (small=1.0D-4)

 900  format (a)

      call mpi_barrier(mpi_comm_world,ierror)

c scan data file to count max topology per atom
c only done by proc 0, then broadcast the results

      if (node == 0) call scan_data(datafile,
     $     idimension,newton_bond,
     $     maxbondper,maxangleper,maxdihedper,maximproper)

      call mpi_bcast(maxbondper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(maxangleper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(maxdihedper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(maximproper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

c read the data file

      call mpi_barrier(mpi_comm_world,ierror)
      if (node == 0) write (6,*) 'Reading data file ...'

      readflag = 1
      ntimestep = 0

      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0

c read data file header - fixed format

      if (node.eq.0) then
        open (unit=2,file=datafile,status='old',iostat=ierror)
        if (ierror /= 0) call abort('No data file')
        read (2,*)
        read (2,*)
        read (2,*) natoms
        read (2,*) nbonds
        read (2,*) nangles
        read (2,*) ndihedrals
        read (2,*) nimpropers
        read (2,*)
        read (2,*) ntypes
        if (nbonds.gt.0) read (2,*) nbondtypes
        if (nangles.gt.0) read (2,*) nangletypes
        if (ndihedrals.gt.0) read (2,*) ndihedtypes
        if (nimpropers.gt.0) read (2,*) nimprotypes
        read (2,*)
        read (2,*) box(1,1),box(2,1)
        read (2,*) box(1,2),box(2,2)
        if (idimension == 3) then
          read (2,*) box(1,3),box(2,3)
        else
          box(1,3) = -0.1
          box(2,3) = 0.1
        endif
      endif

c broadcast values

      call mpi_bcast(natoms,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nbonds,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nangles,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(ndihedrals,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nimpropers,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(ntypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nbondtypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nangletypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(ndihedtypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nimprotypes,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(box,6,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

c initialize ensemble variables

      eta = 0.0
      eta_dot = 0.0

      omega(1) = 0.0
      omega(2) = 0.0
      omega(3) = 0.0
      omega_dot(1) = 0.0
      omega_dot(2) = 0.0
      omega_dot(3) = 0.0

c error check

      if (natoms.le.0) call err('Atoms <= 0')
      if (nbonds.lt.0) call err('Bonds < 0')
      if (nangles.lt.0) call err('Angles < 0')
      if (ndihedrals.lt.0) call err('Dihedrals < 0')
      if (nimpropers.lt.0) call err('Impropers < 0')

      if (ntypes.le.0) call err('Atom types <= 0')
      if (nbonds.gt.0.and.nbondtypes.le.0)
     $     call err('Bond types <= 0')
      if (nangles.gt.0.and.nangletypes.le.0)
     $     call err('Angle types <= 0')
      if (ndihedrals.gt.0.and.ndihedtypes.le.0)
     $     call err('Dihedral types <= 0')
      if (nimpropers.gt.0.and.nimprotypes.le.0)
     $     call err('Improper types <= 0')

      if (box(1,1) >= box(2,1) .or. box(1,2) >= box(2,2) .or.
     $     box(1,3) >= box(2,3) ) call err('Illegal box dimensions')

      if (idimension.eq.2.and.pgrid(1).gt.0) then
        if (pgrid(3).ne.1)
     $       call err('Can only use one z-processor for 2-d')
      endif

c if non-PBC augment box bounds by small
c if user input box is exact, this insures every atom
c   is inside (not on surface) of some proc's sub-domain,
c   else can lose atoms when read-in later in read_data
c if using 2-d Ewald/PPPM, do not adjust z box bounds, even if non-periodic

      if (perflagx == 1) then
        box(1,1) = box(1,1) - small
        box(2,1) = box(2,1) + small
      endif
      if (perflagy == 1) then
        box(1,2) = box(1,2) - small
        box(2,2) = box(2,2) + small
      endif
      if (perflagz == 1 .and. slabflag == 0) then
        box(1,3) = box(1,3) - small
        box(2,3) = box(2,3) + small
      endif

c define simulation box

      xprd = box(2,1) - box(1,1)
      yprd = box(2,2) - box(1,2)
      zprd = box(2,3) - box(1,3)

      xprd_half = 0.5 * xprd
      yprd_half = 0.5 * yprd
      zprd_half = 0.5 * zprd

c define processor grid topology

      call mesh_3d(node,nprocs,idimension,xprd,yprd,zprd,
     $     pgrid,me,mpart)

      if (node == 0) then
        write (6,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
        write (1,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
      endif

c define local sub-box

      border(1,1) = float(me(1))/pgrid(1) * xprd + box(1,1)
      border(2,1) = float(me(1)+1)/pgrid(1) * xprd + box(1,1)
      border(1,2) = float(me(2))/pgrid(2) * yprd + box(1,2)
      border(2,2) = float(me(2)+1)/pgrid(2) * yprd + box(1,2)
      border(1,3) = float(me(3))/pgrid(3) * zprd + box(1,3)
      border(2,3) = float(me(3)+1)/pgrid(3) * zprd + box(1,3)

c allocate atom memory

      call atom_memory

c set atom counters

      nlocal = 0
      nghost = 0
      
c initialize global -> local pointers

      do i = 1,natoms
        localptr(i) = 0
      enddo

c read remainder of data file - free format

      massflag = 0
      natomflag = 0
      nbondflag = 0
      nangleflag = 0
      ndihedflag = 0
      nimproflag = 0

      noncoeffflag = 0
      nbondcoeffflag = 0
      nanglecoeffflag = 0
      ndihedcoeffflag = 0
      nimprocoeffflag = 0

      nbondbondflag = 0
      nbondangleflag = 0
      nmidbondtorsionflag = 0
      nendbondtorsionflag = 0
      nangletorsionflag = 0
      nangleangletorsionflag = 0
      nbondbond13flag = 0
      nangleangleflag = 0

      do

c read identifier string

        if (node.eq.0) then
          read (2,*,end=999,err=999)
          read (2,900,end=999,err=999) str
          read (2,*,end=999,err=999)
        endif
        call mpi_bcast(str,80,mpi_character,0,mpi_comm_world,ierror)

c end of file for all procs except 0

        if (match('All Done',str,m)) then

          goto 999

c read and broadcast atom masses

        else if (match('Masses',str,m)) then

          call read_mass(massflag)

c read and broadcast atoms

        else if (match('Atoms',str,m)) then

          call read_atoms(natomflag)

c read and broadcast velocities

        else if (match('Velocities',str,m)) then

          if (natomflag == 0)
     $         call err('Atoms entry must appear before Velocities')
          call read_vels

c read and broadcast bonds

        else if (match('Bonds',str,m)) then

          if (natomflag == 0)
     $         call err('Atoms entry must appear before Bonds')
          call read_bonds(nbondflag)

c read and broadcast angles

        else if (match('Angles',str,m)) then

          if (natomflag == 0)
     $         call err('Atoms entry must appear before Angles')
          call read_angles(nangleflag)

c read and broadcast dihedrals

        else if (match('Dihedrals',str,m)) then

          if (natomflag == 0)
     $         call err('Atoms entry must appear before Dihedrals')
          call read_dihedrals(ndihedflag)

c read and broadcast impropers

        else if (match('Impropers',str,m)) then

          if (natomflag == 0)
     $         call err('Atoms entry must appear before Impropers')
          call read_impropers(nimproflag)

c read and broadcast nonbond coeffs

        else if (match('Nonbond Coeffs',str,m)) then

          call read_nonbond_coeff(noncoeffflag)

c read and broadcast bond coeffs

        else if (match('Bond Coeffs',str,m)) then

          call read_bond_coeff(nbondcoeffflag)

c read and broadcast angle coeffs

        else if (match('Angle Coeffs',str,m)) then

          call read_angle_coeff(nanglecoeffflag)

c read and broadcast dihedral coeffs

        else if (match('Dihedral Coeffs',str,m)) then

          call read_dihedral_coeff(ndihedcoeffflag)

c read and broadcast improper coeffs

        else if (match('Improper Coeffs',str,m)) then

          call read_improper_coeff(nimprocoeffflag)

c read and broadcast bond-bond coeffs for class 2 force field

        else if (match('BondBond Coeffs',str,m)) then

          call read_bondbond_coeff(nbondbondflag)

c read and broadcast bond-angle coeffs for class 2 force field

        else if (match('BondAngle Coeffs',str,m)) then

          call read_bondangle_coeff(nbondangleflag)

c read and broadcast class 2 middle-bond-torsion coeffs

        else if (match('MiddleBondTorsion Coeffs',str,m)) then

          call read_midbondtorsion_coeff(nmidbondtorsionflag)

c read and broadcast class 2 end-bond-torsion coeffs

        else if (match('EndBondTorsion Coeffs',str,m)) then

          call read_endbondtorsion_coeff(nendbondtorsionflag)

c read and broadcast class 2 angle-torsion coeffs

        else if (match('AngleTorsion Coeffs',str,m)) then

          call read_angletorsion_coeff(nangletorsionflag)

c read and broadcast class 2 angle-angle-torsion coeffs

        else if (match('AngleAngleTorsion Coeffs',str,m)) then

          call read_angleangletorsion_coeff(nangleangletorsionflag)

c read and broadcast class 2 bondbond13 coeffs

        else if (match('BondBond13 Coeffs',str,m)) then

          call read_bondbond13_coeff(nbondbond13flag)

c read and broadcast class 2 angle-angle coeffs

        else if (match('AngleAngle Coeffs',str,m)) then

          call read_angleangle_coeff(nangleangleflag)

c unknown identifier

        else
          if (node == 0) write (6,*) 'UNKNOWN: ',trim(str)
          call err('Unknown identifier in data file')
        endif

      enddo

c close data file

 999  continue

      if (node.eq.0) then
        str = 'All Done'
        call mpi_bcast(str,80,mpi_character,0,mpi_comm_world,ierror)
        close (2)
      endif

c check that all required options were read in from data file

      if (massflag.eq.0) call err('No Masses in data file')
      if (natomflag.eq.0) call err('No Atoms in data file')

      if (nbonds.gt.0.and.nbondflag.eq.0)
     $     call err('No Bonds in data file')
      if (nangles.gt.0.and.nangleflag.eq.0)
     $     call err('No Angles in data file')
      if (ndihedrals.gt.0.and.ndihedflag.eq.0)
     $     call err('No Dihedrals in data file')
      if (nimpropers.gt.0.and.nimproflag.eq.0)
     $     call err('No Impropers in data file')

      if (nbonds.gt.0.and.bondstyle.eq.5.and.nbondcoeffflag.eq.0)
     $     call err('No Bond Coeffs in data file')
      if (nangles.gt.0.and.anglestyle.eq.2.and.nanglecoeffflag.eq.0)
     $     call err('No Angle Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.ndihedcoeffflag.eq.0)
     $     call err('No Dihedral Coeffs in data file')
      if (nimpropers.gt.0.and.improstyle.eq.3.and.nimprocoeffflag.eq.0)
     $     call err('No Improper Coeffs in data file')

      if (nangles.gt.0.and.anglestyle.eq.2.and.nbondbondflag.eq.0)
     $     call err('No BondBond Coeffs in data file')
      if (nangles.gt.0.and.anglestyle.eq.2.and.nbondangleflag.eq.0)
     $     call err('No BondAngle Coeffs in data file')

      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nmidbondtorsionflag.eq.0)
     $     call err('No MiddleBondTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nendbondtorsionflag.eq.0)
     $     call err('No EndBondTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nangletorsionflag.eq.0)
     $     call err('No AngleTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nangleangletorsionflag.eq.0)
     $     call err('No AngleAngleTorsion Coeffs in data file')
      if (ndihedrals.gt.0.and.dihedstyle.eq.2.and.
     $     nbondbond13flag.eq.0)
     $     call err('No BondBond13 Coeffs in data file')

      if (nimpropers.gt.0.and.improstyle.eq.3.and.
     $     nangleangleflag.eq.0)
     $     call err('No AngleAngle Coeffs in data file')

c compute nonbond mixing rules (geometric,arithmetic,sixthpower)
c set nonbond cutoffs to current settings
c  only if nonbond coeffs were input

      if (noncoeffflag.eq.1) then

        do i = 1,ntypes-1
          do j = i+1,ntypes

c mixing of epsilon

            if (mixstyle.eq.1.or.mixstyle.eq.2.or.nonstyle.eq.4.or.
     $           nonstyle.eq.6) then
              noncoeff1(i,j) = sqrt(noncoeff1(i,i)*noncoeff1(j,j))
              if (nonstyle.eq.6) noncoeff14_1(i,j) = 
     $             sqrt(noncoeff14_1(i,i)*noncoeff14_1(j,j))
            else if (mixstyle.eq.3) then
              noncoeff1(i,j) = 2.0 *
     $             sqrt(noncoeff1(i,i)*noncoeff1(j,j)) *
     $             noncoeff2(i,i)**3 * noncoeff2(j,j)**3 /
     $             (noncoeff2(i,i)**6 + noncoeff2(j,j)**6)
            endif

c mixing of sigma

            if (mixstyle.eq.1.or.nonstyle.eq.4) then
              noncoeff2(i,j) = sqrt(noncoeff2(i,i)*noncoeff2(j,j))
            else if (mixstyle.eq.2.or.nonstyle.eq.6) then
              noncoeff2(i,j) = (noncoeff2(i,i)+noncoeff2(j,j))/2.0
              if (nonstyle.eq.6) noncoeff14_2(i,j) = 
     $             (noncoeff14_2(i,i) + noncoeff14_2(j,j))/2.0
            else if (mixstyle.eq.3) then
              noncoeff2(i,j) =
     $             ((noncoeff2(i,i)**6 + noncoeff2(j,j)**6) / 2.0) **
     $             (1.0/6.0)
            endif

c mixing of shift factor

            if (nonstyle.eq.3) noncoeff3(i,j) =
     $           (noncoeff3(i,i)+noncoeff3(j,j))/2.0

          enddo
        enddo

c set cutoffs

        do i = 1,ntypes
          do j = i,ntypes
            if (nonstyle == 1) then
              noncoeff3(i,j) = cutlj
            else if (nonstyle == 2) then
              noncoeff3(i,j) = cutljinterior
              noncoeff4(i,j) = cutlj
            else if (nonstyle == 3) then
              noncoeff4(i,j) = cutlj
            else if (nonstyle == 4) then
              noncoeff3(i,j) = cutlj
            else if (nonstyle == 5) then
              noncoeff3(i,j) = cutlj
            endif
          enddo
        enddo

c mark all type pairs as set

        do i = 1,ntypes
          do j = i,ntypes
            nontypeflag(i,j) = nonstyle
          enddo
        enddo
        
      endif

c perform unit conversions where necessary
c  convert angle theta from degrees to radians
c  convert dihedral phi from degrees to radians
c  convert improper omega from degrees to radians
c  convert class 2 angles from degrees to radians

      if (nanglecoeffflag.eq.1) then
        if (anglestyle.eq.1) then
          do i = 1,nangletypes
            anglecoeff(2,i) = anglecoeff(2,i)/180.0 * 3.1415926
          enddo
        else if (anglestyle.eq.2) then
          do i = 1,nangletypes
            anglecoeff(1,i) = anglecoeff(1,i)/180.0 * 3.1415926
          enddo
        else if (anglestyle.eq.3) then
          do i = 1,nangletypes
            anglecoeff(2,i) = anglecoeff(2,i)/180.0 * 3.1415926
          enddo
        endif
      endif

      if (ndihedcoeffflag.eq.1 .and. dihedstyle == 2) then
        do i = 1,ndihedtypes
          dihedcoeff(2,i) = dihedcoeff(2,i)/180.0 * 3.1415926
          dihedcoeff(4,i) = dihedcoeff(4,i)/180.0 * 3.1415926
          dihedcoeff(6,i) = dihedcoeff(6,i)/180.0 * 3.1415926
        enddo
      endif

      if ((improstyle == 1 .or. improstyle == 3) .and.
     $     nimprocoeffflag.eq.1) then
        do i = 1,nimprotypes
          improcoeff(2,i) = improcoeff(2,i)/180.0 * 3.1415926
        enddo
      endif

      if (dihedstyle.eq.2.and.nangletorsionflag.eq.1) then
        do i = 1,ndihedtypes
          angletorsioncoeff(7,i) =
     $         angletorsioncoeff(7,i)/180.0 * 3.1415926
          angletorsioncoeff(8,i) =
     $         angletorsioncoeff(8,i)/180.0 * 3.1415926
        enddo
      endif

      if (dihedstyle.eq.2.and.nangleangletorsionflag.eq.1) then
        do i = 1,ndihedtypes
          angleangletorsioncoeff(2,i) =
     $         angleangletorsioncoeff(2,i)/180.0 * 3.1415926
          angleangletorsioncoeff(3,i) =
     $         angleangletorsioncoeff(3,i)/180.0 * 3.1415926
        enddo
      endif

      if (improstyle.eq.3.and.nangleangleflag.eq.1) then
        do i = 1,nimprotypes
          angleanglecoeff(4,i) =
     $         angleanglecoeff(4,i)/180.0 * 3.1415926
          angleanglecoeff(5,i) =
     $         angleanglecoeff(5,i)/180.0 * 3.1415926
          angleanglecoeff(6,i) =
     $         angleanglecoeff(6,i)/180.0 * 3.1415926
        enddo
      endif

c compute total mass and charge

      masssum = 0.0
      qsum = 0.0
      qsqsum = 0.0

      do i = 1,nlocal
        masssum = masssum + mass(type(i))
        qsum = qsum + q(i)
        qsqsum = qsqsum + q(i)*q(i)
      enddo

      rtmp = masssum
      call mpi_allreduce(rtmp,masssum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      rtmp = qsum
      call mpi_allreduce(rtmp,qsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      rtmp = qsqsum
      call mpi_allreduce(rtmp,qsqsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read masses from data file

      subroutine read_mass(massflag)
      use global
      use mpi
      implicit none

c argument variables

      integer massflag

c local variables

      integer i,ierror,jtmp

      if (node.eq.0) write (6,*) 'Masses ...'
      massflag = 1

      if (node.eq.0) then
        do i = 1,ntypes
          read (2,*) jtmp,mass(i)
        enddo
      endif

      call mpi_bcast(mass,ntypes,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read atoms from data file

      subroutine read_atoms(natomflag)
      use global
      use mpi
      implicit none

c argument variables

      integer natomflag

c local variables

      integer i,j,itag,itrue,isize,nblock,nleft,nset,ierror
      integer iflagmax,natoms_tmp
      real*8, allocatable :: buf(:)

      if (node.eq.0) write (6,*) 'Atoms ...'
      natomflag = 1

c read nblock atoms at a time, do one broadcast

      isize = 7
      if (trueflag > 1) isize = 10
      nblock = 1000
      nleft = natoms
      allocate(buf(isize*nblock))

      do while (nleft > 0)
        nset = min(nleft,nblock)
        if (node == 0) read (2,*) (buf(j),j=1,isize*nset)
        call mpi_bcast(buf,isize*nset,mpi_double_precision,0,
     $       mpi_comm_world,ierror)

        j = 0
        do i = 1,nset
          if (trueflag <= 1) then
            itrue = 500500500
          else
            itrue = 500+nint(buf(j+8)) +
     $           (500+nint(buf(j+9)))*1000 +
     $           (500+nint(buf(j+10)))*1000000
          endif

c put atom in periodic box no matter how far away it is

          call remap(buf(j+5),buf(j+6),buf(j+7),itrue)

c add atom to my list if in my sub-box

          if (buf(j+5) >= border(1,1) .and.
     $         buf(j+5) < border(2,1).and.
     $         buf(j+6) >= border(1,2) .and.
     $         buf(j+6) < border(2,2).and.
     $         buf(j+7) >= border(1,3) .and.
     $         buf(j+7) < border(2,3)) then
            nlocal = nlocal + 1
            if (nlocal <= maxown) then
              itag = nint(buf(j+1))
              localptr(itag) = nlocal
              tag(nlocal) = itag
              molecule(nlocal) = nint(buf(j+2))
              type(nlocal) = nint(buf(j+3))
              q(nlocal) = buf(j+4)
              x(1,nlocal) = buf(j+5)
              x(2,nlocal) = buf(j+6)
              x(3,nlocal) = buf(j+7)
              true(nlocal) = itrue
              fix(nlocal) = 0
              numbond(nlocal) = 0
              numangle(nlocal) = 0
              numdihed(nlocal) = 0
              numimpro(nlocal) = 0
              velflag(nlocal) = 0
              v(1,nlocal) = 0.0
              v(2,nlocal) = 0.0
              v(3,nlocal) = 0.0
            endif
          endif
          j = j + isize
        enddo

        nleft = nleft - nset
      enddo

      deallocate(buf)

c check if too many atoms on any proc

      call mpi_allreduce(nlocal,iflagmax,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (iflagmax.gt.maxown) then
        if (node.eq.0)
     $       write (6,*) 'Nlocal,maxown =',iflagmax,maxown
        call err('Too many atoms/processor - boost extra_own')
      endif

c check that all atoms were assigned uniquely to a processor

      call mpi_allreduce(nlocal,natoms_tmp,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      if (natoms.ne.natoms_tmp)
     $     call err('Did not assign all atoms correctly')

      return
      end


c -------------------------------------------------------------------------
c read velocities from data file

      subroutine read_vels
      use global
      use mpi
      implicit none

c local variables

      integer i,j,m,isize,nblock,nleft,nset,ierror
      real*8, allocatable :: buf(:)

      if (node.eq.0) write (6,*) 'Velocities ...'

c read nblock vels at a time, do one broadcast

      isize = 4
      nblock = 1000
      nleft = natoms
      allocate(buf(isize*nblock))

      do while (nleft > 0)
        nset = min(nleft,nblock)
        if (node == 0) read (2,*) (buf(j),j=1,isize*nset)
        call mpi_bcast(buf,isize*nset,mpi_double_precision,0,
     $       mpi_comm_world,ierror)

        j = 0
        do i = 1,nset
          m = localptr(nint(buf(j+1)))
          if (m > 0) then
            velflag(m) = 1
            v(1,m) = buf(j+2)*dtfactor
            v(2,m) = buf(j+3)*dtfactor
            v(3,m) = buf(j+4)*dtfactor
          endif
          j = j + isize
        enddo

        nleft = nleft - nset
      enddo

      deallocate(buf)

      return
      end


c -------------------------------------------------------------------------
c read bond topology from data file

      subroutine read_bonds(nbondflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nbondflag

c local variables

      integer i,j,itrue,isize,nblock,nleft,nset,ierror,ipnt
      integer max_bondper,jtmp,nbonds_tmp,ifactor
      integer, allocatable :: ibuf(:)

      if (node.eq.0) write (6,*) 'Bonds ...'
      if (nbonds.eq.0.or.nbondtypes.eq.0)
     $     call err('Should be no Bonds entry')
      nbondflag = 1

c read nblock bonds at a time, do one broadcast

      isize = 4
      nblock = 1000
      nleft = nbonds
      allocate(ibuf(isize*nblock))

      do while (nleft > 0)
        nset = min(nleft,nblock)
        if (node == 0) read (2,*) (ibuf(j),j=1,isize*nset)
        call mpi_bcast(ibuf,isize*nset,mpi_integer,0,
     $       mpi_comm_world,ierror)

        j = 0
        do i = 1,nset

          ipnt = localptr(ibuf(j+3))
          if (ipnt /= 0) then
            numbond(ipnt) = numbond(ipnt) + 1
            if (numbond(ipnt) <= maxbondper) then
              bondtype(numbond(ipnt),ipnt) = ibuf(j+2)
              bondatom1(numbond(ipnt),ipnt) = ibuf(j+3)
              bondatom2(numbond(ipnt),ipnt) = ibuf(j+4)
            endif
          endif
          
          if (newton_bond == 0) then
            ipnt = localptr(ibuf(j+4))
            if (ipnt /= 0) then
              numbond(ipnt) = numbond(ipnt) + 1
              if (numbond(ipnt) <= maxbondper) then
                bondtype(numbond(ipnt),ipnt) = ibuf(j+2)
                bondatom1(numbond(ipnt),ipnt) = ibuf(j+3)
                bondatom2(numbond(ipnt),ipnt) = ibuf(j+4)
              endif
            endif
          endif

          j = j + isize
        enddo

        nleft = nleft - nset
      enddo

      deallocate(ibuf)

c check for mismatched bonds/atom

      max_bondper = 0
      do i = 1,nlocal
        if (numbond(i).gt.max_bondper) max_bondper = numbond(i)
      enddo

      jtmp = max_bondper
      call mpi_allreduce(jtmp,max_bondper,1,mpi_integer,
     $     mpi_max,mpi_comm_world,ierror)
      if (max_bondper /= maxbondper) then
        if (node.eq.0) write (6,*)
     $       'Max bonds/atom =',max_bondper,maxbondper
        call err('Mismatch in max bonds/atom')
      endif

c check that all bonds were assigned uniquely to a processor

      nbonds_tmp = 0
      do i = 1,nlocal
        nbonds_tmp = nbonds_tmp + numbond(i)
      enddo
      
      jtmp = nbonds_tmp
      call mpi_allreduce(jtmp,nbonds_tmp,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      ifactor = 1
      if (newton_bond.eq.0) ifactor = 2
      if (nbonds_tmp.ne.ifactor*nbonds)
     $     call err('Bonds assigned incorrectly')

      return
      end


c -------------------------------------------------------------------------
c read angle topology from data file

      subroutine read_angles(nangleflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nangleflag

c local variables

      integer i,j,isize,nblock,nleft,nset,ierror,ipnt
      integer max_angleper,jtmp,nangles_tmp,ifactor
      integer, allocatable :: ibuf(:)

      if (node.eq.0) write (6,*) 'Angles ...'
      if (nangles.eq.0.or.nangletypes.eq.0)
     $     call err('Should be no Angles entry')
      nangleflag = 1

c read nblock angles at a time, do one broadcast

      isize = 5
      nblock = 1000
      nleft = nangles
      allocate(ibuf(isize*nblock))

      do while (nleft > 0)
        nset = min(nleft,nblock)
        if (node == 0) read (2,*) (ibuf(j),j=1,isize*nset)
        call mpi_bcast(ibuf,isize*nset,mpi_integer,0,
     $       mpi_comm_world,ierror)

        j = 0
        do i = 1,nset

          ipnt = localptr(ibuf(j+4))
          if (ipnt.ne.0) then
            numangle(ipnt) = numangle(ipnt) + 1
            if (numangle(ipnt).le.maxangleper) then
              angletype(numangle(ipnt),ipnt) = ibuf(j+2)
              angleatom1(numangle(ipnt),ipnt) = ibuf(j+3)
              angleatom2(numangle(ipnt),ipnt) = ibuf(j+4)
              angleatom3(numangle(ipnt),ipnt) = ibuf(j+5)
            endif
          endif

          if (newton_bond.eq.0) then
            ipnt = localptr(ibuf(j+3))
            if (ipnt.ne.0) then
              numangle(ipnt) = numangle(ipnt) + 1
              if (numangle(ipnt).le.maxangleper) then
                angletype(numangle(ipnt),ipnt) = ibuf(j+2)
                angleatom1(numangle(ipnt),ipnt) = ibuf(j+3)
                angleatom2(numangle(ipnt),ipnt) = ibuf(j+4)
                angleatom3(numangle(ipnt),ipnt) = ibuf(j+5)
              endif
            endif
            
            ipnt = localptr(ibuf(j+5))
            if (ipnt.ne.0) then
              numangle(ipnt) = numangle(ipnt) + 1
              if (numangle(ipnt).le.maxangleper) then
                angletype(numangle(ipnt),ipnt) = ibuf(j+2)
                angleatom1(numangle(ipnt),ipnt) = ibuf(j+3)
                angleatom2(numangle(ipnt),ipnt) = ibuf(j+4)
                angleatom3(numangle(ipnt),ipnt) = ibuf(j+5)
              endif
            endif
          endif

          j = j + isize
        enddo

        nleft = nleft - nset
      enddo

      deallocate(ibuf)

c check for mismatched angles/atom

      max_angleper = 0
      do i = 1,nlocal
        if (numangle(i).gt.max_angleper) max_angleper = numangle(i)
      enddo

      jtmp = max_angleper
      call mpi_allreduce(jtmp,max_angleper,1,mpi_integer,
     $     mpi_max,mpi_comm_world,ierror)
      if (max_angleper /= maxangleper) then
        if (node.eq.0) write (6,*)
     $       'Max angles/atom =',max_angleper,maxangleper
        call err('Mismatch in max angles/atom')
      endif

c check that all angles were assigned uniquely to a processor

      nangles_tmp = 0
      do i = 1,nlocal
        nangles_tmp = nangles_tmp + numangle(i)
      enddo
      
      jtmp = nangles_tmp
      call mpi_allreduce(jtmp,nangles_tmp,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      ifactor = 1
      if (newton_bond.eq.0) ifactor = 3
      if (nangles_tmp.ne.ifactor*nangles)
     $     call err('Angles assigned incorrectly')

      return
      end


c -------------------------------------------------------------------------
c read dihedral topology from data file

      subroutine read_dihedrals(ndihedflag)
      use global
      use mpi
      implicit none

c argument variables

      integer ndihedflag

c local variables

      integer i,j,isize,nblock,nleft,nset,ierror,ipnt
      integer max_dihedper,jtmp,ndiheds_tmp,ifactor
      integer, allocatable :: ibuf(:)

      if (node.eq.0) write (6,*) 'Dihedrals ...'
      if (ndihedrals.eq.0.or.ndihedtypes.eq.0)
     $     call err('Should be no Dihedrals entry')
      ndihedflag = 1

c read nblock dihedrals at a time, do one broadcast

      isize = 6
      nblock = 1000
      nleft = ndihedrals
      allocate(ibuf(isize*nblock))

      do while (nleft > 0)
        nset = min(nleft,nblock)
        if (node == 0) read (2,*) (ibuf(j),j=1,isize*nset)
        call mpi_bcast(ibuf,isize*nset,mpi_integer,0,
     $       mpi_comm_world,ierror)

        j = 0
        do i = 1,nset

          ipnt = localptr(ibuf(j+4))
          if (ipnt.ne.0) then
            numdihed(ipnt) = numdihed(ipnt) + 1
            if (numdihed(ipnt).le.maxdihedper) then
              dihedtype(numdihed(ipnt),ipnt) = ibuf(j+2)
              dihedatom1(numdihed(ipnt),ipnt) = ibuf(j+3)
              dihedatom2(numdihed(ipnt),ipnt) = ibuf(j+4)
              dihedatom3(numdihed(ipnt),ipnt) = ibuf(j+5)
              dihedatom4(numdihed(ipnt),ipnt) = ibuf(j+6)
            endif
          endif

          if (newton_bond == 0) then
            ipnt = localptr(ibuf(j+3))
            if (ipnt.ne.0) then
              numdihed(ipnt) = numdihed(ipnt) + 1
              if (numdihed(ipnt).le.maxdihedper) then
                dihedtype(numdihed(ipnt),ipnt) = ibuf(j+2)
                dihedatom1(numdihed(ipnt),ipnt) = ibuf(j+3)
                dihedatom2(numdihed(ipnt),ipnt) = ibuf(j+4)
                dihedatom3(numdihed(ipnt),ipnt) = ibuf(j+5)
                dihedatom4(numdihed(ipnt),ipnt) = ibuf(j+6)
              endif
            endif

            ipnt = localptr(ibuf(j+5))
            if (ipnt.ne.0) then
              numdihed(ipnt) = numdihed(ipnt) + 1
              if (numdihed(ipnt).le.maxdihedper) then
                dihedtype(numdihed(ipnt),ipnt) = ibuf(j+2)
                dihedatom1(numdihed(ipnt),ipnt) = ibuf(j+3)
                dihedatom2(numdihed(ipnt),ipnt) = ibuf(j+4)
                dihedatom3(numdihed(ipnt),ipnt) = ibuf(j+5)
                dihedatom4(numdihed(ipnt),ipnt) = ibuf(j+6)
              endif
            endif
            
            ipnt = localptr(ibuf(j+6))
            if (ipnt.ne.0) then
              numdihed(ipnt) = numdihed(ipnt) + 1
              if (numdihed(ipnt).le.maxdihedper) then
                dihedtype(numdihed(ipnt),ipnt) = ibuf(j+2)
                dihedatom1(numdihed(ipnt),ipnt) = ibuf(j+3)
                dihedatom2(numdihed(ipnt),ipnt) = ibuf(j+4)
                dihedatom3(numdihed(ipnt),ipnt) = ibuf(j+5)
                dihedatom4(numdihed(ipnt),ipnt) = ibuf(j+6)
              endif
            endif
          endif

          j = j + isize
        enddo

        nleft = nleft - nset
      enddo

      deallocate(ibuf)

c check for mismatched dihedrals/atom

      max_dihedper = 0
      do i = 1,nlocal
        if (numdihed(i).gt.max_dihedper) max_dihedper = numdihed(i)
      enddo

      jtmp = max_dihedper
      call mpi_allreduce(jtmp,max_dihedper,1,mpi_integer,
     $     mpi_max,mpi_comm_world,ierror)
      if (max_dihedper /= maxdihedper) then
        if (node.eq.0) write (6,*)
     $       'Max dihedrals/atom =',max_dihedper,maxdihedper
        call err('Mismatch in max dihedrals/atom')
      endif

c check that all dihedrals were assigned uniquely to a processor

      ndiheds_tmp = 0
      do i = 1,nlocal
        ndiheds_tmp = ndiheds_tmp + numdihed(i)
      enddo
      
      jtmp = ndiheds_tmp
      call mpi_allreduce(jtmp,ndiheds_tmp,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      ifactor = 1
      if (newton_bond.eq.0) ifactor = 4
      if (ndiheds_tmp.ne.ifactor*ndihedrals)
     $     call err('Dihedrals assigned incorrectly')

      return
      end


c -------------------------------------------------------------------------
c read improper topology from data file

      subroutine read_impropers(nimproflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nimproflag

c local variables

      integer i,j,isize,nblock,nleft,nset,ierror,ipnt
      integer max_improper,jtmp,nimpros_tmp,ifactor
      integer, allocatable :: ibuf(:)

      if (node.eq.0) write (6,*) 'Impropers ...'
      if (nimpropers.eq.0.or.nimprotypes.eq.0)
     $     call err('Should be no Impropers entry')
      nimproflag = 1

c read nblock impropers at a time, do one broadcast

      isize = 6
      nblock = 1000
      nleft = nimpropers
      allocate(ibuf(isize*nblock))

      do while (nleft > 0)
        nset = min(nleft,nblock)
        if (node == 0) read (2,*) (ibuf(j),j=1,isize*nset)
        call mpi_bcast(ibuf,isize*nset,mpi_integer,0,
     $       mpi_comm_world,ierror)

        j = 0
        do i = 1,nset

          ipnt = localptr(ibuf(j+4))
          if (ipnt.ne.0) then
            numimpro(ipnt) = numimpro(ipnt) + 1
            if (numimpro(ipnt).le.maximproper) then
              improtype(numimpro(ipnt),ipnt) = ibuf(j+2)
              improatom1(numimpro(ipnt),ipnt) = ibuf(j+3)
              improatom2(numimpro(ipnt),ipnt) = ibuf(j+4)
              improatom3(numimpro(ipnt),ipnt) = ibuf(j+5)
              improatom4(numimpro(ipnt),ipnt) = ibuf(j+6)
            endif
          endif

          if (newton_bond.eq.0) then
            ipnt = localptr(ibuf(j+3))
            if (ipnt.ne.0) then
              numimpro(ipnt) = numimpro(ipnt) + 1
              if (numimpro(ipnt).le.maximproper) then
                improtype(numimpro(ipnt),ipnt) = ibuf(j+2)
                improatom1(numimpro(ipnt),ipnt) = ibuf(j+3)
                improatom2(numimpro(ipnt),ipnt) = ibuf(j+4)
                improatom3(numimpro(ipnt),ipnt) = ibuf(j+5)
                improatom4(numimpro(ipnt),ipnt) = ibuf(j+6)
              endif
            endif

            ipnt = localptr(ibuf(j+5))
            if (ipnt.ne.0) then
              numimpro(ipnt) = numimpro(ipnt) + 1
              if (numimpro(ipnt).le.maximproper) then
                improtype(numimpro(ipnt),ipnt) = ibuf(j+2)
                improatom1(numimpro(ipnt),ipnt) = ibuf(j+3)
                improatom2(numimpro(ipnt),ipnt) = ibuf(j+4)
                improatom3(numimpro(ipnt),ipnt) = ibuf(j+5)
                improatom4(numimpro(ipnt),ipnt) = ibuf(j+6)
              endif
            endif
            
            ipnt = localptr(ibuf(j+6))
            if (ipnt.ne.0) then
              numimpro(ipnt) = numimpro(ipnt) + 1
              if (numimpro(ipnt).le.maximproper) then
                improtype(numimpro(ipnt),ipnt) = ibuf(j+2)
                improatom1(numimpro(ipnt),ipnt) = ibuf(j+3)
                improatom2(numimpro(ipnt),ipnt) = ibuf(j+4)
                improatom3(numimpro(ipnt),ipnt) = ibuf(j+5)
                improatom4(numimpro(ipnt),ipnt) = ibuf(j+6)
              endif
            endif
          endif

          j = j + isize
        enddo

        nleft = nleft - nset
      enddo

      deallocate(ibuf)

c check for mismatched impropers/atom

      max_improper = 0
      do i = 1,nlocal
        if (numimpro(i).gt.max_improper) max_improper = numimpro(i)
      enddo

      jtmp = max_improper
      call mpi_allreduce(jtmp,max_improper,1,mpi_integer,
     $     mpi_max,mpi_comm_world,ierror)
      if (max_improper /= maximproper) then
        if (node.eq.0) write (6,*)
     $       'Max dihedrals/atom =',max_improper,maximproper
        call err('Mismatch in max impropers/atom')
      endif

c check that all impropers were assigned uniquely to a processor

      nimpros_tmp = 0
      do i = 1,nlocal
        nimpros_tmp = nimpros_tmp + numimpro(i)
      enddo
      
      jtmp = nimpros_tmp
      call mpi_allreduce(jtmp,nimpros_tmp,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      ifactor = 1
      if (newton_bond.eq.0) ifactor = 4
      if (nimpros_tmp.ne.ifactor*nimpropers)
     $     call err('Impropers assigned incorrectly')

      return
      end


c -------------------------------------------------------------------------
c read nonbond coeffs from data file

      subroutine read_nonbond_coeff(noncoeffflag)
      use global
      use mpi
      implicit none

c argument variables

      integer noncoeffflag

c local variables

      integer i,k,jtmp,ierror
      real*8 tmp(4)

      if (node.eq.0) write (6,*) 'Nonbond Coeffs ...'
      if (nonstyle.eq.0)
     $     call err('Nonbond Coeffs entry not allowed')
      noncoeffflag = 1

      if (nonstyle.eq.1) then
        do i = 1,ntypes
          if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
          call mpi_bcast(tmp,2,mpi_double_precision,0,
     $         mpi_comm_world,ierror)
          noncoeff1(i,i) = tmp(1)
          noncoeff2(i,i) = tmp(2)
        enddo
      else if (nonstyle.eq.2) then
        do i = 1,ntypes
          if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
          call mpi_bcast(tmp,2,mpi_double_precision,0,
     $         mpi_comm_world,ierror)
          noncoeff1(i,i) = tmp(1)
          noncoeff2(i,i) = tmp(2)
        enddo
      else if (nonstyle.eq.3) then
        do i = 1,ntypes
          if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,3)
          call mpi_bcast(tmp,3,mpi_double_precision,0,
     $         mpi_comm_world,ierror)
          noncoeff1(i,i) = tmp(1)
          noncoeff2(i,i) = tmp(2)
          noncoeff3(i,i) = tmp(3)
        enddo
      else if (nonstyle.eq.4) then
        do i = 1,ntypes
          if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
          call mpi_bcast(tmp,2,mpi_double_precision,0,
     $         mpi_comm_world,ierror)
          noncoeff1(i,i) = tmp(1)
          noncoeff2(i,i) = tmp(2)
        enddo
      else if (nonstyle.eq.5) then
        do i = 1,ntypes
          if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,2)
          call mpi_bcast(tmp,2,mpi_double_precision,0,
     $         mpi_comm_world,ierror)
          noncoeff1(i,i) = tmp(1)
          noncoeff2(i,i) = tmp(2)
        enddo
      else if (nonstyle.eq.6) then
        do i = 1,ntypes
          if (node.eq.0) read (2,*) jtmp,(tmp(k),k=1,4)
          call mpi_bcast(tmp,4,mpi_double_precision,0,
     $         mpi_comm_world,ierror)
          noncoeff1(i,i) = tmp(1)
          noncoeff2(i,i) = tmp(2)
          noncoeff14_1(i,i) = tmp(3)
          noncoeff14_2(i,i) = tmp(4)
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c read bond coeffs from data file

      subroutine read_bond_coeff(nbondcoeffflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nbondcoeffflag

c local variables

      integer i,j,jtmp,isize,ierror

      if (node.eq.0) write (6,*) 'Bond Coeffs ...'
      if (nbonds.eq.0.or.bondstyle.eq.0)
     $     call err('Bond Coeffs entry not allowed')
      nbondcoeffflag = 1

      if (bondstyle.eq.1) then
        if (node.eq.0) then
          do i = 1,nbondtypes
            read (2,*) jtmp,(bondcoeff(j,i),j=1,2)
          enddo
        endif
      else if (bondstyle.eq.2) then
        if (node.eq.0) then
          do i = 1,nbondtypes
            read (2,*) jtmp,(bondcoeff(j,i),j=1,4)
          enddo
        endif
      else if (bondstyle.eq.3) then
        if (node.eq.0) then
          do i = 1,nbondtypes
            read (2,*) jtmp,(bondcoeff(j,i),j=1,5)
          enddo
        endif
      else if (bondstyle.eq.4) then
        if (node.eq.0) then
          do i = 1,nbondtypes
            read (2,*) jtmp,(bondcoeff(j,i),j=1,3)
          enddo
        endif
      else if (bondstyle.eq.5) then
        if (node.eq.0) then
          do i = 1,nbondtypes
            read (2,*) jtmp,(bondcoeff(j,i),j=1,4)
          enddo
        endif
      endif
      
      isize = nbondtypes*5
      call mpi_bcast(bondcoeff,isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      do i = 1,nbondtypes
        bondtypeflag(i) = bondstyle
      enddo

      return
      end


c -------------------------------------------------------------------------
c read angle coeffs from data file

      subroutine read_angle_coeff(nanglecoeffflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nanglecoeffflag

c local variables

      integer i,j,ierror,isize,jtmp

      if (node.eq.0) write (6,*) 'Angle Coeffs ...'
      if (nangles.eq.0.or.anglestyle.eq.0)
     $     call err('Angle Coeffs entry not allowed')
      nanglecoeffflag = 1

      if (anglestyle.eq.1) then
        if (node.eq.0) then
          do i = 1,nangletypes
            read (2,*) jtmp,(anglecoeff(j,i),j=1,2)
          enddo
        endif
      else if (anglestyle.eq.2) then
        if (node.eq.0) then
          do i = 1,nangletypes
            read (2,*) jtmp,(anglecoeff(j,i),j=1,4)
          enddo
        endif
      else if (anglestyle.eq.3) then
        if (node.eq.0) then
          do i = 1,nangletypes
            read (2,*) jtmp,(anglecoeff(j,i),j=1,4)
          enddo
        endif
      else if (anglestyle.eq.4) then
        if (node.eq.0) then
          do i = 1,nangletypes
            read (2,*) jtmp,(anglecoeff(j,i),j=1,1)
          enddo
        endif
      endif

      isize = nangletypes*4
      call mpi_bcast(anglecoeff,isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      do i = 1,nangletypes
        angletypeflag(i) = anglestyle
      enddo

      return
      end


c -------------------------------------------------------------------------
c read dihedral coeffs from data file

      subroutine read_dihedral_coeff(ndihedcoeffflag)
      use global
      use mpi
      implicit none

c argument variables

      integer ndihedcoeffflag

c local variables

      integer i,isize,ierror,j,jtmp

      if (node.eq.0) write (6,*) 'Dihedral Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.eq.0)
     $     call err('Dihedral Coeffs entry not allowed')
      ndihedcoeffflag = 1

      if (dihedstyle.eq.1) then
        if (node.eq.0) then
          do i = 1,ndihedtypes
            read (2,*) jtmp,(dihedcoeff(j,i),j=1,3)
          enddo
        endif
      else if (dihedstyle.eq.2) then
        if (node.eq.0) then
          do i = 1,ndihedtypes
            read (2,*) jtmp,(dihedcoeff(j,i),j=1,6)
          enddo
        endif
      else if (dihedstyle.eq.3) then
        if (node.eq.0) then
          do i = 1,ndihedtypes
            read (2,*) jtmp,(dihedcoeff(j,i),j=1,5)
          enddo
        endif
      else if (dihedstyle.eq.4) then
        if (node.eq.0) then
          do i = 1,ndihedtypes
            read (2,*) jtmp,(dihedcoeff(j,i),j=1,4)
          enddo
        endif
      endif

      isize = ndihedtypes*6
      call mpi_bcast(dihedcoeff,isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      do i = 1,ndihedtypes
        dihedtypeflag(i) = dihedstyle
      enddo

      return
      end


c -------------------------------------------------------------------------
c read improper coeffs from data file

      subroutine read_improper_coeff(nimprocoeffflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nimprocoeffflag

c local variables

      integer i,isize,ierror,j,jtmp

      if (node.eq.0) write (6,*) 'Improper Coeffs ...'
      if (nimpropers.eq.0.or.improstyle.eq.0)
     $     call err('Improper Coeffs entry not allowed')
      nimprocoeffflag = 1

      if (improstyle.eq.1) then
        if (node.eq.0) then
          do i = 1,nimprotypes
            read (2,*) jtmp,(improcoeff(j,i),j=1,2)
          enddo
        endif
      else if (improstyle.eq.2) then
        if (node.eq.0) then
          do i = 1,nimprotypes
            read (2,*) jtmp,(improcoeff(j,i),j=1,3)
          enddo
        endif
      else if (improstyle.eq.3) then
        if (node.eq.0) then
          do i = 1,nimprotypes
            read (2,*) jtmp,(improcoeff(j,i),j=1,2)
          enddo
        endif
      endif

      isize = nimprotypes*3
      call mpi_bcast(improcoeff,isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      do i = 1,nimprotypes
        improtypeflag(i) = improstyle
      enddo

      return
      end


c -------------------------------------------------------------------------
c read bondbond coeffs from data file

      subroutine read_bondbond_coeff(nbondbondflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nbondbondflag

c local variables

      integer i,isize,ierror,j,jtmp

      if (node.eq.0) write (6,*) 'BondBond Coeffs ...'
      if (nangles.eq.0.or.anglestyle.ne.2)
     $     call err('BondBond Coeffs entry not allowed')
      nbondbondflag = 1

      if (node.eq.0) then
        do i = 1,nangletypes
          read (2,*) jtmp,(bondbondcoeff(j,i),j=1,3)
        enddo
      endif

      isize = nangletypes*3
      call mpi_bcast(bondbondcoeff,isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read bondangle coeffs from data file

      subroutine read_bondangle_coeff(nbondangleflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nbondangleflag

c local variables

      integer i,isize,ierror,j,jtmp

      if (node.eq.0) write (6,*) 'BondAngle Coeffs ...'
      if (nangles.eq.0.or.anglestyle.ne.2)
     $     call err('BondAngle Coeffs entry not allowed')
      nbondangleflag = 1

      if (node.eq.0) then
        do i = 1,nangletypes
          read (2,*) jtmp,(bondanglecoeff(j,i),j=1,4)
        enddo
      endif

      isize = nangletypes*4
      call mpi_bcast(bondanglecoeff,isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read midbondtorsion coeffs from data file

      subroutine read_midbondtorsion_coeff(nmidbondtorsionflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nmidbondtorsionflag

c local variables

      integer i,isize,ierror,j,jtmp

      if (node.eq.0) write (6,*) 'MiddleBondTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) call err(
     $     'MiddleBondTorsion Coeffs entry not allowed')
      nmidbondtorsionflag = 1

      if (node.eq.0) then
        do i = 1,ndihedtypes
          read (2,*) jtmp,(midbondtorsioncoeff(j,i),j=1,4)
        enddo
      endif

      isize = ndihedtypes*4
      call mpi_bcast(midbondtorsioncoeff,
     $     isize,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read endbondtorsion coeffs from data file

      subroutine read_endbondtorsion_coeff(nendbondtorsionflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nendbondtorsionflag

c local variables

      integer i,isize,ierror,jtmp,j

      if (node.eq.0) write (6,*) 'EndBondTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) call err(
     $     'EndBondTorsion Coeffs entry not allowed')
      nendbondtorsionflag = 1

      if (node.eq.0) then
        do i = 1,ndihedtypes
          read (2,*) jtmp,(endbondtorsioncoeff(j,i),j=1,8)
        enddo
      endif

      isize = ndihedtypes*8
      call mpi_bcast(endbondtorsioncoeff,
     $     isize,mpi_double_precision,0,mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read angletorsion coeffs from data file

      subroutine read_angletorsion_coeff(nangletorsionflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nangletorsionflag

c local variables

      integer i,isize,ierror,jtmp,j

      if (node.eq.0) write (6,*) 'AngleTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) call err(
     $     'AngleTorsion Coeffs entry not allowed')
      nangletorsionflag = 1

      if (node.eq.0) then
        do i = 1,ndihedtypes
          read (2,*) jtmp,(angletorsioncoeff(j,i),j=1,8)
        enddo
      endif

      isize = ndihedtypes*8
      call mpi_bcast(angletorsioncoeff,
     $     isize,mpi_double_precision,0,mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read angleangletorsion coeffs from data file

      subroutine read_angleangletorsion_coeff(nangleangletorsionflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nangleangletorsionflag

c local variables

      integer i,isize,ierror,jtmp,j

      if (node.eq.0) write (6,*) 'AngleAngleTorsion Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) call err(
     $     'AngleAngleTorsion Coeffs entry not allowed')
      nangleangletorsionflag = 1

      if (node.eq.0) then
        do i = 1,ndihedtypes
          read (2,*) jtmp,(angleangletorsioncoeff(j,i),j=1,3)
        enddo
      endif

      isize = ndihedtypes*3
      call mpi_bcast(angleangletorsioncoeff,
     $     isize,mpi_double_precision,0,mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read bondbond13 coeffs from data file

      subroutine read_bondbond13_coeff(nbondbond13flag)
      use global
      use mpi
      implicit none

c argument variables

      integer nbondbond13flag

c local variables

      integer i,isize,ierror,jtmp,j

      if (node.eq.0) write (6,*) 'BondBond13 Coeffs ...'
      if (ndihedrals.eq.0.or.dihedstyle.ne.2) call err(
     $     'BondBond13 Coeffs entry not allowed')
      nbondbond13flag = 1

      if (node.eq.0) then
        do i = 1,ndihedtypes
          read (2,*) jtmp,(bondbond13coeff(j,i),j=1,3)
        enddo
      endif

      isize = ndihedtypes*3
      call mpi_bcast(bondbond13coeff,
     $     isize,mpi_double_precision,0,mpi_comm_world,ierror)

      return
      end


c -------------------------------------------------------------------------
c read masses from data file

      subroutine read_angleangle_coeff(nangleangleflag)
      use global
      use mpi
      implicit none

c argument variables

      integer nangleangleflag

c local variables

      integer i,isize,ierror,jtmp,j

      if (node.eq.0) write (6,*) 'AngleAngle Coeffs ...'
      if (nimpropers.eq.0.or.improstyle.ne.3) call err(
     $     'AngleAngle Coeffs entry not allowed')
      nangleangleflag = 1

      if (node.eq.0) then
        do i = 1,nimprotypes
          read (2,*) jtmp,(angleanglecoeff(j,i),j=1,6)
        enddo
      endif

      isize = nimprotypes*6
      call mpi_bcast(angleanglecoeff,
     $     isize,mpi_double_precision,0,mpi_comm_world,ierror)

      return
      end
