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
c scan an input data file to compute max # of bonded terms per atom
c only called by proc 0

      subroutine scan_data(datafile,idimension,newton_bond,
     $     maxbondper,maxangleper,maxdihedper,maximproper)
      implicit none

c argument variables

      character*(*) datafile
      integer idimension,newton_bond
      integer maxbondper,maxangleper,maxdihedper,maximproper

c local variables

      logical match
      integer i,m,atom1,atom2,atom3,atom4,itmp1,itmp2
      integer natoms,nbonds,nangles,ndihedrals,nimpropers
      integer nbondflag,nangleflag,ndihedflag,nimproflag
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
      integer ierror
      real*8 tmp1,tmp2
      character*80 str
      integer, allocatable :: count(:)

 900  format (a)

c initialize the returned values

      maxbondper = 0
      maxangleper = 0
      maxdihedper = 0
      maximproper = 0

c scan the data file

      write (6,*) 'Scanning data file ...'

      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0

c read data file header - fixed format

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
      read (2,*) tmp1,tmp2
      read (2,*) tmp1,tmp2
      if (idimension.eq.3) read (2,*) tmp1,tmp2

c error check

      if (natoms.le.0) call abort('Atoms <= 0')
      if (nbonds.lt.0) call abort('Bonds < 0')
      if (nangles.lt.0) call abort('Angles < 0')
      if (ndihedrals.lt.0) call abort('Dihedrals < 0')
      if (nimpropers.lt.0) call abort('Impropers < 0')

      if (ntypes.le.0) call abort('Atom types <= 0')
      if (nbonds.gt.0.and.nbondtypes.le.0)
     $     call abort('Bond types <= 0')
      if (nangles.gt.0.and.nangletypes.le.0)
     $     call abort('Angle types <= 0')
      if (ndihedrals.gt.0.and.ndihedtypes.le.0)
     $     call abort('Dihedral types <= 0')
      if (nimpropers.gt.0.and.nimprotypes.le.0)
     $     call abort('Improper types <= 0')

c read remainder of data file - free format

      nbondflag = 0
      nangleflag = 0
      ndihedflag = 0
      nimproflag = 0

      do

c read identifier string

        read (2,*,end=999,err=999)
        read (2,900,end=999,err=999) str
        read (2,*,end=999,err=999)

c read past atom masses

        if (match('Masses',str,m)) then

          do i = 1,ntypes
            read (2,*)
          enddo

c read past atoms

        else if (match('Atoms',str,m)) then

          write (6,*) 'Scanning atoms ...'
          do i = 1,natoms
            read (2,*)
          enddo

c read past atom velocities

        else if (match('Velocities',str,m)) then

          write (6,*) 'Scanning velocities ...'
          do i = 1,natoms
            read (2,*)
          enddo

c read bonds and count how many per atom
c depends on newton_bond flag

        else if (match('Bonds',str,m)) then

          write (6,*) 'Scanning bonds ...'
          if (nbonds.eq.0.or.nbondtypes.eq.0)
     $         call abort('Should be no Bonds entry')
          nbondflag = 1

          allocate(count(natoms))
          count = 0

          do i = 1,nbonds
            read (2,*) itmp1,itmp2,atom1,atom2
            count(atom1) = count(atom1) + 1
            if (newton_bond == 0) count(atom2) = count(atom2) + 1
          enddo

          maxbondper = 0
          do i = 1,natoms
            maxbondper = max(maxbondper,count(i))
          enddo
          write (6,*) 'Max # bonds/atom =',maxbondper
          write (1,*) 'Max # bonds/atom =',maxbondper

          deallocate(count)

c read angles and count how many per atom
c depends on newton_bond flag

        else if (match('Angles',str,m)) then

          write (6,*) 'Scanning angles ...'
          if (nangles.eq.0.or.nangletypes.eq.0)
     $         call abort('Should be no Angles entry')
          nangleflag = 1

          allocate(count(natoms))
          count = 0

          do i = 1,nangles
            read (2,*) itmp1,itmp2,atom1,atom2,atom3
            count(atom2) = count(atom2) + 1
            if (newton_bond == 0) then
              count(atom1) = count(atom1) + 1
              count(atom3) = count(atom3) + 1
            endif
          enddo

          maxangleper = 0
          do i = 1,natoms
            maxangleper = max(maxangleper,count(i))
          enddo
          write (6,*) 'Max # angles/atom =',maxangleper
          write (1,*) 'Max # angles/atom =',maxangleper

          deallocate(count)

c read dihedrals and count how many per atom
c depends on newton_bond flag

        else if (match('Dihedrals',str,m)) then

          write (6,*) 'Scanning dihedrals ...'
          if (ndihedrals.eq.0.or.ndihedtypes.eq.0)
     $         call abort('Should be no Dihedrals entry')
          ndihedflag = 1

          allocate(count(natoms))
          count = 0

          do i = 1,ndihedrals
            read (2,*) itmp1,itmp2,atom1,atom2,atom3,atom4
            count(atom2) = count(atom2) + 1
            if (newton_bond == 0) then
              count(atom1) = count(atom1) + 1
              count(atom3) = count(atom3) + 1
              count(atom4) = count(atom4) + 1
            endif
          enddo

          maxdihedper = 0
          do i = 1,natoms
            maxdihedper = max(maxdihedper,count(i))
          enddo
          write (6,*) 'Max # dihedrals/atom =',maxdihedper
          write (1,*) 'Max # dihedrals/atom =',maxdihedper

          deallocate(count)

c read impropers and count how many per atom
c depends on newton_bond flag

        else if (match('Impropers',str,m)) then

          write (6,*) 'Scanning impropers ...'
          if (nimpropers.eq.0.or.nimprotypes.eq.0)
     $         call abort('Should be no Impropers entry')
          nimproflag = 1

          allocate(count(natoms))
          count = 0

          do i = 1,nimpropers
            read (2,*) itmp1,itmp2,atom1,atom2,atom3,atom4
            count(atom2) = count(atom2) + 1
            if (newton_bond == 0) then
              count(atom1) = count(atom1) + 1
              count(atom3) = count(atom3) + 1
              count(atom4) = count(atom4) + 1
            endif
          enddo

          maximproper = 0
          do i = 1,natoms
            maximproper = max(maximproper,count(i))
          enddo
          write (6,*) 'Max # impropers/atom =',maximproper
          write (1,*) 'Max # impropers/atom =',maximproper

          deallocate(count)

c read past nonbond coeffs

        else if (match('Nonbond Coeffs',str,m)) then

          do i = 1,ntypes
            read (2,*)
          enddo

c read past bond coeffs

        else if (match('Bond Coeffs',str,m)) then

          do i = 1,nbondtypes
            read (2,*)
          enddo

c read past angle coeffs

        else if (match('Angle Coeffs',str,m)) then

          do i = 1,nangletypes
            read (2,*)
          enddo

c read past dihedral coeffs

        else if (match('Dihedral Coeffs',str,m)) then

          do i = 1,ndihedtypes
            read (2,*)
          enddo

c read past improper coeffs

        else if (match('Improper Coeffs',str,m)) then

          do i = 1,nimprotypes
            read (2,*)
          enddo

c read past bond-bond coeffs for class 2 force field

        else if (match('BondBond Coeffs',str,m)) then

          do i = 1,nangletypes
            read (2,*)
          enddo

c read past bond-angle coeffs for class 2 force field

        else if (match('BondAngle Coeffs',str,m)) then

          do i = 1,nangletypes
            read (2,*)
          enddo

c read past class 2 middle-bond-torsion coeffs

        else if (match('MiddleBondTorsion Coeffs',str,m)) then

          do i = 1,ndihedtypes
            read (2,*)
          enddo

c read past class 2 end-bond-torsion coeffs

        else if (match('EndBondTorsion Coeffs',str,m)) then

          do i = 1,ndihedtypes
            read (2,*)
          enddo

c read past class 2 angle-torsion coeffs

        else if (match('AngleTorsion Coeffs',str,m)) then

          do i = 1,ndihedtypes
            read (2,*)
          enddo

c read past class 2 angle-angle-torsion coeffs

        else if (match('AngleAngleTorsion Coeffs',str,m)) then

          do i = 1,ndihedtypes
            read (2,*)
          enddo

c read past class 2 bondbond13 coeffs

        else if (match('BondBond13 Coeffs',str,m)) then

          do i = 1,ndihedtypes
            read (2,*)
          enddo

c read past class 2 angle-angle coeffs

        else if (match('AngleAngle Coeffs',str,m)) then

          do i = 1,nimprotypes
            read (2,*)
          enddo

c unknown identifier

        else
          write (6,*) 'UNKNOWN: ',trim(str)
          call abort('Unknown identifier in data file')
        endif

      enddo

c close data file

 999  continue

      close (2)

c check that all required options were read in from data file

      if (nbonds.gt.0.and.nbondflag.eq.0)
     $     call abort('No Bonds in data file')
      if (nangles.gt.0.and.nangleflag.eq.0)
     $     call abort('No Angles in data file')
      if (ndihedrals.gt.0.and.ndihedflag.eq.0)
     $     call abort('No Dihedrals in data file')
      if (nimpropers.gt.0.and.nimproflag.eq.0)
     $     call abort('No Impropers in data file')

      return
      end
