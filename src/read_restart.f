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
c read a restart file

      subroutine read_restart
      use global
      use mpi
      implicit none

c local variables

      integer, parameter :: iperatom = 16
      integer ierror,iversion_tmp,nbyte,ibyte
      integer nprocs_tmp,pgridx,pgridy,pgridz,idimension_tmp
      integer perflagx_tmp,perflagy_tmp,perflagz_tmp
      integer units_tmp,newton_tmp,i,j,k,jj,n
      integer numbond_tmp,numangle_tmp,numdihed_tmp,numimpro_tmp
      integer isize,jtmp,itag,iflagmax,natoms_tmp
      integer nleft,nmin
      real*8 tmp
      real*8, allocatable :: buf(:)
      
c read the restart file

      call mpi_barrier(mpi_comm_world,ierror)
      if (node.eq.0) write (6,*) 'Reading restart file ...'
      
      readflag = 1

      if (node == 0) then
        open (unit=2,file=restart_in,
     $       status='old',form='unformatted',iostat=ierror)
        if (ierror /= 0) call abort('No restart file')
      endif

c read 1st line of restart file and check version #

      if (restart_version == 5 .or. restart_version == 6) then
        if (node == 0) read (2) iversion_tmp,nbyte,ibyte
      else
        if (node == 0) read (2) iversion_tmp
      endif

      call mpi_bcast(iversion_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      if (iversion_tmp /= restart_version) then
        if (node == 0) write (6,*)
     $       'Use restart version command, set to',iversion_tmp
        call err('Version mis-match in restart file')
      endif

      if (iversion_tmp /= 2001 .and. iversion_tmp /= 2000
     $     .and. iversion_tmp /= 6 .and. iversion_tmp /= 5)
     $     call err('Too old a restart file -- cannot read')

c read remainder of restart file header

      nbondtypes = 0
      nangletypes = 0
      ndihedtypes = 0
      nimprotypes = 0

      if (node.eq.0) then
        read (2) ntimestep,nprocs_tmp,pgridx,pgridy,pgridz
        read (2) idimension_tmp,perflagx_tmp,perflagy_tmp,perflagz_tmp
        read (2) units_tmp,newton_tmp
        read (2) natoms,nbonds,nangles,ndihedrals,nimpropers
        read (2) ntypes
        if (nbonds > 0) read (2) nbondtypes
        if (nangles > 0) read (2) nangletypes
        if (ndihedrals > 0) read (2) ndihedtypes
        if (nimpropers > 0) read (2) nimprotypes
        read (2) maxbondper,maxangleper,maxdihedper,maximproper
        read (2) box
        if (restart_version == 2000 .or. restart_version == 2001) then
          read (2) eta,eta_dot
          read (2) omega(1),omega(2),omega(3),
     $         omega_dot(1),omega_dot(2),omega_dot(3)
        else
          read (2) eta,eta_dot,omega_dot(1),omega_dot(2),omega_dot(3)
          omega(1) = 0.0
          omega(2) = 0.0
          omega(3) = 0.0
        endif
      endif

c broadcast info

      call mpi_bcast(ntimestep,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(nprocs_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(pgridx,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(pgridy,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(pgridz,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(idimension_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(perflagx_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(perflagy_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(perflagz_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(units_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(newton_tmp,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

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

      call mpi_bcast(maxbondper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(maxangleper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(maxdihedper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(maximproper,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(box,6,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(eta,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(eta_dot,1,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(omega,3,mpi_double_precision,0,
     $     mpi_comm_world,ierror)
      call mpi_bcast(omega_dot,3,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

c error check
      
      if (nprocs_tmp.ne.nprocs) then
        if (node.eq.0) write (6,*) 'WARNING: Restarting on ',
     $       'different number of processors'
      endif

      if (idimension_tmp.ne.idimension) then
        if (node.eq.0) write (6,*) 'WARNING: Old dimension ',
     $       'does not match current dimension'
      endif

      if (perflagx_tmp.ne.perflagx.or.perflagy_tmp.ne.perflagy.or.
     $     perflagz_tmp.ne.perflagz) then
        if (node.eq.0) write (6,*) 'WARNING: Old periodicity ',
     $       'does not match current periodicity'
      endif

      if (units_tmp.ne.units) then
        if (node.eq.0) write (6,*) 'WARNING: Old units ',
     $       'does not match current units'
      endif

      if (newton_tmp.ne.newton) then
        if (node.eq.0) write (6,*) 'WARNING: Old newton flag ',
     $       'does not match current newton flag'
        if (mod(newton_tmp,2) /= mod(newton,2) .and.
     $       (nbonds > 0 .or. nangles > 0 .or.
     $       ndihedrals > 0 .or. nimpropers > 0))
     $       call err(
     $       'Cannot change bonded newton flag for bonded systems')
      endif

      if (idimension.eq.2.and.pgrid(1).gt.0) then
        if (pgrid(3).ne.1)
     $       call err('Can only use one z-processor for 2-d')
      endif

c define simulation box
c for non-PBC, assume bounds were already augmented by "small" when
c  restart file was written

      xprd = box(2,1) - box(1,1)
      yprd = box(2,2) - box(1,2)
      zprd = box(2,3) - box(1,3)

      xprd_half = 0.5 * xprd
      yprd_half = 0.5 * yprd
      zprd_half = 0.5 * zprd

c define mesh topology and error check against restart file values

      call mesh_3d(node,nprocs,idimension,xprd,yprd,zprd,
     $     pgrid,me,mpart)

      if (node.eq.0) then
        write (6,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
        write (1,*) '3-d grid of procs =',pgrid(1),pgrid(2),pgrid(3)
      endif

      if (pgridx.ne.pgrid(1).or.pgridy.ne.pgrid(2).or.
     $     pgridz.ne.pgrid(3)) then
        if (node.eq.0) write (6,*) 'WARNING: Restarting on ',
     $       'different 3-d grid of processors'
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

c read all coefficients

      if (node.eq.0) then

        write (6,*) 'Coefficients ...'

        read (2) (mass(i),i=1,ntypes)

        read (2) nonstyle
        if (nonstyle.eq.0) then
          write (6,*) 'Nonbond style set to none'
        else if (nonstyle.eq.1) then
          write (6,*) 'Nonbond style set to lj/cutoff'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.2) then
          write (6,*) 'Nonbond style set to lj/smooth'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.3) then
          write (6,*) 'Nonbond style set to lj/shift'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.4) then
          write (6,*) 'Nonbond style set to soft'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.5) then
          write (6,*) 'Nonbond style set to class2/cutoff'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
        else if (nonstyle.eq.6) then
          write (6,*) 'Nonbond style set to lj/charmm'
          read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff14_1(i,j),j=i,ntypes),i=1,ntypes),
     $         ((noncoeff14_2(i,j),j=i,ntypes),i=1,ntypes)
        endif

        if (nbonds.gt.0) then
          read (2) bondstyle
          if (bondstyle.eq.0) then
            write (6,*) 'Bond style set to none'
          else if (bondstyle.eq.1) then
            write (6,*) 'Bond style set to harmonic'
            read (2) ((bondcoeff(i,j),i=1,2),j=1,nbondtypes)
          else if (bondstyle.eq.2) then
            write (6,*) 'Bond style set to fene/standard'
            read (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          else if (bondstyle.eq.3) then
            write (6,*) 'Bond style set to fene/shift'
            read (2) ((bondcoeff(i,j),i=1,5),j=1,nbondtypes)
          else if (bondstyle.eq.4) then
            write (6,*) 'Bond style set to nonlinear'
            read (2) ((bondcoeff(i,j),i=1,3),j=1,nbondtypes)
          else if (bondstyle.eq.5) then
            write (6,*) 'Bond style set to class2'
            read (2) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
          endif
        endif

        if (nangles.gt.0) then
          read (2) anglestyle
          if (anglestyle.eq.0) then
            write (6,*) 'Angle style set to none'
          else if (anglestyle.eq.1) then
            write (6,*) 'Angle style set to harmonic'
            read (2) ((anglecoeff(i,j),i=1,2),j=1,nangletypes)
          else if (anglestyle.eq.2) then
            write (6,*) 'Angle style set to class2'
            read (2) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
            read (2) ((bondbondcoeff(i,j),i=1,3),j=1,nangletypes)
            read (2) ((bondanglecoeff(i,j),i=1,4),j=1,nangletypes)
          else if (anglestyle.eq.3) then
            write (6,*) 'Angle style set to charmm'
            read (2) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
          else if (anglestyle.eq.4) then
            write (6,*) 'Angle style set to cosine'
            read (2) ((anglecoeff(i,j),i=1,1),j=1,nangletypes)
          endif
        endif

        if (ndihedrals.gt.0) then
          read (2) dihedstyle
          if (dihedstyle.eq.0) then
            write (6,*) 'Dihedral style set to none'
          else if (dihedstyle.eq.1) then
            write (6,*) 'Dihedral style set to harmonic'
            read (2) ((dihedcoeff(i,j),i=1,3),j=1,ndihedtypes)
          else if (dihedstyle.eq.2) then
            write (6,*) 'Dihedral style set to class2'
            read (2) ((dihedcoeff(i,j),i=1,6),j=1,ndihedtypes)
            read (2) ((midbondtorsioncoeff(i,j),i=1,4),
     $           j=1,ndihedtypes)
            read (2) ((endbondtorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            read (2) ((angletorsioncoeff(i,j),i=1,8),
     $           j=1,ndihedtypes)
            read (2) ((angleangletorsioncoeff(i,j),i=1,3),
     $           j=1,ndihedtypes)
            if (restart_version == 2000 .or. restart_version == 2001)
     $           read (2) ((bondbond13coeff(i,j),i=1,3),
     $           j=1,ndihedtypes)
          else if (dihedstyle.eq.3) then
            write (6,*) 'Dihedral style set to multiharmonic'
            read (2) ((dihedcoeff(i,j),i=1,5),j=1,ndihedtypes)
          else if (dihedstyle.eq.4) then
            write (6,*) 'Dihedral style set to charmm'
            read (2) ((dihedcoeff(i,j),i=1,4),j=1,ndihedtypes)
          endif
        endif

        if (nimpropers.gt.0) then
          read (2) improstyle
          if (improstyle.eq.0) then
            write (6,*) 'Improper style set to none'
          else if (improstyle.eq.1) then
            write (6,*) 'Improper style set to harmonic'
            read (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
          else if (improstyle.eq.2) then
            write (6,*) 'Improper style set to cvff'
            read (2) ((improcoeff(i,j),i=1,3),j=1,nimprotypes)
          else if (improstyle.eq.3) then
            write (6,*) 'Improper style set to class2'
            read (2) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
            read (2) ((angleanglecoeff(i,j),i=1,6),j=1,nimprotypes)
          endif
        endif

      endif

c broadcast all coefficients

      call mpi_bcast(mass,ntypes,mpi_double_precision,0,
     $     mpi_comm_world,ierror)

      call mpi_bcast(nonstyle,1,mpi_integer,0,
     $     mpi_comm_world,ierror)

      if (nonstyle.eq.1) then
        call mpi_bcast(noncoeff1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.2) then
        call mpi_bcast(noncoeff1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff4,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.3) then
        call mpi_bcast(noncoeff1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff4,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.4) then
        call mpi_bcast(noncoeff1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.5) then
        call mpi_bcast(noncoeff1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      else if (nonstyle.eq.6) then
        call mpi_bcast(noncoeff1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff3,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff14_1,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_bcast(noncoeff14_2,ntypes*ntypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
      endif

      do i = 1,ntypes
        do j = i,ntypes
          nontypeflag(i,j) = nonstyle
        enddo
      enddo

      if (nbonds.gt.0) then
        call mpi_bcast(bondstyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (bondstyle.gt.0)
     $       call mpi_bcast(bondcoeff,5*nbondtypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        do i = 1,nbondtypes
          bondtypeflag(i) = bondstyle
        enddo
      endif

      if (nangles.gt.0) then
        call mpi_bcast(anglestyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (anglestyle.gt.0)
     $       call mpi_bcast(anglecoeff,4*nangletypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        if (anglestyle.eq.2) then
          call mpi_bcast(bondbondcoeff,3*nangletypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(bondanglecoeff,4*nangletypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
        endif
        do i = 1,nangletypes
          angletypeflag(i) = anglestyle
        enddo
      endif

      if (ndihedrals.gt.0) then
        call mpi_bcast(dihedstyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (dihedstyle.gt.0)
     $       call mpi_bcast(dihedcoeff,6*ndihedtypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        if (dihedstyle.eq.2) then
          call mpi_bcast(midbondtorsioncoeff,4*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(endbondtorsioncoeff,8*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(angletorsioncoeff,8*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(angleangletorsioncoeff,3*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
          call mpi_bcast(bondbond13coeff,3*ndihedtypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
        endif
        do i = 1,ndihedtypes
          dihedtypeflag(i) = dihedstyle
        enddo
      endif

      if (nimpropers.gt.0) then
        call mpi_bcast(improstyle,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        if (improstyle.gt.0)
     $       call mpi_bcast(improcoeff,3*nimprotypes,
     $       mpi_double_precision,0,mpi_comm_world,ierror)
        if (improstyle.eq.3) then
          call mpi_bcast(angleanglecoeff,6*nimprotypes,
     $         mpi_double_precision,0,mpi_comm_world,ierror)
        endif
        do i = 1,nimprotypes
          improtypeflag(i) = improstyle
        enddo
      endif

c allocate a temporary buffer for the restart info
c big enough for many atom's info
c nmin = max size one atom's info may require in buffer

      allocate(buf(10000))

      if (node.eq.0) write (6,*) 'Reading',natoms,' atoms ...'

      nleft = natoms
      nmin = 1000

c while some atoms are still left to read

      do while (nleft > 0)

c proc 0 reads info for as many atoms as will fit into buffer

        if (node == 0) then

          j = 0
          do while (nleft > 0 .and. j < 10000-nmin)
            nleft = nleft - 1
            read (2) (buf(j+k),k=1,iperatom)
            j = j + iperatom
            numbond_tmp = nint(buf(j-3))
            numangle_tmp = nint(buf(j-2))
            numdihed_tmp = nint(buf(j-1))
            numimpro_tmp = nint(buf(j))
            do jj = 1,numbond_tmp
              read (2) buf(j+1),buf(j+2),buf(j+3)
              j = j + 3
            enddo
            do jj = 1,numangle_tmp
              read (2) buf(j+1),buf(j+2),buf(j+3),buf(j+4)
              j = j + 4
            enddo
            do jj = 1,numdihed_tmp
              read (2) buf(j+1),buf(j+2),buf(j+3),
     $             buf(j+4),buf(j+5)
              j = j + 5
            enddo
            do jj = 1,numimpro_tmp
              read (2) buf(j+1),buf(j+2),buf(j+3),
     $             buf(j+4),buf(j+5)
              j = j + 5
            enddo
          enddo
          isize = j

        endif

c broadcast the buffer of atom info to all procs along with its size

        call mpi_bcast(nleft,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        call mpi_bcast(isize,1,mpi_integer,0,
     $       mpi_comm_world,ierror)
        call mpi_bcast(buf,isize,mpi_double_precision,0,
     $       mpi_comm_world,ierror)

c loop over all atoms in the buffer
c add atom to my list if in my sub-box
c do not check for PBC or remap since assume restart file is "correct"

        j = 0
        do while (j < isize)

          if (buf(j+1) >= border(1,1) .and.
     $         buf(j+1) < border(2,1) .and.
     $         buf(j+2) >= border(1,2) .and.
     $         buf(j+2) < border(2,2) .and.
     $         buf(j+3) >= border(1,3) .and.
     $         buf(j+3) < border(2,3)) then
            nlocal = nlocal + 1
            if (nlocal <= maxown) then
              x(1,nlocal) = buf(j+1)
              x(2,nlocal) = buf(j+2)
              x(3,nlocal) = buf(j+3)
              itag = nint(buf(j+4))
              localptr(itag) = nlocal
              tag(nlocal) = itag
              molecule(nlocal) = nint(buf(j+5))
              type(nlocal) = nint(buf(j+6))
              q(nlocal) = buf(j+7)
              v(1,nlocal) = buf(j+8)
              v(2,nlocal) = buf(j+9)
              v(3,nlocal) = buf(j+10)
              true(nlocal) = nint(buf(j+11))
              fix(nlocal) = nint(buf(j+12))
              numbond(nlocal) = nint(buf(j+13))
              numangle(nlocal) = nint(buf(j+14))
              numdihed(nlocal) = nint(buf(j+15))
              numimpro(nlocal) = nint(buf(j+16))
              velflag(nlocal) = 1
              n = j + iperatom
              do jj = 1,numbond(nlocal)
                bondatom1(jj,nlocal) = nint(buf(n+1))
                bondatom2(jj,nlocal) = nint(buf(n+2))
                bondtype(jj,nlocal) = nint(buf(n+3))
                n = n + 3
              enddo
              do jj = 1,numangle(nlocal)
                angleatom1(jj,nlocal) = nint(buf(n+1))
                angleatom2(jj,nlocal) = nint(buf(n+2))
                angleatom3(jj,nlocal) = nint(buf(n+3))
                angletype(jj,nlocal) = nint(buf(n+4))
                n = n + 4
              enddo
              do jj = 1,numdihed(nlocal)
                dihedatom1(jj,nlocal) = nint(buf(n+1))
                dihedatom2(jj,nlocal) = nint(buf(n+2))
                dihedatom3(jj,nlocal) = nint(buf(n+3))
                dihedatom4(jj,nlocal) = nint(buf(n+4))
                dihedtype(jj,nlocal) = nint(buf(n+5))
                n = n + 5
              enddo
              do jj = 1,numimpro(nlocal)
                improatom1(jj,nlocal) = nint(buf(n+1))
                improatom2(jj,nlocal) = nint(buf(n+2))
                improatom3(jj,nlocal) = nint(buf(n+3))
                improatom4(jj,nlocal) = nint(buf(n+4))
                improtype(jj,nlocal) = nint(buf(n+5))
                n = n + 5
              enddo
            endif
          endif

c skip over one atom in buffer

          j = j + iperatom +
     $         3*nint(buf(j+iperatom-3)) +
     $         4*nint(buf(j+iperatom-2)) +
     $         5*nint(buf(j+iperatom-1)) +
     $         5*nint(buf(j+iperatom))

        enddo

      enddo

c free temporary buffer

      deallocate(buf)

c close restart file

      if (node == 0) close (2)

c check if too many atoms on any proc

      call mpi_allreduce(nlocal,iflagmax,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (iflagmax > maxown) then
        if (node.eq.0) write (6,*) 'Nlocal,maxown =',iflagmax,maxown
        call err('Too many atoms/processor - boost extra_own')
      endif

c check that all atoms were assigned uniquely to a processor

      call mpi_allreduce(nlocal,natoms_tmp,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (natoms /= natoms_tmp)
     $     call err('Did not assign all atoms correctly')

c compute total mass and charge

      masssum = 0.0
      qsum = 0.0
      qsqsum = 0.0

      do i = 1,nlocal
        masssum = masssum + mass(type(i))
        qsum = qsum + q(i)
        qsqsum = qsqsum + q(i)*q(i)
      enddo

      tmp = masssum
      call mpi_allreduce(tmp,masssum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      tmp = qsum
      call mpi_allreduce(tmp,qsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      tmp = qsqsum
      call mpi_allreduce(tmp,qsqsum,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)

      call mpi_barrier(mpi_comm_world,ierror)

      return
      end
