c Peek at (read) a binary LAMMPS restart file (version 2001 and before)
c MUST compile and run this program on the same machine the 
c    restart file was created on
c Syntax: peek_restart restart.file

      program peek_restart
      implicit real*8 (a-h,o-z)

      character*200 filename
      dimension buf(1000)
      real*8 omega(3),omega_dot(3)
      integer restart_version,iperatom
      integer bondstyle,anglestyle,dihedstyle,improstyle
      real*8, allocatable :: mass(:)
      real*8, allocatable :: noncoeff1(:,:)
      real*8, allocatable :: noncoeff2(:,:)
      real*8, allocatable :: noncoeff3(:,:)
      real*8, allocatable :: noncoeff4(:,:)
      real*8, allocatable :: noncoeff14_1(:,:)
      real*8, allocatable :: noncoeff14_2(:,:)
      real*8, allocatable :: bondcoeff(:,:)
      real*8, allocatable :: anglecoeff(:,:)
      real*8, allocatable :: bondbondcoeff(:,:)
      real*8, allocatable :: bondanglecoeff(:,:)
      real*8, allocatable :: dihedcoeff(:,:)
      real*8, allocatable :: midbondtorsioncoeff(:,:)
      real*8, allocatable :: endbondtorsioncoeff(:,:)
      real*8, allocatable :: angletorsioncoeff(:,:)
      real*8, allocatable :: angleangletorsioncoeff(:,:)
      real*8, allocatable :: bondbond13coeff(:,:)
      real*8, allocatable :: improcoeff(:,:)
      real*8, allocatable :: angleanglecoeff(:,:)

 900  format (a)

c user prompts

      if (iargc() /= 1) then
        write (6,*) 'Syntax: peek_restart restart.file'
        call exit(0)
      else
        call getarg(1,filename)
      endif

      write (6,*) 'LAMMPS version:'
      write (6,*) '  2001 = LAMMPS 2001'
      write (6,*) '  2000 = LAMMPS 2000'
      write (6,*) '  6 = LAMMPS 99'
      write (6,*) '  5 = LAMMPS - Version 5'
      read (5,*) restart_version
      write (6,*) 'Print out every Nth atom (N = 0 -> no atoms)'
      read (5,*) nevery

      if (restart_version /= 2001 .and. restart_version /= 2000 .and.
     $     restart_version /= 6 .and. restart_version /= 5) then
        write (6,*) 'ERROR: Illegal version request'
        stop
      endif

      open (unit=2,file=filename,form='unformatted',status='old')

      if (restart_version == 5 .or. restart_version == 6) then
        read (2) iversion,nbyte,ibyte
      else
        read (2) iversion
      endif

      if (iversion /= restart_version) then
        write (6,*) 'ERROR: File version does not match request'
        stop
      endif

c read header

      read (2) ntimestep,nprocs,igridx,igridy,igridz
      read (2) idim,iflagx,iflagy,iflagz
      read (2) iunits,inewton
      read (2) natoms,nbonds,nangles,ndihedrals,nimpropers
      read (2) ntypes
      if (nbonds.gt.0) read (2) nbondtypes
      if (nangles.gt.0) read (2) nangletypes
      if (ndihedrals.gt.0) read (2) ndihedtypes
      if (nimpropers.gt.0) read (2) nimprotypes
      read (2) maxbondper,maxangleper,maxdihedper,maximproper
      read (2) xboundlo,xboundhi,yboundlo,yboundhi,
     $     zboundlo,zboundhi
      if (restart_version == 2000 .or. restart_version == 2001) then
        read (2) eta,eta_dot
        read (2) omega(1),omega(2),omega(3),
     $       omega_dot(1),omega_dot(2),omega_dot(3)
      else
        read (2) eta,eta_dot,omega_dot(1),omega_dot(2),omega_dot(3)
        omega(1) = 0.0
        omega(2) = 0.0
        omega(3) = 0.0
      endif
      
      write (6,*) 'Version =',iversion
      write (6,*) 'Ntimestep =',ntimestep
      write (6,*) 'Nprocs =',nprocs
      write (6,*) 'Proc grid =',igridx,igridy,igridz
      write (6,*) 'Dimension =',idim
      write (6,*) 'Periodicity =',iflagx,iflagy,iflagz
      write (6,*) 'Units =',iunits
      write (6,*) 'Newton flag =',inewton
      write (6,*) 'Natoms =',natoms
      write (6,*) 'Nbonds =',nbonds
      write (6,*) 'Nangles =',nangles
      write (6,*) 'Ndihedrals =',ndihedrals
      write (6,*) 'Nimpropers =',nimpropers
      write (6,*) 'Ntypes =',ntypes
      if (nbonds.gt.0) write (6,*) 'Nbondtypes =',nbondtypes
      if (nangles.gt.0) write (6,*) 'Nangletypes =',nangletypes
      if (ndihedrals.gt.0) write (6,*) 'Ndihedtypes =',ndihedtypes
      if (nimpropers.gt.0) write (6,*) 'Nimprotypes =',nimprotypes
      write (6,*) 'Max_bondper =',maxbondper
      write (6,*) 'Max_angleper =',maxangleper
      write (6,*) 'Max_dihedper =',maxdihedper
      write (6,*) 'Max_improper =',maximproper
      write (6,*) 'Xlo/hi =',xboundlo,xboundhi
      write (6,*) 'Ylo/hi =',yboundlo,yboundhi
      write (6,*) 'Zlo/hi =',zboundlo,zboundhi
      write (6,*) 'Eta/Eta dot =',eta,eta_dot
      write (6,*) 'Omega =',omega(1),omega(2),omega(3)
      write (6,*) 'Omega dot =',omega_dot(1),omega_dot(2),omega_dot(3)

c allocate memory

      allocate(mass(ntypes))

      allocate(noncoeff1(ntypes,ntypes))
      allocate(noncoeff2(ntypes,ntypes))
      allocate(noncoeff3(ntypes,ntypes))
      allocate(noncoeff4(ntypes,ntypes))
      allocate(noncoeff14_1(ntypes,ntypes))
      allocate(noncoeff14_2(ntypes,ntypes))

      if (nbonds > 0) allocate(bondcoeff(5,nbondtypes))
      if (nangles > 0) then
        allocate(anglecoeff(5,nangletypes))
        allocate(bondbondcoeff(5,nangletypes))
        allocate(bondanglecoeff(5,nangletypes))
      endif
      if (ndihedrals > 0) then
        allocate(dihedcoeff(6,ndihedtypes))
        allocate(midbondtorsioncoeff(4,ndihedtypes))
        allocate(endbondtorsioncoeff(8,ndihedtypes))
        allocate(angletorsioncoeff(8,ndihedtypes))
        allocate(angleangletorsioncoeff(3,ndihedtypes))
        allocate(bondbond13coeff(3,ndihedtypes))
      endif
      if (nimpropers > 0) then
        allocate(improcoeff(3,nimprotypes))
        allocate(angleanglecoeff(6,nimprotypes))
      endif

c read and print coeffs

      read (2) (mass(i),i=1,ntypes)

      write (6,*) 'Masses =',(mass(i),i=1,ntypes)

      read (2) nonstyle
      if (nonstyle.eq.0) then
        write (6,*) 'Nonbond style set to none'
      else if (nonstyle.eq.1) then
        write (6,*) 'Nonbond style set to lj/cutoff'
        read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle.eq.2) then
        write (6,*) 'Nonbond style set to lj/smooth'
        read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle.eq.3) then
        write (6,*) 'Nonbond style set to lj/shift'
        read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle.eq.4) then
        write (6,*) 'Nonbond style set to soft'
        read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle.eq.5) then
        write (6,*) 'Nonbond style set to class2/cutoff'
        read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle.eq.6) then
        write (6,*) 'Nonbond style set to lj/charmm'
        read (2) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff14_1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff14_2(i,j),j=i,ntypes),i=1,ntypes)
      endif

      if (nonstyle == 1) then
        do i = 1,ntypes
          do j = 1,ntypes
            write (6,*) i,j,
     $           noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
          enddo
        enddo
      else if (nonstyle == 2) then
        do i = 1,ntypes
          do j = 1,ntypes
            write (6,*) i,j,
     $           noncoeff1(i,j),noncoeff2(i,j),
     $           noncoeff3(i,j),noncoeff4(i,j)
          enddo
        enddo
      else if (nonstyle == 3) then
        do i = 1,ntypes
          do j = 1,ntypes
            write (6,*) i,j,
     $           noncoeff1(i,j),noncoeff2(i,j),
     $           noncoeff3(i,j),noncoeff4(i,j)
          enddo
        enddo
      else if (nonstyle == 4) then
        do i = 1,ntypes
          do j = 1,ntypes
            write (6,*) i,j,
     $           noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
          enddo
        enddo
      else if (nonstyle == 5) then
        do i = 1,ntypes
          do j = 1,ntypes
            write (6,*) i,j,
     $           noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
          enddo
        enddo
      else if (nonstyle == 6) then
        do i = 1,ntypes
          do j = 1,ntypes
            write (6,*) i,j,
     $           noncoeff1(i,j),noncoeff2(i,j),
     $           noncoeff14_1(i,j),noncoeff14_2(i,j)
          enddo
        enddo
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

      if (bondstyle == 1) then
        do j = 1,nbondtypes
          write (6,*) j,(bondcoeff(i,j),i=1,2)
        enddo
      else if (bondstyle == 2) then
        do j = 1,nbondtypes
          write (6,*) j,(bondcoeff(i,j),i=1,4)
        enddo
      else if (bondstyle == 3) then
        do j = 1,nbondtypes
          write (6,*) j,(bondcoeff(i,j),i=1,5)
        enddo
      else if (bondstyle == 4) then
        do j = 1,nbondtypes
          write (6,*) j,(bondcoeff(i,j),i=1,3)
        enddo
      else if (bondstyle == 5) then
        do j = 1,nbondtypes
          write (6,*) j,(bondcoeff(i,j),i=1,4)
        enddo
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
        endif
      endif

      if (anglestyle == 1) then
        do j = 1,nangletypes
          write (6,*) j,(anglecoeff(i,j),i=1,2)
        enddo
      else if (anglestyle == 2) then
        do j = 1,nangletypes
          write (6,*) j,(anglecoeff(i,j),i=1,4)
          write (6,*) j,(bondbondcoeff(i,j),i=1,3)
          write (6,*) j,(bondanglecoeff(i,j),i=1,4)
        enddo
      else if (anglestyle == 3) then
        do j = 1,nangletypes
          write (6,*) j,(anglecoeff(i,j),i=1,4)
        enddo
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
     $         j=1,ndihedtypes)
          read (2) ((endbondtorsioncoeff(i,j),i=1,8),
     $         j=1,ndihedtypes)
          read (2) ((angletorsioncoeff(i,j),i=1,8),
     $         j=1,ndihedtypes)
          read (2) ((angleangletorsioncoeff(i,j),i=1,3),
     $         j=1,ndihedtypes)
          if (restart_version == 2000 .or. restart_version == 2001)
     $         read (2) ((bondbond13coeff(i,j),i=1,3),
     $         j=1,ndihedtypes)
        else if (dihedstyle.eq.3) then
          write (6,*) 'Dihedral style set to multiharmonic'
          read (2) ((dihedcoeff(i,j),i=1,5),j=1,ndihedtypes)
        else if (dihedstyle.eq.4) then
          write (6,*) 'Dihedral style set to charmm'
          read (2) ((dihedcoeff(i,j),i=1,4),j=1,ndihedtypes)
        endif
      endif

      if (dihedstyle == 1) then
        do j = 1,ndihedtypes
          write (6,*) j,(dihedcoeff(i,j),i=1,3)
        enddo
      else if (dihedstyle == 2) then
        do j = 1,ndihedtypes
          write (6,*) j,(dihedcoeff(i,j),i=1,6)
          write (6,*) j,(midbondtorsioncoeff(i,j),i=1,4)
          write (6,*) j,(endbondtorsioncoeff(i,j),i=1,8)
          write (6,*) j,(angletorsioncoeff(i,j),i=1,8)
          write (6,*) j,(angleangletorsioncoeff(i,j),i=1,3)
          if (restart_version == 2000 .or. restart_version == 2001)
     $         write (6,*) j,(bondbond13coeff(i,j),i=1,3)
        enddo
      else if (dihedstyle == 3) then
        do j = 1,ndihedtypes
          write (6,*) j,(dihedcoeff(i,j),i=1,5)
        enddo
      else if (dihedstyle == 4) then
        do j = 1,ndihedtypes
          write (6,*) j,(dihedcoeff(i,j),i=1,4)
        enddo
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

      if (improstyle == 1) then
        do j = 1,nimprotypes
          write (6,*) j,(improcoeff(i,j),i=1,2)
        enddo
      else if (improstyle == 2) then
        do j = 1,nimprotypes
          write (6,*) j,(improcoeff(i,j),i=1,3)
        enddo
      else if (improstyle == 3) then
        do j = 1,nimprotypes
          write (6,*) j,(improcoeff(i,j),i=1,2)
          write (6,*) j,(angleanglecoeff(i,j),i=1,6)
        enddo
      endif

c read atoms

      if (nevery.eq.0) stop

      iperatom = 16
      do i = 1,natoms

        read (2) (buf(k),k=1,iperatom)

        x = buf(1)
        y = buf(2)
        z = buf(3)
        itag = nint(buf(4))
        molecule = nint(buf(5))
        itype = nint(buf(6))
        q = buf(7)
        vx = buf(8)
        vy = buf(9)
        vz = buf(10)
        itrue = nint(buf(11))
        ifix = nint(buf(12))

        numbond = nint(buf(13))
        numangle = nint(buf(14))
        numdihed = nint(buf(15))
        numimpro = nint(buf(16))

        j = iperatom
        do jj = 1,numbond
          read (2) buf(j+1),buf(j+2),buf(j+3)
          j = j + 3
        enddo
        do jj = 1,numangle
          read (2) buf(j+1),buf(j+2),buf(j+3),buf(j+4)
          j = j + 4
        enddo
        do jj = 1,numdihed
          read (2) buf(j+1),buf(j+2),buf(j+3),
     $         buf(j+4),buf(j+5)
          j = j + 5
        enddo
        do jj = 1,numimpro
          read (2) buf(j+1),buf(j+2),buf(j+3),
     $         buf(j+4),buf(j+5)
          j = j + 5
        enddo

        if (i == 1 .or. mod(i,nevery) == 0) then
          write (6,*) ' Atom:',itag,molecule,itype,
     $         x,y,z,q,vx,vy,vz,itrue,ifix,
     $         numbond,numangle,numdihed,numimpro
          j = iperatom
          do jj = 1,numbond
            write (6,*) '    Bond:',
     $           nint(buf(j+1)),nint(buf(j+2)),nint(buf(j+3))
            j = j + 3
          enddo
          do jj = 1,numangle
            write (6,*) '    Angle:',
     $           nint(buf(j+1)),nint(buf(j+2)),
     $           nint(buf(j+3)),nint(buf(j+4))
            j = j + 4
          enddo
          do jj = 1,numdihed
            write (6,*) '    Dihedral:',
     $           nint(buf(j+1)),nint(buf(j+2)),
     $           nint(buf(j+3)),nint(buf(j+4)),nint(buf(j+5))
            j = j + 5
          enddo
          do jj = 1,numimpro
            write (6,*) '    Improper:',
     $           nint(buf(j+1)),nint(buf(j+2)),
     $           nint(buf(j+3)),nint(buf(j+4)),nint(buf(j+5))
            j = j + 5
          enddo
        endif

      enddo

      end
