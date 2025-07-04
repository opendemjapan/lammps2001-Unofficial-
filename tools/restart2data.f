c Convert a binary restart file into a LAMMPS data file (version 2001)
c MUST compile and run this program on the same machine the 
c   restart file was created on
c Syntax: restart2data restart.file data.file

      program restart2data
      implicit real*8 (a-h,o-z)

      integer restart_version
      integer bondstyle,anglestyle,dihedstyle,improstyle
      real*8 buf(1000)
      real*8 omega(3),omega_dot(3)
      character*200 restartfile,datafile

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
      real*8, allocatable :: x(:,:),q(:),v(:,:)

      integer, allocatable :: tag(:),type(:),molecule(:),true(:)
      integer, allocatable :: numbond(:),numangle(:)
      integer, allocatable :: numdihed(:),numimpro(:)
      integer, allocatable :: bondatom1(:,:),bondatom2(:,:)
      integer, allocatable :: bondtype(:,:)
      integer, allocatable :: angleatom1(:,:),angleatom2(:,:)
      integer, allocatable :: angleatom3(:,:),angletype(:,:)
      integer, allocatable :: dihedatom1(:,:),dihedatom2(:,:)
      integer, allocatable :: dihedatom3(:,:),dihedatom4(:,:)
      integer, allocatable :: dihedtype(:,:)
      integer, allocatable :: improatom1(:,:),improatom2(:,:)
      integer, allocatable :: improatom3(:,:),improatom4(:,:)
      integer, allocatable :: improtype(:,:)
      integer, allocatable :: bondatom1all(:),bondatom2all(:)
      integer, allocatable :: bondtypeall(:)
      integer, allocatable :: angleatom1all(:),angleatom2all(:)
      integer, allocatable :: angleatom3all(:),angletypeall(:)
      integer, allocatable :: dihedatom1all(:),dihedatom2all(:)
      integer, allocatable :: dihedatom3all(:),dihedatom4all(:)
      integer, allocatable :: dihedtypeall(:)
      integer, allocatable :: improatom1all(:),improatom2all(:)
      integer, allocatable :: improatom3all(:),improatom4all(:)
      integer, allocatable :: improtypeall(:)

      integer, allocatable :: indx(:)

c 904 is format statement for atoms
c can make fields larger if needed

 900  format(a)
 901  format(2f10.6,a)
 902  format(i3,f8.4)
 903  format(i3,8f12.6)
 904  format(i7,i7,i4,f8.3,3f10.5,3i4)
 905  format(i6,i4,4i7)
 906  format(i3,f12.6,i4,i4)
 907  format(i3,f12.6,i4,2f12.6)
 908  format(i7,3f12.6)

c command line args and prompts

      if (iargc().ne.2) then
        write (6,*) 'Syntax: restart2data restart.file data.file'
        call exit(0)
      else
        call getarg(1,restartfile)
        call getarg(2,datafile)
      endif

      write (6,*) 'LAMMPS version:'
      write (6,*) '  2001 = LAMMPS 2001'
      write (6,*) '  2000 = LAMMPS 2000'
      write (6,*) '  6 = LAMMPS 99'
      write (6,*) '  5 = LAMMPS - Version 5'
      read (5,*) restart_version

      if (restart_version /= 2001 .and. restart_version /= 2000 .and.
     $     restart_version /= 6 .and. restart_version /= 5) then
        write (6,*) 'ERROR: Illegal version request'
        stop
      endif

      open (unit=1,file=restartfile,form='unformatted',status='old')

      if (restart_version == 5 .or. restart_version == 6) then
        read (1) iversion,nbyte,ibyte
      else
        read (1) iversion
      endif

      if (iversion /= restart_version) then
        write (6,*) 'ERROR: File version does not match request'
        stop
      endif

      write (6,*) 'Include true flags in data file: (0) no (1) yes'
      read (5,*) itrue

c read header

      read (1) ntimestep,nprocs,igridx,igridy,igridz
      read (1) idim,iflagx,iflagy,iflagz
      read (1) iunits,inewton
      read (1) natoms,nbonds,nangles,ndihedrals,nimpropers
      read (1) ntypes
      if (nbonds > 0) read (1) nbondtypes
      if (nangles > 0) read (1) nangletypes
      if (ndihedrals > 0) read (1) ndihedtypes
      if (nimpropers > 0) read (1) nimprotypes
      read (1) maxbondper,maxangleper,maxdihedper,maximproper
      read (1) xboundlo,xboundhi,yboundlo,yboundhi,
     $     zboundlo,zboundhi
      if (restart_version == 2000 .or. restart_version == 2001) then
        read (1) eta,eta_dot
        read (1) omega(1),omega(2),omega(3),
     $       omega_dot(1),omega_dot(2),omega_dot(3)
      else
        read (1) eta,eta_dot,omega_dot(1),omega_dot(2),omega_dot(3)
      endif
      
c allocate coeff memory

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

c read coeffs

      read (1) (mass(i),i=1,ntypes)

      read (1) nonstyle
      if (nonstyle == 0) then
        write (6,*) 'Nonbond style = none'
      else if (nonstyle == 1) then
        write (6,*) 'Nonbond style = lj/cutoff'
        read (1) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle == 2) then
        write (6,*) 'Nonbond style = lj/smooth'
        read (1) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle == 3) then
        write (6,*) 'Nonbond style = lj/shift'
        read (1) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff4(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle == 4) then
        write (6,*) 'Nonbond style = soft'
        read (1) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle == 5) then
        write (6,*) 'Nonbond style = class2/cutoff'
        read (1) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff3(i,j),j=i,ntypes),i=1,ntypes)
      else if (nonstyle == 6) then
        write (6,*) 'Nonbond style = lj/charmm'
        read (1) ((noncoeff1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff2(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff14_1(i,j),j=i,ntypes),i=1,ntypes),
     $       ((noncoeff14_2(i,j),j=i,ntypes),i=1,ntypes)
      endif

      if (nbonds > 0) then
        read (1) bondstyle
        if (bondstyle == 0) then
          write (6,*) 'Bond style = none'
        else if (bondstyle == 1) then
          write (6,*) 'Bond style = harmonic'
          read (1) ((bondcoeff(i,j),i=1,2),j=1,nbondtypes)
        else if (bondstyle == 2) then
          write (6,*) 'Bond style = fene/standard'
          read (1) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
        else if (bondstyle == 3) then
          write (6,*) 'Bond style = fene/shift'
          read (1) ((bondcoeff(i,j),i=1,5),j=1,nbondtypes)
        else if (bondstyle == 4) then
          write (6,*) 'Bond style = nonlinear'
          read (1) ((bondcoeff(i,j),i=1,3),j=1,nbondtypes)
        else if (bondstyle == 5) then
          write (6,*) 'Bond style = class2'
          read (1) ((bondcoeff(i,j),i=1,4),j=1,nbondtypes)
        endif
      endif

      if (nangles > 0) then
        read (1) anglestyle
        if (anglestyle == 0) then
          write (6,*) 'Angle style = none'
        else if (anglestyle == 1) then
          write (6,*) 'Angle style = harmonic'
          read (1) ((anglecoeff(i,j),i=1,2),j=1,nangletypes)
          do i = 1,nangletypes
            anglecoeff(2,i) = anglecoeff(2,i)/3.1415926 * 180.0
          enddo
        else if (anglestyle == 2) then
          write (6,*) 'Angle style = class2'
          read (1) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
          read (1) ((bondbondcoeff(i,j),i=1,3),j=1,nangletypes)
          read (1) ((bondanglecoeff(i,j),i=1,4),j=1,nangletypes)
          do i = 1,nangletypes
            anglecoeff(1,i) = anglecoeff(1,i)/3.1415926 * 180.0
          enddo
        else if (anglestyle == 3) then
          write (6,*) 'Angle style = charmm'
          read (1) ((anglecoeff(i,j),i=1,4),j=1,nangletypes)
          do i = 1,nangletypes
            anglecoeff(2,i) = anglecoeff(2,i)/3.1415926 * 180.0
          enddo
        endif
      endif

      if (ndihedrals > 0) then
        read (1) dihedstyle
        if (dihedstyle == 0) then
          write (6,*) 'Dihedral style = none'
        else if (dihedstyle == 1) then
          write (6,*) 'Dihedral style = harmonic'
          read (1) ((dihedcoeff(i,j),i=1,3),j=1,ndihedtypes)
        else if (dihedstyle == 2) then
          write (6,*) 'Dihedral style = class2'
          read (1) ((dihedcoeff(i,j),i=1,6),j=1,ndihedtypes)
          read (1) ((midbondtorsioncoeff(i,j),i=1,4),
     $         j=1,ndihedtypes)
          read (1) ((endbondtorsioncoeff(i,j),i=1,8),
     $         j=1,ndihedtypes)
          read (1) ((angletorsioncoeff(i,j),i=1,8),
     $         j=1,ndihedtypes)
          read (1) ((angleangletorsioncoeff(i,j),i=1,3),
     $         j=1,ndihedtypes)
          if (restart_version == 2000 .or. restart_version == 2001) then
            read (1) ((bondbond13coeff(i,j),i=1,3),j=1,ndihedtypes)
          else
            do j = 1,ndihedtypes
              do i = 1,3
                bondbond13coeff(i,j) = 0.0
              enddo
            enddo
          endif
          do i = 1,ndihedtypes
            dihedcoeff(2,i) = dihedcoeff(2,i)/3.1415926 * 180.0
            dihedcoeff(4,i) = dihedcoeff(4,i)/3.1415926 * 180.0
            dihedcoeff(6,i) = dihedcoeff(6,i)/3.1415926 * 180.0
          enddo
          do i = 1,ndihedtypes
            angletorsioncoeff(7,i) =
     $           angletorsioncoeff(7,i)/3.1415926 * 180.0
            angletorsioncoeff(8,i) =
     $           angletorsioncoeff(8,i)/3.1415926 * 180.0
          enddo
          do i = 1,ndihedtypes
            angleangletorsioncoeff(2,i) =
     $           angleangletorsioncoeff(2,i)/3.1415926 * 180.0
            angleangletorsioncoeff(3,i) =
     $           angleangletorsioncoeff(3,i)/3.1415926 * 180.0
          enddo
        else if (dihedstyle == 3) then
          write (6,*) 'Dihedral style = multiharmonic'
          read (1) ((dihedcoeff(i,j),i=1,5),j=1,ndihedtypes)
        else if (dihedstyle == 4) then
          write (6,*) 'Dihedral style = charmm'
          read (1) ((dihedcoeff(i,j),i=1,4),j=1,ndihedtypes)
        endif
      endif

      if (nimpropers > 0) then
        read (1) improstyle
        if (improstyle == 0) then
          write (6,*) 'Improper style = none'
        else if (improstyle == 1) then
          write (6,*) 'Improper style = harmonic'
          read (1) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
          do i = 1,nimprotypes
            improcoeff(2,i) = improcoeff(2,i)/3.1415926 * 180.0
          enddo
        else if (improstyle == 2) then
          write (6,*) 'Improper style = cvff'
          read (1) ((improcoeff(i,j),i=1,3),j=1,nimprotypes)
        else if (improstyle == 3) then
          write (6,*) 'Improper style = class2'
          read (1) ((improcoeff(i,j),i=1,2),j=1,nimprotypes)
          read (1) ((angleanglecoeff(i,j),i=1,6),j=1,nimprotypes)
          do i = 1,nimprotypes
            improcoeff(2,i) = improcoeff(2,i)/3.1415926 * 180.0
          enddo
          do i = 1,nimprotypes
            angleanglecoeff(4,i) =
     $           angleanglecoeff(4,i)/3.1415926 * 180.0
            angleanglecoeff(5,i) =
     $           angleanglecoeff(5,i)/3.1415926 * 180.0
            angleanglecoeff(6,i) =
     $           angleanglecoeff(6,i)/3.1415926 * 180.0
          enddo
        endif
      endif

c allocate atom and topology memory

      allocate(x(3,natoms))
      allocate(v(3,natoms))
      allocate(q(natoms))

      allocate(tag(natoms))
      allocate(type(natoms))
      allocate(molecule(natoms))
      allocate(true(natoms))

      allocate(numbond(natoms))
      allocate(numangle(natoms))
      allocate(numdihed(natoms))
      allocate(numimpro(natoms))

      allocate(bondatom1(maxbondper,natoms))
      allocate(bondatom2(maxbondper,natoms))
      allocate(bondtype(maxbondper,natoms))

      allocate(angleatom1(maxangleper,natoms))
      allocate(angleatom2(maxangleper,natoms))
      allocate(angleatom3(maxangleper,natoms))
      allocate(angletype(maxangleper,natoms))

      allocate(dihedatom1(maxdihedper,natoms))
      allocate(dihedatom2(maxdihedper,natoms))
      allocate(dihedatom3(maxdihedper,natoms))
      allocate(dihedatom4(maxdihedper,natoms))
      allocate(dihedtype(maxdihedper,natoms))

      allocate(improatom1(maximproper,natoms))
      allocate(improatom2(maximproper,natoms))
      allocate(improatom3(maximproper,natoms))
      allocate(improatom4(maximproper,natoms))
      allocate(improtype(maximproper,natoms))

c read atoms

      iperatom = 16

      do i = 1,natoms

        read (1) (buf(j),j=1,iperatom)
        j = iperatom
        numbond_tmp = nint(buf(j-3))
        numangle_tmp = nint(buf(j-2))
        numdihed_tmp = nint(buf(j-1))
        numimpro_tmp = nint(buf(j))
        do jj = 1,numbond_tmp
          read (1) buf(j+1),buf(j+2),buf(j+3)
          j = j + 3
        enddo
        do jj = 1,numangle_tmp
          read (1) buf(j+1),buf(j+2),buf(j+3),buf(j+4)
          j = j + 4
        enddo
        do jj = 1,numdihed_tmp
          read (1) buf(j+1),buf(j+2),buf(j+3),
     $         buf(j+4),buf(j+5)
          j = j + 5
        enddo
        do jj = 1,numimpro_tmp
          read (1) buf(j+1),buf(j+2),buf(j+3),
     $         buf(j+4),buf(j+5)
          j = j + 5
        enddo

        x(1,i) = buf(1)
        x(2,i) = buf(2)
        x(3,i) = buf(3)
        tag(i) = nint(buf(4))
        molecule(i) = nint(buf(5))
        type(i) = nint(buf(6))
        q(i) = buf(7)
        v(1,i) = buf(8)
        v(2,i) = buf(9)
        v(3,i) = buf(10)
        true(i) = nint(buf(11))
        ifix = nint(buf(12))
        numbond(i) = nint(buf(13))
        numangle(i) = nint(buf(14))
        numdihed(i) = nint(buf(15))
        numimpro(i) = nint(buf(16))

        n = iperatom
        do jj = 1,numbond(i)
          bondatom1(jj,i) = nint(buf(n+1))
          bondatom2(jj,i) = nint(buf(n+2))
          bondtype(jj,i) = nint(buf(n+3))
          n = n + 3
        enddo
        do jj = 1,numangle(i)
          angleatom1(jj,i) = nint(buf(n+1))
          angleatom2(jj,i) = nint(buf(n+2))
          angleatom3(jj,i) = nint(buf(n+3))
          angletype(jj,i) = nint(buf(n+4))
          n = n + 4
        enddo
        do jj = 1,numdihed(i)
          dihedatom1(jj,i) = nint(buf(n+1))
          dihedatom2(jj,i) = nint(buf(n+2))
          dihedatom3(jj,i) = nint(buf(n+3))
          dihedatom4(jj,i) = nint(buf(n+4))
          dihedtype(jj,i) = nint(buf(n+5))
          n = n + 5
        enddo
        do jj = 1,numimpro(i)
          improatom1(jj,i) = nint(buf(n+1))
          improatom2(jj,i) = nint(buf(n+2))
          improatom3(jj,i) = nint(buf(n+3))
          improatom4(jj,i) = nint(buf(n+4))
          improtype(jj,i) = nint(buf(n+5))
          n = n + 5
        enddo

      enddo

      close (1)

c unscale velocities

      if (iunits == 0) dtfactor = 48.88821
      if (iunits == 1) dtfactor = 1.0

      do i = 1,natoms
        v(1,i) = v(1,i)/dtfactor
        v(2,i) = v(2,i)/dtfactor
        v(3,i) = v(3,i)/dtfactor
      enddo

c create global topology lists from individual atom lists

      if (nbonds > 0) then

        allocate(bondtypeall(nbonds))
        allocate(bondatom1all(nbonds))
        allocate(bondatom2all(nbonds))

        m = 0
        do i = 1,natoms
          do j = 1,numbond(i)
            if (tag(i) == bondatom1(j,i)) then
              m = m + 1
              bondtypeall(m) = bondtype(j,i)
              bondatom1all(m) = bondatom1(j,i)
              bondatom2all(m) = bondatom2(j,i)
            endif
          enddo
        enddo

        if (m /= nbonds) then
          write (6,*) 'ERROR: Mismatch in bond topologies'
          stop
        endif

      endif

      if (nangles > 0) then

        allocate(angletypeall(nangles))
        allocate(angleatom1all(nangles))
        allocate(angleatom2all(nangles))
        allocate(angleatom3all(nangles))

        m = 0
        do i = 1,natoms
          do j = 1,numangle(i)
            if (tag(i) == angleatom2(j,i)) then
              m = m + 1
              angletypeall(m) = angletype(j,i)
              angleatom1all(m) = angleatom1(j,i)
              angleatom2all(m) = angleatom2(j,i)
              angleatom3all(m) = angleatom3(j,i)
            endif
          enddo
        enddo

        if (m /= nangles) then
          write (6,*) 'ERROR: Mismatch in angle topologies'
          stop
        endif

      endif

      if (ndihedrals > 0) then

        allocate(dihedtypeall(ndihedrals))
        allocate(dihedatom1all(ndihedrals))
        allocate(dihedatom2all(ndihedrals))
        allocate(dihedatom3all(ndihedrals))
        allocate(dihedatom4all(ndihedrals))

        m = 0
        do i = 1,natoms
          do j = 1,numdihed(i)
            if (tag(i) == dihedatom2(j,i)) then
              m = m + 1
              dihedtypeall(m) = dihedtype(j,i)
              dihedatom1all(m) = dihedatom1(j,i)
              dihedatom2all(m) = dihedatom2(j,i)
              dihedatom3all(m) = dihedatom3(j,i)
              dihedatom4all(m) = dihedatom4(j,i)
            endif
          enddo
        enddo

        if (m /= ndihedrals) then
          write (6,*) 'ERROR: Mismatch in dihedral topologies'
          stop
        endif

      endif

      if (nimpropers > 0) then

        allocate(improtypeall(nimpropers))
        allocate(improatom1all(nimpropers))
        allocate(improatom2all(nimpropers))
        allocate(improatom3all(nimpropers))
        allocate(improatom4all(nimpropers))

        m = 0
        do i = 1,natoms
          do j = 1,numimpro(i)
            if (tag(i) == improatom2(j,i)) then
              m = m + 1
              improtypeall(m) = improtype(j,i)
              improatom1all(m) = improatom1(j,i)
              improatom2all(m) = improatom2(j,i)
              improatom3all(m) = improatom3(j,i)
              improatom4all(m) = improatom4(j,i)
            endif
          enddo
        enddo

        if (m /= nimpropers) then
          write (6,*) 'ERROR: Mismatch in improper topologies'
          stop
        endif

      endif

c write LAMMPS data file header

      open (1,file=datafile)

      write (1,900) 'LAMMPS Description'
      write (1,*)

      write (1,*) natoms,' atoms'
      write (1,*) nbonds,' bonds'
      write (1,*) nangles,' angles'
      write (1,*) ndihedrals,' dihedrals'
      write (1,*) nimpropers,' impropers'
      write (1,*)

      write (1,*) ntypes,' atom types'
      if (nbonds > 0) write (1,*) nbondtypes,' bond types'
      if (nangles > 0) write (1,*) nangletypes,' angle types'
      if (ndihedrals > 0) write (1,*) ndihedtypes,' dihedral types'
      if (nimpropers > 0) write (1,*) nimprotypes,' improper types'
      write (1,*)

      write (1,901) xboundlo,xboundhi,' xlo xhi'
      write (1,901) yboundlo,yboundhi,' ylo yhi'
      if (idim == 3) write (1,901) zboundlo,zboundhi,' zlo zhi'

c masses

      write (1,*)
      write (1,900) 'Masses'
      write (1,*)

      do i = 1,ntypes
        write (1,902) i,mass(i)
      enddo

c nonbond coeffs

      if (nonstyle > 0) then
        write (1,*)
        write (1,900) 'Nonbond Coeffs'
        write (1,*)
      endif

      if (nonstyle == 1) then
        do i = 1,ntypes
          write (1,903) i,noncoeff1(i,i),noncoeff2(i,i)
        enddo
      else if (nonstyle == 2) then
        do i = 1,ntypes
          write (1,903) i,noncoeff1(i,i),noncoeff2(i,i)
        enddo
      else if (nonstyle == 3) then
        do i = 1,ntypes
          write (1,903) i,noncoeff1(i,i),noncoeff2(i,i),noncoeff3(i,i)
        enddo
      else if (nonstyle == 4) then
        do i = 1,ntypes
          write (1,903) i,noncoeff1(i,i),noncoeff2(i,i)
        enddo
      else if (nonstyle == 5) then
        do i = 1,ntypes
          write (1,903) i,noncoeff1(i,i),noncoeff2(i,i)
        enddo
      else if (nonstyle == 6) then
        do i = 1,ntypes
          write (1,903) i,noncoeff1(i,i),noncoeff2(i,i),
     $         noncoeff14_1(i,i),noncoeff14_2(i,i)
        enddo
      endif

c bond coeffs

      if (bondstyle > 0) then
        write (1,*)
        write (1,900) 'Bond Coeffs'
        write (1,*)
      endif

      if (bondstyle == 1) then
        do j = 1,nbondtypes
          write (1,903) j,(bondcoeff(i,j),i=1,2)
        enddo
      else if (bondstyle == 2) then
        do j = 1,nbondtypes
          write (1,903) j,(bondcoeff(i,j),i=1,4)
        enddo
      else if (bondstyle == 3) then
        do j = 1,nbondtypes
          write (1,903) j,(bondcoeff(i,j),i=1,5)
        enddo
      else if (bondstyle == 4) then
        do j = 1,nbondtypes
          write (1,903) j,(bondcoeff(i,j),i=1,3)
        enddo
      else if (bondstyle == 5) then
        do j = 1,nbondtypes
          write (1,903) j,(bondcoeff(i,j),i=1,4)
        enddo
      endif

c angle coeffs

      if (anglestyle > 0) then
        write (1,*)
        write (1,900) 'Angle Coeffs'
        write (1,*)
      endif

      if (anglestyle == 1) then
        do j = 1,nangletypes
          write (1,903) j,(anglecoeff(i,j),i=1,2)
        enddo
      else if (anglestyle == 2) then
        do j = 1,nangletypes
          write (1,903) j,(anglecoeff(i,j),i=1,4)
        enddo
      else if (anglestyle == 3) then
        do j = 1,nangletypes
          write (1,903) j,(anglecoeff(i,j),i=1,4)
        enddo
      endif

      if (anglestyle == 2) then
        write (1,*)
        write (1,900) 'BondBond Coeffs'
        write (1,*)

        do j = 1,nangletypes
          write (1,903) j,(bondbondcoeff(i,j),i=1,3)
        enddo

        write (1,*)
        write (1,900) 'BondAngle Coeffs'
        write (1,*)

        do j = 1,nangletypes
          write (1,903) j,(bondanglecoeff(i,j),i=1,4)
        enddo
      endif

c dihedral coeffs

      if (dihedstyle > 0) then
        write (1,*)
        write (1,900) 'Dihedral Coeffs'
        write (1,*)
      endif

      if (dihedstyle == 1) then
        do j = 1,ndihedtypes
          write (1,906) j,dihedcoeff(1,j),
     $         nint(dihedcoeff(2,j)),nint(dihedcoeff(3,j))
        enddo
      else if (dihedstyle == 2) then
        do j = 1,ndihedtypes
          write (1,903) j,(dihedcoeff(i,j),i=1,6)
        enddo
      else if (dihedstyle == 3) then
        do j = 1,ndihedtypes
          write (1,903) j,(dihedcoeff(i,j),i=1,5)
        enddo
      else if (dihedstyle == 4) then
        do j = 1,ndihedtypes
          write (1,907) j,dihedcoeff(1,j),nint(dihedcoeff(2,j)),
     $         dihedcoeff(3,j),dihedcoeff(4,j)
        enddo
      endif

      if (dihedstyle == 2) then
        write (1,*)
        write (1,900) 'MiddleBondTorsion Coeffs'
        write (1,*)

        do j = 1,ndihedtypes
          write (1,903) j,(midbondtorsioncoeff(i,j),i=1,4)
        enddo

        write (1,*)
        write (1,900) 'EndBondTorsion Coeffs'
        write (1,*)

        do j = 1,ndihedtypes
          write (1,903) j,(endbondtorsioncoeff(i,j),i=1,8)
        enddo

        write (1,*)
        write (1,900) 'AngleTorsion Coeffs'
        write (1,*)

        do j = 1,ndihedtypes
          write (1,903) j,(angletorsioncoeff(i,j),i=1,8)
        enddo

        write (1,*)
        write (1,900) 'AngleAngleTorsion Coeffs'
        write (1,*)

        do j = 1,ndihedtypes
          write (1,903) j,(angleangletorsioncoeff(i,j),i=1,3)
        enddo

        write (1,*)
        write (1,900) 'BondBond13 Coeffs'
        write (1,*)

        do j = 1,ndihedtypes
          write (1,903) j,(bondbond13coeff(i,j),i=1,3)
        enddo
      endif

c improper coeffs

      if (improstyle > 0) then
        write (1,*)
        write (1,900) 'Improper Coeffs'
        write (1,*)
      endif

      if (improstyle == 1) then
        do j = 1,nimprotypes
          write (1,903) j,(improcoeff(i,j),i=1,2)
        enddo
      else if (improstyle == 2) then
        do j = 1,nimprotypes
          write (1,906) j,improcoeff(1,j),
     $         nint(improcoeff(2,j)),nint(improcoeff(3,j))
        enddo
      else if (improstyle == 3) then
        do j = 1,nimprotypes
          write (1,903) j,(improcoeff(i,j),i=1,2)
        enddo
      endif

      if (improstyle == 3) then
        write (1,*)
        write (1,900) 'AngleAngle Coeffs'
        write (1,*)

        do j = 1,nimprotypes
          write (1,903) j,(angleanglecoeff(i,j),i=1,6)
        enddo
      endif

c write atoms in sorted order via tags

      allocate(indx(natoms))
      call indexx(natoms,tag,indx)

      write (1,*)
      write (1,900) 'Atoms'
      write (1,*)

      if (itrue == 0) then
        do i = 1,natoms
          j = indx(i)
          write (1,904) tag(j),molecule(j),type(j),q(j),
     $         x(1,j),x(2,j),x(3,j)
        enddo
      else
        do i = 1,natoms
          j = indx(i)
          itruex = mod(true(j),1000) - 500
          itruey = mod(true(j)/1000,1000) - 500
          itruez = true(j)/1000000 - 500
          write (1,904) tag(j),molecule(j),type(j),q(j),
     $         x(1,j),x(2,j),x(3,j),itruex,itruey,itruez
        enddo
      endif

c write velocities in same sorted order as atoms

      write (1,*)
      write (1,900) 'Velocities'
      write (1,*)

      do i = 1,natoms
        j = indx(i)
        write (1,908) tag(j),v(1,j),v(2,j),v(3,j)
      enddo
      
      deallocate(indx)

c write bonds in sorted order via bondatom1all, ignore ties

      if (nbonds > 0) then

        allocate(indx(nbonds))
        call indexx(nbonds,bondatom1all,indx)

        write (1,*)
        write (1,900) 'Bonds'
        write (1,*)

        do i = 1,nbonds
          j = indx(i)
          write (1,905) i,bondtypeall(j),
     $         bondatom1all(j),bondatom2all(j)
        enddo

        deallocate(indx)

      endif

c write angles in sorted order via angleatom2all, ignore ties

      if (nangles > 0) then

        allocate(indx(nangles))
        call indexx(nangles,angleatom2all,indx)

        write (1,*)
        write (1,900) 'Angles'
        write (1,*)

        do i = 1,nangles
          j = indx(i)
          write (1,905) i,angletypeall(j),
     $         angleatom1all(j),angleatom2all(j),angleatom3all(j)
        enddo

        deallocate(indx)

      endif

c write dihedrals in sorted order via dihedatom2all, ignore ties

      if (ndihedrals > 0) then

        allocate(indx(ndihedrals))
        call indexx(ndihedrals,dihedatom2all,indx)

        write (1,*)
        write (1,900) 'Dihedrals'
        write (1,*)

        do i = 1,ndihedrals
          j = indx(i)
          write (1,905) i,dihedtypeall(j),
     $         dihedatom1all(j),dihedatom2all(j),
     $         dihedatom3all(j),dihedatom4all(j)
        enddo

        deallocate(indx)

      endif

c write impropers in sorted order via improatom2all, ignore ties

      if (nimpropers > 0) then

        allocate(indx(nimpropers))
        call indexx(nimpropers,improatom2all,indx)

        write (1,*)
        write (1,900) 'Impropers'
        write (1,*)

        do i = 1,nimpropers
          j = indx(i)
          write (1,905) i,improtypeall(j),
     $         improatom1all(j),improatom2all(j),
     $         improatom3all(j),improatom4all(j)
        enddo

        deallocate(indx)

      endif

      close (1)

      end

c ****************
c Subroutines
c ****************

c index sorting routine from Numerical Recipes

      subroutine indexx(n,arr,indx)
      integer n,indx(n),m,nstack
      integer arr(n)
      parameter (m=7,nstack=50)
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      real a

      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.m)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.nstack)pause 'nstack too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1

      end
