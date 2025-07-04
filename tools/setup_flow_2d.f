c Create LAMMPS data file (version 2001) for 2-d LJ simulation of channel flow
c   all flow atoms in channel are type 1
c   lower wall atoms are type 2, upper wall atoms are type 3
c   virtual wall atoms for "warm" wall (if enabled) are type 4,5
c   flow is in x-direction, channel walls are +/- y direction
c Syntax: setup_flow_2d data.file

      program setup_flow_2d

      real*8, allocatable :: x(:,:)
      integer, allocatable :: type(:),bondatom(:,:)
      character*80 datafile
      real*8 epsilon,sigma,mass
      parameter (epsilon=1.0,sigma=1.0,mass=1.0)
 900  format(a)
 901  format(2f10.6,a)
 902  format(i3,f5.1)
 903  format(i4,i4,i4,f4.1,2f10.4,f4.1)
 904  format(i3,i3,2i6)

      if (iargc() /= 1) then
        write (6,*) 'Syntax: setup_flow_2d data.file'
        call exit(0)
      else
        call getarg(1,datafile)
      endif

c user prompts
c lattice is nx by ny, but add an extra layer of atoms in vertical y
c   since system will not be periodic in y-dimension

      write (6,*) '(0) Square or (1) Hex lattice'
      read (5,*) isetup
      write (6,*) '2-d Lattice size: nx,ny'
      read (5,*) nx,ny

      write (6,*) '# of layers of wall atoms on top/bottom of box'
      read (5,*) nwall
      write (6,*) '(0) Cold or (1) Warm wall'
      read (5,*) iwarm

      if (isetup == 0) then
        if (2*nwall > ny+1) then
          write (6,*) 'Too many wall layers'
          call exit(0)
        endif
      else
        if (2*nwall > 2*ny+1) then
          write (6,*) 'Too many wall layers'
          call exit(0)
        endif
      endif

      write (6,*) 'Rhostar'
      read (5,*) rhostar

c allocate memory for all atoms

      if (isetup == 0) then
        natoms = nx*ny + nx
      else
        natoms = 2*nx*ny + nx
      endif
      if (iwarm == 1) natoms = 2*natoms

      allocate(x(2,natoms))
      allocate(type(natoms))
      allocate(bondatom(2,natoms))

c setup rectangular or square box
c set volume (via rhostar) as if extra layer of atoms weren't there
      
      volume = (natoms-nx) * sigma*sigma / rhostar
      alattice = (volume/(nx*ny))**(1.0/2.0)

      xboundlo = -(nx*alattice)/2.0
      xboundhi = -xboundlo
      yboundlo = -(ny*alattice)/2.0
      yboundhi = -yboundlo

c initial lattice - flow atoms are type 1, wall atoms are type 2,3

      if (isetup == 0) then

        m = 0
        do j = 1,ny+1
          do i = 1,nx
            m = m + 1
            x(1,m) = xboundlo + (i-1)*alattice
            x(2,m) = yboundlo + (j-1)*alattice
            if (j <= nwall) then
              type(m) = 2
            else if (ny+1-j < nwall) then
              type(m) = 3
            else
              type(m) = 1
            endif
          enddo
        enddo

      else

        m = 0
        do j = 1,ny*2+1
          do i = 1,nx*2
            if (mod(i+j,2) == 0) then
              m = m + 1
              x(1,m) = xboundlo + (i-1)*alattice/2.0
              x(2,m) = yboundlo + (j-1)*alattice/2.0
              if (j <= nwall) then
                type(m) = 2
              else if (ny*2+1-j < nwall) then
                type(m) = 3
              else
                type(m) = 1
              endif
            endif
          enddo
        enddo
        
      endif

      natoms = m

c set up cold or warm wall
c warm wall atoms do not contribute to volume

      if (iwarm == 0) then
        ntypes = 3
        nbonds = 0
      else
        m = natoms
        nbonds = 0
        do i = 1,natoms
          if (type(i) == 2 .or. type(i) == 3) then
            m = m + 1
            type(m) = type(i) + 2
            x(1,m) = x(1,i)
            x(2,m) = x(2,i)
            nbonds = nbonds + 1
            bondatom(1,nbonds) = i
            bondatom(2,nbonds) = m
          endif
        enddo
        natoms = m
        ntypes = 5
        nbondtypes = 1
      endif

c write out file

      open (1,file=datafile)

      write (1,900) 'LAMMPS Description'
      write (1,*)

      write (1,*) natoms,' atoms'
      write (1,*) nbonds,' bonds'
      write (1,*) 0,' angles'
      write (1,*) 0,' dihedrals'
      write (1,*) 0,' impropers'
      write (1,*)

      write (1,*) ntypes,' atom types'
      if (nbonds > 0) write (1,*) nbondtypes,' bond types'
      write (1,*)

      write (1,901) xboundlo,xboundhi,' xlo xhi'
      write (1,901) yboundlo,yboundhi,' ylo yhi'

      write (1,*)
      write (1,900) 'Masses'
      write (1,*)

      do i = 1,ntypes
        write (1,902) i,1.0
      enddo

      write (1,*)
      write (1,900) 'Atoms'
      write (1,*)

      do i = 1,natoms
        write (1,903) i,0,type(i),0.0,x(1,i),x(2,i),0.0
      enddo

      if (nbonds > 0) then

        write (1,*)
        write (1,900) 'Bonds'
        write (1,*)

        do i = 1,nbonds
          write (1,904) i,1,bondatom(1,i),bondatom(2,i)
        enddo

      endif

      close (1)

      end
