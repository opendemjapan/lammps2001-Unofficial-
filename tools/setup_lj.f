c Create LAMMPS data file (version 2001) for 3-d LJ simulation of monomers
c sets up mixture of LJ atom types
c Syntax: setup_lj data.file

      program setup_lj
      implicit real*8 (a-h,o-z)

      real*8, allocatable :: x(:,:)
      integer, allocatable :: target(:),type(:)
      character*80 datafile
      real*8 random
 900  format(a)
c lo-precision output
 901  format(2f10.6,a)
 902  format(i3,f5.1)
 903  format(i6,i4,i4,f4.1,3f10.4)
c hi-precision output
c 901  format(2f,a)
c 902  format(i3,f5.1)
c 903  format(i6,i4,i4,f4.1,3f)


      real*8 epsilon,sigma,mass
c liquid Argon Lennard-Jones values
      parameter (epsilon=120.0,sigma=3.4,mass=39.948)
c unitless LJ
c      parameter (epsilon=1.0d0,sigma=1.0d0,mass=1.0d0)

c prompts

      if (iargc() /= 1) then
        write (6,*) 'Syntax: setup_lj data.file'
        call exit(0)
      else
        call getarg(1,datafile)
      endif

      write (6,*) '(0) FCC Lattice or (1) Number of atoms'
      read (5,*) isetup
      if (isetup == 0) then
        write (6,*) 'Lattice: nx,ny,nz'
        read (5,*) nx,ny,nz
        natoms = 4*nx*ny*nz
      else
        write (6,*) 'Number of atoms: n'
        read (5,*) natoms
      endif

      write (6,*) 'Rhostar'
      read (5,*) rhostar

      write (6,*) 'Random # seed (8 digit max)'
      read (5,*) iseed

      write (6,*) 'Number of atom types'
      read (5,*) ntypes

      allocate(target(ntypes))

      if (ntypes > 1) then
        nsum = 0
        do i = 1,ntypes
          write (6,*) '# of type',i
          read (5,*) target(i)
          nsum = nsum + target(i)
        enddo
        if (nsum /= natoms) then
          write (6,*) 'Total of types <> total # of atoms'
          call exit(0)
        endif
      else
        target(1) = natoms
      endif

c allocate memory

      allocate(x(3,natoms))
      allocate(type(natoms))

c setup rectangular or cubic box

      volume = natoms * sigma*sigma*sigma / rhostar

      if (isetup == 0) then
        alattice = (volume/(nx*ny*nz))**(1.0d0/3.0d0)
        xboundlo = -(nx*alattice)/2.0
        xboundhi = -xboundlo
        yboundlo = -(ny*alattice)/2.0
        yboundhi = -yboundlo
        zboundlo = -(nz*alattice)/2.0
        zboundhi = -zboundlo
      else
        xboundlo = -(volume**(1.0d0/3.0d0))/2.0
        xboundhi = -xboundlo
        yboundlo = xboundlo
        yboundhi = xboundhi
        zboundlo = xboundlo
        zboundhi = xboundhi
      endif

c initial fcc lattice or random points

      if (isetup.eq.0) then

        m = 0
        do k = 1,nz*2
          do j = 1,ny*2
            do i = 1,nx*2
              if (mod(i+j+k,2).eq.1) then
                m = m + 1
                x(1,m) = xboundlo + (i-1)*alattice/2.0
                x(2,m) = yboundlo + (j-1)*alattice/2.0
                x(3,m) = zboundlo + (k-1)*alattice/2.0
              endif
            enddo
          enddo
        enddo
        
      else

        do m = 1,natoms
          x(1,m) = xboundlo + random(iseed)*(xboundhi-xboundlo)
          x(2,m) = yboundlo + random(iseed)*(yboundhi-yboundlo)
          x(3,m) = zboundlo + random(iseed)*(zboundhi-zboundlo)
        enddo

      endif

c set types of atoms randomly according to user-specified targets

      m = 0
      do i = 1,ntypes
        do j = 1,target(i)
          m = m + 1
          type(m) = i
        enddo
      enddo

      do i = 1,natoms
        j = natoms*random(iseed) + 1
        k = type(i)
        type(i) = type(j)
        type(j) = k
      enddo

c write out file

      open (1,file=datafile)

      write (1,900) 'LAMMPS Description'
      write (1,*)

      write (1,*) natoms,' atoms'
      write (1,*) 0,' bonds'
      write (1,*) 0,' angles'
      write (1,*) 0,' dihedrals'
      write (1,*) 0,' impropers'
      write (1,*)

      write (1,*) ntypes,' atom types'
      write (1,*)

      write (1,901) xboundlo,xboundhi,' xlo xhi'
      write (1,901) yboundlo,yboundhi,' ylo yhi'
      write (1,901) zboundlo,zboundhi,' zlo zhi'

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
        write (1,903) i,0,type(i),0.0,x(1,i),x(2,i),x(3,i)
      enddo

      close (1)

      end


c ************
c Subroutines
c ************

c RNG from Numerical Recipes
      
      real*8 function random(iseed)
      real*8 aa,mm,sseed
      parameter (aa=16807.0D0,mm=2147483647.0D0)
      
      sseed = iseed
      sseed = mod(aa*sseed,mm)
      random = sseed/mm
      iseed = sseed

      return
      end
