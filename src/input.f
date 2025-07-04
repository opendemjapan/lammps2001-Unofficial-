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

c ------------------------------------------------------------------------
c read all info from text input file one line at a time

      subroutine input(iopflag)
      use global
      use mpi
      implicit none

c argument variables

      integer iopflag

c local variables

      logical match
      character*200 str
      integer ierror,m

 900  format (a)
      
      do
        
c read one line from file and broadcast it

        if (node == 0) read (5,900,end=999,err=999) str
        call mpi_bcast(str,200,mpi_character,0,mpi_comm_world,ierror)

c blank line

        if (len_trim(str) == 0) cycle

c comment line

        if (match('#',str,m)) then
          if (node == 0) write (1,900) trim(str)
          cycle
        endif

c if it appears, "units" must be 1st command in input script

        if (match('units',str,m)) then
          if (firstflag == 1) call err(
     $         'Units command must be first command in file')
          call in_units(str,m)
          cycle
        endif

        firstflag = 1

c commands that must occur BEFORE read_data and read_restart

        if (match('extra memory',str,m)) then
          call in_extra_memory(str,m)
        else if (match('dimension',str,m)) then
          call in_dimension(str,m)
        else if (match('processor grid',str,m)) then
          call in_processor_grid(str,m)
        else if (match('periodicity',str,m)) then
          call in_periodicity(str,m)
        else if (match('slab volume',str,m)) then
          call in_slab_volume(str,m)
        else if (match('newton flag',str,m)) then
          call in_newton_flag(str,m)
        else if (match('true flag',str,m)) then
          call in_true_flag(str,m)
        else if (match('maximum cutoff',str,m)) then
          call in_maximum_cutoff(str,m)
        else if (match('mixing style',str,m)) then
          call in_mixing_style(str,m)
        else if (match('restart version',str,m)) then
          call in_restart_version(str,m)

c commands that can appear BEFORE or AFTER read_data and read_restart

        else if (match('neighbor',str,m)) then
          call in_neighbor(str,m)

        else if (match('nonbond style',str,m)) then
          call in_nonbond_style(str,m)
        else if (match('coulomb style',str,m)) then
          call in_coulomb_style(str,m)

        else if (match('bond style',str,m)) then
          call in_bond_style(str,m)
        else if (match('angle style',str,m)) then
          call in_angle_style(str,m)
        else if (match('dihedral style',str,m)) then
          call in_dihedral_style(str,m)
        else if (match('improper style',str,m)) then
          call in_improper_style(str,m)

c read_data and read_restart commands

        else if (match('read data',str,m)) then
          iopflag = 1
          call in_read_data(str,m)
          return

        else if (match('read restart',str,m)) then
          iopflag = 2
          call in_read_restart(str,m)
          return

c commands that must occur AFTER read_data and read_restart

c velocity creation

        else if (match('rotation zero',str,m)) then
          call in_rotation_zero(str,m)
        else if (match('create group',str,m)) then
          call in_create_group(str,m)
        else if (match('create temp',str,m)) then
          iopflag = 3
          call in_create_temp(str,m)
          return

c force field parameters

        else if (match('nonbond coeff',str,m)) then
          call in_nonbond_coeff(str,m)
        else if (match('special bonds',str,m)) then
          call in_special_bonds(str,m)

        else if (match('pppm mesh',str,m)) then
          call in_pppm_mesh(str,m)
        else if (match('pppm order',str,m)) then
          call in_pppm_order(str,m)
        else if (match('dielectric',str,m)) then
          call in_dielectric(str,m)

        else if (match('bond coeff',str,m)) then
          call in_bond_coeff(str,m)
        else if (match('angle coeff',str,m)) then
          call in_angle_coeff(str,m)
        else if (match('dihedral coeff',str,m)) then
          call in_dihedral_coeff(str,m)
        else if (match('improper coeff',str,m)) then
          call in_improper_coeff(str,m)

c constraints

        else if (match('fix style',str,m)) then
          call in_fix_style(str,m)
        else if (match('assign fix',str,m)) then
          iopflag = 4
          call in_assign_fix(str,m)
          return

c ensemble control

        else if (match('temp control',str,m)) then
          call in_temp_control(str,m)
        else if (match('press control',str,m)) then
          call in_press_control(str,m)
        else if (match('volume control',str,m)) then
          call in_volume_control(str,m)

c output control

        else if (match('thermo flag',str,m)) then
          call in_thermo_flag(str,m)
        else if (match('thermo style',str,m)) then
          call in_thermo_style(str,m)
        else if (match('dump atoms',str,m)) then
          call in_dump_atoms(str,m)
        else if (match('dump velocities',str,m)) then
          call in_dump_velocities(str,m)
        else if (match('dump forces',str,m)) then
          call in_dump_forces(str,m)
        else if (match('restart',str,m)) then
          call in_restart(str,m)
        else if (match('diagnostic',str,m)) then
          call in_diagnostic(str,m)

c integrator settings

        else if (match('timestep',str,m)) then
          call in_timestep(str,m)
        else if (match('respa',str,m)) then
          call in_respa(str,m)
        else if (match('reset timestep',str,m)) then
          call in_reset_timestep(str,m)

c minimizer settings

        else if (match('min style',str,m)) then
          call in_min_style(str,m)
        else if (match('min flag',str,m)) then
          call in_min_flag(str,m)

c dynamics or minimization

        else if (match('run',str,m)) then
          iopflag = 5
          call in_run(str,m)
          return

        else if (match('minimize',str,m)) then
          iopflag = 6
          call in_minimize(str,m)
          return

c invalid command

        else
          call err('Command: '//trim(str))
        endif
        
      enddo
      
 999  continue
      
      str = 'All Done'
      call mpi_bcast(str,200,mpi_character,0,mpi_comm_world,ierror)
      call err(str(1:8))
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_units(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer n
      character*16 work
      logical match

      call strread(str,m,work)
      if (len_trim(work) == 0) call err('Bad units param')

      if (match('real',work,n)) then

        units = 0

c set constants

        boltz = 0.001987191
        dtfactor = 48.88821
        pfactor = 68589.796
        efactor = 332.0636

c set default cutoffs

        skin = 2.0
        cutlj = 10.0
        cutcoul = 10.0

      else if (match('lj',work,n)) then

        units = 1

c set constants

        boltz = 1.0
        dtfactor = 1.0
        pfactor = 1.0
        efactor = 1.0

c set default cutoffs

        skin = 0.3
        cutlj = 2.5
        cutcoul = 2.5

      else
        call err('Bad units param')
      endif

      if (node == 0) then
        if (units == 0) write (1,*) 'Units real'
        if (units == 1) write (1,*) 'Units lj'
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_extra_memory(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      real*8 realread

      call input_check(str,0)
      extra_own = realread(str,m)
      extra_ghost = realread(str,m)
      extra_neigh = realread(str,m)
      extra_buf = realread(str,m)
      if (node == 0)
     $     write (1,*) 'Extra memory',
     $     extra_own,extra_ghost,extra_neigh,extra_buf

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dimension(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,0)
      idimension = intread(str,m)
      if (idimension < 2 .or. idimension > 3)
     $     call err('Bad dimension parameter')
      if (node == 0)
     $     write (1,*) 'Dimension',idimension

      return
      end

c ------------------------------------------------------------------------

      subroutine in_processor_grid(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,0)
      pgrid(1) = intread(str,m)
      pgrid(2) = intread(str,m)
      pgrid(3) = intread(str,m)
      if (pgrid(1) <= 0.or.pgrid(2) <= 0.or.pgrid(3) <= 0)
     $     call err('Bad processor grid parameter')
      if (pgrid(1)*pgrid(2)*pgrid(3) /= nprocs)
     $     call err('Specified processor grid not equal Nprocs')
      if (node == 0)
     $     write (1,*) 'Processor grid',pgrid(1),pgrid(2),pgrid(3)

      return
      end

c ------------------------------------------------------------------------

      subroutine in_periodicity(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,0)
      perflagx = intread(str,m)
      perflagy = intread(str,m)
      perflagz = intread(str,m)
      if (perflagx < 0.or.perflagx > 1.or.
     $     perflagy < 0.or.perflagy > 1.or.
     $     perflagz < 0.or.perflagz > 1)
     $     call err('Bad periodicity parameter')
      if (perflagx+perflagy+perflagz > 0) then
        nonperiodic = .TRUE.
      else
        nonperiodic = .FALSE.
      endif
      if (node == 0)
     $     write (1,*) 'Periodicity',perflagx,perflagy,perflagz

      return
      end

c ------------------------------------------------------------------------

      subroutine in_slab_volume(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread
      real*8 realread

      call input_check(str,0)
      slab_volfactor = realread(str,m)
      if (slab_volfactor < 2.0)
     $     call err('Bad slab volume parameter')
      slabflag = 1
      if (node == 0)
     $     write (1,*) 'Slab volume',slab_volfactor

      return
      end

c ------------------------------------------------------------------------

      subroutine in_newton_flag(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,0)
      newton = intread(str,m)
      if (newton < 0.or.newton > 3)
     $     call err('Bad newton flag parameter')
      newton_bond = mod(newton,2)
      newton_nonbond = newton/2
      if (node == 0)
     $     write (1,*) 'Newton flag',newton

      return
      end

c ------------------------------------------------------------------------

      subroutine in_true_flag(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,0)
      trueflag = intread(str,m)
      if (trueflag < 0.or.trueflag > 3)
     $     call err('Bad true flag parameter')
      if (node == 0)
     $     write (1,*) 'True flag',trueflag

      return
      end

c ------------------------------------------------------------------------

      subroutine in_maximum_cutoff(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      real*8 realread

      call input_check(str,0)
      cutmax = realread(str,m)
      if (cutmax < 0.0)
     $     call err('Bad maximum cutoff parameter')
      if (node == 0)
     $     write (1,*) 'Maximum cutoff',cutmax

      return
      end

c ------------------------------------------------------------------------

      subroutine in_mixing_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm

      logical match
      character*16 work

      call input_check(str,0)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad mixing style parameter')
      mixflag = 1
      if (match('geometric',work,mm)) then
        mixstyle = 1
      else if (match('arithmetic',work,mm)) then
        mixstyle = 2
      else if (match('sixthpower',work,mm)) then
        mixstyle = 3
      else
        call err('Bad mixing style parameter')
      endif
      if (node == 0) then
        if (mixstyle == 1) then
          write (1,*) 'Mixing style geometric'
        else if (mixstyle == 2) then
          write (1,*) 'Mixing style arithmetic'
        else if (mixstyle == 3) then
          write (1,*) 'Mixing style sixthpower'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_restart_version(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,0)
      restart_version = intread(str,m)
      if (restart_version < 5)
     $     call err('Bad restart version parameter')
      if (node == 0)
     $     write (1,*) 'Restart version',restart_version

      return
      end

c ------------------------------------------------------------------------

      subroutine in_neighbor(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread
      real*8 realread

      skin = realread(str,m)
      neighstyle = intread(str,m)
      neighfreq = intread(str,m)
      neighdelay = intread(str,m)
      neightrigger = intread(str,m)
      if (skin < 0.0.or.
     $     (neighstyle < 0.or.neighstyle > 1).or.
     $     neighfreq <= 0.or.neighdelay < 0.or.
     $     (neightrigger < 0.or.neightrigger > 1))
     $     call err('Bad neighbor parameter')
      if (node == 0)
     $     write (1,*) 'Neighbor',skin,neighstyle,
     $     neighfreq,neighdelay,neightrigger

      return
      end

c ------------------------------------------------------------------------

      subroutine in_nonbond_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,i,j,intread
      real*8 realread

      logical match
      character*16 work

c if data file already read in:
c  set cutoffs in noncoeff3/noncoeff4 for all type pairs
c  else cutoffs will be set for all type pairs in read_data
c  or subsequent nonbond coeff commands

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad nonbond style parameter')

      if (match('none',work,mm)) then
        nonstyle = 0

      else if (match('lj/cutoff',work,mm)) then
        nonstyle = 1
        cutlj = realread(str,m)
        offsetflag = intread(str,m)
        if (cutlj <= 0.0.or.offsetflag < 0.or.offsetflag > 1)
     $       call err('Bad nonbond style parameter')
        if (readflag == 1) then
          do i = 1,ntypes
            do j = i,ntypes
              noncoeff3(i,j) = cutlj
            enddo
          enddo
        endif
        if (mixflag == 0) mixstyle = 1

      else if (match('lj/smooth',work,mm)) then
        nonstyle = 2
        cutljinterior = realread(str,m)
        cutlj = realread(str,m)
        if (cutljinterior <= 0.0.or.cutljinterior >= cutlj)
     $       call err('Bad nonbond style parameter')
        if (readflag == 1) then
          do i = 1,ntypes
            do j = i,ntypes
              noncoeff3(i,j) = cutljinterior
              noncoeff4(i,j) = cutlj
            enddo
          enddo
        endif
        if (mixflag == 0) mixstyle = 1

      else if (match('lj/shift',work,mm)) then
        nonstyle = 3
        cutlj = realread(str,m)
        offsetflag = intread(str,m)
        if (cutlj <= 0.0.or.offsetflag < 0.or.offsetflag > 1)
     $       call err('Bad nonbond style parameter')
        if (readflag == 1) then
          do i = 1,ntypes
            do j = i,ntypes
              noncoeff4(i,j) = cutlj
            enddo
          enddo
        endif
        if (mixflag == 0) mixstyle = 1

      else if (match('soft',work,mm)) then
        nonstyle = 4
        cutlj = realread(str,m)
        if (cutlj <= 0.0) call err('Bad nonbond style parameter')
        if (readflag == 1) then
          do i = 1,ntypes
            do j = i,ntypes
              noncoeff3(i,j) = cutlj
            enddo
          enddo
        endif
        mixstyle = 1

      else if (match('class2/cutoff',work,mm)) then
        nonstyle = 5
        cutlj = realread(str,m)
        offsetflag = intread(str,m)
        if (cutlj <= 0.0.or.offsetflag < 0.or.offsetflag > 1)
     $       call err('Bad nonbond style parameter')
        if (readflag == 1) then
          do i = 1,ntypes
            do j = i,ntypes
              noncoeff3(i,j) = cutlj
            enddo
          enddo
        endif
        if (mixflag == 0) mixstyle = 3

      else if (match('lj/charmm',work,mm)) then
        nonstyle = 6
        cutljinterior = realread(str,m)
        cutlj = realread(str,m)
        cutljint_sq = cutljinterior*cutljinterior
        cutlj_sq = cutlj*cutlj
        ch_denom_lj = (cutlj_sq - cutljint_sq)**3
        if (cutlj <= 0.0 .or. cutljinterior <= 0)
     $       call err('Bad nonbond style parameter')
        mixstyle = 2

      else
        call err('Bad nonbond style parameter')
      endif

      if (node == 0) then
        if (nonstyle == 0) then
          write (1,*) 'Nonbond style none'
        else if (nonstyle == 1) then
          write (1,*) 'Nonbond style lj/cutoff',
     $         cutlj,offsetflag
        else if (nonstyle == 2) then
          write (1,*) 'Nonbond style lj/smooth',
     $         cutljinterior,cutlj
        else if (nonstyle == 3) then
          write (1,*) 'Nonbond style lj/shift',
     $         cutlj,offsetflag
        else if (nonstyle == 4) then
          write (1,*) 'Nonbond style soft',
     $         cutlj
        else if (nonstyle == 5) then
          write (1,*) 'Nonbond style class2/cutoff',
     $         cutlj,offsetflag
        else if (nonstyle == 6) then
          write (1,*) 'Nonbond style lj/charmm',
     $         cutljinterior,cutlj
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_coulomb_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm
      real*8 realread

      logical match
      character*16 work

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad coulomb style parameter')
      if (match('none',work,mm)) then
        coulstyle = 0
        cutcoul = 0.0
      else if (match('cutoff',work,mm)) then
        coulstyle = 1
        cutcoul = realread(str,m)
        if (cutcoul <= 0.0)
     $       call err('Bad coulomb style parameter')
      else if (match('smooth',work,mm)) then
        coulstyle = 2
        cutcoulint = realread(str,m)
        cutcoul = realread(str,m)
        if (cutcoulint <= 0.0.or.cutcoulint >= cutcoul)
     $       call err('Bad coulomb style parameter')
      else if (match('ewald',work,mm)) then
        coulstyle = 3
        cutcoul = realread(str,m)
        long_prec = realread(str,m)
        if (cutcoul <= 0.0.or.
     $       long_prec <= 0.0.or.long_prec > 1.0)
     $       call err('Bad coulomb style parameter')
      else if (match('pppm',work,mm)) then
        coulstyle = 4
        cutcoul = realread(str,m)
        long_prec = realread(str,m)
        if (cutcoul <= 0.0.or.
     $       long_prec <= 0.0.or.long_prec > 1.0)
     $       call err('Bad coulomb style parameter')
      else if (match('charmm/switch',work,mm)) then
        coulstyle = 5
        cutcoulint = realread(str,m)
        cutcoul = realread(str,m)
        ch_denom_coul = (cutcoul*cutcoul - cutcoulint*cutcoulint)**3
        if (cutcoulint <= 0.0.or.cutcoulint >= cutcoul)
     $       call err('Bad coulomb style parameter')
      else if (match('debye',work,mm)) then
        coulstyle = 6
        cutcoul = realread(str,m)
        kappa = realread(str,m)
        if (cutcoul <= 0.0 .or. kappa < 0.0)
     $       call err('Bad coulomb style parameter')
      else
        call err('Bad coulomb style parameter')
      endif

      if (node == 0) then
        if (coulstyle == 0) then
          write (1,*) 'Coulomb style none'
        else if (coulstyle == 1) then
          write (1,*) 'Coulomb style cutoff',cutcoul
        else if (coulstyle == 2) then
          write (1,*) 'Coulomb style smooth',cutcoulint,cutcoul
        else if (coulstyle == 3) then
          write (1,*) 'Coulomb style ewald',cutcoul,long_prec
        else if (coulstyle == 4) then
          write (1,*) 'Coulomb style pppm',cutcoul,long_prec
        else if (coulstyle == 5) then
          write (1,*) 'Coulomb style charmm/switch',
     $         cutcoulint,cutcoul
        else if (coulstyle == 6) then
          write (1,*) 'Coulomb style debye',cutcoul,kappa
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_bond_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,i
      logical match
      character*16 work

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad bond style parameter')
      if (match('none',work,mm)) then
        bondstyle = 0
      else if (match('harmonic',work,mm)) then
        bondstyle = 1
      else if (match('fene/standard',work,mm)) then
        bondstyle = 2
      else if (match('fene/shift',work,mm)) then
        bondstyle = 3
      else if (match('nonlinear',work,mm)) then
        bondstyle = 4
      else if (match('class2',work,mm)) then
        bondstyle = 5
      else
        call err('Bad bond style parameter')
      endif
      if (node == 0) then
        if (bondstyle == 0) then
          write (1,*) 'Bond style none'
        else if (bondstyle == 1) then
          write (1,*) 'Bond style harmonic'
        else if (bondstyle == 2) then
          write (1,*) 'Bond style fene/standard'
        else if (bondstyle == 3) then
          write (1,*) 'Bond style fene/shift'
        else if (bondstyle == 4) then
          write (1,*) 'Bond style nonlinear'
        else if (bondstyle == 5) then
          write (1,*) 'Bond style class2'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_angle_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,i
      logical match
      character*16 work

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad angle style parameter')
      if (match('none',work,mm)) then
        anglestyle = 0
      else if (match('harmonic',work,mm)) then
        anglestyle = 1
      else if (match('class2',work,mm)) then
        anglestyle = 2
      else if (match('charmm',work,mm)) then
        anglestyle = 3
      else if (match('cosine',work,mm)) then
        anglestyle = 4
      else
        call err('Bad angle style parameter')
      endif

      if (node == 0) then
        if (anglestyle == 0) then
          write (1,*) 'Angle style none'
        else if (anglestyle == 1) then
          write (1,*) 'Angle style harmonic'
        else if (anglestyle == 2) then
          write (1,*) 'Angle style class2'
        else if (anglestyle == 3) then
          write (1,*) 'Angle style charmm'
        else if (anglestyle == 3) then
          write (1,*) 'Angle style cosine'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dihedral_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,i
      logical match
      character*16 work

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad dihedral style parameter')
      if (match('none',work,mm)) then
        dihedstyle = 0
      else if (match('harmonic',work,mm)) then
        dihedstyle = 1
      else if (match('class2',work,mm)) then
        dihedstyle = 2
      else if (match('multiharmonic',work,mm)) then
        dihedstyle = 3
      else if (match('charmm',work,mm)) then
        dihedstyle = 4
      else
        call err('Bad dihedral style parameter')
      endif
      if (node == 0) then
        if (dihedstyle == 0) then
          write (1,*) 'Dihedral style none'
        else if (dihedstyle == 1) then
          write (1,*) 'Dihedral style harmonic'
        else if (dihedstyle == 2) then
          write (1,*) 'Dihedral style class2'
        else if (dihedstyle == 3) then
          write (1,*) 'Dihedral style multiharmonic'
        else if (dihedstyle == 4) then
          write (1,*) 'Dihedral style charmm'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_improper_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,i
      logical match
      character*16 work

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad improper style parameter')
      if (match('none',work,mm)) then
        improstyle = 0
      else if (match('harmonic',work,mm)) then
        improstyle = 1
      else if (match('cvff',work,mm)) then
        improstyle = 2
      else if (match('class2',work,mm)) then
        improstyle = 3
      else
        call err('Bad improper style parameter')
      endif
      if (node == 0) then
        if (improstyle == 0) then
          write (1,*) 'Improper style none'
        else if (improstyle == 1) then
          write (1,*) 'Improper style harmonic'
        else if (improstyle == 2) then
          write (1,*) 'Improper style cvff'
        else if (improstyle == 3) then
          write (1,*) 'Improper style class2'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_read_data(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

      call input_check(str,0)
      call strread(str,m,datafile)
      if (len_trim(datafile) == 0)
     $     call err('Bad read data parameter')
      if (node == 0)
     $     write (1,*) 'Read data ',trim(datafile)

      return
      end

c ------------------------------------------------------------------------

      subroutine in_read_restart(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

      call input_check(str,0)
      call strread(str,m,restart_in)
      if (len_trim(restart_in) == 0)
     $     call err('Bad read restart parameter')
      if (node == 0)
     $     write (1,*) 'Read restart ',trim(restart_in)

      return
      end

c ------------------------------------------------------------------------

      subroutine in_timestep(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      real*8 realread

      call input_check(str,1)
      dt = realread(str,m)
      if (dt <= 0.0) call err('Bad timestep parameter')
      if (node == 0)
     $     write (1,*) 'Timestep',dt
      dt = dt/dtfactor

      return
      end

c ------------------------------------------------------------------------

      subroutine in_respa(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      nstretch = intread(str,m)
      nintra = intread(str,m)
      nshort = intread(str,m)
      if (nstretch <= 0 .or. nintra <= 0 .or. nshort <= 0)
     $     call err('Bad respa parameter')
      nrespa = 1
      if (nstretch == 1 .and. nintra == 1 .and. nshort == 1) nrespa = 0
      if (node == 0)
     $     write (1,*) 'Respa',nstretch,nintra,nshort

      return
      end

c ------------------------------------------------------------------------

      subroutine in_reset_timestep(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      ntimestep = intread(str,m)
      if (node == 0)
     $     write (1,*) 'Reset timestep',ntimestep

      return
      end

c ------------------------------------------------------------------------

      subroutine in_special_bonds(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      real*8 realread
      character*16 work
      logical match
      integer mm

      call input_check(str,1)
      mm = m
      call strread(str,mm,work)
      if (len_trim(work) == 0)
     $     call err('Bad special bonds parameter')
      if (match('charmm',work,mm)) then
        amberflag = 0
        special(1) = 0.0
        special(2) = 0.0
        special(3) = 0.0
        if (node == 0) write (1,*) 'Special bonds charmm'
      else if (match('amber',work,mm)) then
        amberflag = 1
        special(1) = 0.0
        special(2) = 0.0
        special(3) = 1.0d0/1.2d0
        if (node == 0) write (1,*) 'Special bonds amber'
      else
        amberflag = 0
        special(1) = realread(str,m)
        special(2) = realread(str,m)
        special(3) = realread(str,m)
        if (special(1) < 0.0.or.special(1) > 1.0.or.
     $       special(2) < 0.0.or.special(2) > 1.0.or.
     $       special(3) < 0.0.or.special(3) > 1.0)
     $       call err('Bad special bonds parameter')
        if (node == 0) write (1,*) 'Special bonds',
     $       special(1),special(2),special(3)
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_thermo_flag(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      nthermo = intread(str,m)
      if (nthermo < 0)
     $     call err('Bad thermo flag parameter')
      if (node == 0)
     $     write (1,*) 'Thermo flag',nthermo

      return
      end

c ------------------------------------------------------------------------

      subroutine in_thermo_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      thermostyle = intread(str,m)
      if (thermostyle < 0.or.thermostyle > 2)
     $     call err('Bad thermo style parameter')
      if (node == 0)
     $     write (1,*) 'Thermo style',thermostyle
      
      return
      end


c ------------------------------------------------------------------------

      subroutine in_dump_atoms(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      ndumpatom = intread(str,m)
      if (ndumpatom < 0) call err('Bad dump atoms parameter')

      if (ndumpatom == 0) then
        if (dumpatomfileflag > 0) call dumpatomclose
        dumpatomfileflag = 0
        if (node == 0) write (1,*) 'Dump atoms',ndumpatom
      else
        call strread(str,m,dumpatomfile)
        if (len_trim(dumpatomfile) == 0)
     $       call err('Bad dump atoms param')
        if (dumpatomfileflag > 0) call dumpatomclose
        call dumpatomopen
        dumpatomfileflag = 1
        if (node == 0) write (1,*) 'Dump atoms',ndumpatom,' ',
     $       trim(dumpatomfile)
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dump_velocities(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      ndumpvel = intread(str,m)
      if (ndumpvel < 0) call err('Bad dump velocities parameter')

      if (ndumpvel == 0) then
        if (dumpvelfileflag > 0) call dumpvelclose
        dumpvelfileflag = 0
        if (node == 0) write (1,*) 'Dump velocities',ndumpvel
      else
        call strread(str,m,dumpvelfile)
        if (len_trim(dumpvelfile) == 0)
     $       call err('Bad dump velocities param')
        if (dumpvelfileflag > 0) call dumpvelclose
        call dumpvelopen
        dumpvelfileflag = 1
        if (node == 0) write (1,*) 'Dump velocities',ndumpvel,' ',
     $       trim(dumpvelfile)
      endif
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_dump_forces(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      ndumpforce = intread(str,m)
      if (ndumpforce < 0) call err('Bad dump forces parameter')

      if (ndumpforce == 0) then
        if (dumpforcefileflag > 0) call dumpforceclose
        dumpforcefileflag = 0
        if (node == 0) write (1,*) 'Dump forces',ndumpforce
      else
        call strread(str,m,dumpforcefile)
        if (len_trim(dumpforcefile) == 0)
     $       call err('Bad dump forces param')
        if (dumpforcefileflag > 0) call dumpforceclose
        call dumpforceopen
        dumpforcefileflag = 1
        if (node == 0) write (1,*) 'Dump forces',ndumpforce,' ',
     $       trim(dumpforcefile)
      endif
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_restart(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      nrestart = intread(str,m)
      if (nrestart < 0) call err('Bad restart param')

      if (nrestart == 0) then
        if (node == 0) write (1,*) 'Restart',nrestart
      else
        restartstyle = intread(str,m)
        if (restartstyle < 1 .or. restartstyle > 2)
     $       call err('Bad restart param')
        if (restartstyle == 1) then
          call strread(str,m,restart_out)
          if (len_trim(restart_out) == 0)
     $         call err('Bad restart param')
        else
          call strread(str,m,restart_out1)
          call strread(str,m,restart_out2)
          if (len_trim(restart_out1) == 0 .or.
     $         len_trim(restart_out2) == 0)
     $         call err('Bad restart param')
          restartlast = 2
        endif
        if (node == 0) then
          if (restartstyle == 1) then
            write (1,*) 'Restart',nrestart,restartstyle,
     $           trim(restart_out)
          else
            write (1,*) 'Restart',nrestart,restartstyle,
     $           trim(restart_out1),' ',trim(restart_out2)
          endif
        endif
      endif
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_diagnostic(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer idiag,intread,i
      real*8 realread
      character*16 work

      call input_check(str,1)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad diagnostic parameter')

      call diag_lookup(trim(work),idiag)
      if (idiag == 0) call err('Non-existent diagnostic')

      ndiag(idiag) = intread(str,m)
      if (ndiag(idiag) < 0) call err('Bad diagnostic param')
      diagprev(idiag) = -1

c close out this diagnostic

      if (ndiag(idiag) == 0) then
        if (diagfileflag(idiag) > 0) call diag_close(idiag)
        diagfileflag(idiag) = 0
        if (node == 0) write (1,*) 'Diagnostic ',
     $       trim(work),' ',ndiag(idiag)
      else

c initialize this diagnostic

        call strread(str,m,diagfile(idiag))
        if (len_trim(diagfile(idiag)) == 0)
     $       call err('Bad diagnostic param')
        if (diagfileflag(idiag) > 0) call diag_close(idiag)
        if (trim(diagfile(idiag)) /= 'none') then
          call diag_open(idiag)
          diagfileflag(idiag) = 1
        endif
        diagnparams(idiag) = intread(str,m)
        if (diagnparams(idiag) < 0.or.diagnparams(idiag) > 5)
     $       call err('Bad diagnostic param')
        do i = 1,5
          diagparam(i,idiag) = 0.0
          if (diagnparams(idiag) >= i)
     $         diagparam(i,idiag) = realread(str,m)
        enddo
        if (node == 0) write (1,*) 'Diagnostic',
     $       trim(work),' ',ndiag(idiag),
     $       trim(diagfile(idiag)),
     $       diagparam(1,idiag),diagparam(2,idiag),diagparam(3,idiag),
     $       diagparam(4,idiag),diagparam(5,idiag)
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_temp_control(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,intread
      real*8 realread,tmp,ranmars
      logical match
      character*16 work

      call input_check(str,1)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad temp control parameter')
      if (match('none',work,mm)) then
        tempstyle = 0
      else if (match('rescale',work,mm)) then
        tempstyle = 1
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_every = intread(str,m)
        t_window = realread(str,m)
        t_fraction = realread(str,m)
        if (t_start < 0.0 .or. t_stop < 0.0 .or.
     $       t_every <= 0 .or. t_window < 0.0 .or.
     $       t_fraction < 0.0 .or. t_fraction > 1.0)
     $       call err('Bad temp control parameter')
      else if (match('replace',work,mm)) then
        tempstyle = 2
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_every = intread(str,m)
        iseed = realread(str,m)
        if (t_start < 0.0.or.t_stop < 0.0.or.t_every <= 0.or.
     $       iseed <= 0.or.iseed >= 1000000000)
     $       call err('Bad temp control parameter')
        tmp = ranmars(iseed+node)
      else if (match('langevin',work,mm)) then
        tempstyle = 3
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_freq = realread(str,m)
        iseed = intread(str,m)
        if (t_start < 0.0.or.t_stop < 0.0.or.
     $       t_freq < 0.0.or.iseed <= 0.or.iseed >= 1000000000)
     $       call err('Bad temp control parameter')
        tmp = ranmars(iseed+node)
      else if (match('nose/hoover',work,mm)) then
        tempstyle = 4
        t_start = realread(str,m)
        t_stop = realread(str,m)
        t_freq = realread(str,m)
        if (t_start < 0.0 .or. t_stop <= 0.0 .or. t_freq < 0.0)
     $       call err('Bad temp control parameter')
      else
        call err('Bad temp control parameter')
      endif
      if (node == 0) then
        if (tempstyle == 0) then
          write (1,*) 'Temp control none'
        else if (tempstyle == 1) then
          write (1,*) 'Temp control rescale',
     $         t_start,t_stop,t_every,t_window,t_fraction
        else if (tempstyle == 2) then
          write (1,*) 'Temp control replace',
     $         t_start,t_stop,t_every,iseed
        else if (tempstyle == 3) then
          write (1,*) 'Temp control langevin',
     $         t_start,t_stop,t_freq,iseed
        else if (tempstyle == 4) then
          write (1,*) 'Temp control nose/hoover',
     $         t_start,t_stop,t_freq
        endif
      endif
      if (tempstyle == 3) t_freq = t_freq*dtfactor
      if (tempstyle == 4) t_freq = t_freq*dtfactor

      return
      end

c ------------------------------------------------------------------------

      subroutine in_press_control(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm
      real*8 realread
      logical match
      character*16 work
      real*8 freq

      call input_check(str,1)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad press control parameter')

      if (match('none',work,mm)) then
        pressstyle = 0
      else if (match('nose/hoover',work,mm)) then
        pressstyle = 1

c pressure coupling

        call strread(str,m,work)
        if (len_trim(work) == 0)
     $       call err('Bad press control parameter')

        if (match('xyz',work,mm)) then
          presscouple = 0
        else if (match('xy',work,mm)) then
          presscouple = 1
        else if (match('yz',work,mm)) then
          presscouple = 2
        else if (match('xz',work,mm)) then
          presscouple = 3
        else if (match('aniso',work,mm)) then
          presscouple = 4
        else
          call err('Bad press control parameter')
        endif

c desired pressure values
c check for NULL entries

        if (presscouple == 0) then

          p_start(1) = realread(str,m)
          p_stop(1) = realread(str,m)
          p_start(2) = p_start(1)
          p_stop(2) = p_stop(1)
          p_start(3) = p_start(1)
          p_stop(3) = p_stop(1)
          p_freq(1) = 1.0
          p_freq(2) = 1.0
          p_freq(3) = 1.0

        else

          mm = m
          call strread(str,m,work)
          if (work(1:4) == 'NULL') then
            p_start(1) = 0.0
            p_stop(1) = 0.0
            p_freq(1) = 0.0
            call strread(str,m,work)
          else
            m = mm
            p_start(1) = realread(str,m)
            p_stop(1) = realread(str,m)
            p_freq(1) = 1.0
          endif

          mm = m
          call strread(str,m,work)
          if (work(1:4) == 'NULL') then
            p_start(2) = 0.0
            p_stop(2) = 0.0
            p_freq(2) = 0.0
            call strread(str,m,work)
          else
            m = mm
            p_start(2) = realread(str,m)
            p_stop(2) = realread(str,m)
            p_freq(2) = 1.0
          endif

          mm = m
          call strread(str,m,work)
          if (work(1:4) == 'NULL') then
            p_start(3) = 0.0
            p_stop(3) = 0.0
            p_freq(3) = 0.0
            call strread(str,m,work)
          else
            m = mm
            p_start(3) = realread(str,m)
            p_stop(3) = realread(str,m)
            p_freq(3) = 1.0
          endif

        endif

c set p_freq to input for non-NULL components, 0.0 for NULL components

        freq = realread(str,m)
        if (freq <= 0.0) call err('Bad press control parameter')
        if (p_freq(1) /= 0.0) p_freq(1) = freq
        if (p_freq(2) /= 0.0) p_freq(2) = freq
        if (p_freq(3) /= 0.0) p_freq(3) = freq

      else
        call err('Bad press control parameter')
      endif

c error check - coupled dimensions cannot be NULL

      if (presscouple == 1 .and.
     $     (p_freq(1) == 0.0 .or. p_freq(2) == 0.0))
     $     call err('Bad press control parameter')

      if (presscouple == 2 .and.
     $     (p_freq(2) == 0.0 .or. p_freq(3) == 0.0))
     $     call err('Bad press control parameter')

      if (presscouple == 3 .and.
     $     (p_freq(1) == 0.0 .or. p_freq(3) == 0.0))
     $     call err('Bad press control parameter')

c print out

      if (node == 0) then
        if (pressstyle == 0) then
          write (1,*) 'Press control none'
        else
          m = 1
          call strread(str,m,work)
          call strread(str,m,work)
          call strread(str,m,work)
          call strread(str,m,work)
          write (1,*) 'Press control nose/hoover ',
     $         trim(work),trim(str(m:))
        endif
      endif

c unit conversion

      if (pressstyle > 0) then
        p_start(1) = p_start(1)/pfactor
        p_stop(1) = p_stop(1)/pfactor
        p_freq(1) = p_freq(1)*p_freq(1)*dtfactor*dtfactor
        p_start(2) = p_start(2)/pfactor
        p_stop(2) = p_stop(2)/pfactor
        p_freq(2) = p_freq(2)*p_freq(2)*dtfactor*dtfactor
        p_start(3) = p_start(3)/pfactor
        p_stop(3) = p_stop(3)/pfactor
        p_freq(3) = p_freq(3)*p_freq(3)*dtfactor*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_volume_control(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm
      real*8 realread
      logical match
      character*16 work
      real*8 freq

      call input_check(str,1)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad volume control parameter')

      if (match('none',work,mm)) then
        volstyle = 0
      else if (match('linear',work,mm)) then
        volstyle = 1

        call strread(str,m,work)
        if (len_trim(work) == 0)
     $       call err('Bad volume control parameter')

        if (match('x',work,mm)) then
          voldimx = 1
          volstop_xlo = realread(str,m)
          volstop_xhi = realread(str,m)
          if (volstop_xlo >= volstop_xhi)
     $         call err('Bad volume control parameter')
        else if (match('y',work,mm)) then
          voldimy = 1
          volstop_ylo = realread(str,m)
          volstop_yhi = realread(str,m)
          if (volstop_ylo >= volstop_yhi)
     $         call err('Bad volume control parameter')
        else if (match('z',work,mm)) then
          voldimz = 1
          volstop_zlo = realread(str,m)
          volstop_zhi = realread(str,m)
          if (volstop_zlo >= volstop_zhi)
     $         call err('Bad volume control parameter')
        else
          call err('Bad volume control parameter')
        endif
      else
        call err('Bad press control parameter')
      endif

      if (node == 0) then
        if (volstyle == 0) then
          write (1,*) 'Volume control none'
        else
          m = 1
          call strread(str,m,work)
          call strread(str,m,work)
          call strread(str,m,work)
          call strread(str,m,work)
          write (1,*) 'Volume control linear ',
     $         trim(work),trim(str(m:))
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_nonbond_coeff(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer i,intread,j
      real*8 realread

      call input_check(str,1)
      if (nonstyle == 0) call err('Cannot define nonbond '//
     $     'coeffs when nonbond style = none')
      i = intread(str,m)
      j = intread(str,m)
      if (i <= 0.or.i > ntypes.or.j < i.or.j > ntypes)
     $     call err('Bad nonbond coeff parameter')
      nontypeflag(i,j) = nonstyle
      if (nonstyle == 1) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
      else if (nonstyle == 2) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
        noncoeff4(i,j) = realread(str,m)
      else if (nonstyle == 3) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
        noncoeff4(i,j) = realread(str,m)
      else if (nonstyle == 4) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
      else if (nonstyle == 5) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff3(i,j) = realread(str,m)
      else if (nonstyle == 6) then
        noncoeff1(i,j) = realread(str,m)
        noncoeff2(i,j) = realread(str,m)
        noncoeff14_1(i,j) = realread(str,m)
        noncoeff14_2(i,j) = realread(str,m)
      endif

      if (nonstyle == 1) then
        if (noncoeff1(i,j) < 0.0.or.noncoeff2(i,j) < 0.0.or.
     $       noncoeff3(i,j) < 0.0)
     $       call err('Bad nonbond coeff parameter')
      else if (nonstyle == 2) then
        if (noncoeff1(i,j) < 0.0.or.noncoeff2(i,j) < 0.0.or.
     $       noncoeff3(i,j) < 0.0.or.noncoeff4(i,j) < 0.0)
     $       call err('Bad nonbond coeff parameter')
        if (noncoeff3(i,j) == 0.0.and.noncoeff4(i,j) == 0.0) then
          continue
        else
          if (noncoeff3(i,j) == 0.0.or.
     $         noncoeff3(i,j) >= noncoeff4(i,j))
     $         call err('Bad nonbond coeff parameter')
        endif
      else if (nonstyle == 3) then
        if (noncoeff1(i,j) < 0.0.or.noncoeff2(i,j) < 0.0.or.
     $       noncoeff4(i,j) < 0.0)
     $       call err('Bad nonbond coeff parameter')
      else if (nonstyle == 4) then
        if (noncoeff1(i,j) < 0.0.or.noncoeff2(i,j) < 0.0.or.
     $       noncoeff3(i,j) < 0.0)
     $       call err('Bad nonbond coeff parameter')
      else if (nonstyle == 5) then
        if (noncoeff1(i,j) < 0.0.or.noncoeff2(i,j) < 0.0.or.
     $       noncoeff3(i,j) < 0.0)
     $       call err('Bad nonbond coeff parameter')
      else if (nonstyle == 6) then
        if (noncoeff1(i,j) < 0.0.or.noncoeff2(i,j) < 0.0.or.
     $      noncoeff14_1(i,j) < 0.0.or.noncoeff14_2(i,j) < 0.0)
     $       call err('Bad nonbond coeff parameter')
      endif

      if (node == 0) then
        if (nonstyle == 1) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
        else if (nonstyle == 2) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),
     $         noncoeff3(i,j),noncoeff4(i,j)
        else if (nonstyle == 3) then
          write (1,*) 'Nonbond coeff',
     $         noncoeff1(i,j),noncoeff2(i,j),
     $         noncoeff3(i,j),noncoeff4(i,j)
        else if (nonstyle == 4) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
        else if (nonstyle == 5) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),noncoeff3(i,j)
        else if (nonstyle == 6) then
          write (1,*) 'Nonbond coeff',i,j,
     $         noncoeff1(i,j),noncoeff2(i,j),
     $         noncoeff14_1(i,j),noncoeff14_2(i,j)
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_pppm_mesh(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      meshflag = 1
      nx_pppm_input = intread(str,m)
      ny_pppm_input = intread(str,m)
      nz_pppm_input = intread(str,m)
      if (nx_pppm_input <= 0.or.ny_pppm_input <= 0.or.
     $     nz_pppm_input <= 0)
     $     call err('Bad pppm mesh parameter')
      if (node == 0) write (1,*) 'Pppm mesh',
     $     nx_pppm_input,ny_pppm_input,nz_pppm_input

      return
      end

c ------------------------------------------------------------------------

      subroutine in_pppm_order(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      orderflag = intread(str,m)
      if (orderflag <= 0) call err('Bad pppm order parameter')
      if (node == 0)
     $     write (1,*) 'Pppm order',orderflag

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dielectric(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      real*8 realread

      call input_check(str,1)
      dielectric = realread(str,m)
      if (dielectric <= 0.0)
     $     call err('Bad dielectric parameter')
      if (node == 0)
     $     write (1,*) 'Dielectric',dielectric

      return
      end

c ------------------------------------------------------------------------

      subroutine in_bond_coeff(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer i,intread,j
      real*8 realread

      call input_check(str,1)
      if (bondstyle == 0) call err(
     $     'Cannot define bond coeffs when bond style = none')
      if (bondstyle == 5) call err(
     $     'Cannot define bond coeffs when bond style = class2')
      i = intread(str,m)
      if (i <= 0.or.i > nbondtypes)
     $     call err('Bad bond coeff parameter')
      bondtypeflag(i) = bondstyle
      if (bondstyle == 1) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
      else if (bondstyle == 2) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
        bondcoeff(3,i) = realread(str,m)
        bondcoeff(4,i) = realread(str,m)
      else if (bondstyle == 3) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
        bondcoeff(3,i) = realread(str,m)
        bondcoeff(4,i) = realread(str,m)
        bondcoeff(5,i) = realread(str,m)
      else if (bondstyle == 4) then
        bondcoeff(1,i) = realread(str,m)
        bondcoeff(2,i) = realread(str,m)
        bondcoeff(3,i) = realread(str,m)
      endif
      if (bondstyle == 1) then
        if (bondcoeff(1,i) < 0.0.or.
     $       bondcoeff(2,i) < 0.0)
     $       call err('Bad bond coeff parameter')
      else if (bondstyle == 2) then
        if (bondcoeff(1,i) < 0.0.or.
     $       bondcoeff(2,i) <= 0.0.or.
     $       bondcoeff(3,i) < 0.0.or.
     $       bondcoeff(4,i) <= 0.0)
     $       call err('Bad bond coeff parameter')
      else if (bondstyle == 3) then
        if (bondcoeff(1,i) < 0.0.or.
     $       bondcoeff(2,i) <= 0.0.or.
     $       bondcoeff(3,i) < 0.0.or.
     $       bondcoeff(4,i) <= 0.0)
     $       call err('Bad bond coeff parameter')
      else if (bondstyle == 4) then
        if (bondcoeff(1,i) < 0.0.or.
     $       bondcoeff(2,i) <= 0.0.or.
     $       bondcoeff(3,i) <= 0.0)
     $       call err('Bad bond coeff parameter')
      endif
      if (node == 0) then
        if (bondstyle == 1) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,2)
        else if (bondstyle == 2) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,4)
        else if (bondstyle == 3) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,5)
        else if (bondstyle == 4) then
          write (1,*) 'Bond coeff',(bondcoeff(j,i),j=1,3)
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_angle_coeff(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer i,intread,j
      real*8 realread

      call input_check(str,1)
      if (anglestyle == 0) call err(
     $     'Cannot define angle coeffs when angle style = none')
      if (anglestyle == 2) call err(
     $     'Cannot define angle coeffs when angle style = class2')

      i = intread(str,m)
      if (i <= 0 .or. i > nangletypes)
     $     call err('Bad angle coeff parameter')
      angletypeflag(i) = anglestyle

      if (anglestyle == 1) then
        anglecoeff(1,i) = realread(str,m)
        anglecoeff(2,i) = realread(str,m)
      else if (anglestyle == 3) then
        anglecoeff(1,i) = realread(str,m)
        anglecoeff(2,i) = realread(str,m)
        anglecoeff(3,i) = realread(str,m)
        anglecoeff(4,i) = realread(str,m)
      else if (anglestyle == 4) then
        anglecoeff(1,i) = realread(str,m)
      endif

      if (anglestyle == 1) then
        if (anglecoeff(1,i) < 0.0 .or.
     $       anglecoeff(2,i) <= 0.0)
     $       call err('Bad angle coeff parameter')
      else if (anglestyle == 3) then
        if (anglecoeff(1,i) < 0.0 .or. anglecoeff(2,i) <= 0.0 .or.
     $       anglecoeff(3,i) < 0.0 .or. anglecoeff(4,i) < 0.0)
     $       call err('Bad angle coeff parameter')
      else if (anglestyle == 4) then
        if (anglecoeff(1,i) < 0.0)
     $       call err('Bad angle coeff parameter')
      endif

      if (node == 0) then
        if (anglestyle == 1) then
          write (1,*) 'Angle coeff',(anglecoeff(j,i),j=1,2)
        else if (anglestyle == 3) then
          write (1,*) 'Angle coeff',(anglecoeff(j,i),j=1,4)
        else if (anglestyle == 4) then
          write (1,*) 'Angle coeff',(anglecoeff(j,i),j=1,1)
        endif
      endif

      if (anglestyle == 1)
     $     anglecoeff(2,i) = anglecoeff(2,i)/180.0 * 3.1415926
      if (anglestyle == 3)
     $     anglecoeff(2,i) = anglecoeff(2,i)/180.0 * 3.1415926

      return
      end

c ------------------------------------------------------------------------

      subroutine in_dihedral_coeff(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer i,intread,j
      real*8 realread

      call input_check(str,1)
      if (dihedstyle == 0) call err(
     $     'Cannot define dihedral coeffs when dihed style = none')
      if (dihedstyle == 2) call err(
     $     'Cannot define dihedral coeffs when dihed style = class2')

      i = intread(str,m)
      if (i <= 0 .or. i > ndihedtypes)
     $     call err('Bad dihedral coeff parameter')
      dihedtypeflag(i) = dihedstyle

      if (dihedstyle == 1) then
        dihedcoeff(1,i) = realread(str,m)
        dihedcoeff(2,i) = realread(str,m)
        dihedcoeff(3,i) = realread(str,m)
      else if (dihedstyle == 3) then
        dihedcoeff(1,i) = realread(str,m)
        dihedcoeff(2,i) = realread(str,m)
        dihedcoeff(3,i) = realread(str,m)
        dihedcoeff(4,i) = realread(str,m)
        dihedcoeff(5,i) = realread(str,m)
      else if (dihedstyle == 4) then
        dihedcoeff(1,i) = realread(str,m)
        dihedcoeff(2,i) = realread(str,m)
        dihedcoeff(3,i) = realread(str,m)
        dihedcoeff(4,i) = realread(str,m)
      endif

      if (dihedstyle == 1) then
        if (dihedcoeff(1,i) < 0.0 .or.
     $       (nint(dihedcoeff(2,i)) /= 1 .and.
     $       nint(dihedcoeff(2,i)) /= -1) .or.
     $       nint(dihedcoeff(3,i)) < 1 .or. nint(dihedcoeff(3,i)) > 6)
     $       call err('Bad dihedral coeff parameter')
      else if (dihedstyle == 4) then
        if ((nint(dihedcoeff(3,i)) /= 0 .and.
     $       nint(dihedcoeff(3,i)) /= 180) .or.
     $       nint(dihedcoeff(2,i)) < 0 .or. nint(dihedcoeff(2,i)) > 6
     $       .or. dihedcoeff(4,i) < 0.0 .or. dihedcoeff(4,i) > 1.0) 
     $       call err('Bad dihedral coeff parameter')
      endif

      if (node == 0) then
        if (dihedstyle == 1) then
          write (1,*) 'Dihedral coeff',(dihedcoeff(j,i),j=1,3)
        else if (dihedstyle == 3) then
          write (1,*) 'Dihedral coeff',(dihedcoeff(j,i),j=1,5)
        else if (dihedstyle == 4) then
          write (1,*) 'Dihedral coeff',(dihedcoeff(j,i),j=1,4)
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_improper_coeff(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer i,intread,j
      real*8 realread

      call input_check(str,1)
      if (improstyle == 0) call err(
     $     'Cannot define improper coeffs when impro style = none')
      if (improstyle == 3) call err(
     $     'Cannot define improper coeffs when impro style = class2')
      i = intread(str,m)
      if (i <= 0 .or. i > nimprotypes)
     $     call err('Bad improper coeff parameter')
      improtypeflag(i) = improstyle
      if (improstyle == 1) then
        improcoeff(1,i) = realread(str,m)
        improcoeff(2,i) = realread(str,m)
      else if (improstyle == 2) then
        improcoeff(1,i) = realread(str,m)
        improcoeff(2,i) = realread(str,m)
        improcoeff(3,i) = realread(str,m)
      endif
      if (improstyle == 1) then
        if (improcoeff(1,i) < 0.0 .or.
     $       improcoeff(2,i) < 0.0)
     $       call err('Bad improper coeff parameter')
      else if (improstyle == 2) then
        if (improcoeff(1,i) < 0.0 .or.
     $       (nint(improcoeff(2,i)) /= 1 .and.
     $       nint(improcoeff(2,i)) /= -1) .or.
     $       nint(improcoeff(3,i)) < 0 .or. nint(improcoeff(3,i)) > 6)
     $       call err('Bad improper coeff parameter')
      endif
      if (node == 0) then
        if (improstyle == 1) then
          write (1,*) 'Improper coeff',(improcoeff(j,i),j=1,2)
        else if (improstyle == 2) then
          write (1,*) 'Improper coeff',(improcoeff(j,i),j=1,3)
        endif
      endif
      if (improstyle == 1)
     $     improcoeff(2,i) = improcoeff(2,i)/180.0 * 3.1415926

      return
      end

c ------------------------------------------------------------------------

      subroutine in_min_style(str,m)
      use global
      implicit none

c argument variables

      character*(*)  str
      integer m

c local variables

      integer mm
      logical match
      character*16 work

      call input_check(str,1)
      call strread(str,m,work)
      if (match('hftn',work,mm)) then
        optstyle = 1
      else
        call err('Bad min style parameter')
      endif
      if (node == 0) then
        if (optstyle == 1) write (1,*) 'Min style hftn'
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_min_flag(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      optflag = intread(str,m)
      if (optflag < 0) call err('Bad min flag parameter')
      if (node == 0) write (1,*) 'Min flag',optflag

      return
      end

c ------------------------------------------------------------------------

      subroutine in_rotation_zero(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      rotationflag = intread(str,m)
      if (rotationflag < 0 .or. rotationflag > 1)
     $     call err('Bad rotation zero parameter')
      if (node == 0) write (1,*) 'Rotation zero',rotationflag

      return
      end

c ------------------------------------------------------------------------

      subroutine in_create_group(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,intread,i
      real*8 toggle,realread

      logical match
      character*16 work

      call input_check(str,1)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad create group parameter')
      if (match('types',work,mm)) then
        creategroup = 1
        createlo = intread(str,m)
        createhi = intread(str,m)
        if (createlo <= 0 .or. createhi < createlo)
     $       call err('Bad create style parameter')
      else if (match('molecules',work,mm)) then
        creategroup = 2
        createlo = intread(str,m)
        createhi = intread(str,m)
        if (createlo < 0 .or. createhi < createlo)
     $       call err('Bad create style parameter')
      else if (match('region',work,mm)) then
        creategroup = 3
        toggle = -1.0
        do i = 1,6
          mm = m
          call strread(str,m,work)
          if (work(1:3) == 'INF') then
            createregion(i) = toggle*1.0E20
          else
            m = mm
            createregion(i) = realread(str,m)
          endif
          toggle = -toggle
        enddo
        if (createregion(1) > createregion(2).or.
     $       createregion(3) > createregion(4).or.
     $       createregion(5) > createregion(6))
     $       call err('Bad create group parameter')
      else if (match('remainder',work,mm)) then
        creategroup = 4
      else
        call err('Bad create group parameter')
      endif
      if (node == 0) then
        if (creategroup == 1) then
          write (1,*) 'Create group types',createlo,createhi
        else if (creategroup == 1) then
          write (1,*) 'Create group molecules',createlo,createhi
        else if (creategroup == 3) then
          write (1,*) 'Create group region',
     $         createregion(1),createregion(2),createregion(3),
     $         createregion(4),createregion(5),createregion(6)
        else if (creategroup == 4) then
          write (1,*) 'Create group remainder'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_create_temp(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,intread
      real*8 realread

      logical match
      character*16 work

      call input_check(str,1)
      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad create temp parameter')
      if (match('uniform',work,mm)) then
        createstyle = 1
        t_create = realread(str,m)
        iseed = intread(str,m)
        if (t_create < 0.0.or.iseed <= 0.or.iseed >= 1000000000)
     $       call err('Bad create temp parameter')
      else if (match('gaussian',work,mm)) then
        createstyle = 2
        t_create = realread(str,m)
        iseed = intread(str,m)
        if (t_create < 0.0.or.iseed <= 0.or.iseed >= 1000000000)
     $       call err('Bad create temp parameter')
      else if (match('velocity',work,mm)) then
        createstyle = 3
        createvec(1) = realread(str,m)
        createvec(2) = realread(str,m)
        createvec(3) = realread(str,m)
      else
        call err('Bad create temp parameter')
      endif
      if (node == 0) then
        if (createstyle == 1) then
          write (1,*) 'Create temp uniform ',t_create,iseed
        else if (createstyle == 2) then
          write (1,*) 'Create temp gaussian ',t_create,iseed
        else if (createstyle == 3) then
          write (1,*) 'Create temp velocity ',
     $         createvec(1),createvec(2),createvec(3)
        endif
      endif
      if (createstyle == 3) then
        createvec(1) = createvec(1)*dtfactor
        createvec(2) = createvec(2)*dtfactor
        createvec(3) = createvec(3)*dtfactor
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_fix_style(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer mm,i,iwhich,intread
      real*8 realread,tmp,ranmars

      logical match
      character*16 work

      call input_check(str,1)
      mm = m

      call strread(str,mm,work)
      if (len_trim(work) == 0)
     $     call err('Bad fix style parameter')

      if (match('none',work,mm)) then
        nfixes = 0
        nfixes_respa = 0
        nshake = 0
        do i = 1,maxfix
          fixstyle(i) = 0
        enddo
        do i = 1,nlocal
          fix(i) = 0
        enddo

      else
        iwhich = intread(str,m)
        if (iwhich <= 0) call err('Bad fix style parameter')
        if (iwhich >= maxfix)
     $       call err('Too many fixes - boost maxfix')
        nfixes = max(nfixes,iwhich)
        call strread(str,m,work)

        if (match('setforce',work,mm)) then
          fixstyle(iwhich) = 1
          do i = 1,3
            mm = m
            call strread(str,m,work)
            if (work(1:4) == 'NULL') then
              fixflag(i,iwhich) = 1
              fixcoeff(i,iwhich) = 0.0
            else
              fixflag(i,iwhich) = 0
              m = mm
              fixcoeff(i,iwhich) = realread(str,m)
            endif
          enddo

        else if (match('addforce',work,mm)) then
          fixstyle(iwhich) = 2
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)

        else if (match('aveforce',work,mm)) then
          fixstyle(iwhich) = 3
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)

        else if (match('rescale',work,mm)) then
          fixstyle(iwhich) = 4
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = intread(str,m)
          fixcoeff(4,iwhich) = realread(str,m)
          fixcoeff(5,iwhich) = realread(str,m)
          if (fixcoeff(1,iwhich) < 0.0 .or.
     $         fixcoeff(2,iwhich) < 0.0 .or.
     $         nint(fixcoeff(3,iwhich)) <= 0 .or.
     $         fixcoeff(4,iwhich) < 0.0 .or.
     $         fixcoeff(5,iwhich) < 0.0 .or. fixcoeff(5,iwhich) > 1.0)
     $         call err('Bad fix style parameter')

        else if (match('hoover/drag',work,mm)) then
          fixstyle(iwhich) = 5
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)
c set initial drag to 0.0
          fixcoeff(4,iwhich) = 0.0
          if (fixcoeff(1,iwhich) < 0.0 .or.
     $         fixcoeff(2,iwhich) < 0.0 .or.
     $         fixcoeff(3,iwhich) < 0.0 .or. fixcoeff(3,iwhich) > 1.0)
     $         call err('Bad fix style parameter')

        else if (match('langevin',work,mm)) then
          fixstyle(iwhich) = 6
          fixcoeff(1,iwhich) = realread(str,m)
          fixcoeff(2,iwhich) = realread(str,m)
          fixcoeff(3,iwhich) = realread(str,m)
          iseed = intread(str,m)
          fixflag(1,iwhich) = intread(str,m)
          fixflag(2,iwhich) = intread(str,m)
          fixflag(3,iwhich) = intread(str,m)
          if (fixcoeff(1,iwhich) < 0.0.or.
     $         fixcoeff(2,iwhich) < 0.0.or.
     $         fixcoeff(3,iwhich) < 0.0.or.
     $         iseed <= 0.or.iseed >= 1000000000)
     $         call err('Bad fix style parameter')
          if (fixflag(1,iwhich) < 0.or.fixflag(1,iwhich) > 1.or.
     $         fixflag(2,iwhich) < 0.or.fixflag(2,iwhich) > 1.or.
     $         fixflag(3,iwhich) < 0.or.fixflag(3,iwhich) > 1)
     $         call err('Bad fix style parameter')
          tmp = ranmars(iseed+node)

        else if (match('springforce',work,mm)) then
          fixstyle(iwhich) = 7
          do i = 1,3
            mm = m
            call strread(str,m,work)
            if (work(1:4) == 'NULL') then
              fixflag(i,iwhich) = 1
              fixcoeff(i,iwhich) = 0.0
            else
              fixflag(i,iwhich) = 0
              m = mm
              fixcoeff(i,iwhich) = realread(str,m)
            endif
          enddo
          fixcoeff(4,iwhich) = realread(str,m)
          if (fixcoeff(4,iwhich) < 0.0)
     $         call err('Bad fix style parameter')

        else if (match('dragforce',work,mm)) then
          fixstyle(iwhich) = 8
          do i = 1,3
            mm = m
            call strread(str,m,work)
            if (work(1:4) == 'NULL') then
              fixflag(i,iwhich) = 1
              fixcoeff(i,iwhich) = 0.0
            else
              fixflag(i,iwhich) = 0
              m = mm
              fixcoeff(i,iwhich) = realread(str,m)
            endif
          enddo
          fixcoeff(4,iwhich) = realread(str,m)
          fixcoeff(5,iwhich) = realread(str,m)
          if (fixcoeff(4,iwhich) < 0.0.or.fixcoeff(5,iwhich) < 0.0)
     $         call err('Bad fix style parameter')

        else if (match('shake',work,mm)) then
          fixstyle(iwhich) = 9
          shakeiter = intread(str,m)
          shaketol = realread(str,m)
          nshakestats = intread(str,m)
          if (shakeiter <= 0 .or. shaketol < 0.0 .or. nshakestats < 0)
     $         call err('Bad fix style parameter')
          nshake = 1
          do i = 1,nbondtypes
            shakeablebond(i) = 0
          enddo
          shakeableangle = 0

        else
          call err('Bad fix style parameter')
        endif
      endif

      if (node == 0) then
        if (nfixes == 0) then
          write (1,*) 'Fix style none'
        else if (fixstyle(iwhich) == 1) then
          write (1,*) 'Fix style',iwhich,' setforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich) == 2) then
          write (1,*) 'Fix style',iwhich,' addforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich) == 3) then
          write (1,*) 'Fix style',iwhich,' aveforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich) == 4) then
          write (1,*) 'Fix style',iwhich,' rescale',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         nint(fixcoeff(3,iwhich)),fixcoeff(4,iwhich),
     $         fixcoeff(5,iwhich)
        else if (fixstyle(iwhich) == 5) then
          write (1,*) 'Fix style',iwhich,' hoover/drag',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich)
        else if (fixstyle(iwhich) == 6) then
          write (1,*) 'Fix style',iwhich,' langevin',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich),iseed,
     $         fixflag(1,iwhich),fixflag(2,iwhich),fixflag(3,iwhich)
        else if (fixstyle(iwhich) == 7) then
          write (1,*) 'Fix style',iwhich,' springforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich),fixcoeff(4,iwhich)
        else if (fixstyle(iwhich) == 8) then
          write (1,*) 'Fix style',iwhich,' dragforce',
     $         fixcoeff(1,iwhich),fixcoeff(2,iwhich),
     $         fixcoeff(3,iwhich),fixcoeff(4,iwhich),
     $         fixcoeff(5,iwhich)
        else if (fixstyle(iwhich) == 9) then
          write (1,*) 'Fix style',iwhich,' shake',
     $         shakeiter,shaketol,nshakestats
        endif
      endif

      if (nfixes > 0) then
        if (fixstyle(iwhich) == 1 .or. fixstyle(iwhich) == 2 .or.
     $       fixstyle(iwhich) == 3) then
          fixcoeff(1,iwhich) = fixcoeff(1,iwhich)*dtfactor*dtfactor
          fixcoeff(2,iwhich) = fixcoeff(2,iwhich)*dtfactor*dtfactor
          fixcoeff(3,iwhich) = fixcoeff(3,iwhich)*dtfactor*dtfactor
        endif
        if (fixstyle(iwhich) == 5)
     $       fixcoeff(3,iwhich) = fixcoeff(3,iwhich)*dtfactor
        if (fixstyle(iwhich) == 6)
     $       fixcoeff(3,iwhich) = fixcoeff(3,iwhich)*dtfactor
        if (fixstyle(iwhich) == 7)
     $       fixcoeff(4,iwhich) = fixcoeff(4,iwhich)*dtfactor*dtfactor
        if (fixstyle(iwhich) == 8)
     $       fixcoeff(4,iwhich) = fixcoeff(4,iwhich)*dtfactor*dtfactor
      endif
      
      return
      end

c ------------------------------------------------------------------------

      subroutine in_assign_fix(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread,mm,i
      real*8 toggle,realread

      logical match
      character*16 work

      call input_check(str,1)
      fixwhich = intread(str,m)
      if (fixwhich <= 0 .or. fixwhich >= maxfix)
     $     call err('Bad assign fix parameter')

      call strread(str,m,work)
      if (len_trim(work) == 0)
     $     call err('Bad assign fix parameter')

      if (match('atom',work,mm)) then
        fixgroup = 1
        fixatom = intread(str,m)
        if (fixatom <= 0 .or. fixatom > natoms) 
     $       call err('Bad assign fix parameter')
      else if (match('molecule',work,mm)) then
        fixgroup = 2
        fixatom = intread(str,m)
        if (fixatom <= 0) 
     $       call err('Bad assign fix parameter')
      else if (match('type',work,mm)) then
        fixgroup = 3
        fixtype = intread(str,m)
        if (fixtype <= 0 .or. fixtype > ntypes)
     $       call err('Bad assign fix parameter')
      else if (match('region',work,mm)) then
        fixgroup = 4
        toggle = -1.0
        do i = 1,6
          mm = m
          call strread(str,m,work)
          if (work(1:3) == 'INF') then
            fixregion(i) = toggle*1.0E20
          else
            m = mm
            fixregion(i) = realread(str,m)
          endif
          toggle = -toggle
        enddo
        if (fixregion(1) > fixregion(2).or.
     $       fixregion(3) > fixregion(4).or.
     $       fixregion(5) > fixregion(6))
     $       call err('Bad assign fix parameter')
      else if (match('bondtype',work,mm)) then
        fixgroup = 5
        fixbond = intread(str,m)
        if (fixbond <= 0 .or. fixbond > nbondtypes)
     $       call err('Bad assign fix parameter')
        if (fixstyle(fixwhich) /= 9) call err(
     $       'Assign fix bondtype must be used with SHAKE fix')
      else if (match('angletype',work,mm)) then
        fixgroup = 6
        fixangle = intread(str,m)
        fixbond = intread(str,m)
        if (fixangle <= 0 .or. fixangle > nangletypes)
     $       call err('Bad assign fix parameter')
        if (fixbond <= 0 .or. fixbond > nbondtypes)
     $       call err('Bad assign fix parameter')
        if (fixstyle(fixwhich) /= 9) call err(
     $       'Assign fix angletype must be used with SHAKE fix')
      else if (match('remainder',work,mm)) then
        fixgroup = 7
      else
        call err('Bad assign fix parameter')
      endif
      if (node == 0) then
        if (fixgroup == 1) then
          write (1,*) 'Assign fix',fixwhich,' atom',fixatom
        else if (fixgroup == 2) then
          write (1,*) 'Assign fix',fixwhich,' molecule',fixatom
        else if (fixgroup == 3) then
          write (1,*) 'Assign fix',fixwhich,' type',fixtype
        else if (fixgroup == 4) then
          write (1,*) 'Assign fix',fixwhich,' region',
     $         (fixregion(i),i=1,6)
        else if (fixgroup == 5) then
          write (1,*) 'Assign fix',fixwhich,' bondtype',fixbond
        else if (fixgroup == 6) then
          write (1,*) 'Assign fix',fixwhich,' angletype',
     $         fixangle,fixbond
        else if (fixgroup == 7) then
          write (1,*) 'Assign fix',fixwhich,' remainder'
        endif
      endif

      return
      end

c ------------------------------------------------------------------------

      subroutine in_run(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread

      call input_check(str,1)
      nsteps = intread(str,m)
      if (nsteps < 0) call err('Bad run parameter')
      if (node == 0) write (1,*) 'Run',nsteps

      return
      end

c ------------------------------------------------------------------------

      subroutine in_minimize(str,m)
      use global
      implicit none

c argument variables

      character*(*) str
      integer m

c local variables

      integer intread
      real*8 realread

      call input_check(str,1)
      opt_stop_tol = realread(str,m)
      opt_max_iters = intread(str,m)
      opt_max_fns = intread(str,m)
      nsteps = 1
      if (opt_stop_tol <= 0.0.or.
     $     opt_max_iters < 0.or.opt_max_fns < 0)
     $     call err('Bad minimize parameter')
      if (node == 0)
     $     write (1,*) 'Minimize',
     $     opt_stop_tol,opt_max_iters,opt_max_fns

      return
      end

c ------------------------------------------------------------------------
c check for positioning of command
c readflag tells whether "read data" or "read restart" has been seen already
c iflag = 0 = command must be before "read data"
c iflag = 1 = command must be after "read data"

      subroutine input_check(str,iflag)
      use global
      implicit none

c argument variables

      character*(*) str
      integer iflag

      if (iflag == 1 .and. readflag == 0) then
        if (node == 0) write (6,*) 'Command: ',trim(str)
        call err(
     $     'Illegal command before read_data or read_restart')
      endif
      if (iflag == 0 .and. readflag == 1) then
        if (node == 0) write (6,*) 'Command: ',trim(str)
        call err(
     $     'Illegal command after read_data or read_restart')
      endif

      return
      end
