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
c create velocities for atoms according to createstyle and creategroup
c
c   createstyle = 1 -> init atoms by uniform RNG to t_create
c   createstyle = 2 -> init atoms by gaussian RNG to t_create
c   createstyle = 3 -> init atoms to velocity of createvec(3)
c
c   creategroup = 0 -> apply to all atoms
c   creategroup = 1 -> only apply to atoms with createlo <= type <= createhi
c   creategroup = 2 -> only apply to molecules with 
c                      createlo <= molecule # <= createhi
c   creategroup = 3 -> apply to atoms within a spatial region
c   creategroup = 4 -> apply to all unassigned atoms

      subroutine temp_create
      use global
      use mpi
      implicit none

c local variables

      integer ncount,i,iflag,itmp,ierror,j
      real*8 countmass,tmp,vx,ranpark,vy,vz,gausspark
      real*8 factor

      if (node == 0) write (6,*) 'Creating temperature ...'

      if (readflag == 0)
     $     call err('Data file must be read before creating temp')

      if (creategroup == 1) then
        if (createlo.gt.ntypes.or.createhi.gt.ntypes)
     $       call err('Creating temp for too large a type')
      endif

c mark atoms whose velocities will be initialized and count them

      ncount = 0
      countmass = 0.0

      do i = 1,nlocal
        iflag = 0
        if (creategroup == 0) then
          iflag = 1
        else if (creategroup == 1) then
          if (createlo <= type(i) .and.
     $         type(i) <= createhi) iflag = 1
        else if (creategroup == 2) then
          if (createlo <= molecule(i) .and.
     $         molecule(i) <= createhi) iflag = 1
        else if (creategroup == 3) then
          if (x(1,i) >= createregion(1).and.
     $         x(1,i) <= createregion(2).and.
     $         x(2,i) >= createregion(3).and.
     $         x(2,i) <= createregion(4).and.
     $         x(3,i) >= createregion(5).and.
     $         x(3,i) <= createregion(6)) iflag = 1
        else if (creategroup == 4) then
          if (velflag(i) == 0) iflag = 1
        endif
        if (iflag == 1) then
          velflag(i) = -1
          ncount = ncount + 1
          countmass = countmass + mass(type(i))
        endif
      enddo

      itmp = ncount
      call mpi_allreduce(itmp,ncount,1,mpi_integer,
     $     mpi_sum,mpi_comm_world,ierror)
      if (ncount == 0) return

      tmp = countmass
      call mpi_allreduce(tmp,countmass,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      if (createstyle.le.2) then

c for uniform and Gaussian initializers, loop over all atoms to insure
c   RNG usage assigns velocities the same, independent of # of procs
c weight initial velocities by 1/sqrt(mass)
c   so will contribute equally to temp
c old style was independent of mass, factor = 1.0

        do j = 1,natoms
          if (createstyle == 1) then
            vx = ranpark(iseed)
            vy = ranpark(iseed)
            vz = ranpark(iseed)
          else
            vx = gausspark(iseed)
            vy = gausspark(iseed)
            vz = gausspark(iseed)
          endif
          i = localptr(j)
          if (i /= 0) then
            if (velflag(i) == -1) then
              factor = 1.0/sqrt(mass(type(i)))
              v(1,i) = vx * factor
              v(2,i) = vy * factor
              v(3,i) = vz * factor
              if (idimension == 2) v(3,i) = 0.0
            endif
          endif
        enddo

c zero out momentum and angular rotation of initialized atoms only
c don't bother if just one atom

        if (ncount > 1) call momentum_zero(ncount,countmass)
        if (rotationflag == 1 .and. ncount > 1)
     $       call rotation_zero(ncount,countmass)

c compute temp of initialized atoms only

        t_current = 0.0

        do i = 1,nlocal
          if (velflag(i) == -1) t_current = t_current +
     $         (v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i)) *
     $         mass(type(i))
        enddo
        call mpi_allreduce(t_current,tmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        t_current = tmp / (boltz*idimension*ncount)

c rescale only initialized atoms to temp t_create

        if (t_current /= 0.0) then
          factor = sqrt(t_create/t_current)
          do i = 1,nlocal
            if (velflag(i) == -1) then
              v(1,i) = v(1,i)*factor
              v(2,i) = v(2,i)*factor
              v(3,i) = v(3,i)*factor
            endif
          enddo
        endif

      else if (createstyle == 3) then

        do i = 1,nlocal
          if (velflag(i) == -1) then
            v(1,i) = createvec(1)
            v(2,i) = createvec(2)
            v(3,i) = createvec(3)
          endif
        enddo

      endif

c mark atoms as initialized

      do i = 1,nlocal
        if (velflag(i) == -1) velflag(i) = 1
      enddo

c print message on how many atoms initialized

      if (node == 0) then
        write (6,*) '  Velocites set for this many atoms:',ncount
        write (1,*) '  Velocites set for this many atoms:',ncount
      endif

c reset creategroup to all atoms

      creategroup = 0

      return
      end


c -------------------------------------------------------------------------
c rescale temperature to target temp

      subroutine temp_rescale(target)
      use global
      implicit none

c argument variables

      real*8 target

c local variables

      integer i
      real*8 factor

      factor = sqrt(target/t_current)
      do i = 1,nlocal
        v(1,i) = v(1,i)*factor
        v(2,i) = v(2,i)*factor
        v(3,i) = v(3,i)*factor
      enddo
      t_current = target

      return
      end


c -------------------------------------------------------------------------
c Gaussian replacement of all velocities at target temp

      subroutine temp_replace(target)
      use global
      implicit none

c argument variables

      real*8 target

c local variables

      integer i
      real*8 gaussmars

      do i = 1,nlocal
        v(1,i) = gaussmars()
        v(2,i) = gaussmars()
        v(3,i) = gaussmars()
        if (idimension == 2) v(3,i) = 0.0
      enddo
      call momentum_zero(0,masssum)
      if (rotationflag == 1) call rotation_zero(0,masssum)
      call temperature
      if (t_current.ne.0.0) call temp_rescale(target)
      
      return
      end


c -------------------------------------------------------------------------
c subtract center-of-mass momentum from velocities
c  ncount = 0 -> include all atoms
c  ncount > 1 -> include only atoms tagged with velflag(i) = -1

      subroutine momentum_zero(ncount,countmass)
      use global
      use mpi
      implicit none

c argument variables

      integer ncount
      real*8 countmass

c local variables

      integer i,ierror
      real*8 rmass
      real*8 ptot(3),ptotall(3)

      ptot(1) = 0.0
      ptot(2) = 0.0
      ptot(3) = 0.0

      if (ncount == 0) then
        do i = 1,nlocal
          rmass = mass(type(i))
          ptot(1) = ptot(1) + v(1,i)*rmass
          ptot(2) = ptot(2) + v(2,i)*rmass
          ptot(3) = ptot(3) + v(3,i)*rmass
        enddo
      else
        do i = 1,nlocal
          if (velflag(i) == -1) then
            rmass = mass(type(i))
            ptot(1) = ptot(1) + v(1,i)*rmass
            ptot(2) = ptot(2) + v(2,i)*rmass
            ptot(3) = ptot(3) + v(3,i)*rmass
          endif
        enddo
      endif

      call mpi_allreduce(ptot,ptotall,3,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      ptot(1) = ptotall(1) / countmass
      ptot(2) = ptotall(2) / countmass
      ptot(3) = ptotall(3) / countmass

      if (ncount == 0) then
        do i = 1,nlocal
          v(1,i) = v(1,i) - ptot(1)
          v(2,i) = v(2,i) - ptot(2)
          v(3,i) = v(3,i) - ptot(3)
        enddo
      else
        do i = 1,nlocal
          if (velflag(i) == -1) then
            v(1,i) = v(1,i) - ptot(1)
            v(2,i) = v(2,i) - ptot(2)
            v(3,i) = v(3,i) - ptot(3)
          endif
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c subtract aggregate angular rotation from velocities
c  ncount = 0 -> include all atoms
c  ncount > 1 -> include only atoms tagged with velflag(i) = -1

      subroutine rotation_zero(ncount,countmass)
      use global
      use mpi
      implicit none

c argument variables

      integer ncount
      real*8 countmass

c local variables

      real*8 rm(3),rmall(3),r2,rmass,small
      integer i,ierror

      small = 0.001

      rm(1) = 0.0
      rm(2) = 0.0
      rm(3) = 0.0
      
      if (ncount == 0) then
        do i = 1,nlocal
          rmass = mass(type(i))
          rm(1) = rm(1) + rmass * (x(2,i)*v(3,i) - x(3,i)*v(2,i))
          rm(2) = rm(2) + rmass * (x(3,i)*v(1,i) - x(1,i)*v(3,i))
          rm(3) = rm(3) + rmass * (x(1,i)*v(2,i) - x(2,i)*v(1,i))
        enddo
      else
        do i = 1,nlocal
          if (velflag(i) == -1) then
            rmass = mass(type(i))
            rm(1) = rm(1) + rmass * (x(2,i)*v(3,i) - x(3,i)*v(2,i))
            rm(2) = rm(2) + rmass * (x(3,i)*v(1,i) - x(1,i)*v(3,i))
            rm(3) = rm(3) + rmass * (x(1,i)*v(2,i) - x(2,i)*v(1,i))
          endif
        enddo
      endif

      call mpi_allreduce(rm,rmall,3,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      rm(1) = rmall(1) / countmass
      rm(2) = rmall(2) / countmass
      rm(3) = rmall(3) / countmass
      
      if (ncount == 0) then
        do i = 1,nlocal
          r2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i) + small
          v(1,i) = v(1,i) - (rm(2)*x(3,i) - rm(3)*x(2,i)) / r2
          v(2,i) = v(2,i) - (rm(3)*x(1,i) - rm(1)*x(3,i)) / r2
          v(3,i) = v(3,i) - (rm(1)*x(2,i) - rm(2)*x(1,i)) / r2
        enddo
      else
        do i = 1,nlocal
          if (velflag(i) == -1) then
            r2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i) + small
            v(1,i) = v(1,i) - (rm(2)*x(3,i) - rm(3)*x(2,i)) / r2
            v(2,i) = v(2,i) - (rm(3)*x(1,i) - rm(1)*x(3,i)) / r2
            v(3,i) = v(3,i) - (rm(1)*x(2,i) - rm(2)*x(1,i)) / r2
          endif
        enddo
      endif

      return
      end
