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
c zero out force and virial arrays

      subroutine zero_force_virial
      use global
      implicit none

c local variables

      integer i

      do i = 1,6
        virial(i) = 0.0
      enddo

      if (newton == 0) then
        do i = 1,nlocal
          f(1,i) = 0.0
          f(2,i) = 0.0
          f(3,i) = 0.0
        enddo
      else
        do i = 1,nlocal+nghost
          f(1,i) = 0.0
          f(2,i) = 0.0
          f(3,i) = 0.0
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c copy force array to target array

      subroutine copy_force(n,f_target)
      use global
      implicit none

c argument variables

      integer n
      real*8 f_target
      dimension f_target(3,*)

c local variables

      integer i

      do i = 1,n
        f_target(1,i) = f(1,i)
        f_target(2,i) = f(2,i)
        f_target(3,i) = f(3,i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c copy virial array to target array

      subroutine copy_virial(n,vir_target)
      use global
      implicit none

c argument variables

      integer n
      real*8 vir_target
      dimension vir_target(6)

c local variables

      integer i

      do i = 1,6
        vir_target(i) = virial(i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c minimum image convention inside a periodic box
c adjust dx,dy,dz to magnitude and sign of shortest distance in box

      subroutine minimg(dx,dy,dz)
      use global
      implicit none

c argument variables

      real*8 dx,dy,dz

      if (perflagx == 0) then
        if (abs(dx) > xprd_half) then
          if (dx < 0.0) then
            dx = dx + xprd
          else
            dx = dx - xprd
          endif
        endif
      endif

      if (perflagy == 0) then
        if (abs(dy) > yprd_half) then
          if (dy < 0.0) then
            dy = dy + yprd
          else
            dy = dy - yprd
          endif
        endif
      endif

      if (perflagz == 0) then
        if (abs(dz) > zprd_half) then
          if (dz < 0.0) then
            dz = dz + zprd
          else
            dz = dz - zprd
          endif
        endif
      endif
      
      return
      end


c -------------------------------------------------------------------------
c remap the point (xx,yy,zz) into the periodic box 
c  no matter how far away it is
c adjust true flag accordingly
c only do it if periodicity is on

      subroutine remap(xx,yy,zz,itrue)
      use global
      implicit none

c argument variables

      real*8 xx,yy,zz
      integer itrue

      if (perflagx == 0) then
        do while (xx >= box(2,1))
          xx = xx - xprd
          itrue = itrue + 1
        enddo
        do while (xx < box(1,1))
          xx = xx + xprd
          itrue = itrue - 1
        enddo
      endif

      if (perflagy == 0) then
        do while (yy >= box(2,2))
          yy = yy - yprd
          itrue = itrue + 1000
        enddo
        do while (yy < box(1,2))
          yy = yy + yprd
          itrue = itrue - 1000
        enddo
      endif

      if (perflagz == 0) then
        do while (zz >= box(2,3))
          zz = zz - zprd
          itrue = itrue + 1000000
        enddo
        do while (zz < box(1,3))
          zz = zz + zprd
          itrue = itrue - 1000000
        enddo
      endif

      return
      end
