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
c compute SHAKE constraints on appropriate atom groups
c iflag determines whether being called from start (apply for 1/2 step)
c   or from integrate (apply for full step)

      subroutine shake(iflag,vflag)
      use global
      use mpi
      implicit none

c argument variables

      integer iflag,vflag

c local variables

      integer i,ngroup
      real*8 dtsq_factor

c compute new coords without SHAKE constraints
c don't need to communicate if 1 proc, since ghost atoms won't be accessed

      call shake_update(iflag)
      if (nprocs > 1) call shake_comm

c compute SHAKE constraints for each SHAKE group in local neighbor list
c ngroup = 3,4,5 for bond-only SHAKE, -3 for bond/angle SHAKE

      if (iflag == 0) dtsq_factor = 0.5*dt*dt
      if (iflag == 1) dtsq_factor = dt*dt

      do i = 1,nshakelocal

        ngroup = shakesize(i)

        if (ngroup == 2) then
          call shake2(shakeatom(1,i),shakeatom(2,i),
     $         shakebondlen(1,i),dtsq_factor,vflag)
        else if (ngroup == 3) then
          call shake3(shakeatom(1,i),shakeatom(2,i),shakeatom(3,i),
     $         shakebondlen(1,i),shakebondlen(2,i),
     $         shaketol,shakeiter,dtsq_factor,vflag)
        else if (ngroup == 4) then
          call shake4(shakeatom(1,i),shakeatom(2,i),
     $         shakeatom(3,i),shakeatom(4,i),
     $         shakebondlen(1,i),shakebondlen(2,i),shakebondlen(3,i),
     $         shaketol,shakeiter,dtsq_factor,vflag)
        else
          call shake3angle(shakeatom(1,i),shakeatom(2,i),
     $         shakeatom(3,i),
     $         shakebondlen(1,i),shakebondlen(2,i),shakeanglebond,
     $         shaketol,shakeiter,dtsq_factor,vflag)
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c 2-atom SHAKE

      subroutine shake2(i0,i1,bond1,dtsq_factor,vflag)
      use global
      implicit none

c argument variables

      integer i0,i1,ifactor,vflag
      real*8 bond1,dtsq_factor

c local variables

      real*8 r01sq,s01sq
      real*8 a,b,c,determ
      real*8 lamda,lamda1,lamda2
      real*8 r01(3),s01(3)
      real*8 invmass0,invmass1

      invmass0 = 1.0/mass(type(i0))
      invmass1 = 1.0/mass(type(i1))

c r01 = distance vec between atoms, with PBC

      r01(1) = x(1,i0) - x(1,i1)
      r01(2) = x(2,i0) - x(2,i1)
      r01(3) = x(3,i0) - x(3,i1)
      call minimg(r01(1),r01(2),r01(3))

c s01 = distance vec after unconstrained update, with PBC

      s01(1) = xshake(1,i0) - xshake(1,i1)
      s01(2) = xshake(2,i0) - xshake(2,i1)
      s01(3) = xshake(3,i0) - xshake(3,i1)
      call minimg(s01(1),s01(2),s01(3))

c scalar distances between atoms

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)

c a,b,c = coeffs in quadratic equation for lamda

      a = (invmass0+invmass1)**2 * r01sq
      b = 2.0 * (invmass0+invmass1) *
     $     (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      c = s01sq - bond1*bond1

c error check

      determ = b*b - 4.0*a*c
      if (determ < 0.0) then
        write (6,*) 'WARNING: SHAKE determinant < 0.0'
        determ = 0.0
      endif

c exact quadratic solution for lamda

      lamda1 = (-b+sqrt(determ)) / (2.0*a)
      lamda2 = (-b-sqrt(determ)) / (2.0*a)

      if (abs(lamda1) <= abs(lamda2)) then
        lamda = lamda1
      else
        lamda = lamda2
      endif

c update forces if atom is owned by this processor

      lamda = lamda/dtsq_factor

      if (i0 <= nlocal) then
        f(1,i0) = f(1,i0) + lamda*r01(1)
        f(2,i0) = f(2,i0) + lamda*r01(2)
        f(3,i0) = f(3,i0) + lamda*r01(3)
      endif

      if (i1 <= nlocal) then
        f(1,i1) = f(1,i1) - lamda*r01(1)
        f(2,i1) = f(2,i1) - lamda*r01(2)
        f(3,i1) = f(3,i1) - lamda*r01(3)
      endif

      if (vflag >= 1) then
        ifactor = 0
        if (i0 <= nlocal) ifactor = ifactor + 1
        if (i1 <= nlocal) ifactor = ifactor + 1

        virial(1) = virial(1) + ifactor*0.5*r01(1)*r01(1)*lamda
        virial(2) = virial(2) + ifactor*0.5*r01(2)*r01(2)*lamda
        virial(3) = virial(3) + ifactor*0.5*r01(3)*r01(3)*lamda
        virial(4) = virial(4) + ifactor*0.5*r01(1)*r01(2)*lamda
        virial(5) = virial(5) + ifactor*0.5*r01(1)*r01(3)*lamda
        virial(6) = virial(6) + ifactor*0.5*r01(2)*r01(3)*lamda
      endif

      return
      end


c -------------------------------------------------------------------------
c 3-atom SHAKE

      subroutine shake3(i0,i1,i2,bond1,bond2,
     $     tolerance,maxiter,dtsq_factor,vflag)
      use global
      implicit none

c argument variables

      integer i0,i1,i2,maxiter,ifactor,vflag
      real*8 bond1,bond2,tolerance,dtsq_factor

c local variables

      real*8 s01sq,s02sq
      real*8 a11,a12,a21,a22,b1,b2
      real*8 a11inv,a12inv,a21inv,a22inv
      real*8 determ,determinv
      real*8 lamda01,lamda02
      real*8 r01(3),r02(3),s01(3),s02(3)
      real*8 invmass0,invmass1,invmass2
      real*8 quad1,quad2
      real*8 quad1_0101,quad1_0202,quad1_0102
      real*8 quad2_0101,quad2_0202,quad2_0102
      real*8 lamda01_new,lamda02_new
      real*8 r01sq,r02sq,r0102
      integer niter
      logical done

      invmass0 = 1.0/mass(type(i0))
      invmass1 = 1.0/mass(type(i1))
      invmass2 = 1.0/mass(type(i2))

c r01,r02 = distance vecs between atoms, with PBC

      r01(1) = x(1,i0) - x(1,i1)
      r01(2) = x(2,i0) - x(2,i1)
      r01(3) = x(3,i0) - x(3,i1)
      call minimg(r01(1),r01(2),r01(3))

      r02(1) = x(1,i0) - x(1,i2)
      r02(2) = x(2,i0) - x(2,i2)
      r02(3) = x(3,i0) - x(3,i2)
      call minimg(r02(1),r02(2),r02(3))

c s01,s02 = distance vecs after unconstrained update, with PBC

      s01(1) = xshake(1,i0) - xshake(1,i1)
      s01(2) = xshake(2,i0) - xshake(2,i1)
      s01(3) = xshake(3,i0) - xshake(3,i1)
      call minimg(s01(1),s01(2),s01(3))

      s02(1) = xshake(1,i0) - xshake(1,i2)
      s02(2) = xshake(2,i0) - xshake(2,i2)
      s02(3) = xshake(3,i0) - xshake(3,i2)
      call minimg(s02(1),s02(2),s02(3))

c scalar distances between atoms

      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)
      s02sq = s02(1)*s02(1) + s02(2)*s02(2) + s02(3)*s02(3)

c matrix coeffs and rhs for lamda equations

      a11 = 2.0 * (invmass0+invmass1) *
     $     (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      a12 = 2.0 * invmass0 *
     $     (s01(1)*r02(1) + s01(2)*r02(2) + s01(3)*r02(3))
      a21 = 2.0 * invmass0 *
     $     (s02(1)*r01(1) + s02(2)*r01(2) + s02(3)*r01(3))
      a22 = 2.0 * (invmass0+invmass2) *
     $     (s02(1)*r02(1) + s02(2)*r02(2) + s02(3)*r02(3))

      b1 = bond1*bond1 - s01sq
      b2 = bond2*bond2 - s02sq

c inverse of matrix

      determ = a11*a22 - a12*a21
      if (determ == 0.0) call abort('SHAKE determinant = 0.0')
      determinv = 1.0/determ

      a11inv = a22*determinv
      a12inv = -a12*determinv
      a21inv = -a21*determinv
      a22inv = a11*determinv

c quadratic correction coeffs

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r0102 = (r01(1)*r02(1) + r01(2)*r02(2) + r01(3)*r02(3))

      quad1_0101 = (invmass0+invmass1)**2 * r01sq
      quad1_0202 = invmass0*invmass0 * r02sq
      quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102

      quad2_0202 = (invmass0+invmass2)**2 * r02sq
      quad2_0101 = invmass0*invmass0 * r01sq
      quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102

c iterate until converged

      lamda01 = 0.0
      lamda02 = 0.0
      niter = 0
      done = .FALSE.

      do while (.not. done .and. niter < maxiter)

        quad1 = quad1_0101 * lamda01*lamda01 + 
     $       quad1_0202 * lamda02*lamda02 + 
     $       quad1_0102 * lamda01*lamda02

        quad2 = quad2_0101 * lamda01*lamda01 + 
     $       quad2_0202 * lamda02*lamda02 + 
     $       quad2_0102 * lamda01*lamda02
        
        b1 = bond1*bond1 - s01sq - quad1
        b2 = bond2*bond2 - s02sq - quad2
        
        lamda01_new = a11inv*b1 + a12inv*b2
        lamda02_new = a21inv*b1 + a22inv*b2

        done = .TRUE.
        if (abs(lamda01_new-lamda01) .gt. tolerance) done = .FALSE.
        if (abs(lamda02_new-lamda02) .gt. tolerance) done = .FALSE.

        lamda01 = lamda01_new
        lamda02 = lamda02_new
        niter = niter + 1

      enddo

c update forces if atom is owned by this processor

      lamda01 = lamda01/dtsq_factor
      lamda02 = lamda02/dtsq_factor

      if (i0 <= nlocal) then
        f(1,i0) = f(1,i0) + lamda01*r01(1) + lamda02*r02(1)
        f(2,i0) = f(2,i0) + lamda01*r01(2) + lamda02*r02(2)
        f(3,i0) = f(3,i0) + lamda01*r01(3) + lamda02*r02(3)
      endif

      if (i1 <= nlocal) then
        f(1,i1) = f(1,i1) - lamda01*r01(1)
        f(2,i1) = f(2,i1) - lamda01*r01(2)
        f(3,i1) = f(3,i1) - lamda01*r01(3)
      endif

      if (i2 <= nlocal) then
        f(1,i2) = f(1,i2) - lamda02*r02(1)
        f(2,i2) = f(2,i2) - lamda02*r02(2)
        f(3,i2) = f(3,i2) - lamda02*r02(3)
      endif

      if (vflag >= 1) then
        ifactor = 0
        if (i0 <= nlocal) ifactor = ifactor + 1
        if (i1 <= nlocal) ifactor = ifactor + 1
        if (i2 <= nlocal) ifactor = ifactor + 1

        virial(1) = virial(1) + ifactor/3.0*(r01(1)*r01(1)*lamda01)
        virial(2) = virial(2) + ifactor/3.0*(r01(2)*r01(2)*lamda01)
        virial(3) = virial(3) + ifactor/3.0*(r01(3)*r01(3)*lamda01)
        virial(4) = virial(4) + ifactor/3.0*(r01(1)*r01(2)*lamda01)
        virial(5) = virial(5) + ifactor/3.0*(r01(1)*r01(3)*lamda01)
        virial(6) = virial(6) + ifactor/3.0*(r01(2)*r01(3)*lamda01)
        virial(1) = virial(1) + ifactor/3.0*(r02(1)*r02(1)*lamda02)
        virial(2) = virial(2) + ifactor/3.0*(r02(2)*r02(2)*lamda02)
        virial(3) = virial(3) + ifactor/3.0*(r02(3)*r02(3)*lamda02)
        virial(4) = virial(4) + ifactor/3.0*(r02(1)*r02(2)*lamda02)
        virial(5) = virial(5) + ifactor/3.0*(r02(1)*r02(3)*lamda02)
        virial(6) = virial(6) + ifactor/3.0*(r02(2)*r02(3)*lamda02)
      endif

      return
      end


c -------------------------------------------------------------------------
c 4-atom SHAKE

      subroutine shake4(i0,i1,i2,i3,bond1,bond2,bond3,
     $     tolerance,maxiter,dtsq_factor,vflag)
      use global
      implicit none

c argument variables

      integer i0,i1,i2,i3,maxiter,ifactor,vflag
      real*8 bond1,bond2,bond3,tolerance,dtsq_factor

c local variables

      real*8 s01sq,s02sq,s03sq
      real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
      real*8 a11inv,a12inv,a13inv,a21inv,a22inv,a23inv
      real*8 a31inv,a32inv,a33inv
      real*8 ainv,binv,cinv,dinv,einv,finv,ginv,hinv,iinv
      real*8 determ,determinv
      real*8 lamda01,lamda02,lamda03
      real*8 r01(3),r02(3),r03(3),s01(3),s02(3),s03(3)
      real*8 invmass0,invmass1,invmass2,invmass3
      real*8 quad1,quad2,quad3
      real*8 quad1_0101,quad1_0202,quad1_0303
      real*8 quad1_0102,quad1_0103,quad1_0203
      real*8 quad2_0101,quad2_0202,quad2_0303
      real*8 quad2_0102,quad2_0103,quad2_0203
      real*8 quad3_0101,quad3_0202,quad3_0303
      real*8 quad3_0102,quad3_0103,quad3_0203
      real*8 lamda01_new,lamda02_new,lamda03_new
      real*8 r01sq,r02sq,r03sq,r0102,r0103,r0203
      integer niter
      logical done

      invmass0 = 1.0/mass(type(i0))
      invmass1 = 1.0/mass(type(i1))
      invmass2 = 1.0/mass(type(i2))
      invmass3 = 1.0/mass(type(i3))

c r01,r02,r03 = distance vecs between atoms, with PBC

      r01(1) = x(1,i0) - x(1,i1)
      r01(2) = x(2,i0) - x(2,i1)
      r01(3) = x(3,i0) - x(3,i1)
      call minimg(r01(1),r01(2),r01(3))

      r02(1) = x(1,i0) - x(1,i2)
      r02(2) = x(2,i0) - x(2,i2)
      r02(3) = x(3,i0) - x(3,i2)
      call minimg(r02(1),r02(2),r02(3))

      r03(1) = x(1,i0) - x(1,i3)
      r03(2) = x(2,i0) - x(2,i3)
      r03(3) = x(3,i0) - x(3,i3)
      call minimg(r03(1),r03(2),r03(3))

c s01,s02,s03 = distance vecs after unconstrained update, with PBC

      s01(1) = xshake(1,i0) - xshake(1,i1)
      s01(2) = xshake(2,i0) - xshake(2,i1)
      s01(3) = xshake(3,i0) - xshake(3,i1)
      call minimg(s01(1),s01(2),s01(3))

      s02(1) = xshake(1,i0) - xshake(1,i2)
      s02(2) = xshake(2,i0) - xshake(2,i2)
      s02(3) = xshake(3,i0) - xshake(3,i2)
      call minimg(s02(1),s02(2),s02(3))

      s03(1) = xshake(1,i0) - xshake(1,i3)
      s03(2) = xshake(2,i0) - xshake(2,i3)
      s03(3) = xshake(3,i0) - xshake(3,i3)
      call minimg(s03(1),s03(2),s03(3))

c scalar distances between atoms

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r03sq = r03(1)*r03(1) + r03(2)*r03(2) + r03(3)*r03(3)
      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)
      s02sq = s02(1)*s02(1) + s02(2)*s02(2) + s02(3)*s02(3)
      s03sq = s03(1)*s03(1) + s03(2)*s03(2) + s03(3)*s03(3)

c matrix coeffs and rhs for lamda equations

      a11 = 2.0 * (invmass0+invmass1) *
     $     (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      a12 = 2.0 * invmass0 *
     $     (s01(1)*r02(1) + s01(2)*r02(2) + s01(3)*r02(3))
      a13 = 2.0 * invmass0 *
     $     (s01(1)*r03(1) + s01(2)*r03(2) + s01(3)*r03(3))
      a21 = 2.0 * invmass0 *
     $     (s02(1)*r01(1) + s02(2)*r01(2) + s02(3)*r01(3))
      a22 = 2.0 * (invmass0+invmass2) *
     $     (s02(1)*r02(1) + s02(2)*r02(2) + s02(3)*r02(3))
      a23 = 2.0 * invmass0 *
     $     (s02(1)*r03(1) + s02(2)*r03(2) + s02(3)*r03(3))
      a31 = 2.0 * invmass0 *
     $     (s03(1)*r01(1) + s03(2)*r01(2) + s03(3)*r01(3))
      a32 = 2.0 * invmass0 *
     $     (s03(1)*r02(1) + s03(2)*r02(2) + s03(3)*r02(3))
      a33 = 2.0 * (invmass0+invmass3) *
     $     (s03(1)*r03(1) + s03(2)*r03(2) + s03(3)*r03(3))

      b1 = bond1*bond1 - s01sq
      b2 = bond2*bond2 - s02sq
      b3 = bond3*bond3 - s03sq

c inverse of matrix

      determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
     $     a11*a23*a32 - a12*a21*a33 - a13*a22*a31
      if (determ == 0.0) call abort('SHAKE determinant = 0.0')
      determinv = 1.0/determ

      a11inv = determinv * (a22*a33 - a23*a32)
      a12inv = -determinv * (a12*a33 - a13*a32)
      a13inv = determinv * (a12*a23 - a13*a22)
      a21inv = -determinv * (a21*a33 - a23*a31)
      a22inv = determinv * (a11*a33 - a13*a31)
      a23inv = -determinv * (a11*a23 - a13*a21)
      a31inv = determinv * (a21*a32 - a22*a31)
      a32inv = -determinv * (a11*a32 - a12*a31)
      a33inv = determinv * (a11*a22 - a12*a21)

c quadratic correction coeffs

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r03sq = r03(1)*r03(1) + r03(2)*r03(2) + r03(3)*r03(3)
      r0102 = (r01(1)*r02(1) + r01(2)*r02(2) + r01(3)*r02(3))
      r0103 = (r01(1)*r03(1) + r01(2)*r03(2) + r01(3)*r03(3))
      r0203 = (r02(1)*r03(1) + r02(2)*r03(2) + r02(3)*r03(3))

      quad1_0101 = (invmass0+invmass1)**2 * r01sq
      quad1_0202 = invmass0*invmass0 * r02sq
      quad1_0303 = invmass0*invmass0 * r03sq
      quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102
      quad1_0103 = 2.0 * (invmass0+invmass1)*invmass0 * r0103
      quad1_0203 = 2.0 * invmass0*invmass0 * r0203

      quad2_0101 = invmass0*invmass0 * r01sq
      quad2_0202 = (invmass0+invmass2)**2 * r02sq
      quad2_0303 = invmass0*invmass0 * r03sq
      quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102
      quad2_0103 = 2.0 * invmass0*invmass0 * r0103
      quad2_0203 = 2.0 * (invmass0+invmass2)*invmass0 * r0203

      quad3_0101 = invmass0*invmass0 * r01sq
      quad3_0202 = invmass0*invmass0 * r02sq
      quad3_0303 = (invmass0+invmass3)**2 * r03sq
      quad3_0102 = 2.0 * invmass0*invmass0 * r0102
      quad3_0103 = 2.0 * (invmass0+invmass3)*invmass0 * r0103
      quad3_0203 = 2.0 * (invmass0+invmass3)*invmass0 * r0203

c iterate until converged

      lamda01 = 0.0
      lamda02 = 0.0
      lamda03 = 0.0
      niter = 0
      done = .FALSE.

      do while (.not. done .and. niter < maxiter)

        quad1 = quad1_0101 * lamda01*lamda01 + 
     $       quad1_0202 * lamda02*lamda02 +
     $       quad1_0303 * lamda03*lamda03 + 
     $       quad1_0102 * lamda01*lamda02 +
     $       quad1_0103 * lamda01*lamda03 +
     $       quad1_0203 * lamda02*lamda03

        quad2 = quad2_0101 * lamda01*lamda01 + 
     $       quad2_0202 * lamda02*lamda02 +
     $       quad2_0303 * lamda03*lamda03 + 
     $       quad2_0102 * lamda01*lamda02 +
     $       quad2_0103 * lamda01*lamda03 +
     $       quad2_0203 * lamda02*lamda03

        quad3 = quad3_0101 * lamda01*lamda01 + 
     $       quad3_0202 * lamda02*lamda02 +
     $       quad3_0303 * lamda03*lamda03 + 
     $       quad3_0102 * lamda01*lamda02 +
     $       quad3_0103 * lamda01*lamda03 +
     $       quad3_0203 * lamda02*lamda03

        b1 = bond1*bond1 - s01sq - quad1
        b2 = bond2*bond2 - s02sq - quad2
        b3 = bond3*bond3 - s03sq - quad3
        
        lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3
        lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3
        lamda03_new = a31inv*b1 + a32inv*b2 + a33inv*b3
        
        done = .TRUE.
        if (abs(lamda01_new-lamda01) .gt. tolerance) done = .FALSE.
        if (abs(lamda02_new-lamda02) .gt. tolerance) done = .FALSE.
        if (abs(lamda03_new-lamda03) .gt. tolerance) done = .FALSE.

        lamda01 = lamda01_new
        lamda02 = lamda02_new
        lamda03 = lamda03_new
        niter = niter + 1

      enddo

c update forces if atom is owned by this processor

      lamda01 = lamda01/dtsq_factor
      lamda02 = lamda02/dtsq_factor
      lamda03 = lamda03/dtsq_factor

      if (i0 <= nlocal) then
        f(1,i0) = f(1,i0) +
     $       lamda01*r01(1) + lamda02*r02(1) + lamda03*r03(1)
        f(2,i0) = f(2,i0) +
     $       lamda01*r01(2) + lamda02*r02(2) + lamda03*r03(2)
        f(3,i0) = f(3,i0) +
     $       lamda01*r01(3) + lamda02*r02(3) + lamda03*r03(3)
      endif

      if (i1 <= nlocal) then
        f(1,i1) = f(1,i1) - lamda01*r01(1)
        f(2,i1) = f(2,i1) - lamda01*r01(2)
        f(3,i1) = f(3,i1) - lamda01*r01(3)
      endif

      if (i2 <= nlocal) then
        f(1,i2) = f(1,i2) - lamda02*r02(1)
        f(2,i2) = f(2,i2) - lamda02*r02(2)
        f(3,i2) = f(3,i2) - lamda02*r02(3)
      endif

      if (i3 <= nlocal) then
        f(1,i3) = f(1,i3) - lamda03*r03(1)
        f(2,i3) = f(2,i3) - lamda03*r03(2)
        f(3,i3) = f(3,i3) - lamda03*r03(3)
      endif

      if (vflag >= 1) then
        ifactor = 0
        if (i0 <= nlocal) ifactor = ifactor + 1
        if (i1 <= nlocal) ifactor = ifactor + 1
        if (i2 <= nlocal) ifactor = ifactor + 1
        if (i3 <= nlocal) ifactor = ifactor + 1
        virial(1) = virial(1) + ifactor*0.25*r01(1)*r01(1)*lamda01
        virial(2) = virial(2) + ifactor*0.25*r01(2)*r01(2)*lamda01
        virial(3) = virial(3) + ifactor*0.25*r01(3)*r01(3)*lamda01
        virial(4) = virial(4) + ifactor*0.25*r01(1)*r01(2)*lamda01
        virial(5) = virial(5) + ifactor*0.25*r01(1)*r01(3)*lamda01
        virial(6) = virial(6) + ifactor*0.25*r01(2)*r01(3)*lamda01
        virial(1) = virial(1) + ifactor*0.25*r02(1)*r02(1)*lamda02
        virial(2) = virial(2) + ifactor*0.25*r02(2)*r02(2)*lamda02
        virial(3) = virial(3) + ifactor*0.25*r02(3)*r02(3)*lamda02
        virial(4) = virial(4) + ifactor*0.25*r02(1)*r02(2)*lamda02
        virial(5) = virial(5) + ifactor*0.25*r02(1)*r02(3)*lamda02
        virial(6) = virial(6) + ifactor*0.25*r02(2)*r02(3)*lamda02
        virial(1) = virial(1) + ifactor*0.25*r03(1)*r03(1)*lamda03
        virial(2) = virial(2) + ifactor*0.25*r03(2)*r03(2)*lamda03
        virial(3) = virial(3) + ifactor*0.25*r03(3)*r03(3)*lamda03
        virial(4) = virial(4) + ifactor*0.25*r03(1)*r03(2)*lamda03
        virial(5) = virial(5) + ifactor*0.25*r03(1)*r03(3)*lamda03
        virial(6) = virial(6) + ifactor*0.25*r03(2)*r03(3)*lamda03
      endif

      return
      end


c -------------------------------------------------------------------------
c 3-atom SHAKE with fixed angle

      subroutine shake3angle(i0,i1,i2,bond1,bond2,bond12,
     $     tolerance,maxiter,dtsq_factor,vflag)
      use global
      implicit none

c argument variables

      integer i0,i1,i2,maxiter,ifactor,vflag
      real*8 bond1,bond2,bond12,tolerance,dtsq_factor

c local variables

      real*8 dtsq,s01sq,s02sq,s12sq
      real*8 a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
      real*8 a11inv,a12inv,a13inv,a21inv,a22inv,a23inv
      real*8 a31inv,a32inv,a33inv
      real*8 ainv,binv,cinv,dinv,einv,finv,ginv,hinv,iinv
      real*8 determ,determinv
      real*8 lamda01,lamda02,lamda12
      real*8 s0(3),s1(3),s2(3),s3(3)
      real*8 r01(3),r02(3),r12(3),s01(3),s02(3),s12(3)
      real*8 invmass0,invmass1,invmass2
      real*8 quad1,quad2,quad3
      real*8 quad1_0101,quad1_0202,quad1_1212
      real*8 quad1_0102,quad1_0112,quad1_0212
      real*8 quad2_0101,quad2_0202,quad2_1212
      real*8 quad2_0102,quad2_0112,quad2_0212
      real*8 quad3_0101,quad3_0202,quad3_1212
      real*8 quad3_0102,quad3_0112,quad3_0212
      real*8 lamda01_new,lamda02_new,lamda12_new
      real*8 r01sq,r02sq,r12sq,r0102,r0112,r0212
      integer niter
      logical done

      invmass0 = 1.0/mass(type(i0))
      invmass1 = 1.0/mass(type(i1))
      invmass2 = 1.0/mass(type(i2))

c r01,r02,r12 = distance vecs between atoms, with PBC

      r01(1) = x(1,i0) - x(1,i1)
      r01(2) = x(2,i0) - x(2,i1)
      r01(3) = x(3,i0) - x(3,i1)
      call minimg(r01(1),r01(2),r01(3))

      r02(1) = x(1,i0) - x(1,i2)
      r02(2) = x(2,i0) - x(2,i2)
      r02(3) = x(3,i0) - x(3,i2)
      call minimg(r02(1),r02(2),r02(3))

      r12(1) = x(1,i1) - x(1,i2)
      r12(2) = x(2,i1) - x(2,i2)
      r12(3) = x(3,i1) - x(3,i2)
      call minimg(r12(1),r12(2),r12(3))

c s01,s02,s12 = distance vecs after unconstrained update, with PBC

      s01(1) = xshake(1,i0) - xshake(1,i1)
      s01(2) = xshake(2,i0) - xshake(2,i1)
      s01(3) = xshake(3,i0) - xshake(3,i1)
      call minimg(s01(1),s01(2),s01(3))

      s02(1) = xshake(1,i0) - xshake(1,i2)
      s02(2) = xshake(2,i0) - xshake(2,i2)
      s02(3) = xshake(3,i0) - xshake(3,i2)
      call minimg(s02(1),s02(2),s02(3))

      s12(1) = xshake(1,i1) - xshake(1,i2)
      s12(2) = xshake(2,i1) - xshake(2,i2)
      s12(3) = xshake(3,i1) - xshake(3,i2)
      call minimg(s12(1),s12(2),s12(3))

c scalar distances between atoms

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r12sq = r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3)
      s01sq = s01(1)*s01(1) + s01(2)*s01(2) + s01(3)*s01(3)
      s02sq = s02(1)*s02(1) + s02(2)*s02(2) + s02(3)*s02(3)
      s12sq = s12(1)*s12(1) + s12(2)*s12(2) + s12(3)*s12(3)

c matrix coeffs and rhs for lamda equations

      a11 = 2.0 * (invmass0+invmass1) *
     $     (s01(1)*r01(1) + s01(2)*r01(2) + s01(3)*r01(3))
      a12 = 2.0 * invmass0 *
     $     (s01(1)*r02(1) + s01(2)*r02(2) + s01(3)*r02(3))
      a13 = - 2.0 * invmass1 *
     $     (s01(1)*r12(1) + s01(2)*r12(2) + s01(3)*r12(3))
      a21 = 2.0 * invmass0 *
     $     (s02(1)*r01(1) + s02(2)*r01(2) + s02(3)*r01(3))
      a22 = 2.0 * (invmass0+invmass2) *
     $     (s02(1)*r02(1) + s02(2)*r02(2) + s02(3)*r02(3))
      a23 = 2.0 * invmass2 *
     $     (s02(1)*r12(1) + s02(2)*r12(2) + s02(3)*r12(3))
      a31 = - 2.0 * invmass1 *
     $     (s12(1)*r01(1) + s12(2)*r01(2) + s12(3)*r01(3))
      a32 = 2.0 * invmass2 *
     $     (s12(1)*r02(1) + s12(2)*r02(2) + s12(3)*r02(3))
      a33 = 2.0 * (invmass1+invmass2) *
     $     (s12(1)*r12(1) + s12(2)*r12(2) + s12(3)*r12(3))

      b1 = bond1*bond1 - s01sq
      b2 = bond2*bond2 - s02sq
      b3 = bond12*bond12 - s12sq

c inverse of matrix

      determ = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
     $     a11*a23*a32 - a12*a21*a33 - a13*a22*a31
      if (determ == 0.0) call abort('SHAKE determinant = 0.0')
      determinv = 1.0/determ

      a11inv = determinv * (a22*a33 - a23*a32)
      a12inv = -determinv * (a12*a33 - a13*a32)
      a13inv = determinv * (a12*a23 - a13*a22)
      a21inv = -determinv * (a21*a33 - a23*a31)
      a22inv = determinv * (a11*a33 - a13*a31)
      a23inv = -determinv * (a11*a23 - a13*a21)
      a31inv = determinv * (a21*a32 - a22*a31)
      a32inv = -determinv * (a11*a32 - a12*a31)
      a33inv = determinv * (a11*a22 - a12*a21)

c quadratic correction coeffs

      r01sq = r01(1)*r01(1) + r01(2)*r01(2) + r01(3)*r01(3)
      r02sq = r02(1)*r02(1) + r02(2)*r02(2) + r02(3)*r02(3)
      r12sq = r12(1)*r12(1) + r12(2)*r12(2) + r12(3)*r12(3)
      r0102 = (r01(1)*r02(1) + r01(2)*r02(2) + r01(3)*r02(3))
      r0112 = (r01(1)*r12(1) + r01(2)*r12(2) + r01(3)*r12(3))
      r0212 = (r02(1)*r12(1) + r02(2)*r12(2) + r02(3)*r12(3))

      quad1_0101 = (invmass0+invmass1)**2 * r01sq
      quad1_0202 = invmass0*invmass0 * r02sq
      quad1_1212 = invmass1*invmass1 * r12sq
      quad1_0102 = 2.0 * (invmass0+invmass1)*invmass0 * r0102
      quad1_0112 = - 2.0 * (invmass0+invmass1)*invmass1 * r0112
      quad1_0212 = - 2.0 * invmass0*invmass1 * r0212

      quad2_0101 = invmass0*invmass0 * r01sq
      quad2_0202 = (invmass0+invmass2)**2 * r02sq
      quad2_1212 = invmass2*invmass2 * r12sq
      quad2_0102 = 2.0 * (invmass0+invmass2)*invmass0 * r0102
      quad2_0112 = 2.0 * invmass0*invmass2 * r0112
      quad2_0212 = 2.0 * (invmass0+invmass2)*invmass2 * r0212

      quad3_0101 = invmass1*invmass1 * r01sq
      quad3_0202 = invmass2*invmass2 * r02sq
      quad3_1212 = (invmass1+invmass2)**2 * r12sq
      quad3_0102 = - 2.0 * invmass1*invmass2 * r0102
      quad3_0112 = - 2.0 * (invmass1+invmass2)*invmass1 * r0112
      quad3_0212 = 2.0 * (invmass1+invmass2)*invmass2 * r0212

c iterate until converged

      lamda01 = 0.0
      lamda02 = 0.0
      lamda12 = 0.0
      niter = 0
      done = .FALSE.

      do while (.not. done .and. niter < maxiter)

        quad1 = quad1_0101 * lamda01*lamda01 + 
     $       quad1_0202 * lamda02*lamda02 +
     $       quad1_1212 * lamda12*lamda12 + 
     $       quad1_0102 * lamda01*lamda02 +
     $       quad1_0112 * lamda01*lamda12 +
     $       quad1_0212 * lamda02*lamda12

        quad2 = quad2_0101 * lamda01*lamda01 + 
     $       quad2_0202 * lamda02*lamda02 +
     $       quad2_1212 * lamda12*lamda12 + 
     $       quad2_0102 * lamda01*lamda02 +
     $       quad2_0112 * lamda01*lamda12 +
     $       quad2_0212 * lamda02*lamda12

        quad3 = quad3_0101 * lamda01*lamda01 + 
     $       quad3_0202 * lamda02*lamda02 +
     $       quad3_1212 * lamda12*lamda12 + 
     $       quad3_0102 * lamda01*lamda02 +
     $       quad3_0112 * lamda01*lamda12 +
     $       quad3_0212 * lamda02*lamda12

        b1 = bond1*bond1 - s01sq - quad1
        b2 = bond2*bond2 - s02sq - quad2
        b3 = bond12*bond12 - s12sq - quad3
        
        lamda01_new = a11inv*b1 + a12inv*b2 + a13inv*b3
        lamda02_new = a21inv*b1 + a22inv*b2 + a23inv*b3
        lamda12_new = a31inv*b1 + a32inv*b2 + a33inv*b3

        done = .TRUE.
        if (abs(lamda01_new-lamda01) .gt. tolerance) done = .FALSE.
        if (abs(lamda02_new-lamda02) .gt. tolerance) done = .FALSE.
        if (abs(lamda12_new-lamda12) .gt. tolerance) done = .FALSE.

        lamda01 = lamda01_new
        lamda02 = lamda02_new
        lamda12 = lamda12_new
        niter = niter + 1

      enddo

c update forces if atom is owned by this processor

      lamda01 = lamda01/dtsq_factor
      lamda02 = lamda02/dtsq_factor
      lamda12 = lamda12/dtsq_factor

      if (i0 <= nlocal) then
        f(1,i0) = f(1,i0) + lamda01*r01(1) + lamda02*r02(1)
        f(2,i0) = f(2,i0) + lamda01*r01(2) + lamda02*r02(2)
        f(3,i0) = f(3,i0) + lamda01*r01(3) + lamda02*r02(3)
      endif

      if (i1 <= nlocal) then
        f(1,i1) = f(1,i1) - lamda01*r01(1) + lamda12*r12(1)
        f(2,i1) = f(2,i1) - lamda01*r01(2) + lamda12*r12(2)
        f(3,i1) = f(3,i1) - lamda01*r01(3) + lamda12*r12(3)
      endif

      if (i2 <= nlocal) then
        f(1,i2) = f(1,i2) - lamda02*r02(1) - lamda12*r12(1)
        f(2,i2) = f(2,i2) - lamda02*r02(2) - lamda12*r12(2)
        f(3,i2) = f(3,i2) - lamda02*r02(3) - lamda12*r12(3)
      endif

      if (vflag >= 1) then
        ifactor = 0
        if (i0 <= nlocal) ifactor = ifactor + 1
        if (i1 <= nlocal) ifactor = ifactor + 1
        if (i2 <= nlocal) ifactor = ifactor + 1
        virial(1) = virial(1) + ifactor/3.0*(r01(1)*r01(1)*lamda01)
        virial(2) = virial(2) + ifactor/3.0*(r01(2)*r01(2)*lamda01)
        virial(3) = virial(3) + ifactor/3.0*(r01(3)*r01(3)*lamda01)
        virial(4) = virial(4) + ifactor/3.0*(r01(1)*r01(2)*lamda01)
        virial(5) = virial(5) + ifactor/3.0*(r01(1)*r01(3)*lamda01)
        virial(6) = virial(6) + ifactor/3.0*(r01(2)*r01(3)*lamda01)
        virial(1) = virial(1) + ifactor/3.0*(r02(1)*r02(1)*lamda02)
        virial(2) = virial(2) + ifactor/3.0*(r02(2)*r02(2)*lamda02)
        virial(3) = virial(3) + ifactor/3.0*(r02(3)*r02(3)*lamda02)
        virial(4) = virial(4) + ifactor/3.0*(r02(1)*r02(2)*lamda02)
        virial(5) = virial(5) + ifactor/3.0*(r02(1)*r02(3)*lamda02)
        virial(6) = virial(6) + ifactor/3.0*(r02(2)*r02(3)*lamda02)
        virial(1) = virial(1) + ifactor/3.0*(r12(1)*r12(1)*lamda12)
        virial(2) = virial(2) + ifactor/3.0*(r12(2)*r12(2)*lamda12)
        virial(3) = virial(3) + ifactor/3.0*(r12(3)*r12(3)*lamda12)
        virial(4) = virial(4) + ifactor/3.0*(r12(1)*r12(2)*lamda12)
        virial(5) = virial(5) + ifactor/3.0*(r12(1)*r12(3)*lamda12)
        virial(6) = virial(6) + ifactor/3.0*(r12(2)*r12(3)*lamda12)
      endif

      return
      end


c -------------------------------------------------------------------------
c create SHAKE list for each atom

      subroutine shake_create
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,m
      integer inext,iprev,iatom,jatom,nbuf,nbufall,ipnt
      integer mess_size,mess_tag,itmp,itmp2,itmp3,itmp4,itmpangle
      integer ibad,ibadall,num,ierror,irequest,loop,maxshake
      integer iatom_local,jatom_local,ibond,iflag
      integer istatus(mpi_status_size)
      integer, allocatable :: buf(:),bufcopy(:)
      real*8 bond,angle

      call mpi_barrier(mpi_comm_world,ierror)
      if (node == 0) write (6,*) 'Finding SHAKE neighbors ...'

c shakewhichbondcoeff = location in bond coeffs where bond length resides
c FENE bondstyle not allowed

      if (bondstyle == 1) shakewhichbondcoeff = 2
      if (bondstyle == 4) shakewhichbondcoeff = 2
      if (bondstyle == 5) shakewhichbondcoeff = 1
      if (bondstyle == 2 .or. bondstyle == 3)
     $     call err('No current SHAKE support for FENE bonds')

c shakeanglebond = pseudo-bond distance due to angle constraint

      if (shakeableangle > 0) then
        if (anglestyle == 1) angle = anglecoeff(2,shakeableangle)
        if (anglestyle == 2) angle = anglecoeff(1,shakeableangle)
        if (anglestyle == 3) angle = anglecoeff(2,shakeableangle)
        bond = bondcoeff(shakewhichbondcoeff,shakeableanglebond)
        shakeanglebond = 2.0*bond*bond * (1.0-cos(angle))
        shakeanglebond = sqrt(shakeanglebond)
      endif

c processor IDs for cyclic ring passing

      inext = node + 1
      iprev = node - 1
      if (inext == nprocs) inext = 0
      if (iprev < 0) iprev = nprocs - 1

c ---------------------------------------------------------
c make list of SHAKE partner atoms bonded to each atom
c at end of this section:
c   shakegroup = # of SHAKE bonds attached to each atom
c   shakepartner = global ID of SHAKE partners for each atom (no self)
c ---------------------------------------------------------

c for newton off, all bonds are stored with atom

      if (newton_bond == 0) then
        
        do i = 1,nlocal
          shakegroup(i) = 0
          do j = 1,numbond(i)
            if (shakeablebond(bondtype(j,i)) == 1) then
              shakegroup(i) = shakegroup(i) + 1
              if (bondatom1(j,i) /= tag(i)) jatom = bondatom1(j,i)
              if (bondatom2(j,i) /= tag(i)) jatom = bondatom2(j,i)
              if (shakegroup(i) <= 3)
     $             shakepartner(shakegroup(i),i) = jatom
            endif
          enddo
        enddo

      else

c for newton on, only 1/2 of bonds are stored with atom
c increment shakegroup for other atom if I own it
c nbuf = # of off-proc bond partner info I'll need to communicate

        do i = 1,nlocal
          shakegroup(i) = 0
        enddo

        nbuf = 0
        do i = 1,nlocal
          do j = 1,numbond(i)
            if (shakeablebond(bondtype(j,i)) == 1) then
              shakegroup(i) = shakegroup(i) + 1
              jatom = bondatom2(j,i)
              if (shakegroup(i) <= 3)
     $             shakepartner(shakegroup(i),i) = jatom
              m = localptr(jatom)
              if (m /= 0 .and. m <= nlocal) then
                shakegroup(m) = shakegroup(m) + 1
                if (shakegroup(m) <= 3)
     $               shakepartner(shakegroup(m),m) = tag(i)
              else
                nbuf = nbuf + 1
              endif
            endif
          enddo
        enddo

c nbufall = largest amount of info any proc will communicate

        call mpi_allreduce(2*nbuf,nbufall,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)

        if (nbufall > 0) then

          allocate(buf(nbufall))
          allocate(bufcopy(nbufall))

c load buffer with pairs of global IDs of bonded atoms
c 2nd atom is off-processor

          ipnt = 0
          do i = 1,nlocal
            do j = 1,numbond(i)
              if (shakeablebond(bondtype(j,i)) == 1) then
                jatom = bondatom2(j,i)
                m = localptr(jatom)
                if (m == 0 .or. m > nlocal) then
                  buf(ipnt+1) = tag(i)
                  buf(ipnt+2) = jatom
                  ipnt = ipnt + 2
                endif
              endif
            enddo
          enddo

c cycle buffer around to all nodes in handshaked fashion
c do not send/recv buffer to self, no need to process own buf
c if I own 2nd atom, increment its shakegroup

          mess_size = ipnt
          mess_tag = 1

          do loop = 1,nprocs-1
            call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $           mpi_comm_world,irequest,ierror)
            call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
            do i = 1,mess_size
              buf(i) = bufcopy(i)
            enddo
            do i = 2,mess_size,2
              j = localptr(buf(i))
              if (j /= 0 .and. j <= nlocal) then
                shakegroup(j) = shakegroup(j) + 1
                if (shakegroup(j) <= 3)
     $               shakepartner(shakegroup(j),j) = buf(i-1)
              endif
            enddo
          enddo

c free buffers

          deallocate(buf,bufcopy)

        endif

      endif

c -------------------------------------------------
c maxshake = size of largest SHAKE group = shakegroup() + 1
c error check: should not be > 4 for current SHAKE implementation
c if maxshake = 1, just return since SHAKE is effectively non-active
c -------------------------------------------------

      itmp = 0
      do i = 1,nlocal
        itmp = max(itmp,shakegroup(i))
      enddo

      call mpi_allreduce(itmp+1,maxshake,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (maxshake > 4) then
        if (node == 0) write (6,*) 'Largest SHAKE cluster =',maxshake
        call err('Too large a SHAKE cluster')
      endif

      if (maxshake == 1) return

c -------------------------------------------------
c SHAKE groups must be independent
c check that no 2 atoms are connected which both have shakegroup > 1
c -------------------------------------------------

c nbuf = size of buffer needed to hold info from this proc
c info for each atom = shakegroup, shakepartner
c only necessary for groups with shakegroup > 1

      nbuf = 0
      do i = 1,nlocal
        if (shakegroup(i) > 1) nbuf = nbuf + 1 + shakegroup(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      allocate(buf(nbufall))
      allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = # of SHAKE partners
c  (2:N) = global tags of each SHAKE partner

      ipnt = 0
      do i = 1,nlocal
        if (shakegroup(i) > 1) then
          ipnt = ipnt + 1
          buf(ipnt) = shakegroup(i)
          do j = 1,shakegroup(i)
            ipnt = ipnt + 1
            buf(ipnt) = shakepartner(j,i)
          enddo
        endif
      enddo

c cycle buffer around to all nodes in handshaked fashion
c do not send/recv buffer to self, but process own buf
c for each received SHAKE group, check each SHAKE partner (1:num)
c if I own partner, insure its shakegroup <= 1, else error

      mess_size = ipnt
      mess_tag = 2

      ibad = 0
      do loop = 1,nprocs
        i = 1
        do while (i <= mess_size)
          num = buf(i)
          do j = 1,num
            m = localptr(buf(i+j))
            if (m /= 0 .and. m <= nlocal) then
              if (shakegroup(m) > 1) ibad = ibad + 1
            endif
          enddo
          i = i + 1 + num
        enddo
        if (node /= inext) then
          call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
          do i = 1,mess_size
            buf(i) = bufcopy(i)
          enddo
        endif
      enddo

c free buffers

      deallocate(buf,bufcopy)

c error check

      call mpi_allreduce(ibad,ibadall,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)

      if (ibadall > 0) then
        if (node == 0)
     $       write (6,*) '# of SHAKE group connections =',ibadall
        call err('Two or more SHAKE groups are connected')
      endif

c -------------------------------------------------
c add self to beginning of every SHAKE group
c -------------------------------------------------

c increment size of SHAKE group to include self
c put 1st partner at end of list to make room for self

      do i = 1,nlocal
        if (shakegroup(i) > 0) then
          shakegroup(i) = shakegroup(i) + 1
          shakepartner(shakegroup(i),i) = shakepartner(1,i)
          shakepartner(1,i) = tag(i)
        endif
      enddo

c -------------------------------------------------
c synchronize all atoms in each SHAKE group, so they store same info
c each will know size of entire group
c   and store central atom 1st in shakepartner list
c SHAKE groups of size 2 are an exception, both atoms can consider
c   themselves as central atom, OK to store shakepartner in different orders
c -------------------------------------------------

c nbuf = size of buffer needed to hold info from this proc
c info for each atom = shakegroup, shakepartner
c only necessary for groups with shakegroup > 2

      nbuf = 0
      do i = 1,nlocal
        if (shakegroup(i) > 2) nbuf = nbuf + 1 + shakegroup(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (nbufall > 0) then

        allocate(buf(nbufall))
        allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = size of SHAKE group
c  (2:N) = global tags of each atom in SHAKE group

        ipnt = 0
        do i = 1,nlocal
          if (shakegroup(i) > 2) then
            ipnt = ipnt + 1
            buf(ipnt) = shakegroup(i)
            do j = 1,shakegroup(i)
              ipnt = ipnt + 1
              buf(ipnt) = shakepartner(j,i)
            enddo
          endif
        enddo

c cycle buffer around to all nodes in handshaked fashion
c do not send/recv buffer to self, but process own buf
c for each received SHAKE group, check each non-central SHAKE partner (2:num)
c if I own partner, then it is not a central atom so overwrite that
c   atom's SHAKE group size and partner list with received info
        
        mess_size = ipnt
        mess_tag = 3

        do loop = 1,nprocs
          i = 1
          do while (i <= mess_size)
            num = buf(i)
            do j = 2,num
              m = localptr(buf(i+j))
              if (m /= 0 .and. m <= nlocal) then
                shakegroup(m) = num
                do k = 1,num
                  shakepartner(k,m) = buf(i+k)
                enddo
              endif
            enddo
            i = i + 1 + num
          enddo
          if (node /= inext) then
            call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $           mpi_comm_world,irequest,ierror)
            call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $           mpi_comm_world,ierror)
            call mpi_wait(irequest,istatus,ierror)
            call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
            do i = 1,mess_size
              buf(i) = bufcopy(i)
            enddo
          endif
        enddo

c free buffers

        deallocate(buf,bufcopy)

      endif

c -------------------------------------------------
c create shakebondtype = which global bond each SHAKE bond is
c -------------------------------------------------

c zero out the data structure

      do i = 1,nlocal
        do j = 1,shakegroup(i)-1
          shakebondtype(j,i) = 0
        enddo
      enddo

c fill shakebondtype with bonds I know about
c iatom/jatom = 2 global IDs of atoms in a SHAKE bond
c iatom_local/jatom_local = local IDs of those atoms
c if I own iatom or jatom, scan their bond lists
c if the bond to the other atom is there, extract the bondtype

      do i = 1,nlocal
        do j = 2,shakegroup(i)
          iatom = shakepartner(1,i)
          jatom = shakepartner(j,i)
          iatom_local = localptr(iatom)
          jatom_local = localptr(jatom)
          if (iatom_local /= 0 .and. iatom_local <= nlocal) then
            do k = 1,numbond(iatom_local)
              if ((bondatom1(k,iatom_local) == iatom .and.
     $             bondatom2(k,iatom_local) == jatom) .or.
     $             (bondatom1(k,iatom_local) == jatom .and.
     $             bondatom2(k,iatom_local) == iatom)) then
                shakebondtype(j-1,i) = bondtype(k,iatom_local)
                exit
              endif
            enddo
          endif
          if (jatom_local /= 0 .and. jatom_local <= nlocal) then
            do k = 1,numbond(jatom_local)
              if ((bondatom1(k,jatom_local) == iatom .and.
     $             bondatom2(k,jatom_local) == jatom) .or.
     $             (bondatom1(k,jatom_local) == jatom .and.
     $             bondatom2(k,jatom_local) == iatom)) then
                shakebondtype(j-1,i) = bondtype(k,jatom_local)
                exit
              endif
            enddo
          endif
        enddo
      enddo

c nbuf = size of buffer needed to hold info from this proc
c info for each atom = shakegroup, shakepartner, shakebondtype
c only necessary for groups with some atoms off proc
c if nbufall = 0 -> all shakebondtype values are set, skip to end

      nbuf = 0
      do i = 1,nlocal
        iflag = 0
        do j = 1,shakegroup(i)
          m = localptr(shakepartner(j,i))
          if (m == 0 .or. m > nlocal) iflag = 1
        enddo
        if (iflag == 1) nbuf = nbuf + 2*shakegroup(i)
      enddo

      call mpi_allreduce(nbuf,nbufall,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (nbufall > 0) then

        allocate(buf(nbufall))
        allocate(bufcopy(nbufall))

c fill buffer with:
c  (1) = size of SHAKE group
c  (2:N) = global tags of each atom in SHAKE group
c  (N+1:2*N) = bondtype of each bond in SHAKE group
      
        ipnt = 0
        do i = 1,nlocal
          iflag = 0
          do j = 1,shakegroup(i)
            m = localptr(shakepartner(j,i))
            if (m == 0 .or. m > nlocal) iflag = 1
          enddo
          if (iflag == 1) then
            ipnt = ipnt + 1
            buf(ipnt) = shakegroup(i)
            do j = 1,shakegroup(i)
              ipnt = ipnt + 1
              buf(ipnt) = shakepartner(j,i)
            enddo
            do j = 1,shakegroup(i)-1
              ipnt = ipnt + 1
              buf(ipnt) = shakebondtype(j,i)
            enddo
          endif
        enddo

c cycle buffer around to all nodes in handshaked fashion
c do not send/recv buffer to self, no need to process own buf
c for each received SHAKE group, check each SHAKE partner (1:num)
c if I own partner, its shakebondtype list matches received one
c so use any nonzero received shakebondtype values to update my list

        mess_size = ipnt
        mess_tag = 4
        
        do loop = 1,nprocs-1
          call mpi_irecv(bufcopy,nbufall,mpi_integer,iprev,mess_tag,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf,mess_size,mpi_integer,inext,mess_tag,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
          call mpi_get_count(istatus,mpi_integer,mess_size,ierror)
          do i = 1,mess_size
            buf(i) = bufcopy(i)
          enddo
          
          i = 1
          do while (i <= mess_size)
            num = buf(i)
            do j = 1,num
              m = localptr(buf(i+j))
              if (m /= 0 .and. m <= nlocal) then
                do k = 1,num-1
                  ibond = buf(i+num+k)
                  if (ibond /= 0) shakebondtype(k,m) = ibond
                enddo
              endif
            enddo
            i = i + 2*num
          enddo
        enddo
        
c error check that all shakebondtype values are non-zero
        
        ibad = 0
        do i = 1,nlocal
          do j = 1,shakegroup(i)-1
            if (shakebondtype(j,i) == 0) ibad = ibad + 1
          enddo
        enddo
        
        call mpi_allreduce(ibad,ibadall,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)
        
        if (ibadall > 0) then
          if (node == 0)
     $         write (6,*) '# of zero-value SHAKE bondtypes = ',
     $         ibadall
          call err('Invalid SHAKE bondtypes')
        endif

      endif
        
c -------------------------------------------------
c print info on size of SHAKE groups
c -------------------------------------------------

      itmp2 = 0
      itmp3 = 0
      itmp4 = 0
      itmpangle = 0
      do i = 1,nlocal
        if (shakegroup(i) == 2) itmp2 = itmp2 + 1
        if (shakegroup(i) == 3) then
          if (shakeableangle > 0 .and.
     $         shakebondtype(1,i) == shakeableanglebond .and.
     $         shakebondtype(2,i) == shakeableanglebond) then
            itmpangle = itmpangle + 1
          else
            itmp3 = itmp3 + 1
          endif
        endif
        if (shakegroup(i) == 4) itmp4 = itmp4 + 1
      enddo
      
      itmp = itmp2
      call mpi_allreduce(itmp,itmp2,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (node == 0) then
        write (6,*) 'SHAKE clusters of size 2 =',itmp2/2
        write (1,*) 'SHAKE clusters of size 2 =',itmp2/2
      endif

      itmp = itmp3
      call mpi_allreduce(itmp,itmp3,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (node == 0) then
        write (6,*) 'SHAKE clusters of size 3 =',itmp3/3
        write (1,*) 'SHAKE clusters of size 3 =',itmp3/3
      endif

      itmp = itmp4
      call mpi_allreduce(itmp,itmp4,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (node == 0) then
        write (6,*) 'SHAKE clusters of size 4 =',itmp4/4
        write (1,*) 'SHAKE clusters of size 4 =',itmp4/4
      endif

      itmp = itmpangle
      call mpi_allreduce(itmp,itmpangle,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (node == 0) then
        write (6,*) 'SHAKE 3-clusters with frozen angles =',itmpangle/3
        write (1,*) 'SHAKE 3-clusters with frozen angles =',itmpangle/3
      endif

      nshakebonds = itmp2/2 + 2*itmp3/3 + 3*itmp4/4 + itmpangle

      return
      end


c -----------------------------------------------------------------------
c allocate SHAKE memory - done at start of every run

      subroutine shake_allocate
      use global
      implicit none

      allocate(shakegroup(maxown))
      allocate(shakepartner(4,maxown))
      allocate(shakebondtype(3,maxown))
      allocate(shakesize(maxown))
      allocate(shakeatom(4,maxown))
      allocate(shakebondlen(3,maxown))
      allocate(xshake(3,maxatom))

      return
      end


c -----------------------------------------------------------------------
c deallocate SHAKE memory - done at end of every run

      subroutine shake_deallocate
      use global
      implicit none

      deallocate(shakegroup,shakepartner,shakebondtype)
      deallocate(shakesize,shakeatom,shakebondlen)
      deallocate(xshake)

      return
      end


c -----------------------------------------------------------------------
c SHAKE update to unconstrained coords
c xshake = predicted new coords using current velocity/force
c assumes NVE update, seems to be accurate enough for NVT,NPH,NPT as well

      subroutine shake_update(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variablest

      integer i
      real*8 invmass,dtsq_factor

      if (iflag == 0) dtsq_factor = 0.5*dt*dt
      if (iflag == 1) dtsq_factor = dt*dt

      do i = 1,nlocal
        invmass = 1.0/mass(type(i))
        xshake(1,i) = x(1,i) + dt*v(1,i) + dtsq_factor*invmass*f(1,i)
        xshake(2,i) = x(2,i) + dt*v(2,i) + dtsq_factor*invmass*f(2,i)
        xshake(3,i) = x(3,i) + dt*v(3,i) + dtsq_factor*invmass*f(3,i)
      enddo

      return
      end


c -----------------------------------------------------------------------
c SHAKE interprocessor communication of unconstrained coords via 6-way stencil
c swap slabs of atoms with other processors in all 3 directions
c send/receive message always to another proc - never to self
      
      subroutine shake_comm
      use global
      use mpi
      implicit none

c local variables

      integer i,j,k,m,ii,ierror
      integer iflag,jflag,kflag,isize_send,isize_recv
      integer irequest,istatus(mpi_status_size)

      do k = 1,nswap

c sending to another proc, pack up message buffer with current atom coords

        if (spart(k) /= node) then

          j = 0

c not sending across periodic boundary

          if (commflagall(k) == 0) then
            do ii = nsfirst(k),nslast(k)
              i = slist(ii)
              buf1(j+1) = xshake(1,i)
              buf1(j+2) = xshake(2,i)
              buf1(j+3) = xshake(3,i)
              j = j + 3
            enddo

          else

c sending across periodic boundary, add/subtract box length

            iflag = commflag(1,k)
            jflag = commflag(2,k)
            kflag = commflag(3,k)
            do ii = nsfirst(k),nslast(k)
              i = slist(ii)
              buf1(j+1) = xshake(1,i) + iflag*xprd
              buf1(j+2) = xshake(2,i) + jflag*yprd
              buf1(j+3) = xshake(3,i) + kflag*zprd
              j = j + 3
            enddo

          endif

c send/receive message

          isize_send = nslast(k) - nsfirst(k) + 1
          isize_recv = nrlast(k) - nrfirst(k) + 1

#ifdef SENDRECV
          call mpi_sendrecv(buf1,3*isize_send,
     $         mpi_double_precision,spart(k),0,
     $         xshake(1,nrfirst(k)),3*isize_recv,
     $         mpi_double_precision,rpart(k),0,
     $         mpi_comm_world,istatus,ierror)
#endif
#ifndef SENDRECV
          call mpi_irecv(xshake(1,nrfirst(k)),3*isize_recv,
     $         mpi_double_precision,rpart(k),0,
     $         mpi_comm_world,irequest,ierror)
          call mpi_send(buf1,3*isize_send,
     $         mpi_double_precision,spart(k),0,
     $         mpi_comm_world,ierror)
          call mpi_wait(irequest,istatus,ierror)
#endif

        else

c only swapping with self, avoid buffering, do direct copy

          m = nrfirst(k)
          iflag = commflag(1,k)
          jflag = commflag(2,k)
          kflag = commflag(3,k)
          do ii = nsfirst(k),nslast(k)
            i = slist(ii)
            xshake(1,m) = xshake(1,i) + iflag*xprd
            xshake(2,m) = xshake(2,i) + jflag*yprd
            xshake(3,m) = xshake(3,i) + kflag*zprd
            m = m + 1
          enddo

        endif
      enddo
      
      return
      end


c -----------------------------------------------------------------------
c print-out bond length statistics

      subroutine shake_bond_stats
      use global
      use mpi
      implicit none

c local variables

      integer i,j,m,nb1,iatom,jatom,ierror
      real*8 delx,dely,delz,dist,dist1,dist2,dist3,angle
      integer, allocatable :: ncount(:),ncount_all(:)
      real*8, allocatable :: distave(:),distave_all(:)
      real*8, allocatable :: distmax(:),distmax_all(:)
      real*8, allocatable :: distmin(:),distmin_all(:)

c statistics arrays

      nb1 = nbondtypes + 1
      allocate(ncount(nb1),ncount_all(nb1))
      allocate(distave(nb1),distave_all(nb1))
      allocate(distmax(nb1),distmax_all(nb1))
      allocate(distmin(nb1),distmin_all(nb1))

      do i = 1,nb1
        ncount(i) = 0
        distave(i) = 0.0
        distmax(i) = 0.0
        distmin(i) = 1.0E20
      enddo

c log stats for each bond
c OK to double count since are just averaging
c   and both numerator/denominator will be doubled

      do i = 1,nlocal
        if (shakegroup(i) > 0) then
          iatom = localptr(shakepartner(1,i))
          do j = 2,shakegroup(i)
            jatom = localptr(shakepartner(j,i))
            delx = x(1,iatom) - x(1,jatom)
            dely = x(2,iatom) - x(2,jatom)
            delz = x(3,iatom) - x(3,jatom)
            call minimg(delx,dely,delz)
            dist = sqrt(delx*delx + dely*dely + delz*delz)

            m = shakebondtype(j-1,i)
            ncount(m) = ncount(m) + 1
            distave(m) = distave(m) + dist
            if (dist > distmax(m)) distmax(m) = dist
            if (dist < distmin(m)) distmin(m) = dist
          enddo
          if (shakeableangle > 0 .and. shakegroup(i) == 3) then
            if (shakebondtype(1,i) == shakeableanglebond .and.
     $           shakebondtype(2,i) == shakeableanglebond) then
              dist2 = dist

              jatom = localptr(shakepartner(2,i))
              delx = x(1,iatom) - x(1,jatom)
              dely = x(2,iatom) - x(2,jatom)
              delz = x(3,iatom) - x(3,jatom)
              call minimg(delx,dely,delz)
              dist1 = sqrt(delx*delx + dely*dely + delz*delz)

              iatom = localptr(shakepartner(3,i))
              delx = x(1,iatom) - x(1,jatom)
              dely = x(2,iatom) - x(2,jatom)
              delz = x(3,iatom) - x(3,jatom)
              call minimg(delx,dely,delz)
              dist3 = sqrt(delx*delx + dely*dely + delz*delz)
            
              angle = acos((dist1*dist1 + dist2*dist2 - dist3*dist3) /
     $             (2.0*dist1*dist2))
              angle = angle * 180.0/3.1415926
              ncount(nb1) = ncount(nb1) + 1
              distave(nb1) = distave(nb1) + angle
              if (angle > distmax(nb1)) distmax(nb1) = angle
              if (angle < distmin(nb1)) distmin(nb1) = angle
            endif
          endif
        endif
      enddo

      call mpi_allreduce(ncount,ncount_all,nb1,
     $     mpi_integer,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(distave,distave_all,nb1,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(distmax,distmax_all,nb1,
     $     mpi_double_precision,mpi_max,mpi_comm_world,ierror)
      call mpi_allreduce(distmin,distmin_all,nb1,
     $     mpi_double_precision,mpi_min,mpi_comm_world,ierror)

      if (node == 0) then
        write (6,*)
     $       'SHAKE bond distance stats ',
     $       '(type,ave,delta) on timestep',ntimestep
        do i = 1,nbondtypes
          if (ncount_all(i) > 0 .and. shakeablebond(i) == 1)
     $         write (6,*) '  ',i,distave_all(i)/ncount_all(i),
     $         distmax_all(i)-distmin_all(i)
        enddo
        if (shakeableangle > 0 .and. ncount_all(nb1) > 0) then
          write (6,*) '      angle stats ',
     $         '(type,ave,delta) on timestep',ntimestep
          write (6,*) '  ',shakeableangle,
     $         distave_all(nb1)/ncount_all(nb1),
     $         distmax_all(nb1)-distmin_all(nb1)
        endif

        write (1,*)
     $       'SHAKE bond distance stats ',
     $       '(type,ave,delta) on timestep',ntimestep
        do i = 1,nbondtypes
          if (ncount_all(i) > 0 .and. shakeablebond(i) == 1)
     $         write (1,*) '  ',i,distave_all(i)/ncount_all(i),
     $         distmax_all(i)-distmin_all(i)
        enddo
        if (shakeableangle > 0 .and. ncount_all(nb1) > 0) then
          write (1,*) '      angle stats ',
     $         '(type,ave,delta) on timestep',ntimestep
          write (1,*) '  ',shakeableangle,
     $         distave_all(nb1)/ncount_all(nb1),
     $         distmax_all(nb1)-distmin_all(nb1)
        endif
      endif

      deallocate(ncount,distave,distmax,distmin)
      deallocate(ncount_all,distave_all,distmax_all,distmin_all)

      nshake_next = min(ntimestep+nshakestats,ntime_last)

      return
      end
