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
c non-bonded forces
c  force components stored for atoms I own
c  variable force = (-1/r * d(phi)/dr)
c  variable virial = (-r * d(phi)/dr)

c cut LJ and no Coulombic

      subroutine lj_cut(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 xtmp,ytmp,ztmp,spfactor,delx,dely,delz,rsq
      real*8 r2inv,r6inv,forcelj,force

      do i = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then
            r2inv = 1.0/rsq
            r6inv = r2inv*r2inv*r2inv
            forcelj = r6inv*(lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            force = spfactor*forcelj*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and cut Coulombic

      subroutine lj_cut_coul_cut(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 rsq,r2inv,forcecoul,r6inv,forcelj,force
      real*8 spfactor_coul,spfactor_lj

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich)
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else
              forcecoul = 0.0
            endif

            if (rsq < cutljsq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (spfactor_coul*forcecoul +
     $           spfactor_lj*forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and cut Coulombic

      subroutine lj_cut_coul_debye(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 rsq,r2inv,forcecoul,r6inv,forcelj,force
      real*8 r,rinv,screening
      real*8 spfactor_coul,spfactor_lj

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich)
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              r = sqrt(rsq)
              rinv = 1.0/r
              screening = exp(-kappa*r)
              forcecoul = coulpre*qtmp*q(j)*screening*(kappa + rinv)
            else
              forcecoul = 0.0
            endif

            if (rsq < cutljsq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (spfactor_coul*forcecoul +
     $           spfactor_lj*forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and smooth Coulombic

      subroutine lj_cut_coul_smooth(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 spfactor_coul,spfactor_lj
      real*8 rsq,r2inv,forcecoul,r,dd1,switch_r,deriv_sw_r
      real*8 r6inv,forcelj,force

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich)
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulintsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else if (rsq < cutcoulsq) then
              r = sqrt(rsq)
              dd1 = (cutcoul - cutcoulint)**3
              switch_r = (cutcoul - r)**2 *
     $             (cutcoul + (2.0D0*r) - (3.0D0*cutcoulint)) / dd1
              deriv_sw_r = 6.0D0*(r - cutcoul)*(r - cutcoulint) / dd1
              forcecoul = (coulpre *qtmp*q(j)*sqrt(r2inv)*switch_r) -
     $             (coulpre*qtmp*q(j)*deriv_sw_r)
            else
              forcecoul = 0.0D0
            endif

            if (rsq < cutljsq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (spfactor_coul*forcecoul + 
     $           spfactor_lj*forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut LJ and long-range Coulombic

      subroutine lj_cut_coul_long(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4
      real*8 ewald_a5,qtmp,xtmp,ytmp,ztmp
      real*8 spfactor_coul,spfactor_lj
      real*8 delx,dely,delz,rsq,r2inv,r,grij,expm2,t,erfc
      real*8 prefactor,forcecoul,r6inv,forcelj,force

      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich) - 2.0
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              r = sqrt(rsq)
              grij = gewald*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)/r
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor_coul < 1.0) forcecoul = forcecoul -
     $             (1.0-spfactor_coul)*prefactor
            else
              forcecoul = 0.0
            endif

            if (rsq < cutljsq(itype,jtype) .and.
     $           spfactor_lj > 0.0) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = spfactor_lj *
     $             r6inv*(lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and no Coulombic

      subroutine lj_smooth(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 xtmp,ytmp,ztmp,spfactor,delx,dely,delz,rsq
      real*8 r2inv,r6inv,forcelj,r,t,tsq,fskin,force

      do i = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = fskin*r
            endif

            force = spfactor*forcelj*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
          
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and cut Coulombic

      subroutine lj_smooth_coul_cut(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 spfactor_coul,spfactor_lj
      real*8 rsq,r2inv,forcecoul,r6inv,forcelj,r,t,tsq,fskin
      real*8 force

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich)
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else
              forcecoul = 0.0
            endif

            if (rsq < cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else if (rsq < cutljsq(itype,jtype)) then
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = fskin*r
            else
              forcelj = 0.0
            endif

            force = (spfactor_coul*forcecoul +
     $           spfactor_lj*forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
          
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and smooth Coulombic

      subroutine lj_smooth_coul_smooth(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 spfactor_coul,spfactor_lj
      real*8 rsq,r2inv,forcecoul,r,dd1,switch_r,deriv_sw_r
      real*8 r6inv,forcelj,t,tsq,fskin,force

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich)
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulintsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)
            else if (rsq < cutcoulsq) then
              r = sqrt(rsq)
              dd1 = (cutcoul - cutcoulint)**3
              switch_r = (cutcoul - r)**2 *
     $             (cutcoul + (2.0D0*r) - (3.0D0*cutcoulint)) / dd1
              deriv_sw_r = 6.0D0*(r - cutcoul)*(r - cutcoulint) / dd1
              forcecoul = (coulpre *qtmp*q(j)*sqrt(r2inv)*switch_r) -
     $             (coulpre*qtmp*q(j)*deriv_sw_r)
            else
              forcecoul = 0.0D0
            endif

            if (rsq < cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else if (rsq < cutljsq(itype,jtype)) then
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = fskin*r
            else
              forcelj = 0.0
            endif

            force = (spfactor_coul*forcecoul +
     $           spfactor_lj*forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
          
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c smooth LJ and long-range Coulombic

      subroutine lj_smooth_coul_long(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 spfactor_coul,spfactor_lj
      real*8 rsq,r2inv,r,grij,expm2,t,erfc,prefactor
      real*8 forcecoul,r6inv,forcelj,tsq,fskin,force

      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich) - 2.0
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              r = sqrt(rsq)
              grij = gewald*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)/r
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor_coul < 1.0) forcecoul = forcecoul -
     $             (1.0-spfactor_coul)*prefactor
            else
              forcecoul = 0.0
            endif

            if (spfactor_lj == 0.0) then
              forcelj = 0.0
            else if (rsq < cutljinnersq(itype,jtype)) then
              r6inv = r2inv*r2inv*r2inv
              forcelj = spfactor_lj*r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            else if (rsq < cutljsq(itype,jtype)) then
              r = sqrt(rsq) 
              t = r - cutljinner(itype,jtype)
              tsq = t*t
              fskin = ljsw1(itype,jtype) + ljsw2(itype,jtype)*t +
     $             ljsw3(itype,jtype)*tsq + ljsw4(itype,jtype)*tsq*t
              forcelj = spfactor_lj*fskin*r
            else
              forcelj = 0.0
            endif

            force = (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif

        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c shifted LJ and no Coulombic

      subroutine lj_shift(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 xtmp,ytmp,ztmp,spfactor,delx,dely,delz,rsq
      real*8 r,rshift,rshiftsq,r2inv,r6inv,forcelj,force

      do i = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r = sqrt(rsq)
            rshift = r - lj5(itype,jtype)
            rshiftsq = rshift*rshift
            r2inv = 1.0/rshiftsq
            r6inv = r2inv*r2inv*r2inv
            forcelj = r6inv*(lj1(itype,jtype)*r6inv-lj2(itype,jtype))
            force = spfactor*forcelj/rshift/r

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c soft cosine potential and no Coulombic
c set current values of ramped pre-factors before computing forces

      subroutine soft(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,j,itype,k,iwhich,jtype
      real*8 pi,xtmp,ytmp,ztmp,spfactor,delx,dely,delz
      real*8 rsq,r,arg,force

      data pi /3.1415926/

      do i = 1,ntypes
        do j = 1,ntypes
          lj4(i,j) = lj1(i,j) +
     $         float(itime)/nsteps*(lj2(i,j)-lj1(i,j))
        enddo
      enddo
      
      do i = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r = sqrt(rsq)
            arg = pi*r/lj3(itype,jtype)
            force = spfactor * lj4(itype,jtype) *
     $           pi/lj3(itype,jtype)/r * sin(arg)

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo
      
      return
      end


c -------------------------------------------------------------------------
c cut class 2 LJ and no Coulombic 

      subroutine ljclass2_cut(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 xtmp,ytmp,ztmp,spfactor,delx,dely,delz,rsq
      real*8 rinv,r2inv,r3inv,r6inv,forcelj,force

      do i = 1,nlocal
        itype = type(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            rinv = 1.0/sqrt(rsq)
            r2inv = rinv*rinv
            r3inv = r2inv*rinv
            r6inv = r3inv*r3inv
            forcelj = r6inv*(lj1(itype,jtype)*r3inv-lj2(itype,jtype))
            force = spfactor*forcelj*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut class 2 LJ and cut Coulombic

      subroutine ljclass2_cut_coul_cut(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz
      real*8 spfactor_coul,spfactor_lj
      real*8 rsq,rinv,r2inv,forcecoul,r3inv,r6inv,forcelj
      real*8 force

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich)
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            rinv = 1.0/sqrt(rsq)
            r2inv = rinv*rinv

            if (rsq < cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*rinv
            else
              forcecoul = 0.0
            endif

            if (rsq < cutljsq(itype,jtype)) then
              r3inv = r2inv*rinv
              r6inv = r3inv*r3inv
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r3inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (spfactor_coul*forcecoul + 
     $           spfactor_lj*forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c cut class 2 LJ and long-range Coulombic

      subroutine ljclass2_cut_coul_long(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz,rsq
      real*8 spfactor_coul,spfactor_lj
      real*8 r,rinv,r2inv,grij,expm2,t,erfc,prefactor,forcecoul
      real*8 r3inv,r6inv,forcelj,force

      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor_coul = 1.0
            spfactor_lj = 1.0
          else
            iwhich = j/maxatomp1
            spfactor_coul = special(iwhich) - 2.0
            spfactor_lj = spfactor_coul
            if (amberflag == 1) spfactor_lj = 0.5
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz
          jtype = type(j)

          if (rsq < cutforcesq(itype,jtype)) then

            r = sqrt(rsq)
            rinv = 1.0/r
            r2inv = rinv*rinv

            if (rsq < cutcoulsq) then
              grij = gewald*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)*rinv
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor_coul < 1.0) forcecoul = forcecoul -
     $             (1.0-spfactor_coul)*prefactor
            else
              forcecoul = 0.0
            endif

            if (rsq < cutljsq(itype,jtype) .and.
     $           spfactor_lj > 0.0) then
              r3inv = r2inv*rinv
              r6inv = r3inv*r3inv
              forcelj = spfactor_lj *
     $             r6inv*(lj1(itype,jtype)*r3inv-lj2(itype,jtype))
            else
              forcelj = 0.0
            endif

            force = (forcecoul+forcelj) * r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c CHARMM switched LJ and switched Coulombic
c roff and ron must be same for Coulomb and LJ

      subroutine lj_charmm_coul_charmm(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 qtmp,xtmp,ytmp,ztmp,spfactor,delx,dely,delz
      real*8 rsq,r2inv,forcecoul,r6inv,forcelj,force
      real*8 philj,switch,dswitch

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich)
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz

          if (rsq < cutforce_sq) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              forcecoul = coulpre*qtmp*q(j)*sqrt(r2inv)

              if (rsq > cutcoulintsq) then
                switch = (cutcoulsq - rsq)*(cutcoulsq - rsq)*
     $            (cutcoulsq + 2.*rsq - 3.*cutcoulintsq)/ch_denom_coul
                dswitch = 12.0*rsq*(cutcoulsq - rsq)*
     $                        (rsq - cutcoulintsq)/ch_denom_coul
                forcecoul = forcecoul*(switch + dswitch)
              endif

            else
              forcecoul = 0.0
            endif

            if (rsq < cutlj_sq) then
              r6inv = r2inv*r2inv*r2inv
              jtype = type(j)
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))

              if (rsq > cutljint_sq) then
                switch = (cutlj_sq - rsq)*(cutlj_sq - rsq)*                   
     $            (cutlj_sq + 2.*rsq - 3.*cutljint_sq)/ch_denom_lj
                dswitch = 12.0*rsq*(cutlj_sq - rsq)* 
     $                          (rsq - cutljint_sq)/ch_denom_lj
                philj = r6inv *
     $               (lj3(itype,jtype)*r6inv - lj4(itype,jtype))
                forcelj = forcelj * switch + philj * dswitch
              end if

            else
              forcelj = 0.0
            endif

            force = spfactor*(forcecoul + forcelj)*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c CHARMM switched LJ and long-range Coulombic

      subroutine lj_charmm_coul_long(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer i,itype,k,j,iwhich,jtype
      real*8 ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
      real*8 qtmp,xtmp,ytmp,ztmp,delx,dely,delz,spfactor
      real*8 rsq,r2inv,r,grij,expm2,t,erfc,prefactor
      real*8 forcecoul,r6inv,philj,forcelj,force
      real*8 switch,dswitch

      data ewald_p,ewald_a1,ewald_a2,ewald_a3,ewald_a4,ewald_a5
     $     /0.3275911,0.254829592,-0.284496736,
     $     1.421413741,-1.453152027,1.061405429/

      do i = 1,nlocal
        itype = type(i)
        qtmp = q(i)
        xtmp = x(1,i)
        ytmp = x(2,i)
        ztmp = x(3,i)
        do k = nnfirst(i),nnlast(i)

          j = nlist(k)

          if (j <= maxatom) then
            spfactor = 1.0
          else
            iwhich = j/maxatomp1
            spfactor = special(iwhich) - 2.0
            j = mod(j,maxatomp1) + 1
          endif

          delx = xtmp - x(1,j)
          dely = ytmp - x(2,j)
          delz = ztmp - x(3,j)
          rsq = delx*delx + dely*dely + delz*delz

          if (rsq < cutforce_sq) then

            r2inv = 1.0/rsq

            if (rsq < cutcoulsq) then
              r = sqrt(rsq)
              grij = gewald*r
              expm2 = exp(-grij*grij)
              t = 1.0/(1.0+ewald_p*grij)
              erfc = t*(ewald_a1+t*(ewald_a2+t*(ewald_a3+t*
     $             (ewald_a4+t*ewald_a5))))*expm2
              prefactor = coulpre*qtmp*q(j)/r
              forcecoul = prefactor*(erfc+1.12837917*grij*expm2)
              if (spfactor < 1.0) forcecoul = forcecoul -
     $             (1.0-spfactor)*prefactor
            else
              forcecoul = 0.0
            endif

            if (rsq < cutlj_sq) then
              r6inv = r2inv*r2inv*r2inv
              jtype = type(j)
              forcelj = r6inv *
     $             (lj1(itype,jtype)*r6inv-lj2(itype,jtype))

              if (rsq > cutljint_sq) then
                switch = (cutlj_sq - rsq)*(cutlj_sq - rsq)*                   
     $            (cutlj_sq + 2.*rsq - 3.*cutljint_sq)/ch_denom_lj
                dswitch = 12.0*rsq*(cutlj_sq - rsq)* 
     $                          (rsq - cutljint_sq)/ch_denom_lj
                philj = r6inv *
     $               (lj3(itype,jtype)*r6inv - lj4(itype,jtype))
                forcelj = forcelj * switch + philj * dswitch
              end if
            
              forcelj = forcelj*spfactor

            else
              forcelj = 0.0
            endif
            
            force = (forcecoul + forcelj)*r2inv

            f(1,i) = f(1,i) + delx*force
            f(2,i) = f(2,i) + dely*force
            f(3,i) = f(3,i) + delz*force
            if (newton_nonbond == 1 .or. j <= nlocal) then
              f(1,j) = f(1,j) - delx*force
              f(2,j) = f(2,j) - dely*force
              f(3,j) = f(3,j) - delz*force
            endif

            if (iflag == 1) then
              if (newton_nonbond == 1 .or. j <= nlocal) then
                virial(1) = virial(1) + delx*delx*force
                virial(2) = virial(2) + dely*dely*force
                virial(3) = virial(3) + delz*delz*force
                virial(4) = virial(4) + delx*dely*force
                virial(5) = virial(5) + delx*delz*force
                virial(6) = virial(6) + dely*delz*force
              else
                virial(1) = virial(1) + 0.5*delx*delx*force
                virial(2) = virial(2) + 0.5*dely*dely*force
                virial(3) = virial(3) + 0.5*delz*delz*force
                virial(4) = virial(4) + 0.5*delx*dely*force
                virial(5) = virial(5) + 0.5*delx*delz*force
                virial(6) = virial(6) + 0.5*dely*delz*force
              endif
            endif

          endif
        enddo
      enddo

      return
      end
