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

c Contributing authors: Roy Pollock (LLNL), Paul Crozier (Sandia)

c -------------------------------------------------------------------------
c pre-compute Ewald coefficients for potential energy and fields
c called at start of a run
c called during a run in simulation box volume changes
c called at end of run for timing purposes

      subroutine ewald_coeff
      use global
      implicit none

c local variables

      integer kmaxold,nkxmx,nkymx,nkzmx,m,k,l,kmaxall
      real*8 pi,gsqmx,gewaldsqinv,unitkx,unitky,unitkz,preu,sqk,vterm

      pi = 4.0D0*atan(1.0D0)
      gewald = (1.35D0 - 0.15D0*log(long_prec))/cutcoul

c adjustment of z dimension for 2d slab Ewald
c 3d Ewald just uses zprd since volfactor = 1.0

      zprd_slab = zprd*slab_volfactor

      nkxmx = (gewald*xprd/pi)*sqrt(-log(long_prec))
      nkymx = (gewald*yprd/pi)*sqrt(-log(long_prec))
      nkzmx = (gewald*zprd_slab/pi)*sqrt(-log(long_prec))

      kmaxold = kmax
      kmax = max(nkxmx,nkymx,nkzmx)
      kmaxall = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax

c print out size of Ewald problem

      if (kmax /= kmaxold) then
        if (node == 0) then
          write (6,*) 'Ewald G, kmax, max-vecs =',gewald,kmax,kmaxall
          write (1,*) 'Ewald G, kmax, max-vecs =',gewald,kmax,kmaxall
        endif
      endif

c allocate new Ewald memory if current size > old size
c deallocate previous memory if necessary

      if (kmax > kmaxold) then
        if (allocated(kxvecs)) then
          deallocate(kxvecs,kyvecs,kzvecs)
          deallocate(ug,eg,vg,ek)
          deallocate(sfacrl,sfacim,sfacrl_all,sfacim_all)
          deallocate(cs,sn)
        endif

        allocate (kxvecs(kmaxall),kyvecs(kmaxall),kzvecs(kmaxall))
        allocate (ug(kmaxall),eg(3,kmaxall),vg(6,kmaxall))
        allocate (ek(3,maxown))
        allocate (sfacrl(kmaxall),sfacim(kmaxall))
        allocate (sfacrl_all(kmaxall),sfacim_all(kmaxall))
        allocate (cs(maxown,3,-kmax:kmax),sn(maxown,3,-kmax:kmax))
      endif

c pre-compute coefficients for each Ewald vector

      gsqmx = -4.0D0*gewald*gewald*log(long_prec)
      gewaldsqinv = 1.0/(gewald*gewald)

      unitkx = 2.0*pi/xprd
      unitky = 2.0*pi/yprd
      unitkz = 2.0*pi/zprd_slab
      preu = 4.0*pi/(xprd*yprd*zprd_slab) 

      kcount = 0
c
c (k,0,0),(0,l,0),(0,0,m)
c
      do m = 1,kmax
        sqk = (m*unitkx)**2
        if (sqk.le.gsqmx) then
          kcount = kcount + 1
          kxvecs(kcount) = m
          kyvecs(kcount) = 0
          kzvecs(kcount) = 0
          ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
          eg(1,kcount) = 2.0*unitkx*m*ug(kcount)
          eg(2,kcount) = 0.0
          eg(3,kcount) = 0.0
          vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
          vg(1,kcount) = 1.0 + vterm*(unitkx*m)**2
          vg(2,kcount) = 1.0
          vg(3,kcount) = 1.0
          vg(4,kcount) = 0.0
          vg(5,kcount) = 0.0
          vg(6,kcount) = 0.0
        endif
        sqk = (m*unitky)**2
        if (sqk.le.gsqmx) then
          kcount = kcount + 1
          kxvecs(kcount) = 0
          kyvecs(kcount) = m
          kzvecs(kcount) = 0
          ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
          eg(1,kcount) = 0.0
          eg(2,kcount) = 2.0*unitky*m*ug(kcount)
          eg(3,kcount) = 0.0
          vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
          vg(1,kcount) = 1.0
          vg(2,kcount) = 1.0 + vterm*(unitky*m)**2
          vg(3,kcount) = 1.0
          vg(4,kcount) = 0.0
          vg(5,kcount) = 0.0
          vg(6,kcount) = 0.0
        endif
        sqk = (m*unitkz)**2
        if (sqk.le.gsqmx) then
          kcount = kcount + 1
          kxvecs(kcount) = 0
          kyvecs(kcount) = 0
          kzvecs(kcount) = m
          ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
          eg(1,kcount) = 0.0
          eg(2,kcount) = 0.0
          eg(3,kcount) = 2.0*unitkz*m*ug(kcount)
          vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
          vg(1,kcount) = 1.0
          vg(2,kcount) = 1.0
          vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
          vg(4,kcount) = 0.0
          vg(5,kcount) = 0.0
          vg(6,kcount) = 0.0
        endif
      enddo
c
c 1 = (k,l,0), 2 = (k,-l,0)
c
      do k = 1,kmax
        do l = 1,kmax      
          sqk = (unitkx*k)**2 + (unitky*l)**2
          if (sqk.le.gsqmx) then
            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = l
            kzvecs(kcount) = 0
            ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
            eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
            eg(2,kcount) = 2.0*unitky*l*ug(kcount)
            eg(3,kcount) = 0.0
            vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
            vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
            vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
            vg(3,kcount) = 1.0
            vg(4,kcount) = vterm*unitkx*k*unitky*l
            vg(5,kcount) = 0.0
            vg(6,kcount) = 0.0

            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = -l
            kzvecs(kcount) = 0
            ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
            eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
            eg(2,kcount) = -2.0*unitky*l*ug(kcount)
            eg(3,kcount) = 0.0
            vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
            vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
            vg(3,kcount) = 1.0
            vg(4,kcount) = -vterm*unitkx*k*unitky*l
            vg(5,kcount) = 0.0
            vg(6,kcount) = 0.0
          endif
        enddo
      enddo
c
c 1 = (0,l,m), 2 = (0,l,-m)
c
      do l = 1,kmax
        do m = 1,kmax      
          sqk = (unitky*l)**2+(unitkz*m)**2
          if (sqk.le.gsqmx) then
            kcount = kcount + 1
            kxvecs(kcount) = 0
            kyvecs(kcount) = l
            kzvecs(kcount) = m
            ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
            eg(1,kcount) =  0.0
            eg(2,kcount) =  2.0*unitky*l*ug(kcount)
            eg(3,kcount) =  2.0*unitkz*m*ug(kcount)
            vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
            vg(1,kcount) = 1.0
            vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
            vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
            vg(4,kcount) = 0.0
            vg(5,kcount) = 0.0
            vg(6,kcount) = vterm*unitky*l*unitkz*m

            kcount = kcount + 1
            kxvecs(kcount) = 0
            kyvecs(kcount) = l
            kzvecs(kcount) = -m
            ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
            eg(1,kcount) =  0.0
            eg(2,kcount) =  2.0*unitky*l*ug(kcount)
            eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
            vg(1,kcount) = 1.0
            vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
            vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
            vg(4,kcount) = 0.0
            vg(5,kcount) = 0.0
            vg(6,kcount) = -vterm*unitky*l*unitkz*m
          endif
        enddo
      enddo
c
c 1 = (k,0,m), 2 = (k,0,-m)
c
      do k = 1,kmax
        do m = 1,kmax      
          sqk = (unitkx*k)**2 + (unitkz*m)**2
          if (sqk.le.gsqmx) then
            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = 0
            kzvecs(kcount) = m
            ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
            eg(1,kcount) =  2.0*unitkx*k*ug(kcount)
            eg(2,kcount) =  0.0
            eg(3,kcount) =  2.0*unitkz*m*ug(kcount)
            vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
            vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
            vg(2,kcount) = 1.0
            vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
            vg(4,kcount) = 0.0
            vg(5,kcount) = vterm*unitkx*k*unitkz*m
            vg(6,kcount) = 0.0

            kcount = kcount + 1
            kxvecs(kcount) = k
            kyvecs(kcount) = 0
            kzvecs(kcount) = -m
            ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
            eg(1,kcount) =  2.0*unitkx*k*ug(kcount)
            eg(2,kcount) =  0.0
            eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
            vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
            vg(2,kcount) = 1.0
            vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
            vg(4,kcount) = 0.0
            vg(5,kcount) = -vterm*unitkx*k*unitkz*m
            vg(6,kcount) = 0.0
          endif
        enddo
      enddo
c
c 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
c
      do k = 1,kmax
        do l = 1,kmax
          do m = 1,kmax      
            sqk = (unitkx*k)**2 + (unitky*l)**2 + (unitkz*m)**2
            if (sqk.le.gsqmx) then
              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = l
              kzvecs(kcount) = m
              ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = 2.0*unitky*l*ug(kcount)
              eg(3,kcount) = 2.0*unitkz*m*ug(kcount)
              vterm = -2.0*(1.0/sqk + 0.25*gewaldsqinv)
              vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
              vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
              vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
              vg(4,kcount) = vterm*unitkx*k*unitky*l
              vg(5,kcount) = vterm*unitkx*k*unitkz*m
              vg(6,kcount) = vterm*unitky*l*unitkz*m

              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = -l
              kzvecs(kcount) = m
              ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = -2.0*unitky*l*ug(kcount)
              eg(3,kcount) = 2.0*unitkz*m*ug(kcount)
              vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
              vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
              vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
              vg(4,kcount) = -vterm*unitkx*k*unitky*l
              vg(5,kcount) = vterm*unitkx*k*unitkz*m
              vg(6,kcount) = -vterm*unitky*l*unitkz*m

              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = l
              kzvecs(kcount) = -m
              ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = 2.0*unitky*l*ug(kcount)
              eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
              vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
              vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
              vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
              vg(4,kcount) = vterm*unitkx*k*unitky*l
              vg(5,kcount) = -vterm*unitkx*k*unitkz*m
              vg(6,kcount) = -vterm*unitky*l*unitkz*m

              kcount = kcount + 1
              kxvecs(kcount) = k
              kyvecs(kcount) = -l
              kzvecs(kcount) = -m
              ug(kcount) = preu*exp(-0.25*sqk*gewaldsqinv)/sqk
              eg(1,kcount) = 2.0*unitkx*k*ug(kcount)
              eg(2,kcount) = -2.0*unitky*l*ug(kcount)
              eg(3,kcount) = -2.0*unitkz*m*ug(kcount)
              vg(1,kcount) = 1.0 + vterm*(unitkx*k)**2
              vg(2,kcount) = 1.0 + vterm*(unitky*l)**2
              vg(3,kcount) = 1.0 + vterm*(unitkz*m)**2
              vg(4,kcount) = -vterm*unitkx*k*unitky*l
              vg(5,kcount) = -vterm*unitkx*k*unitkz*m
              vg(6,kcount) = vterm*unitky*l*unitkz*m
            endif
          enddo
        enddo
      enddo

c print out actual count of Ewald vectors

      if (kmax /= kmaxold) then
        if (node == 0) then
          write (6,*) 'Actual Ewald vecs =',kcount
          write (1,*) 'Actual Ewald vecs =',kcount
        endif
      endif

      return
      end


c -------------------------------------------------------------------------
c deallocate all Ewald memory, called at end of a run

      subroutine ewald_deallocate
      use global
      implicit none

      if (allocated(kxvecs)) then
        deallocate(kxvecs,kyvecs,kzvecs)
        deallocate(ug,eg,vg,ek)
        deallocate(sfacrl,sfacim,sfacrl_all,sfacim_all)
        deallocate(cs,sn)
      endif

c this will force reallocation of memory on subsequent run

      kmax = 0

      return
      end
