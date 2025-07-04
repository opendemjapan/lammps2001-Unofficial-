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
c compute k-space part of Ewald forces and potential

      subroutine ewald(iflag)
      use global
      use mpi
      implicit none

c argument variables

      integer iflag

c local variables

      integer ierror,i,kk,kx,ky,kz
      real*8 pi,cypz,sypz,exprl,expim,partial,ukk

c compute partial structure factors by summing over particles on processor

      call eikdotr

c compute total structure factor by summing over all processors

      call mpi_allreduce(sfacrl,sfacrl_all,kcount,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(sfacim,sfacim_all,kcount,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

c compute k-space part of electrostatic potential energy

      e_long = 0.0D0

      if (iflag == 0) then

        do kk = 1,kcount
          ukk = ug(kk)*(sfacrl_all(kk)**2+sfacim_all(kk)**2)
          e_long = e_long + ukk
        enddo

      else

        do i = 1,6
          vir_long(i) = 0.0D0
        enddo
        do kk = 1,kcount
          ukk = ug(kk)*(sfacrl_all(kk)**2+sfacim_all(kk)**2)
          e_long = e_long + ukk
          do i = 1,6
            vir_long(i) = vir_long(i) + ukk*vg(i,kk)
          enddo
        enddo
        do i = 1,6
          vir_long(i) = coulpre*vir_long(i)
        enddo

      endif

      pi = 4.0D0*atan(1.0D0)
      e_long = e_long - gewald*qsqsum/1.772453851D0 -
     $     0.5D0*pi*qsum*qsum/(gewald*gewald*xprd*yprd*zprd_slab)
      e_long = coulpre*e_long

c compute k-space part of electric field for particles on processor

      do i = 1,nlocal
        ek(1,i) = 0.0D0
        ek(2,i) = 0.0D0
        ek(3,i) = 0.0D0
      enddo

      do kk = 1,kcount
        kx = kxvecs(kk)
        ky = kyvecs(kk)
        kz = kzvecs(kk)

        do i = 1,nlocal
          cypz = cs(i,2,ky)*cs(i,3,kz) - sn(i,2,ky)*sn(i,3,kz)
          sypz = sn(i,2,ky)*cs(i,3,kz) + cs(i,2,ky)*sn(i,3,kz)
          exprl = cs(i,1,kx)*cypz - sn(i,1,kx)*sypz
          expim = sn(i,1,kx)*cypz + cs(i,1,kx)*sypz
          partial = expim*sfacrl_all(kk) - exprl*sfacim_all(kk)
          ek(1,i) = ek(1,i) + eg(1,kk)*partial
          ek(2,i) = ek(2,i) + eg(2,kk)*partial
          ek(3,i) = ek(3,i) + eg(3,kk)*partial
        enddo
      enddo

c convert e-field to force

      do i = 1,nlocal
        f(1,i) = f(1,i) + q(i)*coulpre*ek(1,i)
        f(2,i) = f(2,i) + q(i)*coulpre*ek(2,i)
        f(3,i) = f(3,i) + q(i)*coulpre*ek(3,i)
      enddo

c correct for slab geometry if needed

      if (slabflag == 1) call slabcorr

      return
      end


c -------------------------------------------------------------------------
c calculate structure factors and exp(k*rr(i)) 

      subroutine eikdotr
      use global

c local variables

      integer ic,i,m,k,l
      real*8 pi,gsqmx,sqk,cstr1,sstr1,cstr2,sstr2
      real*8 cstr3,sstr3,cstr4,sstr4,clpm,slpm
      real*8 unitk(3)

      pi = 4.0*atan(1.0)
      unitk(1) = 2.0*pi/xprd
      unitk(2) = 2.0*pi/yprd
      unitk(3) = 2.0*pi/zprd_slab
      gsqmx = -4.0*gewald*gewald*log(long_prec)

      kcount = 0

c (k,0,0), (0,l,0), (0,0,m)

      do ic = 1,3
        sqk = unitk(ic)**2
        if (sqk.le.gsqmx) then
          cstr1 = 0.0
          sstr1 = 0.0
          do i = 1,nlocal
            cs(i,ic,0) = 1.0
            sn(i,ic,0) = 0.0
            cs(i,ic,1) = cos(unitk(ic)*x(ic,i))
            sn(i,ic,1) = sin(unitk(ic)*x(ic,i))
            cs(i,ic,-1) = cs(i,ic,1)
            sn(i,ic,-1) = -sn(i,ic,1)
            cstr1 = cstr1 + q(i)*cs(i,ic,1)
            sstr1 = sstr1 + q(i)*sn(i,ic,1)
          enddo
          kcount = kcount + 1
          sfacrl(kcount) = cstr1
          sfacim(kcount) = sstr1
        endif
      enddo

      do m = 2,kmax
        do ic = 1,3
          sqk = (m*unitk(ic))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            do i = 1,nlocal
              cs(i,ic,m) = cs(i,ic,m-1)*cs(i,ic,1) -
     $             sn(i,ic,m-1)*sn(i,ic,1)
              sn(i,ic,m) = sn(i,ic,m-1)*cs(i,ic,1) +
     $             cs(i,ic,m-1)*sn(i,ic,1)
              cs(i,ic,-m) = cs(i,ic,m)
              sn(i,ic,-m) = -sn(i,ic,m)
              cstr1 = cstr1 + q(i)*cs(i,ic,m)
              sstr1 = sstr1 + q(i)*sn(i,ic,m)
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1
          endif
        enddo
      enddo

c 1 = (k,l,0), 2 = (k,-l,0)

      do k = 1,kmax
        do l = 1,kmax
          sqk = (k*unitk(1))**2 + (l*unitk(2))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            cstr2 = 0.0
            sstr2 = 0.0
            do i = 1,nlocal
              cstr1 = cstr1 + q(i)*(cs(i,1,k)*cs(i,2,l) -
     $             sn(i,1,k)*sn(i,2,l))
              sstr1 = sstr1 + q(i)*(sn(i,1,k)*cs(i,2,l) +
     $             cs(i,1,k)*sn(i,2,l))
              cstr2 = cstr2 + q(i)*(cs(i,1,k)*cs(i,2,l) +
     $             sn(i,1,k)*sn(i,2,l))
              sstr2 = sstr2 + q(i)*(sn(i,1,k)*cs(i,2,l) -
     $             cs(i,1,k)*sn(i,2,l))
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1

            kcount = kcount + 1
            sfacrl(kcount) = cstr2
            sfacim(kcount) = sstr2
          endif
        enddo
      enddo

c 1 = (0,l,m), 2 = (0,l,-m)

      do l = 1,kmax
        do m = 1,kmax
          sqk = (l*unitk(2))**2 + (m*unitk(3))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            cstr2 = 0.0
            sstr2 = 0.0
            do i = 1,nlocal
              cstr1 = cstr1 + q(i)*(cs(i,2,l)*cs(i,3,m) -
     $             sn(i,2,l)*sn(i,3,m))
              sstr1 = sstr1 + q(i)*(sn(i,2,l)*cs(i,3,m) +
     $             cs(i,2,l)*sn(i,3,m))
              cstr2 = cstr2 + q(i)*(cs(i,2,l)*cs(i,3,m) +
     $             sn(i,2,l)*sn(i,3,m))
              sstr2 = sstr2 + q(i)*(sn(i,2,l)*cs(i,3,m) -
     $             cs(i,2,l)*sn(i,3,m))
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1

            kcount = kcount + 1
            sfacrl(kcount) = cstr2
            sfacim(kcount) = sstr2
          endif
        enddo
      enddo

c 1 = (k,0,m), 2 = (k,0,-m)

      do k = 1,kmax
        do m = 1,kmax
          sqk = (k*unitk(1))**2 + (m*unitk(3))**2
          if (sqk.le.gsqmx) then
            cstr1 = 0.0
            sstr1 = 0.0
            cstr2 = 0.0
            sstr2 = 0.0
            do i = 1,nlocal
              cstr1 = cstr1 + q(i)*(cs(i,1,k)*cs(i,3,m) -
     $             sn(i,1,k)*sn(i,3,m))
              sstr1 = sstr1 + q(i)*(sn(i,1,k)*cs(i,3,m) +
     $             cs(i,1,k)*sn(i,3,m))
              cstr2 = cstr2 + q(i)*(cs(i,1,k)*cs(i,3,m) +
     $             sn(i,1,k)*sn(i,3,m))
              sstr2 = sstr2 + q(i)*(sn(i,1,k)*cs(i,3,m) -
     $             cs(i,1,k)*sn(i,3,m))
            enddo
            kcount = kcount + 1
            sfacrl(kcount) = cstr1
            sfacim(kcount) = sstr1

            kcount = kcount + 1
            sfacrl(kcount) = cstr2
            sfacim(kcount) = sstr2
          endif
        enddo
      enddo

c 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

      do k = 1,kmax
        do l = 1,kmax
          do m = 1,kmax
            sqk = (k*unitk(1))**2 + (l*unitk(2))**2 + (m*unitk(3))**2
            if (sqk.le.gsqmx) then
              cstr1 = 0.0
              sstr1 = 0.0
              cstr2 = 0.0
              sstr2 = 0.0
              cstr3 = 0.0
              sstr3 = 0.0
              cstr4 = 0.0
              sstr4 = 0.0
              do i = 1,nlocal
                clpm = cs(i,2,l)*cs(i,3,m) - sn(i,2,l)*sn(i,3,m)
                slpm = sn(i,2,l)*cs(i,3,m) + cs(i,2,l)*sn(i,3,m)
                cstr1 = cstr1 + q(i)*(cs(i,1,k)*clpm -
     $               sn(i,1,k)*slpm)
                sstr1 = sstr1 + q(i)*(sn(i,1,k)*clpm +
     $               cs(i,1,k)*slpm)

                clpm = cs(i,2,l)*cs(i,3,m) + sn(i,2,l)*sn(i,3,m)
                slpm = -sn(i,2,l)*cs(i,3,m) + cs(i,2,l)*sn(i,3,m)
                cstr2 = cstr2 + q(i)*(cs(i,1,k)*clpm -
     $               sn(i,1,k)*slpm)
                sstr2 = sstr2 + q(i)*(sn(i,1,k)*clpm +
     $               cs(i,1,k)*slpm)

                clpm = cs(i,2,l)*cs(i,3,m) + sn(i,2,l)*sn(i,3,m)
                slpm = sn(i,2,l)*cs(i,3,m) - cs(i,2,l)*sn(i,3,m)
                cstr3 = cstr3 + q(i)*(cs(i,1,k)*clpm -
     $               sn(i,1,k)*slpm)
                sstr3 = sstr3 + q(i)*(sn(i,1,k)*clpm +
     $               cs(i,1,k)*slpm)

                clpm = cs(i,2,l)*cs(i,3,m) - sn(i,2,l)*sn(i,3,m)
                slpm = -sn(i,2,l)*cs(i,3,m) - cs(i,2,l)*sn(i,3,m)
                cstr4 = cstr4 + q(i)*(cs(i,1,k)*clpm -
     $               sn(i,1,k)*slpm)
                sstr4 = sstr4 + q(i)*(sn(i,1,k)*clpm +
     $               cs(i,1,k)*slpm)
              enddo
              kcount = kcount + 1
              sfacrl(kcount) = cstr1
              sfacim(kcount) = sstr1

              kcount = kcount + 1
              sfacrl(kcount) = cstr2
              sfacim(kcount) = sstr2

              kcount = kcount + 1
              sfacrl(kcount) = cstr3
              sfacim(kcount) = sstr3

              kcount = kcount + 1
              sfacrl(kcount) = cstr4
              sfacim(kcount) = sstr4
            endif
          enddo
        enddo
      enddo

      return
      end


c -------------------------------------------------------------------------
c Slab-geometry correction term to dampen inter-slab interactions between
c  periodically repeating slabs. Yields good approximation to 2-D Ewald if 
c  adequate empty space is left between repeating slabs (J. Chem. Phys. 
c  111, 3155). Slabs defined here to be parallel to the 1-2 (or x-y) plane.

      subroutine slabcorr
      use global
      use mpi
      implicit none

c local variables

      integer i,ierror
      real*8 dipole,dipole_all,vol,ffact,e_slabcorr

c compute local contribution to global dipole moment

      dipole = 0.0
      do i = 1,nlocal
         dipole = dipole + q(i)*x(3,i)
      enddo    

c sum local contributions to get global dipole moment

      call mpi_allreduce(dipole,dipole_all,1,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

c compute corrections
  
      vol = xprd*yprd*zprd_slab
      ffact = -4*3.14159265359*dipole_all/vol   
      e_slabcorr = 2*3.14159265359*dipole_all*dipole_all/vol
      e_long = e_long + coulpre*e_slabcorr

c add on force corrections

      do i = 1,nlocal
        f(3,i) = f(3,i) + q(i)*coulpre*ffact
      enddo
 
      return
      end
