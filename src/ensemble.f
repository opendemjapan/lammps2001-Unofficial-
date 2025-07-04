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

c Contributing author: Mark Stevens (Sandia)

c -------------------------------------------------------------------------
c constant NVE ensemble

      subroutine nve_v(f_partial,dt_partial)
      use global
      implicit none

c argument variables

      real*8 f_partial(3,*),dt_partial

c local variables

      integer i
      real*8 dtmass

      do i = 1,nlocal
        dtmass = dt_partial/mass(type(i))
        v(1,i) = v(1,i) + dtmass*f_partial(1,i)
        v(2,i) = v(2,i) + dtmass*f_partial(2,i)
        v(3,i) = v(3,i) + dtmass*f_partial(3,i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c constant NVE ensemble

      subroutine nve_x
      use global
      implicit none

c local variables

      integer i

      do i = 1,nlocal
        x(1,i) = x(1,i) + dt*v(1,i)
        x(2,i) = x(2,i) + dt*v(2,i)
        x(3,i) = x(3,i) + dt*v(3,i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c constant NVT ensemble

      subroutine nvt_eta(dt_partial,nvt_flag)
      use global
      implicit none

c argument variables

      real*8 dt_partial
      integer nvt_flag

c local variables

      real*8 f_eta

      t_target = t_start + float(itime)/nsteps*(t_stop-t_start)
      f_eta = t_freq*t_freq * (t_current/t_target - 1.0)
      eta_dot = eta_dot + f_eta*dt_partial
      if (nvt_flag == 1) eta = eta + 2.0*dt_partial*eta_dot

      return
      end


c -------------------------------------------------------------------------
c constant NVT ensemble

      subroutine nvt_v(f_partial,dt_partial,nvt_flag)
      use global
      implicit none

c argument variables

      integer nvt_flag
      real*8 f_partial(3,*),dt_partial

c local variables

      integer i
      real*8 factor,factor1,dt_4

      factor = exp(-dt_partial*eta_dot)

      if (nvt_flag == 1) then
         factor1 = 1.0
      else if (nvt_flag == 2) then
         factor1 = factor
      endif

      do i = 1,nlocal
	dt_4 = dt_partial/(mass(type(i)))*factor1
        v(1,i) = v(1,i)*factor + dt_4*f_partial(1,i)
        v(2,i) = v(2,i)*factor + dt_4*f_partial(2,i)
        v(3,i) = v(3,i)*factor + dt_4*f_partial(3,i)
      enddo

      return
      end


c -------------------------------------------------------------------------
c constant NPH ensemble
c reuses most routines from NPT ensemble

c set t_target = t_current since npt_omegadot needs t_target

      subroutine nph_eta
      use global
      implicit none

      t_target = t_current

      return
      end


c -------------------------------------------------------------------------
c constant NPT ensemble

      subroutine npt_omegadot(dt_partial,p_omega_dt)
      use global
      implicit none

c argument variables

      real*8 dt_partial,p_omega_dt(3)

c local variables

      real*8 dtime,denskt,f_omega(3)

      dtime = float(itime)/nsteps
      denskt = (natoms*boltz*t_target) / (xprd*yprd*zprd)

c advance omega_dot in time
c for non-varying dimensions:
c   p_freq is 0.0, so f_omega will be 0.0
c   omega_dot is preset to 0.0, so p_omega_dt will also be 0.0

      p_target(1) = p_start(1) + dtime*(p_stop(1)-p_start(1))
      f_omega(1) = p_freq(1) * (p_current(1)-p_target(1))/denskt 
      omega_dot(1) = omega_dot(1) + f_omega(1)*dt_partial
      p_omega_dt(1) = dt_partial*omega_dot(1)

      p_target(2) = p_start(2) + dtime*(p_stop(2)-p_start(2))
      f_omega(2) = p_freq(2) * (p_current(2)-p_target(2))/denskt 
      omega_dot(2) = omega_dot(2) + f_omega(2)*dt_partial
      p_omega_dt(2) = dt_partial*omega_dot(2)

      p_target(3) = p_start(3) + dtime*(p_stop(3)-p_start(3))
      f_omega(3) = p_freq(3) * (p_current(3)-p_target(3))/denskt 
      omega_dot(3) = omega_dot(3) + f_omega(3)*dt_partial
      p_omega_dt(3) = dt_partial*omega_dot(3)

      return
      end


c -------------------------------------------------------------------------
c constant NPT ensemble

      subroutine npt_v(f_partial,dt_partial,npt_flag)
      use global
      implicit none

c argument variables

      integer npt_flag
      real*8 f_partial(3,*),dt_partial

c local variables

      integer i
      real*8 factor(3),factor1(3),dt_4

      factor(1) = exp(-dt_partial*(eta_dot+omega_dot(1)))
      factor(2) = exp(-dt_partial*(eta_dot+omega_dot(2)))
      factor(3) = exp(-dt_partial*(eta_dot+omega_dot(3)))

      if (npt_flag == 1) then
         factor1(1) = 1.0
         factor1(2) = 1.0
         factor1(3) = 1.0
      else if (npt_flag == 2) then
         factor1(1) = factor(1)
         factor1(2) = factor(2)
         factor1(3) = factor(3)
      endif

      do i = 1,nlocal
        dt_4 = dt_partial/(mass(type(i)))
        v(1,i) = v(1,i)*factor(1) + dt_4*f_partial(1,i)*factor1(1)
        v(2,i) = v(2,i)*factor(2) + dt_4*f_partial(2,i)*factor1(2)
        v(3,i) = v(3,i)*factor(3) + dt_4*f_partial(3,i)*factor1(3)
      enddo

      return
      end


c -------------------------------------------------------------------------
c constant NPT ensemble

      subroutine npt_x(dt_partial,p_omega_dt,npt_flag)
      use global
      use mpi
      implicit none

c argument variables

      integer npt_flag
      real*8 dt_partial,p_omega_dt(3)

c local variables

      integer i,ix,iy,iz,ierror
      real*8 dilation(3),volume,rmass
      real*8 cm(3),cmall(3)

c dilate box & particles
c must use same dilation factor for box & particles
c propagate omega even though don't use it for dynamics
c  can be used in diagnostics as constant of motion
c  volume = exp(omega(1)+omega(2)+omega(3))
c dilation factor = 1.0 for non-varying dimensions
c  since omega and p_omega_dt are 0.0

      omega(1) = omega(1) + p_omega_dt(1)
      dilation(1) = exp(p_omega_dt(1))
      omega(2) = omega(2) + p_omega_dt(2)
      dilation(2) = exp(p_omega_dt(2))
      omega(3) = omega(3) + p_omega_dt(3)
      dilation(3) = exp(p_omega_dt(3))

c center of mass using true positions

      cm(1) = 0.0
      cm(2) = 0.0
      cm(3) = 0.0
      do i = 1,nlocal
        rmass = mass(type(i))
        ix = mod(true(i),1000) - 500
        iy = mod(true(i)/1000,1000) - 500
        iz = true(i)/1000000 - 500
        cm(1) = cm(1) + rmass*(x(1,i) + ix * xprd)
        cm(2) = cm(2) + rmass*(x(2,i) + iy * yprd)
        cm(3) = cm(3) + rmass*(x(3,i) + iz * zprd)
      enddo
      call mpi_allreduce(cm,cmall,3,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)
      cm(1) = cmall(1)/masssum
      cm(2) = cmall(2)/masssum
      cm(3) = cmall(3)/masssum

c dilate the simulation box about the center of mass 
c re-define simulation box via xprd/yprd/zprd
c scale atom coords
c don't do anything if non-periodic or press style is constant volume

      if (perflagx == 0 .and. p_freq(1) /= 0.0) then
        box(1,1) = (box(1,1)-cm(1))*dilation(1) + cm(1)
        box(2,1) = (box(2,1)-cm(1))*dilation(1) + cm(1)
        xprd = box(2,1) - box(1,1)
        xprd_half = 0.5 * xprd
        do i = 1,nlocal
          x(1,i) = cm(1) + (x(1,i)-cm(1))*dilation(1)
        enddo
      endif

      if (perflagy == 0 .and. p_freq(2) /= 0.0) then
        box(1,2) = (box(1,2)-cm(2))*dilation(2) + cm(2)
        box(2,2) = (box(2,2)-cm(2))*dilation(2) + cm(2)
        yprd = box(2,2) - box(1,2)
        yprd_half = 0.5 * yprd
        do i = 1,nlocal
          x(2,i) = cm(2) + (x(2,i)-cm(2))*dilation(2)
        enddo
      endif

      if (perflagz == 0 .and. p_freq(3) /= 0.0) then
        box(1,3) = (box(1,3)-cm(3))*dilation(3) + cm(3)
        box(2,3) = (box(2,3)-cm(3))*dilation(3) + cm(3)
        zprd = box(2,3) - box(1,3)
        zprd_half = 0.5 * zprd
        do i = 1,nlocal
          x(3,i) = cm(3) + (x(3,i)-cm(3))*dilation(3)
        enddo
      endif

c redo long-range coeffs because volume has changed
c only done once per timestep, before Ewald/PPPM force calls

      if (npt_flag == 1) then
        if (coulstyle == 3) then
          call ewald_coeff
        else if (coulstyle == 4) then
          call pppm_coeff(1)
        endif
      endif

      return
      end


c -------------------------------------------------------------------------
c Langevin thermostat
c used both on full forces (for non-RESPA) and partial forces (for RESPA)
c don't apply new force to z dimension if 2d simulation

      subroutine langevin(f_partial,dt_partial)
      use global
      implicit none

c argument variables

      real*8 f_partial(3,*),dt_partial

c local variables

      integer i
      real*8 target,gconst,rmass,gamma1,gamma2,ranmars

      target = t_start + float(itime)/nsteps*(t_stop-t_start)
      gconst = sqrt(12.0*boltz*target*t_freq/dt_partial)

      if (idimension == 3) then
        do i = 1,nlocal
          rmass = mass(type(i))
          gamma1 = -t_freq*rmass
          gamma2 = gconst*sqrt(rmass)
          f_partial(1,i) = f_partial(1,i) +
     $         gamma1*v(1,i) + gamma2*(ranmars(0)-0.5)
          f_partial(2,i) = f_partial(2,i) +
     $         gamma1*v(2,i) + gamma2*(ranmars(0)-0.5)
          f_partial(3,i) = f_partial(3,i) +
     $       gamma1*v(3,i) + gamma2*(ranmars(0)-0.5)
        enddo
      else
        do i = 1,nlocal
          rmass = mass(type(i))
          gamma1 = -t_freq*rmass
          gamma2 = gconst*sqrt(rmass)
          f_partial(1,i) = f_partial(1,i) +
     $         gamma1*v(1,i) + gamma2*(ranmars(0)-0.5)
          f_partial(2,i) = f_partial(2,i) +
     $         gamma1*v(2,i) + gamma2*(ranmars(0)-0.5)
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c linear volume expand/contract in each of 3 dimensions
c reset global box and scale all atom coords

      subroutine vol_rescale
      use global

c local variables

      integer i
      real*8 old_lo,old_hi,old_box,ratio

      if (voldimx == 1) then
        old_lo = box(1,1)
        old_hi = box(2,1)
        old_box = xprd
        box(1,1) = volstart_xlo +
     $       float(itime)/nsteps*(volstop_xlo-volstart_xlo)
        box(2,1) = volstart_xhi +
     $       float(itime)/nsteps*(volstop_xhi-volstart_xhi)
        xprd = box(2,1) - box(1,1)
        xprd_half = 0.5 * xprd
        ratio = xprd/old_box

        do i = 1,nlocal
          x(1,i) = box(1,1) + (x(1,i)-old_lo) * ratio
        enddo
      endif

      if (voldimy == 1) then
        old_lo = box(1,2)
        old_hi = box(2,2)
        old_box = yprd
        box(1,2) = volstart_ylo +
     $       float(itime)/nsteps*(volstop_ylo-volstart_ylo)
        box(2,2) = volstart_yhi +
     $       float(itime)/nsteps*(volstop_yhi-volstart_yhi)
        yprd = box(2,2) - box(1,2)
        yprd_half = 0.5 * yprd
        ratio = yprd/old_box

        do i = 1,nlocal
          x(2,i) = box(1,2) + (x(2,i)-old_lo) * ratio
        enddo
      endif

      if (voldimz == 1) then
        old_lo = box(1,3)
        old_hi = box(2,3)
        old_box = zprd
        box(1,3) = volstart_zlo +
     $       float(itime)/nsteps*(volstop_zlo-volstart_zlo)
        box(2,3) = volstart_zhi +
     $       float(itime)/nsteps*(volstop_zhi-volstart_zhi)
        zprd = box(2,3) - box(1,3)
        zprd_half = 0.5 * zprd
        ratio = zprd/old_box

        do i = 1,nlocal
          x(3,i) = box(1,3) + (x(3,i)-old_lo) * ratio
        enddo
      endif

c redo long-range coeffs because volume has changed

      if (coulstyle == 3) then
        call ewald_coeff
      else if (coulstyle == 4) then
        call pppm_coeff(1)
      endif

      return
      end
