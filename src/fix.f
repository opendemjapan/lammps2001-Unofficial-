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
c assign group of atoms to a constraint
c   fixgroup = 1 -> single atom
c   fixgroup = 2 -> molecule
c   fixgroup = 3 -> atom type
c   fixgroup = 4 -> spatial region
c   fixgroup = 5 -> bondtype to SHAKE
c   fixgroup = 6 -> angletype to SHAKE
c   fixgroup = 7 -> all remaining unconstrained atoms

      subroutine assign_fix
      use global
      implicit none

c local variables

      integer i

      if (node == 0) write (6,*) 'Assigning atoms/bonds to a fix ...'

      if (fixgroup == 1) then
        i = localptr(fixatom)
        if (i /= 0 .and. i <= nlocal)
     $       fix(i) = maxfix*fix(i) + fixwhich
      else if (fixgroup == 2) then
        do i = 1,nlocal
          if (molecule(i) == fixatom)
     $         fix(i) = maxfix*fix(i) + fixwhich
        enddo
      else if (fixgroup == 3) then
        do i = 1,nlocal
          if (type(i) == fixtype) fix(i) = maxfix*fix(i) + fixwhich
        enddo
      else if (fixgroup == 4) then
        do i = 1,nlocal
          if (x(1,i) >= fixregion(1).and.x(1,i) <= fixregion(2) .and.
     $         x(2,i) >= fixregion(3).and.x(2,i) <= fixregion(4) .and.
     $         x(3,i) >= fixregion(5).and.x(3,i) <= fixregion(6))
     $         fix(i) = maxfix*fix(i) + fixwhich
        enddo
      else if (fixgroup == 5) then
        shakeablebond(fixbond) = 1
      else if (fixgroup == 6) then
        shakeableangle = fixangle
        shakeableanglebond = fixbond
      else if (fixgroup == 7) then
        do i = 1,nlocal
          if (fix(i) == 0) fix(i) = fixwhich
        enddo
      endif

      return
      end


c -------------------------------------------------------------------------
c figure out how many global quantities each fix will need to accumulate
c called from start routine before a run
c fixptr(k) = 1st loc in fixstore where values associated with fix k
c             will be accumulated
c fixnum = total # of values to be accumulated
c nfixes_respa = # of constraints that will be active at all phases of RESPA,
c   only setforce (1) and aveforce (3) styles

      subroutine fix_setup
      use global
      use mpi
      implicit none

c local variables

      integer i,k,iwhich,istyle,itmp,ierror

      fixnum = 0
      do k = 1,nfixes
        fixptr(k) = fixnum + 1
        istyle = fixstyle(k)
        if (istyle <= 3) then
          fixnum = fixnum + 3
        else if (istyle == 4) then
          fixnum = fixnum + 1
        else if (istyle == 5) then
          fixnum = fixnum + 1
        else if (istyle == 7) then
          fixnum = fixnum + 3
        endif
      enddo

      do k = 1,nfixes
        fixcount(k) = 0
        fixmass(k) = 0.0
      enddo

      do i = 1,nlocal
        itmp = fix(i)
        do while (itmp > 0)
          iwhich = mod(itmp,maxfix)
          fixcount(iwhich) = fixcount(iwhich) + 1
          fixmass(iwhich) = fixmass(iwhich) + mass(type(i))
          itmp = itmp/maxfix
        enddo
      enddo

      call mpi_allreduce(fixcount,fixcounttmp,nfixes,
     $     mpi_integer,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(fixmass,fixmasstmp,nfixes,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
      do i = 1,nfixes
        fixcount(i) = fixcounttmp(i)
        fixmass(i) = fixmasstmp(i)
      enddo

      nfixes_respa = 0
      do k = 1,nfixes
        if (fixstyle(k) == 1 .or. fixstyle(k) == 3)
     $       nfixes_respa = nfixes_respa + 1
      enddo

      return
      end


c -------------------------------------------------------------------------
c fix forces and velocities for each constrained atom at end of timestep
c   applies to all fix styles
c iflag helps determines which fixes are active (see comment below)
c SHAKE constraints are not applied here, but directly from
c   start/integrate routines
c dt_fix is 1/2 a timestep whether called from integrate or start
c   this is so adjusted force at beginning of timestep will be same
c     whether on first timestep (from restart file) or successive timesteps
c   should enable exact restarts

      subroutine fix_apply(iflag,f_fix,dt_fix)
      use global
      use mpi
      implicit none

c argument variables

      integer iflag
      real*8 f_fix(3,*),dt_fix

c local variables

      integer n,i,itmp,iwhich,istyle,iboxx,iboxy,iboxz
      integer ierror,k,ncount,itrue
      real*8 rmass,target,dx,dy,dz,gamma1,gamma2,ranmars
      real*8 massfrac,r,pre
      real*8 factor,fraction,dtmass,vx,vy,vz

c decide which fixes are active:
c thermostatting fixes (4-6) are not active when called from minimizer
c   (itime = 0)
c thermostatting fixes (4-6) are not active when called from start routine
c   (itime = 0) else will be rescaled by 1st half of Verlet update
c   before 1st timestep
c rescale (4) is not active except every N timesteps
c if iflag = 0, being called from RESPA when only setforce (1)
c   and aveforce (3) are active

      do k = 1,nfixes
        fixactive(k) = 1
        istyle = fixstyle(k)
        if (istyle >= 4 .and. istyle <= 6 .and. itime == 0)
     $       fixactive(k) = 0
        if (istyle == 4) then
          if (mod(itime,nint(fixcoeff(3,k))) /= 0) fixactive(k) = 0
        endif
        if (iflag == 0) then
          fixactive(k) = 0
          if (istyle == 1 .or. istyle == 3) fixactive(k) = 1
        endif
      enddo

c zero out storage before accumulating info about each fix

      do n = 1,fixnum
        fixstore(n) = 0.0
      enddo

c accumulate each processor's contribution to quantities
c   associated with each fix
c for rescale and hoover/drag (istyles 4,5),
c   compute predicted velocity/temperature at end of this timestep by
c   adding in Verlet update (dtmass*f)
c for springforce (istyle 7),
c   compute group's center-of-mass in true coords
c   else will get bogus answer if group straddles a periodic boundary

      do i = 1,nlocal
        itmp = fix(i)
        do while (itmp > 0)

          iwhich = mod(itmp,maxfix)
          if (fixactive(iwhich) == 0) then
            itmp = itmp/maxfix
            cycle
          endif

          n = fixptr(iwhich)
          istyle = fixstyle(iwhich)

          if (istyle <= 3) then
            fixstore(n) = fixstore(n) + f_fix(1,i)
            fixstore(n+1) = fixstore(n+1) + f_fix(2,i)
            fixstore(n+2) = fixstore(n+2) + f_fix(3,i)
          else if (istyle == 4) then
            dtmass = dt_fix/mass(type(i))
            vx = v(1,i) + dtmass*f_fix(1,i)
            vy = v(2,i) + dtmass*f_fix(2,i)
            vz = v(3,i) + dtmass*f_fix(3,i)
c            vx = v(1,i)
c            vy = v(2,i)
c            vz = v(3,i)
            fixstore(n) = fixstore(n) +
     $           (vx*vx + vy*vy + vz*vz) * mass(type(i))
          else if (istyle == 5) then
            dtmass = dt_fix/mass(type(i))
            vx = v(1,i) + dtmass*f_fix(1,i)
            vy = v(2,i) + dtmass*f_fix(2,i)
            vz = v(3,i) + dtmass*f_fix(3,i)
            fixstore(n) = fixstore(n) +
     $           (vx*vx + vy*vy + vz*vz) * mass(type(i))
          else if (istyle == 7) then
            iboxx = mod(true(i),1000) - 500
            iboxy = mod(true(i)/1000,1000) - 500
            iboxz = true(i)/1000000 - 500
            rmass = mass(type(i))
            fixstore(n) = fixstore(n) + (x(1,i)+iboxx*xprd)*rmass
            fixstore(n+1) = fixstore(n+1) + (x(2,i)+iboxy*yprd)*rmass
            fixstore(n+2) = fixstore(n+2) + (x(3,i)+iboxz*zprd)*rmass
          endif

          itmp = itmp/maxfix

        enddo
      enddo

c sum fix quantities across all processors

      if (fixnum > 0) then
        call mpi_allreduce(fixstore,fixstoretmp,fixnum,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        do i = 1,fixnum
          fixstore(i) = fixstoretmp(i)
        enddo
      endif

c normalize global quantities and setup fix coeffs
c for hoover/drag (istyle 5), the drag (fixcoeff(4,k)) was initialized
c   in input.f when fix command was specified, accumulates thru all runs
c for springforce (istyle 7), compute distance from center-of-mass to
c   spring center in true coordinates

      do k = 1,nfixes

        if (fixactive(k) == 0) cycle
        n = fixptr(k)
        istyle = fixstyle(k)

        if (istyle <= 2) then
          fixstore(n) = fixstore(n)/fixcount(k)
          fixstore(n+1) = fixstore(n+1)/fixcount(k)
          fixstore(n+2) = fixstore(n+2)/fixcount(k)
        else if (istyle == 3) then
          fixstore(n) = fixstore(n)/fixcount(k)
          fixstore(n+1) = fixstore(n+1)/fixcount(k)
          fixstore(n+2) = fixstore(n+2)/fixcount(k)
          fixcoeff(4,k) = fixstore(n) + fixcoeff(1,k)
          fixcoeff(5,k) = fixstore(n+1) + fixcoeff(2,k)
          fixcoeff(6,k) = fixstore(n+2) + fixcoeff(3,k)
        else if (istyle == 4) then
          fixstore(n) = fixstore(n) / (boltz*idimension*fixcount(k))
          target = fixcoeff(1,k) +
     $         float(itime)/nsteps*(fixcoeff(2,k)-fixcoeff(1,k))
          if (abs(fixstore(n)-target) > fixcoeff(4,k) .and.
     $         fixstore(n) > 0.0) then
            target = fixstore(n) -
     $           fixcoeff(5,k)*(fixstore(n)-target)
            fixcoeff(6,k) = sqrt(target/fixstore(n))
          else
            fixactive(k) = 0
          endif
        else if (istyle == 5) then
          ncount = nint(fixstore(n))
          fixstore(n) = fixstore(n) / (boltz*idimension*fixcount(k))
          target = fixcoeff(1,k) +
     $         float(itime)/nsteps*(fixcoeff(2,k)-fixcoeff(1,k))
          fixcoeff(4,k) = fixcoeff(4,k) +
     $         (fixstore(n)-target)/target * fixcoeff(3,k)
        else if (istyle == 6) then
          target = fixcoeff(1,k) +
     $         float(itime)/nsteps*(fixcoeff(2,k)-fixcoeff(1,k))
          fixcoeff(4,k) = sqrt(12.0*boltz*target*fixcoeff(3,k)/dt_fix)
        else if (istyle == 7) then
          fixstore(n) = fixstore(n)/fixmass(k)
          fixstore(n+1) = fixstore(n+1)/fixmass(k)
          fixstore(n+2) = fixstore(n+2)/fixmass(k)
          dx = fixstore(n) - fixcoeff(1,k)
          dy = fixstore(n+1) - fixcoeff(2,k)
          dz = fixstore(n+2) - fixcoeff(3,k)
          if (fixflag(1,k) == 1) dx = 0.0
          if (fixflag(2,k) == 1) dy = 0.0
          if (fixflag(3,k) == 1) dz = 0.0
          fixcoeff(5,k) = dx*fixcoeff(4,k)
          fixcoeff(6,k) = dy*fixcoeff(4,k)
          fixcoeff(7,k) = dz*fixcoeff(4,k)
        endif

      enddo

c apply constraints to each atom by adjusting force
c 2d and 3d are separate so don't alter v,f of 3rd dimension in 2d case
c loop over all fixes on an atom in reverse of assignment order
c for rescale (istyle 4), are actually adjusting velocities, not force
c for dragforce (istyle 8), apply minimg convention to get direction of force

      if (idimension == 3) then

        do i = 1,nlocal
          itmp = fix(i)
          do while (itmp > 0)

            iwhich = mod(itmp,maxfix)
            if (fixactive(iwhich) == 0) then
              itmp = itmp/maxfix
              cycle
            endif
            istyle = fixstyle(iwhich)

            if (istyle == 1) then
              if (fixflag(1,iwhich) == 0)
     $             f_fix(1,i) = fixcoeff(1,iwhich)
              if (fixflag(2,iwhich) == 0)
     $             f_fix(2,i) = fixcoeff(2,iwhich)
              if (fixflag(3,iwhich) == 0)
     $             f_fix(3,i) = fixcoeff(3,iwhich)
            else if (istyle == 2) then
              f_fix(1,i) = f_fix(1,i) + fixcoeff(1,iwhich)
              f_fix(2,i) = f_fix(2,i) + fixcoeff(2,iwhich)
              f_fix(3,i) = f_fix(3,i) + fixcoeff(3,iwhich)
            else if (istyle == 3) then
              f_fix(1,i) = fixcoeff(4,iwhich)
              f_fix(2,i) = fixcoeff(5,iwhich)
              f_fix(3,i) = fixcoeff(6,iwhich)
            else if (istyle == 4) then
              factor = fixcoeff(6,iwhich)
              v(1,i) = v(1,i)*factor
              v(2,i) = v(2,i)*factor
              v(3,i) = v(3,i)*factor
            else if (istyle == 5) then
              factor = fixcoeff(4,iwhich)
              dtmass = mass(type(i))/dt_fix
              f_fix(1,i) = f_fix(1,i) - dtmass*factor*v(1,i)
              f_fix(2,i) = f_fix(2,i) - dtmass*factor*v(2,i)
              f_fix(3,i) = f_fix(3,i) - dtmass*factor*v(3,i)
            else if (istyle == 6) then
              rmass = mass(type(i))
              gamma1 = -fixcoeff(3,iwhich)*rmass
              gamma2 = fixcoeff(4,iwhich)*sqrt(rmass)
              if (fixflag(1,iwhich) == 1)
     $             f_fix(1,i) = f_fix(1,i) + gamma1*v(1,i) +
     $             gamma2*(ranmars(0)-0.5)
              if (fixflag(2,iwhich) == 1)
     $             f_fix(2,i) = f_fix(2,i) + gamma1*v(2,i) +
     $             gamma2*(ranmars(0)-0.5)
              if (fixflag(3,iwhich) == 1)
     $             f_fix(3,i) = f_fix(3,i) + gamma1*v(3,i) +
     $             gamma2*(ranmars(0)-0.5)
            else if (istyle == 7) then
              massfrac = mass(type(i))/fixmass(iwhich)
              if (fixflag(1,iwhich) == 0) f_fix(1,i) = f_fix(1,i) -
     $             massfrac*fixcoeff(5,iwhich)
              if (fixflag(2,iwhich) == 0) f_fix(2,i) = f_fix(2,i) -
     $             massfrac*fixcoeff(6,iwhich)
              if (fixflag(3,iwhich) == 0) f_fix(3,i) = f_fix(3,i) -
     $             massfrac*fixcoeff(7,iwhich)
            else if (istyle == 8) then
              dx = x(1,i) - fixcoeff(1,iwhich)
              dy = x(2,i) - fixcoeff(2,iwhich)
              dz = x(3,i) - fixcoeff(3,iwhich)
              if (fixflag(1,iwhich) == 1) dx = 0.0
              if (fixflag(2,iwhich) == 1) dy = 0.0
              if (fixflag(3,iwhich) == 1) dz = 0.0
              call minimg(dx,dy,dz)
              r = sqrt(dx*dx + dy*dy + dz*dz)
              if (r > fixcoeff(5,iwhich)) then
                pre = fixcoeff(4,iwhich)/r
                f_fix(1,i) = f_fix(1,i) - pre*dx
                f_fix(2,i) = f_fix(2,i) - pre*dy
                f_fix(3,i) = f_fix(3,i) - pre*dz
              endif
            endif

            itmp = itmp/maxfix

          enddo
        enddo

      else

        do i = 1,nlocal
          itmp = fix(i)
          do while (itmp > 0)

            iwhich = mod(itmp,maxfix)
            if (fixactive(iwhich) == 0) then
              itmp = itmp/maxfix
              cycle
            endif
            istyle = fixstyle(iwhich)

            if (istyle == 1) then
              if (fixflag(1,iwhich) == 0)
     $             f_fix(1,i) = fixcoeff(1,iwhich)
              if (fixflag(2,iwhich) == 0)
     $             f_fix(2,i) = fixcoeff(2,iwhich)
            else if (istyle == 2) then
              f_fix(1,i) = f_fix(1,i) + fixcoeff(1,iwhich)
              f_fix(2,i) = f_fix(2,i) + fixcoeff(2,iwhich)
            else if (istyle == 3) then
              f_fix(1,i) = fixcoeff(4,iwhich)
              f_fix(2,i) = fixcoeff(5,iwhich)
            else if (istyle == 4) then
              factor = fixcoeff(6,iwhich)
              v(1,i) = v(1,i)*factor
              v(2,i) = v(2,i)*factor
            else if (istyle == 5) then
              factor = fixcoeff(4,iwhich)
              dtmass = mass(type(i))/dt_fix
              f_fix(1,i) = f_fix(1,i) - dtmass*factor*v(1,i)
              f_fix(2,i) = f_fix(2,i) - dtmass*factor*v(2,i)
            else if (istyle == 6) then
              rmass = mass(type(i))
              gamma1 = -fixcoeff(3,iwhich)*rmass
              gamma2 = fixcoeff(4,iwhich)*sqrt(rmass)
              if (fixflag(1,iwhich) == 1)
     $             f_fix(1,i) = f_fix(1,i) + gamma1*v(1,i) +
     $             gamma2*(ranmars(0)-0.5)
              if (fixflag(2,iwhich) == 1)
     $             f_fix(2,i) = f_fix(2,i) + gamma1*v(2,i) +
     $             gamma2*(ranmars(0)-0.5)
            else if (istyle == 7) then
              massfrac = mass(type(i))/fixmass(iwhich)
              if (fixflag(1,iwhich) == 0) f_fix(1,i) = f_fix(1,i) -
     $             massfrac*fixcoeff(5,iwhich)
              if (fixflag(2,iwhich) == 0) f_fix(2,i) = f_fix(2,i) -
     $             massfrac*fixcoeff(6,iwhich)
            else if (istyle == 8) then
              dx = x(1,i) - fixcoeff(1,iwhich)
              dy = x(2,i) - fixcoeff(2,iwhich)
              dz = x(3,i) - fixcoeff(3,iwhich)
              if (fixflag(1,iwhich) == 1) dx = 0.0
              if (fixflag(2,iwhich) == 1) dy = 0.0
              if (fixflag(3,iwhich) == 1) dz = 0.0
              call minimg(dx,dy,dz)
              r = sqrt(dx*dx + dy*dy + dz*dz)
              if (r > fixcoeff(5,iwhich)) then
                pre = fixcoeff(4,iwhich)/r
                f_fix(1,i) = f_fix(1,i) - pre*dx
                f_fix(2,i) = f_fix(2,i) - pre*dy
              endif
            endif

            itmp = itmp/maxfix

          enddo
        enddo

      endif

      return
      end
