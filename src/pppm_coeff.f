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
c setup of 2 decompositions (3-d brick and 2-d FFT) and various coeffs
c   for PPPM calculation
c iflag = 0 section is only computed once since those quantities
c   depend only on grid size, it is called at beginning of a run
c NOTE: if box volume changes enough that a different resolution grid
c   is needed then need to put in calls to recompute iflag = 0 section
c iflag = 1 section must be recomputed whenever simulation volume changes,
c   it is called during a run when a volume change occurs

      subroutine pppm_coeff(iflag)
      use global
      use mpi
      implicit none

c argument variables

      integer iflag

c local variables

      integer istatus(mpi_status_size)
      integer nlo,nhi,itotal,ntotal,ierror,nplanes,jflag
      integer kflag,nx,ny,nz,idelx,idely,idelz,nxx,nyy
      integer npey_fft,npez_fft,me_y,me_z,itmp,i,j,k,kper,l
      integer lper,m,mper,nzz,nfft_brick,icount
      real*8 delxinv,delyinv,delzinv,pshift,cuthalf
      real*8 unitk,sqk,vterm

c adjustment of z dimension for 2d slab PPPM
c 3d PPPM just uses zprd since volfactor = 1.0

      zprd_slab = zprd*slab_volfactor

      if (iflag == 0) then

c PPPM grid size and cutoff parameter from estimation routine

        call get_p3m_parms

c nlo_in,nhi_in = lower/upper limits of the 3-d sub-brick of
c   global PPPM grid that I own in each dimension WITHOUT ghost cells
c global indices range from 0 to N-1
c for slab PPPM, assign z grid as if it were not extended

        nxlo_in = me(1)*nx_pppm/pgrid(1)
        nxhi_in = (me(1)+1)*nx_pppm/pgrid(1) - 1
        nylo_in = me(2)*ny_pppm/pgrid(2)
        nyhi_in = (me(2)+1)*ny_pppm/pgrid(2) - 1
        nzlo_in = me(3)*nint(nz_pppm/slab_volfactor)/pgrid(3)
        nzhi_in = (me(3)+1)*nint(nz_pppm/slab_volfactor)/pgrid(3) - 1

c nlower,nupper = stencil size for mapping particles to PPPM grid

        nlower = -(orderflag-1)/2
        nupper = orderflag/2

c nlo_out,nhi_out = lower/upper limits of the 3-d sub-brick of
c   global PPPM grid that I own in each dimension WITH ghost cells
c ghost 3-d brick boundaries = inner + stencil + particles moving skin/2.0
c nlo/nhi = global coords of grid pt to "lower left" of smallest/largest
c           position a particle in my box can be at
c add/subtract 4096 to avoid problem of int(-0.75) = 0 when it needs to be -1

        delxinv = nx_pppm/xprd
        delyinv = ny_pppm/yprd
        delzinv = nz_pppm/zprd_slab

        if (mod(orderflag,2).ne.0) then
          pshift = 4096.5
        else
          pshift = 4096.0
        endif

        cuthalf = skin/2.0

        nlo = int((border(1,1)-cuthalf-box(1,1))*delxinv+pshift) - 4096
        nhi = int((border(2,1)+cuthalf-box(1,1))*delxinv+pshift) - 4096
        nxlo_out = nlo + nlower
        nxhi_out = nhi + nupper

        nlo = int((border(1,2)-cuthalf-box(1,2))*delyinv+pshift) - 4096
        nhi = int((border(2,2)+cuthalf-box(1,2))*delyinv+pshift) - 4096
        nylo_out = nlo + nlower
        nyhi_out = nhi + nupper

        nlo = int((border(1,3)-cuthalf-box(1,3))*delzinv+pshift) - 4096
        nhi = int((border(2,3)+cuthalf-box(1,3))*delzinv+pshift) - 4096
        nzlo_out = nlo + nlower
        nzhi_out = nhi + nupper

c for slab-geometry-corrected calculations, change the grid boundary for
c   processors at +z end to include the empty volume between periodically
c   repeating slabs
c for slab PPPM, want charge data communicated from -z proc to +z proc,
c   but not vice versa, also want field data communicated from +z proc to
c   -z proc, but not vice versa
c this is accomplished by nzhi_in = nzhi_out on +z end (no ghost cells)

        if (slabflag == 1 .and. me(3)+1 == pgrid(3)) then
          nzhi_in =  nz_pppm - 1
          nzhi_out = nz_pppm - 1
        endif

c exchange messages to determine how many planes I will send,recv in each dir
c if no neighbor proc exists, values comes from self
c   since I have ghosts regardless

        nplanes = nxlo_in - nxlo_out
        if (mpart(1,1).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,1),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nxhi_ghost,1,mpi_integer,mpart(2,1),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nxhi_ghost = nplanes
        endif

        nplanes = nxhi_out - nxhi_in
        if (mpart(2,1).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,1),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nxlo_ghost,1,mpi_integer,mpart(1,1),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nxlo_ghost = nplanes
        endif

        nplanes = nylo_in - nylo_out
        if (mpart(1,2).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,2),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nyhi_ghost,1,mpi_integer,mpart(2,2),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nyhi_ghost = nplanes
        endif

        nplanes = nyhi_out - nyhi_in
        if (mpart(2,2).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,2),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nylo_ghost,1,mpi_integer,mpart(1,2),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nylo_ghost = nplanes
        endif

        nplanes = nzlo_in - nzlo_out
        if (mpart(1,3).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(1,3),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nzhi_ghost,1,mpi_integer,mpart(2,3),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nzhi_ghost = nplanes
        endif

        nplanes = nzhi_out - nzhi_in
        if (mpart(2,3).ne.node) then
          call mpi_send(nplanes,1,mpi_integer,mpart(2,3),0,
     $         mpi_comm_world,ierror)
          call mpi_recv(nzlo_ghost,1,mpi_integer,mpart(1,3),0,
     $         mpi_comm_world,istatus,ierror)
        else
          nzlo_ghost = nplanes
        endif

c test that ghost overlap is not bigger than my sub-domain

        jflag = 0
        if (nxlo_ghost.gt.nxhi_in-nxlo_in+1) jflag = 1
        if (nxhi_ghost.gt.nxhi_in-nxlo_in+1) jflag = 1
        if (nylo_ghost.gt.nyhi_in-nylo_in+1) jflag = 1
        if (nyhi_ghost.gt.nyhi_in-nylo_in+1) jflag = 1
        if (nzlo_ghost.gt.nzhi_in-nzlo_in+1) jflag = 1
        if (nzhi_ghost.gt.nzhi_in-nzlo_in+1) jflag = 1

        call mpi_allreduce(jflag,kflag,1,mpi_integer,mpi_sum,
     $       mpi_comm_world,ierror)

        if (kflag.gt.0)
     $       call err('PPPM stencil extends beyond neighbor proc'//
     $       ' - reduce PPPM order')

c decomposition of FFT mesh
c proc owns entire x-dimension, clump of columns in y,z dimensions
c npey_fft,npez_fft = # of procs in y,z dims
c if nprocs is small enough, proc can own 1 or more entire xy planes,
c   else proc owns 2-d sub-blocks of yz plane
c me_y,me_z = which proc (0-npe_fft-1) I am in y,z dimensions
c nlo_fft,nhi_fft = lower/upper limit of the section
c   of the global FFT mesh that I own in each dimension
c indices range from 0 to N-1

        if (nz_pppm >= nprocs) then
          npey_fft = 1
          npez_fft = nprocs
        else
          call proc2grid2d(nprocs,ny_pppm,nz_pppm,npey_fft,npez_fft)
        endif

        me_y = mod(node,npey_fft)
        me_z = node/npey_fft

        nxlo_fft = 0
        nxhi_fft = nx_pppm - 1
        nylo_fft = me_y*ny_pppm/npey_fft
        nyhi_fft = (me_y+1)*ny_pppm/npey_fft - 1
        nzlo_fft = me_z*nz_pppm/npez_fft
        nzhi_fft = (me_z+1)*nz_pppm/npez_fft - 1

c allocate all memory needed for PPPM solution

c allocate grid arrays for PPPM mesh on this proc, including ghosts

        maxgrid = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
     $       (nzhi_out-nzlo_out+1)

        call mpi_allreduce(maxgrid,ntotal,1,mpi_integer,
     $       mpi_max,mpi_comm_world,ierror)
        if (node == 0) then
          write (6,*) 'PPPM max grid size =',ntotal
          write (1,*) 'PPPM max grid size =',ntotal
        endif

        allocate (density_brick(maxgrid))
        allocate (vdx_brick(maxgrid))
        allocate (vdy_brick(maxgrid))
        allocate (vdz_brick(maxgrid))

c allocate FFT arrays for PPPM mesh on this proc, without ghosts
c nfft = FFT points in FFT decomposition on this proc
c nfft_brick = FFT points in 3-d brick-decomposition on this proc

        nfft = (nxhi_fft-nxlo_fft+1) * (nyhi_fft-nylo_fft+1) *
     $       (nzhi_fft-nzlo_fft+1)
        nfft_brick = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) *
     $       (nzhi_in-nzlo_in+1)
        maxfft = max(nfft,nfft_brick)

        call mpi_allreduce(maxfft,ntotal,1,mpi_integer,
     $       mpi_max,mpi_comm_world,ierror)
        if (node == 0) then
          write (6,*) 'PPPM max FFT size =',ntotal
          write (1,*) 'PPPM max FFT size =',ntotal
        endif

        allocate (density_fft(maxfft))
        allocate (greensfn(maxfft))
        allocate (workvec1(maxfft))
        allocate (workvec2(maxfft))
        allocate (vg(6,maxfft))

c allocate per-particle arrays

        allocate (partgrid(3,maxown))
        allocate (ek(3,maxown))

c allocate space for fkvecs_xyz factors

        allocate (fkvecs_x(nxlo_fft:nxhi_fft))
        allocate (fkvecs_y(nylo_fft:nyhi_fft))
        allocate (fkvecs_z(nzlo_fft:nzhi_fft))

c allocate buffer space for use in brick2fft and fillbrick
c idel = max # of ghost planes to send or recv in +/- direction of each dim
c nxx,nyy,nzz = max # of cells to send in each dim
c maxpbuf = max in any dim, factor of 3 for components of vd_xyz in fillbrick

        nx = nxhi_out - nxlo_out + 1
        ny = nyhi_out - nylo_out + 1
        nz = nzhi_out - nzlo_out + 1

        idelx = max(nxlo_ghost,nxhi_ghost)
        idelx = max(idelx,nxhi_out-nxhi_in)
        idelx = max(idelx,nxlo_in-nxlo_out)

        idely = max(nylo_ghost,nyhi_ghost)
        idely = max(idely,nyhi_out-nyhi_in)
        idely = max(idely,nylo_in-nylo_out)

        idelz = max(nzlo_ghost,nzhi_ghost)
        idelz = max(idelz,nzhi_out-nzhi_in)
        idelz = max(idelz,nzlo_in-nzlo_out)

        nxx = idelx * ny * nz
        nyy = idely * nx * nz
        nzz = idelz * nx * ny

        itotal = max(nxx,nyy)
        itotal = max(itotal,nzz)
        maxpbuf = 3*itotal

        call mpi_allreduce(maxpbuf,ntotal,1,mpi_integer,
     $       mpi_max,mpi_comm_world,ierror)
        if (node == 0) then
          write (6,*) 'PPPM max buffer size =',ntotal
          write (1,*) 'PPPM max buffer size =',ntotal
        endif

        allocate (pbuf1(maxpbuf),pbuf2(maxpbuf))

c create 2 FFT plans and remap plan
c 1st FFT plan keeps data in FFT decompostion
c 2nd FFT plan returns data in 3-d brick decomposition
c remap plan takes data from 3-d brick to FFT decomposition
c all indices must be converted to range from 1 to N

        call fft_3d_create_plan(mpi_comm_world,
     $       nx_pppm,ny_pppm,nz_pppm,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       0,0,itmp,plan1_fft)

        call fft_3d_create_plan(mpi_comm_world,
     $       nx_pppm,ny_pppm,nz_pppm,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       nxlo_in+1,nxhi_in+1,nylo_in+1,nyhi_in+1,
     $       nzlo_in+1,nzhi_in+1,
     $       0,0,itmp,plan2_fft)

        call remap_3d_create_plan(mpi_comm_world,
     $       nxlo_in+1,nxhi_in+1,nylo_in+1,nyhi_in+1,
     $       nzlo_in+1,nzhi_in+1,
     $       nxlo_fft+1,nxhi_fft+1,nylo_fft+1,nyhi_fft+1,
     $       nzlo_fft+1,nzhi_fft+1,
     $       1,0,0,2,plan_remap)

      endif

c computations that must be re-done whenever simulation volume changes

c pre-compute fkvecs_xyz(*) for my FFT grid pts

      unitk = (2.0*3.141592654/xprd)
      do k = nxlo_fft,nxhi_fft
        kper = k - nx_pppm*int(2*k/nx_pppm)
        fkvecs_x(k) = unitk*kper
      enddo

      unitk = (2.0*3.141592654/yprd)
      do l = nylo_fft,nyhi_fft
        lper = l - ny_pppm*int(2*l/ny_pppm)
        fkvecs_y(l) = unitk*lper
      enddo

      unitk = (2.0*3.141592654/zprd_slab)
      do m = nzlo_fft,nzhi_fft
        mper = m - nz_pppm*int(2*m/nz_pppm)
        fkvecs_z(m) = unitk*mper
      enddo

c pre-compute virial coefficients

      icount = 0
      do k = nzlo_fft,nzhi_fft
        do j = nylo_fft,nyhi_fft
          do i = nxlo_fft,nxhi_fft
            icount = icount + 1
            sqk = fkvecs_x(i)*fkvecs_x(i) +
     $            fkvecs_y(j)*fkvecs_y(j) +
     $            fkvecs_z(k)*fkvecs_z(k)
            if (sqk == 0.0) then
              vg(1,icount) = 0.0
              vg(2,icount) = 0.0
              vg(3,icount) = 0.0
              vg(4,icount) = 0.0
              vg(5,icount) = 0.0
              vg(6,icount) = 0.0
            else
              vterm = -2.0*(1.0/sqk + 0.25/(gewald*gewald))
              vg(1,icount) = 1.0 + vterm*fkvecs_x(i)*fkvecs_x(i)
              vg(2,icount) = 1.0 + vterm*fkvecs_y(j)*fkvecs_y(j)
              vg(3,icount) = 1.0 + vterm*fkvecs_z(k)*fkvecs_z(k)
              vg(4,icount) = vterm*fkvecs_x(i)*fkvecs_y(j)
              vg(5,icount) = vterm*fkvecs_x(i)*fkvecs_z(k)
              vg(6,icount) = vterm*fkvecs_y(j)*fkvecs_z(k)
            endif
          enddo
        enddo
      enddo

c pre-compute Green's function coeffs

      call make_greensfn(orderflag,gewald,
     $     nx_pppm,ny_pppm,nz_pppm,xprd,yprd,zprd_slab,
     $     nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft,
     $     greensfn)

      return
      end

c -------------------------------------------------------------------------
c get p3m parameters
c   INPUT:
c           orderflag------------------p3m assignment scheme order
c           long_prec------------------desired relative error in forces
c           cutcoul--------------------real space cutoff
c           xprd,yprd,zprd_slab--------edge lengths of periodic cell
c           natoms---------------------total number of atoms
c           qsqsum---------------------sum of the squared charges
c   OUTPUT:
c           gewald---------------------Ewald parameter
c           nx_pppm,ny_pppm,nz_pppm----P3m grid sizes
c

      subroutine get_p3m_parms
      use global
      implicit none

c local variables (old method)

      integer maxorder
      real*8 a,b,gdel
      parameter (maxorder=16)
      dimension a(maxorder),b(maxorder)
      data a /5.0,5.78,6.93,7.76,8.55,9.10,9.66,10.17,10.60,
     $     11.1,14,15,16,17,18,19/
      data b /2.0,2.25,3.41,4.62,5.85,7.05,8.20,9.28,10.34,
     $     11.3,12,13,14,15,16,17/

c local variables (new method)

      logical factorable,usenew
      integer nfacs,factors(3),ncount,m
      real*8 small,large,acons(7,0:6),h,h1,h2,er,er1,er2,gew1,gew2
      real*8 hx,hy,hz,sum,spr,lpr,lprx,lpry,lprz,q2

      parameter (small=0.00001,large=10000)

c comment one of next 2 lines out
c 1st is for power-of-2 FFTs
c 2nd is for non-power-of-2 FFTs

c      data nfacs /1/, factors /2,0,0/
      data nfacs /3/, factors /2,3,5/

c if problems arise with the "new" method, set usenew = .false.

      usenew = .true.

      if (usenew) then

        q2 = qsqsum/dielectric

c see JCP 109, pg. 7698 for derivation of coefficients
c higher order coefficients may be computed if needed

        acons(1,0) = 2./3.
        acons(2,0) = 1./50.
        acons(2,1) = 5./294.
        acons(3,0) = 1./588.
        acons(3,1) = 7./1440.
        acons(3,2) = 21./3872.
        acons(4,0) = 1./4320.
        acons(4,1) = 3./1936.
        acons(4,2) = 7601./2271360.
        acons(4,3) = 143./28800.
        acons(5,0) = 1./23232.
        acons(5,1) = 7601./13628160.
        acons(5,2) = 143./69120.
        acons(5,3) = 517231./106536960.
        acons(5,4) = 106640677./11737571328.
        acons(6,0) = 691./68140800.
        acons(6,1) = 13./57600.
        acons(6,2) = 47021./35512320.
        acons(6,3) = 9694607./2095994880.
        acons(6,4) = 733191589./59609088000.
        acons(6,5) = 326190917./11700633600.
        acons(7,0) = 1./345600.
        acons(7,1) = 3617./35512320.
        acons(7,2) = 745739./838397952.
        acons(7,3) = 56399353./12773376000.
        acons(7,4) = 25091609./1560084480.
        acons(7,5) = 1755948832039./36229939200000.
        acons(7,6) = 4887769399./37838389248.

c error check, if higher order desired, add needed acons values above,
c or set usenew = .false.

        if (orderflag > 7)
     $      call err('PPPM order cannot be greater than 7')

c compute necessary gewald
c based on desired error and real space cutoff
c fluid-occupied volume used to estimate real-space error

        gewald = sqrt(-log(long_prec*
     $    sqrt(natoms*cutcoul*xprd*yprd*zprd)/(2*q2)))/cutcoul

c compute optimal nx_pppm,ny_pppm,nz_pppm
c based on assignment order and desired error

        h = 1.0
        h1 = 2.0
        ncount = 0
        er = large
        er1 = 0
        do while (abs(er) > small)
           sum = 0.0
           do m = 0,orderflag-1
              sum = sum + acons(orderflag,m)*(h*gewald)**(2*m)
           end do
           lpr = q2*(h*gewald)**orderflag *
     $          sqrt(gewald*xprd*sqrt(2.0D0*3.141592654)*sum/natoms)/
     $          (xprd*xprd)
           er = log(lpr) - log(long_prec)
           er2 = er1
           er1 = er
           h2 = h1
           h1 = h
           if ((er1 - er2) == 0.0) then
             h = h1 + er1
           else
             h = h1 + er1*(h2 - h1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call err('can''t get x-mesh spacing')
        end do
        nx_pppm = xprd/h + 1

        ncount = 0
        er = large
        er1 = 0
        do while (abs(er) > small)
           sum = 0.0
           do m = 0,orderflag-1
              sum = sum + acons(orderflag,m)*(h*gewald)**(2*m)
           end do
           lpr = q2*(h*gewald)**orderflag *
     $          sqrt(gewald*yprd*sqrt(2.0D0*3.141592654)*sum/natoms)/
     $          (yprd*yprd)
           er = log(lpr) - log(long_prec)
           er2 = er1
           er1 = er
           h2 = h1
           h1 = h
           if ((er1 - er2) == 0.0) then
             h = h1 + er1
           else
             h = h1 + er1*(h2 - h1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call err('can''t get y-mesh spacing')
        end do
        ny_pppm = yprd/h + 1

        ncount = 0
        er = large
        er1 = 0
        do while (abs(er) > small)
           sum = 0.0
           do m = 0,orderflag-1
              sum = sum + acons(orderflag,m)*(h*gewald)**(2*m)
           end do
           lpr = q2*(h*gewald)**orderflag *
     $          sqrt(gewald*zprd_slab*sqrt(2.0*3.141592654)*sum/natoms)/
     $          (zprd_slab*zprd_slab)
           er = log(lpr) - log(long_prec)
           er2 = er1
           er1 = er
           h2 = h1
           h1 = h
           if ((er1 - er2) == 0.0) then
             h = h1 + er1
           else
             h = h1 + er1*(h2 - h1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call err('can''t get z-mesh spacing')
        end do
        nz_pppm = zprd_slab/h + 1

c use the old method to find the needed parameters

      else

c error check - should be maxorder, but has bug for > 6

c        if (orderflag > maxorder)
c     $       call err('PPPM order cannot be larger than maxorder')
        if (orderflag > 6)
     $       call err('PPPM order cannot be larger than 6')

c compute necessary gewald based on desired error and real space cutoff
c compute optimal nx_pppm,ny_pppm,nz_pppm based on assignment order
c and desired error

        gewald = sqrt(2.0)*(1.35-.15*log(long_prec))/cutcoul
        gdel = (long_prec/(2.5*exp(-a(orderflag))))**(1.0/b(orderflag))
        nx_pppm = gewald*xprd/gdel
        ny_pppm = gewald*yprd/gdel
        nz_pppm = gewald*zprd_slab/gdel
        lpr = long_prec

      endif

c convert PPPM grid size into "good" numbers - factorable by 2 or 2,3,5

      do while (.not.factorable(nx_pppm,nfacs,factors))
        nx_pppm = nx_pppm + 1
      enddo
      do while (.not.factorable(ny_pppm,nfacs,factors))
        ny_pppm = ny_pppm + 1
      enddo
      do while (.not.factorable(nz_pppm,nfacs,factors))
        nz_pppm = nz_pppm + 1
      enddo

c overwrite PPPM grid size with input values

      if (meshflag.eq.1) then
        nx_pppm = nx_pppm_input
        ny_pppm = ny_pppm_input
        nz_pppm = nz_pppm_input
      endif

c print-out actual PPPM grid to be used

      if (node.eq.0) then
        write (6,*) 'Actual PPPM grid =',
     $       nx_pppm,ny_pppm,nz_pppm
        write (1,*) 'Actual PPPM grid =',
     $       nx_pppm,ny_pppm,nz_pppm
      endif

      if (usenew) then

c now that the actual grid has been determined, find best ewald param

        gew1 = gewald + 0.01
        ncount = 0
        er = large
        er1 = 0
        hx = xprd/nx_pppm
        hy = yprd/ny_pppm
        hz = zprd_slab/nz_pppm
        do while (abs(er) > small)
           sum = 0.0
           do m = 0,orderflag-1
              sum = sum + acons(orderflag,m)*(hx*gewald)**(2*m)
           enddo
           lprx = q2*(hx*gewald)**orderflag *
     $          sqrt(gewald*xprd*sqrt(2.0D0*3.141592654)*sum/natoms)/
     $          (xprd*xprd)

           sum = 0.0
           do m = 0,orderflag-1
              sum = sum + acons(orderflag,m)*(hy*gewald)**(2*m)
           enddo
           lpry = q2*(hy*gewald)**orderflag *
     $          sqrt(gewald*yprd*sqrt(2.0D0*3.141592654)*sum/natoms)/
     $          (yprd*yprd)

           sum = 0.0
           do m = 0,orderflag-1
              sum = sum + acons(orderflag,m)*(hz*gewald)**(2*m)
           enddo
           lprz = q2*(hz*gewald)**orderflag *
     $          sqrt(gewald*zprd_slab*sqrt(2.0D0*3.141592654)*
     $          sum/natoms)/(zprd_slab*zprd_slab)
           lpr = sqrt(lprx*lprx +
     $                lpry*lpry + lprz*lprz)/sqrt(3.)
           spr = 2.*q2*exp(-gewald*gewald*cutcoul*cutcoul)/
     $              sqrt(natoms*cutcoul*xprd*yprd*zprd_slab)
           er = log(lpr) - log(spr)
           er2 = er1
           er1 = er
           gew2 = gew1
           gew1 = gewald
           if ((er1 - er2) == 0.0) then
             gewald = gew1 + er1
           else
             gewald = gew1 + er1*(gew2 - gew1)/(er1 - er2)
           endif
           ncount = ncount + 1
           if (ncount > large) call err('Trouble finding Ewald G')
        end do

c h*gewald must be small in order for this "new" error estimation routine
c  to be completely valid. See JCP 109, 7698. However, this routine
c  may produce a reasonable guess even if h*gewald > 1

        if ((hx*gewald > 1) .or. (hy*gewald > 1) .or.
     $       (hz*gewald > 1)) then
          if (node.eq.0) write (6,*) 'WARNING: Invalid assumption ',
     $         'that h*gewald is small was used in choice of PPPM G'
        endif

      endif

      if (node.eq.0) then
        write (6,*) 'PPPM G =',gewald
        write (1,*) 'PPPM G =',gewald
        write (6,*) 'Expected RMS precision =', lpr
        write (1,*) 'Expected RMS precision =', lpr
      endif

      return
      end


c -------------------------------------------------------------------------
c set up modified (Hockney-Eastwood) Coulomb Green's Fn
c        (replaces 4*pi/k**2) for k-vectors on this PE
c     INPUT:
c         na_ord=order of assignment scheme
c         gewald=(G-Ewald) width parameter for Gaussians
c         nx_pppm,ny_pppm,nz_pppm=grid array dimensions
c         cellz,celly,cellz=grid sizes
c         nx1,nx2 = lower, upper x limit of k-space arrays for this PE
c         ny1,ny2 = lower, upper y limit of k-space arrays for this PE
c         nz1,nz2 = lower, upper z limit of k-space arrays for this PE
c     OUTPUT:
c         greensfn=Coulomb Green's Fn for k-vectors on this PE
c     DATA:
c         epshoc=convergence parameter for sums
c     CALLS: gf_denom

      subroutine make_greensfn(na_ord,gewald,nx_pppm,ny_pppm,
     $     nz_pppm,cellx,celly,cellz,nx1,nx2,ny1,ny2,nz1,nz2,greensfn)
      implicit none

c argument variables

      integer na_ord,nx_pppm,ny_pppm,nz_pppm
      integer nx1,nx2,ny1,ny2,nz1,nz2
      real*8 gewald,cellx,celly,cellz,greensfn
      dimension greensfn(*)

c local variables

      integer nbx,nby,nbz,indx,m,mper,l,lper,k,kper,nx,ny
      integer nz
      real*8 epshoc,shpfn,arg,unitkx,unitky,unitkz,form
      real*8 snz2,sny2,snx2,sqk,sqk_perp,sum1,qx,sx,wx
      real*8 argx,qy,sy,wy,argy,qz,sz,wz,argz,dot1,dot2
      real*8 sum2sq
      data epshoc /1.0D-7/

c implicit function

      shpfn(arg) = exp(-.25*(arg/gewald)**2)

      unitkx = (2.*3.141592654/cellx)
      unitky = (2.*3.141592654/celly)
      unitkz = (2.*3.141592654/cellz)

      nbx = (gewald*cellx/(3.141592654*nx_pppm))*((-log(epshoc))**.25)
      nby = (gewald*celly/(3.141592654*ny_pppm))*((-log(epshoc))**.25)
      nbz = (gewald*cellz/(3.141592654*nz_pppm))*((-log(epshoc))**.25)

      form = 1.0
      indx = 1
      do m = nz1,nz2
        mper = m-nz_pppm*int(2*m/nz_pppm)
        snz2 = (sin(.5*unitkz*mper*cellz/nz_pppm))**2
        do l = ny1,ny2
          lper = l-ny_pppm*int(2*l/ny_pppm)
          sny2 = (sin(.5*unitky*lper*celly/ny_pppm))**2
          do k = nx1,nx2

            kper = k-nx_pppm*int(2*k/nx_pppm)
            snx2 = (sin(.5*unitkx*kper*cellx/nx_pppm))**2
            sqk = ((unitkx*kper)**2+(unitky*lper)**2+(unitkz*mper)**2)
            sqk_perp = ((unitkx*kper)**2+(unitky*lper)**2)

            if (sqk.ne.0.0) then
              greensfn(indx) = form*12.5663706/sqk
              sum1 = 0.0
              do nx = -nbx,nbx
                qx = unitkx*(kper+nx_pppm*nx)
                sx = shpfn(qx)
                wx = 1.0
                argx = .5*qx*cellx/nx_pppm
                if (argx.ne.0.0) wx = (sin(argx)/argx)**na_ord
                do ny = -nby,nby
                  qy = unitky*(lper+ny_pppm*ny)
                  sy = shpfn(qy)
                  wy = 1.0
                  argy = .5*qy*celly/ny_pppm
                  if (argy.ne.0.0) wy = (sin(argy)/argy)**na_ord
                  do nz = -nbz,nbz
                    qz = unitkz*(mper+nz_pppm*nz)
                    sz = shpfn(qz)
                    wz = 1.0
                    argz = .5*qz*cellz/nz_pppm
                    if (argz.ne.0.0) wz = (sin(argz)/argz)**na_ord

                    dot1 = unitkx*kper*qx +
     $                   unitky*lper*qy+unitkz*mper*qz
                    dot2 = qx*qx+qy*qy+qz*qz
                    sum1 = sum1+(dot1/dot2)*sx*sy*sz*(wx*wy*wz)**2
                  enddo
                enddo
              enddo
              call gf_denom(na_ord,snx2,sny2,snz2,sum2sq)
              greensfn(indx) = greensfn(indx)*sum1/sum2sq
            else
              greensfn(indx) = 0.0
            endif

            indx = indx+1
          enddo
        enddo
      enddo

      return
      end

c -------------------------------------------------------------------------
c returns denominator for Hockney-Eastwood Green's function
c               inf                 n-1
c      S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
c              j=-inf               l=0
c
c       = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z=sin(x)
c     INPUT:
c           n ====== order of assignment scheme
c           x,y,z=== sin(kx*deltax/2), etc.
c           b ====== denominator expansion coefficients /Gamma(2n)
c     OUTPUT:
c           s ====== denominator of Hockney-Eastwood expression

      subroutine gf_denom(n,x,y,z,s)
      implicit none

c argument variables

      integer n
      real*8 x,y,z,s

c local variables

      integer maxorder,nlast,k,ifact,l,m
      real*8 b,gaminv,sx,sy,sz
      parameter (maxorder=16)
      dimension b(0:maxorder-1)
      save b
      data nlast /0/
      save nlast

      if (n.ne.nlast) then
        ifact = 1
        do k = 1,2*n-1
          ifact = ifact*k
        enddo
        gaminv = 1./float(ifact)
        do l = 0,n-1
          b(l) = 0.0
        enddo
        b(0) = 1.0

        do m = 1,n-1
          do l = m,1,-1
            b(l) = 4.0*(b(l)*(l-m)*(l-m-.5)-b(l-1)*(l-1-m)**2)
          enddo
          b(0) = 4.0*(b(0)*(l-m)*(l-m-.5))
        enddo
        do l = 0,n-1
          b(l) = b(l)*gaminv
        enddo
        nlast = n
      endif

      sx = 0.0
      sy = 0.0
      sz = 0.0
      do l = n-1,0,-1
        sx = b(l) + sx*x
        sy = b(l) + sy*y
        sz = b(l) + sz*z
      enddo
      s = sx*sy*sz
      s = s*s

      return
      end

c -------------------------------------------------------------------------
c check if all the factors of n are in factors(1:nfacs)
c return TRUE if they are
c return FALSE if they are not

      logical function factorable(n,nfacs,factors)
      implicit none

c argument variables

      integer n,nfacs,factors(nfacs)

c local variables

      integer nremain,i,j,iflag
      integer nlist,list(20)

c build list of factors of n
c nlist = # of factors
c after inner while loop, i is a factor of n
c nremain = remaining portion of n that has not been factored

      nlist = 0
      nremain = n

      do while (nremain.gt.1)
        i = 2
        do while (mod(nremain,i).ne.0)
          i = i + 1
        enddo
        nlist = nlist + 1
        list(nlist) = i
        nremain = nremain/i
      enddo

c check if every factor of n in list is also is in input factors
c as soon as one list item is not in factors, return FALSE

      do i = 1,nlist
        iflag = 0
        do j = 1,nfacs
          if (list(i).eq.factors(j)) iflag = 1
        enddo
        if (iflag.eq.0) then
          factorable = .FALSE.
          return
        endif
      enddo

      factorable = .TRUE.

      return
      end

c -------------------------------------------------------------------------
c deallocate all PPPM memory, called at end of a run

      subroutine pppm_deallocate
      use global
      implicit none

      if (allocated(density_brick)) then
        deallocate (density_brick)
        deallocate (vdx_brick,vdy_brick,vdz_brick)
        deallocate (density_fft,greensfn)
        deallocate (workvec1,workvec2,vg)
        deallocate (partgrid,ek)
        deallocate (fkvecs_x,fkvecs_y,fkvecs_z)
        deallocate (pbuf1,pbuf2)
      endif

      return
      end
