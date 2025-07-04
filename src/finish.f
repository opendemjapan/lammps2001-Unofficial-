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
c end of timestepping

      subroutine finish
      use global
      use mpi
      implicit none

c local variables

      integer ierror,i,indx,itmp1,itmp2,itmp,max_bin
      integer ihisto(10),ihistotmp(10)
      real*8 rtmp,ave,xmax,xmin,time1,time2,r8tmp,all
      real*8 opt_own,opt_ghost,opt_neigh,opt_buf

      external fft_3d,fft_3d_destroy_plan
      external remap_3d_destroy_plan

 900  format(' ',a,f15.6,f13.4)
 901  format(' ',a,f13.4,a,f13.4,a,f13.4,a)
 902  format(' ',a,10i5)
 903  format(' ',a,f12.6)
 904  format(4f10.2)

c timings stats

      time_other = time_loop - (time_nonbond+time_long+time_bond+
     $     time_angle+time_dihedral+time_improper+time_neigh1+
     $     time_neigh2+time_shake+time_exch+time_comm+time_fcomm+
     $     time_io)

      rtmp = time_loop
      call mpi_allreduce(rtmp,time_loop,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      time_loop = time_loop/nprocs
      if (node.eq.0) then
        write (6,*) 'Loop time:',time_loop,
     $       ' on',nprocs,' procs for',natoms,' atoms'
        write (1,*) 'Loop time:',time_loop,
     $       ' on',nprocs,' procs for',natoms,' atoms'
      endif

      if (time_loop.eq.0.0) time_loop = 1.0

      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      call mpi_allreduce(time_nonbond,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Nbond time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Nbond time/%:',rtmp,rtmp/time_loop*100
      endif

      if (coulstyle == 3 .or. coulstyle == 4) then
        call mpi_allreduce(time_long,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Long  time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Long  time/%:',rtmp,rtmp/time_loop*100
        endif
      endif

      if (nmolecular == 1) then

        call mpi_allreduce(time_bond,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Bond  time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Bond  time/%:',rtmp,rtmp/time_loop*100
        endif

        call mpi_allreduce(time_angle,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Angle time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Angle time/%:',rtmp,rtmp/time_loop*100
        endif

        call mpi_allreduce(time_dihedral,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Dihed time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Dihed time/%:',rtmp,rtmp/time_loop*100
        endif

        call mpi_allreduce(time_improper,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Impro time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Impro time/%:',rtmp,rtmp/time_loop*100
        endif

      endif

      call mpi_allreduce(time_neigh1,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Nay-1 time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Nay-1 time/%:',rtmp,rtmp/time_loop*100
      endif

      if (nmolecular == 1) then
        call mpi_allreduce(time_neigh2,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Nay-2 time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Nay-2 time/%:',rtmp,rtmp/time_loop*100
        endif
      endif
      
      call mpi_allreduce(time_exch,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Exch  time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Exch  time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_comm,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Comm  time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Comm  time/%:',rtmp,rtmp/time_loop*100
      endif

      if (newton >= 1) then
        call mpi_allreduce(time_fcomm,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Fcomm time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Fcomm time/%:',rtmp,rtmp/time_loop*100
        endif
      endif

      if (nshake == 1) then
        call mpi_allreduce(time_shake,rtmp,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs
        if (node.eq.0) then
          write (6,900) 'Shake time/%:',rtmp,rtmp/time_loop*100
          write (1,900) 'Shake time/%:',rtmp,rtmp/time_loop*100
        endif
      endif

      call mpi_allreduce(time_io,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'I/O   time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'I/O   time/%:',rtmp,rtmp/time_loop*100
      endif

      call mpi_allreduce(time_other,rtmp,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      rtmp = rtmp/nprocs
      if (node.eq.0) then
        write (6,900) 'Other time/%:',rtmp,rtmp/time_loop*100
        write (1,900) 'Other time/%:',rtmp,rtmp/time_loop*100
      endif
      
      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      call stats(time_nonbond,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nbond time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nbond time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (coulstyle == 3 .or. coulstyle == 4) then
        call stats(time_long,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Long  time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Long  time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif
      
      if (nmolecular == 1) then

        call stats(time_bond,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Bond  time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Bond  time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        call stats(time_angle,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Angle time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Angle time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        call stats(time_dihedral,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Dihed time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Dihed time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        call stats(time_improper,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Impro time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Impro time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

      endif

      call stats(time_neigh1,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nay-1 time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nay-1 time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (nmolecular == 1) then
        call stats(time_neigh2,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Nay-2 time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Nay-2 time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      call stats(time_exch,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Exch  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Exch  time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_comm,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Comm  time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Comm  time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (newton >= 1) then
        call stats(time_fcomm,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Fcomm time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Fcomm time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      if (nshake == 1) then
        call stats(time_shake,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Shake time:',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Shake time:',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif
      
      call stats(time_io,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'I/O   time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'I/O   time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      call stats(time_other,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Other time:',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Other time:',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

c timing breakdown for Ewald

      if (coulstyle == 3) then

        if (node.eq.0) then
          write (6,*)
          write (6,*) 'Ewald timing info:'
          write (1,*)
          write (1,*) 'Ewald timing info:'
        endif

        do i = 1,kcount
          sfacrl(i) = 0.0
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,10
          call mpi_allreduce(sfacrl,sfacrl_all,kcount,
     $         mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        call mpi_allreduce(time_long,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        ave = rtmp/nprocs
        if (ave.eq.0.0) ave = 1.0
        rtmp = (time2-time1)/10 * 2

        if (node.eq.0) then
          write (6,*) '  Allreduce time per timestep =',rtmp
          write (6,*) '  Total allreduce time =',rtmp*nsteps
          write (6,*) '  Allreduce % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  Allreduce time per timestep =',rtmp
          write (1,*) '  Total allreduce time =',rtmp*nsteps
          write (1,*) '  Allreduce % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        if (ensemble >= 3) then

          call mpi_barrier(mpi_comm_world,ierror)
          time1 = mpi_wtime()
          do i = 1,10
            call ewald_coeff
          enddo
          call mpi_barrier(mpi_comm_world,ierror)
          time2 = mpi_wtime()

          call mpi_allreduce(time_other,rtmp,1,
     $         mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
          ave = rtmp/nprocs
          if (ave.eq.0.0) ave = 1.0
          rtmp = (time2-time1)/10

          if (node.eq.0) then
            write (6,*) '  Setup time per timestep =',rtmp
            write (6,*) '  Total setup time =',rtmp*nsteps
            write (6,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
            write (1,*) '  Setup time per timestep =',rtmp
            write (1,*) '  Total setup time =',rtmp*nsteps
            write (1,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
          endif

        endif

      endif

c timing breakdown for PPPM

      if (coulstyle == 4) then

        if (node.eq.0) then
          write (6,*)
          write (6,*) 'PPPM timing info:'
          write (1,*)
          write (1,*) 'PPPM timing info:'
        endif

        call mpi_allreduce(time_long,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        ave = rtmp/nprocs
        if (ave.eq.0.0) ave = 1.0

        call mpi_allreduce(time_rho,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs

        if (node.eq.0) then
          write (6,*) '  Make_rho time =',rtmp
          write (6,*) '  Make_rho % of long time =',
     $         rtmp/ave*100.0
          write (1,*) '  Make_rho time =',rtmp
          write (1,*) '  Make_rho % of long time =',
     $         rtmp/ave*100.0
        endif

        call mpi_allreduce(time_poiss,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs

        if (node.eq.0) then
          write (6,*) '  Poisson time =',rtmp
          write (6,*) '  Poisson % of long time =',
     $         rtmp/ave*100.0
          write (1,*) '  Poisson time =',rtmp
          write (1,*) '  Poisson % of long time =',
     $         rtmp/ave*100.0
        endif

        call mpi_allreduce(time_field,rtmp,1,
     $       mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        rtmp = rtmp/nprocs

        if (node.eq.0) then
          write (6,*) '  Electric_field time =',rtmp
          write (6,*) '  Electric_field % of long time =',
     $         rtmp/ave*100.0
          write (1,*) '  Electric_field time =',rtmp
          write (1,*) '  Electric_field % of long time =',
     $         rtmp/ave*100.0
        endif

        do indx = 1,(nxhi_out-nxlo_out+1)*(nyhi_out-nylo_out+1)*
     $       (nzhi_out-nzlo_out+1)
          density_brick(indx) = 0.0
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,10
          call brick2fft(density_brick,density_fft)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        rtmp = (time2-time1)/10

        if (node.eq.0) then
          write (6,*) '  Brick2fft time per timestep =',rtmp
          write (6,*) '  Total Brick2fft time =',rtmp*nsteps
          write (6,*) '  Brick2fft % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  Brick2fft time per timestep =',rtmp
          write (1,*) '  Total brick2fft time =',rtmp*nsteps
          write (1,*) '  Brick2fft % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        do indx = 1,(nxhi_out-nxlo_out+1)*(nyhi_out-nylo_out+1)*
     $       (nzhi_out-nzlo_out+1)
          vdx_brick(indx) = 0.0
          vdy_brick(indx) = 0.0
          vdz_brick(indx) = 0.0
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,10
          call fillbrick(vdx_brick,vdy_brick,vdz_brick)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        rtmp = (time2-time1)/10

        if (node.eq.0) then
          write (6,*) '  Fillbrick time per timestep =',rtmp
          write (6,*) '  Total fillbrick time =',rtmp*nsteps
          write (6,*) '  Fillbrick % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  Fillbrick time per timestep =',rtmp
          write (1,*) '  Total fillbrick time =',rtmp*nsteps
          write (1,*) '  Fillbrick % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        do i = 1,nfft
          workvec1(i) = cmplx(0.0,0.0)
        enddo

        call mpi_barrier(mpi_comm_world,ierror)
        time1 = mpi_wtime()
        do i = 1,5
          call fft_3d(workvec1,workvec1,1,plan1_fft)
          call fft_3d(workvec1,workvec1,-1,plan2_fft)
          call fft_3d(workvec1,workvec1,-1,plan2_fft)
          call fft_3d(workvec1,workvec1,-1,plan2_fft)
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        time2 = mpi_wtime()

        rtmp = (time2-time1)/5

        if (node.eq.0) then
          write (6,*) '  FFT time per timestep =',rtmp
          write (6,*) '  Total FFT time =',rtmp*nsteps
          write (6,*) '  FFT % of long time =',
     $         rtmp*nsteps/ave*100.0
          write (1,*) '  FFT time per timestep =',rtmp
          write (1,*) '  Total FFT time =',rtmp*nsteps
          write (1,*) '  FFT % of long time =',
     $         rtmp*nsteps/ave*100.0
        endif

        if (ensemble >= 3) then

          call mpi_barrier(mpi_comm_world,ierror)
          time1 = mpi_wtime()
          do i = 1,10
            call pppm_coeff(1)
          enddo
          call mpi_barrier(mpi_comm_world,ierror)
          time2 = mpi_wtime()

          call mpi_allreduce(time_other,rtmp,1,
     $         mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
          ave = rtmp/nprocs
          if (ave.eq.0.0) ave = 1.0
          rtmp = (time2-time1)/10

          if (node.eq.0) then
            write (6,*) '  Setup time per timestep =',rtmp
            write (6,*) '  Total setup time =',rtmp*nsteps
            write (6,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
            write (1,*) '  Setup time per timestep =',rtmp
            write (1,*) '  Total setup time =',rtmp*nsteps
            write (1,*) '  Setup % of other time =',
     $           rtmp*nsteps/ave*100.0
          endif

        endif

      endif

c load-balance info

      if (node.eq.0) then
        write (6,*)
        write (1,*)
      endif

      r8tmp = nlocal
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nlocal:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nlocal:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      r8tmp = nghost
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Nghost:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Nghost:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (nbonds > 0) then
        r8tmp = nbondlocal
        call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Nbonds:    ',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Nbonds:    ',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      if (nangles > 0) then
        r8tmp = nanglelocal
        call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Nangle:    ',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Nangle:    ',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      if (ndihedrals > 0) then
        r8tmp = ndihedlocal
        call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Ndihed:    ',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Ndihed:    ',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      if (nimpropers > 0) then
        r8tmp = nimprolocal
        call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Nimpro:    ',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Nimpro:    ',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      if (nlocal > 0) then
        r8tmp = nnlast(nlocal)
      else
        r8tmp = 0.0
      endif
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
      if (node.eq.0) then
        write (6,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
        write (6,902) ' Histogram:',(ihisto(i),i=1,10)
        write (1,901) 'Neighs:    ',ave,' ave',xmax,' max',
     $       xmin,' min'
        write (1,902) ' Histogram:',(ihisto(i),i=1,10)
      endif

      if (nswap > 0) then
        r8tmp = nslast(nswap)
        call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (node.eq.0) then
          write (6,901) 'Nswaps:    ',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Nswaps:    ',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif

      if (nlocal > 0) then
        r8tmp = nnlast(nlocal)
      else
        r8tmp = 0.0
      endif
      call mpi_allreduce(r8tmp,all,1,mpi_double_precision,mpi_sum,
     $     mpi_comm_world,ierror)

      itmp2 = 0
      do i = 1,nlocal
        itmp2 = itmp2 + num3bond(i)
      enddo

      r8tmp = itmp2
      call stats(r8tmp,1,ave,xmax,xmin,ihisto,ihistotmp,10)

      call mpi_allreduce(itmp2,itmp,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      itmp2 = itmp

      if (nmolecular == 1) then
        if (node.eq.0) then
          write (6,901) 'Nspecs:    ',ave,' ave',xmax,' max',xmin,' min'
          write (6,902) ' Histogram:',(ihisto(i),i=1,10)
          write (1,901) 'Nspecs:    ',ave,' ave',xmax,' max',
     $         xmin,' min'
          write (1,902) ' Histogram:',(ihisto(i),i=1,10)
        endif
      endif
      
      if (node.eq.0) then
        write (6,*)
        write (6,*) 'Total # of neighbors =',nint(all)
        write (6,903) 'Ave neighs/atom =',all/natoms
        write (6,903) 'Ave nspecs/atom =',float(itmp2)/natoms
        write (6,*) 'Number of reneighborings =',numneigh
        write (6,*) 'Dangerous reneighborings =',ndanger
        write (6,*)
        write (1,*)
        write (1,*) 'Total # of neighbors =',nint(all)
        write (1,903) 'Ave neighs/atom =',all/natoms
        write (1,903) 'Ave nspecs/atom =',float(itmp2)/natoms
        write (1,*) 'Number of reneighborings =',numneigh
        write (1,*) 'Dangerous reneighborings =',ndanger
        write (1,*)
      endif

      itmp = max_nlocal
      call mpi_allreduce(itmp,max_nlocal,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_bond
      call mpi_allreduce(itmp,max_bond,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_angle
      call mpi_allreduce(itmp,max_angle,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_dihed
      call mpi_allreduce(itmp,max_dihed,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_impro
      call mpi_allreduce(itmp,max_impro,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_neigh
      call mpi_allreduce(itmp,max_neigh,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_exch
      call mpi_allreduce(itmp,max_exch,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_slist
      call mpi_allreduce(itmp,max_slist,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      itmp = max_exch
      call mpi_allreduce(itmp,max_exch,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)
      if (neighstyle == 1) then
        itmp = mbinx*mbiny*mbinz
        call mpi_allreduce(itmp,max_bin,1,mpi_integer,mpi_max,
     $       mpi_comm_world,ierror)
      endif

      if (node == 0) then
        
        opt_own = extra_own * max_nlocal/maxown + 0.01
        opt_ghost = extra_ghost *
     $       max(max_slist,max_nghost)/maxghost + 0.01
        opt_neigh = extra_neigh * max_neigh/maxneigh + 0.01
        opt_buf = extra_buf * max_buf/maxbuf + 0.01

        if (maxown == minown) opt_own = extra_own
        if (maxghost == minghost) opt_ghost = extra_ghost
        if (maxneigh == minneigh) opt_neigh = extra_neigh
        if (maxbuf == minbuf) opt_buf = extra_buf

        write (1,*) 'Actual extra params: Extra own,ghost,neigh,buf'
        write (1,904) extra_own,extra_ghost,extra_neigh,extra_buf
        write (1,*) 'Optimal extra params: Extra own,ghost,neigh,buf'
        write (1,904) opt_own,opt_ghost,opt_neigh,opt_buf
        write (1,*) 'Max # of local atoms =',max_nlocal,
     $       ' out of',maxown
        write (1,*) 'Max # of ghost atoms =',max_nghost,
     $       ' out of',maxghost
        write (1,*) 'Max in neighbor list =',max_neigh,
     $       ' out of',maxneigh
        write (1,*) 'Max in swap list =',max_slist,
     $       ' out of',maxghost
        write (1,*) 'Max atoms exchanged =',max_exch
        write (1,*) 'Max atoms in border =',max_bord,
     $       ' out of',maxbuf/3
        write (1,*) 'Max use of comm buffers =',max_buf,
     $       ' out of',maxbuf
        write (1,*) 'Max # of bonds/atom =',maxbondper
        write (1,*) 'Max # of angles/atom =',maxangleper
        write (1,*) 'Max # of dihedrals/atom =',maxdihedper
        write (1,*) 'Max # of impropers/atom =',maximproper
        write (1,*) 'Max in bond list =',max_bond,
     $       ' out of',maxbondlocal
        write (1,*) 'Max in angle list =',max_angle,
     $       ' out of',maxanglelocal
        write (1,*) 'Max in dihedral list =',max_dihed,
     $       ' out of',maxdihedlocal
        write (1,*) 'Max in improper list =',max_impro,
     $       ' out of',maximprolocal
        if (neighstyle == 1) then
          write (1,*) 'Max # of neighbor bins =',max_bin
          write (1,*) 'Bins in stencil =',nstencil
        endif
        write (1,*) '# of swaps =',nswap,
     $       ' Needs =',need(1),need(2),need(3)
        write (1,*) 'Cutneigh =',cutneigh
        write (1,*) 'Cut/Box =',cutneigh*pgrid(1)/xprd,
     $       cutneigh*pgrid(2)/yprd,cutneigh*pgrid(3)/zprd
      endif

c end-of-run clean-up
c corresponds to reverse-order of how allocation was done in start

c deallocate Ewald & PPPM memory

      if (coulstyle == 3) call ewald_deallocate
      if (coulstyle == 4) call pppm_deallocate

c free up PPPM FFT and remap plans

      if (coulstyle == 4) then
        call fft_3d_destroy_plan(plan1_fft)
        call fft_3d_destroy_plan(plan2_fft)
        call remap_3d_destroy_plan(plan_remap)
      endif

c free neighbor bins and communication arrays

      if (neighstyle == 1) then
        deallocate(binpnt)
        maxbin = 0
      endif

      if (maxswap > 0) then
        deallocate(slablo,slabhi,spart,rpart)
        deallocate(nsfirst,nslast,nrfirst,nrlast)
        deallocate(commflag,commflagall)
        maxswap = 0
      endif

c free SHAKE neighbors and special lists of 1-2, 1-3, 1-4 neighbors

      if (nshake == 1) call shake_deallocate
      deallocate(specbond)

c free up rRESPA arrays

      if (nrespa == 1)
     $     deallocate(f_stretch,f_intra,f_short,f_long)

c free up fhold array

      if (nrespa == 0 .and. newton == 3) deallocate(fhold)

c reset special neighs after Ewald/PPPM

      if (coulstyle == 3 .or. coulstyle == 4) then
        special(1) = special(1) - 2.0
        special(2) = special(2) - 2.0
        special(3) = special(3) - 2.0
      endif

c free up nonbond arrays

      call nonbond_deallocate

      return
      end
