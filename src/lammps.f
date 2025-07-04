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
c main program

      program lammps
      use global
      implicit none

c local variables

      integer iopflag

      call initialize
      
      do
        
        call input(iopflag)

        if (iopflag == 1) then
          call read_data
        else if (iopflag == 2) then
          call read_restart
        else if (iopflag == 3) then
          call temp_create
        else if (iopflag == 4) then
          call assign_fix
        else if (iopflag == 5) then
          call start
          if (nrespa == 0) call integrate
          if (nrespa /= 0) call integrate_respa
          call finish
        else if (iopflag == 6) then
          call start
          call minimize
          call finish
        endif
        
      enddo
      
      end
