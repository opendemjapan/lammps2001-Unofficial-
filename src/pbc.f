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
c enforce PBC on appropriate dims
c order of 2 tests in each dim is important
c   to insure lo-bound <= coord < hi-bound
c   even with round-off errors where (coord +/- epsilon) +/- period = bound

      subroutine pbc
      use global
      implicit none

      integer i

      if (perflagx == 0) then
        do i = 1,nlocal
          if (x(1,i) < box(1,1)) then
            x(1,i) = x(1,i) + xprd
            true(i) = true(i) - 1
          endif
          if (x(1,i) >= box(2,1)) then
            x(1,i) = x(1,i) - xprd
            true(i) = true(i) + 1
          endif
        enddo
      endif

      if (perflagy == 0) then
        do i = 1,nlocal
          if (x(2,i) < box(1,2)) then
            x(2,i) = x(2,i) + yprd
            true(i) = true(i) - 1000
          endif
          if (x(2,i) >= box(2,2)) then
            x(2,i) = x(2,i) - yprd
            true(i) = true(i) + 1000
          endif
        enddo
      endif

      if (perflagz == 0) then
        do i = 1,nlocal
          if (x(3,i) < box(1,3)) then
            x(3,i) = x(3,i) + zprd
            true(i) = true(i) - 1000000
          endif
          if (x(3,i) >= box(2,3)) then
            x(3,i) = x(3,i) - zprd
            true(i) = true(i) + 1000000
          endif
        enddo
      endif

      return
      end
