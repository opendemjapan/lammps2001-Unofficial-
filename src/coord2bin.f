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
c convert xyz atom coords into local bin #
c take special care to insure ghost atoms with
c  coord >= boxhi or coord < boxlo are put in correct bins

      integer function coord2bin(xcoord,ycoord,zcoord)
      use global
      implicit none
      real*8 xcoord,ycoord,zcoord

      integer ix,iy,iz

      if (xcoord.ge.box(2,1)) then
        ix = int((xcoord-box(2,1))*bininvx) + nbinx - mbinxlo
      else if (xcoord.ge.box(1,1)) then
        ix = int((xcoord-box(1,1))*bininvx) - mbinxlo
      else
        ix = int((xcoord-box(1,1))*bininvx) - mbinxlo - 1
      endif

      if (ycoord.ge.box(2,2)) then
        iy = int((ycoord-box(2,2))*bininvy) + nbiny - mbinylo
      else if (ycoord.ge.box(1,2)) then
        iy = int((ycoord-box(1,2))*bininvy) - mbinylo
      else
        iy = int((ycoord-box(1,2))*bininvy) - mbinylo - 1
      endif

      if (zcoord.ge.box(2,3)) then
        iz = int((zcoord-box(2,3))*bininvz) + nbinz - mbinzlo
      else if (zcoord.ge.box(1,3)) then
        iz = int((zcoord-box(1,3))*bininvz) - mbinzlo
      else
        iz = int((zcoord-box(1,3))*bininvz) - mbinzlo - 1
      endif

      coord2bin = iz*mbiny*mbinx + iy*mbinx + ix + 1
      
      return
      end
