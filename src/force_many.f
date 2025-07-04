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
c many-body forces
c  force components stored for atoms I own
c  variable force = (-1/r * d(phi)/dr)
c  variable virial = (-r * d(phi)/dr)

c angle - harmonic spring

      subroutine angle_harmonic(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,itype,ifactor
      real*8 small,delx1,dely1,delz1,rsq1,r1,delx2,dely2
      real*8 delz2,rsq2,r2,c,s,dtheta,tk,a,a11,a12,a22
      real*8 vx1,vx2,vy1,vy2,vz1,vz2
      e_angle = 0.0
      small = 0.001
      
      do m = 1,nanglelocal
        
        i1 = anglelist(1,m)
        i2 = anglelist(2,m)
        i3 = anglelist(3,m)
        itype = anglelist(4,m)

c 1st bond
        
        delx1 = x(1,i1) - x(1,i2)
        dely1 = x(2,i1) - x(2,i2)
        delz1 = x(3,i1) - x(3,i2)
        call minimg(delx1,dely1,delz1)
        
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1
        r1 = sqrt(rsq1)
        
c 2nd bond
        
        delx2 = x(1,i3) - x(1,i2)
        dely2 = x(2,i3) - x(2,i2)
        delz2 = x(3,i3) - x(3,i2)
        call minimg(delx2,dely2,delz2)
        
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2
        r2 = sqrt(rsq2)
        
c angle (cos and sin)
        
        c = delx1*delx2 + dely1*dely2 + delz1*delz2
        c = c / (r1*r2)
        
        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
        s = sqrt(1.0 - c*c)
        if (s.lt.small) s = small
        s = 1.0/s
        
c energy
        
        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 3
        endif

        dtheta = acos(c) - anglecoeff(2,itype)
        tk = anglecoeff(1,itype) * dtheta
        e_angle = e_angle + ifactor/3.0 * tk*dtheta
        
c force on each of 3 atoms
        
        a = 2.0 * tk * s
        
        a11 = a*c / rsq1
        a12 = -a / (r1*r2)
        a22 = a*c / rsq2
        
        vx1 = a11*delx1 + a12*delx2
        vx2 = a22*delx2 + a12*delx1
        vy1 = a11*dely1 + a12*dely2
        vy2 = a22*dely2 + a12*dely1
        vz1 = a11*delz1 + a12*delz2
        vz2 = a22*delz2 + a12*delz1
        
        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - vx1
          f(2,i1) = f(2,i1) - vy1
          f(3,i1) = f(3,i1) - vz1
        endif

        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (vx2+vx1)
          f(2,i2) = f(2,i2) + (vy2+vy1)
          f(3,i2) = f(3,i2) + (vz2+vz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) - vx2
          f(2,i3) = f(2,i3) - vy2
          f(3,i3) = f(3,i3) - vz2
        endif

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor/3.0*(delx1*vx1 + delx2*vx2)
          virial(2) = virial(2) - ifactor/3.0*(dely1*vy1 + dely2*vy2)
          virial(3) = virial(3) - ifactor/3.0*(delz1*vz1 + delz2*vz2)
          virial(4) = virial(4) - ifactor/3.0*(delx1*vy1 + delx2*vy2)
          virial(5) = virial(5) - ifactor/3.0*(delx1*vz1 + delx2*vz2)
          virial(6) = virial(6) - ifactor/3.0*(dely1*vz1 + dely2*vz2)
        endif

      enddo
      
      return
      end


c -------------------------------------------------------------------------
c angle - harmonic spring + Urey-Bradley term

      subroutine angle_charmm(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,itype,ifactor
      real*8 small,delx1,dely1,delz1,rsq1,r1,delx2,dely2
      real*8 delxUB,delyUB,delzUB,force,rsqUB,rUB,dr,rk
      real*8 delz2,rsq2,r2,c,s,dtheta,tk,a,a11,a12,a22
      real*8 vx1,vx2,vy1,vy2,vz1,vz2
      e_angle = 0.0
      small = 0.001
      
      do m = 1,nanglelocal
        
        i1 = anglelist(1,m)
        i2 = anglelist(2,m)
        i3 = anglelist(3,m)
        itype = anglelist(4,m)

c 1st bond
        
        delx1 = x(1,i1) - x(1,i2)
        dely1 = x(2,i1) - x(2,i2)
        delz1 = x(3,i1) - x(3,i2)
        call minimg(delx1,dely1,delz1)
        
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1
        r1 = sqrt(rsq1)
        
c 2nd bond
        
        delx2 = x(1,i3) - x(1,i2)
        dely2 = x(2,i3) - x(2,i2)
        delz2 = x(3,i3) - x(3,i2)
        call minimg(delx2,dely2,delz2)
        
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2
        r2 = sqrt(rsq2)

c Urey-Bradley bond
        
        delxUB = x(1,i3) - x(1,i1)
        delyUB = x(2,i3) - x(2,i1)
        delzUB = x(3,i3) - x(3,i1)
        call minimg(delxUB,delyUB,delzUB)
        
        rsqUB = delxUB*delxUB + delyUB*delyUB + delzUB*delzUB
        rUB = sqrt(rsqUB)
        
c angle (cos and sin)
        
        c = delx1*delx2 + dely1*dely2 + delz1*delz2
        c = c / (r1*r2)
        
        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
        s = sqrt(1.0 - c*c)
        if (s.lt.small) s = small
        s = 1.0/s
        
c energy
        
        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 3
        endif

        dtheta = acos(c) - anglecoeff(2,itype)
        tk = anglecoeff(1,itype) * dtheta
        e_angle = e_angle + ifactor/3.0 * tk*dtheta

c UB portion

        dr = rUB - anglecoeff(4,itype)
        rk = anglecoeff(3,itype) * dr
        e_angle = e_angle + ifactor/3.0 * rk*dr
        
c harmonic force on each of 3 atoms
        
        a = 2.0 * tk * s
        
        a11 = a*c / rsq1
        a12 = -a / (r1*r2)
        a22 = a*c / rsq2
        
        vx1 = a11*delx1 + a12*delx2
        vx2 = a22*delx2 + a12*delx1
        vy1 = a11*dely1 + a12*dely2
        vy2 = a22*dely2 + a12*dely1
        vz1 = a11*delz1 + a12*delz2
        vz2 = a22*delz2 + a12*delz1

c UB force

        if (rUB.gt.0.0) then
          force = -2.0*rk/rUB
        else
          force = 0.0
        endif
        
c total force on each of 3 atoms

        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - vx1 - delxUB*force
          f(2,i1) = f(2,i1) - vy1 - delyUB*force
          f(3,i1) = f(3,i1) - vz1 - delzUB*force
        endif
        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (vx2+vx1)
          f(2,i2) = f(2,i2) + (vy2+vy1)
          f(3,i2) = f(3,i2) + (vz2+vz1)
        endif
        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) - vx2 + delxUB*force
          f(2,i3) = f(2,i3) - vy2 + delyUB*force
          f(3,i3) = f(3,i3) - vz2 + delzUB*force
        endif

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor/3.0*(delx1*vx1 + delx2*vx2 -
     $         delxUB*delxUB*force)
          virial(2) = virial(2) - ifactor/3.0*(dely1*vy1 + dely2*vy2 -
     $         delyUB*delyUB*force)
          virial(3) = virial(3) - ifactor/3.0*(delz1*vz1 + delz2*vz2 -
     $         delzUB*delzUB*force)
          virial(4) = virial(4) - ifactor/3.0*(delx1*vy1 + delx2*vy2 -
     $         delxUB*delyUB*force)
          virial(5) = virial(5) - ifactor/3.0*(delx1*vz1 + delx2*vz2 -
     $         delxUB*delzUB*force)
          virial(6) = virial(6) - ifactor/3.0*(dely1*vz1 + dely2*vz2 -
     $         delyUB*delzUB*force)
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c angle - cosine = K (1 + cos(theta))

      subroutine angle_cosine(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,itype,ifactor
      real*8 small,delx1,dely1,delz1,rsq1,r1,delx2,dely2
      real*8 delz2,rsq2,r2,c,s,dtheta,tk,a,a11,a12,a22
      real*8 vx1,vx2,vy1,vy2,vz1,vz2
      e_angle = 0.0
      small = 0.001
      
      do m = 1,nanglelocal
        
        i1 = anglelist(1,m)
        i2 = anglelist(2,m)
        i3 = anglelist(3,m)
        itype = anglelist(4,m)

c 1st bond
        
        delx1 = x(1,i1) - x(1,i2)
        dely1 = x(2,i1) - x(2,i2)
        delz1 = x(3,i1) - x(3,i2)
        call minimg(delx1,dely1,delz1)
        
        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1
        r1 = sqrt(rsq1)
        
c 2nd bond
        
        delx2 = x(1,i3) - x(1,i2)
        dely2 = x(2,i3) - x(2,i2)
        delz2 = x(3,i3) - x(3,i2)
        call minimg(delx2,dely2,delz2)
        
        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2
        r2 = sqrt(rsq2)
        
c angle (cos and sin)
        
        c = delx1*delx2 + dely1*dely2 + delz1*delz2
        c = c / (r1*r2)
        if (c > 1.0) c = 1.0
        if (c < -1.0) c = -1.0
        
c energy
        
        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 3
        endif

        e_angle = e_angle + ifactor/3.0 * anglecoeff(1,itype)*(1.0+c)
        
c force on each of 3 atoms
        
        a = -anglecoeff(1,itype)
        
        a11 = a*c / rsq1
        a12 = -a / (r1*r2)
        a22 = a*c / rsq2
        
        vx1 = a11*delx1 + a12*delx2
        vx2 = a22*delx2 + a12*delx1
        vy1 = a11*dely1 + a12*dely2
        vy2 = a22*dely2 + a12*dely1
        vz1 = a11*delz1 + a12*delz2
        vz2 = a22*delz2 + a12*delz1
        
        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - vx1
          f(2,i1) = f(2,i1) - vy1
          f(3,i1) = f(3,i1) - vz1
        endif

        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (vx2+vx1)
          f(2,i2) = f(2,i2) + (vy2+vy1)
          f(3,i2) = f(3,i2) + (vz2+vz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) - vx2
          f(2,i3) = f(2,i3) - vy2
          f(3,i3) = f(3,i3) - vz2
        endif

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor/3.0*(delx1*vx1 + delx2*vx2)
          virial(2) = virial(2) - ifactor/3.0*(dely1*vy1 + dely2*vy2)
          virial(3) = virial(3) - ifactor/3.0*(delz1*vz1 + delz2*vz2)
          virial(4) = virial(4) - ifactor/3.0*(delx1*vy1 + delx2*vy2)
          virial(5) = virial(5) - ifactor/3.0*(delx1*vz1 + delx2*vz2)
          virial(6) = virial(6) - ifactor/3.0*(dely1*vz1 + dely2*vz2)
        endif

      enddo
      
      return
      end


c -------------------------------------------------------------------------
c single simple-harmonic dihedral torsion

      subroutine dihedral_harmonic(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,i4,itype,n,ifactor
      real*8 tolerance,small,vb1x,vb1y,vb1z,vb2x,vb2y
      real*8 vb2z,vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z,sb1
      real*8 sb2,sb3,rb1,rb2,rb3,c0,b1mag2,b1mag,b2mag2
      real*8 b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      real*8 c2mag,sc1,sc2,s1,s12,c,p,pd,rc2,a,a11,a22
      real*8 a33,a12,a13,a23,sx1,sx2,sx12,sy1,sy2,sy12
      real*8 sz1,sz2,sz12,s2

      e_dihedral = 0.0
      tolerance = 0.05
      small = 0.001
      
      do m = 1,ndihedlocal
        
        i1 = dihedlist(1,m)
        i2 = dihedlist(2,m)
        i3 = dihedlist(3,m)
        i4 = dihedlist(4,m)
        itype = dihedlist(5,m)
        
c 1st bond in dihedral
        
        vb1x = x(1,i1) - x(1,i2)
        vb1y = x(2,i1) - x(2,i2)
        vb1z = x(3,i1) - x(3,i2)
        call minimg(vb1x,vb1y,vb1z)
        
c 2nd bond in dihedral
        
        vb2x = x(1,i3) - x(1,i2)
        vb2y = x(2,i3) - x(2,i2)
        vb2z = x(3,i3) - x(3,i2)
        call minimg(vb2x,vb2y,vb2z)

        vb2xm = -vb2x
        vb2ym = -vb2y
        vb2zm = -vb2z
        call minimg(vb2xm,vb2ym,vb2zm)

c 3rd bond in dihedral
        
        vb3x = x(1,i4) - x(1,i3)
        vb3y = x(2,i4) - x(2,i3)
        vb3z = x(3,i4) - x(3,i3)
        call minimg(vb3x,vb3y,vb3z)
        
c c0 calculation
        
        sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
        sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
        sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
        
        rb1 = sqrt(sb1)
        rb2 = sqrt(sb2)
        rb3 = sqrt(sb3)
        
        c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
        
c 1st and 2nd angle
        
        b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
        b1mag = sqrt(b1mag2)
        b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
        b2mag = sqrt(b2mag2)
        b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
        b3mag = sqrt(b3mag2)

        ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
        r12c1 = 1.0 / (b1mag*b2mag)
        c1mag = ctmp * r12c1

        ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z
        r12c2 = 1.0 / (b2mag*b3mag)
        c2mag = ctmp * r12c2
        
c cos and sin of 2 angles and final c
        
        sc1 = sqrt(1.0 - c1mag*c1mag)
        if (sc1 < small) sc1 = small
        sc1 = 1.0/sc1
        
        sc2 = sqrt(1.0 - c2mag*c2mag)
        if (sc2 < small) sc2 = small
        sc2 = 1.0/sc2
        
        s1 = sc1 * sc1
        s2 = sc2 * sc2
        s12 = sc1 * sc2
        c = (c0 + c1mag*c2mag) * s12
        
c error check

        if (c > 1.0+tolerance.or.c < (-1.0-tolerance)) then
          write (6,*) 'Dihedral:',node,ntimestep,
     $         tag(i1),tag(i2),tag(i3),tag(i4)
          write (6,*) ' Error1:',node,x(1,i1),x(2,i1),x(3,i1)
          write (6,*) ' Error2:',node,x(1,i2),x(2,i2),x(3,i2)
          write (6,*) ' Error3:',node,x(1,i3),x(2,i3),x(3,i3)
          write (6,*) ' Error4:',node,x(1,i4),x(2,i4),x(3,i4)
        endif

        if (c > 1.0) c = 1.0
        if (c < -1.0) c = -1.0
        
c energy
c  p = 1 + cos(n*phi) for d = 1
c  p = 1 - cos(n*phi) for d = -1
c  pd = dp/dc / 2
        
        n = nint(dihedcoeff(3,itype))

        if (n.eq.2) then
          p = 2.0*c*c
          pd = 2.0*c
        else if (n.eq.3) then
          rc2 = c*c
          p = (4.0*rc2-3.0)*c + 1.0
          pd = 6.0*rc2 - 1.5
        else if (n.eq.4) then
          rc2 = c*c
          p = 8.0*(rc2-1)*rc2 + 2.0
          pd = (16.0*rc2-8.0)*c
        else if (n.eq.6) then
          rc2 = c*c
          p = ((32.0*rc2-48.0)*rc2 + 18.0)*rc2
          pd = (96.0*(rc2-1.0)*rc2 + 18.0)*c
        else if (n.eq.1) then
          p = c + 1.0
          pd = 0.5
        else if (n.eq.5) then
          rc2 = c*c
          p = ((16.0*rc2-20.0)*rc2 + 5.0)*c + 1.0
          pd = (40.0*rc2-30.0)*rc2 + 2.5
        else if (n.eq.0) then
          p = 2.0
          pd = 0.0
        endif

        if (nint(dihedcoeff(2,itype)).eq.-1) then
          p = 2.0 - p
          pd = -pd
        endif

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
          if (i4.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 4
        endif

        e_dihedral = e_dihedral + ifactor*0.25 *
     $       dihedcoeff(1,itype)*p

c force on all 4 atoms
        
        a = 2.0 * dihedcoeff(1,itype) * pd
        c = c * a
        s12 = s12 * a
        a11 = (-c*sb1*s1)
        a22 = sb2*(2.0*c0*s12 - c*(s1+s2))
        a33 = (-c*sb3*s2)
        a12 = r12c1*(c1mag*c*s1 + c2mag*s12)
        a13 = rb1*rb3*s12
        a23 = r12c2*(-c2mag*c*s2 - c1mag*s12)
        
        sx1  = a11*vb1x + a12*vb2x + a13*vb3x
        sx2  = a12*vb1x + a22*vb2x + a23*vb3x
        sx12 = a13*vb1x + a23*vb2x + a33*vb3x
        sy1  = a11*vb1y + a12*vb2y + a13*vb3y
        sy2  = a12*vb1y + a22*vb2y + a23*vb3y
        sy12 = a13*vb1y + a23*vb2y + a33*vb3y
        sz1  = a11*vb1z + a12*vb2z + a13*vb3z
        sz2  = a12*vb1z + a22*vb2z + a23*vb3z
        sz12 = a13*vb1z + a23*vb2z + a33*vb3z

        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - sx1
          f(2,i1) = f(2,i1) - sy1
          f(3,i1) = f(3,i1) - sz1
        endif

        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (sx2 + sx1)
          f(2,i2) = f(2,i2) + (sy2 + sy1)
          f(3,i2) = f(3,i2) + (sz2 + sz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) + sx12 - sx2
          f(2,i3) = f(2,i3) + sy12 - sy2
          f(3,i3) = f(3,i3) + sz12 - sz2
        endif

        if (newton_bond.eq.1.or.i4.le.nlocal) then
          f(1,i4) = f(1,i4) - sx12
          f(2,i4) = f(2,i4) - sy12
          f(3,i4) = f(3,i4) - sz12
        endif

c virial contribution: Sum r_ij * F_ij

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor*0.25 *
     $         (vb1x*sx1 + vb2x*sx2 + vb3x*sx12)
          virial(2) = virial(2) - ifactor*0.25 *
     $         (vb1y*sy1 + vb2y*sy2 + vb3y*sy12)
          virial(3) = virial(3) - ifactor*0.25 *
     $         (vb1z*sz1 + vb2z*sz2 + vb3z*sz12)
          virial(4) = virial(4) - ifactor*0.25 *
     $         (vb1x*sy1 + vb2x*sy2 + vb3x*sy12)
          virial(5) = virial(5) - ifactor*0.25 *
     $         (vb1x*sz1 + vb2x*sz2 + vb3x*sz12)
          virial(6) = virial(6) - ifactor*0.25 *
     $         (vb1y*sz1 + vb2y*sz2 + vb3y*sz12)
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c CHARMM dihedral torsion, includes LJ & coul 1-4 terms

      subroutine dihedral_charmm(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,i4,itype,n,ifactor
      integer itype1, itype4
      real*8 rsq14,forcelj,force,r2inv,r6inv,delx,dely,delz
      real*8 philj,phicoul,forcecoul
      real*8 tolerance,small,vb1x,vb1y,vb1z,vb2x,vb2y
      real*8 vb2z,vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z,sb1
      real*8 sb2,sb3,rb1,rb2,rb3,c0,b1mag2,b1mag,b2mag2
      real*8 b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      real*8 c2mag,sc1,sc2,s1,s12,c,p,pd,rc2,a,a11,a22
      real*8 a33,a12,a13,a23,sx1,sx2,sx12,sy1,sy2,sy12
      real*8 sz1,sz2,sz12,s2
      real*8 tmp

      e_dihedral = 0.0
      e_14_coul = 0.0
      e_14_vdwl = 0.0
      tolerance = 0.05
      small = 0.001
      
      do m = 1,ndihedlocal
        
        i1 = dihedlist(1,m)
        i2 = dihedlist(2,m)
        i3 = dihedlist(3,m)
        i4 = dihedlist(4,m)
        itype = dihedlist(5,m)
        
c 1st bond in dihedral
        
        vb1x = x(1,i1) - x(1,i2)
        vb1y = x(2,i1) - x(2,i2)
        vb1z = x(3,i1) - x(3,i2)
        call minimg(vb1x,vb1y,vb1z)
        
c 2nd bond in dihedral
        
        vb2x = x(1,i3) - x(1,i2)
        vb2y = x(2,i3) - x(2,i2)
        vb2z = x(3,i3) - x(3,i2)
        call minimg(vb2x,vb2y,vb2z)

        vb2xm = -vb2x
        vb2ym = -vb2y
        vb2zm = -vb2z
        call minimg(vb2xm,vb2ym,vb2zm)

c 3rd bond in dihedral
        
        vb3x = x(1,i4) - x(1,i3)
        vb3y = x(2,i4) - x(2,i3)
        vb3z = x(3,i4) - x(3,i3)
        call minimg(vb3x,vb3y,vb3z)

c 1-4 interaction
        
        delx = x(1,i1) - x(1,i4)
        dely = x(2,i1) - x(2,i4)
        delz = x(3,i1) - x(3,i4)
        call minimg(delx,dely,delz)
        rsq14 = delx*delx + dely*dely + delz*delz
        itype1 = type(i1)
        itype4 = type(i4)
        
c c0 calculation
        
        sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
        sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
        sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
        
        rb1 = sqrt(sb1)
        rb2 = sqrt(sb2)
        rb3 = sqrt(sb3)
        
        c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
        
c 1st and 2nd angle
        
        b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
        b1mag = sqrt(b1mag2)
        b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
        b2mag = sqrt(b2mag2)
        b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
        b3mag = sqrt(b3mag2)

        ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
        r12c1 = 1.0 / (b1mag*b2mag)
        c1mag = ctmp * r12c1

        ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z
        r12c2 = 1.0 / (b2mag*b3mag)
        c2mag = ctmp * r12c2
        
c cos and sin of 2 angles and final c
        
        sc1 = sqrt(1.0 - c1mag*c1mag)
        if (sc1.lt.small) sc1 = small
        sc1 = 1.0/sc1
        
        sc2 = sqrt(1.0 - c2mag*c2mag)
        if (sc2.lt.small) sc2 = small
        sc2 = 1.0/sc2
        
        s1 = sc1 * sc1
        s2 = sc2 * sc2
        s12 = sc1 * sc2
        c = (c0 + c1mag*c2mag) * s12
        
c error check

        if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
          write (6,*) 'Dihedral:',node,ntimestep,
     $         tag(i1),tag(i2),tag(i3),tag(i4)
          write (6,*) ' Error1:',node,x(1,i1),x(2,i1),x(3,i1)
          write (6,*) ' Error2:',node,x(1,i2),x(2,i2),x(3,i2)
          write (6,*) ' Error3:',node,x(1,i3),x(2,i3),x(3,i3)
          write (6,*) ' Error4:',node,x(1,i4),x(2,i4),x(3,i4)
        endif

        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
c energy
c  p = 1 + cos(n*phi) for d = 1
c  p = 1 - cos(n*phi) for d = -1
c  pd = dp/dc / 2
        
        n = nint(dihedcoeff(2,itype))

        if (n.eq.2) then
          p = 2.0*c*c
          pd = 2.0*c
        else if (n.eq.3) then
          rc2 = c*c
          p = (4.0*rc2-3.0)*c + 1.0
          pd = 6.0*rc2 - 1.5
        else if (n.eq.4) then
          rc2 = c*c
          p = 8.0*(rc2-1)*rc2 + 2.0
          pd = (16.0*rc2-8.0)*c
        else if (n.eq.6) then
          rc2 = c*c
          p = ((32.0*rc2-48.0)*rc2 + 18.0)*rc2
          pd = (96.0*(rc2-1.0)*rc2 + 18.0)*c
        else if (n.eq.1) then
          p = c + 1.0
          pd = 0.5
        else if (n.eq.5) then
          rc2 = c*c
          p = ((16.0*rc2-20.0)*rc2 + 5.0)*c + 1.0
          pd = (40.0*rc2-30.0)*rc2 + 2.5
        else if (n.eq.0) then
          p = 2.0
          pd = 0.0
        endif

        if (dihedcoeff(3,itype).eq.180.) then
          p = 2.0 - p
          pd = -pd
        endif

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
          if (i4.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 4
        endif

        e_dihedral = e_dihedral + ifactor/4.0 *
     $       dihedcoeff(1,itype)*p

c force on all 4 atoms
        
        a = 2.0 * dihedcoeff(1,itype) * pd
        c = c * a
        s12 = s12 * a
        a11 = (-c*sb1*s1)
        a22 = sb2*(2.0*c0*s12 - c*(s1+s2))
        a33 = (-c*sb3*s2)
        a12 = r12c1*(c1mag*c*s1 + c2mag*s12)
        a13 = rb1*rb3*s12
        a23 = r12c2*(-c2mag*c*s2 - c1mag*s12)
        
        sx1  = a11*vb1x + a12*vb2x + a13*vb3x
        sx2  = a12*vb1x + a22*vb2x + a23*vb3x
        sx12 = a13*vb1x + a23*vb2x + a33*vb3x
        sy1  = a11*vb1y + a12*vb2y + a13*vb3y
        sy2  = a12*vb1y + a22*vb2y + a23*vb3y
        sy12 = a13*vb1y + a23*vb2y + a33*vb3y
        sz1  = a11*vb1z + a12*vb2z + a13*vb3z
        sz2  = a12*vb1z + a22*vb2z + a23*vb3z
        sz12 = a13*vb1z + a23*vb2z + a33*vb3z

        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - sx1
          f(2,i1) = f(2,i1) - sy1
          f(3,i1) = f(3,i1) - sz1
        endif

        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (sx2 + sx1)
          f(2,i2) = f(2,i2) + (sy2 + sy1)
          f(3,i2) = f(3,i2) + (sz2 + sz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) + sx12 - sx2
          f(2,i3) = f(2,i3) + sy12 - sy2
          f(3,i3) = f(3,i3) + sz12 - sz2
        endif

        if (newton_bond.eq.1.or.i4.le.nlocal) then
          f(1,i4) = f(1,i4) - sx12
          f(2,i4) = f(2,i4) - sy12
          f(3,i4) = f(3,i4) - sz12
        endif

c virial contribution: Sum r_ij * F_ij

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor/4.0 *
     $         (vb1x*sx1 + vb2x*sx2 + vb3x*sx12)
          virial(2) = virial(2) - ifactor/4.0 *
     $         (vb1y*sy1 + vb2y*sy2 + vb3y*sy12)
          virial(3) = virial(3) - ifactor/4.0 *
     $         (vb1z*sz1 + vb2z*sz2 + vb3z*sz12)
          virial(4) = virial(4) - ifactor/4.0 *
     $         (vb1x*sy1 + vb2x*sy2 + vb3x*sy12)
          virial(5) = virial(5) - ifactor/4.0 *
     $         (vb1x*sz1 + vb2x*sz2 + vb3x*sz12)
          virial(6) = virial(6) - ifactor/4.0 *
     $         (vb1y*sz1 + vb2y*sz2 + vb3y*sz12)
        endif

c 1-4 LJ & Coulomb interactions

c Warning: are assuming here that all nonbond 1-4 interactions (generally 
c   less that 6 angstroms apart) are within the inner cutoffs,
c   cutljinterior and cutcoulint (generally at least 8 angstroms)

        if (dihedcoeff(4,itype) /= 0) then

          r2inv = 1.0/rsq14
          r6inv = r2inv*r2inv*r2inv

          forcecoul = coulpre*q(i1)*q(i4)*sqrt(r2inv)
          forcelj = r6inv*(lj14_1(itype1,itype4)*r6inv-
     $                     lj14_2(itype1,itype4))
          force = dihedcoeff(4,itype)*(forcelj + forcecoul)*r2inv

c compute energy if this is thermodynamic timestep

          if (iflag >= 1) then  
            phicoul = forcecoul
            phicoul = ifactor/4.0 * phicoul * dihedcoeff(4,itype)

            philj = r6inv *
     $         (lj14_3(itype1,itype4)*r6inv - lj14_4(itype1,itype4))
            philj = ifactor/4.0 * philj * dihedcoeff(4,itype)

            e_14_coul = e_14_coul + phicoul
            e_14_vdwl = e_14_vdwl + philj
          endif

c add nonbond force contribution to total on each atom

          if (newton_nonbond == 1 .or. i1 <= nlocal) then
            f(1,i1) = f(1,i1) + delx*force
            f(2,i1) = f(2,i1) + dely*force
            f(3,i1) = f(3,i1) + delz*force
          endif
          if (newton_nonbond == 1 .or. i4 <= nlocal) then
            f(1,i4) = f(1,i4) - delx*force
            f(2,i4) = f(2,i4) - dely*force
            f(3,i4) = f(3,i4) - delz*force
          endif

          if (iflag >= 1) then
            force = force*ifactor/4.0
            virial(1) = virial(1) + delx*delx*force
            virial(2) = virial(2) + dely*dely*force
            virial(3) = virial(3) + delz*delz*force
            virial(4) = virial(4) + delx*dely*force
            virial(5) = virial(5) + delx*delz*force
            virial(6) = virial(6) + dely*delz*force
          endif

        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c multiple simple-harmonic dihedral torsions

      subroutine dihedral_multiharmonic(iflag)
      use global
      implicit none

c arguemnt variables

      integer iflag

c local variables

      integer m,i1,i2,i3,i4,itype,ifactor
      real*8 tolerance,small,vb1x,vb1y,vb1z,vb2x,vb2y
      real*8 vb2z,vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z,sb1
      real*8 sb2,sb3,rb1,rb2,rb3,c0,b1mag2,b1mag,b2mag2
      real*8 b2mag,b3mag2,b3mag,ctmp,r12c1,c1mag,r12c2
      real*8 c2mag,sc1,sc2,s1,s12,c,p,pd,a,a11,a22
      real*8 a33,a12,a13,a23,sx1,sx2,sx12,sy1,sy2,sy12
      real*8 sz1,sz2,sz12,s2

      e_dihedral = 0.0
      tolerance = 0.05
      small = 0.001

      do m = 1,ndihedlocal
        
        i1 = dihedlist(1,m)
        i2 = dihedlist(2,m)
        i3 = dihedlist(3,m)
        i4 = dihedlist(4,m)
        itype = dihedlist(5,m)
        
c 1st bond in dihedral
        
        vb1x = x(1,i1) - x(1,i2)
        vb1y = x(2,i1) - x(2,i2)
        vb1z = x(3,i1) - x(3,i2)
        call minimg(vb1x,vb1y,vb1z)
        
c 2nd bond in dihedral
        
        vb2x = x(1,i3) - x(1,i2)
        vb2y = x(2,i3) - x(2,i2)
        vb2z = x(3,i3) - x(3,i2)
        call minimg(vb2x,vb2y,vb2z)

c 3rd bond in dihedral
        
        vb3x = x(1,i4) - x(1,i3)
        vb3y = x(2,i4) - x(2,i3)
        vb3z = x(3,i4) - x(3,i3)
        call minimg(vb3x,vb3y,vb3z)
        
c c0 calculation
        
        sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
        sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)
        sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)
        
        rb1 = sqrt(sb1)
        rb2 = sqrt(sb2)
        rb3 = sqrt(sb3)
        
        c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3
        
c 1st and 2nd angle
        
        b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
        b1mag = sqrt(b1mag2)
        b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
        b2mag = sqrt(b2mag2)
        b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
        b3mag = sqrt(b3mag2)

        ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
        r12c1 = 1.0 / (b1mag*b2mag)
        c1mag = ctmp * r12c1

        ctmp = vb2x*vb3x + vb2y*vb3y + vb2z*vb3z
        r12c2 = 1.0 / (b2mag*b3mag)
        c2mag = - ctmp * r12c2
        
c cos and sin of 2 angles and final c
        
        sc1 = sqrt(1.0 - c1mag*c1mag)
        if (sc1 < small) sc1 = small
        sc1 = 1.0/sc1
        
        sc2 = sqrt(1.0 - c2mag*c2mag)
        if (sc2 < small) sc2 = small
        sc2 = 1.0/sc2
        
        s1 = sc1 * sc1
        s2 = sc2 * sc2
        s12 = sc1 * sc2
        c = (c0 + c1mag*c2mag) * s12
        
c error check

        if (c > 1.0+tolerance.or.c < (-1.0-tolerance)) then
          write (6,*) 'Dihedral:',node,ntimestep,
     $         tag(i1),tag(i2),tag(i3),tag(i4)
          write (6,*) ' Error1:',node,x(1,i1),x(2,i1),x(3,i1)
          write (6,*) ' Error2:',node,x(1,i2),x(2,i2),x(3,i2)
          write (6,*) ' Error3:',node,x(1,i3),x(2,i3),x(3,i3)
          write (6,*) ' Error4:',node,x(1,i4),x(2,i4),x(3,i4)
        endif

        if (c > 1.0) c = 1.0
        if (c < -1.0) c = -1.0
        
c energy
c  p = sum (i=1,5) a_i * c**(i-1)
c  pd = dp/dc
        
        p = dihedcoeff(1,itype) +
     $       c*(dihedcoeff(2,itype) +
     $       c*(dihedcoeff(3,itype) +
     $       c*(dihedcoeff(4,itype) +
     $       c*dihedcoeff(5,itype))))

        pd = dihedcoeff(2,itype) +
     $       c*(2.0*dihedcoeff(3,itype) +
     $       c*(3.0*dihedcoeff(4,itype) +
     $       c*4.0*dihedcoeff(5,itype)))

        if (newton_bond.eq.0) then
           ifactor = 0
           if (i1.le.nlocal) ifactor = ifactor + 1
           if (i2.le.nlocal) ifactor = ifactor + 1
           if (i3.le.nlocal) ifactor = ifactor + 1
           if (i4.le.nlocal) ifactor = ifactor + 1
        else
           ifactor = 4
        endif
        
        e_dihedral = e_dihedral + ifactor*0.25 * p

c force on all 4 atoms
        
        a = pd
        c = c * a
        s12 = s12 * a
        a11 = (-c*sb1*s1)
        a22 = sb2*(2.0*c0*s12 - c*(s1+s2))
        a33 = (-c*sb3*s2)
        a12 = r12c1*(c1mag*c*s1 + c2mag*s12)
        a13 = rb1*rb3*s12
        a23 = r12c2*(-c2mag*c*s2 - c1mag*s12)
        
        sx1  = a11*vb1x + a12*vb2x + a13*vb3x
        sx2  = a12*vb1x + a22*vb2x + a23*vb3x
        sx12 = a13*vb1x + a23*vb2x + a33*vb3x
        sy1  = a11*vb1y + a12*vb2y + a13*vb3y
        sy2  = a12*vb1y + a22*vb2y + a23*vb3y
        sy12 = a13*vb1y + a23*vb2y + a33*vb3y
        sz1  = a11*vb1z + a12*vb2z + a13*vb3z
        sz2  = a12*vb1z + a22*vb2z + a23*vb3z
        sz12 = a13*vb1z + a23*vb2z + a33*vb3z

        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - sx1
          f(2,i1) = f(2,i1) - sy1
          f(3,i1) = f(3,i1) - sz1
        endif

        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (sx2 + sx1)
          f(2,i2) = f(2,i2) + (sy2 + sy1)
          f(3,i2) = f(3,i2) + (sz2 + sz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) + sx12 - sx2
          f(2,i3) = f(2,i3) + sy12 - sy2
          f(3,i3) = f(3,i3) + sz12 - sz2
        endif

        if (newton_bond.eq.1.or.i4.le.nlocal) then
          f(1,i4) = f(1,i4) - sx12
          f(2,i4) = f(2,i4) - sy12
          f(3,i4) = f(3,i4) - sz12
        endif

c virial contribution: Sum r_ij * F_ij

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor*0.25 *
     $         (vb1x*sx1 + vb2x*sx2 + vb3x*sx12)
          virial(2) = virial(2) - ifactor*0.25 *
     $         (vb1y*sy1 + vb2y*sy2 + vb3y*sy12)
          virial(3) = virial(3) - ifactor*0.25 *
     $         (vb1z*sz1 + vb2z*sz2 + vb3z*sz12)
          virial(4) = virial(4) - ifactor*0.25 *
     $         (vb1x*sy1 + vb2x*sy2 + vb3x*sy12)
          virial(5) = virial(5) - ifactor*0.25 *
     $         (vb1x*sz1 + vb2x*sz2 + vb3x*sz12)
          virial(6) = virial(6) - ifactor*0.25 *
     $         (vb1y*sz1 + vb2y*sz2 + vb3y*sz12)
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c harmonic out-of-plane improper

      subroutine improper_harmonic(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,i4,itype,ifactor
      real*8 tolerance,small,v1x,v1y,v1z,v2x,v2y,v2z,v3x
      real*8 v3y,v3z,ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2
      real*8 s12,c,s,domega,a,a11,a22,a33,a12,a13,a23,sx1
      real*8 sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12

      e_improper = 0.0
      tolerance = 0.05
      small = 0.001
      
      do m = 1,nimprolocal
        
        i1 = improlist(1,m)
        i2 = improlist(2,m)
        i3 = improlist(3,m)
        i4 = improlist(4,m)
        itype = improlist(5,m)
        
c geometry of 4-body
        
        v1x = x(1,i2) - x(1,i1)
        v1y = x(2,i2) - x(2,i1)
        v1z = x(3,i2) - x(3,i1)
        call minimg(v1x,v1y,v1z)
        
        v2x = x(1,i3) - x(1,i2)
        v2y = x(2,i3) - x(2,i2)
        v2z = x(3,i3) - x(3,i2)
        call minimg(v2x,v2y,v2z)
        
        v3x = x(1,i4) - x(1,i3)
        v3y = x(2,i4) - x(2,i3)
        v3z = x(3,i4) - x(3,i3)
        call minimg(v3x,v3y,v3z)
        
        ss1 = 1.0 / (v1x*v1x + v1y*v1y + v1z*v1z)
        ss2 = 1.0 / (v2x*v2x + v2y*v2y + v2z*v2z)
        ss3 = 1.0 / (v3x*v3x + v3y*v3y + v3z*v3z)
        
        r1 = sqrt(ss1)
        r2 = sqrt(ss2)
        r3 = sqrt(ss3)
        
c sin and cos of angle
        
        c0 = -(v1x * v3x + v1y * v3y + v1z * v3z) * r1 * r3
        c1 = -(v1x * v2x + v1y * v2y + v1z * v2z) * r1 * r2
        c2 = -(v3x * v2x + v3y * v2y + v3z * v2z) * r3 * r2
        
        s1 = 1.0 - c1*c1
        if (s1.lt.small) s1 = small
        s1 = 1.0 / s1
        
        s2 = 1.0 - c2*c2
        if (s2.lt.small) s2 = small
        s2 = 1.0 / s2
        
        s12 = sqrt(s1*s2)
        c = (c1*c2 + c0) * s12
        
c error check
        
        if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
          write (6,*) 'Improper:',node,ntimestep,
     $         tag(i1),tag(i2),tag(i3),tag(i4)
          write (6,*) ' Error1:',node,x(1,i1),x(2,i1),x(3,i1)
          write (6,*) ' Error2:',node,x(1,i2),x(2,i2),x(3,i2)
          write (6,*) ' Error3:',node,x(1,i3),x(2,i3),x(3,i3)
          write (6,*) ' Error4:',node,x(1,i4),x(2,i4),x(3,i4)
        endif

        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0
        
        s = sqrt(1.0 - c*c)
        if (s.lt.small) s = small
        
c energy
        
        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
          if (i4.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 4
        endif

        domega = acos(c) - improcoeff(2,itype)
        a = improcoeff(1,itype) * domega

        e_improper = e_improper + ifactor*0.25 * a*domega
        
c force on all 4 atoms
        
        a = -a * 2.0/s
        c = c * a
        
        s12 = s12 * a
        a11 = (-c * ss1 * s1)
        a22 = ss2 * (2.0 * c0 * s12 - c * (s1 + s2))
        a33 = (-c * ss3 * s2)
        a12 = r1 * r2 * (c1 * c * s1 + c2 * s12)
        a13 = r1 * r3 * s12
        a23 = r2 * r3 * (-c2 * c * s2 - c1 * s12)
        
        sx1  = a12*v2x + a13*v3x - a11*v1x
        sx2  = a22*v2x + a23*v3x - a12*v1x
        sx12 = a23*v2x + a33*v3x - a13*v1x
        sy1  = a12*v2y + a13*v3y - a11*v1y
        sy2  = a22*v2y + a23*v3y - a12*v1y
        sy12 = a23*v2y + a33*v3y - a13*v1y
        sz1  = a12*v2z + a13*v3z - a11*v1z
        sz2  = a22*v2z + a23*v3z - a12*v1z
        sz12 = a23*v2z + a33*v3z - a13*v1z
        
        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - sx1
          f(2,i1) = f(2,i1) - sy1
          f(3,i1) = f(3,i1) - sz1
        endif
        
        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (sx2 + sx1)
          f(2,i2) = f(2,i2) + (sy2 + sy1)
          f(3,i2) = f(3,i2) + (sz2 + sz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) + sx12 - sx2
          f(2,i3) = f(2,i3) + sy12 - sy2
          f(3,i3) = f(3,i3) + sz12 - sz2
        endif

        if (newton_bond.eq.1.or.i4.le.nlocal) then
          f(1,i4) = f(1,i4) - sx12
          f(2,i4) = f(2,i4) - sy12
          f(3,i4) = f(3,i4) - sz12
        endif

c virial contribution: Sum r_ij * F_ij

        if (iflag >= 1) then
          virial(1) = virial(1) + ifactor*0.25 *
     $         (v1x*sx1 - v2x*sx2 - v3x*sx12)
          virial(2) = virial(2) + ifactor*0.25 *
     $         (v1y*sy1 - v2y*sy2 - v3y*sy12)
          virial(3) = virial(3) + ifactor*0.25 *
     $         (v1z*sz1 - v2z*sz2 - v3z*sz12)
          virial(4) = virial(4) + ifactor*0.25 *
     $         (v1x*sy1 - v2x*sy2 - v3x*sy12)
          virial(5) = virial(5) + ifactor*0.25 *
     $         (v1x*sz1 - v2x*sz2 - v3x*sz12)
          virial(6) = virial(6) + ifactor*0.25 *
     $         (v1y*sz1 - v2y*sz2 - v3y*sz12)
        endif

      enddo

      return
      end


c -------------------------------------------------------------------------
c cvff out-of-plane improper

      subroutine improper_cvff(iflag)
      use global
      implicit none

c argument variables

      integer iflag

c local variables

      integer m,i1,i2,i3,i4,itype,n,ifactor
      real*8 tolerance,small,vb1x,vb1y,vb1z,vb2x,vb2y,vb2z
      real*8 vb2xm,vb2ym,vb2zm,vb3x,vb3y,vb3z,sb1,sb2,sb3
      real*8 rb1,rb2,rb3,c0,b1mag2,b1mag,b2mag2,b2mag,b3mag2
      real*8 b3mag,ctmp,r12c1,c1mag,r12c2,c2mag,sc1,sc2,s1,s2
      real*8 s12,c,p,pd,rc2,a,a11,a22,a33,a12,a13,a23,sx1
      real*8 sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12

      e_improper = 0.0
      tolerance = 0.05
      small = 0.001

      do m = 1,nimprolocal

        i1 = improlist(1,m)
        i2 = improlist(2,m)
        i3 = improlist(3,m)
        i4 = improlist(4,m)
        itype = improlist(5,m)

c 1st bond in improper

        vb1x = x(1,i1) - x(1,i2)
        vb1y = x(2,i1) - x(2,i2)
        vb1z = x(3,i1) - x(3,i2)
        call minimg(vb1x,vb1y,vb1z)

c 2nd bond in improper

        vb2x = x(1,i3) - x(1,i2)
        vb2y = x(2,i3) - x(2,i2)
        vb2z = x(3,i3) - x(3,i2)
        call minimg(vb2x,vb2y,vb2z)

        vb2xm = -vb2x
        vb2ym = -vb2y
        vb2zm = -vb2z
        call minimg(vb2xm,vb2ym,vb2zm)

c 3rd bond in improper

        vb3x = x(1,i4) - x(1,i3)
        vb3y = x(2,i4) - x(2,i3)
        vb3z = x(3,i4) - x(3,i3)
        call minimg(vb3x,vb3y,vb3z)

c c0 calculation

        sb1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z)
        sb2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z)

        sb3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z)

        rb1 = sqrt(sb1)
        rb2 = sqrt(sb2)
        rb3 = sqrt(sb3)

        c0 = (vb1x*vb3x + vb1y*vb3y + vb1z*vb3z) * rb1*rb3

c 1st and 2nd angle

        b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z
        b1mag = sqrt(b1mag2)
        b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z
        b2mag = sqrt(b2mag2)
        b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z
        b3mag = sqrt(b3mag2)

        ctmp = vb1x*vb2x + vb1y*vb2y + vb1z*vb2z
        r12c1 = 1.0 / (b1mag*b2mag)
        c1mag = ctmp * r12c1

        ctmp = vb2xm*vb3x + vb2ym*vb3y + vb2zm*vb3z
        r12c2 = 1.0 / (b2mag*b3mag)
        c2mag = ctmp * r12c2

c cos and sin of 2 angles and final c

        sc1 = sqrt(1.0 - c1mag*c1mag)
        if (sc1.lt.small) sc1 = small
        sc1 = 1.0/sc1

        sc2 = sqrt(1.0 - c2mag*c2mag)
        if (sc2.lt.small) sc2 = small
        sc2 = 1.0/sc2

        s1 = sc1 * sc1
        s2 = sc2 * sc2
        s12 = sc1 * sc2
        c = (c0 + c1mag*c2mag) * s12

c error check

        if (c.gt.1.0+tolerance.or.c.lt.(-1.0-tolerance)) then
          write (6,*) 'Improper:',node,ntimestep,
     $         tag(i1),tag(i2),tag(i3),tag(i4)
          write (6,*) ' Error1:',node,x(1,i1),x(2,i1),x(3,i1)
          write (6,*) ' Error2:',node,x(1,i2),x(2,i2),x(3,i2)
          write (6,*) ' Error3:',node,x(1,i3),x(2,i3),x(3,i3)
          write (6,*) ' Error4:',node,x(1,i4),x(2,i4),x(3,i4)
        endif

        if (c.gt.1.0) c = 1.0
        if (c.lt.-1.0) c = -1.0

c energy
c  p = 1 + cos(n*phi) for d = 1
c  p = 1 - cos(n*phi) for d = -1
c  pd = dp/dc / 2

        n = nint(improcoeff(3,itype))

        if (n.eq.2) then
          p = 2.0*c*c
          pd = 2.0*c
        else if (n.eq.3) then
          rc2 = c*c
          p = (4.0*rc2-3.0)*c + 1.0
          pd = 6.0*rc2 - 1.5
        else if (n.eq.4) then
          rc2 = c*c
          p = 8.0*(rc2-1)*rc2 + 2.0
          pd = (16.0*rc2-8.0)*c
        else if (n.eq.6) then
          rc2 = c*c
          p = ((32.0*rc2-48.0)*rc2 + 18.0)*rc2
          pd = (96.0*(rc2-1.0)*rc2 + 18.0)*c
        else if (n.eq.1) then
          p = c + 1.0
          pd = 0.5
        else if (n.eq.5) then
          rc2 = c*c
          p = ((16.0*rc2-20.0)*rc2 + 5.0)*c + 1.0
          pd = (40.0*rc2-30.0)*rc2 + 2.5
        else if (n.eq.0) then
          p = 2.0
          pd = 0.0
        endif

        if (nint(improcoeff(2,itype)).eq.-1) then
          p = 2.0 - p
          pd = -pd
        endif

        if (newton_bond.eq.0) then
          ifactor = 0
          if (i1.le.nlocal) ifactor = ifactor + 1
          if (i2.le.nlocal) ifactor = ifactor + 1
          if (i3.le.nlocal) ifactor = ifactor + 1
          if (i4.le.nlocal) ifactor = ifactor + 1
        else
          ifactor = 4
        endif

        e_improper = e_improper + ifactor*0.25 *
     $       improcoeff(1,itype)*p

c force on all 4 atoms

        a = 2.0 * improcoeff(1,itype) * pd
        c = c * a
        s12 = s12 * a
        a11 = (-c*sb1*s1)
        a22 = sb2*(2.0*c0*s12 - c*(s1+s2))
        a33 = (-c*sb3*s2)
        a12 = r12c1*(c1mag*c*s1 + c2mag*s12)
        a13 = rb1*rb3*s12
        a23 = r12c2*(-c2mag*c*s2 - c1mag*s12)

        sx1  = a11*vb1x + a12*vb2x + a13*vb3x
        sx2  = a12*vb1x + a22*vb2x + a23*vb3x
        sx12 = a13*vb1x + a23*vb2x + a33*vb3x
        sy1  = a11*vb1y + a12*vb2y + a13*vb3y
        sy2  = a12*vb1y + a22*vb2y + a23*vb3y
        sy12 = a13*vb1y + a23*vb2y + a33*vb3y
        sz1  = a11*vb1z + a12*vb2z + a13*vb3z
        sz2  = a12*vb1z + a22*vb2z + a23*vb3z
        sz12 = a13*vb1z + a23*vb2z + a33*vb3z

        if (newton_bond.eq.1.or.i1.le.nlocal) then
          f(1,i1) = f(1,i1) - sx1
          f(2,i1) = f(2,i1) - sy1
          f(3,i1) = f(3,i1) - sz1
        endif

        if (newton_bond.eq.1.or.i2.le.nlocal) then
          f(1,i2) = f(1,i2) + (sx2 + sx1)
          f(2,i2) = f(2,i2) + (sy2 + sy1)
          f(3,i2) = f(3,i2) + (sz2 + sz1)
        endif

        if (newton_bond.eq.1.or.i3.le.nlocal) then
          f(1,i3) = f(1,i3) + sx12 - sx2
          f(2,i3) = f(2,i3) + sy12 - sy2
          f(3,i3) = f(3,i3) + sz12 - sz2
        endif

        if (newton_bond.eq.1.or.i4.le.nlocal) then
          f(1,i4) = f(1,i4) - sx12
          f(2,i4) = f(2,i4) - sy12
          f(3,i4) = f(3,i4) - sz12
        endif

c virial contribution: Sum r_ij * F_ij

        if (iflag >= 1) then
          virial(1) = virial(1) - ifactor*0.25 *
     $         (vb1x*sx1 + vb2x*sx2 + vb3x*sx12)
          virial(2) = virial(2) - ifactor*0.25 *
     $         (vb1y*sy1 + vb2y*sy2 + vb3y*sy12)
          virial(3) = virial(3) - ifactor*0.25 *
     $         (vb1z*sz1 + vb2z*sz2 + vb3z*sz12)
          virial(4) = virial(4) - ifactor*0.25 *
     $         (vb1x*sy1 + vb2x*sy2 + vb3x*sy12)
          virial(5) = virial(5) - ifactor*0.25 *
     $         (vb1x*sz1 + vb2x*sz2 + vb3x*sz12)
          virial(6) = virial(6) - ifactor*0.25 *
     $         (vb1y*sz1 + vb2y*sz2 + vb3y*sz12)
        endif

      enddo

      return
      end
