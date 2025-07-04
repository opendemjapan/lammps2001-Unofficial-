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

c Contributing author: Todd Plantenga (Sandia)

c -------------------------------------------------------------------------
C   FILE min_support.f
C
C   Externally callable routines:
C      minimize
C      Compute_f_g
C      Get_Smallest_Box_Dim
C      Get_Atom_Positions
C      Put_Atom_Positions
C      Compute_Inf_Norm_PAR
C      Compute_L2_Norm_PAR
C      Compute_ddot_PAR
C      Normalize_ddot
C      Find_Biggest_Jump
C      Repartitn_and_Reneighbor
C
C   Local routines (in order)
C      Write_Energy_Vals
C
C   Code development:
C      21 Oct 96 - Originated by T. Plantenga, Sandia National Labs.

c -------------------------------------------------------------------------
      SUBROUTINE minimize
      use global
      use mpi
      implicit none

C       orig:  T. Plantenga   4 Aug 97
C
C   This routine is called from "lammps.f" to find a molecular
C   configuration that gives a local minimum of the potential energy.
C
C   The following LAMMPS input command lines pass in relevant parameters:
C   "neighbor      real int int int int" --> skin, neighstyle, neighfreq,
C                                             neighdelay, neightrigger
C   "dump atom  int file"             --> ndumpatom
C   "dump force int file"             --> ndumpforce
C   "restart    int file"             --> nrestart
C   "min style  string"               --> optstyle (always 1 for hftn)
C   "min flag   int"                  --> freq of iteration output
C   "minimize   real int int"         --> opt_stop_tol, opt_max_iters,
C                                             opt_max_fns
C
C   The passed algorithm control parameters have the following effects:
C        neighfreq     = exchange & reneighbor after this many acc steps
C        neighdelay    = delay reneighboring until after this many acc steps
C        nrestart      = write restart file after this many evals (0=none)
C        opt_stop_tol  = applied to ||grad||_inf
C        opt_max_iters = maximum number of outer iterations
C        opt_max_fns   = maximum number of force/energy evaluations
C
C   At the end of execution, LAMMPS dumps some timing information.
C   For minimization, these mean:
C        time_loop      - total time in minimize
C        time_bond      - computing forces due to covalent bonds
C        time_angle     - computing forces due to angles
C        time_dihedral  - computing forces due to dihedral angles
C        time_improper  - computing forces due to improper angles
C        time_nonbond   - computing forces due to long-range interactions
C        time_fcomm     - comm for updating force calculations
C        time_comm      - comm for updating atom positions
C        time_exch      - comm for shifting atoms to new processors
C        time_neigh1    - computing new neighbor lists
C        time_neigh2    - verifying new neighbor lists
C
C**********************************************************************

C     ---- LOCAL VARIABLES.
      integer           i,iflag
      INTEGER           procNum
      INTEGER           totalN, localN
      INTEGER           num_fns
      DOUBLE PRECISION  endTime

C----------------------------------------------------------------------

c do not allow thermostatting fixes when doing a minimization

      iflag = 0
      do i = 1,nfixes
        if (fixstyle(i) == 4 .or. fixstyle(i) == 5 .or. 
     $       fixstyle(i) == 6) iflag = 1
      enddo
      if (iflag > 0) call err(
     $     'Minimiziation incompatible with temperature fix')

      procNum = node

c write minimizer header

      if (procnum == 0) then
        write (6,*)
        write (6,*) 'Starting Minimizer ...'
      endif

C check smoothness of nonbond terms

      if (coulstyle == 1) THEN
        if (procNum == 0) then
          write (6,*) 'WARNING: Coulomb potential is cutoff'
        endif
      endif

      if (nonstyle == 1 .or. nonstyle == 5) then
        if (procNum == 0) then
          write (6,*) 'WARNING: LJ potential is cutoff'
        endif
      endif

C     ---- DUMP THE INITIAL POTENTIAL ENERGY (COMPUTED IN "start.f").

      CALL Write_Energy_Vals (procNum)

C     ---- CREATE VARIABLES FOR THE MINIMIZATION ALGORITHM.

      totalN = 3 * natoms
      localN = 3 * nlocal

C     ---- CALL THE MINIMIZATION ALGORITHM.

      time_loop = mpi_wtime()
      IF (optstyle .EQ. 1) THEN
        CALL Min_HFTN (procNum, totalN, localN,
     $       opt_stop_tol, opt_max_fns, opt_max_iters,
     $       neighfreq, neighdelay, num_fns, optflag,
     $       ndumpatom, ndumpforce, maxown)
      ENDIF
      endTime = mpi_wtime()
      time_loop = endTime - time_loop

C     ---- DUMP THE FINAL POTENTIAL ENERGY.

      CALL Write_Energy_Vals (procNum)

c write final output files if requested

      if (ndumpatom > 0) call dump_atom(1)
      if (ndumpforce > 0) call dump_force(1)
      if (nrestart > 0) call write_restart(1)

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Compute_f_g(obj_f, g_vec_size, g_vec)
      use global
      use mpi
      implicit none
C
C       orig: 22 Oct 96  T. Plantenga
C       mod:   8 May 97  (TDP)  ignore constrained forces
C       mod:  27 Jul 97  (TDP)  adapted for LAMMPS V4.0
C
C   Input/output variable definitions:
C       obj_f            O  potential energy for minimization
C       g_vec_size      I
C       g_vec            O  gradient vector
C
      integer ierror
      real*8 tmp
      DOUBLE PRECISION  obj_f
      INTEGER           g_vec_size
      DOUBLE PRECISION  g_vec(g_vec_size)
C
C   This routine consists of calls extracted from "integrate" and
C   "thermo" that compute interatomic forces and potential energy.
C   Calculations are done at the current atom positions stored in
C   the global variable "x".
C
C   LAMMPS V4.0 combines the computation of energy and forces for
C   covalent bond terms, but not for the non-bonded interactions.
C   This subroutine should support all bonded and non-bonded model
C   "styles", but only the following have been tested:
C        harmonic (bond, angle, dihedral, improper)
C        lj/cutoff, lj/smooth
C        coulomb/cutoff, coul/smooth
C
C   Force computations are done in parallel, but without interprocessor
C   communication until the very end of this routine.
C   
C   External subroutines called:
C        zero_force_virial       "misc.f"
C        error                   "misc.f"
C        bond_harmonic           "force_bond.f"
C        bond_fene_standard      "force_bond.f"
C        bond_fene_shift         "force_bond.f"
C        bond_nonlinear          "force_bond.f"
C        bond_class2             "force_bond.f"
C        angle_harmonic          "force_many.f"
C        angle_charmm            "force_many.f"
C        angle_cosine            "force_many.f"
C        angle_class2            "force_class2.f"
C        dihedral_harmonic       "force_many.f"
C        dihedral_class2         "force_class2.f"
C        dihedral_multiharmonic  "force_many.f"
C        dihedral_charmm         "force_many.f"
C        improper_harmonic       "force_many.f"
C        improper_class2         "force_class2.f"
C        angleangle_class2       "force_class2.f"
C        lj_cut_coul_cut         "force.f"
C        lj_cut_coul_long        "force.f"
C        lj_cut                  "force.f"
C        lj_smooth               "force.f"
C        lj_smooth_coul_cut      "force.f"
C        lj_smooth_coul_long     "force.f"
C        lj_shift                "force.f"
C        lj_charmm_coul_long     "force.f"
C        lj_charmm_coul_charmm   "force.f"
C        soft                    "force.f"
C        ljclass2_cut            "force.f"
C        ljclass2_cut_coul_cut   "force.f"
C        ljclass2_cut_coul_long  "force.f"
C        lj_cut_coul_smooth      "force.f"
C        lj_smooth_coul_smooth   "force.f"
C        ewald                   "ewald.f"
C        pppm                    "pppm.f"
C        reverse_comm            "communicate.f"
C        fix_apply               "fix.f"
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER           i, j, k
      INTEGER           iflag
      DOUBLE PRECISION  time1, time2, time3, time4, time5
      DOUBLE PRECISION  e_potential, rnorm

C----------------------------------------------------------------------

C     ---- SET IFLAG TO ZERO TO DISABLE VIRIAL CALCULATIONS.
      iflag = 0

C     ---- ZERO THE GLOBAL FORCE VECTOR (SEE "integrate").
      CALL zero_force_virial ()

C----------------------------------------------------------------------
C   COMPUTE COVALENT BOND FORCES AND ENERGIES SIMULTANEOUSLY (TAKEN
C   FROM "integrate").  RESULTS ARE STORED AS SIDE EFFECTS IN THE
C   GLOBAL VARIABLES:
C        FORCES   - f(1..3,1..maxatom)
C        ENERGIES - e_bond
C                   e_angle
C                   e_dihedral
C                   e_improper
C----------------------------------------------------------------------

C     ---- THE barrier CALL GIVES MORE ACCURATE TIMING INFORMATION.

      if (nmolecular == 1) then

#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()

C     ---- BOND LENGTH TERMS.
        IF (nbonds .GT. 0) THEN
          IF (bondstyle .EQ. 1) THEN
            CALL bond_harmonic (iflag)
          ELSE IF (bondstyle .EQ. 2) THEN
            CALL bond_fene_standard (iflag)
          ELSE IF (bondstyle .EQ. 3) THEN
            CALL bond_fene_shift (iflag)
          ELSE IF (bondstyle .EQ. 4) THEN
            CALL bond_nonlinear (iflag)
          ELSE
            CALL bond_class2 (iflag)
          ENDIF
        ENDIF

        time2 = mpi_wtime()

C     ---- BOND ANGLE TERMS.
        IF (nangles .GT. 0) THEN
          IF (anglestyle .EQ. 1) THEN
            CALL angle_harmonic (iflag)
          ELSE IF (anglestyle .EQ. 2) THEN
            CALL angle_class2 (iflag)
          ELSE IF (anglestyle .EQ. 3) THEN
            CALL angle_charmm (iflag)
          ELSE IF (anglestyle .EQ. 4) THEN
            CALL angle_cosine (iflag)
          ENDIF
        ENDIF

        time3 = mpi_wtime()

C     ---- DIHEDRAL ANGLE TERMS.
        IF (ndihedrals .GT. 0) THEN
          IF (dihedstyle .EQ. 1) THEN
            CALL dihedral_harmonic (iflag)
          ELSE IF (dihedstyle .EQ. 2) THEN
            CALL dihedral_class2 (iflag)
          ELSE IF (dihedstyle .EQ. 3) THEN
            CALL dihedral_multiharmonic (iflag)
          ELSE IF (dihedstyle .EQ. 4) THEN
            CALL dihedral_charmm (iflag)
          ENDIF
        ENDIF

        time4 = mpi_wtime()

C     ---- IMPROPER ANGLE TERMS.
        IF (nimpropers .GT. 0) THEN
          IF (improstyle .EQ. 1) THEN
            CALL improper_harmonic (iflag)
          ELSE IF (improstyle .EQ. 2) then
            call improper_cvff(iflag)
          ELSE IF (improstyle .EQ. 3) THEN
            CALL improper_class2 (iflag)
            CALL angleangle_class2 (iflag)
          ENDIF
        ENDIF

        time5 = mpi_wtime()

        time_bond = time_bond + time2-time1
        time_angle = time_angle + time3-time2
        time_dihedral = time_dihedral + time4-time3
        time_improper = time_improper + time5-time4

      endif

C----------------------------------------------------------------------
C   COMPUTE NON-BOND FORCES (TAKEN FROM "integrate").  THESE ROUTINES
C   DO NOT SIMULTANEOUSLY CALCULATE ENERGY.  RESULTS APPEAR IN THE
C   GLOBAL VARIABLE "f(1..3,1..maxatom)".
C----------------------------------------------------------------------

#ifdef SYNC
      call mpi_barrier(mpi_comm_world,ierror)
#endif
      time1 = mpi_wtime()

C     ---- LENNARD-JONES (AKA VAN DER WAALS) FORCES.

      IF (nonstyle .EQ. 0) THEN
        IF (coulstyle .EQ. 1) THEN
          CALL lj_cut_coul_cut (iflag)
        ELSE IF (coulstyle .EQ. 2) THEN
          CALL lj_cut_coul_smooth (iflag)
        ELSE IF (coulstyle .EQ. 3 .or. coulstyle == 4) THEN
          CALL lj_cut_coul_long (iflag)
        ENDIF
      ELSE IF (nonstyle .EQ. 1) THEN
        IF (coulstyle .EQ. 0) THEN
          CALL lj_cut (iflag)
        ELSE IF (coulstyle .EQ. 1) THEN
          CALL lj_cut_coul_cut (iflag)
        ELSE IF (coulstyle .EQ. 2) THEN
          CALL lj_cut_coul_smooth (iflag)
        ELSE IF (coulstyle .EQ. 3 .or. coulstyle == 4) THEN
          CALL lj_cut_coul_long (iflag)
        ENDIF
      ELSE IF (nonstyle .EQ. 2) THEN
        IF (coulstyle .EQ. 0) THEN
          CALL lj_smooth (iflag)
        ELSE IF (coulstyle .EQ. 1) THEN
          CALL lj_smooth_coul_cut (iflag)
        ELSE IF (coulstyle .EQ. 2) THEN
          CALL lj_smooth_coul_smooth (iflag)
        ELSE IF (coulstyle .EQ. 3 .or. coulstyle == 4) THEN
          CALL lj_smooth_coul_long (iflag)
        ENDIF
      ELSE IF (nonstyle .EQ. 3) THEN
        CALL lj_shift (iflag)
      ELSE IF (nonstyle .EQ. 4) THEN
        CALL soft (iflag)
      ELSE IF (nonstyle .EQ. 5) THEN
        IF (coulstyle .EQ. 0) THEN
          CALL ljclass2_cut (iflag)
        ELSE IF (coulstyle .EQ. 1) THEN
          CALL ljclass2_cut_coul_cut (iflag)
        ELSE IF (coulstyle .EQ. 3 .or. coulstyle == 4) THEN
          CALL ljclass2_cut_coul_long (iflag)
        ENDIF
      ELSE IF (nonstyle .EQ. 6) THEN
        IF (coulstyle .EQ. 3 .or. coulstyle == 4) THEN
          CALL lj_charmm_coul_long(iflag)
        ELSE IF (coulstyle .EQ. 5) THEN
          CALL lj_charmm_coul_charmm(iflag)
        ENDIF
      ENDIF

      time2 = mpi_wtime()
      time_nonbond = time_nonbond + time2-time1

C     ---- LONG-RANGE COULOMB FORCES

      IF (coulstyle == 3 .OR. coulstyle == 4) THEN
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
        IF (coulstyle .EQ. 3) THEN
          CALL ewald (iflag)
        ELSE
          CALL pppm (iflag)
        ENDIF
        time2 = mpi_wtime()
        time_long = time_long + time2-time1
      endif

C     ---- SWAP FORCE INFORMATION BETWEEN PROCESSORS AND GET THE
C     ---- TOTAL FORCE ON EACH LOCAL ATOM.

      IF (newton .GE. 1) THEN
#ifdef SYNC
        call mpi_barrier(mpi_comm_world,ierror)
#endif
        time1 = mpi_wtime()
        CALL reverse_comm ()
        time2 = mpi_wtime()
        time_fcomm = time_fcomm + (time2 - time1)
      ENDIF

c     ---- CONSTRAINTS

      if (nfixes > 0) call fix_apply(1,f,0.0)

C     ---- COMPUTE NON-BONDED ENERGIES (SEE "energy" IN "thermo.f").
C     ---- RESULTS APPEAR IN THE GLOBAL VARIABLES "e_vdwl" AND "e_coul".

      CALL energy ()

C----------------------------------------------------------------------
C   EXTRACT INFORMATION FROM GLOBAL VARIABLES TO GET f AND g.
C----------------------------------------------------------------------

C     ---- SWAP AND ACCUMULATE EACH ENERGY COMPONENT (THESE CALLS
C     ---- ARE COPIED FROM SUBROUTINE "thermo").  AT THE END, EACH
C     ---- PROCESSOR KNOWS THE TOTAL ENERGIES.

      tmp = e_vdwl
      call mpi_allreduce(tmp,e_vdwl,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)
      if (ncharge == 1) then
        tmp = e_coul
        call mpi_allreduce(tmp,e_coul,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif

      if (nmolecular == 1) then
        tmp = e_bond
        call mpi_allreduce(tmp,e_bond,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_angle
        call mpi_allreduce(tmp,e_angle,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_dihedral
        call mpi_allreduce(tmp,e_dihedral,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_improper
        call mpi_allreduce(tmp,e_improper,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif

      if (anglestyle.eq.2) then
        tmp = e_bondbond
        call mpi_allreduce(tmp,e_bondbond,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_bondangle
        call mpi_allreduce(tmp,e_bondangle,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif
      if (dihedstyle.eq.2) then
        tmp = e_midbondtorsion
        call mpi_allreduce(tmp,e_midbondtorsion,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_endbondtorsion
        call mpi_allreduce(tmp,e_endbondtorsion,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_angletorsion
        call mpi_allreduce(tmp,e_angletorsion,1
     &       ,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_angleangletorsion
        call mpi_allreduce(tmp,e_angleangletorsion,1
     &       ,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
        tmp = e_bondbond13
        call mpi_allreduce(tmp,e_bondbond13,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif
      if (improstyle.eq.3) then
        tmp = e_angleangle
        call mpi_allreduce(tmp,e_angleangle,1,mpi_double_precision,
     $       mpi_sum,mpi_comm_world,ierror)
      endif

C     ---- LUMP IN CLASS 2 TERMS

      IF (thermostyle .LE. 1) THEN
        IF (anglestyle .EQ. 2)
     +       e_angle = e_angle + e_bondbond + e_bondangle
        IF (dihedstyle .EQ. 2)
     +       e_dihedral = e_dihedral + e_midbondtorsion +
     +       e_endbondtorsion + e_angletorsion +
     +       e_angleangletorsion
        IF (improstyle .EQ. 3)
     +       e_improper = e_improper + e_angleangle
      ENDIF

C     ---- ADD EVERYTHING TO GET THE POTENTIAL ENERGY.

      IF (thermostyle .LE. 1) THEN
        e_potential = e_long + e_vdwl + e_coul +
     +       e_bond + e_angle + e_dihedral + e_improper
      ELSE
        e_potential = e_long + e_vdwl + e_coul +
     +       e_bond + e_angle + e_dihedral + e_improper +
     +       e_bondbond + e_bondangle + e_midbondtorsion +
     +       e_endbondtorsion + e_angletorsion +
     +       e_angleangletorsion + e_angleangle
      ENDIF

C     ---- NORMALIZE ENERGIES TO THE NUMBER OF ATOMS IF UNITS ARE LJ

      IF (units .EQ. 1) THEN
        rnorm = 1.0D0 / (DBLE (natoms))
        e_long = e_long*rnorm
        e_vdwl = e_vdwl*rnorm
        e_coul = e_coul*rnorm
        e_bond = e_bond*rnorm
        e_angle = e_angle*rnorm
        e_dihedral = e_dihedral*rnorm
        e_improper = e_improper*rnorm
        IF (anglestyle .EQ. 2) THEN
          e_bondbond = e_bondbond*rnorm
          e_bondangle = e_bondangle*rnorm
        ENDIF
        IF (dihedstyle .EQ. 2) THEN
          e_midbondtorsion = e_midbondtorsion*rnorm
          e_endbondtorsion = e_endbondtorsion*rnorm
          e_angletorsion = e_angletorsion*rnorm
          e_angleangletorsion = e_angleangletorsion*rnorm
        ENDIF
        IF (improstyle .EQ. 3) e_angleangle = e_angleangle*rnorm
        e_potential = e_potential*rnorm
        e_total = e_total*rnorm
      ENDIF

C----------------------------------------------------------------------
C   COPY INFORMATION FROM GLOBAL VARIABLES.
C----------------------------------------------------------------------

      obj_f = e_potential

C     ---- COPY THE FORCE VECTOR.  BY CONVENTION, THE LAMMPS FORCE IS
C     ---- STEEPEST DESCENT, SO NEGATE IT INTO A GRADIENT.

      IF (g_vec_size .LT. (3 * nlocal)) THEN
        call err('g_vec_size in min_support must be larger')
      ENDIF

      j = 0
      DO i = 1, nlocal
        g_vec(j+1) = (-1.0D0) * f(1,i)
        g_vec(j+2) = (-1.0D0) * f(2,i)
        g_vec(j+3) = (-1.0D0) * f(3,i)
        j = j + 3
      ENDDO

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Get_Smallest_Box_Dim(smallest_dim)
      use global
      implicit none
C
C       orig: 30 Jan 97  T. Plantenga
C
C   Input/output variable definitions:
C       pgrid(3)         I  (global variable)
C       box(2,3)         I  (global simulation box)
C       smallest_dim     O MIN (x length, y length, z length)
C
      DOUBLE PRECISION  smallest_dim
C
C   This routine accesses global variables to compute the dimensions
C   of the box size used on each processor.
C**********************************************************************

C     ---- LOCAL VARIABLES
      DOUBLE PRECISION  x_dim, y_dim, z_dim

C----------------------------------------------------------------------

C     ---- PROCESSOR BOX SIZE IS DETERMINED BY SIMPLY DIVIDING THE
C     ---- TOTAL EXTENT OF THE MOLECULE IN EACH DIRECTION BY THE NUMBER
C     ---- OF PROCESSORS IN THE DIRECTION.  LAMMPS DETERMINES THE
C     ---- PROCESSOR LAYOUT ONCE AT STARTUP, STORING THE VALUES IN THE
C     ---- ARRAY pgrid.  THE EXTENT OF THE SIMULATION BOX IS STORED IN
C     ---- box(2,3).

      x_dim = box(2,1) - box(1,1)
      x_dim = x_dim / (DBLE (pgrid(1)) )

      y_dim = box(2,2) - box(1,2)
      y_dim = y_dim / (DBLE (pgrid(2)) )

      z_dim = box(2,3) - box(1,3)
      z_dim = z_dim / (DBLE (pgrid(3)) )

      smallest_dim = MIN (x_dim, y_dim)
      smallest_dim = MIN (z_dim, smallest_dim)

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Get_Atom_Positions(localN, x_vec)
      use global
      implicit none
C
C       orig: 25 Oct 96  T. Plantenga
C
C   Input/output variable definitions:
C       nlocal           I  number local atoms (global variable)
C       localN            O 3 * nlocal
C       x                I  atom positions (global variable)
C       x_vec             O atom positions
C
C
      INTEGER           localN
      DOUBLE PRECISION  x_vec(3*maxown)
C
C
C   This routine copies the atom positions vector from the global
C   LAMMPS array to the local x_vec variable.  It also updates n
C   from nlocal, which could change after Repartitn_and_Reneighbor.
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER  i, j, k

C----------------------------------------------------------------------

      j = 0
      DO i = 1, nlocal
         x_vec(j+1) = x(1,i)
         x_vec(j+2) = x(2,i)
         x_vec(j+3) = x(3,i)
         j = j + 3
      ENDDO

      localN = 3 * nlocal

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Put_Atom_Positions(n, x_vec, update_per_trues)
      use global
      use mpi
      implicit none
C
C       orig: 25 Oct 96  T. Plantenga
C
C   Input/output variable definitions:
C       x_vec            I  atom positions
C       update_per_trues I  boolean
C       x                 O atom positions (global variable)
C
C
      INTEGER           n
      DOUBLE PRECISION  x_vec(n)
      INTEGER           update_per_trues
C
C
C   This routine copies the atom positions x_vec into the global
C   LAMMPS array.  Position information is then swapped to other
C   processors that need to know non-resident atom positions.
C
C   External subroutines called:
C        communicate         "communicate.f"
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER           i, j, k
      DOUBLE PRECISION  time1, time2

C----------------------------------------------------------------------

      j = 0
      DO i = 1, nlocal
         x(1,i) = x_vec(j+1)
         x(2,i) = x_vec(j+2)
         x(3,i) = x_vec(j+3)
         j = j + 3
      ENDDO

      time1 = mpi_wtime()
      CALL communicate ()
      time2 = mpi_wtime()

C     ---- UPDATE GLOBAL TIMER.
      time_comm = time_comm + (time2 - time1)

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Compute_Inf_Norm_PAR(n, x_vec, inf_norm)
      use global
      use mpi
      implicit none
C
C       orig: 23 Oct 96  T. Plantenga
C       mod:  21 May 97  TDP - changed from fn to sub for C
C
      integer ierror
      INTEGER           n
      DOUBLE PRECISION  x_vec(n)
      DOUBLE PRECISION  inf_norm
C
C   This routine returns the infinity norm of x_vec, which is a vector
C   strung out across all processors.  Each processor computes the
C   infinity norm of its locally stored portion (in parallel), then
C   merges this with all the other processor norms (requiring
C   interprocessor communication).  At the end, all processors know
C   ||x_vec||_inf.
C
C**********************************************************************

C     ---- LOCAL VARIABLES
C     ---- (VARIABLES COMMUNICATED VIA MPI MUST USE REAL*8).
      INTEGER           i
      DOUBLE PRECISION  dd1, local_max

C----------------------------------------------------------------------

      local_max = 0.0D0
      DO i = 1, n
         dd1 = DABS (x_vec(i))
         IF (dd1 .GT. local_max)  local_max = dd1
      ENDDO

      call mpi_allreduce(local_max,inf_norm,1,mpi_double_precision,
     $     mpi_max,mpi_comm_world,ierror)

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Compute_L2_Norm_PAR(n, x_vec, L2_norm)
      use global
      use mpi
      implicit none
C
C       orig: 25 Oct 96  T. Plantenga
C       mod:  21 May 97  TDP - changed from fn to sub for C
C
      integer ierror
      INTEGER           n
      DOUBLE PRECISION  x_vec(n)
      DOUBLE PRECISION  L2_norm
      real*8 tmp
C
C   This routine returns the L2 norm of a distributed vector x_vec.
C   Each processor computes the L2 norm of its locally stored portion
C   using BLAS (in parallel), then merges this with all the other
C   processor norms (requiring interprocessor communication).  At the
C   end, all processors know ||x_vec||_2.
C
C   External subroutines called:
C        dnrm2               "blas/dnrm2.f"
C**********************************************************************

      DOUBLE PRECISION  dnrm2
      EXTERNAL          dnrm2

C     ---- LOCAL VARIABLES
C     ---- (VARIABLES COMMUNICATED VIA MPI MUST USE REAL*8).
      DOUBLE PRECISION  x_2norm
      REAL*8            global_real

C----------------------------------------------------------------------

      x_2norm = dnrm2 (n, x_vec, 1)

      tmp = x_2norm * x_2norm
      call mpi_allreduce(tmp,global_real,1,mpi_double_precision,
     $     mpi_sum,mpi_comm_world,ierror)

      L2_norm = SQRT (global_real)

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Compute_ddot_PAR(n, x1_vec, x2_vec, dot_product)
      use global
      use mpi
      implicit none
C
C       orig: 31 Oct 96  T. Plantenga
C       mod:  21 May 97  TDP - changed from fn to sub for C
C
      integer ierror
      INTEGER           n
      DOUBLE PRECISION  x1_vec(n), x2_vec(n)
      DOUBLE PRECISION  dot_product
C
C   This routine returns the dot product of two distributed vectors.
C   Each processor computes the dot product of its locally stored
C   portions using BLAS (in parallel), then adds this with all the
C   other processor norms (requiring interprocessor communication).
C   At the end, all processors know x1_vec'x1_vec
C
C   External subroutines called:
C        ddot                "blas/ddot.f"
C**********************************************************************

      DOUBLE PRECISION  ddot
      EXTERNAL          ddot

C     ---- LOCAL VARIABLES
C     ---- (VARIABLES COMMUNICATED VIA MPI MUST USE REAL*8).
      DOUBLE PRECISION  x1_dot_x2

C----------------------------------------------------------------------

      x1_dot_x2 = ddot (n, x1_vec, 1, x2_vec, 1)

      call mpi_allreduce(x1_dot_x2,dot_product,1,
     $     mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Normalize_ddot(dot_product,n)
      use global
      implicit none
C
C     orig: 15 Jan 01  S. Plimpton
C
      DOUBLE PRECISION dot_product
      INTEGER n
C
C   This routine normalizes a previously computed dot-product
C     by the global length of the vector n
C   This is only called after energy estimation performed
C     by the minimizer via Approximate_Bv
C   Only necessary when units = LJ

C**********************************************************************

      if (units == 1) dot_product = dot_product/n

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Find_Biggest_Jump
     +   (boxes_jumped, grid_num_procs, periodic_flag)
      use global
      use mpi
      implicit none
C
C       orig: 31 Jan 97  T. Plantenga
C
C   Input/output variable definitions:
C       x                I  proposed atom positions (global variable)
C       pgrid            I  num procs in each direction (global variable)
C       box              I  (global simulation box)
C       perflagx         I  (global variable)
C       perflagy         I  (global variable)
C       perflagz         I  (global variable)
C       border           I  borders of local processor (global variable)
C       boxes_jumped      O biggest jump in terms of boxes
C       grid_num_procs    O number processors in direction of biggest jump
C       periodic_flag     O zero if all boundaries are periodic
C
      integer ierror
      DOUBLE PRECISION  boxes_jumped
      INTEGER           grid_num_procs, periodic_flag
      real*8 tmp
C
C   This routine has each processor check the current position of their
C   local atoms to see if any have left the local "box" domain.  The
C   number of boxes away (in any one direction) that an atom has moved
C   is computed.  The single largest jump made by any atom in the whole
C   set of processors is returned, along with the number of processors
C   in that direction.
C
C   Parallel communication is used at the end.
C   
C**********************************************************************

C     ---- LOCAL VARIABLES.
      INTEGER           i, j, k
      REAL*8            dd1, dd2, dd3
      REAL*8            box_dim(3)
      REAL*8            border_lo, border_hi
      INTEGER           per_flags(3)
      DOUBLE PRECISION  biggest_jump(3)

C----------------------------------------------------------------------

C     ---- CALCULATE THE SPATIAL DIMENSIONS OF A GENERIC PROCESSOR BOX.
      box_dim(1) = (box(2,1) - box(1,1)) / ( FLOAT (pgrid(1)) )
      box_dim(2) = (box(2,2) - box(1,2)) / ( FLOAT (pgrid(2)) )
      box_dim(3) = (box(2,3) - box(1,3)) / ( FLOAT (pgrid(3)) )

C     ---- LOOK AT PERIODICITY OF THE BOUNDARIES.
      per_flags(1) = perflagx
      per_flags(2) = perflagy
      per_flags(3) = perflagz
      IF ( (per_flags(1) .EQ. 0) .AND.
     +     (per_flags(2) .EQ. 0) .AND.
     +     (per_flags(3) .EQ. 0)       ) THEN
         periodic_flag = 0
      ELSE
         periodic_flag = 1
      ENDIF

C     ---- THIS LOOP IS ADAPTED FROM SUBROUTINE "exchange".
      DO k = 1, 3
         dd1 = 0.0
         IF (pgrid(k) .GT. 1) THEN

            border_lo = border(1,k)
            border_hi = border(2,k)

            DO i = 1, nlocal

C              ---- CALCULATE HOW MANY BOXES AWAY THE ATOM HAS MOVED.
C              ---- COMPENSATE FOR PERIODIC BOUNDARIES IF NECESSARY.

               IF (x(k,i) .LT. border_lo) THEN
                  dd2 = (border_lo - x(k,i)) / box_dim(k)
                  IF (per_flags(k) .EQ. 0) THEN
                     dd3 = (x(k,i) - box(1,k)) / box_dim(k)
                     dd2 = MIN (dd2, dd3)
                  ENDIF
                  dd1 = MAX (dd1, dd2)
               ENDIF

               IF (x(k,i) .GE. border_hi) THEN
                  dd2 = (x(k,i) - border_hi) / box_dim(k)
                  IF (per_flags(k) .EQ. 0) THEN
                     dd3 = (box(2,k) - x(k,i)) / box_dim(k)
                     dd2 = MIN (dd2, dd3)
                  ENDIF
                  dd1 = MAX (dd1, dd2)
               ENDIF

            ENDDO

C           ---- MERGE EACH PROCESSOR'S BIGGEST JUMP TO FIND OUT THE
C           ---- BIGGEST JUMP MADE BY ANY ATOM IN THE PROBLEM.

            tmp = dd1
            call mpi_allreduce(tmp,dd1,1,
     $           mpi_double_precision,mpi_max,mpi_comm_world,ierror)

         ENDIF
         biggest_jump(k) = dd1
      ENDDO

C     ---- FIND THE BIGGEST JUMP OVER THE X, Y, AND Z DIMENSIONS.
      boxes_jumped   = 0.0D0
      grid_num_procs = 1
      DO k = 1, 3
        IF (biggest_jump(k) .GT. boxes_jumped) THEN
           boxes_jumped   = biggest_jump(k)
           grid_num_procs = pgrid(k)
        ENDIF
      ENDDO

      RETURN
      END


c -------------------------------------------------------------------------

      SUBROUTINE Repartitn_and_Reneighbor
      use global
      use mpi
      implicit none
C
C       orig: 20 Jan 97  T. Plantenga
C
C   This routine consists of calls extracted from "integrate" that
C   redefine each processor's domain, move atoms into the correct domains,
C   and then build new neighbor lists for nonbonded interactions.
C
C   Operations are performed in parallel, with substantial interprocessor
C   communication.
C   
C   External subroutines called:
C        pbc                 "pbc.f"
C        setup_box           "setup.f"
C        setup_comm          "setup.f"
C        setup_neigh         "setup.f"
C        exchange            "communicate.f"
C        borders             "communicate.f"
C        neighbor            "neighbor.f"
C        neighbor_bond       "neighbor.f"
C        neighbor_angle      "neighbor.f"
C        neighbor_dihedral   "neighbor.f"
C        neighbor_improper   "neighbor.f"
C**********************************************************************

C     ---- LOCAL VARIABLES.
      DOUBLE PRECISION  time1, time2, time3, time4, time5

C----------------------------------------------------------------------

c same logic as in integrator
c except don't worry about pressstyle or volstyle

      time1 = mpi_wtime()
      call pbc
      if (nonperiodic) call setup_box
      if (nonperiodic) then
        call setup_comm
        if (neighstyle == 1) call setup_neigh
      endif
      call exchange()
      time2 = mpi_wtime()
      call borders ()
      time3 = mpi_wtime()
      call neighbor()
      time4 = mpi_wtime()
      if (nbonds > 0) call neighbor_bond
      if (nangles > 0) call neighbor_angle
      if (ndihedrals > 0) call neighbor_dihedral
      if (nimpropers > 0) call neighbor_improper
      time5 = mpi_wtime()

      time_exch = time_exch + (time2 - time1)
      time_comm = time_comm + (time2 - time1)
      time_neigh1 = time_neigh1 + (time2 - time1)
      time_neigh2 = time_neigh2 + (time2 - time1)

      return
      end


c -------------------------------------------------------------------------

      SUBROUTINE Write_Energy_Vals(procNum)
      use global
      implicit none
C
C       orig: 14 Nov 96  T. Plantenga
C
      INTEGER  procNum
C
C   This routine obtains and writes out the many components of
C   potential energy.  It is assumed energy has been previously
C   calculated by a call to "Compute_f" or "Compute_f_g".
C**********************************************************************

C     ---- LOCAL VARIABLES.
      DOUBLE PRECISION  e_potential

C----------------------------------------------------------------------

      IF (procNum .NE. 0)  RETURN

C     ---- CALCULATION OF e_potential MATCHES CODE IN SUBROUTINES
C     ---- "Compute_f_g" AND "thermo".
      IF (thermostyle .LE. 1) THEN
         e_potential = e_long + e_vdwl + e_coul +
     +                 e_bond + e_angle + e_dihedral + e_improper
      ELSE
         e_potential = e_long + e_vdwl + e_coul +
     +                 e_bond + e_angle + e_dihedral + e_improper +
     +                 e_bondbond + e_bondangle + e_midbondtorsion +
     +                 e_endbondtorsion + e_angletorsion +
     +                 e_angleangletorsion + e_angleangle
      ENDIF

      WRITE (6,600)
      WRITE (6,601) e_potential
      WRITE (6,602) e_bond, e_coul
      WRITE (6,603) e_angle, e_vdwl
      WRITE (6,604) e_dihedral, e_long
      WRITE (6,605) e_improper
      WRITE (6,606)
      write (6,*)

      WRITE (1,600)
      WRITE (1,601) e_potential
      WRITE (1,602) e_bond, e_coul
      WRITE (1,603) e_angle, e_vdwl
      WRITE (1,604) e_dihedral, e_long
      WRITE (1,605) e_improper
      WRITE (1,606)
      write (1,*)

      RETURN

C----------------------------------------------------------------------
C   THE REMAINING LINES OF CODE ARE I/O STATEMENTS REFERENCED EARLIER.
C----------------------------------------------------------------------

 600  FORMAT (/'======================================================',
     +        '==============')
 601  FORMAT ('  Total PE  = ', F14.4)
 602  FORMAT ('  E (bond)  = ', F14.4, 10X, 'E (coulomb) = ', F14.4)
 603  FORMAT ('  E (angle) = ', F14.4, 10X, 'E (VDW)     = ', F14.4)
 604  FORMAT ('  E (dih)   = ', F14.4, 10X, 'E (long)    = ', F14.4)
 605  FORMAT ('  E (imprp) = ', F14.4)
 606  FORMAT ('======================================================',
     +        '==============')

      END
