LAMMPS = Large-scale Atomic/Molecular Massively Parallel Simulator

This is the documentation for the LAMMPS 2001 version, written in F90,
which has been superceded by more current versions.  See the LAMMPS
WWW Site at http://www.cs.sandia.gov/~sjplimp/lammps.html for more
information.

LAMMPS is a classical molecular dynamics code designed for simulating
molecular and atomic systems on parallel computers using
spatial-decomposition techniques. It runs on any parallel platform
that supports F90 and the MPI message-passing library or on
single-processor workstations.

LAMMPS 2001 is copyrighted code that is distributed freely as
open-source software under the GNU Public License (GPL).  See the
LICENSE file or www.gnu.org for more details.  Basically the GPL
allows you as a user to use, modify, or distribute LAMMPS however you
wish, so long as any software you distribute remains under the GPL.

Features of LAMMPS 2001 include:
    
    short-range pairwise Lennard-Jones and Coulombic interactions 
    long-range Coulombic interactions via Ewald or PPPM (particle-mesh 
        Ewald) 
    short-range harmonic bond potentials (bond, angle, torsion, improper) 
    short-range class II (cross-term) molecular potentials 
    NVE, NVT, NPT dynamics 
    constraints on atoms or groups of atoms 
    rRESPA long-timescale integrator 
    energy minimizer (Hessian-free truncated Newton method) 

For users of LAMMPS 99, this version is written in F90 to take
advantage of dynamic memory allocation.  This means the user does not
have to fiddle with parameter settings and re-compile the code so
often for different problems.  This enhancment means there are new
rules for the ordering of commands in a LAMMPS input script, as well
as a few new commands to guide the memory allocator. Users should read
the beginning sections of the input_commands file for an explanation.

More details about the code can be found here, in the HTML- or
text-based documentation.  The LAMMPS Web page is at
www.cs.sandia.gov/~sjplimp/lammps.html, which includes benchmark
timings and a list of papers written using LAMMPS results.  They
illustrate the kinds of scientific problems that can be modeled with
LAMMPS.  The following papers describe the parallel algorithms used in
the code.  Please cite these if you incorporate LAMMPS results in your
work.  And if you send me citations for your papers, I'll be pleased
to add them to the LAMMPS WWW page.

S. J. Plimpton, R. Pollock, M. Stevens, &quot;Particle-Mesh Ewald and
rRESPA for Parallel Molecular Dynamics Simulations&quot;, in Proc of
the Eighth SIAM Conference on Parallel Processing for Scientific
Computing, Minneapolis, MN, March 1997.

S. J. Plimpton, "Fast Parallel Algorithms for Short-Range Molecular
Dynamics", J Comp Phys, 117, 1-19 (1995).

LAMMPS was originally developed as part of a 5-way CRADA collaboration
between 3 industrial partners (Cray Research, Bristol-Myers Squibb,
and Dupont) and 2 DoE laboratories (Sandia National Laboratories and
Lawrence Livermore National Laboratories).

The primary author of LAMMPS is Steve Plimpton, but others have written 
or worked on significant portions of the code:
    
    Roy Pollock (LLNL): Ewald, PPPM solvers 
    Mark Stevens (Sandia): rRESPA, NPT integrators 
    Eric Simon (Cray Research): class II force fields 
    Todd Plantenga (Sandia): energy minimizer 
    Steve Lustig (Dupont): msi2lmp tool 
    Mike Peachey (Cray Research): msi2lmp tool 

Other CRADA partners involved in the design and testing of LAMMPS are 
    
    John Carpenter (Cray Research) 
    Terry Stouch (Bristol-Myers Squibb) 
    Jim Belak (LLNL) 

If you have questions about LAMMPS, please contact me:

    Steve Plimpton
    sjplimp@sandia.gov 
    www.cs.sandia.gov/~sjplimp 
    Sandia National Labs 
    Albuquerque, NM 87185

More Information about LAMMPS

    Basics 
        how to make, run, and test LAMMPS with the example problems 
    Input Commands 
        a complete listing of input commands used by LAMMPS 
    Data Format 
        the data file format used by LAMMPS 
    Force Fields 
        the equations LAMMPS uses to compute force-fields 
    Units 
        the input/output and internal units for LAMMPS variables 
    History 
        a brief timeline of features added to LAMMPS 
    Deficiencies 
        features LAMMPS does not (yet) have
