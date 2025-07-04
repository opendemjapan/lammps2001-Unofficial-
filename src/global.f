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
c all global variables

      module global

c -------------------------------------------------------------------------
c *** parameters unlikely to need changing

c maxfix = max # of fixes
c maxdiag = max # of user defined diagnostics

c minown = min size of atom arrays on one proc
c minghost = min size of ghost atom arrays on one proc
c minneigh = min size of neighbor list on one proc
c minbuf = min size of communication bufs on one proc

      integer, parameter :: maxfix = 100
      integer, parameter :: maxdiag = 10

      integer, parameter :: minown = 1000
      integer, parameter :: minghost = 1000
      integer, parameter :: minneigh = 10000
      integer, parameter :: minbuf = 10000

c -------------------------------------------------------------------------
c *** constants

c boltz = Boltzmann factor
c dtfactor = multiplier on dt into time units
c efactor = conversion factor for Coulombic to energy units
c pfactor = conversion factor for pressure units
c two_1_3 = 2^(1/3)

      real*8 boltz
      real*8 dtfactor
      real*8 efactor
      real*8 pfactor
      real*8 two_1_3

c -------------------------------------------------------------------------
c *** global settings and flags

c iversion = internal version # of code
c units = 0 for conventional, 1 for LJ units
c idimension = 2 or 3 for 2d/3d
c nsteps = # of timesteps to run
c itime = current timestep from 1 to nsteps
c ntimestep = current global timestep #
c ntime_last = what the global timestep will be on last step of run
c firstflag = 0 if haven't read any input command yet, 1 if have
c readflag = 0 if haven't read-in atoms yet, 1 if have
c newton = user input combination (0-3) of nonbond and bonded Newton flag
c newton_nonbond = 0 if Newton's law off for nonbond forces, 1 if on
c newton_bond = 0 if Newton's law off for bonded forces, 1 if on

c dt = timestep
c dthalf = 1/2 the timestep

      integer iversion
      integer units
      integer idimension
      integer nsteps
      integer itime
      integer ntimestep
      integer ntime_last
      integer firstflag,readflag
      integer newton,newton_nonbond,newton_bond

      real*8 dt,dthalf

c -------------------------------------------------------------------------
c *** domain

c xprd,yprd,zprd = size of global simulation box
c xprd_half,yprd_half,zprd_half = 1/2 the box lengths
c box(2,3) = lower/upper boundary of 3 dims of global box
c border(2,3) = lower/upper boundary of 3 dims of my sub-domain
c perflagx,perflagy,perflagz = 0 for periodic, 1 for non-periodic in 3 dims
c slabflag = 0 for 3-D periodicity of long-range forces, 1 for slab
c nonperiodic = TRUE if any dim is non-periodic, FALSE otherwise
c slab_volfactor = ratio of total volume to occupied volume for slab geometry
c zprd_slab = size of global box for slab geometry long-range interactions

      real*8 xprd,yprd,zprd
      real*8 xprd_half,yprd_half,zprd_half

      real*8 box(2,3)
      real*8 border(2,3)

      integer perflagx,perflagy,perflagz
      integer slabflag
      logical nonperiodic

      real*8 slab_volfactor,zprd_slab

c -------------------------------------------------------------------------
c *** atoms

c natoms = total # of atoms in simulations
c ntypes = # of atom types
c maxown = most # of atoms I can own
c maxghost = most # of ghost atoms I can acquire
c extra_own = multiplier on allocation of maxown
c extra_ghost = multiplier on allocation of maxghost
c maxatom = maxown + maxghost
c maxatomp1 = maxatom + 1, used for storing a flag for special interactions

c nlocal = # of atoms I own
c nghost = # of ghost atoms I have
c max_nlocal = most atoms I ever own during simulation
c max_ghost = most ghost atoms I ever acquire during simulation

c mass(ntypes) = mass of each type of atom
c x(3,n) = coords of my atoms, owned and ghost
c v(3,n) = velocities of owned atoms
c f(3,n) = forces on my atoms, owned and ghost
c q(n) = charge of my atoms, owned and ghost
c xhold(3,n) = coords of owned atoms at last neighbor list build
c fhold(3,n) = force storage in integrate.f for newton=3 virial computation
c tag(n) = global ID tags of my atoms, owned and ghost
c type(n) = type of my atoms, owned and ghost
c molecule(n) = molecule ID # of owned atoms
c true(n) = true flag of owned atoms, which 3-d image of box they are in

      integer natoms,ntypes
      integer maxown,maxghost
      real*8 extra_own,extra_ghost
      integer maxatom,maxatomp1

      integer nlocal,nghost
      integer max_nlocal,max_nghost

      real*8, allocatable :: mass(:)

      real*8, allocatable :: x(:,:),v(:,:),f(:,:)
      real*8, allocatable :: q(:)
      real*8, allocatable :: xhold(:,:),fhold(:,:)
      integer, allocatable :: tag(:),type(:),molecule(:),true(:)

c -------------------------------------------------------------------------
c *** bond connectivity for each atom

c nmolecular = 0 if an atomic system, 1 if molecular
c maxbondper,maxangleper,maxdihedper,maximproper = max # of bonds, etc
c       that any atom in simulation must store, depends on newton_bond
c maxspec = max # of 1-2, 1-3, 1-4 neighbors any atom must store

c numbond(n) = # of bonds stored by each of my atoms
c numangle(n) = # of angles stored by each of my atoms
c numdihed(n) = # of dihedrals stored by each of my atoms
c numimpro(n) = # of impropers stored by each of my atoms

c bondtype,angletype,dihedtype,improtype = 
c  type of each bond,angle,dihedral,improper for each of my atoms
c bondatom,angleatom,dihedatom,improatom 1234 = 
c  global tags of 2/3/4 atoms in each bond,angle,dihedral,improper
c  for each of my atoms

c num 123 bond = how many 1-2, 1-3, 1-4 neighbors this atom has in specbond
c                is cummulative: num2 includes 1-2 count, num3 has 1-2 & 1-3
c specbond = list of global ID tags of 
c            1-2, 1-3, 1-4 neighbors of each of my owned atoms

      integer nmolecular
      integer maxbondper,maxangleper,maxdihedper,maximproper
      integer maxspec

      integer, allocatable :: numbond(:),numangle(:)
      integer, allocatable :: numdihed(:),numimpro(:)

      integer, allocatable :: bondtype(:,:)
      integer, allocatable :: bondatom1(:,:),bondatom2(:,:)
      integer, allocatable :: angletype(:,:)
      integer, allocatable :: angleatom1(:,:),angleatom2(:,:)
      integer, allocatable :: angleatom3(:,:)
      integer, allocatable :: dihedtype(:,:)
      integer, allocatable :: dihedatom1(:,:),dihedatom2(:,:)
      integer, allocatable :: dihedatom3(:,:),dihedatom4(:,:)
      integer, allocatable :: improtype(:,:)
      integer, allocatable :: improatom1(:,:),improatom2(:,:)
      integer, allocatable :: improatom3(:,:),improatom4(:,:)

      integer, allocatable :: num1bond(:),num2bond(:),num3bond(:)
      integer, allocatable :: specbond(:,:)

c -------------------------------------------------------------------------
c *** global ptr

c localptr(i) = local index of atom with global tag of i
c               0 if this proc doesn't own the atom or have it as a ghost

      integer, allocatable :: localptr(:)

c -------------------------------------------------------------------------
c *** nonbond LJ and Coulombic interactions

c nonstyle = style of nonbond VanderWaal interactions
c            0 = none, 1 = cutoff LJ, 2 = smoothed LJ, 3 = shifted LJ,
c            4 = soft potential, 5 = class 2
c offsetflag = whether to add in shifted LJ energy at cutoff distance
c mixflag = 0 if no user input of mixing style, 1 if has been specified
c mixstyle = nonbond mixing style for epsilon,sigma
c            1 = geometric, 2 = arithmetic, 3 = sixthpower (class 2)
c ncharge = 0 if no charges defined in system, 1 if there are charges
c coulstyle = style of Coulombic interactions
c           0 = none, 1 = cutoff, 2 = smoothed, 3 = Ewald, 4 = PPPM
c amberflag = 0/1 if special bonds are set to AMBER force field settings
c coulpre = prefactor on Coulombic energy (efactor/dielectric)
c dielectric = dielectric constant settable by user
c special(3) = nonbond weighting factors on 1-2, 1-3, 1-4 neighbors

c cutforce,cutforce_sq = longest force cutoff of any nonbond interaction
c cutlj,cutlj_sq = cutoff for LJ
c cutljinterior,cutljint_sq = interior cutoff for smoothed LJ
c cutcoul = cutoff for Coulombic
c cutcoulsq = square of Coulombic cutoff
c cutcoulint = interior cutoff for smoothed Coulombic
c cutcoulintsq = square of interior cutoff for smoothed Coulombic
c cutmax = cutoff set by "maximum cutoff" command for memory allocator
c cutmemory = cutoff used by memory allocator, includes neighbor skin
c ch_denom_lj = charmm switching function denominator for LJ interactions
c ch_denom_coul = charmm switching fnct. denom. for Coulombic interactions
c kappa = damping factor for Debye/Huckel interactions

c noncoeff 1234 = nonbond coefficients as input for each atom type pair
c noncoeff14 12 = 1-4 interaction nonbond coefficients as input
c nontypeflag = whether coeffs have been specified for each type pair
c cutforcesq = force cutoff for each type pair
c cutneighsq = neighbor cutoff for each type pair
c cutljsq = LJ cutoff for each type pair
c cutljinner, cutljinnersq = interior LJ cutoff for each type pair
c lj 12345 = nonbond coeffs as used in force routines for each type pair
c lj14 1234 = 1-4 interaction nonbond coeffs as used in charmm force routine
c ljsw 01234 = smoothed LJ coeffs as used in force routines for each type pair
c offset = offset LJ energy for each type pair

      integer nonstyle
      integer offsetflag
      integer mixflag,mixstyle
      integer ncharge
      integer coulstyle
      integer amberflag
      real*8 coulpre
      real*8 dielectric
      real*8 special(3)

      real*8 cutforce,cutforce_sq
      real*8 cutlj,cutljinterior,cutlj_sq,cutljint_sq
      real*8 cutcoul,cutcoulsq,cutcoulint,cutcoulintsq
      real*8 cutmax,cutmemory
      real*8 ch_denom_lj,ch_denom_coul
      real*8 kappa

      real*8, allocatable :: noncoeff1(:,:),noncoeff2(:,:)
      real*8, allocatable :: noncoeff3(:,:),noncoeff4(:,:)
      real*8, allocatable :: noncoeff14_1(:,:),noncoeff14_2(:,:)
      integer, allocatable :: nontypeflag(:,:)

      real*8, allocatable :: cutforcesq(:,:),cutneighsq(:,:)
      real*8, allocatable :: cutljsq(:,:)
      real*8, allocatable :: cutljinner(:,:),cutljinnersq(:,:)
      real*8, allocatable :: lj1(:,:),lj2(:,:),lj3(:,:)
      real*8, allocatable :: lj4(:,:),lj5(:,:)
      real*8, allocatable :: lj14_1(:,:),lj14_2(:,:),lj14_3(:,:)
      real*8, allocatable :: lj14_4(:,:)
      real*8, allocatable :: ljsw0(:,:),ljsw1(:,:),ljsw2(:,:)
      real*8, allocatable :: ljsw3(:,:),ljsw4(:,:)
      real*8, allocatable :: offset(:,:)

c -------------------------------------------------------------------------
c *** bonded interactions

c nbonds = # of bonds in entire simulation
c nangles = # of angles in entire sim
c ndihedrals = # of dihedrals in entire sim
c nimpropers = # of impropers in entier sim
c nbondtypes,nangletypes,ndihedtypes,nimprotypes = 
c              # of types of each interaction
c bondstyle = style of bond interactions
c             0 = none, 1 = harmonic, 2 = FENE, 3 = shifted FENE,
c             4 = nonlinear spring, 5 = class 2
c anglestyle = style of angle interactions
c              0 = none, 1 = harmonic, 2 = class 2
c dihedstyle = style of dihedral interactions
c              0 = none, 1 = harmonic, 2 = class 2
c improstyle = style of improper interactions
c              0 = none, 1 = harmonic, 2 = cvff, 3 = class 2

c bondcoeff = coeffs for all bond types
c bondtypeflag = whether coeffs have been set for each type
c anglecoeff = coeffs for all angle types
c angletypeflag = whether coeffs have been set
c dihedcoeff = coeffs for all dihedral types
c dihedtypeflag = whether coeffs have been set
c improcoeff = coeffs for all improper types
c improtypeflag = whether coeffs have been set

      integer nbonds,nangles,ndihedrals,nimpropers
      integer nbondtypes,nangletypes,ndihedtypes,nimprotypes
      integer bondstyle,anglestyle,dihedstyle,improstyle

      real*8, allocatable :: bondcoeff(:,:)
      integer, allocatable :: bondtypeflag(:)
      real*8, allocatable :: anglecoeff(:,:)
      integer, allocatable :: angletypeflag(:)
      real*8, allocatable :: dihedcoeff(:,:)
      integer, allocatable :: dihedtypeflag(:)
      real*8, allocatable :: improcoeff(:,:)
      integer, allocatable :: improtypeflag(:)

c -------------------------------------------------------------------------
c *** class 2 force field

c coeffs for all class 2 interactions

      real*8, allocatable :: bondbondcoeff(:,:)
      real*8, allocatable :: bondanglecoeff(:,:)
      real*8, allocatable :: bondbond13coeff(:,:)
      real*8, allocatable :: angleanglecoeff(:,:)
      real*8, allocatable :: angleangletorsioncoeff(:,:)
      real*8, allocatable :: angletorsioncoeff(:,:)
      real*8, allocatable :: midbondtorsioncoeff(:,:)
      real*8, allocatable :: endbondtorsioncoeff(:,:)

c -------------------------------------------------------------------------
c *** nonbond neighbor lists

c maxneigh = most # of neighbors that can be stored
c max_neigh = = most # ever stored during simulation
c extra_neigh = multiplier on allocation of maxneigh
c numneigh = # of times neighor lists are rebuilt
c ndanger = # of times a "dangerous" rebuild is done, where an atom may
c           have moved within a force-cutoff distance earlier
c maxbin = # of neighbor bins I own
c nbinx,nbiny,nbinz = # of bins in each dim of global simulation box
c mbinx,mbiny,mbinz = # of bins in each dim of my sub-domain, including ghosts
c mbinxlo,mbinylo,mbinzlo = lowest global bin one of my atoms could be in
c nstencil = # of bins in stencil for checking neighbor interactions
c stencil = list of offsets into set of bins that comprise the stencil
c binsizex,binsizey,binsizez = size of a neighbor bin in 3-d
c bininvx,bininvy,bininvz = inverse sizes of neighbor bins
c cutneigh = neighbor cutoff
c skin = distance beyond force cutoff that neighbor cutoff includes
c triggersq = trigger distance (1/2 of skin) for rebuilding neighbor lists
c neighago = how many steps ago neighbor list was rebuilt
c neighdelay = delay for at least this may steps before rebuilding lists
c neighfreq = check rebuild criterion or rebuild every this many steps
c neighstyle = 0 for brute-force N^2 search, 1 for binning
c neightrigger = 0 if don't rebuild based on distance moved, 1 if do
c nlist = 1-d list of all neighors of all my atoms
c nnfirst,nnlast = ptrs into nlist where neighbor of each atom start/stop
c bin = ptr from each atom to the next atom in its neighbor bin
c binpnt = ptr to 1st atom in each neighbor bin

      integer maxneigh,max_neigh
      real*8 extra_neigh
      integer numneigh
      integer ndanger

      integer maxbin
      integer nbinx,nbiny,nbinz
      integer mbinx,mbiny,mbinz
      integer mbinxlo,mbinylo,mbinzlo

      integer nstencil
      integer stencil(1000)

      real*8 binsizex,binsizey,binsizez
      real*8 bininvx,bininvy,bininvz

      real*8 cutneigh
      real*8 skin
      real*8 triggersq

      integer neighago,neighdelay,neighfreq,neighstyle,neightrigger

      integer, allocatable :: nlist(:)
      integer, allocatable :: nnfirst(:),nnlast(:)

      integer, allocatable :: bin(:)
      integer, allocatable :: binpnt(:)

c -------------------------------------------------------------------------
c *** bonded lists

c nbondlocal,nanglelocal,ndihedlocal,nimprolocal =
c   # of bond,angle,dihedral,improper in current bonded lists
c maxbondlocal,maxanglelocal,maxdihedlocal,maximprolocal
c   max # of bond,angle,dihedral,improper my allocated lists can store
c max_bond,max_angle,max_dihed,max_impro =
c   most # of bond,angle,dihedral,improper my lists every store
c bondlist,anglelist,dihedlist,improlist =
c   lists for bond,angle,dihedral,improper interactions to compute,
c   each entry in list stores type and global tag IDs of atoms involved

      integer nbondlocal,nanglelocal,ndihedlocal,nimprolocal
      integer maxbondlocal,maxanglelocal,maxdihedlocal,maximprolocal
      integer max_bond,max_angle,max_dihed,max_impro

      integer, allocatable :: bondlist(:,:),anglelist(:,:)
      integer, allocatable :: dihedlist(:,:),improlist(:,:)

c -------------------------------------------------------------------------
c *** communication

c node = proc ID of me
c nprocs = total # of procs
c nswap = # of atom swaps each proc does with surrounding procs
c maxswap = current size of swap arrays
c max_exch = most # of atoms I've ever exchanged
c max_bord = most # of atoms sent in one border swap
c max_slist = most # of atoms I send in all my border swaps
c             (max_ghost is most I ever receive)
c maxbuf = size of allocated communication buffers
c max_buf = largest portion of buf ever used during simulation
c extra_buf = multiplier on allocation of maxbuf

c pgrid(3) = # of procs assigned to each dim on simulation box
c me(3) = which sub-box I own in each of 3 dims (0 to n-1)
c mpart(2,3) = my 6 neighboring procs
c need(3) = how many boxes away I need ghost info from in each dim

c slablo,slabhi = boundary inside which I need ghost info for in each swap
c spart = proc to send to in each swap
c rpart = proc to recv from in each swap
c nsfirst,nslast = ptrs into slist of send atoms for each swap
c nrfirst,nrlast = ptrs into ghost list where to put recv info in each swap
c commflag(3,n) = flags for PBC treatment of each dim in each swap
c commflagall(n) = flag for whether PBC treatment is needed in each swap
c slist = list of local atoms to send in all my swaps
c ibuf1,ibuf2,buf1,buf2 = buffers to use in communicating ghost info

      integer node,nprocs
      integer nswap,maxswap
      integer max_exch,max_bord,max_slist
      integer maxbuf,max_buf
      real*8 extra_buf

      integer pgrid(3),me(3),mpart(2,3),need(3)

      real*8, allocatable :: slablo(:),slabhi(:)
      integer, allocatable :: spart(:),rpart(:)
      integer, allocatable :: nsfirst(:),nslast(:)
      integer, allocatable :: nrfirst(:),nrlast(:)
      integer, allocatable :: commflag(:,:),commflagall(:)

      integer, allocatable :: slist(:)

      integer, allocatable :: ibuf1(:),ibuf2(:)
      real*8, allocatable :: buf1(:),buf2(:)

c -------------------------------------------------------------------------
c *** thermodynamics

c t_current = current temperature
c e_total = potential energy
c p_total = scalar pressure
c p_current(3) = pressure trace
c virial(6) = diagonal and off-diagonal virial components
c virialhold(6) = temporary copy of virial for ghost-atom virial computation
c vir_long(6) = virial for long-range Ewald/PPPM forces

      real*8 t_current
      real*8 e_total
      real*8 p_total
      real*8 p_current(3)
      real*8 virial(6),virialhold(6),vir_long(6)

c nonbond and bonded potential energies
c e_14 terms are for nonbond portion computed in CHARMM dihedrals

      real*8 e_vdwl
      real*8 e_coul
      real*8 e_bond
      real*8 e_angle
      real*8 e_dihedral
      real*8 e_improper
      real*8 e_14_coul,e_14_vdwl

c long-range Coulombic energy

      real*8 e_long

c class 2 potential energies

      real*8 e_bondbond
      real*8 e_bondbond13
      real*8 e_bondangle
      real*8 e_endbondtorsion
      real*8 e_midbondtorsion
      real*8 e_angletorsion
      real*8 e_angleangletorsion
      real*8 e_angleangle

c -------------------------------------------------------------------------
c *** ensemble variables: temp, pressure, volume control

c ensemble = which ensemble is being simulated: 1=NVE, 2=NVT, 3=NPH, 4=NPT

c tempstyle = style of temperature control
c             0 = none, 1 = rescale, 2 = replace, 3 = Langevin, 4 = Nose/Hoover
c t_every = check for temp rescaling every this many steps
c t_start,t_stop = desired temperature at start/end of run
c t_window = rescale temperature if it is outside this window
c t_fraction = amount (0% to 100%) of rescaling to perform
c t_freq = drag/mass parameter in Langevin and Nose/Hoover temp control
c t_target = target temperature on this timestep
c eta,eta_dot = evolving Nose variable for temp control

c pressstyle = style of pressure control: 0 = none, 1 = Nose/Hoover
c presscouple = style of coupling:
c               0 = xyz isotropic, 1 = xy, 2 = yz, 3 = xz, 4 = anisotropic
c p_freq(3) = piston mass parameter in Nose/Hoover pressure control
c p_start(3),p_stop(3) = desired pressure at start/end of run
c p_target(3) = target pressure on this timestep
c omega(3),omega_dot(3) = evolving Nose variables for pressure control
c masssum = total mass in system

c volstyle = style of volume contfol: 0 = none, 1 = linear expand/contract
c voldimx,voldimy,voldimz = disabled (0) or enabled (1) for each dimension
c volstart/stop xyz lo/hi = global box boundary at begin/end of run in each dim

      integer ensemble

      integer tempstyle
      integer t_every
      real*8 t_start,t_stop,t_window,t_fraction
      real*8 t_freq,t_target
      real*8 eta,eta_dot

      integer pressstyle,presscouple
      real*8 p_freq(3)
      real*8 p_start(3),p_stop(3),p_target(3)
      real*8 omega(3),omega_dot(3)
      real*8 masssum

      integer volstyle,voldimx,voldimy,voldimz
      real*8 volstart_xlo,volstart_xhi,volstop_xlo,volstop_xhi
      real*8 volstart_ylo,volstart_yhi,volstop_ylo,volstop_yhi
      real*8 volstart_zlo,volstart_zhi,volstop_zlo,volstop_zhi

c -------------------------------------------------------------------------
c *** temperature creation

c createstyle = style of velocity creation
c               1 = uniform, 2 = gaussian, 3 = explicit velocity
c creategroup = kind of group of atoms to create vels for
c               1 = types, 2 = region, 3 = remainder of unset atoms
c createlo,createhi = range of types/molecules to create vels for
c iseed = random # seed, used for velocity creation and temp control
c rotationflag = 0 if don't zero angular momentum, 1 if do
c t_create = desired temperature to create vels for
c createregion(6) = geometric bounds of atoms to create vels for
c createvec(3) = explicit velocity vector to apply to atoms in group

c velflag(:) = used to mark atoms with previously initialized vels

      integer createstyle,creategroup,createlo,createhi
      integer iseed,rotationflag
      real*8 t_create
      real*8 createregion(6)
      real*8 createvec(3)

      integer, allocatable :: velflag(:)

c -------------------------------------------------------------------------
c *** fixes

c nfixes = # of fixes specified by user
c nfixes_respa = # of setforce and aveforce fixes that must be applied
c                at each level of rRESPA

c fixstyle = style of each fix, 0 = none, 1 = setforce, 2 = addforce, 
c            3 = aveforce, 4 = rescale, 5 = Langevin, 6 = Nose/Hoover,
c            7 = springforce, 8 = dragforce, 9 = shake
c fixflag = flags associated with each fix
c fixptr = ptr into fixstore where values associated with a fix are stored
c fixcount = # of atoms assigned to each fix
c fixactive = flag for whether each fix is active on this timestep
c fixcoeff = coeffs/values/constants assocated with each fix
c fixmass = total mass of all atoms assigned to each fix

c fixwhich = which fix (1-N) to assign atoms to
c fixgroup = kind of group that the fix is to be applied to
c            1 = single atom, 2 = molecule, 3 = type, 4 = region, 5 = remainder
c fixatom = atom tag ID or molecule ID to apply fix to
c fixtype,fixbond,fixangle = atom/bond/angle type to apply fix to
c fixregion(6) = geometric region to apply fix to

c fixnum = how many fix quantities need to be summed up across procs
c fixstore = vector of fix values computed on this timestep
c fixstoretmp,fixmasstmp,fixcounttmp = scratch vectors for summing fix info

c fix = fix #'s associated with each of my atoms

      integer nfixes,nfixes_respa

      integer fixstyle(maxfix),fixflag(3,maxfix),fixptr(maxfix)
      integer fixcount(maxfix),fixactive(maxfix)
      real*8 fixcoeff(7,maxfix),fixmass(maxfix)

      integer fixwhich,fixgroup,fixatom,fixtype,fixbond,fixangle
      real*8 fixregion(6)

      integer fixnum
      real*8 fixstore(3*maxfix)

      real*8 fixstoretmp(3*maxfix),fixmasstmp(maxfix)
      integer fixcounttmp(maxfix)

      integer, allocatable :: fix(:)

c -------------------------------------------------------------------------
c *** Ewald (some of these variables are also used by PPPM)

c kcount = actual # of Ewald vectors
c kmax = max dimensionality of Ewald vectors
c gewald = G vector for Ewald/PPPM as function of precision/cutoff
c qsum,qsqsum = total charge, square of total charge
c long_prec = user-specified precision for long-range approximation

c kxvecs,kyvecs,kzvecs = pre-computed indices of Ewald vectors
c ug = pre-computed exponential factor for each Ewald vector
c eg(3,n) = pre-computed Ewald vector components of ug
c vg(6,n) = pre-computed virial components for each Ewald/PPPM vector
c sfacrl,sfacim = real/imag structure factors for each Ewald vector
c sfacrl_all,sfacim_all = scratch vectors for summing sfacs across procs
c cs,sn = atom contributions to each Ewald level

c ek(3,n) = components of electric field at each atom

      integer kcount,kmax
      real*8 gewald
      real*8 qsum,qsqsum
      real*8 long_prec

      integer, allocatable :: kxvecs(:),kyvecs(:),kzvecs(:)
      real*8, allocatable :: ug(:),eg(:,:),vg(:,:)
      real*8, allocatable :: sfacrl(:),sfacim(:)
      real*8, allocatable :: sfacrl_all(:),sfacim_all(:)
      real*8, allocatable :: cs(:,:,:),sn(:,:,:)

      real*8, allocatable :: ek(:,:)

c -------------------------------------------------------------------------
c *** PPPM

c orderflag = order of PPPM, how far into grid the charge overlaps
c meshflag = 1 if user sets PPPM mesh, 0 otherwise
c nfft = # of FFT points I own in FFT decomp
c nlower,nupper = # of grid pts a charge extends to the left,right
c nx_pppm,ny_pppm,nz_pppm = global PPPM grid
c nx_pppm_input,ny_pppm_input,nz_pppm_input = global PPPM grid chosen by user
c nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in = 
c   portion of global grid I own in brick decomp
c nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out = 
c   portion of global grid I own including ghost cells
c nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost =
c   # of planes of ghost cells I receive from neighbor in each direction
c nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft =
c   portion of global grid I own in FFT decomp
c plan1_fft,plan2_fft,plan_remap = FFT and remap plans

c maxgrid = size of local grid array including ghost cells in brick decomp
c maxfft = size of local FFT array, bigger of brick or FFT decomp
c maxpbuf = size of buffers for exchanging ghost cells

c density_brick = stores particle charge as mapped to grid, brick decomp
c vdx_brick,vdy_brick,vdz_brick = stores potential gradient as mapped to grid
c density_fft = particle charge on grid in FFT decomp
c greensfn = Green's function for each point on grid
c workvec1,workvec2 = complex work vectors for FFTs
c partgrid = which grid cell (nx,ny,nz) a particle is centered at
c fkvecs_x,fkvecs_y,fkvecs_z = pre-computed FFT coefficients
c pbuf1,pbuf2 = buffers for exchanging ghost cells

      integer orderflag
      integer meshflag
      integer nfft
      integer nlower,nupper
      integer nx_pppm,ny_pppm,nz_pppm
      integer nx_pppm_input,ny_pppm_input,nz_pppm_input
      integer nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in
      integer nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out
      integer nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost
      integer nzlo_ghost,nzhi_ghost
      integer nxlo_fft,nxhi_fft,nylo_fft,nyhi_fft,nzlo_fft,nzhi_fft
      real*8 plan1_fft,plan2_fft,plan_remap

      integer maxgrid
      integer maxfft
      integer maxpbuf

      real*8, allocatable :: density_brick(:)
      real*8, allocatable :: vdx_brick(:),vdy_brick(:),vdz_brick(:)
      real*8, allocatable :: density_fft(:)
      real*8, allocatable :: greensfn(:)
      complex*16, allocatable :: workvec1(:),workvec2(:)
      integer, allocatable :: partgrid(:,:)
      real*8, allocatable :: fkvecs_x(:),fkvecs_y(:),fkvecs_z(:)
      real*8, allocatable :: pbuf1(:),pbuf2(:)

c -------------------------------------------------------------------------
c *** rRESPA

c nrespa = 1 if rRESPA is on, 0 if off
c nstretch,nintra,nshort = dilation factors between 4 sets of timesteps
c vir_stretch(6),vir_intra(6),vir_short(6) = virials for 3 forces
c dthalf_intra,dthalf_short,dthalf_long = timesteps for 3 scales besides dt
c f_stretch,f_intra,f_short,f_long = force arrays for 4 forces

      integer nrespa
      integer nstretch,nintra,nshort
      real*8 vir_stretch(6),vir_intra(6),vir_short(6)
      real*8 dthalf_intra,dthalf_short,dthalf_long

      real*8, allocatable :: f_stretch(:,:),f_intra(:,:)
      real*8, allocatable :: f_short(:,:),f_long(:,:)

c -------------------------------------------------------------------------
c *** SHAKE

c nshake = 1 if SHAKE is on, 0 if off
c shakeiter = max # of iterations SHAKE will attempt
c nshakestats = print bond statistics every this many timesteps (0 = never)
c shakewhichbondcoeff = which bond coeff is bond length (for bondstyle)
c shakeableangle = 0/1 if angle constraint is not set or set
c shakeableanglebond = which bond type the angle constraint applies to
c nshake_next = next timestep to call bond statistics routine
c shaketol = tolerance for SHAKE bonds
c shakeanglebond = psuedo-bond distance across constrained angle
c shakeablebond = 0/1 for each bondtype, whether (0) not SHAKE or (1) SHAKE
c nshakebonds = total # of bonds & pseudo-bonds constrained by SHAKE

c shakegroup = 0,2,3,4 = size of SHAKE group this atom is part of
c shakepartner = global ID tag of all atoms (including self) in SHAKE group,
c                central atom is listed 1st
c shakebondtype = which type of bond each SHAKE bond is

c nshakelocal = # of SHAKE groups this proc must currently compute (neigh list)
c shakesize = size of each SHAKE group (2,3,4) (in neigh list)
c shakeatom = local IDs for each atom in SHAKE group (in neigh list)
c shakebondlen = bond length of each bond in SHAKE group (in neigh list)
c xshake = updated unconstrained atom coords for owned and ghost atoms

      integer nshake,shakeiter,nshakestats,nshake_next,nshakebonds
      integer shakewhichbondcoeff,shakeableangle,shakeableanglebond
      real*8 shaketol,shakeanglebond
      integer, allocatable :: shakeablebond(:)

      integer, allocatable :: shakegroup(:)
      integer, allocatable :: shakepartner(:,:),shakebondtype(:,:)

      integer nshakelocal
      integer, allocatable :: shakesize(:),shakeatom(:,:)
      real*8, allocatable :: shakebondlen(:,:)
      real*8, allocatable :: xshake(:,:)

c -------------------------------------------------------------------------
c *** minimizer

c optstyle = style of minimizer, 1 = Hessian-free truncated Newton
c optflag = output minimizer iteration info every this many steps
c opt_max_iters,opt_max_fns = max # of iterations and function evals for min
c opt_stop_tol = stopping tolerance

      integer optstyle,optflag
      integer opt_max_iters,opt_max_fns
      real*8 opt_stop_tol

c -------------------------------------------------------------------------
c *** output

c noutput_next = next timestep to do any kind of output on
c trueflag = flag for whether true box flags are included on input/output

c thermostyle = style of thermo output, little to lots
c nthermo,nthermo_next = how-often,when-next to do thermo output
c ndumpatom,ndumpatom_prev,ndumpatom_next = how-often,when-next to dump atoms
c ndumpvel,ndumpvel_prev,ndumpvel_next = how-often,when-next to dump vels
c ndumpforce,ndumpforce_prev,ndumpforce_next = how-often,when-next dump forces
c dumpatomfileflag,dumpvelfileflag,dumpforcefileflag = flags for whether
c   dump files are currently open or closed
c nrestart,nrestart_next = how-often,when-next to write a restart file
c restartstyle = (1) append timestep to file, (2) toggle between 2 files
c restartlast = which file (1 or 2) was used on last restart output
c restart_version = what version to expect when reading in a restart file

      integer noutput_next
      integer trueflag

      integer thermostyle
      integer nthermo,nthermo_next
      integer ndumpatom,ndumpatom_prev,ndumpatom_next
      integer ndumpvel,ndumpvel_prev,ndumpvel_next
      integer ndumpforce,ndumpforce_prev,ndumpforce_next
      integer dumpatomfileflag,dumpvelfileflag,dumpforcefileflag
      integer nrestart,nrestart_next,restartstyle,restartlast
      integer restart_version

c files for various I/O operations

      character*80 datafile
      character*80 dumpatomfile
      character*80 dumpvelfile
      character*80 dumpforcefile
      character*80 restart_in
      character*80 restart_out,restart_out1,restart_out2

c -------------------------------------------------------------------------
c *** diagnostics

c numdiag = # of diagnostic routines defined by user in diagnostic.f
c diagnames = name of each diag routine as specified in input script

c ndiag_next = next timestep to call any diagnostic routine on
c ndiag = how often to call each diag routine
c diagprev,diagnext = previous/next timestep to call each diag routine
c diagcall = flag for whether each diag routine should be called this timestep

c diagnparams = # of params to pass to each diag routine
c diagparam(5,n) = parameters to pass to each diag routine
c diagfileflag = whether file for each diag routine is open or closed
c diagfile = filename for each diag routine

      integer numdiag
      character*16 diagnames(maxdiag)

      integer ndiag_next
      integer ndiag(maxdiag)
      integer diagprev(maxdiag),diagnext(maxdiag)
      integer diagcall(maxdiag)

      integer diagnparams(maxdiag)
      real*8 diagparam(5,maxdiag)
      integer diagfileflag(maxdiag)
      character*80 diagfile(maxdiag)

c -------------------------------------------------------------------------
c *** timers

c CPU timers for various operations within code

      real*8 time_total
      real*8 time_loop
      real*8 time_current
      real*8 time_nonbond
      real*8 time_bond
      real*8 time_angle
      real*8 time_dihedral
      real*8 time_improper
      real*8 time_comm
      real*8 time_fcomm
      real*8 time_exch
      real*8 time_io
      real*8 time_shake
      real*8 time_neigh1
      real*8 time_neigh2
      real*8 time_long
      real*8 time_rho
      real*8 time_poiss
      real*8 time_field
      real*8 time_other

c -------------------------------------------------------------------------

      end module
