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

c ----------------------------------------------------------------------
c allocate atom memory by estimating sizes then adding extra

      subroutine atom_memory
      use global
      implicit none

c local variables

      real*8 sizex,sizey,sizez,volume,density,rneigh
      integer ntotal,iperatom

 900  format(4f10.2)
 901  format(5i10)

c cutmemory = cutoff distance used for allocating ghost memory
c           = cutmax + skin (if "maximum cutoff" command has been used)
c           = max(cutlj,cutcoul) + skin (otherwise)

      if (cutmax > 0.0) then
        cutmemory = cutmax + skin
      else
        cutmemory = max(cutlj,cutcoul) + skin
      endif

c nlocal = my fraction of total atoms
c          just an estimate, will be set explicitly by read_data/read_restart
c volume = my sub-box extended by neighbor cutoff in all dimensions
c nghost = assume constant density in extra neighbor volume
c rneigh = # of neighbors/atom with Newton 1/2 factor

      nlocal = natoms/nprocs
      density = natoms / (xprd*yprd*zprd)

      sizex = xprd/pgrid(1) + 2.0*cutmemory
      sizey = yprd/pgrid(2) + 2.0*cutmemory
      sizez = zprd/pgrid(3) + 2.0*cutmemory
      volume = sizex * sizey * sizez
      ntotal = volume * density
      nghost = ntotal - nlocal + 1

      volume = 4.0/3.0 * 3.1415926 * cutmemory**3
      rneigh = volume * density
      rneigh = 0.5 * rneigh

c pad amounts by extra factors
c insure at least min amounts are allocated
c this enables small problems to be easily run without overflow

      maxown = nlocal * extra_own
      maxghost = nghost * extra_ghost
      maxneigh = rneigh * nlocal * extra_neigh
      maxbuf = nghost * extra_buf

      maxown = max(minown,maxown)
      maxghost = max(minghost,maxghost)
      maxneigh = max(minneigh,maxneigh)
      maxbuf = max(minbuf,maxbuf)

      maxatom = maxown + maxghost
      maxatomp1 = maxatom + 1

c set bond neighbor lists to maximum possible size
c on one proc, set list sizes exactly

      maxbondlocal = maxbondper * maxown
      maxanglelocal = maxangleper * maxown
      maxdihedlocal = maxdihedper * maxown
      maximprolocal = maximproper * maxown

      if (nprocs == 1) then
        maxbondlocal = nbonds
        maxanglelocal = nangles
        maxdihedlocal = ndihedrals
        maximprolocal = nimpropers
      endif

c print out results

      if (node.eq.0) then
        write (6,*) 'Extra params: Extra own,ghost,neigh,buf'
        write (6,900) extra_own,extra_ghost,extra_neigh,extra_buf
        write (6,*) 'Atom arrays: Max own,ghost,atom,buf'
        write (6,901) maxown,maxghost,maxatom,maxbuf
        write (6,*) 'Neighbor arrays: ',
     $       'Max nonbond,bond,angle,dihedral,improper'
        write (6,901) maxneigh,maxbondlocal,maxanglelocal,
     $       maxdihedlocal,maximprolocal
        write (1,*) 'Extra params: Extra own,ghost,neigh,buf'
        write (1,900) extra_own,extra_ghost,extra_neigh,extra_buf
        write (1,*) 'Atom arrays: Max own,ghost,atom,buf'
        write (1,901) maxown,maxghost,maxatom,maxbuf
        write (1,*) 'Neighbor arrays: ',
     $       'Max nonbond,bond,angle,dihedral,improper'
        write (1,901) maxneigh,maxbondlocal,maxanglelocal,
     $       maxdihedlocal,maximprolocal
      endif

c allocate atom arrays

      allocate(x(3,maxatom))
      allocate(v(3,maxown))
      allocate(f(3,maxatom))

      allocate(q(maxatom))
      allocate(xhold(3,maxown))

      allocate(tag(maxatom))
      allocate(type(maxatom))
      allocate(molecule(maxown))
      allocate(true(maxown))

      allocate(fix(maxown))
      allocate(velflag(maxown))

      allocate(numbond(maxown))
      allocate(numangle(maxown))
      allocate(numdihed(maxown))
      allocate(numimpro(maxown))

      allocate(num1bond(maxown))
      allocate(num2bond(maxown))
      allocate(num3bond(maxown))

c allocate nonbond neighbor lists

      allocate(nlist(maxneigh+1000))
      allocate(nnfirst(maxown))
      allocate(nnlast(maxown))
      allocate(bin(maxatom))

c allocate communication arrays and buffers

      allocate(slist(maxghost))

      allocate(ibuf1(maxbuf+1000))
      allocate(ibuf2(maxbuf+1000))
      allocate(buf1(maxbuf+1000))
      allocate(buf2(maxbuf+1000))

c allocate global->local ptr = size of global problem

      allocate(localptr(natoms))

c allocate atom masses

      allocate(mass(ntypes))

c allocate bonded topology arrays stored per atom

      allocate(bondtype(maxbondper,maxown))
      allocate(bondatom1(maxbondper,maxown))
      allocate(bondatom2(maxbondper,maxown))

      allocate(angletype(maxangleper,maxown))
      allocate(angleatom1(maxangleper,maxown))
      allocate(angleatom2(maxangleper,maxown))
      allocate(angleatom3(maxangleper,maxown))

      allocate(dihedtype(maxdihedper,maxown))
      allocate(dihedatom1(maxdihedper,maxown))
      allocate(dihedatom2(maxdihedper,maxown))
      allocate(dihedatom3(maxdihedper,maxown))
      allocate(dihedatom4(maxdihedper,maxown))

      allocate(improtype(maximproper,maxown))
      allocate(improatom1(maximproper,maxown))
      allocate(improatom2(maximproper,maxown))
      allocate(improatom3(maximproper,maxown))
      allocate(improatom4(maximproper,maxown))

c allocate bonded neighbor lists

      allocate(bondlist(3,maxbondlocal))
      allocate(anglelist(4,maxanglelocal))
      allocate(dihedlist(5,maxdihedlocal))
      allocate(improlist(5,maximprolocal))

c allocate force field arrays

      allocate(noncoeff1(ntypes,ntypes))
      allocate(noncoeff2(ntypes,ntypes))
      allocate(noncoeff3(ntypes,ntypes))
      allocate(noncoeff4(ntypes,ntypes))
      allocate(noncoeff14_1(ntypes,ntypes))
      allocate(noncoeff14_2(ntypes,ntypes))
      allocate(nontypeflag(ntypes,ntypes))

      allocate(bondcoeff(5,nbondtypes))
      allocate(bondtypeflag(nbondtypes))
      allocate(shakeablebond(nbondtypes))

      allocate(anglecoeff(4,nangletypes))
      allocate(angletypeflag(nangletypes))

      allocate(dihedcoeff(6,ndihedtypes))
      allocate(dihedtypeflag(ndihedtypes))

      allocate(improcoeff(3,nimprotypes))
      allocate(improtypeflag(nimprotypes))

      allocate(bondbondcoeff(3,nangletypes))
      allocate(bondanglecoeff(4,nangletypes))
      allocate(midbondtorsioncoeff(4,ndihedtypes))
      allocate(endbondtorsioncoeff(8,ndihedtypes))
      allocate(angletorsioncoeff(8,ndihedtypes))
      allocate(angleangletorsioncoeff(3,ndihedtypes))
      allocate(bondbond13coeff(3,ndihedtypes))
      allocate(angleanglecoeff(6,nimprotypes))

      return
      end


c ----------------------------------------------------------------------
c allocate arrays for all nonbond interactions

      subroutine nonbond_memory
      use global
      implicit none

      allocate(cutforcesq(ntypes,ntypes))
      allocate(cutneighsq(ntypes,ntypes))

      if (nonstyle == 0) then
        allocate(cutljsq(ntypes,ntypes))
      else if (nonstyle == 1) then
        allocate(cutljsq(ntypes,ntypes))
        allocate(lj1(ntypes,ntypes))
        allocate(lj2(ntypes,ntypes))
        allocate(lj3(ntypes,ntypes))
        allocate(lj4(ntypes,ntypes))
        allocate(offset(ntypes,ntypes))
      else if (nonstyle == 2) then
        allocate(cutljsq(ntypes,ntypes))
        allocate(lj1(ntypes,ntypes))
        allocate(lj2(ntypes,ntypes))
        allocate(lj3(ntypes,ntypes))
        allocate(lj4(ntypes,ntypes))
        allocate(cutljinner(ntypes,ntypes))
        allocate(cutljinnersq(ntypes,ntypes))
        allocate(ljsw0(ntypes,ntypes))
        allocate(ljsw1(ntypes,ntypes))
        allocate(ljsw2(ntypes,ntypes))
        allocate(ljsw3(ntypes,ntypes))
        allocate(ljsw4(ntypes,ntypes))
        allocate(offset(ntypes,ntypes))
      else if (nonstyle == 3) then
        allocate(cutljsq(ntypes,ntypes))
        allocate(lj1(ntypes,ntypes))
        allocate(lj2(ntypes,ntypes))
        allocate(lj3(ntypes,ntypes))
        allocate(lj4(ntypes,ntypes))
        allocate(lj5(ntypes,ntypes))
        allocate(offset(ntypes,ntypes))
      else if (nonstyle == 4) then
        allocate(cutljsq(ntypes,ntypes))
        allocate(lj1(ntypes,ntypes))
        allocate(lj2(ntypes,ntypes))
        allocate(lj3(ntypes,ntypes))
        allocate(lj4(ntypes,ntypes))
      else if (nonstyle == 5) then
        allocate(cutljsq(ntypes,ntypes))
        allocate(lj1(ntypes,ntypes))
        allocate(lj2(ntypes,ntypes))
        allocate(lj3(ntypes,ntypes))
        allocate(lj4(ntypes,ntypes))
        allocate(offset(ntypes,ntypes))
      else if (nonstyle == 6) then
        allocate(cutljsq(ntypes,ntypes))
        allocate(lj1(ntypes,ntypes))
        allocate(lj2(ntypes,ntypes))
        allocate(lj3(ntypes,ntypes))
        allocate(lj4(ntypes,ntypes))
        allocate(lj14_1(ntypes,ntypes))
        allocate(lj14_2(ntypes,ntypes))
        allocate(lj14_3(ntypes,ntypes))
        allocate(lj14_4(ntypes,ntypes))
      endif

      return
      end


c ----------------------------------------------------------------------
c free arrays for all nonbond interactions

      subroutine nonbond_deallocate
      use global
      implicit none

      deallocate(cutforcesq,cutneighsq)

      if (nonstyle == 0) then
        deallocate(cutljsq)
      else if (nonstyle == 1) then
        deallocate(cutljsq)
        deallocate(lj1,lj2,lj3,lj4)
        deallocate(offset)
      else if (nonstyle == 2) then
        deallocate(cutljsq)
        deallocate(lj1,lj2,lj3,lj4)
        deallocate(cutljinner,cutljinnersq)
        deallocate(ljsw0,ljsw1,ljsw2,ljsw3,ljsw4)
        deallocate(offset)
      else if (nonstyle == 3) then
        deallocate(cutljsq)
        deallocate(lj1,lj2,lj3,lj4,lj5)
        deallocate(offset)
      else if (nonstyle == 4) then
        deallocate(cutljsq)
        deallocate(lj1,lj2,lj3,lj4)
      else if (nonstyle == 5) then
        deallocate(cutljsq)
        deallocate(lj1,lj2,lj3,lj4)
        deallocate(offset)
      else if (nonstyle == 6) then
        deallocate(cutljsq)
        deallocate(lj1,lj2,lj3,lj4)
        deallocate(lj14_1,lj14_2,lj14_3,lj14_4)
      endif

      return
      end

c ----------------------------------------------------------------------
c tally up total memory used for all dynamically-allocated arrays
c done just before a run to report total memory usage
c sizes of arrays are the same on all procs

      subroutine memory_usage
      use global
      use mpi
      implicit none

c local variables

      integer bytes,maxbytes,kmaxall,ierror
      real*8 mbytes

 900  format(' ',a,f8.3)

c assume integer is 4 bytes, real is 8 bytes, complex is 16 bytes

      integer, parameter :: ibyte = 4
      integer, parameter :: rbyte = 8
      integer, parameter :: cbyte = 16

      bytes = 0

c x,v,f arrays - allocated in atom_memory
c q,xhold arrays
c tag,type,molecule,true arrays
c fix,velflag arrays
c numbond,numangle,numdihed,numimpro arrays
c num123bond arrays

      bytes = bytes + (3*maxatom + 3*maxown + 3*maxatom) * rbyte
      bytes = bytes + (maxatom + 3*maxown) * rbyte
      bytes = bytes + (maxatom + maxatom + maxown + maxown) * ibyte
      bytes = bytes + (maxown + maxown) * ibyte
      bytes = bytes + (4*maxown) * ibyte
      bytes = bytes + (3*maxown) * ibyte
      
c nlist,nnfirst,nnlast,bin arrays - allocated in atom_memory
c slist array
c ibuf1,ibuf2,buf1,buf2 arrays
c localptr array

      bytes = bytes + (maxneigh + maxown + maxown + maxatom) * ibyte
      bytes = bytes + (maxghost) * ibyte
      bytes = bytes + (2*maxbuf) * ibyte + (2*maxbuf) * rbyte
      bytes = bytes + (natoms) * ibyte

c mass array - allocated in atom_memory
c bondtype,bondatom12 arrays
c angletype,angleatom123 arrays
c dihedtype,dihedatom1234 arrays
c improtype,improatom1234 arrays

      bytes = bytes + (ntypes) * rbyte
      bytes = bytes + (3*maxbondper*maxown) * ibyte
      bytes = bytes + (4*maxangleper*maxown) * ibyte
      bytes = bytes + (5*maxdihedper*maxown) * ibyte
      bytes = bytes + (5*maximproper*maxown) * ibyte

c bondlist,anglelist,dihedlist,improlist arrays - allocated in atom_memory

      bytes = bytes + (3*maxbondlocal + 4*maxanglelocal) * ibyte
      bytes = bytes + (5*maxdihedlocal + 5*maximprolocal) * ibyte

c noncoeff1234,noncoeff14_12,nontypeflag arrays - allocated in atom_memory
c bondcoeff,bondtypeflag,shakeablebond arrays
c anglecoeff,angletypeflag arrays
c dihedcoeff,dihedtypeflag arrays
c improcoeff,improtypeflag arrays

      bytes = bytes + (6*ntypes*ntypes) * rbyte +
     $     (ntypes*ntypes) * ibyte
      bytes = bytes + (5*nbondtypes) * rbyte + (2*nbondtypes) * ibyte
      bytes = bytes + (4*nangletypes) * rbyte + (nangletypes) * ibyte
      bytes = bytes + (6*ndihedtypes) * rbyte + (ndihedtypes) * ibyte
      bytes = bytes + (3*nimprotypes) * rbyte + (nimprotypes) * ibyte

c class 2 force field arrays - allocated in atom_memory

      bytes = bytes + (3*nangletypes + 4*nangletypes +
     $     4*ndihedtypes + 8*ndihedtypes + 8*ndihedtypes +
     $     3*ndihedtypes + 3*ndihedtypes + 6*nimprotypes) * rbyte

c swap arrays - allocated in setup_comm

      bytes = bytes + (2*nswap) * rbyte + (10*nswap) * ibyte

c neighbor bin array - allocated in setup_neigh

      if (neighstyle == 1)
     $     bytes = bytes + (mbinx*mbiny*mbinz) * ibyte

c nonbond interaction arrays - allocated in nonbond_memory

      bytes = bytes + (2*ntypes*ntypes) * rbyte
      if (nonstyle == 0) then
        bytes = bytes + (ntypes*ntypes) * rbyte
      else if (nonstyle == 1) then
        bytes = bytes + (6*ntypes*ntypes) * rbyte
      else if (nonstyle == 2) then
        bytes = bytes + (13*ntypes*ntypes) * rbyte
      else if (nonstyle == 3) then
        bytes = bytes + (7*ntypes*ntypes) * rbyte
      else if (nonstyle == 4) then
        bytes = bytes + (5*ntypes*ntypes) * rbyte
      else if (nonstyle == 5) then
        bytes = bytes + (6*ntypes*ntypes) * rbyte
      else if (nonstyle == 6) then
        bytes = bytes + (9*ntypes*ntypes) * rbyte
      endif

c fhold array - allocated in start

      if (nrespa == 0 .and. newton == 3)
     $     bytes = bytes + (3*maxatom) * rbyte

c rRESPA arrays - allocated in start

      if (nrespa == 1) bytes = bytes + (4*3*maxown) * rbyte

c specbond array - allocated in special_create

      bytes = bytes + (maxspec*maxown) * ibyte

c SHAKE arrays - allocated in shake_allocate

      if (nshake == 1) then
        bytes = bytes + (8*maxown) * ibyte
        bytes = bytes + (5*maxown) * ibyte
        bytes = bytes + (3*maxown) * rbyte
        bytes = bytes + (3*maxatom) * rbyte
      endif

c Ewald arrays - allocated in ewald_coeff

      if (coulstyle == 3) then
        kmaxall = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax
        bytes = bytes + (3*kmaxall) * ibyte
        bytes = bytes + (10*kmaxall + 3*maxown) * rbyte
        bytes = bytes + (2*kmaxall + 2*kmaxall) * rbyte
        bytes = bytes + (2*(2*kmax+1)*3*maxown) * rbyte
      endif

c PPPM arrays - allocated in pppm_coeff
c assume maxgrid = maxfft

      if (coulstyle == 4) then
        bytes = bytes + (4*maxgrid) * rbyte
        bytes = bytes + (8*maxgrid) * rbyte + (2*maxgrid) * cbyte
        bytes = bytes + (3*maxown) * ibyte + (3*maxown) * rbyte
        bytes = bytes + (3*(nxhi_fft-nxlo_fft+1)) * rbyte
        bytes = bytes + (2*maxpbuf) * rbyte
      endif

c print-out max result across procs

      call mpi_allreduce(bytes,maxbytes,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      mbytes = maxbytes/1024.0/1024.0
      if (node == 0) then
        write (6,900) 'Memory use per processor (MBytes) =',mbytes
        write (1,900) 'Memory use per processor (MBytes) =',mbytes
      endif

      return
      end
