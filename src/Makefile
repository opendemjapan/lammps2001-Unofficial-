# Multiple-machine Makefile

SHELL = /bin/sh
.IGNORE:

# Files

SRCF =	lammps.f communicate.f coord2bin.f diagnostic.f dump_atom.f \
	dump_force.f dump_vel.f ensemble.f error.f ewald.f ewald_coeff.f \
	finish.f fix.f force.f force_bond.f force_class2.f force_many.f \
	force_respa.f global.f initialize.f input.f integrate.f \
	integrate_respa.f lapack.f memory.f mesh_3d.f min_support.f \
	minimize.f misc.f mpi.f neighbor.f pbc.f pppm.f pppm_coeff.f \
	pppm_remap.f random.f read_data.f read_restart.f \
	scan_data.f setup.f shake.f special.f start.f stats.f string.f \
	thermo.f velocity.f write_restart.f

SRCC = 	fft_3d.c

# Definitions

ROOT =		lmp
EXE =		$(ROOT)_$@
SRC =		$(SRCF) $(SRCC)
OBJ = 		$(SRCF:.f=.o) $(SRCC:.c=.o)

# Help

help:
	@echo 'Type "make target" where target is one of:'
	@echo '      serial	(for desktop w/out MPI)'
	@echo '      linux	(for desktop Linux box)'
	@echo '      tflop	(for Intel Tflop)'
	@echo '      alaska	(for CPlant - alaska)'
	@echo '      ross	(for CPlant - ross)'
	@echo '      liberty	(for Institutional cluster)'

# Targets

serial linux tflop alaska ross liberty:
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) Makefile.depend Obj_$@
	@cp Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE) "OBJ = $(OBJ)" "EXE = ../$(EXE)" ../$(EXE)
	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) Makefile*; fi

depend:
	Depend.pl $(SRC)

# Cleans

clean_serial:
	rm -r Obj_serial

clean_linux:
	rm -r Obj_linux

clean_tflop:
	rm -r Obj_tflop

clean_alaska:
	rm -r Obj_alaska

clean_ross:
	rm -r Obj_ross

clean_liberty:
	rm -r Obj_liberty
