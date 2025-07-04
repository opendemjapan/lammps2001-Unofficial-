# Makefile for serial tools - modify for your compilers

.IGNORE:

#
# Targets
#

all:	peek_restart replicate restart2data setup_chain setup_flow_2d setup_lj

peek_restart:	peek_restart.o
	f90 -n32 peek_restart.o -o peek_restart

replicate:	replicate.o
	cc replicate.o -o replicate

restart2data:	restart2data.o
	f90 -n32 restart2data.o -o restart2data

setup_chain:	setup_chain.o
	f90 -n32 setup_chain.o -o setup_chain

setup_flow_2d:	setup_flow_2d.o
	f90 -n32 setup_flow_2d.o -o setup_flow_2d 

setup_lj:	setup_lj.o
	f90 -n32 setup_lj.o -o setup_lj 

clean:
	rm peek_restart replicate restart2data
	rm setup_chain setup_flow_2d setup_lj
	rm *.o

#
# Rules
#

.f.o:	
	f90 -O -n32 -c $<

.c.o:
	cc -O -n32 -c $<
