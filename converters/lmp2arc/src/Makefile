
OBJS = lmp2arc.o ReadCarFile.o ProcessPosFile.o unwrap_molecules.o WriteArcFrame.o

lmp2arc.exe: $(OBJS)
	cc -o lmp2arc.exe $(OBJS)

lmp2arc.o: lmp2arc.c lmp2.h
ReadCarFile.o: ReadCarFile.c lmp2.h
ProcessPosFile.o: ProcessPosFile.c lmp2.h
WriteArcFrame.o: WriteArcFrame.c lmp2.h
unwrap_molecules.o: unwrap_molecules.c lmp2.h

.c.o:
	cc -O2 -c $<

clean: 
	rm lmp2arc.exe *.o
