SHELL=/bin/sh

.SILENT:

MESHOPT3D = ../../../src/programs/meshopt3d 

test: rand1.vmesh.out

rand1.vmesh.out:
	$(MESHOPT3D) -i rand1.vmesh -O wiwsmf5mf5rmf5mf5

rand2.vmesh.out:
	$(MESHOPT3D) -i rand2.vmesh -O wiwsmf5mf5rmf5mf5

all: rand1.vmesh.out rand2.vmesh.out

clean:
	-rm -f *.msg *.qual *.out *.vtk.surf

distclean:
	-rm -f *~ *.rej *.orig \#* *.qual *.msg *.out *.vtk *.vmesh.*

regress:
	rm $(TARGET).vmesh
	$(MAKE) $(TARGET).vmesh
	diff $(TARGET).vmesh ref-output/$(TARGET).vmesh
