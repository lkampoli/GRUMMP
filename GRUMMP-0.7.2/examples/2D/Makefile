.SUFFIXES: .poly .mesh .bdry

.SILENT:

BDRY_FILES=0012-spline.bdry \
0012.bdry \
6IMR.test.bdry \
6IMR.test2.bdry \
af4-spline.bdry \
af4.bdry \
aniso1.bdry \
aniso2.bdry \
aniso3.bdry \
aniso4.bdry \
aniso5.bdry \
aniso6.bdry \
annulus.bdry \
anslab-round.bdry \
blob.bdry \
branch.bdry \
circle.bdry \
cubicbezier.bdry \
curvestar.bdry \
cyl.bdry \
cyl2.bdry \
destruct.bdry \
easyspline.bdry \
face.bdry \
front-step.bdry \
habs.bdry \
hexagon.bdry \
hexagon2.bdry \
hex-flow.bdry \
longarc.bdry \
mapleleaf-circle.bdry \
multi-region.bdry \
nozzle.bdry \
plate-hole.bdry \
quad.bdry \
quarter-sphere-peel.bdry \
saddle.bdry \
scram.bdry \
sector.bdry \
shortarc.bdry \
square-circhole-arcside.bdry \
square-circhole-bezier.bdry \
square-hole.bdry \
square.bdry \
square2.bdry

MESH_FILES=$(BDRY_FILES:.bdry=.mesh)
COARSEN_FILES=$(MESH_FILES:.mesh=.mesh.out.mesh)
REGRESS_FILES=$(MESH_FILES)

all: $(MESH_FILES)

test: square-circhole-bezier.mesh shortarc.mesh longarc.mesh saddle.mesh

clean:
	-rm -f *~ *.rej *.orig \#* 

distclean: clean
	-rm -f *.msg *.qual *.mesh *.bpatch *.vtk *.fea *.grmesh

.bdry.mesh:
	../../src/programs/tri -i $* -r 8 -g 8

%.mesh.out.mesh:
	../../src/programs/tri -i $* -r 4 -g 4 -c 1
	../../src/programs/coarsen2d -i $*.mesh -b $*

coarsen: 
	-rm -f $(MESH_FILES) $(MESH_FILES:.mesh=.mesh.*);
	$(MAKE) -k -i $(COARSEN_FILES);
	for i in $(BDRY_FILES) ; do \
		if ! test -e `basename $$i .bdry`.mesh; then \
			echo "tri: No output for $$i."; fi; done 
	for i in $(MESH_FILES) ; do \
		if ! test -e `basename $$i .mesh`.mesh.out.mesh; then \
			echo "coarsen2d: No output for $$i."; fi; done 
regress:
	for i in $(REGRESS_FILES) ; do \
	  $(MAKE) regress-single TARGET=$$i OUTPUT="> /dev/null" ; \
	done

regress-single:
	echo -n "Regression testing 2D/$(TARGET)..."
	rm -f $(TARGET)
	$(MAKE) $(TARGET) $(OUTPUT) || (echo "---MESHING FAILED---" && false)
	diff $(TARGET) ref-output/$(TARGET) $(OUTPUT) || (echo "---DIFF FAILED---" && false)
	echo "success"
