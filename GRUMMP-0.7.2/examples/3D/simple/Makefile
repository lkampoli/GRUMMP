SHELL=/bin/sh

.SUFFIXES: .smesh .bdry .vmesh .grvmesh

.SILENT:

all: smesh bdry multi-domain

test: hex3.vmesh cube-clip.vmesh hex-dual.vmesh

clean:
	-rm -f *~ *.rej *.orig \#* 

distclean: clean
	-rm -f *.msg *.qual *.vmesh *pie* *vtk* *.grvmesh

SMESH_FILES = badprism.smesh \
doorroom.smesh \
doubtet.smesh \
faceknife.smesh \
goodprism.smesh \
hex.smesh \
hex1.smesh \
hex2.smesh \
hex3.smesh \
indpint.smesh \
knife.smesh \
knife2.smesh \
knife3.smesh \
normprism.smesh \
octa.smesh \
patho1.smesh \
pint.smesh \
skinny.smesh \
sliver.smesh \
tet.smesh \
tet2.smesh \
tint.smesh \
tint2.smesh 

BDRY_FILES = cube-clip.bdry \
cutout-cube.bdry \
hex-bdry.bdry \
hex-cavity.bdry \
hex-cavity2.bdry \
hex-tunnel.bdry \
pyramid.bdry \
tet3.bdry \
tet4.bdry \
tunnel2.bdry \
tunnel3.bdry \
seg-split.bdry 

MULTI_DOMAIN = hex-dual.bdry \
hex-embed.bdry \
hex-embed2.bdry \
wagner.bdry

TETRA = ../../../src/programs/tetra 

smesh: $(SMESH_FILES:.smesh=.vmesh)

multi-domain: $(MULTI_DOMAIN:.bdry=.grvmesh)

bdry: $(BDRY_FILES:.bdry=.vmesh)

.smesh.vmesh:
	./../../../src/programs/tetra -i $* 

.bdry.vmesh:
	./../../../src/programs/tetra -i $* 

.bdry.grvmesh:
	./../../../src/programs/tetra -i $* 


MESH_FILES=$(SMESH_FILES:.smesh=.vmesh) $(MULTI_DOMAIN:.bdry=.vmesh) $(BDRY_FILES:.bdry=.vmesh)
REGRESS_FILES=$(MESH_FILES)
#REGRESS_FILES=badprism.vmesh cube-clip.vmesh cutout-cube.vmesh doorroom.vmesh doubtet.vmesh faceknife.vmesh goodprism.vmesh hex1.vmesh hex2.vmesh hex3.vmesh hex-bdry.vmesh hex-cavity2.vmesh hex-cavity.vmesh hex-dual.vmesh hex-embed2.vmesh hex-embed.vmesh hex-tunnel.vmesh hex.vmesh indpint.vmesh knife2.vmesh knife3.vmesh knife.vmesh normprism.vmesh octa.vmesh patho1.vmesh pint.vmesh pyramid.vmesh seg-split.vmesh skinny.vmesh sliver.vmesh sphere2.vmesh sphere.vmesh tet2.vmesh tet3.vmesh tet4.vmesh tet.vmesh tint2.vmesh tint.vmesh tunnel2.vmesh tunnel3.vmesh wagner.vmesh


regress:
	for i in $(REGRESS_FILES) ; do \
	  $(MAKE) regress-single TARGET=$$i OUTPUT="> /dev/null" ; \
	done

regress-single:
	echo -n "Regression testing 3D/simple/$(TARGET)..."
	rm -f $(TARGET)
	$(MAKE) $(TARGET) $(OUTPUT) || (echo "---MESHING FAILED---" && false)
	diff $(TARGET) ref-output/$(TARGET) $(OUTPUT) || (echo "---DIFF FAILED---" && false)
	echo "success"
