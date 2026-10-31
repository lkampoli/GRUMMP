include ../../../src/conf/default-rules
include /home/cfog/GRUMMP/iMesh-Defs.inc

SRC=vTestMeshInterface.cc \
    vTestEntityInterface.cc \
    vTestArrInterface.cc \
    vTestModifyInterface.cc \
    vTestEntSetInterface.cc \
    vTestMeshWithSets.cc \
    vTestEntityWithSets.cc \
    vTestArrayWithSets.cc \
    vTestEntTagInterface.cc \
    vTestArrTagInterface.cc \
    vTestSetRelationInterface.cc \
    vTestBoolOpsInterface.cc \
    vTestSetTagInterface.cc \
    vTestLowerDimensionedCoordinates.cc

OBJ=$(SRC:.cc=.o)

# Temporarily maintain backward compat. with older iMesh-Defs.inc varnames
# Eventually, we should remove reference to IMESH_INCLUDES
TEST_CPPFLAGS=${IMESH_CPPFLAGS} ${IMESH_INCLUDES}

default: iMesh_unitTest

iMesh_unitTest: iMesh_unitTest.o $(OBJ) $(IMESH_DEPENDENCIES)
	echo "Linking iMesh_unitTest"
	mpic++ -fPIC -fPIC  -fopenmp -O2 -fopenmp $(EXTRAFLAGS) -fopenmp -o iMesh_unitTest iMesh_unitTest.o $(OBJ) ${IMESH_LIBS} ${SYS_LIBS}

%.o: %.cc iMesh_unitTest.hh iMesh_unitTest_config.h
	mpic++ -fPIC  -fPIC  -fopenmp -O2 -fopenmp $(EXTRAFLAGS) ${TEST_CPPFLAGS} -c $<

check: test

test: iMesh_unitTest
	./iMesh_unitTest --silent annulus.vtk
	./iMesh_unitTest --silent globe_mats.vtk

distclean: localdistclean

localdistclean:
	-rm -f Makefile *~ config.* iMesh_unitTest iMesh_unitTest_config.h

# DO NOT DELETE

iMesh_unitTest.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
iMesh_unitTest.o: /home/cfog/GRUMMP/include/iMesh.h
iMesh_unitTest.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestArrInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestArrInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestArrInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestArrTagInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestArrTagInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestArrTagInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestArrayWithSets.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestArrayWithSets.o: /home/cfog/GRUMMP/include/iMesh.h
vTestArrayWithSets.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestBoolOpsInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestBoolOpsInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestBoolOpsInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestEntSetInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestEntSetInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestEntSetInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestEntTagInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestEntTagInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestEntTagInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestEntityInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestEntityInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestEntityInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestEntityWithSets.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestEntityWithSets.o: /home/cfog/GRUMMP/include/iMesh.h
vTestEntityWithSets.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestLowerDimensionedCoordinates.o: iMesh_unitTest.hh
vTestLowerDimensionedCoordinates.o: /home/cfog/GRUMMP/include/iBase.h
vTestLowerDimensionedCoordinates.o: /home/cfog/GRUMMP/include/iMesh.h
vTestLowerDimensionedCoordinates.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestMeshInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestMeshInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestMeshInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestMeshWithSets.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestMeshWithSets.o: /home/cfog/GRUMMP/include/iMesh.h
vTestMeshWithSets.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestModifyInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestModifyInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestModifyInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestSetRelationInterface.o: iMesh_unitTest.hh
vTestSetRelationInterface.o: /home/cfog/GRUMMP/include/iBase.h
vTestSetRelationInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestSetRelationInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
vTestSetTagInterface.o: iMesh_unitTest.hh /home/cfog/GRUMMP/include/iBase.h
vTestSetTagInterface.o: /home/cfog/GRUMMP/include/iMesh.h
vTestSetTagInterface.o: /home/cfog/GRUMMP/include/iMesh_protos.h
