#ifndef IMESH_F_H
#define IMESH_F_H

#define IMESH_VERSION_MAJOR ITAPS_VERSION_MAJOR
#define IMESH_VERSION_MINOR ITAPS_VERSION_MINOR
#define IMESH_VERSION_PATCH ITAPS_VERSION_PATCH
#define IMESH_MAJOR_VERSION IMESH_VERSION_MAJOR
#define IMESH_MINOR_VERSION IMESH_VERSION_MINOR
#define IMESH_PATCH_VERSION IMESH_VERSION_PATCH
#define IMESH_VERSION_GE(Maj,Min,Pat) ITAPS_VERSION_GE(Maj,Min,Pat)
#define IMESH_VERSION_STRING ITAPS_VERSION_STRING_(iMesh)
#define IMESH_VERSION_TAG ITAPS_VERSION_TAG_(iMesh)
#define IMESH_NEW_MESH_NAME__(A,B,C) A##_##B##_##C
#define IMESH_NEW_MESH_NAME_(A,B,C) IMESH_NEW_MESH_NAME__(A,B,C)
#define IMESH_NEW_MESH_NAME(A) IMESH_NEW_MESH_NAME_(A,IMESH_VERSION_MAJOR,IMESH_VERSION_MINOR)

#define iMesh_Instance IBASE_HANDLE_T

#endif

integer iMesh_EntityTopology_MIN
integer iMesh_POINT
integer iMesh_LINE_SEGMENT
integer iMesh_POLYGON
integer iMesh_TRIANGLE
integer iMesh_QUADRILATERAL
integer iMesh_POLYHEDRON
integer iMesh_TETRAHEDRON
integer iMesh_HEXAHEDRON
integer iMesh_PRISM
integer iMesh_PYRAMID
integer iMesh_SEPTAHEDRON
integer iMesh_ALL_TOPOLOGIES
integer iMesh_EntityTopology_MAX

parameter(iMesh_EntityTopology_MIN=0)
parameter(iMesh_POINT=0)
parameter(iMesh_LINE_SEGMENT=1)
parameter(iMesh_POLYGON=2)
parameter(iMesh_TRIANGLE=3)
parameter(iMesh_QUADRILATERAL=4)
parameter(iMesh_POLYHEDRON=5)
parameter(iMesh_TETRAHEDRON=6)
parameter(iMesh_HEXAHEDRON=7)
parameter(iMesh_PRISM=8)
parameter(iMesh_PYRAMID=9)
parameter(iMesh_SEPTAHEDRON=10)
parameter(iMesh_ALL_TOPOLOGIES=11)
parameter(iMesh_EntityTopology_MAX=11)

