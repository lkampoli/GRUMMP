#include <set>
#include <queue>
#include <utility>
#include "GR_LengthAniso3D.h"
#include "GR_misc.h"
#include "GR_VolMesh.h"
#include "GR_Entity.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_GRPoint.h"
#include "GR_GRCurve.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"
#include "RTree.hpp"
#include "GR_ADT.h"

using std::queue;
using std::set;
using std::pair;
using std::make_pair;

// requires simplicial mesh
void
LengthAniso3D::provideMesh(Mesh* /*pM3D*/, char* /*strLenFileName*/)
{

}

void
LengthAniso3D::evaluateLengthFromNeighbors(
    Vert* const /*vertex*/, double* /*adMetric*/,
    const std::set<Vert*>& /*neigh_verts*/) const
{

}

/*
 * Currently copied from scat3d. Should refactor so that this code exists only in one place.
 */
void
LengthAniso3D::evaluateLengthFromMesh(Vert* const /*vertex*/,
				      double* /*adMetric*/) const
{

}
