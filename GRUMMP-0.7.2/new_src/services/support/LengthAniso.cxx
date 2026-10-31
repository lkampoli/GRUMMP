#include <set>
#include <queue>
#include <utility>
#include "GR_LengthAniso.h"
#include "GR_misc.h"
#include "GR_Mesh.h"
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

LengthAniso::LengthAniso(const eLengthType& lengthType, double refine,
			 double grading, double minLength) :
    m_lengthType(lengthType), m_refine(refine), m_grade(grading), m_minLength(
	minLength), m_lengthFunc(NULL), m_lengthMesh(NULL), m_BBoxTree(
    NULL), m_lengthCheckQueue(), m_neigh()
{
  assert(m_refine > 0.);
  assert(m_grade > 0.);
  if (m_minLength <= 0.)
    m_minLength = LARGE_DBL;

  switch (m_lengthType)
    {
    case eFunc:
      break;
    case eNeigh:
      break;
    case eMesh:
      break;
    default:
      assert(0);
      break;
    }
}

LengthAniso::~LengthAniso()
{
  if (m_lengthMesh)
    delete m_lengthMesh;
  m_lengthMesh = NULL;
}

void
LengthAniso::initLength(const Mesh* const mesh)
{
  int num_verts = mesh->getNumVerts();
  Vert* vertex;

  for (int i = 0; i < num_verts; i++) {
    vertex = mesh->getVert(i);
    if (vertex->isDeleted())
      continue;
    assert(vertex->isValid());

    switch (m_lengthType)
      {
      case eFunc:
	initLengthFunc(vertex);
	break;
      case eNeigh:
	initLengthNeigh(vertex);
	break;
      case eMesh:
	initLengthMesh(vertex);
	break;
      default:
	assert(0);
	break;
      }

    m_lengthCheckQueue.push(vertex);
  }
  setGradedLength();
}

void
LengthAniso::initLengthFunc(Vert* const vertex)
{
  setLengthScaleFunc(vertex);
}

void
LengthAniso::initLengthNeigh(Vert* const vertex)
{
  setLengthScaleNeigh(vertex);
}

void
LengthAniso::initLengthMesh(Vert* const vertex)
{
  setLengthScaleMesh(vertex);
}

void
LengthAniso::setLengthScale(const Mesh* const mesh)
{
  int num_verts = mesh->getNumVerts();
  Vert* vertex;

  for (int i = 0; i < num_verts; i++) {
    vertex = mesh->getVert(i);
    if (vertex->isDeleted())
      continue;
    assert(vertex->isValid());

    switch (m_lengthType)
      {
      case eFunc:
	setLengthScaleFunc(vertex);
	break;
      case eNeigh:
	setLengthScaleNeigh(vertex);
	break;
      case eMesh:
	setLengthScaleMesh(vertex);
	break;
      default:
	assert(0);
	break;
      }

    m_lengthCheckQueue.push(vertex);
  }
  setGradedLength();
}

void
LengthAniso::setLengthScale(Vert* const vertex)
{
  assert(vertex->isValid());
  assert(!vertex->isDeleted());

  switch (m_lengthType)
    {
    case eFunc:
      setLengthScaleFunc(vertex);
      break;
    case eNeigh:
      setLengthScaleNeigh(vertex);
      break;
    case eMesh:
      setLengthScaleMesh(vertex);
      break;
    default:
      assert(0);
      break;
    }

  m_lengthCheckQueue.push(vertex);
  setGradedLength();
}

void
LengthAniso::setLengthScaleFunc(Vert* const vertex)
{
  double adMetric[5];
  evaluateLengthFromFunc(vertex, adMetric);

  vertex->setMetric(adMetric);
}

void
LengthAniso::setLengthScaleNeigh(Vert* const vertex)
{
  if (vertex->getNumFaces() == 0) {
    // This vertex is not connected into the mesh.  Set a dummy value.
    vertex->setLengthScale(LARGE_DBL);
    return;
  }
  getNeighborhood(vertex);

  assert(!m_neigh.verts.empty());
  assert(!m_neigh.cells.empty());

  double adMetric[5];
  evaluateLengthFromNeighbors(vertex, adMetric, m_neigh.verts);

  vertex->setMetric(adMetric);
}

void
LengthAniso::setLengthScaleMesh(Vert* const vertex)
{
  double adMetric[5];
  evaluateLengthFromMesh(vertex, adMetric);

  vertex->setMetric(adMetric);
}

// Still need to figure out if/how we should grade metric.
void
LengthAniso::setGradedLength()
{
  /*
   Vert *vertex = NULL, *neigh_vertex = NULL;
   double length_scale, lip_length, neigh_ls, neigh_dist;

   while (!m_lengthCheckQueue.empty()) {
   lip_length = LARGE_DBL;

   vertex = m_lengthCheckQueue.front();
   m_lengthCheckQueue.pop();

   getNeighborhood(vertex);

   length_scale = vertex->getLengthScale();

   set<Vert*>::iterator it = m_neigh.verts.begin(), it_end =
   m_neigh.verts.end();
   for (; it != it_end; ++it) {
   neigh_vertex = *it;
   neigh_ls = neigh_vertex->getLengthScale();
   neigh_dist = calcDistanceBetween(vertex, neigh_vertex);

   assert(iFuzzyComp(neigh_ls, 0.) == 1);
   assert(iFuzzyComp(neigh_dist, 0.) == 1);

   lip_length = std::min(lip_length, neigh_ls + (neigh_dist / m_grade));
   }

   assert(iFuzzyComp(m_minLength, 0.) == 1);
   assert(iFuzzyComp(lip_length, 0.) == 1);

   length_scale =
   MIN3( m_minLength < 0. ? LARGE_DBL : m_minLength, lip_length, length_scale);
   vertex->setLengthScale(length_scale);
   }
   */
}

void
LengthAniso::provideFunc(void
(*func)(const double adLoc[], double adRetVal[]))
{
  assert(m_lengthType == eFunc);
  m_lengthFunc = func;
}

void
LengthAniso::evaluateLengthFromFunc(Vert* const vertex, double* adMetric) const
{
  m_lengthFunc(vertex->getCoords(), adMetric);
}

void
LengthAniso::getNeighborhood(Vert* const vertex)
{

  m_neigh.my_vert = vertex;
  findNeighborhoodInfo(vertex, m_neigh.cells, m_neigh.verts);
}

