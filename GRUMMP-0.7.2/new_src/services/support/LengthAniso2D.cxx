#include <set>
#include <queue>
#include <utility>
#include "GR_LengthAniso2D.h"
#include "GR_misc.h"
#include "GR_Mesh2D.h"
#include "GR_Entity.h"
#include "GR_Vertex.h"
#include "GR_Face.h"
#include "GR_Cell.h"
#include "GR_GRPoint.h"
#include "GR_GRCurve.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "CubitMatrix.hpp"
#include "DLIList.hpp"
#include "RTree.hpp"
#include "GR_ADT.h"

using std::queue;
using std::set;
using std::pair;
using std::make_pair;

// requires simplicial mesh
void
LengthAniso2D::provideMesh(Mesh* pM2D, char* /*strLenFileName*/)
{
  assert(m_lengthType == eMesh);

  // Copy mesh verts and cells into new mesh object.
  m_lengthMesh = new Mesh2D(1);
  for (GR_index_t iV = 0; iV < pM2D->getNumVerts(); iV++) {
    Vert* pV = pM2D->getVert(iV);
    m_lengthMesh->createVert(pV->getCoords());
  }
  for (GR_index_t iC = 0; iC < pM2D->getNumCells(); iC++) {
    Cell* pC = pM2D->getCell(iC);
    int iV0 = pM2D->getVertIndex(pC->getVert(0));
    int iV1 = pM2D->getVertIndex(pC->getVert(1));
    int iV2 = pM2D->getVertIndex(pC->getVert(2));
    dynamic_cast<Mesh2D*>(m_lengthMesh)->createTriCell(
	m_lengthMesh->getVert(iV0), m_lengthMesh->getVert(iV1),
	m_lengthMesh->getVert(iV2));
  }
  m_lengthMesh->purgeAllEntities();

  // Initialize length metric for m_lengthMesh using distances to neighboring vertices.
  LengthAniso2D* pL2D = new LengthAniso2D(LengthAniso2D::eNeigh, 1, 1,
					  -LARGE_DBL);
  pL2D->initLength(m_lengthMesh);
  delete pL2D;

  /*
   // Overwrite length metric if initialized in the given mesh.
   for(int iV = 0; iV < pM2D->getNumVerts(); iV++){
   Vert* pVnew = m_lengthMesh->getVert(iV);
   double adMetric[5];
   for(int iMet = 0; iMet < 5; iMet++){
   adMetric[iMet] = pM2D->getVert(iV)->getMetric(iMet);
   }
   if(adMetric[0] < LARGE_DBL*0.9999 && adMetric[0] > 1e-10){
   pVnew->setMetric(adMetric);
   }
   }
   */

  // Create a tree of the cells from that mesh to make interpolation
  // efficient
  createADTFromAllCells(m_lengthMesh, m_BBoxTree);
  writeVTK(*(dynamic_cast<Mesh2D*>(m_lengthMesh)), "lengthMesh",
		" ",
	   Mesh::eMetric);
  writeVTK(*(dynamic_cast<Mesh2D*>(m_lengthMesh)), "lengthMesh2");
}

void
LengthAniso2D::evaluateLengthFromNeighbors(
    Vert* const vertex, double* adMetric,
    const std::set<Vert*>& neigh_verts) const
{
  double r_max = 0, r_min = LARGE_DBL, angle;
  bool qEllipseFound = false;

  if (neigh_verts.size() >= 4) {
    // Try to fit an ellipse to neighbor positions.
    CubitMatrix B(neigh_verts.size() + 1, 1);
    for (unsigned i = 0; i < neigh_verts.size(); i++) {
      B.set(i, 0, 0);
    }
    B.set(neigh_verts.size(), 0, -1);

    CubitMatrix A(neigh_verts.size() + 1, 4);
    int iRow = 0;
    for (Vert* pVneigh : neigh_verts) {
      double x = pVneigh->x() - vertex->x();
      double y = pVneigh->y() - vertex->y();

      A.set(iRow, 0, x * x);
      A.set(iRow, 1, x * y);
      A.set(iRow, 2, y * y);
      A.set(iRow, 3, 1);

      iRow++;
    }
    iRow = neigh_verts.size();
    A.set(iRow, 0, 0);
    A.set(iRow, 1, 0);
    A.set(iRow, 2, 0);
    A.set(iRow, 3, 1);

    CubitMatrix X = (A.transpose() * A).inverse() * A.transpose() * B;
    assert(X.num_rows() == 4 && X.num_cols() == 1);

    double a = X.get(0, 0);
    double b = X.get(1, 0);
    double c = X.get(2, 0);
    double f = X.get(3, 0);
    if (b * b - 4 * a * c < 0) {
      qEllipseFound = true;

      // http://math.stackexchange.com/questions/616645/determining-the-major-minor-axes-of-an-ellipse-from-general-form
      double q = 64 * f / (4 * a * c - b * b);
      double s = 0.25 * sqrt(fabs(q) * sqrt(b * b + (a - c) * (a - c)));
      r_max = 0.125
	  * sqrt(
	      2 * fabs(q) * sqrt(b * b + (a - c) * (a - c)) - 2 * q * (a + c));
      r_min = sqrt(r_max * r_max - s * s);

      double p1 = q * (a - c);
      double p2 = q * b;
      if (p1 == 0) {
	angle = M_PI * (p2 > 0 ? 0.25 : 0.75);
      }
      else if (p1 > 0) {
	angle = 0.5 * atan(b / (a - c)) + (p2 < 0 ? M_PI : 0);
      }
      else {
	angle = 0.5 * atan(b / (a - c)) + 0.5 * M_PI;
      }
    }
  }

  if (!qEllipseFound) {
    // Conic fit wasn't an ellipse, or less than four neighbors.
    // Take longest neighbor edge as semi-major axis. Use length of shortest neighbor edge as semi_minor.

    for (Vert* pVneigh : neigh_verts) {
      double x = pVneigh->x() - vertex->x();
      double y = pVneigh->y() - vertex->y();

      double neigh_length = sqrt(x * x + y * y);

      if (neigh_length > r_max) {
	r_max = neigh_length;
	angle = atan2(y, x);
      }
      if (neigh_length < r_min) {
	r_min = neigh_length;
      }
    }
  }

  if (!finite(r_max) || r_min == 0) {

    // Take longest neighbor edge as semi-major axis. Use length of shortest neighbor edge as semi_minor.
    r_max = 0;
    r_min = LARGE_DBL;

    for (Vert* pVneigh : neigh_verts) {
      double x = pVneigh->x() - vertex->x();
      double y = pVneigh->y() - vertex->y();
      double neigh_length = sqrt(x * x + y * y);

      if (neigh_length > r_max) {
	r_max = neigh_length;
	angle = atan2(y, x);
      }
      if (neigh_length < r_min) {
	r_min = neigh_length;
      }
    }
  }

  double c = cos(angle);
  double s = sin(angle);
  adMetric[0] = (c * c) / (r_max * r_max) + (s * s) / (r_min * r_min);
  adMetric[1] = (c * s) / (r_max * r_max) - (c * s) / (r_min * r_min);
  adMetric[2] = (s * s) / (r_max * r_max) + (c * c) / (r_min * r_min);

  assert(finite(adMetric[0]));
}

static void
diagonalizeMetric(double* adMetric, CubitMatrix &R, CubitMatrix &L)
{
  double a = adMetric[0];
  double b = adMetric[1];
  double c = adMetric[2];

  double l1 = 0.5 * (c + a + sqrt(4 * b * b - 2 * a * c + a * a + c * c));
  double l2 = 0.5 * (c + a - sqrt(4 * b * b - 2 * a * c + a * a + c * c));

  assert(finite(l1) && finite(l2));

  L.set(0, 0, l1);
  L.set(1, 1, l2);
  L.set(0, 1, 0);
  L.set(1, 0, 0);

  double bl1 = b / (l1 - a);
  double bl2 = b / (l2 - a);

  R.set(0, 0, bl1 / sqrt(1 + bl1 * bl1));
  R.set(1, 0, 1.0 / sqrt(1 + bl1 * bl1));
  R.set(0, 1, bl2 / sqrt(1 + bl2 * bl2));
  R.set(1, 1, 1.0 / sqrt(1 + bl2 * bl2));

  if (fabs(b) / fmax(fabs(a), fabs(c)) < 1e-7) {
    R.set(0, 0, 1);
    R.set(1, 0, 0);

    R.set(0, 1, 0);
    R.set(1, 1, 1);
  }
  else {
    assert(finite(bl1) && finite(bl2));
  }
}

static CubitMatrix
logDiagonal(CubitMatrix diagonal)
{
  CubitMatrix ret = diagonal;
  ret.set(0, 0, log(ret.get(0, 0)));
  ret.set(1, 1, log(ret.get(1, 1)));
  return ret;
}
static CubitMatrix
expDiagonal(CubitMatrix diagonal)
{
  CubitMatrix ret = diagonal;
  ret.set(0, 0, exp(ret.get(0, 0)));
  ret.set(1, 1, exp(ret.get(1, 1)));
  return ret;
}

/*
 * Currently copied from scat2d. Should refactor so that this code exists only in one place.
 */
void
LengthAniso2D::evaluateLengthFromMesh(Vert* const vertex,
				      double* adMetric) const
{
  assert(m_lengthType == eMesh);
  assert(m_lengthMesh);
  assert(m_BBoxTree);

  double adCoord[2] =
    { vertex->x(), vertex->y() };
  std::vector<GR_index_t> veciQueryResult;

  // Find cell (nearly?) containing the target point
  double dEps = 1.e-8;
  do {
    double adBB[] =
      { adCoord[0] - dEps, adCoord[0] + dEps, adCoord[1] - dEps, adCoord[1]
	  + dEps };
    veciQueryResult = m_BBoxTree->veciRangeQuery(adBB);
    dEps *= 10;
  }
  while (veciQueryResult.empty());

  //Pick the cell that contains the vertex point by evaluating barycentric coordinates.
  double dBestValue = 1.e50;
  int iBestCell = -1;
  double adBestBary[3];
  Cell *pC;
  int iQ;
  GR_index_t iCell;
  for (iQ = veciQueryResult.size() - 1; iQ >= 0 && dBestValue > 0.5; iQ--) {
    // Get barycentric coords.  If best_value <(=) 0.5 on completion
    // then the vert lies within (on an edge of) the cell in list
    iCell = veciQueryResult[iQ];
    pC = m_lengthMesh->getCell(iCell);
    double adBary[3];
    dynamic_cast<TriCell*>(pC)->calcBarycentricCoords(adCoord, adBary);

    double dThisValue = max(max(fabs(adBary[0] - 0.5), fabs(adBary[1] - 0.5)),
			    fabs(adBary[2] - 0.5));

    if (dThisValue < dBestValue) {
      iBestCell = iCell;
      dBestValue = dThisValue;
      adBestBary[0] = adBary[0];
      adBestBary[1] = adBary[1];
      adBestBary[2] = adBary[2];
    }
  } // Done checking all cells which have bboxes overlapping with the
    // vertex
  assert(iBestCell >= 0);
  pC = m_lengthMesh->getCell(iBestCell);

  dEps = 1e-8;
  if (dBestValue < (0.5 + dEps)) {
    // If a cell contains the vertex point, interpolate metric using Log-Euclidean framework.
    for (int iMet = 0; iMet < 5; iMet++) {
      adMetric[iMet] = 0;
    }

    for (int ii = 0; ii < 3; ii++) {
      Vert* pV = pC->getVert(ii);

//			for(int iMet = 0; iMet < 5; iMet++){
//				adMetric[iMet] += pV->getMetric(iMet) * adBestBary[ii];
//			}

      double adPVMet[3] =
	{ pV->getMetric(0), pV->getMetric(1), pV->getMetric(2) };
      CubitMatrix R(2, 2);
      CubitMatrix L(2, 2);
      diagonalizeMetric(adPVMet, R, L);

      CubitMatrix logM = R * logDiagonal(L) * R.transpose();

      adMetric[0] += logM.get(0, 0) * adBestBary[ii];
      adMetric[1] += logM.get(1, 0) * adBestBary[ii];
      adMetric[2] += logM.get(1, 1) * adBestBary[ii];
    }

    CubitMatrix R(2, 2);
    CubitMatrix L(2, 2);
    diagonalizeMetric(adMetric, R, L);

    CubitMatrix M = R * expDiagonal(L) * R.transpose();
    adMetric[0] = M.get(0, 0);
    adMetric[1] = M.get(1, 0);
    adMetric[2] = M.get(1, 1);

  }
  else {
    // Currently require m_lengthMesh to cover the domain.
    assert(0);
  }
}
