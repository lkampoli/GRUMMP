/*
 * VTK2D.cxx
 *
 *  Created on: Nov 6, 2017
 *      Author: cfog
 */

#include <cstdio>

#include "GR_Mesh2D.h"
#include "GR_Aniso.h"

#ifdef HAVE_MESQUITE
#include "CubitMatrix.hpp"
#endif

#define DISCARD_LINE do {			\
     lineNum++;                                 \
    int _res_ = fscanf(pFInFile, "%*[^\n]\n");	\
    _res_++;					\
  } while(0)

#define CHECK_READ(b, a) do {						\
    int _res_ = (a);							\
    if (_res_ != b) {							\
      logMessage(0, "File %s, line %d: \n", strFileName, lineNum);	\
      vFatalError("File format error", __func__);			\
    }									\
  } while(0)

#define CHECK_RANGE(val_, min_, max_) do {				\
    if (val_ < min_ || val_ > sGR_index_t(max_)) {				\
      logMessage(0, "File %s, line %d: \n", strFileName, lineNum);	\
      vFatalError("Index range error", __func__);			\
    }									\
  } while(0)

#define CHECK_POS_RANGE(val_, max_) do {					\
    if (val_ > max_) {							\
      logMessage(0, "File %s, line %d: \n", strFileName, lineNum);	\
      vFatalError("Index range error", __func__);			\
    }									\
  } while(0)

#define CHECK_VERTEX_RANGE(val_) CHECK_POS_RANGE(val_, nVerts)

void
readVTK(Mesh2D& M2D, const char strFileName[])
{
//  Mesh::eMeshType eMT = M2D.getType();
  FILE *pFInFile = fopen(strFileName, "r");
  if (NULL == pFInFile) {
    vFatalError("Couldn't open input file for reading", __func__);
  }
  logMessage(1,
	     "Reading VTK legacy input (triangles, quads, bdry edges only)\n");

  // Discard the four line header
  GR_index_t lineNum = 1;
  DISCARD_LINE
  ;
  DISCARD_LINE
  ;
  DISCARD_LINE
  ;
  DISCARD_LINE
  ;
  GR_index_t nVerts, nCells;
  CHECK_READ(1, (fscanf(pFInFile, " POINTS %u float\n", &nVerts)));
  lineNum++;
  if (nVerts == 0)
    vFatalError("Number of vertices is specified to be 0!", __func__);

  for (GR_index_t i = 0; i < nVerts; i++) {
    double x, y, z;
    CHECK_READ(3, (fscanf(pFInFile, "%lf%lf%lf\n", &x, &y, &z)));
    lineNum++;
    (void) M2D.createVert(x, y, z);
  }

  CHECK_READ(1, (fscanf(pFInFile, " CELLS %u %*d\n", &nCells)));
  lineNum++;
  if (nCells == 0)
    vFatalError("Number of cells is specified to be 0!", __func__);

  // These next three only really make sense if there is bdry data in the file
  int (*edgeVerts)[2] = new int[nCells][2];

  // bdryEdgeType is -1 for cells, 0 for domain bdry, 1 for internal bdry,
  // 2 for part bdry (owned), 3 for part bdry (ghost)
  int *bdryEdgeType = new int[nCells];

  // bdryData is the BC for domain and internal bdrys, and the opposite part
  // number for part bdrys
  int *bdryData = new int[nCells];

  // With this initialization, it's easy to tell later which data is valid for
  // edges and which isn't
  for (GR_index_t ii = 0; ii < nCells; ii++) {
    edgeVerts[ii][0] = edgeVerts[ii][1] = bdryEdgeType[ii] = bdryData[ii] = -1;
  }

  GR_index_t maxVert = 0;
  GR_index_t numBdryEdges = 0;
  for (GR_index_t c = 0; c < nCells; c++) {
    int nVertsThisCell;
    lineNum++;
    CHECK_READ(1, (fscanf(pFInFile, "%d", &nVertsThisCell)));
    switch (nVertsThisCell)
      {
      case 2:
	{
	  numBdryEdges++;
	  GR_index_t vertA, vertB;
	  CHECK_READ(2, (fscanf(pFInFile, "%u%u", &vertA, &vertB)));
	  CHECK_VERTEX_RANGE(vertA);
	  CHECK_VERTEX_RANGE(vertB);
	  maxVert = std::max(maxVert, std::max(vertA, vertB));
	  edgeVerts[c][0] = vertA;
	  edgeVerts[c][1] = vertB;
	  break;
	}
      case 3:
	{
	  GR_index_t vertA, vertB, vertC;
	  CHECK_READ(3, (fscanf(pFInFile, "%u%u%u", &vertA, &vertB, &vertC)));
	  CHECK_VERTEX_RANGE(vertA);
	  CHECK_VERTEX_RANGE(vertB);
	  CHECK_VERTEX_RANGE(vertC);
	  maxVert = std::max(std::max(maxVert, vertA), std::max(vertB, vertC));
	  (void) M2D.createTriCell(M2D.getVert(vertA), M2D.getVert(vertB),
				   M2D.getVert(vertC), iDefaultRegion, true);
	  break;
	}
      case 4:
	{
	  M2D.allowNonSimplicial();
	  GR_index_t vertA, vertB, vertC, vertD;
	  CHECK_READ(
	      4,
	      (fscanf(pFInFile, "%u%u%u%u", &vertA, &vertB, &vertC, &vertD)));
	  CHECK_VERTEX_RANGE(vertA);
	  CHECK_VERTEX_RANGE(vertB);
	  CHECK_VERTEX_RANGE(vertC);
	  CHECK_VERTEX_RANGE(vertD);
	  maxVert = std::max(
	      maxVert,
	      std::max(std::max(vertA, vertB), std::max(vertC, vertD)));
	  (void) M2D.createQuadCell(M2D.getVert(vertA), M2D.getVert(vertB),
				    M2D.getVert(vertC), M2D.getVert(vertD),
				    iDefaultRegion);
	  break;
	}
      default:
	vFatalError("Invalid number of verts for cell.\n", __func__);
	break;
      }
  }

  if (numBdryEdges != 0) {
    // Read the cell type data primarily to discard it.
    GR_index_t check;
    CHECK_READ(1, (fscanf(pFInFile, " CELL_TYPES %u\n", &check)));
    lineNum++;
    assert(check == nCells);
    for (GR_index_t ii = 0; ii < nCells; ii++) {
      unsigned int cellType;
      CHECK_READ(1, fscanf(pFInFile, "%u ", &cellType));
      lineNum++;
    }

    // Now read the cell attribute data
    CHECK_READ(1, (fscanf(pFInFile, " CELL_DATA %u ", &check)));
    lineNum++;
    assert(check == nCells);

    // Region data gets discarded presently.
    // TODO: Actually use this properly.  This will likely require reading
    // all the connectivity and region data up front, and -then- starting to
    // build cells.  This is similar to what's done in 3D.

    DISCARD_LINE
    ;
    DISCARD_LINE
    ;
    for (GR_index_t ii = 0; ii < nCells; ii++) {
      unsigned int region;
      CHECK_READ(1, fscanf(pFInFile, "%d ", &region));
      lineNum++;
    }

    // Bdry face type
    DISCARD_LINE
    ;
    DISCARD_LINE
    ;
    for (GR_index_t ii = 0; ii < nCells; ii++) {
      CHECK_READ(1, fscanf(pFInFile, "%d ", &bdryEdgeType[ii]));
      lineNum++;
      assert(bdryEdgeType[ii] >= -1 && bdryEdgeType[ii] <= 3);
    }

    DISCARD_LINE
    ;
    DISCARD_LINE
    ;
    for (GR_index_t ii = 0; ii < nCells; ii++) {
      CHECK_READ(1, fscanf(pFInFile, "%d ", &bdryData[ii]));
      lineNum++;
      assert(bdryData[ii] >= -1);
    }

    // Now work through the edge data.
    for (GR_index_t ii = 0; ii < nCells; ii++) {
      if (edgeVerts[ii][0] == -1) {
	assert(edgeVerts[ii][1] == -1);
	assert(bdryEdgeType[ii] == -1);
	assert(bdryData[ii] == -1);
      }
      else {
	assert(edgeVerts[ii][1] >= 0);
	assert(bdryEdgeType[ii] >= 0);
	assert(bdryData[ii] >= 0);

	Vert *pVA = M2D.getVert(edgeVerts[ii][0]);
	Vert *pVB = M2D.getVert(edgeVerts[ii][1]);
	Face *pF = findCommonFace(pVA, pVB);

	switch (bdryEdgeType[ii])
	  {
	  case 0: // Domain bdry edge
	    assert(
		pF->getLeftCell() == pCInvalidCell || pF->getRightCell() == pCInvalidCell);
	    M2D.createBFace(pF, bdryData[ii]);
	    break;
	  case 1: // Internal bdry edge
	    assert(
		pF->getLeftCell() != pCInvalidCell || pF->getRightCell() != pCInvalidCell);
	    assert(0); // We're not reading regions yet, so this is
	    // going to fail anyway.
	    break;
#ifdef HAVE_MPI
	  case 2: // Part bdry edge, owned by this part
	    {
	      assert(
		  pF->getLeftCell() == pCInvalidCell || pF->getRightCell() == pCInvalidCell);
	      pF->setOwnership(true); // Should be a no-op.
	      pF->setEntStatus(Entity::ePartBdry);
	      M2D.addGhostEntry(pF, M2D.getPartID(), 1);
	      bool success = false;
	      GRUMMP::EntityGID egid(bdryData[ii]);
	      success = M2D.addCopyToGhostEntry(pF, egid);
	      assert(success);
	      GRUMMP::Ghost g = M2D.getGhostEntry(pF, success);
	      assert(success);
	      g.setNumPartBdryCopies(1);
	      assert(g.getNumCopies() == 1 && g.getMaxNumCopies() == 1);
	      break;
	    }
	  case 3: // Part bdry edge, owned by the other part
	    {
	      assert(
		  pF->getLeftCell() == pCInvalidCell || pF->getRightCell() == pCInvalidCell);
	      pF->setOwnership(false);
	      GRUMMP::PartGID pgid(bdryData[ii]);
	      pF->setEntStatus(Entity::ePartBdry);
	      M2D.addGhostEntry(pF, pgid, 1);
	      bool success = false;
	      GRUMMP::Ghost g = M2D.getGhostEntry(pF, success);
	      assert(success);
	      g.setNumPartBdryCopies(1);
	      assert(g.getNumCopies() == 0 && g.getMaxNumCopies() == 1);
	      break;
	    }
#endif
	  default:
	    // Should never get here.
	    assert(0);
	    break;
	  } // switch on bdry edge type
      } // this is actually a bdry edge
    } // loop over cells

#ifdef HAVE_MPI
    // Finally, read info about vertex ownership.
    CHECK_READ(1, (fscanf(pFInFile, " POINT_DATA %u ", &check)));
    lineNum++;
    assert(check == nVerts);
    DISCARD_LINE
    ;
    DISCARD_LINE
    ;
    for (GR_index_t ii = 0; ii < nVerts; ii++) {
      int ownerInfo;
      // Values:
      //   <0: I own it; value indicates total number of copies, including mine
      //   >=0:  Part ID of owner.
      CHECK_READ(1, fscanf(pFInFile, "%d ", &ownerInfo));
      lineNum++;
      if (ownerInfo == -1)
	continue;
      Vert *pV = M2D.getVert(ii);
      pV->setEntStatus(Entity::ePartBdry);
      if (ownerInfo < 0) {
	// I own it, and I now know how many ghost copies there are
	M2D.addGhostEntry(pV, M2D.getPartID(), (-ownerInfo) - 1);
	bool success = false;
	GRUMMP::Ghost g = M2D.getGhostEntry(pV, success);
	assert(success);
	g.setNumPartBdryCopies((-ownerInfo) - 1);
      }
      else {
	// Someone else owns it, but I have a copy.
	GRUMMP::PartGID pgid(ownerInfo);
	M2D.addGhostEntry(pV, pgid, 1);
	pV->setOwnership(false);
	bool success = false;
	GRUMMP::Ghost g = M2D.getGhostEntry(pV, success);
	assert(success);
	g.setNumPartBdryCopies(1);
      }
    }
#endif
  } // bdry edges exist in file

  if (nVerts != maxVert + 1)
    vFatalError("Vertex numbering in mesh doesn't match number of verts given.",
		__func__);

  // Now create some bdry faces, in case VTK doesn't provide them all.

  for (GR_index_t f = 0; f < M2D.getNumFaces(); f++) {
    Face *pF = M2D.getFace(f);
    if (pF->isDeleted() || pF->getEntStatus() == Entity::ePartBdry)
      continue;
    if (pF->getLeftCell() == pCInvalidCell
	|| pF->getRightCell() == pCInvalidCell) {
      (void) M2D.createBFace(pF, iDefaultBC);
    }
  }
  M2D.identifyVertexTypesGeometrically();
}

//-----------------------------------------------------------------------------
void
writeVTK(Mesh2D& OutMesh, const char strBaseFileName[],
	 const char strFileNameExt[], const enum Mesh::ExtraData ED)
//-----------------------------------------------------------------------------
{
  FILE *pFOutFile = NULL;
  char strFileName[1024];
  OutMesh.purgeAllEntities();
  if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
      || strlen(strBaseFileName) < 4) {
    // File name doesn't end in .vtk
    sprintf(strFileName, "%s%s.vtk", strBaseFileName, strFileNameExt);
  }
  else {
    // File name ends in .vtk
    sprintf(strFileName, "%s%s", strBaseFileName, strFileNameExt);
  }
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open Vtk output file for writing", "2D mesh output");

  GR_index_t nno = OutMesh.iNumVerts();   // number of nodes
  //int nel = OutMesh.iNumCells();
  GR_index_t ntri = OutMesh.getNumTriCells();	// number of Tri cells
  GR_index_t nquad = OutMesh.getNumQuadCells(); // numbers of quad cells
  GR_index_t nel = ntri + nquad;

  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0");
  fprintf(pFOutFile, "\nGRUMMP topologically 2D mesh");
  fprintf(pFOutFile, "\nASCII\n");
  fprintf(pFOutFile, "\nDATASET UNSTRUCTURED_GRID");
  fprintf(pFOutFile, "\nPOINTS %d float", nno);

  //-------------------------------------
  // write 2d vertex data
  //-------------------------------------
  Vert* pV = NULL;
  for (GR_index_t iV = 0; iV < nno; ++iV) {
    /*verts: coords*/
    pV = OutMesh.getVert(iV);
    fprintf(pFOutFile, "\n%16.8f %16.8f %16.8f", pV->x(), pV->y(), pV->z());
  }

  //-------------------------------------
  // write cell data as:
  //
  //    CELLS  n  size
  //    3          v1  v2  v3
  //    3          v1  v2  v3
  //    3          v1  v2  v3
  //
  // (#verts)     (vertex list)
  //
  //-------------------------------------

  //                      n size         no. ints to build connectivity
  fprintf(pFOutFile, "\n\nCELLS %d %d", nel, (nel + ntri * 3 + nquad * 4));

  int h, globalNodeNum;
  Cell* pC = NULL;
  for (GR_index_t e = 0; e < nel; ++e) {
    pC = OutMesh.getCell(e);
    h = pC->getNumVerts();
    fprintf(pFOutFile, "\n%d  ", h);
    for (int idx = 0; idx < h; ++idx) {
      globalNodeNum = OutMesh.iVertIndex(pC->getVert(idx));
      fprintf(pFOutFile, " %10d", globalNodeNum);
    }
  }
  fprintf(pFOutFile, "\n");

  //-------------------------------------
  // write cell type (VTK_TRIANGLE == 5, VTK_QUAD == 9)
  //-------------------------------------
  fprintf(pFOutFile, "\nCELL_TYPES %d\n", nel);
  for (GR_index_t iC = 0; iC < nel; ++iC) {
    int iNV = OutMesh.getCell(iC)->getNumVerts();
    fprintf(pFOutFile, "%d\n", iNV == 3 ? 5 : 9);
  }
  fprintf(pFOutFile, "\n");

  if (ED == Mesh::eStandard) {
    fprintf(pFOutFile, "POINT_DATA %u\n", nno);
    fprintf(pFOutFile, "SCALARS LengthScale float \n");
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    for (GR_index_t iV = 0; iV < nno; iV++) {
      fprintf(pFOutFile, "%f\n", OutMesh.getVert(iV)->getLengthScale());
    }
  }
  else {
#ifdef HAVE_MESQUITE
    if (ED == Mesh::eMetric) {
      // Now write metric data as tensors
      //----------------------------------------
      // Compute the following quantities implied by the metric at each vertex:
      //	- Local aspect ratio.
      //	- Angle of semi-major axis.
      //	- Scale of semi-major axis.
      //----------------------------------------

      double adAspect[nno];
      double adAngle[nno];
      double adScale[nno];

      CubitMatrix R(2, 2);
      CubitMatrix L(2, 2);
      for (GR_index_t iV = 0; iV < nno; iV++) {
	pV = OutMesh.getVert(iV);

	double a = pV->getMetric(0);
	double b = pV->getMetric(1);
	double c = pV->getMetric(2);

	double l1 = 0.5
	    * (c + a + sqrt(fmax(4 * b * b - 2 * a * c + a * a + c * c, 0.)));
	double l2 = 0.5
	    * (c + a - sqrt(fmax(4 * b * b - 2 * a * c + a * a + c * c, 0.)));

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

	if (fabs(b) < 1e-10) {
	  R.set(0, 0, 1);
	  R.set(1, 0, 0);

	  R.set(0, 1, 0);
	  R.set(1, 1, 1);
	}
	else {
	  assert(finite(bl1) && finite(bl2));
	}

	if (l1 > l2) {
	  adAspect[iV] = sqrt(l1 / l2);
	  adAngle[iV] = atan(R.get(1, 1) / R.get(0, 1));
	  adScale[iV] = 1.0 / sqrt(l2);
	}
	else {
	  adAspect[iV] = sqrt(l2 / l1);
	  adAngle[iV] = atan(R.get(1, 0) / R.get(0, 0));
	  adScale[iV] = 1.0 / sqrt(l1);
	}
      }

      fprintf(pFOutFile, "POINT_DATA %u\n", nno);

      fprintf(pFOutFile, "SCALARS Aspect float 1\n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iV = 0; iV < nno; iV++) {
	fprintf(pFOutFile, "%.15g\n", adAspect[iV]);
      }

      fprintf(pFOutFile, "SCALARS Angle float 1\n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iV = 0; iV < nno; iV++) {
	fprintf(pFOutFile, "%.15g\n", adAngle[iV]);
      }

      fprintf(pFOutFile, "SCALARS Scale float 1\n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iV = 0; iV < nno; iV++) {
	fprintf(pFOutFile, "%.15g\n", adScale[iV]);
      }

      fprintf(pFOutFile, "TENSORS Metric float\n");
      for (GR_index_t iV = 0; iV < nno; iV++) {
	pV = OutMesh.getVert(iV);
	fprintf(pFOutFile, "%.15g %.15g %.15g\n", pV->getMetric(0),
		pV->getMetric(1), 0.);
	fprintf(pFOutFile, "%.15g %.15g %.15g\n", pV->getMetric(1),
		pV->getMetric(2), 0.);
	fprintf(pFOutFile, "%.15g %.15g %.15g\n", 0., 0., 0.);
      }
    }
#endif
    if (ED == Mesh::eTMOPQual) {
      // Without Mesquite, this request silently does nothing
#ifdef HAVE_MESQUITE
      // Now compute and write TMOP quality information
      //-------------------------------------------------------
      // Computes TMOP quality for mesh (Size, Shape, ShapeOrient, ShapeSize, ShapeSizeOrient)
      //-------------------------------------------------------

      fprintf(pFOutFile, "CELL_DATA %d\n", nel);
      fprintf(pFOutFile, "SCALARS ShapeOrient float \n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
	pC = OutMesh.getCell(iter);

	double dT[9];
	double dq;
	GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeOrient;
	GRUMMP::TMOPQual TQ;

	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pC->getVert(0),
			       pC->getVert(1), pC->getVert(2));
	dq = TQ.dEval(QM, CellSkel::eTriCell, dT, pC->getVert(0)->getCoords(),
		      pC->getVert(1)->getCoords(), pC->getVert(2)->getCoords());
	fprintf(pFOutFile, "%f\n", dq);
      }
      fprintf(pFOutFile, "\n");

      fprintf(pFOutFile, "SCALARS ShapeSizeOrient float \n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
	pC = OutMesh.getCell(iter);

	double dT[9];
	double dq;
	GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeSizeOrient;
	GRUMMP::TMOPQual TQ;

	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pC->getVert(0),
			       pC->getVert(1), pC->getVert(2));
	dq = TQ.dEval(QM, CellSkel::eTriCell, dT, pC->getVert(0)->getCoords(),
		      pC->getVert(1)->getCoords(), pC->getVert(2)->getCoords());
	fprintf(pFOutFile, "%f\n", dq);
      }
      fprintf(pFOutFile, "\n");

      fprintf(pFOutFile, "SCALARS Shape float \n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
	pC = OutMesh.getCell(iter);

	double dT[9];
	double dq;
	GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShape;
	GRUMMP::TMOPQual TQ;

	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pC->getVert(0),
			       pC->getVert(1), pC->getVert(2));
	dq = TQ.dEval(QM, CellSkel::eTriCell, dT, pC->getVert(0)->getCoords(),
		      pC->getVert(1)->getCoords(), pC->getVert(2)->getCoords());
	fprintf(pFOutFile, "%f\n", dq);
      }
      fprintf(pFOutFile, "\n");

      fprintf(pFOutFile, "SCALARS Size float \n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
	pC = OutMesh.getCell(iter);

	double dT[9];
	double dq;
	GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eSize;
	GRUMMP::TMOPQual TQ;

	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pC->getVert(0),
			       pC->getVert(1), pC->getVert(2));
	dq = TQ.dEval(QM, CellSkel::eTriCell, dT, pC->getVert(0)->getCoords(),
		      pC->getVert(1)->getCoords(), pC->getVert(2)->getCoords());
	fprintf(pFOutFile, "%f\n", dq);
      }
      fprintf(pFOutFile, "\n");

      fprintf(pFOutFile, "SCALARS ShapeSize float \n");
      fprintf(pFOutFile, "LOOKUP_TABLE default\n");
      for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
	pC = OutMesh.getCell(iter);

	double dT[9];
	double dq;
	GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeSize;
	GRUMMP::TMOPQual TQ;

	TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pC->getVert(0),
			       pC->getVert(1), pC->getVert(2));
	dq = TQ.dEval(QM, CellSkel::eTriCell, dT, pC->getVert(0)->getCoords(),
		      pC->getVert(1)->getCoords(), pC->getVert(2)->getCoords());
	fprintf(pFOutFile, "%f\n", dq);
      }
      fprintf(pFOutFile, "\n");
#endif
    }
  }
  fclose(pFOutFile);

  logMessage(MSG_MANAGER, "Finished writing VTK dataFile, %s\n", strFileName);
  logMessage(MSG_SERVICE, "The mesh has %7u vertices\n", OutMesh.iNumVerts());
  logMessage(MSG_SERVICE, "             %7u faces\n", OutMesh.getNumFaces());
  logMessage(MSG_SERVICE, "             %7u bdry faces\n",
	     OutMesh.getNumBdryFaces());
  logMessage(MSG_SERVICE, "        and: %7u cells\n", OutMesh.getNumCells());
  logMessage(MSG_SERVICE, "   of which: %7u are triangles\n",
	     OutMesh.getNumTriCells());
  logMessage(MSG_SERVICE, "             %7u are quads\n",
	     OutMesh.getNumQuadCells());
}

