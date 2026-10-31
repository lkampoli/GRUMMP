#include <stdio.h>
#include "GR_Mesh.h"
#include "GR_Mesh2D.h"
#include "GR_misc.h"
#include "GR_TMOPQual.h"
#include "GR_VolMesh.h"

//-----------------------------------------------------------------------------
void
writeVTKLegacy(VolMesh& OutMesh, const char strBaseFileName[],
	       const char strFileNameExt[], const int *aCellMatch)
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
    sprintf(strFileName, "%s", strBaseFileName);
  }
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d Vtk mesh output");

  GR_index_t nno = OutMesh.getNumVerts();   // number of nodes

  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "GRUMMP Tetra example\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", nno);

  //-------------------------------------
  // write 3d vertex data
  //-------------------------------------
  for (GR_index_t i = 0; i < nno; ++i) {
    Vert *pV = OutMesh.getVert(i);
    fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
  }

  // How much topology data is there?

  // Tri and Quad data not currently output in 3D. TODO: But could
  // use this for bdry faces.
  // const int numTris = OutMesh.getNumTriFaces();
  // const int numQuads = OutMesh.getNumQuadFaces();
  const int numTets = OutMesh.getNumTetCells();
  const int numPyrs = OutMesh.getNumPyrCells();
  const int numPrisms = OutMesh.getNumPrismCells();
  const int numHexes = OutMesh.getNumHexCells();
  const int numIBs = OutMesh.getNumIntBdryFaces();
  const int numEnts = /*numTris + numQuads +*/numTets + numPyrs + numPrisms
      + numHexes + numIBs;
  const int dataSize = /*4*numTris*/+5 * (/*numQuads +*/numTets) + 6 * numPyrs
      + 7 * numPrisms + 9 * numHexes + 4 * numIBs;

  fprintf(pFOutFile, "CELLS %d %d\n", numEnts, dataSize);

  int h, globalNodeNum;
  // First all the faces.
  // Not writing faces at this point
  //  for (int iF=0; iF < numTris + numQuads; iF++) {
  //    Face *pF = OutMesh.getFace(iF);
  //    h  = pF->getNumVerts();
  //    fprintf(pFOutFile, "%d  ", h);
  //    for (int idx=0; idx<h; ++idx)
  //    {
  //      globalNodeNum = OutMesh.getVertIndex(pF->getVert(idx));
  //      fprintf(pFOutFile, " %10d", globalNodeNum);
  //    }
  //    fprintf(pFOutFile, "\n");
  //  }

  // Now all the cells.
  for (int iC = 0; iC < numTets + numPyrs + numPrisms + numHexes; iC++) {
    Cell *pC = OutMesh.getCell(iC);
    h = pC->getNumVerts();
    fprintf(pFOutFile, "%d  ", h);
    for (int idx = 0; idx < h; ++idx) {
      globalNodeNum = OutMesh.getVertIndex(pC->getVert(idx));
      fprintf(pFOutFile, " %10d", globalNodeNum);
    }
    fprintf(pFOutFile, "\n");
  }
  for (int iC = 0; iC < numIBs; iC++) {
    Cell *pC = OutMesh.getIntBFace(iC);
    h = pC->getNumVerts();
    fprintf(pFOutFile, "%d  ", h);
    for (int idx = 0; idx < h; ++idx) {
      globalNodeNum = OutMesh.getVertIndex(pC->getVert(idx));
      fprintf(pFOutFile, " %10d", globalNodeNum);
    }
    fprintf(pFOutFile, "\n");
  }

  //-------------------------------------
  // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  // VTK_PYRAMID = 14)
  //-------------------------------------
  fprintf(pFOutFile, "CELL_TYPES %d\n", numEnts);
  // Not writing faces at this point
  //  for (int ct = 0; ct < numTris; ++ct)
  //    fprintf(pFOutFile, "5\n");
  //  for (int ct = 0; ct < numQuads; ++ct)
  //    fprintf(pFOutFile, "9\n");
  for (int ct = 0; ct < numTets; ++ct)
    fprintf(pFOutFile, "10\n");
  for (int ct = 0; ct < numPyrs; ++ct)
    fprintf(pFOutFile, "14\n");
  for (int ct = 0; ct < numPrisms; ++ct)
    fprintf(pFOutFile, "13\n");
  for (int ct = 0; ct < numHexes; ++ct)
    fprintf(pFOutFile, "12\n");
  for (int ct = 0; ct < numIBs; ++ct)
    fprintf(pFOutFile, "5\n");
  //----------------------------------------
  // Write vertex length scale
  //----------------------------------------
  if (aCellMatch == NULL) {
    //----------------------------------------
    // Write vertex length scale
    //----------------------------------------
    fprintf(pFOutFile, "POINT_DATA %u\n", nno);
    fprintf(pFOutFile, "SCALARS Vert-Type float 1\n");
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    fflush(pFOutFile);
    for (GR_index_t iV = 0; iV < nno; iV++) {
      fprintf(pFOutFile, "%d\n", OutMesh.getVert(iV)->getVertType());
    }
    fprintf(pFOutFile, "SCALARS LengthScale float 1\n");
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    fflush(pFOutFile);
    for (GR_index_t iV = 0; iV < nno; iV++) {
      fprintf(pFOutFile, "%12G\n",
	      std::min(double(FLT_MAX), OutMesh.getVert(iV)->getLengthScale()));
    }

  }
  else {

    //-------------------------------------
    // Highlight cells that do not match
    //-------------------------------------

    fprintf(pFOutFile, "CELL_DATA %d\n", numEnts);
    fprintf(pFOutFile, "SCALARS diffCells int %d\n", 1);
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    for (int iter = 0; iter < numTets; ++iter) {
      fprintf(pFOutFile, "%d\n", aCellMatch[iter]);
    }
    for (int ct = 0; ct < numIBs; ++ct)
      fprintf(pFOutFile, "3\n");
  }

  fprintf(pFOutFile, "CELL_DATA %u\n", numEnts);
  fprintf(pFOutFile, "SCALARS CellType float 1\n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (int ct = 0; ct < numTets; ++ct)
    fprintf(pFOutFile, "%d\n", OutMesh.getCell(ct)->getRegion());
  for (int ct = 0; ct < numIBs; ++ct)
    fprintf(pFOutFile, "-1\n");

  fclose(pFOutFile);
  logMessage(2, "Finished writing VTK dataFfile, %s\n", strFileName);
}

#ifdef HAVE_MESQUITE
//-----------------------------------------------------------------------------
void
writeVTKLegacy3DTMOPQual(VolMesh& OutMesh, const char strBaseFileName[],
			 const char strFileNameExt[], const int *aCellMatch)
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
    sprintf(strFileName, "%s", strBaseFileName);
  }
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d Vtk mesh output");

  GR_index_t nno = OutMesh.getNumVerts();   // number of nodes

  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "GRUMMP Tetra example\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", nno);

  //-------------------------------------
  // write 3d vertex data
  //-------------------------------------
  for (GR_index_t i = 0; i < nno; ++i) {
    Vert *pV = OutMesh.getVert(i);
    fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
  }

  // How much topology data is there?

  // Tri and Quad data not currently output in 3D. TODO: But could
  // use this for bdry faces.
  // const int numTris = OutMesh.getNumTriFaces();
  // const int numQuads = OutMesh.getNumQuadFaces();
  const int numTets = OutMesh.getNumTetCells();
  const int numPyrs = OutMesh.getNumPyrCells();
  const int numPrisms = OutMesh.getNumPrismCells();
  const int numHexes = OutMesh.getNumHexCells();
  const int numIBs = OutMesh.getNumIntBdryFaces();
  const int numEnts = /*numTris + numQuads +*/numTets + numPyrs + numPrisms
      + numHexes + numIBs;
  const int dataSize = /*4*numTris*/+5 * (/*numQuads +*/numTets) + 6 * numPyrs
      + 7 * numPrisms + 9 * numHexes + 4 * numIBs;

  fprintf(pFOutFile, "CELLS %d %d\n", numEnts, dataSize);

  int h, globalNodeNum;
  // First all the faces.
  // Not writing faces at this point
  //  for (int iF=0; iF < numTris + numQuads; iF++) {
  //    Face *pF = OutMesh.getFace(iF);
  //    h  = pF->getNumVerts();
  //    fprintf(pFOutFile, "%d  ", h);
  //    for (int idx=0; idx<h; ++idx)
  //    {
  //      globalNodeNum = OutMesh.getVertIndex(pF->getVert(idx));
  //      fprintf(pFOutFile, " %10d", globalNodeNum);
  //    }
  //    fprintf(pFOutFile, "\n");
  //  }

  // Now all the cells.
  for (int iC = 0; iC < numTets + numPyrs + numPrisms + numHexes; iC++) {
    Cell *pC = OutMesh.getCell(iC);
    h = pC->getNumVerts();
    fprintf(pFOutFile, "%d  ", h);
    for (int idx = 0; idx < h; ++idx) {
      globalNodeNum = OutMesh.getVertIndex(pC->getVert(idx));
      fprintf(pFOutFile, " %10d", globalNodeNum);
    }
    fprintf(pFOutFile, "\n");
  }
  for (int iC = 0; iC < numIBs; iC++) {
    Cell *pC = OutMesh.getIntBFace(iC);
    h = pC->getNumVerts();
    fprintf(pFOutFile, "%d  ", h);
    for (int idx = 0; idx < h; ++idx) {
      globalNodeNum = OutMesh.getVertIndex(pC->getVert(idx));
      fprintf(pFOutFile, " %10d", globalNodeNum);
    }
    fprintf(pFOutFile, "\n");
  }

  //-------------------------------------
  // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  // VTK_PYRAMID = 14)
  //-------------------------------------
  fprintf(pFOutFile, "CELL_TYPES %d\n", numEnts);
  // Not writing faces at this point
  //  for (int ct = 0; ct < numTris; ++ct)
  //    fprintf(pFOutFile, "5\n");
  //  for (int ct = 0; ct < numQuads; ++ct)
  //    fprintf(pFOutFile, "9\n");
  for (int ct = 0; ct < numTets; ++ct)
    fprintf(pFOutFile, "10\n");
  for (int ct = 0; ct < numPyrs; ++ct)
    fprintf(pFOutFile, "14\n");
  for (int ct = 0; ct < numPrisms; ++ct)
    fprintf(pFOutFile, "13\n");
  for (int ct = 0; ct < numHexes; ++ct)
    fprintf(pFOutFile, "12\n");
  for (int ct = 0; ct < numIBs; ++ct)
    fprintf(pFOutFile, "5\n");
  //----------------------------------------
  // Write vertex length scale
  //----------------------------------------
  if (aCellMatch == NULL) {
    //----------------------------------------
    // Write vertex length scale
    //----------------------------------------
    fprintf(pFOutFile, "POINT_DATA %u\n", nno);
    fprintf(pFOutFile, "SCALARS Vert-Type float 1\n");
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    for (GR_index_t iV = 0; iV < nno; iV++) {
      fprintf(pFOutFile, "%d\n", OutMesh.getVert(iV)->getVertType());
    }
  }
  else {

    //-------------------------------------
    // Highlight cells that do not match
    //-------------------------------------

    fprintf(pFOutFile, "CELL_DATA %d\n", numEnts);
    fprintf(pFOutFile, "SCALARS diffCells int %d\n", 1);
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    for (int iter = 0; iter < numTets; ++iter) {
      fprintf(pFOutFile, "%d\n", aCellMatch[iter]);
    }
    for (int ct = 0; ct < numIBs; ++ct)
      fprintf(pFOutFile, "3\n");
  }

  fprintf(pFOutFile, "CELL_DATA %u\n", numEnts);
  fprintf(pFOutFile, "SCALARS CellType float 1\n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (int ct = 0; ct < numTets; ++ct)
    fprintf(pFOutFile, "%d\n", OutMesh.getCell(ct)->getRegion());
  for (int ct = 0; ct < numIBs; ++ct)
    fprintf(pFOutFile, "-1\n");

  fprintf(pFOutFile, "SCALARS Shape float \n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
    Cell *pC;
    pC = OutMesh.getCell(iter);
    Vert *pVA = pC->getVert(0);
    Vert *pVB = pC->getVert(1);
    Vert *pVC = pC->getVert(2);
    Vert *pVD = pC->getVert(3);

    double dT[12];
    double dq;
    GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShape;
    GRUMMP::TMOPQual TQ;

    if (checkOrient3D(pVA, pVB, pVC, pVD) == -1) {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVC, pVB, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVC->getCoords(),
		    pVB->getCoords(), pVD->getCoords());
    }
    else {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVB, pVC, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVB->getCoords(),
		    pVC->getCoords(), pVD->getCoords());
    }

    fprintf(pFOutFile, "%f\n", dq);
  }
  fprintf(pFOutFile, "\n");

  fprintf(pFOutFile, "SCALARS Size float \n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
    Cell *pC;
    pC = OutMesh.getCell(iter);
    Vert *pVA = pC->getVert(0);
    Vert *pVB = pC->getVert(1);
    Vert *pVC = pC->getVert(2);
    Vert *pVD = pC->getVert(3);

    double dT[12];
    double dq;
    GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eSize;
    GRUMMP::TMOPQual TQ;

    if (checkOrient3D(pVA, pVB, pVC, pVD) == -1) {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVC, pVB, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVC->getCoords(),
		    pVB->getCoords(), pVD->getCoords());
    }
    else {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVB, pVC, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVB->getCoords(),
		    pVC->getCoords(), pVD->getCoords());
    }

    fprintf(pFOutFile, "%f\n", dq);
  }
  fprintf(pFOutFile, "\n");

  fprintf(pFOutFile, "SCALARS ShapeSize float \n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
    Cell *pC;
    pC = OutMesh.getCell(iter);
    Vert *pVA = pC->getVert(0);
    Vert *pVB = pC->getVert(1);
    Vert *pVC = pC->getVert(2);
    Vert *pVD = pC->getVert(3);

    double dT[12];
    double dq;
    GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeSize;
    GRUMMP::TMOPQual TQ;

    if (checkOrient3D(pVA, pVB, pVC, pVD) == -1) {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVC, pVB, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVC->getCoords(),
		    pVB->getCoords(), pVD->getCoords());
    }
    else {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVB, pVC, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVB->getCoords(),
		    pVC->getCoords(), pVD->getCoords());
    }

    fprintf(pFOutFile, "%f\n", dq);
  }
  fprintf(pFOutFile, "\n");

  fprintf(pFOutFile, "SCALARS ShapeOrient float \n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
    Cell *pC;
    pC = OutMesh.getCell(iter);
    Vert *pVA = pC->getVert(0);
    Vert *pVB = pC->getVert(1);
    Vert *pVC = pC->getVert(2);
    Vert *pVD = pC->getVert(3);

    double dT[12];
    double dq;
    GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeOrient;
    GRUMMP::TMOPQual TQ;

    if (checkOrient3D(pVA, pVB, pVC, pVD) == -1) {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVC, pVB, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVC->getCoords(),
		    pVB->getCoords(), pVD->getCoords());
    }
    else {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVB, pVC, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVB->getCoords(),
		    pVC->getCoords(), pVD->getCoords());
    }

    fprintf(pFOutFile, "%f\n", dq);
  }
  fprintf(pFOutFile, "\n");

  fprintf(pFOutFile, "SCALARS ShapeSizeOrient float \n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iter = 0; iter < OutMesh.getNumCells(); ++iter) {
    Cell *pC;
    pC = OutMesh.getCell(iter);
    Vert *pVA = pC->getVert(0);
    Vert *pVB = pC->getVert(1);
    Vert *pVC = pC->getVert(2);
    Vert *pVD = pC->getVert(3);

    double dT[12];
    double dq;
    GRUMMP::TMOPQual::eQMetric QM = GRUMMP::TMOPQual::eShapeSizeOrient;
    GRUMMP::TMOPQual TQ;

    if (checkOrient3D(pVA, pVB, pVC, pVD) == -1) {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVC, pVB, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVC->getCoords(),
		    pVB->getCoords(), pVD->getCoords());
    }
    else {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVB, pVC, pVD);
      dq = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(), pVB->getCoords(),
		    pVC->getCoords(), pVD->getCoords());
    }

    fprintf(pFOutFile, "%f\n", dq);
  }
  fprintf(pFOutFile, "\n");

  fprintf(pFOutFile, "POINT_DATA %d\n", nno);
  fprintf(pFOutFile, "SCALARS VertIndex float \n");
  fprintf(pFOutFile, "LOOKUP_TABLE default\n");
  for (GR_index_t iter = 0; iter < OutMesh.getNumVerts(); ++iter) {
    fprintf(pFOutFile, "%i\n", iter);
  }
  fprintf(pFOutFile, "\n");

  fclose(pFOutFile);
  logMessage(2, "Finished writing VTK dataFfile, %s\n", strFileName);
}
#endif

void
writeVTKLegacyWithoutPurge(const Mesh& OutMesh, const char strBaseFileName[])
{

//#ifndef NDEBUG

  const char* const strFileNameExt = "";

  FILE *pFOutFile = NULL;
  char strFileName[1024];
  if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
      || strlen(strBaseFileName) < 4) {
    // File name doesn't end in .vtk
    sprintf(strFileName, "%s%s.vtk", strBaseFileName, strFileNameExt);
  }
  else {
    // File name ends in .vtk
    sprintf(strFileName, "%s%s", strBaseFileName, strFileNameExt);
  }

  /*newfile mesh*/
  //logMessage(2, "Opening output file %s\n\n", strFileName);
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open VTK output file for writing",
		"unpurged mesh output");

  GR_index_t nVerts = OutMesh.getNumVerts();   // number of nodes
  GR_index_t nCells = OutMesh.getNumCells();

  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "Produced by GRUMMP\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", nVerts);

  //-------------------------------------
  // write vertex data
  //-------------------------------------
  Vert* pV = NULL;
  for (GR_index_t iV = 0; iV < nVerts; ++iV) {
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
  // Compute total connectivity size.
  GR_index_t connectSize = 0;
  GR_index_t realCells = 0;
  for (GR_index_t iC = 0; iC < nCells; iC++) {
    Cell *pC = OutMesh.getCell(iC);
    if (!pC->isValid() || pC->isDeleted())
      continue;
    connectSize += 1 + OutMesh.getCell(iC)->getNumVerts();
    realCells++;
  }

  fprintf(pFOutFile, "\n\nCELLS %d %d", realCells, connectSize);

  for (GR_index_t iC = 0; iC < nCells; ++iC) {
    Cell *pC = OutMesh.getCell(iC);
    if (!pC->isValid() || pC->isDeleted())
      continue;
    int h = pC->getNumVerts();
    fprintf(pFOutFile, "\n%d  ", h);
    for (int idx = 0; idx < h; ++idx) {
      GR_index_t globalNodeNum = OutMesh.getVertIndex(pC->getVert(idx));
      fprintf(pFOutFile, " %10d", globalNodeNum);
    }
  }
  fprintf(pFOutFile, "\n");

  //-------------------------------------
  // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  // VTK_PYRAMID = 14)
  //-------------------------------------
  fprintf(pFOutFile, "\nCELL_TYPES %d\n", realCells);
  for (GR_index_t iC = 0; iC < nCells; ++iC) {
    Cell *pC = OutMesh.getCell(iC);
    if (!pC->isValid() || pC->isDeleted())
      continue;
    int CT = OutMesh.getCell(iC)->getType();
    int VTK_type = 1;
    switch (CT)
      {
      case CellSkel::eTriCell:
	VTK_type = 5;
	break;
      case CellSkel::eQuadCell:
	VTK_type = 9;
	break;
      case CellSkel::eTet:
	VTK_type = 10;
	break;
      case CellSkel::ePyr:
	VTK_type = 14;
	break;
      case CellSkel::ePrism:
	VTK_type = 13;
	break;
      case CellSkel::eHex:
	VTK_type = 12;
	break;
      default:
	break;
      }
    fprintf(pFOutFile, "%d\n", VTK_type);
  }
  fprintf(pFOutFile, "\n");

  fclose(pFOutFile);

  logMessage(
      2,
      "-----------------------------------------------------------------------------\n");
  logMessage(2, "Finished writing VTK dataFile, %s\n", strFileName);
  logMessage(2, "The mesh has %9u vertices\n", OutMesh.getNumVerts());
  logMessage(2, "             %9u faces\n", OutMesh.getNumFaces());
  logMessage(2, "             %9u bdry faces\n", OutMesh.getNumBdryFaces());
  logMessage(2, "        and: %9u cells\n", OutMesh.getNumCells());
}

void
writeVTK_CellSet(std::vector<Cell*> & cell_list, const char strBaseFileName[])
{
  //-----------------------------------------------------------------------------

  FILE *pFOutFile = NULL;
  char strFileName[1024];
  if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
      || strlen(strBaseFileName) < 4) {
    // File name doesn't end in .vtk
    sprintf(strFileName, "%s.vtk", strBaseFileName);
  }
  else {
    // File name ends in .vtk
    sprintf(strFileName, "%s", strBaseFileName);
  }
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d Vtk mesh output");

  GR_index_t nno = 0;	// number of nodes
  GR_index_t ncells = cell_list.size();
  for (std::vector<Cell*>::iterator itc = cell_list.begin();
      itc != cell_list.end(); itc++) {
    nno += (*itc)->getNumVerts();
  }

  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "GRUMMP Tetra example\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", nno);

  //-------------------------------------
  // write 3d vertex data
  //-------------------------------------
  for (std::vector<Cell*>::iterator itc = cell_list.begin();
      itc != cell_list.end(); itc++) {
    const Cell * cell = *itc;
    for (int i = 0; i < cell->getNumVerts(); i++) {
      const Vert *pV = cell->getVert(i);
      fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
    }
  }

  int counter = -1;
  fprintf(pFOutFile, "CELLS %d %d\n", ncells, nno + ncells);
  for (std::vector<Cell*>::iterator itc = cell_list.begin();
      itc != cell_list.end(); itc++) {
    const Cell * cell = *itc;

    if (cell->getNumVerts() == 4)
      fprintf(pFOutFile, "4 %d %d %d %d\n", counter + 1, counter + 2,
	      counter + 3, counter + 4);
    else
      fprintf(pFOutFile, "3 %d %d %d\n", counter + 1, counter + 2, counter + 3);
    counter += cell->getNumVerts();
  }

  //-------------------------------------
  // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  // VTK_PYRAMID = 14)
  //-------------------------------------
  fprintf(pFOutFile, "CELL_TYPES %d\n", ncells);
  for (std::vector<Cell*>::iterator itc = cell_list.begin();
      itc != cell_list.end(); itc++) {

    if ((*itc)->getNumVerts() == 3)
      fprintf(pFOutFile, "5\n");
    else
      fprintf(pFOutFile, "10\n");
  }
  fclose(pFOutFile);
}
void
writeVTK_SingleFace(const Face * face, const char strBaseFileName[])
{
  FILE *pFOutFile = NULL;
  char strFileName[1024];
  if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
      || strlen(strBaseFileName) < 4) {
    // File name doesn't end in .vtk
    sprintf(strFileName, "%s.vtk", strBaseFileName);
  }
  else {
    // File name ends in .vtk
    sprintf(strFileName, "%s", strBaseFileName);
  }
  pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "3d Vtk mesh output");

  GR_index_t nno = 3;   // number of nodes

  //-------------------------------------------------------
  // Write the VTK header details
  //-------------------------------------------------------

  fprintf(pFOutFile, "# vtk DataFile Version 1.0\n");
  fprintf(pFOutFile, "GRUMMP Tetra example\n");
  fprintf(pFOutFile, "ASCII\n");
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(pFOutFile, "POINTS %d float\n", nno);

  //-------------------------------------
  // write 3d vertex data
  //-------------------------------------
  for (GR_index_t i = 0; i < nno; ++i) {
    const Vert *pV = face->getVert(i);
    fprintf(pFOutFile, "%16.8f %16.8f %16.8f\n", pV->x(), pV->y(), pV->z());
  }

  fprintf(pFOutFile, "CELLS 1 4\n");
  fprintf(pFOutFile, "3 0 1 2 \n");

  //-------------------------------------
  // write cell type (VTK_TRIANGLE = 5, VTK_QUAD = 9,
  // VTK_TETRA = 10, VTK_HEXAHEDRON = 12, VTK_WEDGE = 13,
  // VTK_PYRAMID = 14)
  //-------------------------------------
  fprintf(pFOutFile, "CELL_TYPES %d\n", 1);

  fprintf(pFOutFile, "5\n");

  fclose(pFOutFile);
}

// this is only for data viz, otherwise no purpose
void
writeVTKPointSet(std::vector<CubitVector> &cv_list, double * data,
		 const char strBaseFileName[], const char strExtraFileSuffix[])
{

  char strFileName[1024];
  if (strcasecmp(strBaseFileName + strlen(strBaseFileName) - 4, ".vtk")
      || strlen(strBaseFileName) < 4) {
    // File name doesn't end in .vtk
    sprintf(strFileName, "%s%s.vtk", strBaseFileName, strExtraFileSuffix);
  }
  else {
    // File name ends in .vtk
    sprintf(strFileName, "%s%s", strBaseFileName, strExtraFileSuffix);
  }
  FILE * pFOutFile = fopen(strFileName, "w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open output file for writing", "surface mesh output");

  /*"# vtk DataFile Version" discard_float*/
  fprintf(pFOutFile, "# vtk DataFile Version 1.1.1\n");
  /*discard_line*/
  fprintf(pFOutFile, "GRUMMP\n");
  /*ASCII*/
  fprintf(pFOutFile, "ASCII\n");
  /*DATASET UNSTRUCTURED_GRID*/
  fprintf(pFOutFile, "DATASET UNSTRUCTURED_GRID\n");
  /*POINTS nverts float*/
  fprintf(pFOutFile, "POINTS %u float\n",
	  static_cast<GR_index_t>(cv_list.size()));
  GR_index_t iPoints = cv_list.size();

  for (GR_index_t iV = 0; iV < iPoints; iV++)
    /*verts: coords*/
    fprintf(pFOutFile, "%.14g %.14g %.14g\n", cv_list[iV].x(), cv_list[iV].y(),
	    cv_list[iV].z());

  /*"CELLS" ncells */
  fprintf(pFOutFile, "CELLS %u %u\n", iPoints, iPoints * 2);

  for (GR_index_t iP = 0; iP < iPoints; iP++) {
    fprintf(pFOutFile, "1 %u\n", iP);
  }

  /*CELL_TYPES ncells*/
  fprintf(pFOutFile, "CELL_TYPES %u\n", iPoints);

  for (GR_index_t iC = 0; iC < iPoints; iC++) {
    /*cells: 5*/
    fprintf(pFOutFile, "1\n");
  }
  if (data) {
    fprintf(pFOutFile, "POINT_DATA %u\n", iPoints);
    fprintf(pFOutFile, "SCALARS Data_set float 1\n");
    fprintf(pFOutFile, "LOOKUP_TABLE default\n");
    for (GR_index_t iV = 0; iV < iPoints; iV++) {
      fprintf(pFOutFile, "%f\n", data[iV]);
    }
  }

  fclose(pFOutFile);
}
