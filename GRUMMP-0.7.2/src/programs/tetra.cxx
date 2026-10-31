#ifndef WIN32
#include <unistd.h>
#else /* WIN32 defined */
int getopt(int argc, char *argv[], char *opstring);
#endif

#include <stdlib.h>

#include <math.h>
#include <string.h>
#include "GR_AdaptPred.h"
#include "GR_Bdry3D.h"
#include "GR_InsertionQueue.h"
#include <GR_Length.h>
#include "GR_SurfMesh.h"
#include "GR_VolMesh.h"
#include "GR_misc.h"
#include "GR_SurfMeshBuilder.h"
#include "GR_TetMeshBuilder.h"
#include "GR_RefinementManager3D.h"
#include "GR_CoarsenManager.h"
#include "GR_SmoothingManager.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_SurfaceInsertion.h"
#include "GR_GraphicsOut.h"

#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "DLIList.hpp"
#include "GeometryQueryTool.hpp"
#include "FacetEvalTool.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetSurface.hpp"
#include "RefFace.hpp"

#ifdef SIM_ANNEAL_TEST
extern double simAnnealTemp;
#endif

static char GR_tetra_Usage[][80] =
  { "    -i    Base file name for input file", " ",
      "  Boolean options (argument of 1 for true, 0 for false):",
      "    -A    Use adaptive precision geometric predicates [0]",
      "    -b    Allow boundary changes [1]",
      "    -e    Allow edge swapping [1]",
      "    -m    Merge co-planar boundary patches as appropriate [1]", "  ",
      "  Numeric options",
      "    -a #  Max allowable angle for surface reconfiguration [5]",
      "    -c    Use diametral lenses instead of diametral circles [0]",
      "          0 -> circumradius < shortest edge length * 2 [good]",
      "          1 -> circumradius < shortest edge length * sqrt(2) [better]",
      "          2 -> circumradius < shortest edge length [best known result]",
      " ", "  Length scale options (values < 1 don't make sense):",
      "    -g #  Grading (rate of cell size change); larger is slower [1]",
      "    -r #  Resolution (feature size / cell size); larger is finer [1]",
      " ", "  Adaptation option: (meshopt2d only)",
      "    -l filename  Refine to length scale specified (per vertex) in file",
      " ", "    -o    Base file name for output file [input file base name]",
      " ", "  Mesh optimization string, in shorthand; see manual:",
      "    -O    [wsmf5]", " ", "  Surface Insertion Parameters",
      "    -s    [string]", " ", "    -q #  Mesh quality measure used [2]",
      "         0  Maximum dihedral angle",
      "         1  Minimum dihedral angle", "         2  All dihedral angles",
      "         3  Maximum solid angle", "         4  Minimum solid angle",
      "         5  All solid angles",
      "         6  Aspect ratio I (inscribed to circum radius)",
      "         7  Sliver quality measure",
      "         8  Skewness quality measure",
      "         9  Aspect ratio II (based on volume versus surface area)",
      "        10  Aspect ratio III (shortest edge to circumradius)", " ",
      "    -W (0,1) Output Temp Meshes - Useful for Debugging", "" // Must end with a null string.
    };

static void
vUsage(const char strCanonName[])
{
  logMessage(0, "Usage:  %s -i basefilename [options] \n", strCanonName);
  for (int i = 0; strlen(GR_tetra_Usage[i]) > 0; i++) {
    logMessage(0, "%s\n", GR_tetra_Usage[i]);
  }
  exit(1);
}

int
main(int iNArg, char *apcArgs[])
{
  char strBaseFileName[FILE_NAME_LEN], strQualFileName[FILE_NAME_LEN];
  char strBaseOutFileName[FILE_NAME_LEN];
  char strBdryFileName[FILE_NAME_LEN];
  bool qReadOutFileName = false;
  bool qReadAuxBdryFileName = false;
  bool qLengthScaleGiven = false;
  bool isSTLfile = false;
  bool qWriteTempMeshes = false;
  bool OMP_Flag = false; //********* added for enabling/disabling openmp parallel
  int NumOctLayers = 3; //********* added for setting the number of layers for Octtree
  char strLengthFileName[FILE_NAME_LEN];
  char strSTLFileName[FILE_NAME_LEN];
  int iBdryOptLevel = 1;

  char strOptParams[256], strSurfParams[256];
  char strExecName[FILE_NAME_LEN], strCanonicalName[20];
  enum {
    UNKNOWN, TETRA, MESHOPT, COARSEN, TRANSLATE, SURFINS, TMOPADAPT3D
  } iExec = UNKNOWN;

  sprintf(strExecName, "%s", apcArgs[0]);
  if (strstr(strExecName, "coarsen3d")) {
    iExec = COARSEN;
    sprintf(strCanonicalName, "coarsen3d");
  }
  else if (strstr(strExecName, "meshopt3d")) {
    iExec = MESHOPT;
    sprintf(strCanonicalName, "meshopt3d");
  }
  else if (strstr(strExecName, "trans3d")) {
    iExec = TRANSLATE;
    sprintf(strCanonicalName, "trans3d");
  }
  else if (strstr(strExecName, "surfins3d")) {
    iExec = SURFINS;
    sprintf(strCanonicalName, "surfins3d");
  }
  else if (strstr(strExecName, "tetra")) {
    iExec = TETRA;
    sprintf(strCanonicalName, "tetra");
  }
  else if (strstr(strExecName, "tmop_adapt3d")) {
    iExec = TMOPADAPT3D;
    sprintf(strCanonicalName, "tmop_adapt3d");
#ifndef HAVE_MESQUITE
    vFatalError(
	"TMOP-based adaptation requires Mesquite (version >= 2.99)\n",
	"tri.cxx:main()");
#endif
  }
  else {
    vFatalError("Invalid executable name.  Be less creative renaming files.\n",
		"tetra.cxx:main()");
  }

#ifdef SIM_ANNEAL_TEST
  simAnnealTemp = 10;
#endif

  // Set default values for arguments
  bool qDoEdgeSwapping = true;
  bool qAllowBdryChanges = true;
  double dGrading = 1;
  double dResolution = 1;
  double dMaxAngle = 5; // In degrees
  int iQualMeasure = 2;
  // Currently, only spheres are allowed for encroachment; lenses
  // weren't robust enough.
//	eEncroachType eET = eBall;
  // The following options are currently inactive in 3D.
  //    bool qAnisoCoarsening = false;
  // This one is global, declared in GR_AdaptPred.h
  qAdaptPred = false;
  bool qFile = false, qError = false, qParams = false, qSurfParams = false;
  int iPChar;
  while ((iPChar = getopt(iNArg, apcArgs,
			  "A:a:B:b:e:g:i:j:l:m:o:O:p:q:r:s:W:z:N:")) != EOF) {
    switch (iPChar)
      {
      case 'a':
	sscanf(optarg, "%lf", &dMaxAngle);
	if (dMaxAngle < 0)
	  dMaxAngle = 0;
	if (dMaxAngle > 180)
	  dMaxAngle = 180;
	break;
      case 'A':
	qAdaptPred = (optarg[0] == '1');
	break;
      case 'B':
	strncpy(strBdryFileName, optarg, FILE_NAME_LEN - 1);
	strBdryFileName[FILE_NAME_LEN - 1] = 0;
	qReadAuxBdryFileName = true;
	{
	  int iDum;
	  for (iDum = 0; iDum < FILE_NAME_LEN && strBdryFileName[iDum];
	      iDum++) {
	  }
	  if (iDum == FILE_NAME_LEN) // No null; file name too long
	    vFatalError("Boundary file name too long", "command line");
	}
	break;
      case 'b':
	qAllowBdryChanges = (optarg[0] == '1');
	break;
	//     case 'c':
	//       switch (optarg[0]) {
	//       case '0':
	// 	eET = eBall;
	// 	break;
	//       case '1':
	//       default:
	// 	eET = eLens;
	// 	break;
	//       case '2':
	// 	eET = eNewLens;
	// 	break;
	//       }
	//       break;
      case 'e':
	qDoEdgeSwapping = (optarg[0] == '1');
	break;
      case 'g':
	sscanf(optarg, "%lf", &dGrading);
	if (dGrading < 1)
	  dGrading = 1;
	break;
      case 'i':
	if (qFile)
	  qError = true;
	strncpy(strBaseFileName, optarg, FILE_NAME_LEN - 1);
	strBaseFileName[FILE_NAME_LEN - 1] = 0;
	char *baseNameEnd;
	baseNameEnd = strstr(strBaseFileName, ".stl");
	if (baseNameEnd != NULL) {
	  isSTLfile = true;
	  strcpy(strSTLFileName, strBaseFileName);
	  int len = baseNameEnd - strBaseFileName;
	  if (len > 1023)
	    len = 1023;
	  strncpy(strBaseFileName, strSTLFileName, len);
	  strBaseFileName[len] = 0;

	}

	{
	  GR_index_t iDum;
	  for (iDum = 0; iDum < FILE_NAME_LEN && strBaseFileName[iDum];
	      iDum++) {
	  }
	  if (iDum == FILE_NAME_LEN) // No null; file name too long
	    vFatalError("File name too long", "command line");
	}
	qFile = true;
	break;
      case 'j':
	//        qAnisoCoarsening = (optarg[0] == '1');
	break;
      case 'l':
	if (iExec == MESHOPT || iExec == TETRA) {
	  strncpy(strLengthFileName, optarg, FILE_NAME_LEN - 1);
	  strLengthFileName[FILE_NAME_LEN - 1] = 0;
	  {
	    int iDum;
	    for (iDum = 0; iDum < FILE_NAME_LEN && strLengthFileName[iDum];
		iDum++) {
	    }
	    if (iDum == FILE_NAME_LEN) // No null; file name too long
	      vFatalError("Lengthscale file name too long", "command line");
	  }
	  qLengthScaleGiven = true;
	}
	else {
	  qError = true; // Be sure to crash out.
	}
	break;
      case 'm':
	sscanf(optarg, "%d", &iBdryOptLevel);
	if (iBdryOptLevel != 0 && iBdryOptLevel != 1)
	  qError = true;
	break;
      case 'o':
	strncpy(strBaseOutFileName, optarg, FILE_NAME_LEN - 1);
	strBaseOutFileName[FILE_NAME_LEN - 1] = 0;
	{
	  int iDum;
	  for (iDum = 0; iDum < FILE_NAME_LEN && strBaseFileName[iDum];
	      iDum++) {
	  }
	  if (iDum == FILE_NAME_LEN) // No null; file name too long
	    vFatalError("Output file name too long", "command line");
	}
	qReadOutFileName = true;
	break;
      case 'O':
	sscanf(optarg, "%255s", strOptParams);
	qParams = true;
	break;
      case 'p':
	int ithreadE;
	sscanf(optarg, "%d", &ithreadE);
	if (ithreadE == 0)
	  OMP_Flag = false;
	else
	  OMP_Flag = true;
	break;
      case 'N':
	int NOLayer;
	sscanf(optarg, "%d", &NOLayer);
	if (NOLayer < 2) {
	  vFatalError("Number of layers should be at least 2", "command line");
	  assert(0);
	}
	NumOctLayers = NOLayer;
	break;
      case 'q':
	sscanf(optarg, "%d", &iQualMeasure);
	break;
      case 'r':
	sscanf(optarg, "%lf", &dResolution);
	break;
      case 's':
	sscanf(optarg, "%255s", strSurfParams);
	qSurfParams = true;
	break;
      case 'W':
	qWriteTempMeshes = (optarg[0] == '1');
	break;
      case 'z':
	strncpy(strSTLFileName, optarg, FILE_NAME_LEN - 1);
	strSTLFileName[FILE_NAME_LEN - 1] = 0;
	int iDum;
	for (iDum = 0; iDum < FILE_NAME_LEN && strSTLFileName[iDum]; iDum++) {
	}
	if (iDum == FILE_NAME_LEN) // No null; file name too long
	  vFatalError("Surface file name too long", "command line");

	break;
      default:
	qError = true;
	break;
      }
  }

  if (qError || !qFile) {
    vUsage(strCanonicalName);
  }
  if (!qReadOutFileName) {
    strncpy(strBaseOutFileName, strBaseFileName, FILE_NAME_LEN - 1);
  }

  // Do this after option parsing so that the usage message doesn't have
  // other stuff around it.
  GRUMMPInit(strCanonicalName);
  // Open log file and name quality file
  openMessageFile(strBaseFileName);
  makeFileName(strQualFileName, "%s.qual", strBaseOutFileName,
	       "main() [tetra.C]");
  logMessage(1, "Mesh quality output file: %s\n", strQualFileName);
  // small change, lets log the call to generate this function, to make repeating results possible
  for (int iArg = 0; iArg < iNArg; iArg++) {
    logMessage(1, "%s ", apcArgs[iArg]);
  }
  logMessage(1, "\n");
  VolMesh *pVMesh;
  SurfMeshBuilder* SMB = NULL;
  TetMeshBuilder* TMB = NULL;
  Bdry3D *pB3D = NULL;
  if (qReadAuxBdryFileName) {
    pB3D = new Bdry3D(strBaseFileName, iBdryOptLevel);
  }
  std::set<Subseg*> * mesh_subsegs = NULL;

  if (TETRA == iExec) {
    // tetra initialization code

    if (!isSTLfile) {
      // Read, validate, and optimize boundary rep
			pVMesh = createVolMeshFromBdryFile(strBaseFileName, iBdryOptLevel,
					iQualMeasure, dResolution, dGrading, pB3D);
    }
    else {
      CGMApp::instance()->startup(iNArg, apcArgs);
      CubitBoolean use_feature_angle = CUBIT_TRUE;
      double feature_angle = 135;
      double tolerance = 1.e-6;
      int interp_order = 4;
      CubitBoolean smooth_non_manifold = CUBIT_TRUE;
      CubitBoolean split_surfaces = CUBIT_FALSE;
      CubitBoolean stitch = CUBIT_TRUE;
      CubitBoolean improve = CUBIT_TRUE;
      DLIList<CubitQuadFacet*> quad_facet_list;
      DLIList<CubitFacet*> tri_facet_list;
      DLIList<Surface*> surface_list;

      //This is the call to CGM to build the geometry from
      //a facetted surface representation.
      GeometryQueryTool *gqt = GeometryQueryTool::instance();
      FacetQueryEngine *fqe = FacetQueryEngine::instance();
      fqe->import_facets(strSTLFileName, use_feature_angle, feature_angle,
			 tolerance, interp_order, smooth_non_manifold,
			 split_surfaces, stitch, improve, quad_facet_list,
			 tri_facet_list, surface_list, STL_FILE);
      GraphicsOut go_surf("CGMsurf.vtk");
      GraphicsOut go_curv("CGMcurv.vtk");
      go_surf.output_faces();
      go_curv.output_edges();
      //Set the compare tolerance to a lower value for all FacetSurface.
      //This is needed to project points onto the surface interpolated from facets.

      // this is here for now, because its not set up to handle multiple volumes,
      // but the capability now exists

      assert(gqt->num_ref_volumes() == 1);
      RefVolume * refVolume = gqt->get_first_ref_volume();
      DLIList<RefFace*> all_surfaces;
      refVolume->ref_faces(all_surfaces);
      int num_surfaces = all_surfaces.size();
      RefFace* surface;
      FacetSurface* facet_surface;
      for (int i = 0; i < num_surfaces; ++i) {
	surface = all_surfaces.get_and_step();
	facet_surface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());
	assert(facet_surface);
	facet_surface->get_eval_tool()->compare_tol(1.e-6);
      }

      double allowed_TVT = 0.45 * M_PI;

      SMB = new SurfMeshBuilder(refVolume, SurfMeshBuilder::FACET, allowed_TVT);
      SMB->sample_all_surfaces();
      SMB->output_stitched_surface("final-sampled.vtk");
      mesh_subsegs = new std::set<Subseg*>();
      TMB = new TetMeshBuilder(SMB);
      TMB->buildMesh(mesh_subsegs);
      pVMesh = TMB->get_mesh();
//			std::vector<CubitVector> pointList;
//			for(int iB = 0; iB < pVMesh->getNumBdryFaces(); iB++){
//				TriBFace * bface = dynamic_cast<TriBFace*>(pVMesh->getBFace(iB));
//				double split[3];
//				bface->calcSplitPoint(split);
//				pointList.push_back(split);
//
//			}
//			writeVTKPointSet(pointList,NULL,"projections","");
//			abort();

    }
  }
  else {
    // meshopt3d, trans3d and coarsen3d initialization code
    if (strstr(strBaseFileName, ".hin")) {
      pB3D = new Bdry3D(strBaseFileName);
      SurfMesh SM(pB3D, iQualMeasure);
      writeVTKLegacyWithoutPurge(SM, "tempSurf");
      assert(SM.isValid());
      // Create an initial tetrahedralization of that surface mesh
      pVMesh = new VolMesh(SM, iQualMeasure, 1. / dResolution, 1. / dGrading);

    }
    else {
      // Read volume mesh
      pVMesh = new VolMesh(iQualMeasure, dMaxAngle);
      pVMesh->setBdryData(pB3D); // Arg could be null
      pVMesh->readFromFile(strBaseFileName);
    }
    writeVTKLegacy(*pVMesh, strBaseFileName);
//		pVMesh->allowBdryChanges();
//		pVMesh->allowEdgeSwapping();
//		GRUMMP::MinMaxDihedSwapDecider3D * ttt = new GRUMMP::MinMaxDihedSwapDecider3D(true);
//		GRUMMP::SwapManager3D *sss = new GRUMMP::SwapManager3D(ttt,pVMesh);
//		for(int i = 0; i < 50; i++){
//			sss->swapAllFaces();
//			sss->swapAllBdryEdges();
//			pVMesh->purgeAllEntities();
//		}
//		delete ttt;
//		delete sss;
//		writeVTKLegacy(*pVMesh, "initialMeshSwap");

    logMessage(1, "Number of cells read: %u\n", pVMesh->getNumCells());

    pVMesh->setWriteTempMesh(qWriteTempMeshes);
  }
  // safer to start with
  if (SURFINS == iExec) {
    qDoEdgeSwapping = false;
    qAllowBdryChanges = false;

  }

  if (qDoEdgeSwapping)
    pVMesh->allowEdgeSwapping();
  else
    pVMesh->disallowEdgeSwapping();

  if (qAllowBdryChanges)
    pVMesh->allowBdryChanges();
  else
    pVMesh->disallowBdryChanges();

  switch (iExec)
    {

    case TMOPADAPT3D:
      {

	/*
	 * Extension of the 2D Aniso edge adaptation with TMOP. This is still an unfinished code.
	 * The main problem of the algorithm was in swapping cells. Two different quality criterias are
	 * used for swapping: SineDihedralAngles and VolumeLength3D. VolumeLentgh is used to force swapping cells
	 * by improving the best cell into a better quality, which also creates a lot of bad cells. SineDihedral angles
	 * is used to fix as much of the bad cells created by the previous swapper. This gave good alignment but only at surfaces,
	 * and at interior cells it created a lot of bad cells.
	 * For this reason both swapping criteria are set to improve the worst quality (as it should be... assuring bad cells are not created), but
	 * as a result, good alignment isn't obtained even at the surface. A solution to fix this is still needed.
	 */

	vFatalError("Status: undergoing, not yet a reliable approach for 3D...",
		    "TMOPadapt3D");

#ifdef HAVE_MESQUITE

// Setting metric for test cases. Only supports analytic metric test cases at the moment.
	AnisoRefinement3D *pAR = new AnisoRefinement3D(pVMesh);
	pAR->vAssignMetricsToVertsTest(pAR->e3DTestQuarterCircle);
	double dMinMLength = sqrt(2) / 2.;
	double dMaxMLength = sqrt(2);

// Setting both criteria for swapping
	GRUMMP::QualMeasure *pQM = new GRUMMP::VolumeLength3D();
	pQM->vEvalInMetric(true);
	GRUMMP::SwapDecider3D* pSD3D = new GRUMMP::MetricSwapDecider3D(true,
								       true,
								       pQM);
	GRUMMP::SwapManager3D *pSM3D = new GRUMMP::SwapManager3D(pSD3D, pVMesh);

	GRUMMP::QualMeasure *pQM_2 = new GRUMMP::SineDihedralAngles();
	pQM_2->vEvalInMetric(true);
	GRUMMP::SwapDecider3D* pSD3D_2 = new GRUMMP::MetricSwapDecider3D(true,
									 true,
									 pQM_2);
	GRUMMP::SwapManager3D *pSM3D_2 = new GRUMMP::SwapManager3D(pSD3D_2,
								   pVMesh);

// Used for insertion, one is for computing the location of the point to be inserted and the other one call the insertin function
	GRUMMP::InsertionPointCalculator *pIPC = new GRUMMP::AnisoEdgeIPC3D();
	GRUMMP::SwappingInserter3D *pSI3D = new GRUMMP::SwappingInserter3D(
	    pVMesh, NULL);

// For smoothing
	GRUMMP::SmoothingManagerMesquite3D *pSMM3D =
	    new GRUMMP::SmoothingManagerMesquite3D(pVMesh, false);

	GRUMMP::AnisoEdgeRefinementManager3D ARM3D(pVMesh, dMaxMLength,
						   dMinMLength, pIPC, pSMM3D,
						   pAR, pSM3D, pSM3D_2, pSI3D);

	char strAnisoFileName[60];
	sprintf(strAnisoFileName, "Aniso_Edge_Adapted_Mesh_3D");

	ARM3D.AnisoRefineMesh();

	writeVTKLegacy(*pVMesh, strAnisoFileName);
	writeNativeMESH(*pVMesh, strAnisoFileName, "");
	writeNativeGRMESH(*pVMesh, strAnisoFileName, "");

	delete pAR;
	delete pSI3D;
	delete pSMM3D;
	delete pIPC;
	delete pSM3D_2;
	delete pSD3D_2;
	delete pQM_2;
	delete pSM3D;
	delete pSD3D;
	delete pQM;

#endif
	break;
      }

    case MESHOPT:
      dResolution /= sqrt(2.);
      if (!qParams) {
	sscanf("wsmf5mf5rmf5mf5", "%255s", strOptParams);
	qParams = true;
      }
      if (qLengthScaleGiven) {
	/*
	 pVMesh->initLengthScale(1./dResolution, 1./dGrading);
	 // Read it and store it.
	 FILE *pFLenFile = fopen(strLengthFileName, "r");
	 if (NULL == pFLenFile)
	 vFatalError("Couldn't open file for reading",
	 "length scale input");
	 double dLen;
	 int iV, iDataRead = fscanf(pFLenFile, "%d%lf", &iV, &dLen);
	 int iNV = pVMesh->getNumVerts();
	 while (iDataRead == 2) {
	 if (iV < 0 || iV > iNV) {
	 vFatalError("Vertex index out of range",
	 "length scale input");
	 }
	 Vert *pV = pVMesh->getVert(iV);
	 pV->setLengthScale(dLen);
	 pVMesh->markForLengthScaleCheck(pV);
	 logMessage(4, "Set length scale for vert %d.\n", iV);
	 iDataRead = fscanf(pFLenFile, "%d%lf", &iV, &dLen);
	 };
	 // Now use it.
	 pVMesh->updateLengthScale();
	 */

	pVMesh->setInterpolatedLengthScale();

	std::shared_ptr<GRUMMP::Length> pL3D = GRUMMP::Length::Create(
	    pVMesh, dResolution, dGrading);
	// TODO: Restore the read-length-scale-from-file functionality.
//			pL3D->provideMesh(pVMesh, strLengthFileName);

	//L3D.initLength(pVMesh);
	//pVMesh->adaptToLengthScale(eET);

	// This stuff is from Adapt.cxx, should put it back there with the CoarsenManager
	GRUMMP::CoarsenManager CM3D(pVMesh, pL3D);

	bool qResult = pVMesh->makeDelaunay();
	if (!qResult) {
	  vWarning("Couldn't make mesh fully Delaunay before coarsening!\n");
	  vWarning("This MIGHT still turn out okay...\n");
	}
	CM3D.coarsen(1);
	qResult = pVMesh->makeDelaunay();
	if (!qResult) {
	  vWarning("Couldn't make mesh fully Delaunay before coarsening!\n");
	  vWarning("This MIGHT still turn out okay...\n");
	}

	GRUMMP::OptMSSmoothingManager* OMS =
	    GRUMMP::OptMSSmoothingManager::Create(pVMesh);
	if (OMS) {
	  OMS->smoothAllVerts(2);
	  delete OMS;
	}
      }
      pVMesh->classifyTets();
      pVMesh->evaluateQuality();
      break;
    case TETRA:
      {
	assert(pVMesh->isValid());

	std::shared_ptr<GRUMMP::Length> sizing_field = GRUMMP::Length::Create(
	    pVMesh, dResolution, dGrading);
	if (qLengthScaleGiven) {

	  FILE *pFLenFile = fopen(strLengthFileName, "r");
	  if (NULL == pFLenFile)
	    vFatalError("Couldn't open file for reading", "length scale input");

	  double x_loc = -LARGE_DBL, y_loc = -LARGE_DBL, z_loc = -LARGE_DBL,
	      length = -LARGE_DBL;

	  while (1) {
	    int iDataRead = fscanf(pFLenFile, "%lf%lf%lf%lf", &(x_loc),
				   &(y_loc), &(z_loc), &(length));
	    if (iDataRead != 4)
	      break;
	    sizing_field->addToLengthMesh(x_loc, y_loc, z_loc, length);
	  }
	  fclose(pFLenFile);
	}

	else {
	  pVMesh->purgeAllEntities();
	  //Refines the mesh with an automatic sizing field.

	  // assert(pVMesh->getSubsegMap());
	  // sizing_field->set_subseg_map(pVMesh->getSubsegMap());
	}

	GRUMMP::CircumcenterIPC IPC;
	GRUMMP::RefinementManager3D RM3D(pVMesh, &IPC, mesh_subsegs, eBall,
					 sizing_field,
					 OldPriorityCalculator::EDGE_RADIUS);
	//		RM3D.onlyRefineForQuality();
	RM3D.Set_Refine3D_Timer();
	RM3D.refineMesh(OMP_Flag, NumOctLayers);
	//RM3D.refineForQuality(OMP_Flag, NumOctLayers);
	RM3D.Print_Refine3D_Timer(1);

	if (mesh_subsegs)
	  delete mesh_subsegs;
      }
      pVMesh->classifyTets();
      pVMesh->evaluateQuality();
      if (!qParams)
	sscanf("wsmf5mf5", "%255s", strOptParams);
      break;
    case COARSEN:
      {
	// give COARSEN its own smoother
	// previously smoother defined earlier, which let to a lot
	// of observation
	GRUMMP::OptMSSmoothingManager3D *coarsenSmoother =
	    new GRUMMP::OptMSSmoothingManager3D(pVMesh);
	// The constant in the length scale initialization is an empirical
	// constant, chosen to approximately match the edge-length
	// statistics and to get about a factor of 8 reduction in mesh size.
	logMessage(2, "Improving mesh quality before coarsening...\n");
	{
	  GRUMMP::SwapDecider3D *pSD3D = new GRUMMP::MaxMinSineSwapDecider3D(
	      qAllowBdryChanges);
	  GRUMMP::SwapManager3D Swapper(pSD3D, pVMesh);
	  Swapper.swapAllFaces();
	  delete pSD3D;
	}
	coarsenSmoother->smoothAllVerts(2);

	pVMesh->evaluateQuality();
	pVMesh->iSwap_deprecated(1);

	delete coarsenSmoother;

	//pVMesh->initLengthScale(1.2/dResolution, 1./dGrading);

	// The following line leads to minor disasters in terms of (for
	//    instance) applying T22 swaps on the bdry.
	//    pVMesh->vSetStrictPatchChecking();

	// Set the initial optimization-based smoothing threshold fairly
	// low, since the mesh is likely to be fairly poor at that point.
	//		coarsenSmoother->setSmoothingThreshold(10.);

	//dResolution doesn't actually do anything here
	std::shared_ptr<GRUMMP::Length> sizing_field = GRUMMP::Length::Create(
	    pVMesh, dResolution, dGrading);
	//Coarsening constant 3 chosen as similar to behaviour of old coarsen.
//		sizing_field.coarsen(3);
//		sizing_field.initLength(pVMesh);
//		pVMesh->coarsenToLengthScale();

	GRUMMP::CoarsenManager CM3D(pVMesh, sizing_field);
	CM3D.coarsen(2.6);

	logMessage(1, "Number of tets after coarsening: %u\n",
		   pVMesh->getNumCells());

	pVMesh->classifyTets();
	pVMesh->evaluateQuality();

	if (!qParams) {
	  sscanf("wsmf5mf5wsrmf5mf5", "%255s", strOptParams);
	  qParams = true;
	}
      }
      break;
    case SURFINS:
      logMessage(1, "\n\nComputing Surface Insertion \n \n");
      {
	CGMApp::instance()->startup(iNArg, apcArgs);
	if (!qParams) {
	  sscanf("wirwsmf5mf5wsmf5mf5", "%255s", strOptParams);
	}
	if (!qSurfParams) {
	  // This gives a properly null-terminated string of
	  // length 0.
	  strSurfParams[0] = 0;
	}

	std::shared_ptr<GRUMMP::Length> len = GRUMMP::Length::Create(
	    pVMesh, dGrading, dResolution);
	GRUMMP::SurfaceInserter3D *pSurfIns3D = new GRUMMP::SurfaceInserter3D(
	    pVMesh, len, dGrading, dResolution, strSurfParams);
	pSurfIns3D->insertSurfaceInMesh(strSTLFileName, strOptParams);

	assert(pSurfIns3D != NULL);
	delete pSurfIns3D;

	//			GRUMMP::DelaunaySwapDecider3D * ttt = new GRUMMP::DelaunaySwapDecider3D(true);
	//			GRUMMP::SwapManager3D *sss = new GRUMMP::SwapManager3D(ttt,pVMesh);
	//		for(int i = 0; i < 50; i++){
	//			sss->swapAllFaces();
	//
	//			pVMesh->purgeAllEntities();
	//		}
	//			 bool *aqDelaunay = new bool[pVMesh->getNumCells()];
	//int num = identifyNonDelaunayCells(*pVMesh,aqDelaunay);
	//			printf("mesh is %d delaunay\n", num);
	//
	//					delete aqDelaunay;
      }
      break;
    case TRANSLATE:
      // Do nothing
      break;
    default:
      break;
    }

  GRUMMP::OptMSSmoothingManager3D *pSmoother =
      new GRUMMP::OptMSSmoothingManager3D(pVMesh);

  if (iExec != TRANSLATE && iExec != SURFINS) {
    GR_index_t iLen = strlen(strOptParams), i = 0;

    if (qDoEdgeSwapping)
      pVMesh->allowEdgeSwapping();
    else
      pVMesh->disallowEdgeSwapping();

    while (i < iLen) {
      char cType = strOptParams[i++];
      assert(cType);
      switch (cType)
	{
	case 'W':
	case 'w':  // Do some swapping
	  {
	    char cMeasure = strOptParams[i++];
	    pVMesh->allowSwapRecursion();
	    GRUMMP::SwapDecider3D* pSD3D;
	    switch (cMeasure)
	      {
	      case 'I':
	      case 'i': // Insphere criterion
		logMessage(1, "Doing Delaunay face swapping...\n");
		pVMesh->setSwapType(eDelaunay);
		pSD3D = new GRUMMP::DelaunaySwapDecider3D(qAllowBdryChanges);
		break;
	      case 'M':
	      case 'm': // Minmax dihedral criterion
		logMessage(
		    1,
		    "Doing face swapping to reduce maximum dihedral angle...\n");
		pVMesh->setSwapType(eMaxDihed);
		pSD3D = new GRUMMP::MinMaxDihedSwapDecider3D(qAllowBdryChanges);
		break;
	      case 'S':
	      case 's': // Maxmin sine criterion
		logMessage(
		    1,
		    "Doing face swapping to increase minimum sine of dihedral angle...\n");
		pVMesh->setSwapType(eMinSine);
		pSD3D = new GRUMMP::MaxMinSineSwapDecider3D(qAllowBdryChanges);
		break;
	      case 'U':
	      case 'u': // uniform degree
		logMessage(
		    1, "Doing face swapping to equalize vertex degree...\n");
		pSD3D = new GRUMMP::UniformDegreeSwapDecider3D(
		    qAllowBdryChanges);
		break;
	      case 'q':
	      case 'Q':
		logMessage(
		    1,
		    "Doing face swapping based on user-specified quality criterion...\n");
		pVMesh->setSwapType(eFromQual);
		break;
	      default:
		vFatalError("Bad option for swapping criterion", "main()");
	      }
	    GRUMMP::SwapManager3D *pSwapper = new GRUMMP::SwapManager3D(pSD3D,
									pVMesh);
	    //	clock_t clBegin = clock();
	    pSwapper->Set_Swap3D_Timer();
	    int iNSwaps = pSwapper->swapAllFaces(OMP_Flag);
	    pSwapper->Print_Swap3D_Timer(1);

	    //    clock_t clEnd = clock();

#ifdef SIM_ANNEAL_TEST
	    simAnnealTemp = 1.e-10;
	    iNSwaps += pSwapper->swapAllFaces();
	    simAnnealTemp = 1;
#endif
	    delete pSwapper;
	    delete pSD3D;

	    logMessage(1, "Number of swaps: %d\n", iNSwaps);
	    logMessage(1, "Number of tets after optimization:  %u\n",
		       pVMesh->getNumCells());
	    //	vMessage(1, "Number of clock ticks: %d\n", clEnd - clBegin);
	    break;
	  }
	case 'M':
	case 'm': // Do some smoothing
	  {
	    pVMesh->disallowBdryChanges();
	    char cTech = strOptParams[i++];
	    int iFunc = atoi(&strOptParams[i++]);
	    pSmoother->setSmoothingTechnique(cTech);
	    if (iFunc == 0)
	      iFunc = 10;
	    pSmoother->setSmoothingGoal(iFunc);
	    logMessage(1, "Smoothing vertices...\n");
	    pSmoother->smoothAllVerts(1);
	    pVMesh->allowBdryChanges();

	  }
	  break;
	case 'N':
	case 'n':
	  // No-op
	  break;
	case 'R':
	  pVMesh->repairBadCells(true);
	  break;
	case 'r':
#ifdef SIM_ANNEAL_TEST
	  simAnnealTemp = 1;
#endif
	  pVMesh->repairBadCells(false);
	  break;
	case 'B':
	  pVMesh->breakBadBdryCells(true);
	  break;
	case 'b':
	  pVMesh->breakBadBdryCells(false);
	  break;
	default:
	  vFatalError("Bad option for mesh optimization operation", "main()");
	  break; // Previous line kills code anyway; this makes the compiler happy.
	}
      pVMesh->classifyTets();
      pVMesh->evaluateQuality();
    }
  }

  pVMesh->purgeAllEntities();
  pVMesh->writeQualityToFile(strQualFileName);

  //   if ((MESHOPT != iExec) || qLengthScaleGiven) {
  //     pVMesh->vCheckEdgeLengths();
  //   }
  assert(pVMesh->isValid());
  logMessage(1, "Writing volume mesh file...\n");
  /////
//    int nnn = 0;
//    for (GR_index_t iV = 0; iV < pVMesh->getNumVerts();iV++){
//        Vert * v = pVMesh->getVert(iV);
//                if(v->getNumFaces() == 0) {
//                   nnn++;
//        }
//    }
//    printf("Number of verts without connectivity %d\n",nnn);

  ////
  pVMesh->reorderAllEntities();

  writeNative(*pVMesh, strBaseOutFileName);
  writeVTKLegacy(*pVMesh, strBaseOutFileName);
  // The following line is useful for setting up a mesh into which the
  // shock tracking app can insert surfaces.
//  writeFileTetGen(*pVMesh, strBaseOutFileName, "", 0, nullptr);

  //   // Output: distribution of face degree of vertices.  (Interior only
  //   // for now).
  //   {
  //     int aiCount[100];
  //     for (int ii = 0; ii < 100; ii++) {
  //       aiCount[ii] = 0;
  //     }
  //     for (int iV = 0; iV < pVMesh->iNumVerts(); iV++) {
  //       Vert *pV = pVMesh->pVVert(iV);
  //       if (pV->qIsBdryVert()) continue;
  //       int iDeg = pV->iNumFaces();
  //       if (iDeg >= 100) iDeg = 99;
  //       aiCount[iDeg]++;
  //     }
  //     for (int ii = 0; ii < 100; ii++) {
  //       if (aiCount[ii] != 0) {
  // 	vMessage(2, "Degree: %2d.  Count: %6d\n", (ii+6)/3, aiCount[ii]);
  //       }
  //     }
  //   }
  //   {
  //     GRUMMP::SwapDecider3D *pSD3D =
  //       new GRUMMP::UniformDegreeSwapDecider3D(qAllowBdryChanges);
  //     GRUMMP::SwapManager3D Swapper(pSD3D, pVMesh);
  //     int iNSwaps = pSwapper->swapAllFaces();
  //     vMessage(1, "Number of swaps: %d\n", iNSwaps);
  //     iNSwaps = pSwapper->swapAllFaces();
  //     vMessage(1, "Number of swaps: %d\n", iNSwaps);
  //     iNSwaps = pSwapper->swapAllFaces();
  //     vMessage(1, "Number of swaps: %d\n", iNSwaps);
  //     vMessage(1, "Number of tets after optimization:  %u\n",
  // 	     pVMesh->iNumCells());
  //   }
  //   // Output: distribution of face degree of vertices.  (Interior only
  //   // for now).
  {
    int aiCount[100];
    for (int ii = 0; ii < 100; ii++) {
      aiCount[ii] = 0;
    }
    for (size_t iV = 0; iV < pVMesh->getNumVerts(); iV++) {
      Vert *pV = pVMesh->getVert(iV);
      if (pV->isBdryVert())
	continue;
      int iDeg = pV->getNumFaces();
      if (iDeg >= 100)
	iDeg = 99;
      aiCount[iDeg]++;
    }
    for (int ii = 0; ii < 100; ii++) {
      if (aiCount[ii] != 0) {
	logMessage(2, "Degree: %2d.  Count: %6d\n", (ii + 6) / 3, aiCount[ii]);
      }
    }
  }
  //   // Extra diagnostic stuff, to test the theory that bad tets cluster at
  //   // points with extremal degree.
  //   for (int ii = 0; ii < pVMesh->iNumCells(); ii++) {
  //     Cell *pC = pVMesh->pCCell(ii);
  //     double dQual = pVMesh->dEvaluateQuality(pC);
  //     if ((pC->pVVert(0)->iVertType() != Vert::eInterior) ||
  // 	(pC->pVVert(1)->iVertType() != Vert::eInterior) ||
  // 	(pC->pVVert(2)->iVertType() != Vert::eInterior) ||
  // 	(pC->pVVert(3)->iVertType() != Vert::eInterior)) continue;
  //     int minDeg = 100, maxDeg = 0;
  //     for (int iV = 0; iV < 4; iV++) {
  //       int deg = pC->pVVert(iV)->iNumFaces();
  //       if (deg < minDeg) minDeg = deg;
  //       if (deg > maxDeg) maxDeg = deg;
  //     }
  //     vMessage(2, "Cell %5d.  Qual = %10.6f. Vert degrees: %2d %2d.\n",
  // 	     ii, dQual, minDeg, maxDeg);
  //   }

  if (iExec == TETRA) {
//		 this is a bit messed up
    logMessage(1, "  Final verts: %u. Final tets: %u.\n", pVMesh->getNumVerts(),
	       pVMesh->getNumCells());
  }

  delete pSmoother;
  delete pVMesh;
  if (isSTLfile) {
    delete TMB;
    delete SMB;
//    FacetQueryEngine *fqe = FacetQueryEngine::instance();
//    fqe->delete_instance();
    GeometryQueryTool *gqt = GeometryQueryTool::instance();
    gqt->delete_geometry();
    gqt->delete_instance();
    CGMApp::instance()->shutdown();
    CGMApp::delete_instance();
  }
  if (pB3D)
    delete pB3D;

  logMessage(1, "Done\n");

  return (0);
}
