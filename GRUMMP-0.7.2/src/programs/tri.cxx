#ifndef WIN32
#include <unistd.h>
#else /* WIN32 defined */
int getopt(int argc, char *argv[], char *opstring);
extern char * optarg;
#endif

#include <cstdlib>
#include <math.h>
#include <cstring>
#include <ctime>
#include <sys/time.h>

#include "GR_config.h"
#include "GR_misc.h"
#include "GR_AdaptPred.h"
#include "GR_GRGeom2D.h"
#include "GR_GRPoint.h"
#include "GR_GRCurve.h"
#include "GR_GRCurveGeom.h"
#include "GR_InsertionManager.h"
#include "GR_CoarsenManager.h"
#include "GR_Length.h"
#include "GR_Mesh2D.h"
#include "GR_Mesh.h"
#include "GR_TriMeshBuilder.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_SmoothingManager.h"
#include "GR_SurfaceInsertion.h"
typedef struct _refinementInfo_ {
	double dX, dY, dLen;
} RefinementInfo;

static char GR_tri_Usage[][80] =
		{ "    -i    Base file name for input file", " ",
				"  Boolean options (argument of 1 for true, 0 for false):",
				"    -A    Use adaptive precision geometric predicates [0]",
				"    -b    Allow boundary changes [1]",
				"    -c    Use diametral lenses instead of diametral circles [1]",
				"          This will push the minimum angle to 30 degrees (instead of 20)",
				"    -j    Coarsen quasi-structured mesh anisotropically [0]",
				"    -t    Use Lex reordering of vertices instead of default RCM",
				" ", "  Length scale options (values < 1 don't make sense):",
				"    -a #  Minimum Angle used for Refinement (Degrees) [1]",
				"    -g #  Grading (rate of cell size change); larger is slower [1]",
				"    -G    Use global smoothing with TMOP anisotropic adaptation",
				"    -r #  Resolution (feature size / cell size); larger is finer [1]",
				"    -m #  Specifies a minimum length scale [none]", " ",
				"  Adaptation option: (meshopt2d and tri only)",
				"    -l filename  Refine to length scale specified (per vertex) in file",
				" ",
				"    -o filename  Base file name for output file [input file base name]",
				"  Mesh optimization string, in shorthand; used for surface Insertion",
				"    -O    [wsmf5]", " ",
				"    -q #  Mesh quality measure used [2]",
				"         0  Maximum angle", "         1  Minimum angle",
				"         2  All angles",
				"         6  Aspect ratio I (inscribed to circum radius)",
				"         9  Aspect ratio II (based on area versus perimeter)",
				"        10  Aspect ratio III (shortest edge to circumradius)",
				"    -x # Point insertion method used",
				"         0  Circumcenter",
				"         1  Offcenter (allows for threshold right after)",
				"         2  Sink              ",
				"         3  Optimized Offcenter     ",
				"         4  Engwirda's Offcenter",
				"    -X # Point insertion method, extra parameter",
				"    -Q Face-based queue[1] or cell-based queue[0]",
				"    -D Parameter for the termination criteria[1.]", " ",
				"  Surface Insertion Options , pass as string",
				"    -s    [string]",
				"    -I   Melt Internal boundary faces during surface insertion",
				"    -P    Output boundary patch file (0 = no, 1 = yes) [0]",
				"    -W    Write temporary meshes at each unit operation ", " ",
				"" // Must end with a null string.
		};

static void vUsage(const char strCanonName[]) {
	logMessage(0, "Usage:  %s -i basefilename [options] \n", strCanonName);
	for (int i = 0; strlen(GR_tri_Usage[i]) > 0; i++) {
		logMessage(0, "%s\n", GR_tri_Usage[i]);
	}
	exit(1);
}

int main(int iNArg, char *apcArgs[]) {
	clock_t timecheck;
	char strBaseFileName[FILE_NAME_LEN], strQualFileName[FILE_NAME_LEN],
			strBdryFileName[FILE_NAME_LEN], strBaseOutFileName[FILE_NAME_LEN],
			strLengthFileName[FILE_NAME_LEN], strExecName[FILE_NAME_LEN],
			strMetricFileName[FILE_NAME_LEN], strCanonicalName[20],
			strCurveFileName[FILE_NAME_LEN], strOptParams[256],
			strSurfInsParams[256];
	bool qReadOutFileName = false, qReadAuxBdryFileName = false,
			qLengthScaleGiven = false, qAllowBdryChanges = true,
			qAnisoCoarsening = false, qOutputBPatchFile = false,
			qWriteTempMeshes = false, qParams = false, qSurfInsParams = false;
	bool OMP_Flag = false;
	bool doGlobalSmoothing = false;
	bool qMeltIntBdryFaces = false;
	bool qIsAniso = false;
	int iQualMeasure = 2;
	int iPointInsertion = 0;
	int Queueoption = 0;
	int bdrefine = 0;
	int QuadMerge = 0;
	int keyBdryConditionNum = 1;
	double termipara = 1.;
	double dScale = 1., dGrading = 1., dMinLength = 1.e-12;
	double dMinAngle = 0.447832396928932;
	double dOffcenterAngle = 2.0 * dMinAngle;
	enum Mesh::eReorderingType type = Mesh::RCM;
	enum eEncroachType eET = eLens;
	enum {
		UNKNOWN,
		TRI,
		MESHOPT,
		COARSEN,
		HEAL,
		ANISOADAPT,
		SURFINS,
		TMOPADAPT,
		MERGE
	} iExec = UNKNOWN;
	qAdaptPred = false;
	int NumQuadLayers = 5;

	sprintf(strExecName, "%s", apcArgs[0]);
	if (strstr(strExecName, "coarsen2d")) {
		iExec = COARSEN;
		sprintf(strCanonicalName, "coarsen2d");
	} else if (strstr(strExecName, "meshopt2d")) {
		iExec = MESHOPT;
		sprintf(strCanonicalName, "meshopt2d");
	} else if (strstr(strExecName, "tri")) {
		iExec = TRI;
		sprintf(strCanonicalName, "tri");
	} else if (strstr(strExecName, "merge")) {
		iExec = MERGE;
		sprintf(strCanonicalName, "merge");
	} else if (strstr(strExecName, "heal2d")) {
		iExec = HEAL;
		sprintf(strCanonicalName, "heal2d");
	} else if (strstr(strExecName, "surfins2d")) {
		iExec = SURFINS;
		sprintf(strCanonicalName, "surfins2d");
	} else if (strstr(strExecName, "aniso_adapt2d")) {
		iExec = ANISOADAPT;
		qIsAniso = true;
		sprintf(strCanonicalName, "aniso_adapt2d");
	} else if (strstr(strExecName, "tmop_adapt2d")) {
		iExec = TMOPADAPT;
		qIsAniso = true;
		sprintf(strCanonicalName, "tmop_adapt2d");
#ifndef HAVE_MESQUITE
		vFatalError("TMOP-based adaptation requires Mesquite (version >= 2.99)\n",
				"tri.cxx:main()");
#endif
	} else {
		vFatalError(
				"Invalid executable name.  Be less creative renaming files.\n",
				"tri.cxx:main()");
	}
	GRUMMPInit(strCanonicalName);

	int iPChar, iError = 1;

	// Reading command line args.
	while ((iPChar = getopt(iNArg, apcArgs,
			"a:A:b:B:c:d:D:g:Gi:Ii:j:l:o:O:m:p:P:q:r:s:t:W:x:X:z:N:Q:R:M:k:"))
			!= EOF) {
		switch (iPChar) {
		case 'a':
			sscanf(optarg, "%lf", &dMinAngle);
			if (dMinAngle <= 0 || dMinAngle > 60.0)
				dMinAngle = 0.447832396928932;
			else
				dMinAngle *= M_PI / 180.;
			break;
		case 'A':
			qAdaptPred = (optarg[0] == '1');
			break;

		case 'b':
			qAllowBdryChanges = (optarg[0] == '1');
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

		case 'c':
			switch (optarg[0]) {
			case '0':
				eET = eBall;
				break;
			case '1':
			default:
				eET = eLens;
				break;
			case '2':
				eET = eNewLens;
				break;
			}
			break;
		case 'g':
			sscanf(optarg, "%lf", &dGrading);
			if (dGrading < 1)
				dGrading = 1;
			break;
		case 'G':
			doGlobalSmoothing = true;
			break;
		case 'I':
			qMeltIntBdryFaces = true;
			break;
		case 'i':
			strncpy(strBaseFileName, optarg, FILE_NAME_LEN - 1);
			sprintf(strMetricFileName, "%s.metric", strBaseFileName);
			strBaseFileName[FILE_NAME_LEN - 1] = 0;
			{
				int iDum;
				for (iDum = 0; iDum < FILE_NAME_LEN && strBaseFileName[iDum];
						iDum++) {
				}
				if (iDum == FILE_NAME_LEN) // No null; file name too long
					vFatalError("Input file name too long", "command line");
			}
			iError--;
			break;

		case 'j':
			qAnisoCoarsening = (optarg[0] == '1');
			break;
		case 'k':
			sscanf(optarg, "%d", &keyBdryConditionNum);
			break;

		case 'l':
			if (iExec == MESHOPT || iExec == TRI) {
				strncpy(strLengthFileName, optarg, FILE_NAME_LEN - 1);
				strLengthFileName[FILE_NAME_LEN - 1] = 0;
				{
					int iDum;
					for (iDum = 0;
							iDum < FILE_NAME_LEN && strLengthFileName[iDum];
							iDum++) {
					}
					if (iDum == FILE_NAME_LEN) // No null; file name too long
						vFatalError("Lengthscale file name too long",
								"command line");
				}
				qLengthScaleGiven = true;
			} else {
				iError += 10; // Be sure to crash out.
			}
			break;

		case 'm':
			sscanf(optarg, "%lf", &dMinLength);
			if (dMinLength <= 0)
				dMinLength = -LARGE_DBL;
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

		case 'P':
			qOutputBPatchFile = (optarg[0] == '1');
			break;

		case 'q': {
			int iTmp;
			sscanf(optarg, "%d", &iTmp);
			if (iTmp == 0 || iTmp == 1 || iTmp == 2 || iTmp == 6 || iTmp == 9
					|| iTmp == 10)
				iQualMeasure = iTmp;
			break;
		}
		case 'N':
			int NQLayer;
			sscanf(optarg, "%d", &NQLayer);
			if (NQLayer < 2) {
				vFatalError("Number of layers should be at least 2",
						"command line");
				assert(0);
			}
			NumQuadLayers = NQLayer;
			break;

		case 'r':
			sscanf(optarg, "%lf", &dScale);
			break;
		case 's':
			sscanf(optarg, "%255s", strSurfInsParams);
			qSurfInsParams = true;
			break;
		case 't':
			type = Mesh::Lexicographic;
			break;
		case 'W':
			qWriteTempMeshes = (optarg[0] == '1');
			break;
		case 'x':
			int iTmp;
			sscanf(optarg, "%d", &iTmp);
			if (iTmp == 0 || iTmp == 1 || iTmp == 3 || iTmp == 2 || iTmp == 4)
				iPointInsertion = iTmp;
			break;
		case 'X':
			if (iPointInsertion == 1) {
				sscanf(optarg, "%lf", &dOffcenterAngle);
				if (dOffcenterAngle <= 0 || dOffcenterAngle > 60.0)
					dOffcenterAngle = 0.447832396928932;
				else
					dOffcenterAngle *= M_PI / 180.;
			}
			break;
		case 'Q':
			int faceorcell;
			sscanf(optarg, "%d", &faceorcell);
			if (faceorcell == 0 || faceorcell == 1)
				Queueoption = faceorcell;
			break;
		case 'D':
			sscanf(optarg, "%lf", &termipara);
			if (termipara <= 0)
				termipara = -termipara;
			break;
		case 'R':
			sscanf(optarg, "%d", &bdrefine);
			break;
		case 'z':
			strncpy(strCurveFileName, optarg, FILE_NAME_LEN - 1);
			strCurveFileName[FILE_NAME_LEN - 1] = 0;
			{
				int iDum;
				for (iDum = 0; iDum < FILE_NAME_LEN && strCurveFileName[iDum];
						iDum++) {
				}
				if (iDum == FILE_NAME_LEN) // No null; file name too long
					vFatalError("Curve file name too long", "command line");
			}
			break;
		case 'M':
			sscanf(optarg, "%d", &QuadMerge);
			break;

		default:
			iError += 10;
			break;
		}
	}

	// Error message if the command line is too badly mangled.
	if (iError)
		vUsage(strCanonicalName);

	// File name manipulation
	makeFileName(strQualFileName, "%s.qual", strBaseFileName, "main() [tri.C]");

	if (!qReadOutFileName)
		strncpy(strBaseOutFileName, strBaseFileName, FILE_NAME_LEN - 1);

	// Open input and log files
	openMessageFile(strBaseFileName);
	logMessage(1, "Mesh quality output file: %s\n", strQualFileName);

	// small change, lets log the call to generate this function, to make repeating results possible
	for (int iArg = 0; iArg < iNArg; iArg++) {
		logMessage(1, "%s ", apcArgs[iArg]);
	}
	logMessage(1, "\n");

	///////////////////////////////////////////////////////////////////////
	// This section of the code initializes the mesh and geometry for all
	// execution paths.
	///////////////////////////////////////////////////////////////////////

	// Create a mesh object, shared by all execution paths.
	Mesh2D M2D(iQualMeasure);
	M2D.setWriteTempMesh(qWriteTempMeshes);
	if (qAllowBdryChanges)
		M2D.allowBdryChanges();
	else
		M2D.disallowBdryChanges();

	//Pointers to the geometry and mesh objects.
	GRGeom2D* pG2D = NULL;

	if (iExec == TRI) {
		pG2D = new GRGeom2D(strBaseFileName);
	}
	else if (qReadAuxBdryFileName) {
		pG2D = new GRGeom2D(strBdryFileName);
	}
	else {
		pG2D = new GRGeom2D();
	}
#ifndef NDEBUG
	pG2D->writeGeometryVTK("geometry");
#endif

	//Build objects for triangular mesh generation.
	if (iExec == TRI) {
		// FIXME: These public member variables are an abomination introduced by
		// Hongliang (?).  At the very least, they need to be private; more likely,
		// they're redundant.
		M2D.gradeparameter = dGrading;
		M2D.determineparameter = termipara;
		M2D.minangle = dMinAngle;

		// Initialize the mesh by computing the constrained Delaunay triangulation.
		M2D.setEncroachmentType(eET);
		{
			TriMeshBuilder* TMB = new TriMeshBuilder(pG2D, &M2D);
			TMB->buildConstrainedDelaunayTriangulation();
			delete TMB;
		}
		// This purge is currently redundant, because writeVTK does it anyway.
		// But purging explicitly removes that hidden dependence on behavior
		// that could change (including the act of writing the file itself).
		M2D.purgeAllEntities();
		writeVTK(M2D, "constrained.vtk");
	}

	// Read objects for other executables
	else {
		// Reading an existing mesh file

		M2D.readFromFile(strBaseFileName);
		writeVTK(M2D, "initialMesh");
		//    writeTriangle(M2D,"initialMesh");
		M2D.relateMeshToGeometry(pG2D);

		logMessage(1, "Number of cells read: %u\n", M2D.getNumCells());

		M2D.evaluateQuality(); // 0th Quality Eval
	}

	assert(M2D.isValid());

	///////////////////////////////////////////////////////////////////////
	// Set up a length scale, either scalar or metric.
	// Also, for isotropic meshes, set up a refinement manager, since
	// they're all going to need one.
	///////////////////////////////////////////////////////////////////////

	std::shared_ptr<GRUMMP::Length> pSizingField;
	if (qIsAniso) {

	}
	else {
		//Set up the mesh with an automatic sizing field.
		 pSizingField = GRUMMP::Length::Create(&M2D, dScale, dGrading,
				 dMinLength);
		if (qLengthScaleGiven) {
			// Read a length-scale file and update the sizing field.

			FILE *pFLenFile = fopen(strLengthFileName, "r");
			if (NULL == pFLenFile)
				vFatalError("Couldn't open file for reading",
						"length scale input");

			double x_loc = -LARGE_DBL, y_loc = -LARGE_DBL, length = -LARGE_DBL;

			while (1) {
				int iDataRead = fscanf(pFLenFile, "%lf%lf%lf", &(x_loc),
						&(y_loc), &(length));
				if (iDataRead != 3)
					break;
				pSizingField->addToLengthMesh(length, x_loc, y_loc);
			}
			fclose(pFLenFile);
		}
		pSizingField->setLengthScale(&M2D);

		double minFaceLength = LARGE_DBL;

		for (GR_index_t ii = 0; ii < M2D.getNumFaces(); ii++)
			minFaceLength = min(M2D.getFace(ii)->calcSize(), minFaceLength);
		M2D.checkEdgeLengths();
		printf("Minimum face length %.16f\n", minFaceLength);
	}


	int iNumPasses = 3;
	switch (iExec) {
	case ANISOADAPT: { // Doug's adaptation driver, circa 2008.
		// TODO Mahkame Is there any user documentation for aniso adapt?
		const double alpha = 0.0;
		const double angleTolerance = 28 * M_PI / 180.0;
		const double refRatio = 1.5; //1.3

		// TODO Mahkame A hard restriction on max mesh size is a very bad idea.
#define MAXVERTS 20000

		GRUMMP::SwapDecider2D *pSD_Degree = new GRUMMP::AnisoSwapDecider2D();
		// GRUMMP::SwapDecider2D *pSD_Degree = new GRUMMP::UniformDegreeSwapDecider2D();

		AnisoRefinement ARMesh2D(&M2D, strMetricFileName, alpha);
//	AnisoRefinement ARMesh2D(&M2D, 0);
//	ARMesh2D.vAssignMetricsToVertsTest(
//	    AnisoRefinement::e2DTestQuarterCircle);

		//    AnisoRefinement ARMesh2D(&M2D,pvfMetric,2,alpha);
		M2D.attachAnisoRef(&ARMesh2D);
//		ARMesh2D.vHealBdry();
		// if (iNumPasses > 0) {
		for (int iPass = 0; iPass < iNumPasses; iPass++) {
			// Do some smoothing if requested
#ifdef HAVE_MESQUITE
			//			logMessage(1, "Smoothing vertices anisotropically...\n");
			GRUMMP::SmoothingManagerMesquite2D *pSMM2D =
					new GRUMMP::SmoothingManagerMesquite2D(&M2D);
			GRUMMP::SwapManager2D SM_Degree(pSD_Degree, &M2D);
			pSMM2D->vSetTargetfromMetric();
			pSMM2D->smoothAllVerts(1);
			writeVTK(M2D, "AfterSmoothing");
			//	      writeNative(M2D,"AfterSmoothing");
			// TODO Mahkame Another hard-coded file name.
			//			FILE *pFOutFile2;
			//			pFOutFile2=fopen("/home/mahkame/GRUMMP/mahkame/src/programs/RefinedMetricSmoothed","w");
			//			fprintf(pFOutFile2, "2");
			//			for (unsigned int iV = 0; iV < M2D.iNumVerts(); iV++)
			//			{
			//				Vert *pV = M2D.getVert (iV);
			//				fprintf(pFOutFile2,"%16.15f\t%16.15f\t%16.15f\t%16.15f\t%16.15f\n",pV->x(),pV->y(),pV->getMetric(0),pV->getMetric(1),pV->getMetric(2));
			//			}
			//			fclose(pFOutFile2);
#endif

			for (unsigned int iF = 0; iF < M2D.getNumFaces(); iF++) {
				Face *pF = M2D.getFace(iF);
				M2D.pAnisoRef->iAniSwap(pF);
			}
			//		     SM_Degree.swapAllFaces();
			writeVTK(M2D, "AfterSwapping");
			//    writeNative(M2D,"AfterSwapping");
			double dRho = refRatio;
			double dThetaMin = angleTolerance;
			double dAdes = ARMesh2D.dAvgArea() / dRho;
			double dIdeal = 2 * sqrt(dAdes / sqrt(3));
			double dTau = dIdeal
					* (7.0 - 4.0 * sqrt(3.0) * sin(dThetaMin)
							- 4.0 * cos(dThetaMin) * cos(dThetaMin))
					/ (4.0 * cos(dThetaMin) * cos(dThetaMin) - 1.0);
			double dMaxCircRad = 1.0 / sqrt(3) * (dIdeal + (dTau));
			double MinLength = dIdeal - dTau;

			assert(dMaxCircRad > 0);
			assert(MinLength > 0);
			CellQueue* pCellQueue = new CellQueue(dynamic_cast<Mesh*>(&M2D),
					dMaxCircRad, MinLength, true);
			int iMaxVerts = MAXVERTS;
			ARMesh2D.vImproveMesh(pCellQueue, iMaxVerts, 0.4, false); //0.2 works well
			//			FILE *pFOutFile;
			//			pFOutFile=fopen("/home/mahkame/GRUMMP/mahkame/src/programs/RefinedMetric","w");
			//			fprintf(pFOutFile, "2");
			//			for (unsigned int iV = 0; iV < M2D.getNumVerts(); iV++)
			//			{
			//				Vert *pV = M2D.getVert (iV);
			//				fprintf(pFOutFile,"%16.15f\t%16.15f\t%16.15f\t%16.15f\t%16.15f\n",pV->x(),pV->y(),pV->getMetric(0),pV->getMetric(1),pV->getMetric(2));
			//			}
			//			fclose(pFOutFile);
			//
			//			std::cout << "Wrote \"smoothed_mesh.vtk\"" << std::endl;

			//M2D.vEvaluateAnisoQuality();
		}
		writeVTK(M2D, "smoothed_mesh");
		writeNative(M2D, "smoothed_mesh");

		delete pSD_Degree;
		break;
	}
	case TMOPADAPT: { // David's adaptation driver, circa 2015.

#ifdef HAVE_MESQUITE
		std::vector<Face*> faces_to_lock;
		if (qMeltIntBdryFaces) {
			M2D.meltIntBdryFaces(faces_to_lock);
			M2D.lock_faces(faces_to_lock);
		}

		double dMinMLength;
		double dMaxMLength;
		double alpha = 0;

		AnisoRefinement *pAR = new AnisoRefinement(&M2D, alpha);

		/*
//		 * Set a metric for every point in the mesh
//		 */
//		FILE *pFMetricFile = fopen(strMetricFileName, "r");
//		pAR->vAssignMetricsToVerts(pFMetricFile);
//		/*
//		 * Computing max and min metric length allowed, as described by Pagnutti "Anisotropic Adaptation: Metric and Meshes"
//		 */
//		const double angleTolerance = 28 * M_PI / 180.0;
//		const double refRatio = 1.4;
//		double dRho = refRatio;
//		double dThetaMin = angleTolerance;
//		double dAdes = pAR->dAvgArea() / dRho;
//		double dIdeal = 2 * sqrt(dAdes / sqrt(3));
//		double dTau = dIdeal * (7.0 - 4.0 * sqrt(3.0) * sin(dThetaMin)- 4.0 * cos(dThetaMin) * cos(dThetaMin))/ (4.0 * cos(dThetaMin) * cos(dThetaMin) - 1.0);
//		dMinMLength = (dIdeal - dTau);
//		dMaxMLength = (dIdeal + dTau);
		/*
		 * Used for setting analytic metric functions. Functions must be set in a unit square mesh
		 */
		pAR->vAssignMetricsToVertsTest(pAR->e2DTestCross);    //vAssignMetricsToVertsTest is for test cases such as e2DTestQuarterCircle and e2DTestCross in case of 2D
		dMinMLength = sqrt(2) / 2.;
		dMaxMLength = sqrt(2);
		GRUMMP::QualMeasure *pQM = new GRUMMP::Angle2D;
		GRUMMP::SwapDecider2D* pSD2D = new GRUMMP::MetricSwapDecider2D(pQM);
		GRUMMP::InsertionPointCalculator * pIPC = new GRUMMP::AnisoEdgeIPC2D(
				dMinMLength, dMaxMLength);
		GRUMMP::SmoothingManagerMesquite2D *pSMM2D =
				new GRUMMP::SmoothingManagerMesquite2D(&M2D);
		GRUMMP::AnisoEdgeRefinementManager2D ARM2D(&M2D, dMaxMLength,
				dMinMLength, doGlobalSmoothing, pIPC, pSMM2D, pAR, pSD2D);

		ARM2D.AnisoRefineMesh();

//	sprintf(strAnisoFileName, "Aniso_Edge_Adapted_Mesh");
//
//	writeVTKLegacy(M2D, strAnisoFileName);
//	writeNativeMESH(M2D, strAnisoFileName);
//	writeNativeGRMESH(M2D, strAnisoFileName);
//
//	sprintf(strAnisoFileName, "Aniso_Edge_Adapted_Mesh_Qual");
//	writeVTK(M2D, strAnisoFileName, "", Mesh::TMOPQual);

		delete pAR;
		delete pSMM2D;
		delete pIPC;
		delete pSD2D;
		delete pQM;

//////////////////////////////////////////////////////////////
		/*
		 * Merge tris into quads according to quality
		 * Needs a %s.metric file
		 * Qual = qual range for merging 2 tri
		 * bMaxQual = if true, merges a quad if both tri are best quality for merging
		 */
//
//		GRUMMP::TMOPQual pQM2D;
//		double dQual = 1e6;
//		bool bMaxQual = 0;
//		M2D.allowNonSimplicial();
//		M2D.allowBdryChanges();
//		pQM2D.MergeTriMeshToQuadMesh(bMaxQual,&M2D,&dQual);
////////////////////////////////////////////////////////////////
#endif
		break;
	}

	case MESHOPT:
	case COARSEN:
	case TRI:
	{
		// These three cases are such close variations that it makes sense to combine them.
		M2D.evaluateQuality();

		// Set up the refinement manager here; all three executables will use it.
		GRUMMP::InsertionPointCalculator *pIPC = nullptr;

		if (iPointInsertion == 1) {
			pIPC = new GRUMMP::OffcenterIPC2D(dOffcenterAngle);
		} else if (iPointInsertion == 3) {
			pIPC = new GRUMMP::OptimalSteinerIPC2D(dMinAngle);
		} else if (iPointInsertion == 2) {
			pIPC = new GRUMMP::OffcenterEngwirda2D(dGrading, dOffcenterAngle);
		} else if (iPointInsertion == 4) {
			pIPC = new GRUMMP::AF2D(dGrading, dOffcenterAngle);
		} else
			pIPC = new GRUMMP::CircumcenterIPC();

		writeVTK(M2D, "initial-length.vtk");
		logMessage(1, "Refining with point insertion using %s\n",
				pIPC->getName().c_str());
		GRUMMP::RefinementManager2D RM2D(&M2D, pIPC, pSizingField, NULL);
		RM2D.queuetype(Queueoption);
		RM2D.setDifferentQueue(iPointInsertion);
		RM2D.setQualityTarget(2.0 * sin(dMinAngle) / sqrt(3));
		RM2D.setRefineForSize(true);
		RM2D.setKeyBdryCondition(keyBdryConditionNum);

		// For tri, need to refine the bdry
		if (iExec == TRI) {
			if (bdrefine == 1) {
				std::cout << "using new bdry refinement" << std::endl;
				double dBGrading;
				double dBScale;
				std::cout
						<< "Enter Bdry Refinement control parameter, scale and grading"
						<< std::endl;
				std::cin >> dBScale >> dBGrading;
				std::cout << "Bdry Scale: " << dBScale << "\tBdry Grading: "
						<< dBGrading << std::endl;
				std::shared_ptr<GRUMMP::Length> len2d = GRUMMP::Length::Create(&M2D,
						dBScale, dBGrading);
				GRUMMP::SurfaceInserter2D SurfIns2D(&M2D, len2d, dBGrading, dBScale,
						strSurfInsParams);
	//			GRUMMP::SurfaceInserter2D SurfIns2D(&M2D, dGrading, dScale, strSurfInsParams);
				SurfIns2D.refineBoundaryForLength(pG2D);
				RM2D.Set_Refine2D_Timer();
			} else {
				M2D.writeTempMesh();
				RM2D.refineBoundaryForLength();
				RM2D.Set_Refine2D_Timer();
			}
			timecheck = clock();
			writeVTK(M2D, "prerefine.vtk");
		}

		// For the other two, coarsen appropriately
		else if (iExec == MESHOPT) {
			// For meshopt, eliminate edges that are too short according to the
			// sizing implied by the mesh.  Shouldn't do very much, if anything.
			GRUMMP::CoarsenManager CM2D(&M2D, pSizingField);
			CM2D.coarsen(1);
			pSizingField->scale(1.25); // This factor of 1.25 is a heuristic that
									   // preserves mesh size (approximately) when
									   // running meshopt2d.
			pSizingField->setLengthScale(&M2D);
		}
		else {
			// For coarsen, this bit does the main job of coarsening the mesh.
			GRUMMP::CoarsenManager CM2D(&M2D, pSizingField);
			if (qAnisoCoarsening) {
				CM2D.setInteriorSkip(3);
				CM2D.setFillBoundary(true);
			} else {
				CM2D.setInteriorSkip(1);
				CM2D.setFillBoundary(false);
			}

			CM2D.coarsen(2);
			pSizingField->scale(1.25); // This is the same factor of 1.25 as above.
			pSizingField->setLengthScale(&M2D);
			writeVTK(M2D, "after-coarsen.vtk");

			logMessage(1, "Number of cells after coarsening: %d\n",
					M2D.getNumCells());
		}
		bool qResult = M2D.makeDelaunay();
		if (!qResult) {
			vWarning(
					"Couldn't make mesh fully Delaunay before refining!\n");
			vWarning("This MIGHT still turn out okay...\n");
		}

		if (iExec == COARSEN && qAnisoCoarsening) {
			// Iso refinement for quality fixes is a terrible idea, so dont.
		}
		else {
			// For all three, refine to match length scale and quality measure
			if (Queueoption == 1) { //face based queue
				RM2D.refineMeshonfaces(INT_MAX, OMP_Flag, NumQuadLayers);
			}
			if (Queueoption == 0) { //cell based queue
				RM2D.refineMesh(INT_MAX, OMP_Flag, NumQuadLayers);
			}
			RM2D.Print_Refine2D_Timer(1);
		}
		RM2D.detachFromObservable();
		delete pIPC;

		logMessage(1, "Number of vertices after refinement:  %u\n",
				M2D.getNumVerts());
		M2D.purgeAllEntities();
		M2D.evaluateQuality();

		writeVTK(M2D, "after-refine.vtk");

		if (QuadMerge == 1) {
			/////need to be fixed!!!////
			M2D.allowNonSimplicial();
			M2D.TriMergetoQuad();
			writeVTK(M2D, "quadmesh.vtk");
		}

		//////////////////////////////////////////////////////////////////
		// Post-process the mesh with a prescribed recipe of swapping and
		// smoothing.
		//////////////////////////////////////////////////////////////////

		break;
	}
	case MERGE: { // Merge triangles to quads

		M2D.allowNonSimplicial();
		M2D.TriMergetoQuad2();
		writeVTK(M2D, "Quadmesh2.vtk");
		M2D.purgeAllEntities();
		//    M2D.vReorder();
		writeNative(M2D, strBaseOutFileName, qOutputBPatchFile, "");
		writeVTK(M2D, strBaseOutFileName);
		iNumPasses = 0;

	}
		break;
	case HEAL: // More or less isotropic adaptation, all commented out.
		for (GR_index_t iV = 0; iV < M2D.getNumVerts(); iV++) {
			printf("Before - Vert Number %d X: %f Y %f\n", iV,
					M2D.getVert(iV)->x(), M2D.getVert(iV)->y());
		}
		M2D.reorderAllEntities(type);
		//     M2D.vSetLengthScale(7.5);

		// //     {
		// //       // Total hack job for simulating mesh adaptation.
		// //       for (int ii = 0; ii < M2D.iNumVerts(); ii++) {
		// //	Vert *pV = M2D.pVVert(ii);
		// //	double dX = pV->dX();
		// //	double dY = pV->dY();
		// //	double dX_star = 250 - (dY-541)*(dY-541) / 130;
		// //	double dLen;
		// //	if (dX < dX_star) {
		// //	  dLen = 1. / (0.125 + 10*exp(-0.5*(dX_star - dX)));
		// //	}
		// //	else {
		// //	  dLen = 1. / (0.05 + 10*exp(-(dX-dX_star)/3.));
		// //	}
		// //	pV->vSetLS(dLen);
		// //       }
		// //     }

		//     vMessage(1, "Read mesh, preparing to heal it\n");
		//     M2D.vCheckEdgeLengths();

		//     M2D.vAdaptToLengthScale(eET);
		//		pG2D->OutputCurveVTK("currentGeom");

		//     M2D.vEvaluateQuality();
		//     M2D.vCheckEdgeLengths();
		//     // vWriteFile_Mesh2D(M2D, strBaseOutFileName, ".insert");

		break;
	case SURFINS: { // Insert a curve from a file.
		M2D.setEncroachmentType(eET);
		// Moved down here due to errors with receiving events..

		logMessage(1, "\n\nComputing Surface Insertion \n \n");

		if (!qParams) {
			sscanf("wirmf5wimf5mf5", "%255s", strOptParams);
		}
		if (!qSurfInsParams) {
			sscanf("dd", "%255s", strSurfInsParams);
		}
		iNumPasses = 0;

		std::shared_ptr<GRUMMP::Length> len2d = GRUMMP::Length::Create(&M2D,
				dScale, dGrading);
		GRUMMP::SurfaceInserter2D *pSurfIns2D = new GRUMMP::SurfaceInserter2D(
				&M2D, len2d, dScale, dGrading, strSurfInsParams);
		pSurfIns2D->insertCurveInMesh(strCurveFileName, strOptParams);
		assert(pSurfIns2D != NULL);
		delete pSurfIns2D;

		if (qMeltIntBdryFaces) {
			std::vector<Face*> faces_to_lock;
			M2D.meltIntBdryFaces(faces_to_lock);
		}
//		GRUMMP::SwapDecider2D *pSD_Delaunay = new GRUMMP::DelaunaySwapDecider2D();
//		GRUMMP::SwapManager2D *pSM = new GRUMMP::SwapManager2D(pSD_Delaunay, &M2D);
//
//		pSM->swapAllFaces(false);
//		delete pSM;
//		delete pSD_Delaunay;
	}

		break;
	default:
		assert(0);
		break;
	}

	// moved down here to prevent sendEvents error (D. Zaide Nov 2013)
	GRUMMP::OptMSSmoothingManager2D *smoothOptMS =
			new GRUMMP::OptMSSmoothingManager2D(&M2D);
	GRUMMP::OptMSSmoothingManager2D *smoothOptMSSecond =
			new GRUMMP::OptMSSmoothingManager2D(&M2D);

	if (QuadMerge == 1 or iExec == MERGE) { // Set smoothing for quads
		smoothOptMS->setFunctionID(11);		//MAX_MIN_Interior_Sin
//			smoothOptMS->setFunctionID(15);
		smoothOptMSSecond->setFunctionID(16);
	}

	if (iExec != SURFINS && iExec != ANISOADAPT && iExec != TMOPADAPT) { // quality improvement
		/////////////////
			GRUMMP::SwapDecider2D *pSD;
			if (Queueoption == 1 or iPointInsertion == 2) {
				//AFLR or Engwirda's frontal-Delaunay
				pSD = new GRUMMP::UniformDegreeSwapDecider2D();
			} else {
				pSD = new GRUMMP::DelaunaySwapDecider2D();
			}
//		GRUMMP::SwapDecider2D *pSD = new GRUMMP::MinMaxAngleSwapDecider2D();
			GRUMMP::SwapManager2D *pSM = new GRUMMP::SwapManager2D(pSD,
					&M2D);
			pSM->Set_Swap2D_Timer();
			pSM->swapAllFaces(OMP_Flag);
			if (iNumPasses > 0) {
				// Do some smoothing if requested
				logMessage(1, "Smoothing vertices...\n");
				if (Queueoption == 1 or iPointInsertion == 2) {
					//AFLR or Engwirda's frontal-Delaunay
					smoothOptMS->setSmoothingTechnique('o');
				} else {
					smoothOptMS->setSmoothingThreshold(30);
				}

				for (int ii = 0; ii < iNumPasses; ii++) {
					logMessage(2, "Pass %d of %d...\n", ii + 1, iNumPasses);
					smoothOptMS->smoothAllVerts(1);
					// TODO  Can we actually smooth a mesh with quads?
					// Not sure.  Depends on both smoothing qual funcs
					// and on data management in OptMS.
//				smoothOptMSSecond->smoothAllVerts(1);
					// 	iMesh_Plane *iMP = new iMesh_Plane(&M2D);
					// if Mesquite smoothing {
					//   Create a Mesquite::MsqIMesh from iMP
					//   Set up a smoothing context, etc, and smooth in Mesquite.
					// }
					if ((ii + 1) % 3 == 0) {
						pSM->swapAllFaces(OMP_Flag);
					}
				}

				M2D.evaluateQuality();
			}
			pSM->Print_Swap2D_Timer(1);
			delete pSM;
			delete pSD;
		}

	delete smoothOptMS;
	delete smoothOptMSSecond;
//

	// Mesh output
	logMessage(1, "Writing mesh to file(s)...\n");
	// Delete any extraneous verts
	for (unsigned ii = 0; ii < M2D.getNumVerts(); ii++) {
		Vert *pV = M2D.getVert(ii);
		if (pV->getNumFaces() == 0) {
			M2D.deleteVert(pV);
		}
	}

	if (iExec == TRI) { // Re-ordering and output for mesh gen cases
		if (QuadMerge == 0) {
			M2D.reorderAllEntities(type);
		}
		assert(M2D.isValid());
		writeVTK(M2D, strBaseOutFileName);
		writeNative(M2D, strBaseOutFileName, qOutputBPatchFile);
		// writeFEA(M2D, strBaseOutFileName, 2);
		{
			int intDegree[10], bdryDegree[10];
			for (int ii = 0; ii < 10; ii++) {
				intDegree[ii] = bdryDegree[ii] = 0;
			}
			for (GR_index_t ii = 0; ii < M2D.getNumVerts(); ii++) {
				Vert *pV = M2D.getVert(ii);
				int deg = pV->getNumFaces();
				if (pV->isBdryVert()) {
					bdryDegree[deg]++;
				} else {
					intDegree[deg]++;
				}
			}
			for (int ii = 2; ii < 10; ii++) {
				logMessage(2, "  Degree: %d  Int count: %5d  Bdry count: %5d\n",
						ii, intDegree[ii], bdryDegree[ii]);
			}
		}
	} else { // Output without re-ordering because of quads
		// For now, don't do a full re-ordering, because that doesn't
		// handle quads gracefully.
		M2D.purgeAllEntities();
		//    M2D.vReorder();

		writeNative(M2D, strBaseOutFileName, qOutputBPatchFile, "");
		writeVTK(M2D, strBaseOutFileName);

		//     writeFEA(M2D, strBaseOutFileName, 5, ".out");
	}
	M2D.writeQualityToFile(strQualFileName);
	logMessage(1, "Done\n");
	double minFaceLength = LARGE_DBL;
	for (GR_index_t ii = 0; ii < M2D.getNumFaces(); ii++)
		minFaceLength = min(M2D.getFace(ii)->calcSize(), minFaceLength);
	//   M2D.vCheckEdgeLengths();
	printf("Minimum face length %.16f\n", minFaceLength);

	assert(pG2D != NULL);
	delete pG2D;

	timecheck = clock() - timecheck;
	std::cout << "CPU time: " << timecheck / (double) CLOCKS_PER_SEC
			<< std::endl;
	return 0;

}
