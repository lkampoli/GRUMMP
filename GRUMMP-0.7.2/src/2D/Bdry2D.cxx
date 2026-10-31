#include <ctype.h>
#include <string.h>
#include "GR_Bdry2D.h"
#include "GR_BdryPatch2D.h"
#include "GR_List.h"

// Bdry2D::Bdry2D(const Bdry2D& B2DIn)
//   : BdryRep(), iSz(0), LpBPatches(), LdPoints(), iSizePoints(0), iCount(0),
//     iNumPointsUsed(0), LqBdryPoint(), iNum2DPatchesUsed(0),
//     LsVertConnectivity(), iNumVertsConnectFound(0), iMaxVertsConnectFound(0),
//     LiNumPatchesToPoint(), LiStartVert(), LqVertSmallAngle(), LqVertUsed()
// {
// //  vMessage(2, "In Bdry2D copy constructor\n");

//   // Make sure pointers are null to start out, just as in the default
//   // constructor...
//   LpBPatches.vClear();
//   LqVertUsed.vClear();
// //    LsPatchUsed.vClear();
//   iNPatches = iMaxNPatches = 0; // Test
//   iNumVertsConnectFound = 0;

//   // Copy the Bdry2D info in the internal structure...
//   int i;
//   int iNumSegs = 0;
//   BdryPatch2D *pPatch;
//   structVertConnect sVConnect;

//   vMessage (4, "Copying the Bdry2D information\n");

//   // Initialize size..
//   //    B2D.vSetSize(B2DInfo.iNumBdryPatches(), B2DInfo.iNumVerts());
//   // Should set the size of the num of patches to the num of patches used!
// //  vMessage (3, "Initializing the new Bdry2D class with %d patches and %d verts\n",
// //	    B2DInfo.iNumBdryPatchesUsed (), B2DInfo.iNumPoints ());
//   vSetSize (B2DIn.iNumPatchesUsed (), B2DIn.iNumPoints ());

//   double adC[2];

//   // Add each vertex...
//   for (i = 0; i < B2DIn.iNumPoints (); i++) {
// //    vMessage(2, "   Adding vertex # %d ...", i);
//     B2DIn.vGetPoint(i, adC);
//     vAddPoint(adC);
//     vSetBdryPoint (i, B2DIn.qIsBdryPoint (i));
//     vSetNumPatchesConnectedToPoint(i, B2DIn.iNumPatchesConnectedToPoint(i));
// //    vMessage(2, " (%f, %f) - Used: %d\n", adC[0], adC[1], B2D.qIsBdryPoint (i));
//   }

//   // Copy vert connectivity stuff too
//   for (i = 0; i < B2DIn.iConnectivityListLength(); i++) {
//     sVConnect = B2DIn.sVCConnectivityListContent(i);
//     structVertConnect sVC2;
//     //		sVC2 = new structVertConnect;
//     sVC2.iVert  = sVConnect.iVert;
//     sVC2.iPatch = sVConnect.iPatch;
//     sVC2.dAngle = sVConnect.dAngle;
//     vAddToConnectivityList(sVC2);
//   }

//   // Copy small angles info
//   for (i = 0; i < B2DIn.iNumPoints(); i++) {
//     vSetSmallAngle(i, B2DIn.qVertSmallAngle(i));
//     int iStart = B2DIn.iStartVert(i);
//     vSetStartVert(i, iStart);
//   }

//   // Then add each patch...
//   for (i = 0; i < B2DIn.iNumPatches (); i++) {
//     // Must make sure we only copy USED patches...
//     pPatch = B2DIn.pBdry (i);
//       // Should really check the type of patch first..!!!!!
//     switch (pPatch->eGeom2DType) {
//     case BdryPatch2D::ePolyLine:
//       {
// 	vMessage(4, "   Creating BdrySeg #%d (was %d before) ...", iNumSegs, i);

// 	BdryPatch2D *pBChk = new BdrySeg (this, pPatch->iLeftBdryCond (),
// 					  pPatch->iRightBdryCond (),
// 					  pPatch->iLeftRegion (),
// 					  pPatch->iRightRegion (),
// 					  pPatch->iBeginPointIndex (),
// 					  pPatch->iEndPointIndex ()
// 					  );
// 	vAddPatchToList (pBChk);
// 	iNumSegs++;
// 	/*
// 	  B2D.vAddPatchToList(BdryPatch2D::ePolyLine,
// 	  B2DInfo.pBdry(i)->iLeftBdryCond(),
// 	  B2DInfo.pBdry(i)->iRightBdryCond(),
// 	  B2DInfo.pBdry(i)->iLeftRegion(),
// 	  B2DInfo.pBdry(i)->iRightRegion(),
// 	  B2DInfo.iBeginVert(i),
// 	  B2DInfo.iEndVert(i)
// 	  );
// 	*/
// 	vMessage(4, " and added\n");
//       }
//       break;
//     case BdryPatch2D::eCircle:
//       {
// 				// For circles...
// 	vMessage(4, "   Creating circle BdryArc #%d (was %d before) ...", iNumSegs, i);

// 	dynamic_cast<BdryArc*>(pPatch)->vGetCenter(adC);

// 	BdryPatch2D *pBChk = new BdryArc
// 	  (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	   pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	   dynamic_cast<BdryArc*>(pPatch)->dGetRadius(),
// 	   dynamic_cast<BdryArc*>(pPatch)->iCenterPointIndex());
// 	vAddPatchToList (pBChk);
// 	iNumSegs++;

// 	vMessage(4, " and added\n");
//       }
//       break;
//     case BdryPatch2D::eCircularArc:
//       {
// 				// For arcs....
// 	vMessage(4, "   Creating arc BdryArc #%d (was %d before) ...", iNumSegs, i);

// 	dynamic_cast<BdryArc*>(pPatch)->vGetCenter(adC);

// 	BdryPatch2D *pBChk = new BdryArc
// 	  (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	   pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	   pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	   dynamic_cast<BdryArc*>(pPatch)->dGetRadius());
// 	vAddPatchToList (pBChk);
// 	iNumSegs++;

// 	vMessage(4, " and added\n");
//       }
//       break;
//     case BdryPatch2D::eLongArc:
//       {
// 				// For arcs....
// 	vMessage(4, "   Creating long arc BdryArc #%d (was %d before) ...", iNumSegs, i);

// 	dynamic_cast<BdryArc*>(pPatch)->vGetCenter(adC);

// 	BdryPatch2D *pBChk = new BdryArc
// 	  (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	   pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	   pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	   dynamic_cast<BdryArc*>(pPatch)->dGetRadius(), true);
// 	vAddPatchToList (pBChk);
// 	iNumSegs++;

// 	vMessage(4, " and added\n");
//       }
//       break;
//     case BdryPatch2D::eCubicParam:
//       {
// 				// Cubic parametric curves
// 	vMessage(4, "   Creating cubic param curve #%d (was %d before) ...", iNumSegs, i);

// 	double dax, dbx, dcx, ddx, day, dby, dcy, ddy;
// 				// Get the parameters
// 	dynamic_cast<BdryCubicParam2D*>(pPatch)
// 	  ->vGetAllCoefficients(&dax, &dbx, &dcx, &ddx,
// 				&day, &dby, &dcy, &ddy);
// 				// Create the new patch
// 	BdryPatch2D *pBChk = new BdryCubicParam2D
// 	  (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	   pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	   pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	   dax, dbx, dcx, ddx, day, dby, dcy, ddy);

// 	vAddPatchToList(pBChk);
// 	iNumSegs++;

// 	vMessage(4, " and added\n");
//       }
//       break;

//     case BdryPatch2D::eBezier:
//       {
// 	vMessage(4, "   Creating bezier curve #%d (was %d before) ...", iNumSegs, i);

// 				// Create the new patch
// 	BdryPatch2D *pBChk = new BdryBezier2D
// 	  (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	   pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	   pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	   dynamic_cast<BdryBezier2D*>(pPatch)->iControlPointIndex(1),
// 	   dynamic_cast<BdryBezier2D*>(pPatch)->iControlPointIndex(2));

// 	vAddPatchToList(pBChk);
// 	iNumSegs++;

// 	vMessage(4, " and added\n");

//       }

//       break;
//     case BdryPatch2D::eInterpolated:
//       {
// 	vMessage(4, "   Creating interpolated spline #%d (was %d before) ...", iNumSegs, i);

// 				// Get all the points in the list
// 	List<int> LiPts;
// 	for (int ii = 0;
// 	     ii < dynamic_cast<BdryInterpolated2D*>(pPatch)->iNumPoints();
// 	     ii++) {
// 	  int iDum = dynamic_cast<BdryInterpolated2D*>(pPatch)->iGetPoint(ii);
// 	  LiPts.vAppendItem( iDum );
// 	}

// 				// Create the new patch
// 	BdryPatch2D *pBChk =
// 	  new BdryInterpolated2D(this, pPatch->iLeftBdryCond (),
// 				 pPatch->iRightBdryCond (),
// 				 pPatch->iLeftRegion (),
// 				 pPatch->iRightRegion (),
// 				 LiPts);

// 	vAddPatchToList(pBChk);
// 	iNumSegs++;

// 	// Now must set the coefficients manually in case they were changed
// 	// because of a split or something
// 	double d1, d2, d3, d4, d5, d6, d7, d8;

// 	assert( dynamic_cast<BdryInterpolated2D*>(pPatch)->iNumSplines() ==
// 		dynamic_cast<BdryInterpolated2D*>(pBChk)->iNumSplines());

// 	for (int iCP = 0;
// 	     iCP < dynamic_cast<BdryInterpolated2D*>(pBChk)->iNumSplines();
// 	     iCP++) {
// 	  dynamic_cast<BdryInterpolated2D*>(pPatch)->
// 	    vGetSplineCoefficients(iCP, &d1, &d2, &d3, &d4,
// 				   &d5, &d6, &d7, &d8);
// 	  dynamic_cast<BdryInterpolated2D*>(pBChk)->
// 	    vSetSplineCoefficients(iCP, d1, d2, d3, d4, d5, d6, d7, d8);
// 	}

// 	dynamic_cast<BdryInterpolated2D*>(pBChk)->vRefresh();

// 	vMessage(4, " and added\n");

//       }

//       break;
//     case BdryPatch2D::eUnknown:
//       assert(0);
//       break;
//     } // end switch..
//   }

//   vMessage(4, "Finished copying the Bdry2D info\n");
// }

// Bdry2D& Bdry2D::operator=(const Bdry2D & B2DIn)
// {
// //  vMessage(2, "In Bdry2D operator=\n");

//   if (this != &B2DIn) {
//     // Copy the Bdry2D info in the internal structure...
//     int i;
//     int iNumSegs = 0;
//     BdryPatch2D *pPatch;
//     structVertConnect sVConnect;
//     double adC[2];

//     vMessage (4, "Copying the Bdry2D information\n");

//     // Initialize size..
//     //    B2D.vSetSize(B2DIn.iNumBdryPatches(), B2DIn.iNumVerts());
//     // Should set the size of the num of patches to the num of patches used!
//     //  vMessage (3, "Initializing the new Bbdry2D class with %d patches and %d verts\n",
//     //	    B2DIn.iNumBdryPatchesUsed (), B2DIn.iNumPoints ());
//     vSetSize (B2DIn.iNumPatchesUsed (), B2DIn.iNumPoints ());

//     // Add each vertex...
//     for (i = 0; i < B2DIn.iNumPoints (); i++) {
//       //    vMessage(2, "   Adding vertex # %d ...", i);
//       B2DIn.vGetPoint(i, adC);
//       vAddPoint(adC);
//       vSetBdryPoint (i, B2DIn.qIsBdryPoint (i));
//       vSetNumPatchesConnectedToPoint(i, B2DIn.iNumPatchesConnectedToPoint(i));
//       //    vMessage(2, " (%f, %f) - Used: %d\n", adC[0], adC[1], B2D.qIsBdryPoint (i));
//     }

//     // Copy vert connectivity stuff too
//     for (i = 0; i < B2DIn.iConnectivityListLength(); i++) {
//       sVConnect = B2DIn.sVCConnectivityListContent(i);
//       structVertConnect sVC2;
//       //		sVC2 = new structVertConnect;
//       sVC2.iVert  = sVConnect.iVert;
//       sVC2.iPatch = sVConnect.iPatch;
//       sVC2.dAngle = sVConnect.dAngle;
//       vAddToConnectivityList(sVC2);
//     }

//     // Copy small angles info
//     for (i = 0; i < B2DIn.iNumPoints(); i++) {
//       vSetSmallAngle(i, B2DIn.qVertSmallAngle(i));
//       int iStart = B2DIn.iStartVert(i);
//       vSetStartVert(i, iStart);
//     }

//     // Then add each patch...
//     for (i = 0; i < B2DIn.iNumPatches (); i++) {
//       // Must make sure we only copy USED patches...
//       pPatch = B2DIn.pBdry (i);
//       // Should really check the type of patch first..!!!!!
//       switch (pPatch->eGeom2DType) {
//       case BdryPatch2D::ePolyLine:
// 	{
// 	  vMessage(4, "   Creating BdrySeg #%d (was %d before) ...", iNumSegs, i);

// 	  BdryPatch2D *pBChk = new BdrySeg (this, pPatch->iLeftBdryCond (),
// 					    pPatch->iRightBdryCond (),
// 					    pPatch->iLeftRegion (),
// 					    pPatch->iRightRegion (),
// 					    pPatch->iBeginPointIndex (),
// 					    pPatch->iEndPointIndex ()
// 					    );
// 	  vAddPatchToList (pBChk);
// 	  iNumSegs++;
// 	  vMessage(4, " and added\n");
// 	}
// 	break;
//       case BdryPatch2D::eCircle:
// 	{
// 	  // For circles...
// 	  vMessage(4, "   Creating circle BdryArc #%d (was %d before) ...", iNumSegs, i);

// 	  dynamic_cast<BdryArc *>(pPatch)->vGetCenter(adC);

// 	  BdryPatch2D *pBChk = new BdryArc
// 	    (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	     pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	     dynamic_cast<BdryArc*>(pPatch)->dGetRadius(),
// 	     dynamic_cast<BdryArc*>(pPatch)->iCenterPointIndex ());
// 	  vAddPatchToList (pBChk);
// 	  iNumSegs++;

// 	  vMessage(4, " and added\n");
// 	}
// 	break;
//       case BdryPatch2D::eCircularArc:
// 	{
// 	  // For arcs....
// 	  vMessage(4, "   Creating arc BdryArc #%d (was %d before) ...", iNumSegs, i);

// 	  dynamic_cast<BdryArc*>(pPatch)->vGetCenter(adC);

// 	  BdryPatch2D *pBChk = new BdryArc
// 	    (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	     pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	     pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	     dynamic_cast<BdryArc*>(pPatch)->dGetRadius());
// 	  vAddPatchToList (pBChk);
// 	  iNumSegs++;

// 	  vMessage(4, " and added\n");
// 	}
// 	break;
//       case BdryPatch2D::eLongArc:
// 	{
// 	  // For arcs....
// 	  vMessage(4, "   Creating long arc BdryArc #%d (was %d before) ...", iNumSegs, i);

// 	  dynamic_cast<BdryArc*>(pPatch)->vGetCenter(adC);

// 	  BdryPatch2D *pBChk = new BdryArc
// 	    (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	     pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	     pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	     dynamic_cast<BdryArc*>(pPatch)->dGetRadius(), true);
// 	  vAddPatchToList (pBChk);
// 	  iNumSegs++;

// 	  vMessage(4, " and added\n");
// 	}
// 	break;
//       case BdryPatch2D::eCubicParam:
// 	{
// 	  // Cubic parametric curves
// 	  vMessage(4, "   Creating cubic param curve #%d (was %d before) ...", iNumSegs, i);

// 	  double dax, dbx, dcx, ddx, day, dby, dcy, ddy;
// 	  // Get the parameters
// 	  dynamic_cast<BdryCubicParam2D*>(pPatch)->
// 	    vGetAllCoefficients(&dax, &dbx, &dcx, &ddx, &day, &dby, &dcy, &ddy);
// 	  // Create the new patch
// 	  BdryPatch2D *pBChk = new BdryCubicParam2D
// 	    (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	     pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	     pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	     dax, dbx, dcx, ddx, day, dby, dcy, ddy);
// 	  vAddPatchToList(pBChk);
// 	  iNumSegs++;

// 	  vMessage(4, " and added\n");
// 	}
// 	break;

//       case BdryPatch2D::eBezier:
// 	{
// 	  vMessage(4, "   Creating bezier curve #%d (was %d before) ...",
// 		   iNumSegs, i);

// 	  // Create the new patch
// 	  BdryPatch2D *pBChk = new BdryBezier2D
// 	    (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	     pPatch->iLeftRegion (), pPatch->iRightRegion (),
// 	     pPatch->iBeginPointIndex (), pPatch->iEndPointIndex (),
// 	     dynamic_cast<BdryBezier2D*>(pPatch)->iControlPointIndex(1),
// 	     dynamic_cast<BdryBezier2D*>(pPatch)->iControlPointIndex(2));

// 	  vAddPatchToList(pBChk);
// 	  iNumSegs++;

// 	  vMessage(4, " and added\n");

// 	}

// 	break;
//       case BdryPatch2D::eInterpolated:
// 	{
// 	  vMessage(4, "   Creating interpolated spline #%d (was %d before) ...",
// 		   iNumSegs, i);

// 	  // Get all the points in the list
// 	  List<int> LiPts;
// 	  for (int ii = 0;
// 	       ii < dynamic_cast<BdryInterpolated2D*>(pPatch)->iNumPoints();
// 	       ii++) {
// 	    int iDum = dynamic_cast<BdryInterpolated2D*>(pPatch)->iGetPoint(ii);
// 	    LiPts.vAppendItem( iDum );
// 	  }

// 	  // Create the new patch
// 	  BdryPatch2D *pBChk = new BdryInterpolated2D
// 	    (this, pPatch->iLeftBdryCond (), pPatch->iRightBdryCond (),
// 	     pPatch->iLeftRegion (), pPatch->iRightRegion (), LiPts);
// 	  vAddPatchToList(pBChk);
// 	  iNumSegs++;

// 	  // Now must set the coefficients manually in case they were changed
// 	  // because of a split or something
// 	  double d1, d2, d3, d4, d5, d6, d7, d8;
// 	  assert( dynamic_cast<BdryInterpolated2D*>(pPatch)->iNumSplines() ==
// 		  dynamic_cast<BdryInterpolated2D*>(pBChk)->iNumSplines());

// 	  for (int iCP = 0;
// 	       iCP < dynamic_cast<BdryInterpolated2D*>(pBChk)->iNumSplines();
// 	       iCP++) {
// 	    dynamic_cast<BdryInterpolated2D*>(pPatch)->
// 	      vGetSplineCoefficients(iCP, &d1, &d2, &d3, &d4,
// 				     &d5, &d6, &d7, &d8);
// 	    dynamic_cast<BdryInterpolated2D*>(pBChk)->
// 	      vSetSplineCoefficients(iCP, d1, d2, d3, d4,
// 				     d5, d6, d7, d8);
// 	  }

// 	  dynamic_cast<BdryInterpolated2D*>(pBChk)->vRefresh();
// 	  vMessage(4, " and added\n");
// 	}

// 	break;
//       case BdryPatch2D::eUnknown:
// 	assert(0);
// 	break;
//       } // end switch..
//     }

//   }

//   return (*this);
// }

// // Second arg not currently used, but intended to be used to merge
// // patches as appropriate.
// Bdry2D::
// Bdry2D (const char strBaseFileName[], const bool /*qOptimizePatches*/)
//   : BdryRep(), iSz(0), LpBPatches(), LdPoints(), iSizePoints(0), iCount(0),
//     iNumPointsUsed(0), LqBdryPoint(), iNum2DPatchesUsed(0),
//     LsVertConnectivity(), iNumVertsConnectFound(0), iMaxVertsConnectFound(0),
//     LiNumPatchesToPoint(), LiStartVert(), LqVertSmallAngle(), LqVertUsed()
// {
//   FILE *pFileInput = NULL;
//   char strInFileName[FILE_NAME_LEN];

//   // The following is a hack, although apparently necessary...
//   LiNumPatchesToPoint.vReset();
//   LiStartVert.vReset();

//   vMakeFileName (strInFileName, "%s.bdry", strBaseFileName,
// 		 "Bdry2D::Bdry2D(const string)");

//   pFileInput = fopen (strInFileName, "r");
//   if (pFileInput == NULL) {
//     char strErrorMessage[FILE_NAME_LEN + 80];
//     int iLen = snprintf (strErrorMessage, FILE_NAME_LEN + 80,
// 			 "Couldn't open geometry file for reading: %s\n",
// 			 strInFileName);
//     if ((-1 == iLen) || ((FILE_NAME_LEN + 80) == iLen)) {
//       vFatalError ("Couldn't open geometry file for reading\n",
// 		   "Bdry2D::Bdry2D(const string)");
//     }
//     else {
//       vFatalError (strErrorMessage, "Bdry2D::Bdry2D(const string)");
//     }
//   }
//   vMessage (0, "Reading 2D boundary geometry file %s\n", strInFileName);

//   int iNPoints = -1, iNBdryPatches = -1;
//   const int iBufSize = 2048;
//   char acBuffer[iBufSize];

//   //@@ Read size data
//   vGetLineOrAbort (acBuffer, iBufSize, pFileInput);
//   int iDataRead = sscanf (acBuffer, "%d%d%d", &iNPoints, &iNBdryPatches);
//   if (iDataRead != 2) {
//     // Didn't even get 2 pieces of data.
//     vFatalError ("Problem reading total number of verts and boundary patches",
// 		 "2d boundary geometry input");
//   }
//   if (iNPoints <= 0) {
//     vFatalError("Input file claims number of points <= 0!",
// 		"2d boundary geometry input");
//   }
//   if (iNBdryPatches <= 0) {
//     vFatalError("Input file claims number of patches <= 0!",
// 		"2d boundary geometry input");
//   }

//   //  apBPatches = NULL;
//   vSetSize (iNBdryPatches, iNPoints);

//   //@@ Read coordinate data
//   int i;
//   for (i = 0; i < iNPoints; i++) {
//     vGetLineOrAbort (acBuffer, iBufSize, pFileInput);
//     double adCoord[2];
//     if (2 != sscanf (acBuffer, "%lf %lf", &(adCoord[0]), &(adCoord[1]))) {
//       vMessage (0, "Problem with input line: %s\n", acBuffer);
//       vFatalError ("Insufficient data on line while reading vertex coordinates",
// 		   "2D boundary geometry input");
//     }
//     //    VAVerts[i].vSetCoords (2, adCoord);
//     vAddPoint(adCoord);
//     vSetBdryPoint(i, false);
//   }
//   /*
//     for (i = 0; i < iNumPoints(); i++) {
//     double adLoc[2];
//     vGetPoint(i, adLoc);
//     vMessage(0, "Point %d: (%f, %f)\n", i, adLoc[0], adLoc[1]);
//     }
//   */
//   // Switch to reading the file directly, so that line wrap will happen
//   // automatically.
//   int  iBdryPatch = 0;
//   bool qFoundType;
//   BdryPatch2D::eGeom2D eGeomType = BdryPatch2D::eUnknown;
//   //@@ Read boundary geometry data.
//   do {
//     //@@@ Identify the type of geometry entity
//     char strLineType[33], strLowCase[33];
//     char strDummyLine[80];
//     if (1 != fscanf (pFileInput, "%32s", strLineType)) {
//       vFatalError ("Error reading boundary geometry data",
// 		   "Bdry2D::Bdry2D(const string)");
//     }

//     qFoundType = false;

//     if (strLineType[0] == '#') {
//       // Skip the line
//       int iRes = fscanf (pFileInput, "%s\n", strDummyLine);
//       // Don't really care about the value...
//       if (iRes < 0) 
// 	vFatalError ("Hit EOF reading boundary geometry data",
// 		     "Bdry2D::Bdry2D(const string)");
//     }
//     else {
//       // Only take lines that are not comments..!
//       for (i = 0; strLineType[i] != '\000' && i < 33; i++)
// 	strLowCase[i] = tolower (strLineType[i]);

//       if (strncmp (strLowCase, "polyline", 8) == 0) {
// 	eGeomType = BdryPatch2D::ePolyLine;
// 	qFoundType = true;
//       }
//       else if (strncmp (strLowCase, "arc", 3) == 0) {
// 	eGeomType = BdryPatch2D::eCircularArc;
// 	qFoundType = true;
//       }
//       else if (strncmp (strLowCase, "longarc", 3) == 0) {
// 	eGeomType = BdryPatch2D::eLongArc;
// 	qFoundType = true;
//       }
//       else if (strncmp (strLowCase, "circle", 6) == 0) {
// 	eGeomType = BdryPatch2D::eCircle;
// 	qFoundType = true;
//       }
//       else if (strncmp (strLowCase, "cubicparam", 10) == 0) {
// 	eGeomType = BdryPatch2D::eCubicParam;
// 	qFoundType = true;
//       }
//       else if (strncmp (strLowCase, "bezier", 6) == 0) {
// 	eGeomType = BdryPatch2D::eBezier;
// 	qFoundType = true;
//       }
//       else if (strncmp (strLowCase, "spline", 6) == 0) {
// 	eGeomType = BdryPatch2D::eInterpolated;
// 	qFoundType = true;
//       }

//       if (!qFoundType) {
// 	vMessage (0, "Unknown geometric entity name: %s\n", strLineType);
// 	vFatalError ("Error reading boundary geometry data",
// 		     "Bdry2D::Bdry2D(const string)");
//       }

//       vMessage(1, "Reading patch #%d\n", iBdryPatch);

//       //@@@ Regardless, read the BC and/or region data.
//       char cTypeL, cTypeR;
//       int iValueL, iValueR;
//       if (4 != fscanf (pFileInput, " %1c%d %1c%d", &cTypeL,
// 		       &iValueL, &cTypeR, &iValueR)) {
// 	vMessage (0, "Bad boundary condition or region specification in:\n%s\n",
// 		  acBuffer);
// 	vFatalError ("Error reading boundary geometry data",
// 		     "Bdry2D::Bdry2D(const string)");
//       }

//       int iBCL, iBCR, iRegL, iRegR;
//       if ((cTypeL == 'b') || (cTypeL == 'B')) {
// 	iBCL = iValueL;
// 	iRegL = iInvalidRegion;
//       }
//       else {
// 	assert (cTypeL == 'r' || cTypeL == 'R');
// 	iBCL = iInvalidBC;
// 	iRegL = iValueL;
//       }
//       if ((cTypeR == 'b') || (cTypeR == 'B')) {
// 	iBCR = iValueR;
// 	iRegR = iInvalidRegion;
//       }
//       else {
// 	assert (cTypeR == 'r' || cTypeR == 'R');
// 	iBCR = iInvalidBC;
// 	iRegR = iValueR;
//       }

//       // Here, should really make the patches themselves know
//       // how to read the line describing them.. Or should we?

//       // Or at least push that down to the vAddPatchToList..
//       // Now it justs duplicates code.

//       switch (eGeomType) {
//       case BdryPatch2D::ePolyLine:
// 	{
// 	  int iSeg, iNumSegs;
// 	  int iV0, iV1;

// 	  if (1 != fscanf (pFileInput, "%d", &iNumSegs)) {
// 	    vFatalError ("Error reading boundary geometry data for polylines",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }
// 	  iNumSegs--;
// 	  if (1 != fscanf (pFileInput, "%d", &iV1))
// 	    vFatalError ("Error reading boundary geometry data for polylines", "Bdry2D::Bdry2D(const string)");

// 	  for (iSeg = 0; iSeg < iNumSegs; iSeg++) {
// 	    iV0 = iV1;
// 	    if (1 != fscanf (pFileInput, "%d", &iV1))
// 	      vFatalError ("Error reading boundary geometry data for polylines",
// 			   "Bdry2D::Bdry2D(const string)");

// 	    BdrySeg *pBS = new BdrySeg (this, iBCL, iBCR, iRegL, iRegR,
// 					iV0, iV1);
// 	    vAddPatchToList (pBS);

// 	    // Must still increment..
// 	    iBdryPatch++;

// 	    double adPt[2];
// 	    double adVertPt[2];
// 	    double adNorm;
// 	    double dTanX, dTanY;
// 	    structVertConnect sVConnect;

// 	    // Might as well build the connectivity list as we go along, right?
// 	    sVConnect.iPatch = iNumPatches()-1;
// 	    sVConnect.iVert  = iV0;
// 	    vSetBdryPoint(iV0, true);
// 	    // Get the tangent vector for iV0
// 	    pBS->vPointAtParameterRatio(0.01, adPt);
// 	    // Get iV0...
// 	    vGetPoint(iV0, adVertPt);
// 	    dTanX = adPt[0] - adVertPt[0];
// 	    dTanY = adPt[1] - adVertPt[1];
// 	    adNorm = sqrt(dTanX*dTanX + dTanY*dTanY);
// 	    dTanX /= adNorm;
// 	    dTanY /= adNorm;
// 	    sVConnect.dAngle = atan2(dTanY, dTanX);
// 	    // Add it to list
// 	    vAddToConnectivityList(sVConnect);

// 	    // Do the same for iV1 too
// 	    sVConnect.iPatch = iBdryPatch-1;
// 	    sVConnect.iVert  = iV1;
// 	    vSetBdryPoint(iV1, true);
// 	    // Get the tangent vector for iV0
// 	    pBS->vPointAtParameterRatio(0.99, adPt);
// 	    vGetPoint(iV1, adVertPt);
// 	    dTanX = adPt[0] - adVertPt[0];
// 	    dTanY = adPt[1] - adVertPt[1];
// 	    adNorm = sqrt(dTanX*dTanX + dTanY*dTanY);
// 	    dTanX /= adNorm;
// 	    dTanY /= adNorm;
// 	    sVConnect.dAngle = atan2(dTanY, dTanX);
// 	    // Add it to list
// 	    vAddToConnectivityList(sVConnect);
// 	  }
// 	  break;
// 	}			// End case for polylines

//       case BdryPatch2D::eCircularArc:
// 	{
// 	  // Arcs are defined like this:
// 	  // arc b 1 r 1 Radius VertIndexBegin VertIndexEnd
// 	  double dRadius;
// 	  int    iV0, iV1;

// 	  // First item is always the radius
// 	  if (1 != fscanf (pFileInput, "%lf", &dRadius)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }
// 	  // 2nd & 3rd Items are the indices...
// 	  if (1 != fscanf (pFileInput, "%d", &iV0)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%d", &iV1)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // Initialize the arc...
// 	  BdryArc *pBA = new BdryArc (this, iBCL, iBCR, iRegL, iRegR,
// 				      iV0, iV1, dRadius);
// 	  vAddPatchToList (pBA);

// 	  // Must still increment..
// 	  iBdryPatch++;

// 	  double adPt[2], adVertPt[2];
// 	  double adNorm;
// 	  double dTanX, dTanY;
// 	  structVertConnect sVConnect;

// 	  // Might as well build the connectivity list as we go along, right?
// 	  sVConnect.iPatch = iNumPatches()-1;
// 	  sVConnect.iVert  = iV0;
// 	  vSetBdryPoint(iV0, true);
// 	  // Get the tangent vector for iV0
// 	  pBA->vPointAtParameterRatio(0.01, adPt);
// 	  vGetPoint(iV0, adVertPt);
// 	  dTanX = adPt[0] - adVertPt[0];
// 	  dTanY = adPt[1] - adVertPt[1];
// 	  adNorm = sqrt(dTanX*dTanX + dTanY*dTanY);
// 	  dTanX /= adNorm;
// 	  dTanY /= adNorm;
// 	  sVConnect.dAngle = atan2(dTanY, dTanX);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  // Do the same for iV1 too
// 	  sVConnect.iPatch = iBdryPatch-1;
// 	  sVConnect.iVert  = iV1;
// 	  vSetBdryPoint(iV1, true);
// 	  // Get the tangent vector for iV0
// 	  pBA->vPointAtParameterRatio(0.99, adPt);
// 	  vGetPoint(iV1, adVertPt);
// 	  dTanX = adPt[0] - adVertPt[0];
// 	  dTanY = adPt[1] - adVertPt[1];
// 	  adNorm = sqrt(dTanX*dTanX + dTanY*dTanY);
// 	  dTanX /= adNorm;
// 	  dTanY /= adNorm;
// 	  sVConnect.dAngle = atan2(dTanY, dTanX);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  break;
// 	} // End case for circular arcs

//       case BdryPatch2D::eLongArc:
// 	{
// 	  // Long arcs are defined like this:
// 	  // longarc b 1 r 1 Radius VertIndexBegin VertIndexEnd
// 	  // Unlike regular arcs, these span more than 180 degrees.
// 	  double dRadius;
// 	  int    iV0, iV1;

// 	  // First item is always the radius
// 	  if (1 != fscanf (pFileInput, "%lf", &dRadius)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }
// 	  // 2nd & 3rd Items are the indices...
// 	  if (1 != fscanf (pFileInput, "%d", &iV0)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%d", &iV1)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // Initialize the arc...
// 	  BdryArc *pBA = new BdryArc (this, iBCL, iBCR, iRegL, iRegR,
// 				      iV0, iV1, dRadius, true);
// 	  vAddPatchToList (pBA);

// 	  // Must still increment..
// 	  iBdryPatch++;

// 	  double adPt[2], adVertPt[2];
// 	  double adNorm;
// 	  double dTanX, dTanY;
// 	  structVertConnect sVConnect;

// 	  // Might as well build the connectivity list as we go along, right?
// 	  sVConnect.iPatch = iNumPatches()-1;
// 	  sVConnect.iVert  = iV0;
// 	  vSetBdryPoint(iV0, true);
// 	  // Get the tangent vector for iV0
// 	  pBA->vPointAtParameterRatio(0.01, adPt);
// 	  vGetPoint(iV0, adVertPt);
// 	  dTanX = adPt[0] - adVertPt[0];
// 	  dTanY = adPt[1] - adVertPt[1];
// 	  adNorm = sqrt(dTanX*dTanX + dTanY*dTanY);
// 	  dTanX /= adNorm;
// 	  dTanY /= adNorm;
// 	  sVConnect.dAngle = atan2(dTanY, dTanX);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  // Do the same for iV1 too
// 	  sVConnect.iPatch = iBdryPatch-1;
// 	  sVConnect.iVert  = iV1;
// 	  vSetBdryPoint(iV1, true);
// 	  // Get the tangent vector for iV0
// 	  pBA->vPointAtParameterRatio(0.99, adPt);
// 	  vGetPoint(iV1, adVertPt);
// 	  dTanX = adPt[0] - adVertPt[0];
// 	  dTanY = adPt[1] - adVertPt[1];
// 	  adNorm = sqrt(dTanX*dTanX + dTanY*dTanY);
// 	  dTanX /= adNorm;
// 	  dTanY /= adNorm;
// 	  sVConnect.dAngle = atan2(dTanY, dTanX);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  break;
// 	} // End case for circular arcs

//       case BdryPatch2D::eCircle:
// 	{
// 	  // Circles are defined like this:
// 	  // arc b 1 r 1 radius center_x center_y
// 	  double dRadius;
// 	  int iC;

// 	  // First item is always the radius
// 	  if (1 != fscanf (pFileInput, "%lf", &dRadius)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // Second Items is the center index
// 	  if (1 != fscanf (pFileInput, "%d", &iC)) {
// 	    vFatalError ("Error reading boundary geometry data for circular arcs",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

//  	  // A magic value so that this point isn't added to the mesh.
//  	  LiNumPatchesToPoint[iC] = -1;

// 	  // Initialize the arc...
// 	  BdryArc *pBA = new BdryArc (this, iBCL, iBCR, iRegL, iRegR,
// 				      dRadius, iC);

// 	  vAddPatchToList (pBA);

// 	  // Must still increment..
// 	  iBdryPatch++;

// 	  break;
// 	}

//       case BdryPatch2D::eCubicParam:
// 	{
// 	  // Generic cubic parametric curve. Expect 10 arguments,
// 	  // i.e. beginning and end vertices index and then the 8
// 	  // parameters Ax Bx Cx Dx, Ay By Cy Dy.
// 	  // Dx must be equal to pV0->dX()
// 	  // Dy must be equal to pV1->dY()

// 	  double dAx, dAy, dBx, dBy, dCx, dCy, dDx, dDy;
// 	  int    iV0, iV1;

// 	  // Definition looks like this
// 	  // cubicparam b 1 r 1 iV0 iV1 ax bx cx dx ay by cy dy

// 	  // First item is iV0
// 	  if (1 != fscanf (pFileInput, "%d", &iV0)) {
// 	    vFatalError ("Error reading boundary geometry data (iV0) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // 2nd item is iV1
// 	  if (1 != fscanf (pFileInput, "%d", &iV1)) {
// 	    vFatalError ("Error reading boundary geometry data (iV1) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // The 8 next items are the parameters
// 	  if (1 != fscanf (pFileInput, "%lf", &dAx)) {
// 	    vFatalError ("Error reading boundary geometry data (Ax) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dBx)) {
// 	    vFatalError ("Error reading boundary geometry data (Bx) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dCx)) {
// 	    vFatalError ("Error reading boundary geometry data (Cx) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dDx)) {
// 	    vFatalError ("Error reading boundary geometry data (Dx) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dAy)) {
// 	    vFatalError ("Error reading boundary geometry data (Ay) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dBy)) {
// 	    vFatalError ("Error reading boundary geometry data (By) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dCy)) {
// 	    vFatalError ("Error reading boundary geometry data (Cy) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%lf", &dDy)) {
// 	    vFatalError ("Error reading boundary geometry data (Dy) for cubic parametric curves",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // Initialize the curve...
// 	  BdryCubicParam2D *pBCP =
// 	    new BdryCubicParam2D (this, iBCL, iBCR, iRegL, iRegR, iV0, iV1,
// 				  dAx, dBx, dCx, dDx, dAy, dBy, dCy, dDy);

// 	  vAddPatchToList (pBCP);

// 	  // Must still increment..
// 	  iBdryPatch++;

// 	  // Make sure connectivity is stored too...
// 	  double adTan[2];
// 	  structVertConnect sVConnect;

// 	  sVConnect.iPatch = iNumPatches()-1;
// 	  sVConnect.iVert  = iV0;
// 	  vSetBdryPoint(iV0, true);
// 	  // Get the tangent vector for iV0
// 	  pBCP->vTangentAtParameterRatio(0.0, adTan);
// 	  sVConnect.dAngle = atan2(adTan[1], adTan[0]);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  // Do the same for iV1
// 	  sVConnect.iVert  = iV1;
// 	  vSetBdryPoint(iV1, true);
// 	  // Get the tangent vector for iV1
// 	  pBCP->vTangentAtParameterRatio(1.0, adTan);
// 	  sVConnect.dAngle = atan2(-adTan[1], -adTan[0]);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  break;
// 	}

//       case BdryPatch2D::eBezier:
// 	{
// 	  // Bezier curves. Arguments are in the following format:
// 	  // bezier r 1 b 1 iV0 iV1 iC1 iC2

// 	  int    iV0, iV1, iC1, iC2;

// 	  // First item is iV0
// 	  if (1 != fscanf (pFileInput, "%d", &iV0)) {
// 	    vFatalError ("Error reading boundary geometry data (iV0) for bezier curve",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // 2nd item is iV1
// 	  if (1 != fscanf (pFileInput, "%d", &iV1)) {
// 	    vFatalError ("Error reading boundary geometry data (iV1) for bezier curve",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // The next 2 items are the control points
// 	  if (1 != fscanf (pFileInput, "%d", &iC1)) {
// 	    vFatalError ("Error reading boundary geometry data (C1x) for bezier curve",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  if (1 != fscanf (pFileInput, "%d", &iC2)) {
// 	    vFatalError ("Error reading boundary geometry data (C1y) for bezier curve",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // A magic value so that the control points aren't added to the mesh.
// 	  LiNumPatchesToPoint[iC1] = -1;
// 	  LiNumPatchesToPoint[iC2] = -1;

// 	  // Initialize the curve...
// 	  BdryBezier2D *pBB = new BdryBezier2D (this, iBCL, iBCR, iRegL, iRegR,
// 				      iV0, iV1, iC1, iC2);

// 	  vAddPatchToList (pBB);

// 	  // Must still increment..
// 	  iBdryPatch++;

// 	  // Make sure connectivity is stored too...
// 	  double adTan[2];
// 	  structVertConnect sVConnect;

// 	  sVConnect.iPatch = iNumPatches()-1;
// 	  sVConnect.iVert  = iV0;
// 	  vSetBdryPoint(iV0, true);
// 	  // Get the tangent vector for iV0
// 	  pBB->vTangentAtParameterRatio(0.0, adTan);
// 	  sVConnect.dAngle = atan2(adTan[1], adTan[0]);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  // Do the same for iV1
// 	  sVConnect.iVert  = iV1;
// 	  vSetBdryPoint(iV1, true);
// 	  // Get the tangent vector for iV1
// 	  pBB->vTangentAtParameterRatio(1.0, adTan);
// 	  sVConnect.dAngle = atan2(-adTan[1], -adTan[0]);
// 	  // Add it to list
// 	  vAddToConnectivityList(sVConnect);

// 	  break;
// 	}

//       case BdryPatch2D::eInterpolated:
// 	{
// 	  // Interpolated splines...
// 	  // format is spline r 1 b 1 iNumPoints iPt0 iPt1 iPt2 ... iPtNumPoints-1

// 	  int		iNP, iPt;
// 	  List<int> LiAllPoints;

// 	  // First item is iNP
// 	  if (1 != fscanf (pFileInput, "%d", &iNP)) {
// 	    vFatalError ("Error reading boundary geometry data (iNP) for spline curve",
// 			 "Bdry2D::Bdry2D(const string)");
// 	  }

// 	  // now get all the points
// 	  int j;
// 	  for (j = 0; j < iNP; j++) {
// 	    // Read the next point
// 	    if (1 != fscanf (pFileInput, "%d", &iPt)) {
// 	      vFatalError ("Error reading points for spline curve",
// 			   "Bdry2D::Bdry2D(const string)");
// 	    }
// 	    LiAllPoints.vAppendItem(iPt);
// 	    if (j > 0 && j < iNP-1) {
// 	      // A magic value so that the control points aren't added to the mesh.
// 	      LiNumPatchesToPoint[iPt] = -1;
// 	    }
// 	  }

// 	  // Initialize the curve...
// 	  BdryInterpolated2D *pBS = new BdryInterpolated2D (this, iBCL,
// 							    iBCR, iRegL,
// 							    iRegR,
// 							    LiAllPoints);

// 	  // Don't add it to the list just yet...
// 	  //	  vAddPatchToList (pBS);

// 	  // Must still increment..
// 	  iBdryPatch++;

// 	  // Splines are special because it is possible they will be
// 	  // transformed into several splines. So we have to do this
// 	  // first...

// 	  double dTotalChange = pBS->dTotalOrientationChange();
// 	  int iNumSplines = static_cast<int>(floor(dTotalChange)) + 1;
// 	  double dChange = dTotalChange/iNumSplines;
// 	  double dt = 0.0;
// 	  List<double> LdNewPoints;

// 	  LdNewPoints.vClear();

// 	  vMessage(4, "Spline splitting info\n--------------------------\n");
// 	  vMessage(4, "New Point list:\n");
// 	  // Find these points..
// 	  int iS;
// 	  for (iS = 1; iS < iNumSplines; iS++) {
// 	    double dEnd = pBS->dParameterAtOrientationChange(dChange, dt);
// 	    dt = dEnd;
// 	    LdNewPoints.vAppendItem(dt);
// 	    vMessage(4, "   %f\n", dt);
// 	  }

// 	  // LiAllPoints contains all the control points...
// 	  // first one is at parameter 0.0, next one at 1.0,
// 	  // and so on.

// 	  // This contains a list of the parameters used for the new
// 	  // splines i.e. LLdPointsInSplines[0] contains a list of
// 	  // parameters used for the new spline 0.
// 	  // LLdPointsInSplines[1][3] contains the 4th parameter for
// 	  // spline 1, etc..

// 	  List < List <double> > LLdPointsInSplines;
// 	  List < List <double> > LLdCoeffsInSplines;
// 	  List < List <int>    > LLiPtIdxInSplines;
// 	  LLdPointsInSplines.vClear();

// 	  List <double> LdSplinePoints;
// 	  List <double> LdCoeffs;
// 	  List <int>    LiIndex;
// 	  double dLastItem = 0.0;

// 	  // Fill the list of splines...
// 	  for (iS = 0; iS < iNumSplines; iS++) {
// 	    // Clear the list
// 	    LdSplinePoints.vClear();
// 	    // Add the first item
// 	    LdSplinePoints.vAppendItem(dLastItem);
// 	    double dNextPoint;
// 	    j = 0;
// 	    if (iS < iNumSplines-1){
// 	      dNextPoint = LdNewPoints[iS];
// 	    }
// 	    else {
// 	      dNextPoint = pBS->dMaxParameter();
// 	    }
// 	    if (floor(dLastItem) < floor(dNextPoint)) {
// 				// insert all the control points in between..
// 	      for (j = static_cast<int>(floor(dLastItem))+1;
// 		   j <= static_cast<int>(floor(dNextPoint)); j++) {
// 		LdSplinePoints.vAppendItem(j);
// 	      }
// 	    }
// 	    dLastItem = dNextPoint;
// 	    if (iFuzzyComp(dNextPoint, (j-1)) != 0) {
// 				// Add the next point
// 	      LdSplinePoints.vAppendItem(dNextPoint);
// 	    }
// 	    // Add that to the list
// 	    LLdPointsInSplines.vAppendItem(LdSplinePoints);
// 	  }

// 	  // Now build the list of coefficients..
// 	  for (iS = 0; iS < LLdPointsInSplines.iLength(); iS++) {
// 	    // Clear the list...
// 	    LdCoeffs.vClear();
// 	    for (int iP = 0; iP < LLdPointsInSplines[iS].iLength()-1; iP++) {
// 	      double dFirst = LLdPointsInSplines[iS][iP];
// 	      int iFirst = static_cast<int>(floor(dFirst));
// 	      // Now add the coefficients..
// 	      double d1, d2, d3, d4, d5, d6, d7, d8;
// 	      pBS->vGetSplineCoefficients(iFirst, &d1, &d2, &d3, &d4,
// 					  &d5, &d6, &d7, &d8);
// 	      // Add them...
// 	      LdCoeffs.vAppendItem(d1);
// 	      LdCoeffs.vAppendItem(d2);
// 	      LdCoeffs.vAppendItem(d3);
// 	      LdCoeffs.vAppendItem(d4);
// 	      LdCoeffs.vAppendItem(d5);
// 	      LdCoeffs.vAppendItem(d6);
// 	      LdCoeffs.vAppendItem(d7);
// 	      LdCoeffs.vAppendItem(d8);
// 	    }
// 	    LLdCoeffsInSplines.vAppendItem(LdCoeffs);
// 	  }

// 	  // Now build the list of points indices..
// 	  for (iS = 0; iS < LLdPointsInSplines.iLength(); iS++) {
// 	    LiIndex.vClear();
// 	    for (int iP = 0; iP < LLdPointsInSplines[iS].iLength(); iP++) {
// 	      double dAL = LLdPointsInSplines[iS][iP];
// 	      double dMantissa = dAL - (floor(dAL));
// 	      if (iFuzzyComp(dMantissa, 0.0) == 0) {
// 		// That point exists..
// 		int iAL = static_cast<int>(floor(dAL));
// 		LiIndex.vAppendItem(LiAllPoints[iAL]);
// 	      }
// 	      else {
// 		if (iP != 0) {
// 		  // We need to create it.
// 		  // Get the coords
// 		  double adC[2];
// 		  pBS->vPointAtParameter(dAL, adC);
// 		  int iIndex = iNumPoints();
// 		  vAddPoint(adC);
// 		  LiIndex.vAppendItem(iIndex);
// 		}
// 		else {
// 		  // That point has already been created
// 		  int iL = LLiPtIdxInSplines[iS-1].iLength();
// 		  LiIndex.vAppendItem(LLiPtIdxInSplines[iS-1][iL-1]);
// 		}
// 	      }
// 	    }
// 	    LLiPtIdxInSplines.vAppendItem(LiIndex);
// 	  }

// 	  vMessage(4, "Contents of spline list\n");
// 	  // Check to see what we got..
// 	  for (int ii = 0; ii < LLdPointsInSplines.iLength(); ii++) {
// 	    vMessage(4, "Spline %d:\n", ii);
// 	    List <double> LdP = LLdPointsInSplines[ii];
// 	    for (int jj = 0; jj < LdP.iLength(); jj++) {
// 	      double adPt[2];
// 	      vGetPoint(LLiPtIdxInSplines[ii][jj], adPt);
// 	      vMessage(4, "   %f -- Index: %d (%f, %f)\n", LdP[jj],
// 		       LLiPtIdxInSplines[ii][jj], adPt[0], adPt[1]);
// 	    }
// 	  }

// 	  int iLastSpline = -1;
// 	  double dLastCut = 0.0;

// 	  // Now build each spline...
// 	  for (int iP = 0; iP < LdNewPoints.iLength(); iP++) {
// 	    // Split it..
// 	    double d11, d12, d13, d14, d15, d16, d17, d18;
// 	    double d21, d22, d23, d24, d25, d26, d27, d28;
// 	    int iSpline;
// 	    double dLeft, dNewLastCut;
// 	    iSpline = static_cast<int>(floor(LdNewPoints[iP]));
// 	    dLeft = LdNewPoints[iP] - iSpline;
// 	    dNewLastCut = dLeft;
// 	    if (iLastSpline == iSpline) {
// 	      dLeft = (dLeft - dLastCut)/(1.0 - dLastCut);
// 	    }
// 	    iLastSpline = iSpline;
// 	    dLastCut = dNewLastCut;
// 	    // Spline index...
// 	    int iSi = iP;
// 	    if (iFuzzyComp(dLeft, 0.0) == 0) {
// 				dLeft = 1.0;
// 				iSi--;
// 	    }
// 	    // Get the last coefficients from Spline iP..
// 	    int iBase = LLdCoeffsInSplines[iSi].iLength()/8;
// 	    iBase--;
// 	    iBase *= 8;
// 	    d11 = LLdCoeffsInSplines[iSi][iBase+0];
// 	    d12 = LLdCoeffsInSplines[iSi][iBase+1];
// 	    d13 = LLdCoeffsInSplines[iSi][iBase+2];
// 	    d14 = LLdCoeffsInSplines[iSi][iBase+3];
// 	    d15 = LLdCoeffsInSplines[iSi][iBase+4];
// 	    d16 = LLdCoeffsInSplines[iSi][iBase+5];
// 	    d17 = LLdCoeffsInSplines[iSi][iBase+6];
// 	    d18 = LLdCoeffsInSplines[iSi][iBase+7];
// 	    // Get also the first coefficients from spline iP+1
// 	    d21 = LLdCoeffsInSplines[iSi+1][0];
// 	    d22 = LLdCoeffsInSplines[iSi+1][1];
// 	    d23 = LLdCoeffsInSplines[iSi+1][2];
// 	    d24 = LLdCoeffsInSplines[iSi+1][3];
// 	    d25 = LLdCoeffsInSplines[iSi+1][4];
// 	    d26 = LLdCoeffsInSplines[iSi+1][5];
// 	    d27 = LLdCoeffsInSplines[iSi+1][6];
// 	    d28 = LLdCoeffsInSplines[iSi+1][7];
// 	    // Split..
// 	    pBS->vSplitSplineAtParameter(iSpline, dLeft, &d11, &d12,
// 					 &d13, &d14, &d15, &d16, &d17,
// 					 &d18, &d21, &d22, &d23, &d24,
// 					 &d25, &d26, &d27, &d28);
// 	    // And put that back in the lists...
// 	    LLdCoeffsInSplines[iSi][iBase+0] = d11;
// 	    LLdCoeffsInSplines[iSi][iBase+1] = d12;
// 	    LLdCoeffsInSplines[iSi][iBase+2] = d13;
// 	    LLdCoeffsInSplines[iSi][iBase+3] = d14;
// 	    LLdCoeffsInSplines[iSi][iBase+4] = d15;
// 	    LLdCoeffsInSplines[iSi][iBase+5] = d16;
// 	    LLdCoeffsInSplines[iSi][iBase+6] = d17;
// 	    LLdCoeffsInSplines[iSi][iBase+7] = d18;
// 	    // And this one too....
// 	    LLdCoeffsInSplines[iSi+1][0] = d21;
// 	    LLdCoeffsInSplines[iSi+1][1] = d22;
// 	    LLdCoeffsInSplines[iSi+1][2] = d23;
// 	    LLdCoeffsInSplines[iSi+1][3] = d24;
// 	    LLdCoeffsInSplines[iSi+1][4] = d25;
// 	    LLdCoeffsInSplines[iSi+1][5] = d26;
// 	    LLdCoeffsInSplines[iSi+1][6] = d27;
// 	    LLdCoeffsInSplines[iSi+1][7] = d28;
// 	  }

// 	  // I think we are ready to create them now...

// 	  for (iS = 0; iS < iNumSplines; iS++) {
// 	    BdryInterpolated2D *pNewS;
// 	    List<int> LiSplinePoints = LLiPtIdxInSplines[iS];
// 	    LdCoeffs = LLdCoeffsInSplines[iS];

// 	    // Create the patch
// 	    pNewS = new BdryInterpolated2D(this, iBCL, iBCR, iRegL,
// 					   iRegR, LiSplinePoints);

// 	    vAddPatchToList(pNewS);

// 	    // Change the coefficients..
// 	    for (int iNewS = 0; iNewS < pNewS->iNumSplines(); iNewS++) {
// 	      int iBase = 8*iNewS;
// 	      pNewS->vSetSplineCoefficients(iNewS, LdCoeffs[iBase],
// 					    LdCoeffs[iBase+1],
// 					    LdCoeffs[iBase+2],
// 					    LdCoeffs[iBase+3],
// 					    LdCoeffs[iBase+4],
// 					    LdCoeffs[iBase+5],
// 					    LdCoeffs[iBase+6],
// 					    LdCoeffs[iBase+7]);
// 	    }

// 	    pNewS->vRefresh();
// 	    double adB[2], adE[2];
// 	    pNewS->vGetBeginPoint(adB);
// 	    pNewS->vGetEndPoint(adE);

// 	    vMessage(4, "Total orientation change: %f\n",
// 		     pNewS->dTotalOrientationChange());
// 	    vMessage(4, "Starts at: (%f, %f), ends at (%f, %f)\n",
// 		     adB[0], adB[1], adE[0], adE[1]);
// 		  vMessage(4, "Roots:\n");
// 	    for (int iRt = 0; iRt < pNewS->iNumOrderedRoots(); iRt++) {
// 	      double dRoot = pNewS->dGetOrderedRoot(iRt);
// 	      double adC[2];
// 	      pNewS->vPointAtParameter(dRoot, adC);
// 	      vMessage(4, "%d: %3.16e -- (%f, %f)\n", iRt, dRoot,
// 		       adC[0], adC[1]);
// 	    }

// 	    // Make sure connectivity is stored too...
// 	    double adTan[2];
// 	    structVertConnect sVConnect;

// 	    sVConnect.iPatch = iNumPatches()-1;
// 	    sVConnect.iVert  = LiSplinePoints[0];
// 	    vSetBdryPoint(LiSplinePoints[0], true);
// 	    // Get the tangent vector for iV0
// 	    pNewS->vTangentAtParameterRatio(0.0, adTan);
// 	    sVConnect.dAngle = atan2(adTan[1], adTan[0]);
// 	    // Add it to list
// 	    vAddToConnectivityList(sVConnect);
// 	    // Make sure that the spline end points aren't messed up.
// 	    LiNumPatchesToPoint[sVConnect.iVert] =
// 	      max(LiNumPatchesToPoint[sVConnect.iVert], 0);

// 	    // Do the same for iV1
// 	    sVConnect.iVert  = LiSplinePoints[LiSplinePoints.iLength()-1];
// 	    vSetBdryPoint(LiSplinePoints[LiSplinePoints.iLength()-1], true);
// 	    // Get the tangent vector for iV1
// 	    pNewS->vTangentAtParameterRatio(1.0, adTan);
// 	    sVConnect.dAngle = atan2(-adTan[1], -adTan[0]);
// 	    // Add it to list
// 	    vAddToConnectivityList(sVConnect);
// 	    // Make sure that the spline end points aren't messed up.
// 	    LiNumPatchesToPoint[sVConnect.iVert] =
// 	      max(LiNumPatchesToPoint[sVConnect.iVert], 0);
// 	  }

// 	  delete pBS;

// 	  break;
// 	}
//       case BdryPatch2D::eUnknown:
// 	assert(0);
// 	break;
//       }				// End switch to read all geometry types.
//     }				// End if to make sure it's not a comment line
//   }  while (iBdryPatch < iNBdryPatches);	// End loop to read all
// 						// boundary patches
//   assert (iBdryPatch == iNBdryPatches);

//   // Optimize them if the user asked for it...
//   //  if (qOptimizePatches)
//   //    vOptimizePatches ();

//   // Order stuff correctly
//   vOrderConnectivityInfo();

//   // Find small angles
//   vFindSmallAngles();

//   // Print out the info
//   for (int kk = 0; kk < iNumPoints(); kk++) {
//     if (qIsBdryPoint(kk)) {
//       vMessage(3, "   Point #%d is a boundary point\n", kk);
//       if (LqVertSmallAngle[kk])
// 	vMessage(3, "      It has a small angle\n");
//       else
// 	vMessage(3, "      It does not have a small angle\n");

//       vMessage(3, "      It has %d patches connected to it\n",
// 	       iNumPatchesConnectedToPoint(kk));

//       for (int nn = 0; nn < iNumPatchesConnectedToPoint(kk); nn++) {
// 	vMessage(3, "         Patch %p\n", pB2DPatchConnectedToPoint(kk, nn));
//       }
//     }
//   }

//   fclose(pFileInput);
// }

Bdry2D::~Bdry2D()
{
  // The bdry patches are created by new, so they must be deleted.
  for (int i = 0; i < iNumPatchesUsed(); i++) {
    delete pBdry(i);
  }
}

void
Bdry2D::vAddPatchToList(BdryPatch2D * pPatch)
{

  LpBPatches.vAppendItem(pPatch);
//  apBPatches[iNPatches] = pPatch;

//  // Tag verts as used.
//  aqVertUsed[iBeginVert(iNPatches)] = true;
//  aqVertUsed[iEndVert(iNPatches)] = true;

// Obsolete
//  asPatchUsed[iCount].qIsUsed = true;
//  asPatchUsed[iCount].iWasReplacedBy = 0;

// increment counter
  iNPatches++;
  iNum2DPatchesUsed++;
}

void
Bdry2D::vFindSmallAngles()
{
  // Now that the connectivity at the vertices has been determined,
  // we can easily determine which vertex has a small angle...

  for (int iV = 0; iV < iNumPoints(); iV++) {
    if (qIsBdryPoint(iV)) {
      int iMaxP = iNumPatchesConnectedToPoint(iV);
      int iP = 0;
      LqVertSmallAngle[iV] = false;
      while ((iP < iMaxP) && (LqVertSmallAngle[iV] == false)) {
	BdryPatch2D *pP1, *pP2;
	int iP1, iP2;
	double dA1, dA2;

	iP1 = iP;
	if (iP < iMaxP - 1) {
	  iP2 = iP + 1;
	}
	else {
	  iP2 = 0;
	}

	dA1 = LsVertConnectivity[LiStartVert[iV] + iP1].dAngle;
	dA2 = LsVertConnectivity[LiStartVert[iV] + iP2].dAngle;

	// Have the two patches to look at, now determine if they form
	// a small angle

	double dMinAngle = M_PI / 3.0;

	if (dA1 < 0.0)
	  dA1 += 2.0 * M_PI;
	if (dA2 < 0.0)
	  dA2 += 2.0 * M_PI;

	//			vMessage(2, "dA1; %f -- dA2: %f -- dA2-dA1: %f\n",
	//									dA1, dA2, dA2-dA1);

	if (fabs(dA2 - dA1) <= dMinAngle) {
	  // Possibly a small angle -- but is it in the interior?
	  // Always interior if the patches have two regions..

	  pP1 = pB2DPatchConnectedToPoint(iV, iP1);
	  pP2 = pB2DPatchConnectedToPoint(iV, iP2);

	  if (pP1->iRegionWhichSide() != REGION_SIDE_BOTH) {
	    double adPt1[2], adPt2[2];
	    double adPt[2];
	    double dAL1, dAL2;

	    vGetPoint(iV, adPt);

	    pP1 = pB2DPatchConnectedToPoint(iV, iP1);
	    pP2 = pB2DPatchConnectedToPoint(iV, iP2);

	    dAL1 = pP1->dParameterAtPoint(adPt);
	    dAL2 = pP2->dParameterAtPoint(adPt);

	    // Get the proper points

	    if (dAL1 > (pP1->dMaxParameter() / 2.0)) {
	      pP1->vPointAtParameterRatio(0.99, adPt1);
	    }
	    else {
	      pP1->vPointAtParameterRatio(0.01, adPt1);
	    }

	    if (pP1 != pP2) {
	      if (dAL2 > (pP2->dMaxParameter() / 2.0)) {
		pP2->vPointAtParameterRatio(0.99, adPt2);
	      }
	      else {
		pP2->vPointAtParameterRatio(0.01, adPt2);
	      }
	    }
	    else {
	      if (dAL1 > (pP1->dMaxParameter() / 2.0)) {
		pP2->vPointAtParameterRatio(0.01, adPt2);
	      }
	      else {
		pP2->vPointAtParameterRatio(0.99, adPt2);
	      }
	    }

	    // We have points close to the corner.
	    // Now get a vector going from point 1 to point 2
	    double adVec[2], adNorm[2];
	    adVec[0] = adPt2[0] - adPt1[0];
	    adVec[1] = adPt2[1] - adPt1[1];
	    // Get also the normal at point 1...
	    pP1->vUnitNormal(adPt1, adNorm);

	    // Find dot product of the two..
	    double dDot = adVec[0] * adNorm[0] + adVec[1] * adNorm[1];
	    if (iFuzzyComp(dDot, 0.0) > 0) {
	      //						vMessage(2, "NEW METHOD: Vert %d has a small angle\n", iV);
	      LqVertSmallAngle[iV] = true;
	    }
	    else {
	      LqVertSmallAngle[iV] = false;
	    }
	  }
	  else {
	    LqVertSmallAngle[iV] = true;
	  }

	  /*
	   // Now have the points. Check the distances from the other patch.
	   dDist1 = pP1->dDistFromBdry(adPt2);
	   dDist2 = pP2->dDistFromBdry(adPt1);

	   //				vMessage(2, "dDist1: %f, dDist2: %f\n", dDist1, dDist2);

	   if ( (dDist1 > 0.0) && (dDist2 > 0.0)) {
	   // Really a small angle
	   aqVertSmallAngle[iV] = true;
	   vMessage(2, "Vert %d has a small angle\n", iV);
	   } else {
	   aqVertSmallAngle[iV] = false;
	   }
	   */

	}
	else {
	  // Large angle
	  LqVertSmallAngle[iV] = false;
	} // if angle between the two is smaller than minangle

	iP++;

      } // while

    } // if it is a boundary point

  } // for every vert
}

void
Bdry2D::vOrderConnectivityInfo()
{
  // This will order the connectivity list using these criteria:
  //  1) Vertex index
  //  2) Angle of patch from that vertex
  // Also keep track of where the vertex index changes

  structVertConnect sVConnect;
  int iVertNow = 0;
  /*
   vMessage(2, "Here is the unordered list\n");
   for (int kk = 0; kk < iNumVertsConnectFound; kk++) {
   sVConnect = LsVertConnectivity[kk];
   vMessage(2, "iV: %d -- iP: %d -- Angle: %f\n",
   sVConnect.iVert, sVConnect.iPatch, sVConnect.dAngle);
   }
   */

  // sort the list now
  LsVertConnectivity.vOrderList();

  // Save the indices...
  // First one is 0..
  LiStartVert[0] = 0;
  int iNumP = -1;
  bool qHit = false;
  for (int i = 0; i < iNumVertsConnectFound; i++) {
    iNumP++;
    sVConnect = LsVertConnectivity[i];
    if (sVConnect.iVert != iVertNow) {
      // new one...
      qHit = true;
      LiNumPatchesToPoint[iVertNow] = iNumP;
      iNumP = 0;
      LiStartVert[sVConnect.iVert] = i;
      iVertNow = sVConnect.iVert;
    }
  }
  // The if takes care of the single circle case... what a hack!
  if (qHit)
    LiNumPatchesToPoint[iVertNow] = iNumP + 1;
  /*
   vMessage(2, "Here is the ordered list\n");
   for (int kk = 0; kk < iNumVertsConnectFound; kk++) {
   sVConnect = LsVertConnectivity[kk];
   vMessage(2, "iV: %d -- iP: %d -- Angle: %f\n",
   sVConnect.iVert, sVConnect.iPatch, sVConnect.dAngle);
   }
   */
}

void
Bdry2D::vSetNumPatchesConnectedToPoint(int iVertIndex, int iNPatch)
{
  assert((iVertIndex >= 0) && (iVertIndex < iNumPoints()));
  LiNumPatchesToPoint[iVertIndex] = iNPatch;
}

BdryPatch2D*
Bdry2D::pB2DPatchConnectedToPoint(int iVertIndex, int iPI)
{
  assert((iVertIndex >= 0) && (iVertIndex < iNumPoints()));
  assert(iPI < iNumPatchesConnectedToPoint(iVertIndex));
  return pBdry(LsVertConnectivity[LiStartVert[iVertIndex] + iPI].iPatch);
}

bool
Bdry2D::qIsBdryPoint(int iV) const
{
  return LqBdryPoint[iV];
}

void
Bdry2D::vSetBdryPoint(int iV, bool qBdry)
{
  LqBdryPoint[iV] = qBdry;
}

//  bool Bdry2D::
//  qIsPatchUsed (int iChk)
//    const {
//    return LsPatchUsed[iChk].qIsUsed;
//  }

//  int Bdry2D::
//  iPatchWasReplacedBy (int iChk)
//  {
//    return LsPatchUsed[iChk].iWasReplacedBy;
//  }

bool
BdryPatch2D::qDisjointFrom(const BdryPatch* const pBP) const
{
  // For 2D bdry patches, the only legal connections are at the ends, so
  // patches are disjoint if (and only if) they do not have a common
  // endpoint.  However, since we don't -know- that the BdryPatch* arg
  // is a BdryPatch2D, we have to use subterfuge here, with recursive
  // calls.
  switch (eGeom2DType)
    {
    case eCircle:
      return true;
    default:
      return (dynamic_cast<const BdryPatch2D*>(pBP)->qDisjointFrom(
	  iBeginPointIndex())
	  && dynamic_cast<const BdryPatch2D*>(pBP)->qDisjointFrom(
	      iEndPointIndex()));

    }
}

bool
BdryPatch2D::qDisjointFrom(const int iPt) const
{
  // If the vert is an endpoint, we're done.  Otherwise, compare its
  // surface distance to (approximately) zero.
  // Get the coords
  if (iPt == iBeginPointIndex() || iPt == iEndPointIndex())
    return false;
  double adC[2];
  pB2D->vGetPoint(iPt, adC);
  double dDist = dDistFromBdry(adC);
  // Can't exactly use iFuzzyComp(dDist, 0), because the splines don't
  // always pass this test.
  if (fabs(dDist) < 3.e-12)
    return false; //DOUG: trying to debug
  else
    return true;
}

bool
BdryPatch2D::qDisjointFrom(const Vert* const pV) const
{
  // No choice but to check it the hard way...
  double adB[2], adE[2];
  bool qB, qE;
  qB = false;
  qE = false;

  //	switch (eGeom2DType) {
  //		case eCircle :
  //			return true;
  //		default:
  if (eGeom2DType != eCircle) {
    vGetBeginPoint(adB);
    vGetEndPoint(adE);
    if (pV->isValid()) {
      qB = ((iFuzzyComp(adB[0], pV->x()) == 0)
	  && (iFuzzyComp(adB[1], pV->y()) == 0));
      qE = ((iFuzzyComp(adE[0], pV->x()) == 0)
	  && (iFuzzyComp(adE[1], pV->y()) == 0));
    }
  }

  if (qB || qE || !pV->isValid())
    return false;
  // Check the distance then...
  double dDist = dDistFromBdry(pV->getCoords());

  //if(dDist < 1.e-8)
  //  printf("distance = %e\n", dDist);

  // make sure we give some slack..
  if (fabs(dDist) < 1.e-6)
    return false;  //DOUG: trying to debug
  else
    return true;
//	}
}

double
BdryPatch2D::dGeodesicDistance(const double dParam1, const double dParam2) const
{
  // This routine computes the geodesic distance along a curve by
  // integrating sqrt((dx/dt)^2 + (dy/dt)^2) dt from point 1 to point
  // 2.  This is done numerically using Simpson's rule with three levels
  // of refinement, plus Richardson extrapolaton to work out the result
  // for an infinitely large number of integration intervals.

  double dDeltaT = dParam2 - dParam1;
  double adSpeeds[9];
  int i;
  for (i = 0; i < 9; i++) {
    adSpeeds[i] = dSpeed(dParam1 + (dDeltaT * i) / 8);
  }
  logMessage(4, "Speeds:\n");
  logMessage(4, "  %12.6g %12.6g %12.6g %12.6g %12.6g\n", adSpeeds[0],
	     adSpeeds[1], adSpeeds[2], adSpeeds[3], adSpeeds[4]);
  logMessage(4, "  %12.6g %12.6g %12.6g %12.6g\n", adSpeeds[5], adSpeeds[6],
	     adSpeeds[7], adSpeeds[8]);

  // Three point Simpson's rule:
  double dCoarseDist = (adSpeeds[0] + 4 * adSpeeds[4] + adSpeeds[8]) * dDeltaT
      / 6;

  // Five point (three-point twice):
  double dMediumDist = (adSpeeds[0] + 4 * adSpeeds[2] + 2 * adSpeeds[4]
      + 4 * adSpeeds[6] + adSpeeds[8]) * dDeltaT / 12;

  // Nine point (three-point four times:
  double dFineDist = (adSpeeds[0] + 4 * adSpeeds[1] + 2 * adSpeeds[2]
      + 4 * adSpeeds[3] + 2 * adSpeeds[4] + 4 * adSpeeds[5] + 2 * adSpeeds[6]
      + 4 * adSpeeds[7] + adSpeeds[8]) * dDeltaT / 24;
  logMessage(4, "%12.6g %12.6g %12.6g\n", dCoarseDist, dMediumDist, dFineDist);

  // Now for the Richardson extrapolation.  This assumes implicitly that
  // the error is converging like O(n^k) for some k.  In any event, the
  // result is likely to be more accurate than the fine estimate, so go
  // with it.
  if (iFuzzyComp(dFineDist, dMediumDist) == 0) {
    return (dFineDist);
  }
  else {
    double dDiff1 = dMediumDist - dFineDist;
    double dDiff2 = dCoarseDist - dMediumDist;
    double dRatio = dDiff2 / dDiff1;
    double dFineError = dDiff1 / (dRatio - 1);
    return (dFineDist - dFineError);
  }
}

void
BdryPatch2D::vGeodesicDistanceRatioBetweenTwoPoints(double dRatio,
						    const double adLoc1[2],
						    const double adLoc2[2],
						    double adPoint[2]) const
{
  double dParam1 = dParameterAtPoint(adLoc1);
  double dParam2 = dParameterAtPoint(adLoc2);
  double dTotalDistance = fabs(dGeodesicDistance(dParam1, dParam2));
  double dTargetDistance = dRatio * dTotalDistance;

  double dLowerParam = dParam1;
  double dUpperParam = dParam2;
  double dThisParam;

  while (true) {
    dThisParam = (dLowerParam + dUpperParam) / 2.;   // try the middle

    double dThisDistance = fabs(dGeodesicDistance(dParam1, dThisParam));

    if (fabs(dLowerParam - dUpperParam) < 1e-12) {
      assert(fabs(dThisDistance - dTargetDistance) < 1e-3);
      break;
    }
    else if (iFuzzyComp(dThisDistance, dTargetDistance) > 0) { // solution is in lower half
      dUpperParam = dThisParam;
    }
    else { // solution is in upper half
      dLowerParam = dThisParam;
    }
  }
  vPointAtParameter(dThisParam, adPoint);
}
;
