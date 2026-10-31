#include "GR_config.h"
#include "GR_assert.h"
#include "GR_InsertionQueueEntry.h"
#include "GR_InsertionQueue.h"
#include "GR_Mesh.h"
#include "GR_Quality.h"
#include "GR_Geometry.h"
#include "GR_QualMeasure.h"
#include "GR_InsertionManager.h"

#undef EXPERIMENTAL

// Major items:
//  o InsertionQueue creation and updates (this class is expected to be a
//    thin wrapper around the STL multimap container)
//  o InsertionQueueEntry stuff, most notably insertion for entries that
//    are legally allowed to insert their circumcenter
//  o A new driver for Shewchuk insertion based on the InsertionQueue
//    model.

// Set up an initial insertion queue depending on what type the mesh
// is.  This approach has the strong advantage that the initialization
// can be mesh-type aware, but the queue processing for insertion can be
// type-blind.
InsertionQueue::InsertionQueue(Mesh * const pMIn, const double dQualTarget,
			       const bool bdry_unencroached) :
    MMap(), dWorstAllowable(0), iMaxCellCount(INT_MAX), qIsInitialEncroachmentFixed(
	bdry_unencroached), qShape(true), qSize(true), m_prio(eWorstQuality), pM(
	NULL)

#ifdef HAVE_MESQUITE
	, m_eIsoQM(GRUMMP::TMOPQual::eDefault), m_eAnisoQM(
	GRUMMP::TMOPQual::eDefault), dMinMLength(0), dMaxMLength(0), dMTol(0)
#endif

{
  logMessage(4, "Initializing insertion priority queue...\n");
#ifdef EXPERIMENTAL
  vWarning("  Experimental insertion priority code in use!\n");
#endif
  assert(pMIn != NULL);
  pM = pMIn;

  if (iFuzzyComp(dQualTarget, InsertionQueueEntry::dInvalidPriority) == 0) {
    dWorstAllowable = pM->getWorstAllowedCellShape();
  }
  else {
    dWorstAllowable = dQualTarget;
  }
  logMessage(4, "Worst allowed cell shape measure is %f.\n", dWorstAllowable);

//  vBuildQueue(bdry_unencroached);
}

//// Lets me pick my own cells to attempt to refine.
//InsertionQueue::InsertionQueue(Mesh * const pMIn, const double dQualTarget,
//		const bool bdry_unencroached, std::set<Cell*> &cellQueue) :
//		MMap(), dWorstAllowable(0), iMaxCellCount(INT_MAX), qIsInitialEncroachmentFixed(
//				bdry_unencroached), qShape(true), qSize(true), m_prio(
//				eWorstQuality), pM(NULL) {
//	logMessage(1, "Initializing insertion priority queue...\n");
//#ifdef EXPERIMENTAL
//	vWarning("  Experimental insertion priority code in use!\n");
//#endif
//	assert(pMIn != NULL);
//	pM = pMIn;
//
//	if (iFuzzyComp(dQualTarget, InsertionQueueEntry::dInvalidPriority) == 0) {
//		dWorstAllowable = pM->getWorstAllowedCellShape();
//	} else {
//		dWorstAllowable = dQualTarget;
//	}
//	logMessage(2, "Worst allowed cell shape measure is %f.\n", dWorstAllowable);
//
//	MMap.clear();
//	logMessage(2,
//			"Queueing simplex cells that don't meet quality criterion...");
//
//	std::set<Cell*>::iterator iter = cellQueue.begin(), iterEnd =
//			cellQueue.end();
//
//	for (; iter != iterEnd; ++iter) {
//		Cell *pC = *iter;
//		if (!pC->isDeleted() && dynamic_cast<SimplexCell*>(pC) != NULL) {
//			qAddEntry(pC);
//		}
//	}
//
//	logMessage(2, "queued %d.\n", int(MMap.size()));
//
//	// Now, identify all encroached bdry faces and add them to the queue.
//	if (!bdry_unencroached && pM->areBdryChangesAllowed())
//		pM->queueEncroachedBdryEntities(*this);
//	logMessage(2, "Total queue length is %u.\n", iQueueLength());
//}

InsertionQueueEntry::InsertionQueueEntry(Cell * const pCIn) :
    dPri(dInvalidPriority), pC(pCIn), eET(eInvalid), filterPT(0)
{
  assert(pC->isValid());
  switch (pC->getType())
    {
    case Cell::eTet:
      eET = eTetCell;
      break;
    case Cell::eTriCell:
      eET = eTriCell;
      break;
    default:
      // Shouldn't get here, but we might for boundary faces.
      assert(0);
      break;
    }
  if (pC->isDeleted()) {
    return;
    // FIX ME: New bdry curve insertion will make this unnecessary,
    // because no deleted cells can ever be queued.
  }
  for (int i = 0; i < pC->getNumVerts(); i++) {
    apV[i] = pC->getVert(i);
  }
  logMessage(4, "Queueing a bad tet: %p\n", pC);
}

InsertionQueueEntry::InsertionQueueEntry(BFace * const pBFIn) :
    dPri(dInvalidPriority), pBF(pBFIn), eET(eInvalid), filterPT(0)
{
  assert(pBF->isValid());
  switch (pBF->getType())
    {
    case Cell::eTriBFace:
    case Cell::eIntTriBFace:
      eET = eTriBFace;
      break;
    case Cell::eBdryEdge:
    case Cell::eIntBdryEdge:
      eET = eBdryEdge;
      break;
    default:
      // Shouldn't get here
      assert(0);
      break;
    }
  if (pBF->isDeleted()) {
    return;
    // FIX ME: New bdry curve insertion will make this unnecessary,
    // because no deleted BFaces can ever be queued.
  }
  for (int i = 0; i < pBF->getNumVerts(); i++) {
    apV[i] = pBF->getVert(i);
  }
}

InsertionQueueEntry::InsertionQueueEntry(Vert * const pV0, Vert * const pV1) :
    dPri(dInvalidPriority), pC(pCInvalidCell), filterPT(0)
{

  if ((pV0->getVertType() == Vert::eBdry
      || pV0->getVertType() == Vert::eBdryApex
      || pV0->getVertType() == Vert::eBdryCurve)
      && (pV1->getVertType() == Vert::eBdry
	  || pV1->getVertType() == Vert::eBdryApex
	  || pV1->getVertType() == Vert::eBdryCurve)) {
    eET = eSubSeg;
  }
  assert(pV0->isValid() && pV1->isValid());
  apV[0] = pV0;
  apV[1] = pV1;
}

//Next two for Aniso Queue Entry
InsertionQueueEntry::InsertionQueueEntry(Face * const pFace) :
    dPri(dInvalidPriority), pF(pFace), filterPT(0)
{
  switch (pFace->getNumVerts())
    {
    case 2:
      eET = e2DFace;
      assert(pF->getVert(0)->isValid() && pF->getVert(1)->isValid());
      assert(pFace->doFullCheck());
      apV[0] = pF->getVert(0);
      apV[1] = pF->getVert(1);
      assert(pFace->getLeftCell()->doFullCheck());
      assert(pFace->getRightCell()->doFullCheck());
      break;

    case 3:
      eET = e3DFace;
      apV[0] = pF->getVert(0);
      apV[1] = pF->getVert(1);
      apV[2] = pF->getVert(2);
      assert(pFace->getLeftCell()->doFullCheck());
      assert(pFace->getRightCell()->doFullCheck());

      break;
    }
}

bool
InsertionQueueEntry::qStillInMesh() const
// Determine whether the entity described by this entry is still in
// the mesh in its original form.
{
  switch (eET)
    {
    case eTriCell:
    case eTetCell:
      {
	// If it's been deleted from the mesh, check no farther.
	if (pC->isDeleted()) {
	  return false;
	}
	bool qRetVal = true;
	for (int i = 0; i < pC->getNumVerts() && qRetVal; i++) {
	  qRetVal = qRetVal && pC->hasVert(apV[i]);
	}
	return qRetVal;
      }
    case eTriBFace:
    case eBdryEdge:
      // Like the cell case.
      {
	// If it's been deleted from the mesh, check no farther.
	if (pBF->isDeleted()) {
	  return false;
	}
	bool qRetVal = true;
	for (int i = 0; i < pBF->getNumVerts() && qRetVal; i++) {
	  qRetVal = qRetVal && pBF->hasVert(apV[i]);
	}
	return qRetVal;
      }

    case e3DFace:
      Face *pFt;
      pFt = findCommonFace(apV[0], apV[1], apV[2], NULL, true);
      if (pFt == pF)
	return true;
      else
	return false;

    case eSubSeg:
    case e2DFace:

      // A subsegment can only be removed by splitting that
      // subsegment, not some other one, so this is always false unless
      // the subseg got into the queue twice (possible).
      {
	if (!apV[0] || apV[0] == NULL || apV[0]->isDeleted()
	    || !apV[0]->isValid())
	  return false;
	std::set<Cell*> spCJunk;
	std::set<Vert*> spVNeigh;
	findNeighborhoodInfo(apV[0], spCJunk, spVNeigh);
	if (spVNeigh.find(apV[1]) != spVNeigh.end())
	  return true;
	else
	  return false;
      }

    default:
      // Should never get here.
      assert(0);
      return false;
    }
}

double
InsertionQueue::dEvaluatePriority(const InsertionQueueEntry IQE) const
{
  InsertionQueueEntry::eEntryType eET = IQE.eType();
  // Compute the priority the first time.  After that, return a stored value.
  double dPri = 0;    //IQE.dPriority();
//  if (iFuzzyComp(dPri, IQE.dInvalidPriority) == 0) {
  if (1) {

    switch (eET)
      {

      case InsertionQueueEntry::eSubSeg:
	dPri = IQE.dSubSegPriority;
	break;

      case InsertionQueueEntry::eTriBFace:
      case InsertionQueueEntry::eBdryEdge:
	dPri = IQE.pBFBFace()->calcInsertionPriority();
	break;

      case InsertionQueueEntry::eTetCell:
      case InsertionQueueEntry::eTriCell:
	{
	  // -D/LS (if bigger than 1) + min(angle, shape), where shape is
	  // what's shown below and angle is min dihed (if min dihed is
	  // small enough for us to care; otherwise, angle = 1)
	  Cell * pC = IQE.pCCell();
	  double dSizeQual = 0, dShapeQual = 0;
	  int iNVerts = pC->getNumVerts();
	  int iDim = pC->getVert(0)->getSpaceDimen();

	  double dRadius = (dynamic_cast<SimplexCell*>(pC))->calcCircumradius();
	  double dSize = 2 / sqrt(double(iDim)) * dRadius;

	  // Find the correct size based on length scale.  Also, if the
	  // vertex is at a small angle in the input, make sure that you
	  // don't try to split this cell.

	  double dCorrectSize = 0;

	  for (int ii = 0; ii < iNVerts; ii++) {

	    Vert *pV = pC->getVert(ii);
	    double length_scale = pV->getLengthScale();
	    assert(iFuzzyComp(length_scale, 0.) == 1);

	    dCorrectSize += length_scale;
	    //double dThisSize = pV->dLS();
	    //dCorrectSize = MIN(dCorrectSize, dThisSize);

	    // Now, only want to do this with Tets.
//	  if (pV->isSmallAngleVert() && eET != InsertionQueueEntry::eTriCell)
//	    dShapeQual = IQE.dNeverQueuedPriority;

	  }

	  //Hack to enforce Miller et al. small angle
	  //condition on triangle meshes. SG Sept 13 2006.

//	if(eET == InsertionQueueEntry::eTriCell) {

	  TriCell* triangle = dynamic_cast<TriCell*>(pC);
	  assert(triangle);

	  if (triangle->isSmallAngleAffected()) {
	    dShapeQual = IQE.dNeverQueuedPriority;
	    dSizeQual = IQE.dNeverQueuedPriority;
	  }

//	}

	  dCorrectSize /= iNVerts;
	  logMessage(
	      4, "Correct size: %10.6f  Actual size: %10.6f  Ratio:  %10.7f\n",
	      dCorrectSize, dSize, dSize / dCorrectSize);
	  if (dSize > dCorrectSize && qSize) {
	    dSizeQual = -dSize / dCorrectSize;
	  }

	  else {

#ifndef EXPERIMENTAL
	    // This call works properly in both 2D and 3D.
	    double dShort = pC->calcShortestEdgeLength();
	    double dShortestToRadius;

	    // The following scale factor is used to normalize the
	    // shortest-edge-to-circumradius ratio to the range 0 to 1.

	    if (iDim == 2)
	      dShortestToRadius = dShort / (sqrt(3.) * dRadius);
	    else
	      dShortestToRadius = dShort * sqrt(.375) / dRadius;
	    if (dShortestToRadius < dWorstAllowable) {
//	  	if(m_prio == eWorstQuality)
//	  		dShapeQual += dShortestToRadius;
//	  	else if(m_prio == eShortestEdge)
	      dShapeQual += dShort;
//	  	else if(m_prio == eLongestShortest)
//	  		dShapeQual += 1./dShort;

	    }
	    else {
	      dShapeQual += IQE.dNeverQueuedPriority;
	    }
#else
	    // 2D same as ever.
	    if (iDim == 2 && qShape) {
	      double dShort = pC->calcShortestEdgeLength();
	      dShapeQual += dShort / (sqrt(3.)*dRadius);
	    }

	    else {
	      double dSliverQual;
	      int ii;
	      vInToCircumSphere(pC, &ii, &dSliverQual);
	      logMessage(4, "This tet has in- to circum- ratio of %f\n",
		  dSliverQual);
	      dShapeQual += dSliverQual;
	    }
#endif
	  }

	  // The following code for including dihedral angle in choosing
	  // cells to split turns out not to be very effective in practice
	  // in eliminating bad tetrahedra.
//	if (dPri > 0.2 && eET == eTetCell) {
//	  // If the cell is a tet and the aspect ratio is reasonably
//	  // good, consider dihedral angle as well.
//	  int iN;
//	  double adRes[1];
//	  vMinDihed(pC, &iN, adRes);
//	  assert(iN == 1);
//	  // Convert from degrees back to radians.
//	  double dAnggleQual = adRes[0] * M_PI/180;
//	  dShapeQual = MIN(dShapeQual, dAngleQual*2);
//	}

	  //printf("size qual = %lf, shape qual = %lf \n", dSizeQual, dShapeQual);
	  if (qSize && qShape)
	    dPri = dSizeQual + dShapeQual;
	  else if (qSize)
	    dPri = dSizeQual;
	  else if (qShape)
	    dPri = dShapeQual;

	  if (dPri < IQE.dWorstCellPriority)
	    dPri = IQE.dWorstCellPriority;
	  logMessage(4, "Priority of a cell is %f (%f + %f), worst angle %f\n",
		     dPri, dSizeQual, dShapeQual, pC->calcMinSolid());
	}
	break;
      case InsertionQueueEntry::eInvalid:
	assert(0);
	break;

      default:
	// Should never get here.
	break;
      }
  }

  IQE.setPriority(dPri);
  return (dPri);
}

#ifdef HAVE_MESQUITE
GRUMMP::TMOPQual::eQMetric
InsertionQueue::SetQMetricType(int iNMetricValues, double dTolerance,
			       std::set<Vert*> setpV)
{

  GRUMMP::TMOPQual::eQMetric QMetric;
  double adMetric[iNMetricValues];
  std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();

  for (int i = 0; i < iNMetricValues; i++)
    adMetric[i] = 0;

  for (; iter != iterEnd; ++iter) {
    Vert* psetV = (*iter);
    for (int i = 0; i < iNMetricValues; i++)
      adMetric[i] += (1. / setpV.size()) * psetV->getMetric(i);
  }

  QMetric = m_eAnisoQM;

  if (iNMetricValues == 3) {
    if (adMetric[1] < dTolerance && adMetric[1] > -dTolerance) {
      if (fabs(adMetric[0] - adMetric[2]) < dTolerance) {
	QMetric = m_eIsoQM;
      }
    }
  }

  if (iNMetricValues == 6) {
    if (adMetric[1] < dTolerance && adMetric[1] > -dTolerance
	&& adMetric[3] < dTolerance && adMetric[3] > -dTolerance
	&& adMetric[4] < dTolerance && adMetric[4] > -dTolerance) {
      if (fabs(adMetric[0] - adMetric[2]) < dTolerance
	  && fabs(adMetric[0] - adMetric[5]) < dTolerance
	  && fabs(adMetric[2] - adMetric[5]) < dTolerance) {
	QMetric = m_eIsoQM;
      }
    }
  }

  return QMetric;
}

bool
InsertionQueue::qAddEntryAnisoEdge(InsertionQueueEntry IQE)
{

  double dPriority = dEvaluateAnisoEdgePriority(IQE);

  if (dPriority < 0) {
    MMap.insert(std::pair<const double, InsertionQueueEntry>(dPriority, IQE));
    return true;
  }
  return false;
}

double
InsertionQueue::d2DFaceBestImprove(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual::eQMetric QM;
  GRUMMP::TMOPQual TQ;

  Face *pF = IQE.pFFace();
  Cell*pCL = pF->getLeftCell();
  Cell*pCR = pF->getRightCell();

  assert(pCL->doFullCheck());
  assert(pCR->doFullCheck());
  assert(pF->doFullCheck());

  Vert *pVA = pF->getVert(0);
  Vert *pVB = pF->getVert(1);
  Vert *pVC = pCL->getOppositeVert(pF);
  Vert *pVD = pCR->getOppositeVert(pF);

  double QualBef[2];
  double QualAft[4];
  double dIP[2];
  double dT[9];
  double adMetric[3];

  std::set<Vert*> setpV;
  setpV.insert(pF->getVert(0));
  setpV.insert(pF->getVert(1));
  QM = SetQMetricType(3, dMTol, setpV);

  dIP[0] = (pVA->x() + pVB->x()) / 2.;
  dIP[1] = (pVA->y() + pVB->y()) / 2.;

  TQ.LogEuclideanInterpMetric(dIP, adMetric, setpV);

  TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
  QualBef[0] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			pVB->getCoords(), pVC->getCoords());

  TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVA, pVC);
  QualAft[0] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(), dIP,
			pVC->getCoords());

  TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVB, pVC);
  QualAft[1] = TQ.dEval(QM, CellSkel::eTriCell, dT, dIP, pVB->getCoords(),
			pVC->getCoords());

  TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVD, pVB);
  QualBef[1] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			pVD->getCoords(), pVB->getCoords());

  TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVA, pVD);
  QualAft[2] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
			pVD->getCoords(), dIP);

  TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVD, pVB);
  QualAft[3] = TQ.dEval(QM, CellSkel::eTriCell, dT, dIP, pVD->getCoords(),
			pVB->getCoords());

  double dQualMaxAft;
  dQualMaxAft = std::max(std::max(QualAft[0], QualAft[1]),
			 std::max(QualAft[2], QualAft[3]));
  dPri = dQualMaxAft < std::max(QualBef[0], QualBef[1]) ? -dQualMaxAft : 1;

  double dML = 1e10;
  double dMinML = 1e10;
  double adM[3];
  double adC[2];
  const Vert *apV[4] =
    { pVA, pVB, pVC, pVD };

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++)
      adM[j] = apV[i]->getMetric(j);

    adC[0] = apV[i]->x();
    adC[1] = apV[i]->y();

    dML = TQ.dLogMetricLength(dIP, adMetric, adC, adM);
    dMinML = std::min(dML, dMinML);
  }

  dPri = dMinML < dMinMLength ? 1 : dPri;

  return dPri;
}

double
InsertionQueue::d2DFaceLongestEdge(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual TQ;

  Face *pF = IQE.pFFace();
  Vert *pVA = pF->getVert(0);
  Vert *pVB = pF->getVert(1);
  Vert *pVC = pF->getLeftCell()->getOppositeVert(pF);
  Vert *pVD = pF->getRightCell()->getOppositeVert(pF);

  double dMLength;
  dMLength = TQ.dLogMetricLength(pF);
  dPri = dMaxMLength - dMLength;

  double dML = 1e10;
  double dMinML = 1e10;
  double adMetric[3];
  double adM[3];
  double adC[2];
  double dIP[2];
  const Vert *apV[4] =
    { pVA, pVB, pVC, pVD };

  dIP[0] = (pVA->x() + pVB->x()) / 2.;
  dIP[1] = (pVA->y() + pVB->y()) / 2.;

  std::set<Vert*> setpV;
  for (int i = 0; i < pF->getNumVerts(); i++)
    setpV.insert(pF->getVert(i));

  TQ.LogEuclideanInterpMetric(dIP, adMetric, setpV);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++)
      adM[j] = apV[i]->getMetric(j);

    adC[0] = apV[i]->x();
    adC[1] = apV[i]->y();

    dML = TQ.dLogMetricLength(dIP, adMetric, adC, adM);
    dMinML = std::min(dML, dMinML);
  }
  dPri = dMinML < dMinMLength ? 1 : dPri;

  return dPri;
}

double
InsertionQueue::d2DFaceShortestEdge(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual TQ;

  Face *pF = IQE.pFFace();
  Vert* pV0 = pF->getVert(0);
  Vert* pV1 = pF->getVert(1);
  double dMLength;

  dMLength = TQ.dLogMetricLength(pF);
  dPri = dMLength - dMinMLength;

  std::set<Vert*>::iterator iter, iterE;
  std::set<Cell*>::iterator iterC, iterEC;
  std::set<Cell*> spC0;
  std::set<Cell*> spC1;
  std::set<Vert*> spV0;
  std::set<Vert*> spV1;

  findNeighborhoodInfo(pV0, spC0, spV0);

  iter = spV0.begin();
  iterE = spV0.end();
  iterC = spC0.begin();
  iterEC = spC0.end();

  bool NFMinimum = false;
  int icount = 0;

  for (; iter != iterE; ++iter) {
    Vert* psV = (*iter);
    Face* pFN = findCommonFace(pV0, psV, true);
    if (dMinMLength > TQ.dLogMetricLength(pFN)) {
      icount++;
    }
  }

  if (icount >= 2) {
    if (pV0->getVertType() != Vert::eBdry) {
      pV0->markToDelete();
      NFMinimum = true;
    }
  }

  icount = 0;

  findNeighborhoodInfo(pV1, spC1, spV1);

  iter = spV1.begin();
  iterE = spV1.end();
  iterC = spC1.begin();
  iterEC = spC1.end();

  for (; iter != iterE; ++iter) {
    Vert* psV = (*iter);
    Face* pFN = findCommonFace(pV1, psV, true);
    if (dMinMLength > TQ.dLogMetricLength(pFN)) {
      icount++;
    }
  }

  if (icount >= 2) {
    if (pV1->getVertType() != Vert::eBdry) {
      pV1->markToDelete();
      NFMinimum = true;
    }
  }

  if (NFMinimum == false) {
    dPri = 1;
    pV0->markToKeep();
    pV1->markToKeep();
  }

  return dPri;
}

double
InsertionQueue::d2DBdryBestImprove(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual::eQMetric QM;
  GRUMMP::TMOPQual TQ;

  Face *pF = IQE.pBFBFace()->getFace(0);

  Vert *pVA;
  Vert *pVB;
  Vert *pVC;
  Cell *pC = pF->getOppositeCell(IQE.pBFBFace());

  assert(pF->doFullCheck());

  pVA = pF->getVert(0);
  pVB = pF->getVert(1);
  pVC = pC->getOppositeVert(pF);

  assert(pC->isValid());
  assert(pVA->isValid());
  assert(pVB->isValid());
  assert(pVC->isValid());

  double adMetric[3];
  double QualBef;
  double QualAft[2];
  double dIP[2];
  double dT[9];

  std::set<Vert*> setpV;
  for (int i = 0; i < pF->getNumVerts(); i++)
    setpV.insert(pF->getVert(i));

  QM = SetQMetricType(3, dMTol, setpV);

  dIP[0] = (pVA->x() + pVB->x()) / 2.;
  dIP[1] = (pVA->y() + pVB->y()) / 2.;

  TQ.LogEuclideanInterpMetric(dIP, adMetric, setpV);

  TQ.SetTargetfromMetric(dT, CellSkel::eTriCell, pVA, pVB, pVC);
  QualBef = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(),
		     pVB->getCoords(), pVC->getCoords());

  TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVA, pVC);
  QualAft[0] = TQ.dEval(QM, CellSkel::eTriCell, dT, pVA->getCoords(), dIP,
			pVC->getCoords());

  TQ.SetTargetfromMetricNewVertice(dT, CellSkel::eTriCell, adMetric, pVB, pVC);
  QualAft[1] = TQ.dEval(QM, CellSkel::eTriCell, dT, dIP, pVB->getCoords(),
			pVC->getCoords());

  dPri =
      std::max(QualAft[0], QualAft[1]) < QualBef ?
	  -std::max(QualAft[0], QualAft[1]) : 1;

  return dPri;
}

double
InsertionQueue::d2DBdryLongestEdge(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual TQ;

  Face *pF = IQE.pBFBFace()->getFace(0);
  Vert *pVA = pF->getVert(0);
  Vert *pVB = pF->getVert(1);

  double dMLength;
  dMLength = TQ.dLogMetricLength(pF);
  dPri = dMaxMLength - dMLength;

  double dML = 1e10;
  double dMinML = 1e10;
  double adMetric[3];
  double adM[3];
  double adC[2];
  double dIP[2];
  const Vert *apV[2] =
    { pVA, pVB };

  dIP[0] = (pVA->x() + pVB->x()) / 2.;
  dIP[1] = (pVA->y() + pVB->y()) / 2.;

  std::set<Vert*> setpV;
  for (int i = 0; i < pF->getNumVerts(); i++)
    setpV.insert(pF->getVert(i));

  TQ.LogEuclideanInterpMetric(dIP, adMetric, setpV);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++)
      adM[j] = apV[i]->getMetric(j);

    adC[0] = apV[i]->x();
    adC[1] = apV[i]->y();

    dML = TQ.dLogMetricLength(dIP, adMetric, adC, adM);
    dMinML = std::min(dML, dMinML);
  }
  dPri = dMinML < dMinMLength ? 1 : dPri;

  return dPri;
}

double
InsertionQueue::d2DBdryShortestEdge(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual TQ;

  Face *pF = IQE.pBFBFace()->getFace(0);
  double dMLength;
  dMLength = TQ.dLogMetricLength(pF);
  dPri = dMLength - dMinMLength;
  if (dPri < 0) {
    pF->getVert(0)->markToDelete();
    pF->getVert(1)->markToDelete();
  }

  return dPri;
}

double
InsertionQueue::d3DSubSegBestImprove(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual::eQMetric QM;
  GRUMMP::TMOPQual TQ;
  TQ.SetMetricValue(6);

  bool bMinLenght = false;

  Vert *pV0 = IQE.pVVert(0);
  Vert *pV1 = IQE.pVVert(1);

  int iNCells = 0;

  std::set<Cell*> spCV0;
  std::set<Cell*> spC;
  std::set<Vert*> spV0Neigh;
  std::set<Cell*>::iterator iterC, iterCEnd;

  findNeighborhoodInfo(pV0, spCV0, spV0Neigh);

  iterC = spCV0.begin();
  iterCEnd = spCV0.end();

  int iNCount = 0;
  for (; iterC != iterCEnd; ++iterC) {
    Cell* psetC = (*iterC);
    Vert *pVV[4];
    for (int i = 0; i < 4; i++)
      pVV[i] = psetC->getVert(i);

    if ((pVV[0] == pV0 || pVV[1] == pV0 || pVV[2] == pV0 || pVV[3] == pV0)
	&& (pVV[0] == pV1 || pVV[1] == pV1 || pVV[2] == pV1 || pVV[3] == pV1)) {
      spC.insert(psetC);
      iNCells++;
    }
  }

  double dQualBef[iNCells];
  double dQualAft1[iNCells];
  double dQualAft2[iNCells];
  const double dIP[3] =
    { (pV0->x() + pV1->x()) / 2., (pV0->y() + pV1->y()) / 2., (pV0->z()
	+ pV1->z()) / 2. };

  iterC = spC.begin();
  iterCEnd = spC.end();

  for (; iterC != iterCEnd; ++iterC) {
    Cell* psetC = (*iterC);

    Vert *pVA;
    Vert *pVB;
    Vert *pVC;
    Vert *pVD;

    if (psetC->getVert(0) == pV0 && psetC->getVert(1) == pV1) {
      pVA = psetC->getVert(0);
      pVB = psetC->getVert(1);
      pVC = psetC->getVert(2);
      pVD = psetC->getVert(3);
    }

    else if (psetC->getVert(0) == pV0 && psetC->getVert(2) == pV1) {
      pVA = psetC->getVert(0);
      pVB = psetC->getVert(2);
      pVC = psetC->getVert(1);
      pVD = psetC->getVert(3);
    }

    else if (psetC->getVert(0) == pV0 && psetC->getVert(3) == pV1) {
      pVA = psetC->getVert(0);
      pVB = psetC->getVert(3);
      pVC = psetC->getVert(1);
      pVD = psetC->getVert(2);
    }

    else if (psetC->getVert(1) == pV0 && psetC->getVert(0) == pV1) {
      pVA = psetC->getVert(1);
      pVB = psetC->getVert(0);
      pVC = psetC->getVert(2);
      pVD = psetC->getVert(3);
    }

    else if (psetC->getVert(1) == pV0 && psetC->getVert(2) == pV1) {
      pVA = psetC->getVert(1);
      pVB = psetC->getVert(2);
      pVC = psetC->getVert(0);
      pVD = psetC->getVert(3);
    }

    else if (psetC->getVert(1) == pV0 && psetC->getVert(3) == pV1) {
      pVA = psetC->getVert(1);
      pVB = psetC->getVert(3);
      pVC = psetC->getVert(0);
      pVD = psetC->getVert(2);
    }

    else if (psetC->getVert(2) == pV0 && psetC->getVert(0) == pV1) {
      pVA = psetC->getVert(2);
      pVB = psetC->getVert(0);
      pVC = psetC->getVert(1);
      pVD = psetC->getVert(3);
    }

    else if (psetC->getVert(2) == pV0 && psetC->getVert(1) == pV1) {
      pVA = psetC->getVert(2);
      pVB = psetC->getVert(1);
      pVC = psetC->getVert(0);
      pVD = psetC->getVert(3);
    }

    else if (psetC->getVert(2) == pV0 && psetC->getVert(3) == pV1) {
      pVA = psetC->getVert(2);
      pVB = psetC->getVert(3);
      pVC = psetC->getVert(1);
      pVD = psetC->getVert(0);
    }

    else if (psetC->getVert(3) == pV0 && psetC->getVert(0) == pV1) {
      pVA = psetC->getVert(3);
      pVB = psetC->getVert(0);
      pVC = psetC->getVert(1);
      pVD = psetC->getVert(2);
    }

    else if (psetC->getVert(3) == pV0 && psetC->getVert(1) == pV1) {
      pVA = psetC->getVert(3);
      pVB = psetC->getVert(1);
      pVC = psetC->getVert(0);
      pVD = psetC->getVert(2);
    }

    else {
      assert(psetC->getVert(3) == pV0 && psetC->getVert(2) == pV1);
      pVA = psetC->getVert(3);
      pVB = psetC->getVert(2);
      pVC = psetC->getVert(0);
      pVD = psetC->getVert(1);
    }

    double dT[12];

    std::set<Vert*> setpV;
    for (int i = 0; i < psetC->getNumVerts(); i++)
      setpV.insert(psetC->getVert(i));
    QM = SetQMetricType(6, dMTol, setpV);

    if (checkOrient3D(pVA, pVB, pVC, pVD) == 1) {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVB, pVC, pVD);

      dQualBef[iNCount] = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(),
				   pVB->getCoords(), pVC->getCoords(),
				   pVD->getCoords());
      dQualAft1[iNCount] = TQ.dEval(QM, CellSkel::eTet, dT, dIP,
				    pVB->getCoords(), pVC->getCoords(),
				    pVD->getCoords());
      dQualAft2[iNCount] = TQ.dEval(QM, CellSkel::eTet, dT, dIP,
				    pVC->getCoords(), pVA->getCoords(),
				    pVD->getCoords());
    }
    else {
      TQ.SetTargetfromMetric(dT, CellSkel::eTet, pVA, pVC, pVB, pVD);

      dQualBef[iNCount] = TQ.dEval(QM, CellSkel::eTet, dT, pVA->getCoords(),
				   pVC->getCoords(), pVB->getCoords(),
				   pVD->getCoords());
      dQualAft1[iNCount] = TQ.dEval(QM, CellSkel::eTet, dT, dIP,
				    pVA->getCoords(), pVC->getCoords(),
				    pVD->getCoords());
      dQualAft2[iNCount] = TQ.dEval(QM, CellSkel::eTet, dT, dIP,
				    pVC->getCoords(), pVB->getCoords(),
				    pVD->getCoords());
    }
    iNCount++;

    double dMLs;
    double dMLsnew;
    double adMid[6];
    double adMTest[6];
    int iNM = 6;

    for (int i = 0; i < iNM; i++) {
      adMid[i] = (pVA->getMetric(i) + pVB->getMetric(i) + pVC->getMetric(i)
	  + pVD->getMetric(i)) / 4.;
    }

    for (int i = 0; i < iNM; i++)
      adMTest[i] = (pVA->getMetric(i) + adMid[i]) / 2.;

    dMLs = TQ.MLength3D(pVA->getCoords(), dIP, adMTest);

    for (int i = 0; i < iNM; i++)
      adMTest[i] = (pVB->getMetric(i) + adMid[i]) / 2.;

    dMLsnew = TQ.MLength3D(pVB->getCoords(), dIP, adMTest);
    dMLs = std::min(dMLs, dMLsnew);

    for (int i = 0; i < iNM; i++)
      adMTest[i] = (pVC->getMetric(i) + adMid[i]) / 2.;

    dMLsnew = TQ.MLength3D(pVC->getCoords(), dIP, adMTest);
    dMLs = std::min(dMLs, dMLsnew);

    for (int i = 0; i < iNM; i++)
      adMTest[i] = (pVD->getMetric(i) + adMid[i]) / 2.;

    dMLsnew = TQ.MLength3D(pVD->getCoords(), dIP, adMTest);
    dMLs = std::min(dMLs, dMLsnew);

    if (dMLs < dMinMLength)
      bMinLenght = true;
  }

  double dMaxAft = 0;
  double dMaxAft1 = 0;
  double dMaxAft2 = 0;
  double dMaxBef = 0;

  int iCountBetter = 0;
  int iCountWorst = 0;

  for (int i = 0; i < iNCount; i++) {
    dMaxAft1 = std::max(dQualAft1[i], dMaxAft1);
    dMaxAft2 = std::max(dQualAft2[i], dMaxAft2);
    dMaxAft = std::max(dMaxAft1, dMaxAft2);
    dMaxBef = std::max(dQualBef[i], dMaxBef);

    if ((dQualAft1[i] < dQualBef[i]) && (dQualAft2[i] < dQualBef[i]))
      iCountBetter++;
    else
      iCountWorst++;
  }

  dPri = iCountWorst < iCountBetter ? -dMaxAft : 1;
  if (bMinLenght == true)
    dPri = 1;

  return dPri;
}

double
InsertionQueue::d3DSubSegLongestEdge(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual TQ;
  TQ.SetMetricValue(6);

  Vert *pV0 = IQE.pVVert(0);
  Vert *pV1 = IQE.pVVert(1);
  double dMLength;
  int iNMetricValues = 6;
  double adMet[iNMetricValues];

  for (int i = 0; i < iNMetricValues; i++) {
    adMet[i] = pV0->getMetric(i) + pV1->getMetric(i);
    adMet[i] /= 2.0;
  }

  dMLength = TQ.MLength3D(pV0->getCoords(), pV1->getCoords(), adMet);
  dPri = dMaxMLength - dMLength;

  std::set<Cell*> spCV0;
  std::set<Cell*> spC;
  std::set<Vert*> spV0Neigh;
  std::set<Cell*>::iterator iterC, iterCEnd;

  findNeighborhoodInfo(pV0, spCV0, spV0Neigh);

  iterC = spCV0.begin();
  iterCEnd = spCV0.end();

  for (; iterC != iterCEnd; ++iterC) {
    Cell* psetC = (*iterC);
    Vert *pVV[4];
    for (int i = 0; i < 4; i++) {
      pVV[i] = psetC->getVert(i);
    }

    if ((pVV[0] == pV0 || pVV[1] == pV0 || pVV[2] == pV0 || pVV[3] == pV0)
	&& (pVV[0] == pV1 || pVV[1] == pV1 || pVV[2] == pV1 || pVV[3] == pV1)) {
      spC.insert(psetC);
    }
  }

  const double dNewCoords[3] =
    { (pV0->x() + pV1->x()) / 2., (pV0->y() + pV1->y()) / 2., (pV0->z()
	+ pV1->z()) / 2. };
  double dNewMetric[6];
  double dMLsnew;
  double dMLs = 1e30;

  for (int i = 0; i < 6; i++) {
    dNewMetric[i] = (pV0->getMetric(i) + pV1->getMetric(i)) / 2.;
  }

  iterC = spC.begin();
  iterCEnd = spC.end();

  for (; iterC != iterCEnd; ++iterC) {
    double adMTest[iNMetricValues];

    Cell* psetC = (*iterC);
    for (int i = 0; i < 4; i++) {

      Vert *pVTest = psetC->getVert(i);
      for (int ii = 0; ii < iNMetricValues; ii++) {
	adMTest[ii] = (pVTest->getMetric(ii) + dNewMetric[ii]) / 2.;
      }

      dMLsnew = TQ.MLength3D(pVTest->getCoords(), dNewCoords, adMTest);
      dMLs = std::min(dMLs, dMLsnew);
    }
  }

  if (dMLs < dMinMLength)
    dPri = 1;

  return dPri;
}

double
InsertionQueue::d3DSubSegShortestEdge(const InsertionQueueEntry& IQE)
{

  double dPri = 0;
  GRUMMP::TMOPQual TQ;
  TQ.SetMetricValue(6);

  Vert *pV0 = IQE.pVVert(0);
  Vert *pV1 = IQE.pVVert(1);
  double dMLength;
  int iNMetricValues = 6;
  double adMet[iNMetricValues];

  for (int i = 0; i < iNMetricValues; i++) {
    adMet[i] = pV0->getMetric(i) + pV1->getMetric(i);
    adMet[i] /= 2.0;
  }

  dMLength = TQ.MLength3D(pV0->getCoords(), pV1->getCoords(), adMet);
  dPri = dMLength - dMinMLength;

  if (dPri < 0) {
    std::set<Cell*> spCJunk0;
    std::set<Cell*> spCJunk1;
    std::set<Vert*> spV0;
    std::set<Vert*> spV1;

    std::set<Vert*>::iterator iterV, iterVEnd;

    findNeighborhoodInfo(pV0, spCJunk0, spV0);
    findNeighborhoodInfo(pV1, spCJunk1, spV1);

    int iSmallEdgeCount;

    iSmallEdgeCount = 0;

    iterV = spV0.begin();
    iterVEnd = spV0.end();

    for (; iterV != iterVEnd; ++iterV) {

      Vert *pVTest = (*iterV);
      double adMTest[iNMetricValues];

      for (int i = 0; i < iNMetricValues; i++)
	adMTest[i] = (pV0->getMetric(i) + pVTest->getMetric(i)) / 2.;

      if (TQ.MLength3D(pVTest->getCoords(), pV0->getCoords(), adMTest)
	  - dMinMLength < 0)
	iSmallEdgeCount++;
    }

    if (iSmallEdgeCount >= 2 && pV0->getVertType() != Vert::eBdryApex)
      pV0->markToDelete();

    iSmallEdgeCount = 0;

    iterV = spV1.begin();
    iterVEnd = spV1.end();

    for (; iterV != iterVEnd; ++iterV) {

      Vert *pVTest = (*iterV);
      double adMTest[iNMetricValues];

      for (int i = 0; i < iNMetricValues; i++)
	adMTest[i] = (pV1->getMetric(i) + pVTest->getMetric(i)) / 2.;

      if (TQ.MLength3D(pVTest->getCoords(), pV1->getCoords(), adMTest)
	  - dMinMLength < 0)
	iSmallEdgeCount++;
    }

    if (iSmallEdgeCount >= 2 && pV1->getVertType() != Vert::eBdryApex)
      pV1->markToDelete();

    if (!pV0->isDeletionRequested() && !pV1->isDeletionRequested())
      dPri = 1;
  }
  else
    dPri = 1;

  return dPri;
}

double
InsertionQueue::dEvaluateAnisoEdgePriority(const InsertionQueueEntry& IQE)
{

  InsertionQueueEntry::eEntryType eET = IQE.eType();
  GRUMMP::TMOPQual TQ;
  double dPri = 0;

  switch (eET)
    {
    case InsertionQueueEntry::e2DFace:

      switch (m_prio)
	{
	case eBestImpr:
	  dPri = d2DFaceBestImprove(IQE);
	  break;
	case eShortestEdge:
	  dPri = d2DFaceShortestEdge(IQE);
	  break;
	case eLongestEdge:
	  dPri = d2DFaceLongestEdge(IQE);
	  break;
	default:
	  assert(0);
	  break;
	}
      break;

    case InsertionQueueEntry::eBdryEdge:

      switch (m_prio)
	{
	case eBestImpr:
	  dPri = d2DBdryBestImprove(IQE);
	  break;
	case eShortestEdge:
	  dPri = d2DBdryShortestEdge(IQE);
	  break;
	case eLongestEdge:
	  dPri = d2DBdryLongestEdge(IQE);
	  break;
	default:
	  assert(0);
	  break;
	}
      break;

    case InsertionQueueEntry::eSubSeg:
//    	TQ.SetMetricValue(6);

      switch (m_prio)
	{
	case eBestImpr:
	  dPri = d3DSubSegBestImprove(IQE);
	  break;
	case eShortestEdge:
	  dPri = d3DSubSegShortestEdge(IQE);
	  break;
	case eLongestEdge:
	  dPri = d3DSubSegLongestEdge(IQE);
	  break;
	default:
	  assert(0);
	  break;
	}
      break;

    case InsertionQueueEntry::eInvalid:
      dPri = 1;
      break;

    default:
      break;
    }
  return (dPri);
}
#endif

////for frontal scheme
double
InsertionQueue::dPriorityFrontal(const InsertionQueueEntry IQE) const
{
  InsertionQueueEntry::eEntryType eET = IQE.eType();
  // Engwirda's Frontal Queue
  // Only triangle whose minimum edge length satisfys the lls is added to the queue
  // The queue priority is based on radius edge length.
  double dPri = 0;	//IQE.dPriority();
//  if (iFuzzyComp(dPri, IQE.dInvalidPriority) == 0) {
  if (1) {

    switch (eET)
      {

      case InsertionQueueEntry::eSubSeg:
	dPri = IQE.dSubSegPriority;
	break;

      case InsertionQueueEntry::eTriBFace:
      case InsertionQueueEntry::eBdryEdge:
	dPri = IQE.pBFBFace()->calcInsertionPriority();
	break;

      case InsertionQueueEntry::eTetCell:
      case InsertionQueueEntry::eTriCell:
	{
	  Cell * pC = IQE.pCCell();
//			int iNVerts = pC->getNumVerts();
//			int iDim = pC->getVert(0)->getSpaceDimen();

	  double dRadius = (dynamic_cast<SimplexCell*>(pC))->calcCircumradius();

	  double dShort = pC->calcShortestEdgeLength();
	  double dShortestToRadius;

	  // The following scale factor is used to normalize the
	  // shortest-edge-to-circumradius ratio to the range 0 to 1.
	  dShortestToRadius = dShort / (sqrt(3.) * dRadius);

	  /////
	  ///Find the shortest edge of the cell
	  Face * pFShortest = pC->getFace(0);
	  pFShortest =
	      pFShortest->calcSize() < pC->getFace(1)->calcSize() ?
		  pFShortest : pC->getFace(1);
	  pFShortest =
	      pFShortest->calcSize() < pC->getFace(2)->calcSize() ?
		  pFShortest : pC->getFace(2);
	  const Vert * pVA = pFShortest->getVert(0);
	  const Vert * pVB = pFShortest->getVert(1);
	  //////
	  //Only triangle whose minimum edge length satisfys the lls is added to the queue

	  if (iFuzzyComp(
	      pC->calcShortestEdgeLength(),
	      4. / 3. * 0.5 * (pVA->getLengthScale() + pVB->getLengthScale()))
	      < 0.) {
	    dPri = dShortestToRadius;
	  }
	  else {
	    dPri = IQE.dNeverQueuedPriority / 3.;
	  }

	  ///Get rid of good triangles
	  if (dynamic_cast<Mesh2D*>(pM)->DeterGoodTri(
	      dynamic_cast<TriCell*>(pC), pFShortest) == true) {
	    dPri = IQE.dNeverQueuedPriority;
	  }

	  if (dPri < IQE.dWorstCellPriority)
	    dPri = IQE.dWorstCellPriority;
	  logMessage(4, "Priority of a cell is %f, worst angle %f\n", dPri,
		     pC->calcMinSolid());

	}

	break;
      case InsertionQueueEntry::eFace:
      case InsertionQueueEntry::e2DFace:
      case InsertionQueueEntry::e3DFace:
      case InsertionQueueEntry::eInvalid:
	assert(0);
	break;
      }
  }

  IQE.setPriority(dPri);
  return (dPri);

}

double
InsertionQueue::dSizePriority(const InsertionQueueEntry IQE) const
{
  InsertionQueueEntry::eEntryType eET = IQE.eType();
  // Compute the priority the first time.  After that, return a stored value.
  double dPri = 0;    	//IQE.dPriority();
//  if (iFuzzyComp(dPri, IQE.dInvalidPriority) == 0) {
  if (1) {

    switch (eET)
      {

      case InsertionQueueEntry::eSubSeg:
	dPri = IQE.dSubSegPriority;
	break;

      case InsertionQueueEntry::eTriBFace:
      case InsertionQueueEntry::eBdryEdge:
	dPri = IQE.pBFBFace()->calcInsertionPriority();
	break;

      case InsertionQueueEntry::eTetCell:
      case InsertionQueueEntry::eTriCell:
	{
	  Cell * pC = IQE.pCCell();
	  double dSizeQual = 0, dShapeQual = 0;
	  int iNVerts = pC->getNumVerts();
	  int iDim = pC->getVert(0)->getSpaceDimen();

	  double dRadius = (dynamic_cast<SimplexCell*>(pC))->calcCircumradius();
	  double dSize = 2 / sqrt(double(iDim)) * dRadius;

	  // Find the correct size based on length scale.  Also, if the
	  // vertex is at a small angle in the input, make sure that you
	  // don't try to split this cell.

	  double dCorrectSize = 0;

	  for (int ii = 0; ii < iNVerts; ii++) {

	    Vert *pV = pC->getVert(ii);
	    double length_scale = pV->getLengthScale();
	    assert(iFuzzyComp(length_scale, 0.) == 1);

	    dCorrectSize += length_scale;
	    //double dThisSize = pV->dLS();
	    //dCorrectSize = MIN(dCorrectSize, dThisSize);

	    // Now, only want to do this with Tets.
	    //	  if (pV->isSmallAngleVert() && eET != InsertionQueueEntry::eTriCell)
	    //	    dShapeQual = IQE.dNeverQueuedPriority;

	  }

	  //Hack to enforce Miller et al. small angle
	  //condition on triangle meshes. SG Sept 13 2006.

	  //	if(eET == InsertionQueueEntry::eTriCell) {

	  TriCell* triangle = dynamic_cast<TriCell*>(pC);
	  assert(triangle);

	  if (triangle->isSmallAngleAffected()) {
	    dShapeQual = IQE.dNeverQueuedPriority;
	    dSizeQual = IQE.dNeverQueuedPriority;
	  }

	  //	}

	  dCorrectSize /= iNVerts;
	  logMessage(
	      4, "Correct size: %10.6f  Actual size: %10.6f  Ratio:  %10.7f\n",
	      dCorrectSize, dSize, dSize / dCorrectSize);
	  if (dSize > dCorrectSize && qSize) {
//    				dSizeQual = -dSize / dCorrectSize;
	    dSizeQual = -1. / pC->calcShortestEdgeLength();
	  }

	  else {

#ifndef EXPERIMENTAL
	    // This call works properly in both 2D and 3D.
	    double dShort = pC->calcShortestEdgeLength();
	    double dShortestToRadius;

	    // The following scale factor is used to normalize the
	    // shortest-edge-to-circumradius ratio to the range 0 to 1.

	    if (iDim == 2)
	      dShortestToRadius = dShort / (sqrt(3.) * dRadius);
	    else
	      dShortestToRadius = dShort * sqrt(.375) / dRadius;
	    if (dShortestToRadius < dWorstAllowable) {
	      //	  	if(m_prio == eWorstQuality)
	      //	  		dShapeQual += dShortestToRadius;
	      //		  dSizeQual = -1./pC->calcShortestEdgeLength();
	      //	  	else if(m_prio == eShortestEdge)
	      //	  		dShapeQual += dShort;
	      //	  	else if(m_prio == eLongestShortest)
	      //	  		dShapeQual += 1./dShort;

	    }
	    else {
	      dShapeQual += IQE.dNeverQueuedPriority;
	    }
#else
	    // 2D same as ever.
	    if (iDim == 2 && qShape) {
	      double dShort = pC->calcShortestEdgeLength();
	      dShapeQual += dShort / (sqrt(3.)*dRadius);
	    }

	    else {
	      double dSliverQual;
	      int ii;
	      vInToCircumSphere(pC, &ii, &dSliverQual);
	      logMessage(4, "This tet has in- to circum- ratio of %f\n",
		  dSliverQual);
	      dShapeQual += dSliverQual;
	    }
#endif
	  }

	  // The following code for including dihedral angle in choosing
	  // cells to split turns out not to be very effective in practice
	  // in eliminating bad tetrahedra.
	  //	if (dPri > 0.2 && eET == eTetCell) {
	  //	  // If the cell is a tet and the aspect ratio is reasonably
	  //	  // good, consider dihedral angle as well.
	  //	  int iN;
	  //	  double adRes[1];
	  //	  vMinDihed(pC, &iN, adRes);
	  //	  assert(iN == 1);
	  //	  // Convert from degrees back to radians.
	  //	  double dAnggleQual = adRes[0] * M_PI/180;
	  //	  dShapeQual = MIN(dShapeQual, dAngleQual*2);
	  //	}

	  //printf("size qual = %lf, shape qual = %lf \n", dSizeQual, dShapeQual);

	  if (qSize && qShape)
	    dPri = dSizeQual + dShapeQual;
	  else if (qSize)
	    dPri = dSizeQual;
	  else if (qShape)
	    dPri = dShapeQual;
//    		dSizeQual = dRadius/(pC->calcShortestEdgeLength());
//    		if(triangle->isSmallAngleAffected() or iFuzzyComp(triangle->calcMinSolid(), 59.) == 1){
//    			dPri = IQE.dNeverQueuedPriority;
//    		}
//    		else
//    		{
//    			dPri = -dSizeQual;
//    		}

	  if (dPri < IQE.dWorstCellPriority)
	    dPri = IQE.dWorstCellPriority;
	  logMessage(4, "Priority of a cell is %f (%f + %f), worst angle %f\n",
		     dPri, dSizeQual, dShapeQual, pC->calcMinSolid());

	  ///Find the shortest edge of the cell
	  const Face * pFShortest = pC->getFace(0);
	  pFShortest =
	      pFShortest->calcSize() < pC->getFace(1)->calcSize() ?
		  pFShortest : pC->getFace(1);
	  pFShortest =
	      pFShortest->calcSize() < pC->getFace(2)->calcSize() ?
		  pFShortest : pC->getFace(2);
//    		const Vert * pVA = pFShortest->getVert(0);
//    		const Vert * pVB = pFShortest->getVert(1);
//    		const Vert * pVC = pC->getOppositeVert(pFShortest);

//    		if(iFuzzyComp(pC->calcShortestEdgeLength(),0.5*(pVA->getLengthScale()+pVB->getLengthScale())) == 1){
//    			dPri = IQE.dNeverQueuedPriority/3.;
//    		}

//    		if (pC->getFace(0)->isBdryFace() or pC->getFace(1)->isBdryFace() or pC->getFace(2)->isBdryFace()){
//    			dPri += IQE.dWorstCellPriority/3.;
//    		}
//    		else if (iFuzzyComp(pC->calcShortestEdgeLength(), dCorrectSize) != 1
//    				and
//    				iFuzzyComp(dDist(pVA, pVC), dCorrectSize) != -1
//    				//iFuzzyComp(dDist(pVA, pVC), pVA->getLengthScale()/sqrt(2.)) != -1
//    				and
//    				iFuzzyComp(dDist(pVB, pVC), dCorrectSize) != -1)
//    				//iFuzzyComp(dDist(pVB, pVC), pVB->getLengthScale()/sqrt(2.)) != -1)
//    		{
//    			dPri += IQE.dWorstCellPriority/6.;
//    		}
//    		else if(iFuzzyComp(pC->calcShortestEdgeLength(), dCorrectSize) != 1
//    				and
//    				(iFuzzyComp(dDist(pVA, pVC), dCorrectSize) != -1
//    				//(iFuzzyComp(dDist(pVA, pVC), pVA->getLengthScale()/sqrt(2.)) != -1
//    				or
//    				iFuzzyComp(dDist(pVB, pVC), dCorrectSize) != -1))
//    				//iFuzzyComp(dDist(pVB, pVC), pVB->getLengthScale()/sqrt(2.)) != -1))
//    		{
//    			dPri += IQE.dWorstCellPriority/6.;
//    		}
//    		else
//    			dPri += IQE.dNeverQueuedPriority/3.;

	}

	break;
      default:
      case InsertionQueueEntry::eInvalid:
	assert(0);
	break;
      }
  }

  IQE.setPriority(dPri);
  return (dPri);

}

double
InsertionQueue::dFrontalModified(const InsertionQueueEntry IQE1) const
{
  InsertionQueueEntry::eEntryType eET = IQE1.eType();
  assert(eET == InsertionQueueEntry::e2DFace);
  // Compute the priority the first time.  After that, return a stored value.
  double dPri = 0;    	//IQE.dPriority();

  Face * IQE = IQE1.pFFace();
  Cell* Cell0 = IQE->getCell(0);
  Cell* Cell1 = IQE->getCell(1);
  TriCell* TriCell0 = NULL;
  TriCell* TriCell1 = NULL;
  Mesh2D* pM2D = dynamic_cast<Mesh2D*>(pM);
  double dShape = 2.;
//	int vertsindex = static_cast<int>( pM->getVertIndex(IQE->getVert(0))+ pM->getVertIndex(IQE->getVert(1)));
  if (IQE->isBdryFace()) {
    //Bdry Faces have higher priority.
    if (Cell0->isBdryCell()) {
      BFace* BF = (dynamic_cast<BFace*>(Cell0));
      if (BF->getBdryCondition() == this->keybdryconditionNum) {
	TriCell1 = dynamic_cast<TriCell*>(Cell1);
	if (TriCell1->isSmallAngleAffected()) {
	  dPri = IQE1.dNeverQueuedPriority;
	}
	else if (pM2D->DeterGoodTri(TriCell1, IQE) == false) {
	  dPri = -1e2;
	  dShape = TriCell1->calcShortestEdgeLength()
	      / (sqrt(3.) * TriCell1->calcCircumradius());
	}
	else {
	  dPri = IQE1.dNeverQueuedPriority;
	}
      }
      else {
	dPri = IQE1.dNeverQueuedPriority / 3.;
      }
    }
    else if (Cell1->isBdryCell()) {
      BFace* BF = (dynamic_cast<BFace*>(Cell1));
      if (BF->getBdryCondition() == this->keybdryconditionNum) {
	TriCell0 = dynamic_cast<TriCell*>(Cell0);
	if (TriCell0->isSmallAngleAffected()) {
	  dPri = IQE1.dNeverQueuedPriority;
	}
	if (pM2D->DeterGoodTri(TriCell0, IQE) == false) {
	  dPri = -1e2;
	  dShape = TriCell0->calcShortestEdgeLength()
	      / (sqrt(3.) * TriCell0->calcCircumradius());
	}
	else {
	  dPri = IQE1.dNeverQueuedPriority;
	}
      }
      else {
	dPri = IQE1.dNeverQueuedPriority / 3.;
      }
    }
    else {
      assert(false);
    }
  }
  //Interior
  else {
    TriCell1 = dynamic_cast<TriCell*>(Cell1);
    TriCell0 = dynamic_cast<TriCell*>(Cell0);
    bool bool0 = pM2D->DeterGoodTri(TriCell0, IQE)
	and (!TriCell0->isSmallAngleAffected());
    bool bool1 = pM2D->DeterGoodTri(TriCell1, IQE)
	and (!TriCell0->isSmallAngleAffected());

    if (bool0 xor bool1) {
      dPri = 0.;
      dShape = std::max(
	  TriCell1->calcShortestEdgeLength()
	      / (sqrt(3.) * TriCell1->calcCircumradius()),
	  TriCell0->calcShortestEdgeLength()
	      / (sqrt(3.) * TriCell0->calcCircumradius()));
    }
    else if (bool0 and bool1) {
      dPri = IQE1.dNeverQueuedPriority / 3.;
    }
    else {
      dPri = IQE1.dNeverQueuedPriority;
    }
  }

  dPri += IQE->calcSize() + IQE->layerindex + dShape;
//	dPri += IQE->calcSize()  + IQE->layerindex + pM->getFaceIndex(IQE)*1.e-10 ;

  IQE1.setPriority(dPri);
  return (dPri);

}

bool
InsertionQueueEntry::operator==(const InsertionQueueEntry& IQE) const
{
  if (eType() != IQE.eType())
    return false;
  switch (eType())
    {
    case InsertionQueueEntry::eTriCell:
    case InsertionQueueEntry::eTetCell:
      {
	if (pCCell() == IQE.pCCell())
	  return true;
	else
	  return false;
      }
    case InsertionQueueEntry::eTriBFace:
    case InsertionQueueEntry::eBdryEdge:
      {
	if (pBFBFace() == IQE.pBFBFace())
	  return true;
	else
	  return false;
      }
    case InsertionQueueEntry::eSubSeg:
      {
	if ((pVVert(0) == IQE.pVVert(0) && pVVert(1) == IQE.pVVert(1))
	    || (pVVert(0) == IQE.pVVert(1) && pVVert(1) == IQE.pVVert(0)))
	  return true;
	else
	  return false;
      }
    default:
      {
	// Should never get here.
	assert(0);
	break;
      }
    }
  return false;
}

bool
operator<(const InsertionQueueEntry& IQE1, const InsertionQueueEntry& IQE2)
{
  if (IQE1.eType() != IQE2.eType())
    return IQE1.eType() < IQE2.eType();

  // if different types sort arbitrarily
  //if same type, sort by pointer UNLESS they are equal in pointers
  switch (IQE1.eType())
    {
    case InsertionQueueEntry::eTriCell:
    case InsertionQueueEntry::eTetCell:
      {
	if (IQE1.pCCell() != IQE2.pCCell())
	  return IQE1.pCCell() < IQE2.pCCell();
	break;
      }
    case InsertionQueueEntry::eTriBFace:
    case InsertionQueueEntry::eBdryEdge:
      {
	if (IQE1.pBFBFace() != IQE2.pBFBFace())
	  return IQE1.pBFBFace() < IQE2.pBFBFace();
	break;
      }
    default:
      break;
    }
  // if one is in the mesh and the other isn't then return comparison there
//  if (IQE1.qStillInMesh() < IQE2.qStillInMesh())
  return IQE1.qStillInMesh() < IQE2.qStillInMesh();

//	int numVerts = 0;
//	// get num verts.. for now do it this way
//	switch (IQE1.eType()) {
//	// subsegs not used
//	case InsertionQueueEntry::eSubSeg:
//	{
//		numVerts = 2;
//		break;
//	}
//	case InsertionQueueEntry::eBdryEdge:
//	{
//		numVerts = 2;
//		break;
//	}
//	case InsertionQueueEntry::eTriBFace:
//	case InsertionQueueEntry::eTriCell:
//	{
//		numVerts = 3;
//		break;
//	}
//
//	case InsertionQueueEntry::eTetCell:
//	{
//		numVerts = 4;
//		break;
//	}
//	default:
//	{
//		// Should never get here.
//		assert(0);
//		break;
//	}
//	}

}

