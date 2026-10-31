#ifndef GR_InsertionQueue
#define GR_InsertionQueue 1

#include <iostream>
#include <set>
#include <map>
#include <limits.h>
#include "GR_config.h"
#include "GR_misc.h"
#include "GR_Cell.h"
#include "GR_InsertionQueueEntry.h"
#include "GR_Mesh.h"
#include "GR_TMOPQual.h"

using namespace std;

struct InsertionCriteria {
  bool
  operator()(const Face* pF1, const Face* pF2) const
  {
    if (pF1 == pF2)
      return false;

    double dPF1, dPF2;
    dPF1 = (pF1->getdQualInc());
    dPF2 = (pF2->getdQualInc());

    return (dPF1 < dPF2);
  }
};

class InsertionQueue {
public:

  enum eQueuePriority {
    eWorstQuality,
    eShortestEdge,
    eLongestShortest,
    eBestImpr,
    eLongestEdge,
    eInvalid
  };
private:

  std::multimap<double, InsertionQueueEntry> MMap;
  double dWorstAllowable;
  GR_index_t iMaxCellCount;
  // If we've cleaned up initial encroachment in this queue, then
  // there's no longer any need to check whether new bdry entities are
  // encroached by points already in the mesh --- they can't be, because
  // no points are ever inserted where they will encroach anything.
  bool qIsInitialEncroachmentFixed;
  // refine for one or the other, or both
  bool qShape;
  bool qSize;
  eQueuePriority m_prio;
  Mesh *pM;
  int keybdryconditionNum = 1;

#ifdef HAVE_MESQUITE
  GRUMMP::TMOPQual::eQMetric m_eIsoQM;
  GRUMMP::TMOPQual::eQMetric m_eAnisoQM;
  double dMinMLength;
  double dMaxMLength;
  double dMTol;
#endif

  // Don't see any good reason to allow default or copy construction, or
  // operator=.
//	InsertionQueue(); // Never defined // Reza commented this
  InsertionQueue(const InsertionQueue&); // Never defined
//	InsertionQueue operator=(const InsertionQueue&); // Never defined // Reza commented this
public:
  InsertionQueue() :
      MMap(), dWorstAllowable(0), iMaxCellCount(INT_MAX), qIsInitialEncroachmentFixed(
	  true), qShape(true), qSize(true), m_prio(), pM(NULL)
#ifdef HAVE_MESQUITE
							 , m_eIsoQM(), m_eAnisoQM(), dMinMLength(
	  0), dMaxMLength(LARGE_DBL), dMTol(0)
#endif
  {
  } // Reza added this for parallel refinement sake
    //******************************************************
    // This constructor has been added to be called by parallel InsertionQueue.
    // it creats an empty container on contrary to other constructors.
    // it also does not call "queueEncroachedBdryEntities" function which
    // frankly I don't know if it causes any problem or not.
  InsertionQueue(Mesh * const pMIn, const bool JustToBeDifferent,
		 const double dQualTarget, const bool bdry_unencroached) :
      MMap(), dWorstAllowable(0), iMaxCellCount(INT_MAX), qIsInitialEncroachmentFixed(
	  bdry_unencroached), qShape(true), qSize(true), m_prio(), pM(NULL)
#ifdef HAVE_MESQUITE
								      , m_eIsoQM(), m_eAnisoQM(), dMinMLength(
	  0), dMaxMLength(LARGE_DBL), dMTol(0)
#endif
  {
    logMessage(1, "Initializing insertion priority queue...\n");
    logMessage(4, "Running Parallel Insertion Queue: %d\n",
	       int(JustToBeDifferent));
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
    logMessage(2, "Worst allowed cell shape measure is %f.\n", dWorstAllowable);
    MMap.clear();
  }

  //**********************************************
  InsertionQueue(Mesh * const pMIn, const double dTarget =
		     InsertionQueueEntry::dInvalidPriority,
		 const bool bdry_unencroached = false);
//	InsertionQueue(Mesh * const pMIn,
//			const double dQualTarget,
//			const bool bdry_unencroached, std::set<Cell*> &cellQueue);

  InsertionQueueEntry
  IQETopValidEntry()
  {
    // Grab the first entry in the multimap, and return the data.
    // MMap.begin() returns an iterator which (through overloading of
    // operator ->) can access the contents of the pair<double,
    // InsertionQueueEntry> that constitutes the first mmap entry.
    assert(!qIsQueueEmpty());
    InsertionQueueEntry IQETop = MMap.begin()->second;
    logMessage(3, "Grabbed queue entry with quality measure %f\n",
	       MMap.begin()->first);
    if (IQETop.eType() == InsertionQueueEntry::eTetCell
	|| IQETop.eType() == InsertionQueueEntry::eTriCell)
      qIsInitialEncroachmentFixed = true;
    return IQETop;
  }

//	bool SearchInQueue(Face *pFaceToSearch);

public:
  void
  clearQueue()
  {
    MMap.clear();
  }
  void
  setKeyBdryCondition(int KBC)
  {
    keybdryconditionNum = KBC;
  }
  void
  vCleanQueue()
  // Remove all null entries from the top of the queue.
  {
    if (qIsQueueEmpty())
      return;
    logMessage(4, "Cleaning queue of length %u...", iQueueLength());
    InsertionQueueEntry IQETop = MMap.begin()->second;
    while (!qIsQueueEmpty() && !IQETop.qStillInMesh()) {
      logMessage(5, "Popping entry for entity eliminated from mesh.\n");
      MMap.erase(MMap.begin());
      if (!MMap.empty())
	IQETop = MMap.begin()->second;
    }
    logMessage(4, "now of length %u.\n", iQueueLength());
  }

  double
  dTopPriority()
  {
    assert(!qIsQueueEmpty());
    // Grab the first entry in the multimap, and return its priority.
    // MMap.begin() returns an iterator which (through overloading of
    // operator ->) can access the contents of the pair<double,

    // InsertionQueueEntry> that constitutes the first mmap entry.
    return MMap.begin()->first;
  }
//	void print_queue() {
//		std::multimap<double, InsertionQueueEntry>::iterator it;
//		for(it = MMap.begin(); it != MMap.end(); it++) {
//			logMessage(1,"%p: Still in mesh = %d, deleted = %d, priority = %lf\n",
//					it->second.pCCell(), it->second.qStillInMesh(),
//					it->second.pCCell()->isDeleted(), it->first);
//		}
//	}

  void
  print_queue_to_Screen()
  {
    std::multimap<double, InsertionQueueEntry>::iterator it;
    for (it = MMap.begin(); it != MMap.end(); it++) {
      printf("%p: Still in mesh = %d, deleted = %d, priority = %lf\n",
	     it->second.pCCell(), it->second.qStillInMesh(),
	     it->second.pCCell()->isDeleted(), it->first);
    }
  }
  void
  vPopTopEntry()
  {
    assert(!qIsQueueEmpty());
    MMap.erase(MMap.begin());
//		vCleanQueue();
  }

  GR_index_t
  iQueueLength() const
  {
    return (MMap.size());
  }
  double
  dEvaluatePriority(const InsertionQueueEntry IQE) const;
  double
  dPriorityFrontal(const InsertionQueueEntry IQE) const;
  double
  dSizePriority(const InsertionQueueEntry IQE) const;
  double
  dFrontalModified(const InsertionQueueEntry IQE) const;
  int iDelaunayorEngwirda = 0;		//Queue: Delaunay: 0 Engwirda: 1

#ifdef HAVE_MESQUITE
////	Fucntions for Aniso Edge Adaptation queue
  double
  dEvaluateAnisoEdgePriority(const InsertionQueueEntry &IQE);
  void
  vSetMinMLength(double dML)
  {
    dMinMLength = dML;
  }
  void
  vSetMaxMLength(double dML)
  {
    dMaxMLength = dML;
  }
  void
  vSetMetricTol(double dMT)
  {
    dMTol = dMT;
  }
//  double
//  getMinMLength()
//  {
//    return dMinMLength;
//  }
//  double
//  getMaxMLength()
//  {
//    return dMaxMLength;
//  }
//  double
//  getMetricTol()
//  {
//    return dMTol;
//  }
  bool
  qAddEntryAnisoEdge(InsertionQueueEntry IQE);
  void
  setIsoQMetric(GRUMMP::TMOPQual::eQMetric eQM)
  {
    m_eIsoQM = eQM;
  }
  ;
  void
  setAnisoQMetric(GRUMMP::TMOPQual::eQMetric eQM)
  {
    m_eAnisoQM = eQM;
  }
  ;
  GRUMMP::TMOPQual::eQMetric
  SetQMetricType(int iNMetricValues, double dTolerance, std::set<Vert*> setpV);
  GRUMMP::TMOPQual::eQMetric
  getIsoQMetric()
  {
    return m_eIsoQM;
  }
  ;
  GRUMMP::TMOPQual::eQMetric
  getAnisoQMetric()
  {
    return m_eAnisoQM;
  }
  ;

  double
  d2DFaceBestImprove(const InsertionQueueEntry& IQE);
  double
  d2DFaceLongestEdge(const InsertionQueueEntry& IQE);
  double
  d2DFaceShortestEdge(const InsertionQueueEntry& IQE);

  double
  d2DBdryBestImprove(const InsertionQueueEntry& IQE);
  double
  d2DBdryLongestEdge(const InsertionQueueEntry& IQE);
  double
  d2DBdryShortestEdge(const InsertionQueueEntry& IQE);

  double
  d3DSubSegBestImprove(const InsertionQueueEntry& IQE);
  double
  d3DSubSegLongestEdge(const InsertionQueueEntry& IQE);
  double
  d3DSubSegShortestEdge(const InsertionQueueEntry& IQE);

#endif

public:
  bool
  qIsQueueEmpty()
  {
    return MMap.empty();

  }

  bool
  qAddFrontal(InsertionQueueEntry IQE1, const double dOffset = 0, int filter =
		  -1)
  {
    // Return value of insert is an iterator; not currently checked.

//		double dKey = dFrontal(IQE1) + dOffset;
    double dKey = dFrontalModified(IQE1) + dOffset;
//	  double dKey = dNewFrontal(IQE1) + dOffset;
    IQE1.updateFilter(filter);

    if (dKey < IQE1.dNeverQueuedPriority / 2) {

      logMessage(4, "Worst allowed %f., Priority %f\n", dWorstAllowable, dKey);

      MMap.insert(std::pair<const double, InsertionQueueEntry>(dKey, IQE1));
      return true;
    }
    else {
      return false;
    }

  }
  bool
  qAddEntrywithPriority(InsertionQueueEntry IQE, double priorityvalue)
  {
    if (priorityvalue == 0)
      return false;
    else {
      MMap.insert(
	  std::pair<const double, InsertionQueueEntry>(priorityvalue, IQE));
      return true;
    }
  }

  // Reza dropped the "const" before "InsertionQueueEntry" so that he can modify the IQE
  bool
  qAddEntry(InsertionQueueEntry IQE, const double dOffset = 0, int filter = -1)
  {
    // Return value of insert is an iterator; not currently checked.
    double dKey = dEvaluatePriority(IQE) + dOffset;
    if (iDelaunayorEngwirda == 1) {
      dKey = dPriorityFrontal(IQE) + dOffset;
    }
    else {
      dKey = dEvaluatePriority(IQE) + dOffset;
    }
    //double dKey = dSizePriority(IQE) + dOffset;
    //double dKey = dAdjacentPriority(IQE) + dOffset;
//		double dKey = dFrontal(IQE) + dOffset;
    IQE.updateFilter(filter);
    //      if (IQE.eType() == InsertionQueueEntry::eTriCell) {
    //	TriCell *pTC = (dynamic_cast<TriCell*>(IQE.pCCell()));
    //	assert(pTC);
    //	double adCircCent[2];
    //	pTC->calcCircumcenter(adCircCent);
    //      }

    if (dKey < IQE.dNeverQueuedPriority / 2) {
      /* This hack is designed to fix interior tets.  It works.  But
       most of the bad ones tend to be at the boundaries, so it
       doesn't have much net effect. */
      /* || */
      /* 	  (IQE.eType() == InsertionQueueEntry::eTetCell && */
      /* 	   (!IQE.pCCell()->pVVert(0)->qIsBdryVert() && */
      /* 	    !IQE.pCCell()->pVVert(1)->qIsBdryVert() && */
      /* 	    !IQE.pCCell()->pVVert(2)->qIsBdryVert() && */
      /* 	    !IQE.pCCell()->pVVert(3)->qIsBdryVert()) && */
      /* 	   dKey < dWorstAllowable*2) */
      /* 	if (IQE.eType() == InsertionQueueEntry::eTetCell || */
      /* 	    IQE.eType() == InsertionQueueEntry::eTriCell) { */
      /* 	  double adCent[3]; adCent[2] = 0; */
      /* 	  IQE.pCCell()->vCentroid(adCent); */
      /* 	  vMessage(3, "Key: %10.6f.  Centroid: (%10.6f, %10.6f, %10.6f)\n", */
      /* 		   dKey, adCent[0], adCent[1], adCent[2]); */
      /* 	} */
      logMessage(4, "Worst allowed %f., Priority %f\n", dWorstAllowable, dKey);
      MMap.insert(std::pair<const double, InsertionQueueEntry>(dKey, IQE));
      return true;
    }
    else {
      return false;
    }
  }

  // Can be used both to add all initially encroached BFaces and
  // subsegments to the insertion queue and to subsequently add
  // newly-created entities that are still encroached.
  void
  vSetQualityTarget(const double dNewTarget)
  {
    dWorstAllowable = dNewTarget;
  }
  double
  dQualityTarget() const
  {
    return dWorstAllowable;
  }
  void
  vSetMaxCells(const int iMax)
  {
    iMaxCellCount = iMax;
  }
  GR_index_t
  iMaxCellsAllowed() const
  {
    return iMaxCellCount;
  }
  /// Resets the Queue, used when observing for a long time.
  void
  emptyQueue()
  {
    while (!qIsQueueEmpty())
      MMap.erase(MMap.begin());
  }
  void
  qSetShape(bool shape)
  {
    qShape = shape;
  }
  void
  qSetSize(bool size)
  {
    qSize = size;
  }
  void
  setPriority(eQueuePriority prio)
  {
    m_prio = prio;
  }
  eQueuePriority
  getPriority()
  {
    return m_prio;
  }
  int
  getSize()
  {
    return MMap.size();
  }

};

#endif

