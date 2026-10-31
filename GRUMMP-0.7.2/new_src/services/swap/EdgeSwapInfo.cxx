#include "GR_BFace.h"
#include "GR_EdgeSwapInfo.h"
#include "GR_Face.h"
#include "GR_FaceSwapInfo.h"
#include "GR_Util.h"
#include "GR_Vertex.h"

namespace GRUMMP
{
  EdgeSwapInfo::EdgeSwapInfo(const FaceSwapInfo3D& FC) :
      m_origBFace0(NULL), m_origBFace1(NULL), m_northVert(FC.getVertNorth()), m_southVert(
	  FC.getVertSouth()), m_numFaces(0), m_bestCanon(0), m_bestPerm(0), m_tooManyFaces(
	  false), m_isBdryEdge(false)
  {
    for (int i = 0; i < m_maxFaces + 1; i++) {
      m_outerVerts[i] = NULL;
      m_origCells[i] = NULL;
    }
    Face *face = FC.getFace();
    Vert *otherVert = FC.getVertOther();
    assert(face->isValid());
    assert(otherVert->isValid());
    gatherEdgeInfo(face, otherVert);

    for (GR_index_t i = 0; i < m_maxFaces + 1; i++) {
      if (m_outerVerts[i] != NULL)
	m_TotalVerts.insert(m_outerVerts[i]);
    }
    assert(m_northVert != NULL);
    assert(m_southVert != NULL);
    m_TotalVerts.insert(m_northVert);
    m_TotalVerts.insert(m_southVert);
  }

//@ Locally reconfigure tets incident on an edge to improve mesh quality
// The edge given by pVNorth, pVSouth is a candidate for removal.  If
// the edge doesn't lie on a boundary and a reconfiguration which
// improves mesh quality exists, then a swap is done, replacing 4 tets
// with 4, 5 with 6, 6 with 8, or 7 with 10.  Generally N tets are
// replaced by 2(N-2).  At this point, the largest transformation
// considered is 7 tets to 10.  The number of possible configurations
// for 8 to 12 swapping is somewhere around 150, so it's not
// considered, at least for now.
//
// After the edge swap is complete, all new faces are face swapped.
// The return value is the total number of face swaps done, with edge
// swaps counting as N-2 face swaps.

// FIX ME: Edge swapping could result in botching internal bdrys
// DAN: It does not, but don't swap Int bdry edges unless you have two ES objects
// and two corresponding objects.

/// Gather info about the neighborhood of the edge that is a candidate
/// for swapping.
  void
  EdgeSwapInfo::gatherEdgeInfo(Face* const face, Vert* const otherVert)
  {
    VALIDATE_INPUT(face->isValid());
    VALIDATE_INPUT(otherVert->isValid());
    VALIDATE_INPUT(face->hasVert(m_northVert));
    VALIDATE_INPUT(face->hasVert(m_southVert));
    VALIDATE_INPUT(face->hasVert(otherVert));

    for (unsigned i = 0; i < m_maxFaces; i++) {
      m_outerVerts[i] = pVInvalidVert;
      m_origCells[i] = pCInvalidCell;
    }

    // Determine the size of the orbit around the pivot edge and list
    // all the verts, faces, and cells involved.
    Face *oldFace = face;
    m_outerVerts[0] = otherVert;
    m_origCells[0] = face->getLeftCell();

    // Forget it if there's a non-tet involved.
    if (m_origCells[0]->getType() == Cell::eTriBFace) {
      gatherBdryEdgeInfo(face, otherVert);
      return;
    }
    if (m_origCells[0]->getType() == Cell::eIntTriBFace) {
      gatherIntBdryEdgeInfo(face, otherVert);
      return;
    }
    if (m_origCells[0]->getType() != Cell::eTet) {
      m_numFaces = 0;
      return;
    }

    int iReg = m_origCells[0]->getRegion();
    Face *newFace = m_origCells[0]->getOppositeFace(m_outerVerts[0]);
    int i = 1;
    bool bLoopAround = false;
    while (bLoopAround == false && i < 10) {
      // Abort if there is more than one region involved.
      //      if (newFace->getFaceLoc() == Face::eBdryTwoSide) {
      //	m_numFaces = 0;
      //	return;
      //      }
      if (newFace == face)
	bLoopAround = true;
      m_outerVerts[i] = m_origCells[i - 1]->getOppositeVert(oldFace);
      m_origCells[i] = newFace->getOppositeCell(m_origCells[i - 1]);
      // Hit a boundary.
      if (m_origCells[i]->getType() == Cell::eTriBFace) {
	gatherBdryEdgeInfo(newFace, m_outerVerts[i]);
	return;
      }
      if (m_origCells[i]->getType() == Cell::eIntTriBFace) {
	gatherIntBdryEdgeInfo(newFace, m_outerVerts[i]);
	return;
      }
      // Give up if more than one region is incident on the edge or if
      // we hit a non-tet.
      if ((m_origCells[i]->getType() != Cell::eTet)
	  || (iReg != m_origCells[i]->getRegion())) {
	m_numFaces = 0;
	return;
      }

      if (bLoopAround == false) {
	oldFace = newFace;
	newFace = m_origCells[i]->getOppositeFace(m_outerVerts[i]);
	i++;
      }

    }

    m_numFaces = i;

    if (newFace != face) {

      m_tooManyFaces = true;

      return; // Ran out of space
    }

    // TODO: Apparently edge swapping breaks something, either in
    // parallel or for aniso.  Whichever it is, that needs to be fixed.
//	////////////////////////////////
//	// erase this to edge swap
//	m_numFaces = 10;
//	////////////////////////////////

    if (m_numFaces > 7) {
      logMessage(4, "Too many tets around edge: %d\n", m_numFaces);
      m_tooManyFaces = true;
      return;
    }
  } // Done scanning around the edge

//@ Locally reconfigure tets incident on a bdry edge to improve mesh quality
// The bdry edge given by northVert, southVert is a candidate for removal.
// pF is required to be a boundary face, and otherVert is the third vertex
// on that face.
//
// If the bdry edge has two incident bdry faces that are coplanar and
// have the same boundary condition, and if a reconfiguration which
// improves mesh quality exists, then a swap is done, replacing 2 tets
// with 2, 3 with 4, ... or 6 with 10.  Generally N tets are replaced by
// 2(N-1).  At this point, the largest transformation considered is 6
// tets to 10.  The number of possible configurations for 7 to 12
// swapping is somewhere around 150, so it's not considered, at least
// for now.
//
// All of these reconfigurations are close cousins to interior edge
// swaps.  A 2-to-2 boundary swap is just a 3-to-2 swap with one tet
// flattened into oblivion.  This means that all of the canonical
// configuration info from interior can be re-used here.
//
// After the edge swap is complete, all new interior faces are face
// swapped. The return value is the total number of face swaps done,
// with boundary edge swaps counting as N-3 face swaps.
  void
  EdgeSwapInfo::gatherBdryEdgeInfo(Face* const pF, Vert* const otherVert)
  {
    // Note the optional args used for forced edge swapping during surface
    // recovery aren't used here because they aren't needed.
    VALIDATE_INPUT(pF->getType() == Face::eTriFace);
    VALIDATE_INPUT(pF->isValid());
    VALIDATE_INPUT(otherVert->isValid());
    VALIDATE_INPUT(pF->hasVert(m_northVert));
    VALIDATE_INPUT(pF->hasVert(m_southVert));
    VALIDATE_INPUT(pF->hasVert(otherVert));
    // The face must be on the boundary.
    VALIDATE_INPUT(
	pF->getLeftCell()->getType() == Cell::eTriBFace
	    || pF->getRightCell()->getType() == Cell::eTriBFace);

    m_isBdryEdge = true;
    for (unsigned i = 0; i < m_maxFaces + 1; i++) {
      m_outerVerts[i] = pVInvalidVert;
      m_origCells[i] = pCInvalidCell;
    }
    m_origBFace0 = pBFInvalidBFace;
    m_origBFace1 = pBFInvalidBFace;

    // Determine the size of the orbit around the pivot edge and list
    // all the verts, faces, and cells involved.
    Face *oldFace = pF;
    m_outerVerts[0] = otherVert;
    // ID the boundary associated with the original face *pF
    {
      Cell *pCTemp = pF->getLeftCell();
      if (pCTemp->getType() == Cell::eTriBFace) {
	m_origBFace0 = dynamic_cast<BFace*>(pCTemp);
	assert(m_origBFace0->isValid());
	m_origCells[0] = pF->getRightCell();
      }
      else {
	m_origBFace0 = dynamic_cast<BFace*>(pF->getRightCell());
	assert(m_origBFace0->isValid());
	m_origCells[0] = pCTemp;
      }
    }

    Face *newFace = m_origCells[0]->getOppositeFace(m_outerVerts[0]);
    int i = 1;
    int iReg = m_origCells[0]->getRegion();
    bool bLoopAround = false;
    while (bLoopAround == false && i && i <= 10) {
      // Abort if there is more than one region involved.

      if (newFace == pF)
	bLoopAround = true;

      if (newFace->getFaceLoc() == Face::eBdryTwoSide) {
	m_numFaces = 0;
	return;
      }
      m_outerVerts[i] = m_origCells[i - 1]->getOppositeVert(oldFace);
      m_origCells[i] = newFace->getOppositeCell(m_origCells[i - 1]);
      // Hit a boundary or non-tet cell.
      if (m_origCells[i]->getType() == Cell::eTriBFace)
	break;
      assert(m_origCells[i]->getType() == Cell::eTet);
      // Give up now if there is more than one region incident on the edge.
      if (m_origCells[i]->getRegion() != iReg) {
	m_numFaces = 0;
	return;
      }

      if (bLoopAround == false) {
	oldFace = newFace;
	newFace = m_origCells[i]->getOppositeFace(m_outerVerts[i]);
	i++;
      }

//		oldFace = newFace;
//		newFace = m_origCells[i]->getOppositeFace(m_outerVerts[i]);
//		i++;
    }
    int numOrigCells = i;
    m_numFaces = numOrigCells + 1;
    //    aiEdgeReq[numOrigCells]++;

    ////////////////////////////////
    // erase this to edge swap
    numOrigCells = 10;
    ////////////////////////////////

    // If we ran out of space while identifying the entities around the
    // edge we'd like to swap, the original count will be too high.
    if (numOrigCells > 6) {
      logMessage(4, "Too many tets around bdry edge: %d\n", numOrigCells);
      m_tooManyFaces = true;

      return;
    }

    // Assign data for the second boundary face
    m_origBFace1 = dynamic_cast<BFace*>(m_origCells[i]);
    assert(m_origBFace1->isValid());
  } // Done gathering info for bdry edge swaps.

// lets basically pretend we are doing a boundary edge swap, but do two of them..
// not sure if this function is really any different, may merge up.
  void
  EdgeSwapInfo::gatherIntBdryEdgeInfo(Face* const pF, Vert* const otherVert)
  {
    // Note the optional args used for forced edge swapping during surface
    // recovery aren't used here because they aren't needed.
    VALIDATE_INPUT(pF->getType() == Face::eTriFace);
    VALIDATE_INPUT(pF->isValid());
    VALIDATE_INPUT(otherVert->isValid());
    VALIDATE_INPUT(pF->hasVert(m_northVert));
    VALIDATE_INPUT(pF->hasVert(m_southVert));
    VALIDATE_INPUT(pF->hasVert(otherVert));
    // The face must be on the boundary.
    VALIDATE_INPUT(
	pF->getLeftCell()->getType() == Cell::eIntTriBFace
	    || pF->getRightCell()->getType() == Cell::eIntTriBFace);

    m_isBdryEdge = true;
    for (unsigned i = 0; i < m_maxFaces + 1; i++) {
      m_outerVerts[i] = pVInvalidVert;
      m_origCells[i] = pCInvalidCell;
    }
    m_origBFace0 = pBFInvalidBFace;
    m_origBFace1 = pBFInvalidBFace;

    // Determine the size of the orbit around the pivot edge and list
    // all the verts, faces, and cells involved.
    Face *oldFace = pF;
    m_outerVerts[0] = otherVert;
    // ID the boundary associated with the original face *pF
    {
      Cell *pCTemp = pF->getLeftCell();
      if (pCTemp->getType() == Cell::eIntTriBFace) {
	m_origBFace0 = dynamic_cast<BFace*>(pCTemp);
	assert(m_origBFace0->isValid());
	m_origCells[0] = pF->getRightCell();
      }
      else {
	m_origBFace0 = dynamic_cast<BFace*>(pF->getRightCell());
	assert(m_origBFace0->isValid());
	m_origCells[0] = pCTemp;
      }
    }

    Face *newFace = m_origCells[0]->getOppositeFace(m_outerVerts[0]);
    int i = 1;
    int iReg = m_origCells[0]->getRegion();
    bool bLoopAround = false;
    while (bLoopAround == false && i <= 10) {

      if (newFace == pF)
	bLoopAround = true;

      m_outerVerts[i] = m_origCells[i - 1]->getOppositeVert(oldFace);
      m_origCells[i] = newFace->getOppositeCell(m_origCells[i - 1]);
      // Hit a boundary or non-tet cell.
      if (m_origCells[i]->getType() == Cell::eIntTriBFace)
	break;
      // Give up now if there is more than one region incident on the edge.
      if (m_origCells[i]->getRegion() != iReg
	  || m_origCells[i]->getType() != Cell::eTet) {
	m_numFaces = 0;
	return;
      }
      assert(m_origCells[i]->getType() == Cell::eTet);

      if (bLoopAround == false) {
	oldFace = newFace;
	newFace = m_origCells[i]->getOppositeFace(m_outerVerts[i]);
	i++;
      }

//		oldFace = newFace;
//		newFace = m_origCells[i]->getOppositeFace(m_outerVerts[i]);
//		i++;
    }
    int numOrigCells = i;
    m_numFaces = numOrigCells + 1;

    ////////////////////////////////
    // erase this to edge swap
    numOrigCells = 10;
    ////////////////////////////////

    //    aiEdgeReq[numOrigCells]++;
    // If we ran out of space while identifying the entities around the
    // edge we'd like to swap, the original count will be too high.
    if (numOrigCells > 6) {

      logMessage(4, "Too many tets around bdry edge: %d\n", numOrigCells);
      m_tooManyFaces = true;
      return;
    }

    // Assign data for the second boundary face
    m_origBFace1 = dynamic_cast<BFace*>(m_origCells[i]);
    assert(m_origBFace1->isValid());
  } // Done gathering info for int bdry edge swaps.

}
