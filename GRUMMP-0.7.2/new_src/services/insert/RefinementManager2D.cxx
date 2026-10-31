/*
 * TriRefinementManager.cxx
 *
 *  Created on: 2013-05-15
 *      Author: cfog
 */

#include "GR_Geometry.h"
#include "GR_Aniso.h"
#include "GR_InsertionManager.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_Mesh2D.h"
#include "GR_Face.h"
#include <utility>

namespace GRUMMP {
void RefinementManager2D::receiveDeletedCells(
		std::vector<Cell*>& deletedCells) {
	if (!IsInParallel) {
		logMessage(4,
				"RefinementManager2D received word about %zu deleted cells.\n",
				deletedCells.size());
		// Can't actually recompute the priority to find this thing in the queue,
		// and at this point there's no other way to find it, so do nothing.
	} else {
		logMessage(4,
				"RefinementManager2D received word about %zu deleted cells in Parallel.\n",
				deletedCells.size());
		// Here I'll strip all the cells that have been deleted from the TrisinQuad vectors in buffer leaves.
		// However, I'm not sure whether I have to remove it from the one in the main leaf if it happened to fall there.
	}
}

void RefinementManager2D::receiveCreatedCells(
		std::vector<Cell*>& createdCells) {
	if (!IsInParallel) {
		logMessage(3,
				"RefinementManager2D received word about %zu newly created cells.\n",
				createdCells.size());
		std::vector<Cell*>::iterator iter = createdCells.begin(), iterEnd =
				createdCells.end();
		for (; iter != iterEnd; ++iter) {
			if (m_mesh->isSimplicial()
					|| (*iter)->getType() == Cell::eTriCell) {
				TriCell *tri = dynamic_cast<TriCell*>(*iter);
				// Some get added, some are too good to bother with.
				if (Queueonface == 0) {
					addToPriorityQueue(tri);
				} else {
					if (m_sizingField) {
						if (tri->getVert(0)->getVertType() != Vert::eBdryApex)
							m_sizingField->setLengthScale(tri->getVert(0));
						if (tri->getVert(1)->getVertType() != Vert::eBdryApex)
							m_sizingField->setLengthScale(tri->getVert(1));
						if (tri->getVert(2)->getVertType() != Vert::eBdryApex)
							m_sizingField->setLengthScale(tri->getVert(2));
					}
					addToFrontal(tri->getFace(0));
					addToFrontal(tri->getFace(1));
					addToFrontal(tri->getFace(2));

				}

//
	}
      }
    }
    else {
      double RCC = omp_get_wtime();
      logMessage(
	  3,
	  "RefinementManager2D received word about %zu newly created cells in Parallel.\n",
	  createdCells.size());

#pragma omp parallel shared(createdCells)
		{
		  std::vector<Cell*>::iterator iter;
#pragma omp for schedule(guided)
	for (iter = createdCells.begin(); iter < createdCells.end(); ++iter) {
	  TriCell *tri = dynamic_cast<TriCell*>(*iter);
	  // Some get added, some are too good to bother with.
	  if (!tri->isValid() || tri->isDeleted())
	    continue;
	  Rect* LeafB = QuadT->WhichRect(tri);
	  LeafB->Push_Back_Cell(tri);
//				if (!OkToGo)
//				{
//					int filter = QuadT->FilterTriangle(tri,LeafB);
//					if (filter <= Filter)
//						OkToGo = true;
//				}
//				LeafB->Push_Back_Cell(tri);
			}
		}		// End of parallel
//		Remaining_Cells.clear();
		RCC = omp_get_wtime() - RCC;
#pragma omp atomic
		SendEventsCCTime += RCC;
	}
}

void RefinementManager2D::receiveCreatedBFaces(
		std::vector<BFace*>& createdBFaces) {
	if (!IsInParallel) {
		logMessage(3,
				"RefinementManager2D received word about %zu newly created bfaces.\n",
				createdBFaces.size());
		std::vector<BFace*>::iterator iter = createdBFaces.begin(), iterEnd =
				createdBFaces.end();
		for (; iter != iterEnd; ++iter) {
			BdryEdgeBase* bdryEdge = dynamic_cast<BdryEdgeBase*>(*iter);
			if (bdryEdge->isDeleted())
				continue;
			if (m_refineBdryForSize
					|| bdryEdge->isEncroached(m_mesh->getEncroachmentType())
							!= eClean)
				m_bdryFacesToSplit.push_back(bdryEdge);
		}
	} else {
		double RCB = omp_get_wtime();
		std::deque<BdryEdgeBase*> bufdeq;
		logMessage(3,
				"RefinementManager2D received word about %zu newly created bfaces in Parallel.\n",
				createdBFaces.size());
		std::vector<BFace*>::iterator iter = createdBFaces.begin(), iterEnd =
				createdBFaces.end();
		for (; iter != iterEnd; ++iter) {
			BFace * BF = *iter;
			if (BF->isValid() && !BF->isDeleted()) {
				Face * pF = BF->getFace();
				TriCell* pTC = dynamic_cast<TriCell*>(pF->getOppositeCell(BF));
				if (pTC->isValid() && !pTC->isDeleted()) {
					Rect* LeafB = QuadT->WhichRect(pTC);
					if (BF->isEncroached(m_mesh->getEncroachmentType())
							!= eClean)
						LeafB->m_newBdryFaces.push_back(BF);
				}
			}
		}
		RCB = omp_get_wtime() - RCB;
#pragma omp atomic
		SendEventsCBTime += RCB;
	}
}

void RefinementManager2D::fixBoundaryEncroachment() {
	logMessage(1, "Fixing boundary encroachment\n");
	queueEncroachedBdryEdges();

	while (!m_bdryFacesToSplit.empty()) {

		InsertionQueueEntry IQE = m_bdryFacesToSplit.front();
		BdryEdgeBase* bdryEdge = dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace());
		if (!IQE.qStillInMesh()) {
			m_bdryFacesToSplit.pop_front();
		} else if (splitBdryEdge(bdryEdge)) {
			// If this succeeded, that means that we split the edge.  If it failed,
			// then we tried to insert on the wrong side of some existing edge,
			// which needs to be split first.  In the latter case, we *don't*
			// pop the front entry, because it's now a different bdry edge.
			m_bdryFacesToSplit.pop_front();
			m_mesh->sendEvents();
#ifndef NDEBUG
			m_mesh->writeTempMesh();
#endif
		}
	}

}

void RefinementManager2D::updateEncroachedBFaceQueue() {
	//			if(pointType != InsertionPointCalculator::CIRCUMCENTER){
	//				logMessage(1,"Rejecting point of type %d at %f %f due to encroachment\n",pointType,insertionPoint[0],insertionPoint[1]);
	//			}
	std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
			m_WI->getEncroachedBFaces();
	for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
			++iter) {
		BdryEdgeBase* pBEB = dynamic_cast<BdryEdgeBase*>(*iter);
		//			m_bdryFacesToSplit.push_back(*range.first);
		m_bdryFacesToSplit.push_back(pBEB);
	}
}

bool RefinementManager2D::splitTriangle(TriCell* const cellToSplit,
		bool Para_Flag, int ID) {
	if (!Para_Flag) {
		assert(cellToSplit);
		assert(m_bdryFacesToSplit.empty());

		double insertionPoint[3];
		//	writeVTK_SingleCell(cellToSplit,"tempcell.vtk");
		//	writeVTKLegacyWithoutPurge(*m_mesh,"tempmesh");
		//  Calculate the Insertion Point
		InsertionPointCalculator::ePointType pointType =
				m_IPC->calcInsertionPoint(cellToSplit, insertionPoint);
		double circCent[2];
		cellToSplit->calcCircumcenter(circCent);
		//			double delta2 = dDIST2D(insertionPoint,circCent)/cellToSplit->calcCircumradius();
		//			if(delta2 >0.000001)
		//			printf("delta2 %.16f\n",delta2);
		//Compute watson hull.
		m_WI->computeHull(insertionPoint, cellToSplit);
		//	logMessage(1,"inserting point of type %d into triangle at %f %f\n",pointType,insertionPoint[0],insertionPoint[1]);
		//	logMessage(1,"( %.2f, %.2f) ( %.2f, %.2f) ( %.2f, %.2f)\n",cellToSplit->getVert(0)->x(),cellToSplit->getVert(0)->y(),
		//			cellToSplit->getVert(1)->x(),cellToSplit->getVert(1)->y(),
		//			cellToSplit->getVert(2)->x(),cellToSplit->getVert(2)->y());
		// if a non-circumcenter encroaches, then recalculate
		//	if(pointType != InsertionPointCalculator::CIRCUMCENTER && m_WI->areBFacesEncroached()){
		//		cellToSplit->calcCircumcenter(insertionPoint);
		//		m_WI->computeHull(insertionPoint, cellToSplit);
		//		pointType = InsertionPointCalculator::CIRCUMCENTER;
		//	}
		//If point is encroaching, queue all encroached edges and return.
		if (m_WI->areBFacesEncroached()) {
			//			if(pointType != InsertionPointCalculator::CIRCUMCENTER){
			//				logMessage(1,"Rejecting point of type %d at %f %f due to encroachment\n",pointType,insertionPoint[0],insertionPoint[1]);
			//			}
			updateEncroachedBFaceQueue();
			return false;
		}

		if (!m_WI->isHullValid()) {
			cellToSplit->calcCentroid(insertionPoint);
			m_WI->computeHull(insertionPoint, cellToSplit);
			if (m_WI->areBFacesEncroached()) {
				updateEncroachedBFaceQueue();
				return false;
			}
		}

		m_IPC->updateStats(pointType);

		//Else, insert the point in m_mesh.
		Vert* newVert = m_WI->insertPointInHull();

		newVert->setType(Vert::eInterior);
		newVert->setParentEntity(static_cast<GeometryEntity*>(NULL));
		// TODO Once a sizing field or equivalent has been added, it'll need to be
		// updated here.
		if (m_sizingField)
			m_sizingField->setLengthScale(newVert);

		// The priority queue is updated elsewhere.
		return true;
	} else {
		assert(cellToSplit);
		assert(!cellToSplit->isDeleted());
		assert(Para_m_bdryFacesToSplit[ID].empty());
		// TODO Generalize this to call a function that computes the split point.
		//Find circumcenter.
		double circCenter[3];
		cellToSplit->calcCircumcenter(circCenter);

		// TODO Generalize this to allow an arbitrary inserter.
		//Compute watson hull.
//		double TimeFIfST = omp_get_wtime();
		Para_m_WI[ID]->computeHull(circCenter, cellToSplit);
		//If point is encroaching, queue all encroached edges and return.
		if (Para_m_WI[ID]->areBFacesEncroached()) {
			std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
					Para_m_WI[ID]->getEncroachedBFaces();
			for (std::set<BFace*>::iterator iter = range.first;
					iter != range.second; ++iter) {
				BdryEdgeBase *pBEB = dynamic_cast<BdryEdgeBase*>(*iter);
				Para_m_bdryFacesToSplit[ID].push_back(pBEB);
			}
//			TimeFIfST = omp_get_wtime() - TimeFIfST;
//			ThreadTime[ID][50] += TimeFIfST;
			return false;
		}

//		double TimeSIfST = omp_get_wtime();
		if (!Para_m_WI[ID]->isHullValid()) {
			cellToSplit->calcCentroid(circCenter);
			Para_m_WI[ID]->computeHull(circCenter, cellToSplit);
			if (Para_m_WI[ID]->areBFacesEncroached()) {
				std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
						Para_m_WI[ID]->getEncroachedBFaces();
				for (std::set<BFace*>::iterator iter = range.first;
						iter != range.second; ++iter) {
					BdryEdgeBase *pBEB = dynamic_cast<BdryEdgeBase*>(*iter);
					Para_m_bdryFacesToSplit[ID].push_back(pBEB);
				}
//				TimeSIfST = omp_get_wtime() - TimeSIfST;
//				ThreadTime[ID][51] += TimeSIfST;
				return false;
			}
		}
		logMessage(4, "inserting point into triangle at %f %f\n", circCenter[0],
				circCenter[1]);

		//Else, insert the point in m_mesh.
		Vert* newVert = Para_m_WI[ID]->insertPointInHull();

		newVert->setType(Vert::eInterior);
		newVert->setParentEntity(static_cast<GeometryEntity*>(NULL));
		// TODO Once a sizing field or equivalent has been added, it'll need to be
		// updated here.
#pragma omp critical(LengthScale)
		{
			if (m_sizingField)
				m_sizingField->setLengthScale(newVert);
		}
//		VertsToSetLengthScale[ID].push_back(newVert);
		//TODO: There's no guaranty that threads will not go further
		// to collide with each other in calling "setLengthScale" simultaneously.
//		if (Para_m_sizing_field[ID]) Para_m_sizing_field[ID]->setLengthScale(newVert);

		// The priority queue is updated elsewhere.
		return true;
	}

}

bool RefinementManager2D::splitBdryEdge(BdryEdgeBase* bdryEdgeToSplit,
		bool Para_Flag, int ID) {
	if (!Para_Flag) {
		assert(bdryEdgeToSplit);
		assert(!bdryEdgeToSplit->isDeleted());

		// TODO This next clause assumes that we're doing a Ruppert/Shewchuk style
		// insertion, which won't always be true.

		{
			std::vector<Cell*> new_cells;
			deleteVertsInEdgeBall(bdryEdgeToSplit);
		}

		//Find the split location and the corresponding parameter.
		bool concentric_shell_split;
		double newLoc[3];
		double splitParam;
		{
			CubitVector split_location;
			splitParam = bdryEdgeToSplit->getSplitData(
					//BdryEdgeBase::MID_TVT, split_location, concentric_shell_split);
					BdryEdgeBase::MID_PARAM, split_location,
					concentric_shell_split);
			newLoc[0] = split_location.x();
			newLoc[1] = split_location.y();
			newLoc[2] = split_location.z();

		}

		Cell* seedCell = bdryEdgeToSplit->getFace(0)->getOppositeCell(
				bdryEdgeToSplit);
		logMessage(3, "Splitting bdry edge %p (adj cell %p).\n",
				bdryEdgeToSplit, seedCell);
		logMessage(3, "Bdry edge coords: ");
		logMessage(3, " (%9f, %9f), (%9f, %9f)\n",
				bdryEdgeToSplit->getVert(0)->x(),
				bdryEdgeToSplit->getVert(0)->y(),
				bdryEdgeToSplit->getVert(1)->x(),
				bdryEdgeToSplit->getVert(1)->y());
		logMessage(3, "Seed cell coords: \n");
		logMessage(3, " (%9f, %9f), (%9f, %9f), (%9f, %9f)\n",
				seedCell->getVert(0)->x(), seedCell->getVert(0)->y(),
				seedCell->getVert(1)->x(), seedCell->getVert(1)->y(),
				seedCell->getVert(2)->x(), seedCell->getVert(2)->y());

		m_WI->computeHull(newLoc, seedCell, bdryEdgeToSplit, splitParam);
		// Before we go any further, let's confirm that we're not inserting on the
		// wrong side of some other encroached edge associated with curved geometry.

		double midCoords[] = { 0.5
				* (bdryEdgeToSplit->getVert(0)->x()
						+ bdryEdgeToSplit->getVert(1)->x()), 0.5
				* (bdryEdgeToSplit->getVert(0)->y()
						+ bdryEdgeToSplit->getVert(1)->y()) };

		std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
				m_WI->getEncroachedBFaces();
		for (std::set<BFace*>::iterator iter = range.first;
				iter != range.second; ++iter) {
			BdryEdgeBase *pBEB = dynamic_cast<BdryEdgeBase*>(*iter);
			if (pBEB == bdryEdgeToSplit)
				continue;
			int orientMid = checkOrient2D(pBEB->getVert(0)->getCoords(),
					pBEB->getVert(1)->getCoords(), midCoords);
			assert(orientMid != 0);
			int orientNew = checkOrient2D(pBEB->getVert(0)->getCoords(),
					pBEB->getVert(1)->getCoords(), newLoc);
			if (orientMid != orientNew) {
				m_bdryFacesToSplit.push_front(pBEB);
				return false;
			}
			m_bdryFacesToSplit.push_back(pBEB);
		}
		//	logMessage(1,"Inserting %f %f\n",newLoc[0],newLoc[1]);

		Vert* newVert = m_WI->insertPointInHull();

		if (bdryEdgeToSplit->getType() == Cell::eIntBdryEdge)
			newVert->setType(Vert::eBdryTwoSide);
		else
			newVert->setType(Vert::eBdry);
		newVert->setParentEntity(bdryEdgeToSplit->getCurve());
		if (concentric_shell_split)
			newVert->markAsShellVert(true);

		if (m_sizingField)
			m_sizingField->setLengthScale(newVert);

		return true;
	}
	//**********************************************************************
	//************************* Parallel section ***************************
	//**********************************************************************
	else {
		assert(bdryEdgeToSplit);
		assert(!bdryEdgeToSplit->isDeleted());

		// TODO This next clause assumes that we're doing a Ruppert/Shewchuk style
		// insertion, which won't always be true.
		deleteVertsInEdgeBall(bdryEdgeToSplit);

		//Find the split location and the corresponding parameter.
		bool concentric_shell_split;
		double newLoc[3];
		double splitParam;
		{
			CubitVector split_location;
			splitParam = bdryEdgeToSplit->getSplitData(BdryEdgeBase::MID_TVT,
					split_location, concentric_shell_split);
			newLoc[0] = split_location.x();
			newLoc[1] = split_location.y();
			newLoc[2] = split_location.z();

		}
		//		printf("point %f %f %d\n",newLoc[0],newLoc[1],concentric_shell_split);
		//		bdryEdgeToSplit->getFace(0)->printFaceInfo();
		Cell* seedCell = bdryEdgeToSplit->getFace(0)->getOppositeCell(
				bdryEdgeToSplit);
		logMessage(3, "Splitting bdry edge %p (adj cell %p).\n",
				bdryEdgeToSplit, seedCell);
		logMessage(3, "Bdry edge coords: \n");
		logMessage(3, " (%9f, %9f), (%9f, %9f)\n",
				bdryEdgeToSplit->getVert(0)->x(),
				bdryEdgeToSplit->getVert(0)->y(),
				bdryEdgeToSplit->getVert(1)->x(),
				bdryEdgeToSplit->getVert(1)->y());
		logMessage(3, "Seed cell coords: \n");
		logMessage(3, " (%9f, %9f), (%9f, %9f), (%9f, %9f)\n",
				seedCell->getVert(0)->x(), seedCell->getVert(0)->y(),
				seedCell->getVert(1)->x(), seedCell->getVert(1)->y(),
				seedCell->getVert(2)->x(), seedCell->getVert(2)->y());

		Para_m_WI[ID]->computeHull(newLoc, seedCell, bdryEdgeToSplit,
				splitParam);
		// Before we go any further, let's confirm that we're not inserting on the
		// wrong side of some other encroached edge associated with curved geometry.

		double midCoords[] = { 0.5
				* (bdryEdgeToSplit->getVert(0)->x()
						+ bdryEdgeToSplit->getVert(1)->x()), 0.5
				* (bdryEdgeToSplit->getVert(0)->y()
						+ bdryEdgeToSplit->getVert(1)->y()) };

		std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
				Para_m_WI[ID]->getEncroachedBFaces();
		for (std::set<BFace*>::iterator iter = range.first;
				iter != range.second; ++iter) {
			BdryEdgeBase *pBEB = dynamic_cast<BdryEdgeBase*>(*iter);
			if (pBEB == bdryEdgeToSplit)
				continue;
			int orientMid = checkOrient2D(pBEB->getVert(0)->getCoords(),
					pBEB->getVert(1)->getCoords(), midCoords);
			assert(orientMid != 0);
			int orientNew = checkOrient2D(pBEB->getVert(0)->getCoords(),
					pBEB->getVert(1)->getCoords(), newLoc);
			if (orientMid != orientNew) {
				Para_m_bdryFacesToSplit[ID].push_front(pBEB);
				return false;
			}
			//		if (bdryEdgeToSplit->getType() == Cell::eIntBdryEdge) {
			//			if(pBEB->getType() == Cell::eIntBdryEdge ) {
			Para_m_bdryFacesToSplit[ID].push_back(pBEB);
			//			} else {
			//				m_bdryFacesToSplit.push_front(pBEB);
			//				return false;
			//			}
			//		}

		}
		//	m_WI->writeCavityVTK("cavity.vtk");
		//	m_WI->writeHullVTK("hull.vtk");
		//
		//	writeVTKLegacyWithoutPurge(*m_mesh,"tempmesh");
		//

		Vert* newVert = Para_m_WI[ID]->insertPointInHull();

		if (bdryEdgeToSplit->getType() == Cell::eIntBdryEdge)
			newVert->setType(Vert::eBdryTwoSide);
		else
			newVert->setType(Vert::eBdry);
		newVert->setParentEntity(bdryEdgeToSplit->getCurve());
		if (concentric_shell_split)
			newVert->markAsShellVert(true);

#pragma omp critical(LengthScale)
		{
			if (m_sizingField)
				m_sizingField->setLengthScale(newVert);
		}
//		VertsToSetLengthScale[ID].push_back(newVert);
//		if (Para_m_sizing_field[ID]) Para_m_sizing_field[ID]->setLengthScale(newVert);

		return true;
	}
}

bool RefinementManager2D::addToFrontal(Face * const face) {
	assert(!IsInParallel);

	if (face->isDeleted())
		return false;

	if (face->isBdryFace()) {
		if (BdryConditionFromFace(face) == this->keybdryconditionNum) {
//				if ((face->getVert(0)->x()*face->getVert(0)->x() + face->getVert(0)->y()*face->getVert(0)->y()) < (1.5*1.5))
	face->layerindex = 0;
      }
    }
    ///test if face is already in queue.
//		if (m_IQ.SearchInQueue(face)){
//			return false;
//		}
	return m_IQ.qAddFrontal(face);

}

bool RefinementManager2D::addToPriorityQueue(Cell * const cell) {
	if (!IsInParallel) {
		if (cell->isDeleted())
			return false;
		// really hacky...
		//	double quality = 	cell->calcShortestEdgeLength()/(sqrt(3.)*cell->calcCircumradius());
		//	if( quality < m_IQ.dQualityTarget()){
		return m_IQ.qAddEntry(cell);
		//	}
		//	return false;
	} else {
		if (cell->isDeleted())
			return false;
		TriCell * pTC = dynamic_cast<TriCell*>(cell);
		Rect* LeafB = QuadT->WhichRect(pTC);
		if (!OkToGo) {
			int filter = QuadT->FilterTriangle(pTC, LeafB);
			if (filter <= Filter)
				OkToGo = true;
		}
/// The following commands are for the case when you want the filtering in effect.
//				if (filter <= Filter)
//				{
//					assert(QuadT->QContainsCell(pTC,ThreadLeaves[ID]));
//					assert(QuadT->QContainsCell(pTC,LeafB));
//					ThreadLeaves[ID]->TrisInQuad.push_back(pTC);
//					return Para_m_IQ[ID].qAddEntry(cell);
//					return true;
//				}
//				else
//				{
//					assert(QuadT->QContainsCell(pTC,ThreadLeaves[ID]));
//					assert(QuadT->QContainsCell(pTC,LeafB));
//					ThreadLeaves[ID]->TrisInQuad.push_back(pTC);
//					return true;
//				}
//			}
#pragma omp critical(receiveCreatedCells)
		{
			LeafB->TrisInQuad.push_back(pTC);
		}
		return true;
	}
	return false;
}

void RefinementManager2D::queueEncroachedBdryEdges() {

	// This loop handles both regular and internal bdry faces
	for (GR_index_t i = 0; i < m_mesh->getNumTotalBdryFaces(); i++) {
		BdryEdgeBase* bdryEdge =
				dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		assert(bdryEdge);
		if (bdryEdge->isDeleted())
			continue;
		if (bdryEdge->isEncroached(m_mesh->getEncroachmentType()) != eClean)
			m_bdryFacesToSplit.push_back(bdryEdge);
	}
}

GR_index_t RefinementManager2D::queueAllCells(bool Para_Flag, int NQLayers) {
	if (!Para_Flag) {
		GR_index_t cellsAdded = 0;
		for (GR_index_t i = 0; i < m_mesh->getNumTriCells(); i++) {
			TriCell* tri = dynamic_cast<TriCell*>(m_mesh->getCell(i));
			assert(tri);
			if (tri->isDeleted())
				continue;
			if (addToPriorityQueue(tri)) {
				cellsAdded++;
			}
		}

		assert(cellsAdded <= this->m_IQ.iQueueLength());
		return cellsAdded;
	} else {
		double MaxCircR = 0;
		double XCirB[2] = { 0 };
		double YCirB[2] = { 0 };
		GR_index_t cellsAdded = 0;
		GR_index_t MaxCells = m_mesh->getNumTriCells();
		omp_set_dynamic(0);
		omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(MaxCells, NQLayers)
		{
			double MaxCircRT = 0;
			double XCirBT[2] = { 0 };
			double YCirBT[2] = { 0 };
			double CirPos[2] = { 0 };
#pragma omp for nowait schedule(static)
			for (GR_index_t i = 0; i < MaxCells; i++) {
				TriCell* tri = dynamic_cast<TriCell*>(m_mesh->getCell(i));
				assert(tri);
				if (tri->isDeleted())
					continue;
				if (tri->calcCircumradius() > MaxCircRT)
					MaxCircRT = tri->calcCircumradius();
				tri->calcCircumcenter(CirPos);
				if (CirPos[0] >= XCirBT[1])
					XCirBT[1] = CirPos[0];
				if (CirPos[0] <= XCirBT[0])
					XCirBT[0] = CirPos[0];
				if (CirPos[1] >= YCirBT[1])
					YCirBT[1] = CirPos[1];
				if (CirPos[1] <= YCirBT[0])
					YCirBT[0] = CirPos[1];
			}
#pragma omp critical(queueAllCells)
			{
				if (XCirBT[1] >= XCirB[1])
					XCirB[1] = XCirBT[1];
				if (XCirBT[0] <= XCirB[0])
					XCirB[0] = XCirBT[0];
				if (YCirBT[1] >= YCirB[1])
					YCirB[1] = YCirBT[1];
				if (YCirBT[0] <= YCirB[0])
					YCirB[0] = YCirBT[0];
				if (MaxCircRT >= MaxCircR)
					MaxCircR = MaxCircRT;
			}
#pragma omp barrier
#pragma omp master
			{
				QuadT->UpdateMaxCirR(MaxCircR);
				QuadT->CalcBoundingBox();
				QuadT->UpdateBBound(XCirB, YCirB);
				QuadT->Initial_Quadtree_Creator(true, NQLayers);
				for (GR_index_t i = 0; i < MaxCells; i++) {
					TriCell* tri = dynamic_cast<TriCell*>(m_mesh->getCell(i));
					assert(tri);
					if (tri->isDeleted())
						continue;
					QuadT->AddTriToRect(tri);
					cellsAdded++;
				}
			}
		}
		QuadT->FillWBZ(true);
		QuadT->FindIndSets(IndependentSet, true);
		QuadT->writeVTKquadtree();
		assert(cellsAdded <= this->m_IQ.iQueueLength());
		return cellsAdded;
	}
}

int RefinementManager2D::BdryConditionFromFace(Face* face) {
	///Codes here determine if the boundary condition of the bdry face is 1
	///Starts***
	assert(face->isBdryFace());
	BFace* bdryF = NULL;
	Cell* bdcell = NULL;
	if (face->getLeftCell()->isBdryCell()) {
		bdcell = face->getLeftCell();
	} else {
		bdcell = face->getRightCell();
	}
	bdryF = (dynamic_cast<BFace*>(bdcell));
	assert(bdryF != NULL);
	int bdryConditionNumber = bdryF->getBdryCondition();
	return bdryConditionNumber;
}

/// This function queue all the faces ******Queue all boundary faces for now*******
GR_index_t RefinementManager2D::queueAllFaces(bool Para_Flag, int NQLayers) {

	(void) Para_Flag;
	(void) NQLayers;	//They are not used. Silence a warning
	GR_index_t facesAdded = 0;
	for (GR_index_t i = 0; i < m_mesh->getNumFaces(); i++) {
		Face* face = (dynamic_cast<Face*>(m_mesh->getFace(i)));
		assert(face);

		if (face->isDeleted())
			continue;
		if (face->isBdryFace()) {

			///Only Queue Faces with Bdry Condition 1
//				if( BdryConditionFromFace(face) != 1){
//					continue;
//				}

			if (addToFrontal(face)) {
				facesAdded++;
			}
		}
	}

	return facesAdded;

}

/// This function updates the insertionQueue for each thread meeting the filter criteria.
GR_index_t RefinementManager2D::queueAllCellsThreadWise(int ID, int FLevel) {
	assert(Para_m_IQ[ID].qIsQueueEmpty());
	assert(Para_m_bdryFacesToSplit[ID].empty());
	std::list<TriCell*> BL;
	GR_index_t NumAdded = 0;
	if (ThreadLeaves[ID]->TrisInQuad.empty())
		return 0;
	std::list<TriCell*>::iterator itS = ThreadLeaves[ID]->TrisInQuad.begin();
	for (; itS != ThreadLeaves[ID]->TrisInQuad.end(); itS++) {
		TriCell* pTC = *itS;
		if (!pTC->isValid() || pTC->isDeleted()) {
			continue;
		}
		bool checkQ = QuadT->QContainsCell(pTC, ThreadLeaves[ID]);
		if (checkQ) {
//			int ff = QuadT->FilterTriangle(pTC,ThreadLeaves[ID]);
//			if (ff <= FLevel)
//			{
//				NumAdded++;
			Para_m_IQ[ID].qAddEntry(pTC);
//			}
//			else
//				BL.push_back(pTC);
		} else
			continue;
	}
	logMessage(4, "Level of refinement is: %d\n", FLevel);
	ThreadLeaves[ID]->TrisInQuad.clear();
	ThreadLeaves[ID]->TrisInQuad.swap(BL);

	std::vector<BFace*>::iterator ItBS =
			ThreadLeaves[ID]->m_newBdryFaces.begin(), ItBE =
			ThreadLeaves[ID]->m_newBdryFaces.end();
	for (; ItBS != ItBE; ItBS++) {
		BFace* BF = *ItBS;
		if (!BF->isValid() || BF->isDeleted())
			continue;
		Face * pF = BF->getFace();
		TriCell* pTC = dynamic_cast<TriCell*>(pF->getOppositeCell(BF));
		if (!pTC->isValid() || pTC->isDeleted())
			continue;
		if (BF->isEncroached(m_mesh->getEncroachmentType()) != eClean)
			Para_m_bdryFacesToSplit[ID].push_back(
					dynamic_cast<BdryEdgeBase*>(BF));
	}
	ThreadLeaves[ID]->m_newBdryFaces.clear();

	return NumAdded;
}

GR_index_t RefinementManager2D::refineMeshonfaces(const GR_index_t maxMeshVerts,
		bool ParaFlag, int NQLayers) {
	logMessage(1, "In refineMesh; parallel flag is %s\n",
			ParaFlag ? "true" : "false");
	IsInParallel = ParaFlag;
	if (!IsInParallel) {
		m_bdryFacesToSplit.clear();
		m_IQ.clearQueue();
		queueAllFaces();
		printf("the length of queue: %d\n", m_IQ.iQueueLength());
//		queueAllCells();
		return refineForQualityonfaces(maxMeshVerts);
	} else
		return refineForQuality(maxMeshVerts, ParaFlag, NQLayers);
}

GR_index_t RefinementManager2D::refineMesh(const GR_index_t maxMeshVerts,
		bool ParaFlag, int NQLayers) {
	logMessage(1, "In refineMesh; parallel flag is %s\n",
			ParaFlag ? "true" : "false");
	IsInParallel = ParaFlag;
	if (!IsInParallel) {
		m_bdryFacesToSplit.clear();
		m_IQ.clearQueue();
//		queueAllFaces();
		queueAllCells();
		return refineForQuality(maxMeshVerts);
	} else
		return refineForQuality(maxMeshVerts, ParaFlag, NQLayers);
}

GR_index_t RefinementManager2D::refineAllQueuedCells(
		const GR_index_t maxMeshVerts) {

	logMessage(1, "%d cells queued for refinement.\n", m_IQ.iQueueLength());
	return refineForQuality(maxMeshVerts);
}

///This subroutine has the same functionality as "refineForQuality" only within the given thread.
GR_index_t RefinementManager2D::refinForQualityFilteredbyQuadtree(int ID,
		int Fil, GR_index_t &RVerts, const GR_index_t maxMeshVerts) {
	double TTimeST = 0, TTimeSBE = 0;
	queueAllCellsThreadWise(ID, Fil);

	RVerts = Para_m_IQ[ID].iQueueLength();
	logMessage(3, "Refining for quality in parallel for leaf %d\n",
			ThreadLeaves[ID]->Identifier);
	assert(m_refineBdryForSize == false);
	GR_index_t vertsAdded = 0;
	GR_index_t maxVertsAdded = maxMeshVerts - m_mesh->getNumVerts();
	// This next line handles the case of unsigned integer wraparound from
	// subtracting to get a negative.  Using GR_sindex_t might overflow instead
	// of underflowing.
	if (maxMeshVerts < m_mesh->getNumVerts())
		maxVertsAdded = 0;

	Para_m_IQ[ID].vCleanQueue();

	while (!Para_m_IQ[ID].qIsQueueEmpty() && vertsAdded < maxVertsAdded) {
		while (!Para_m_bdryFacesToSplit[ID].empty()) {
			InsertionQueueEntry IQE = Para_m_bdryFacesToSplit[ID].front();

			BdryEdgeBase* bdryEdge = dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace());
			if (!IQE.qStillInMesh()) {
				Para_m_bdryFacesToSplit[ID].pop_front();
				continue;
			}

			// This fails only if there's another edge we should split first.
			double TimeSBE = omp_get_wtime();
			bool success = splitBdryEdge(bdryEdge, true, ID);
			TimeSBE = omp_get_wtime() - TimeSBE;
			TTimeSBE += TimeSBE;
			if (success) {
				++vertsAdded;
				Para_m_bdryFacesToSplit[ID].pop_front();
#ifndef NDEBUG
				m_mesh->writeTempMesh();
#endif
				if (vertsAdded >= maxVertsAdded)
					return vertsAdded;
			}

		}

		Para_m_bdryFacesToSplit[ID].clear();
		Para_m_IQ[ID].vCleanQueue();

		if (Para_m_IQ[ID].qIsQueueEmpty())
			break;
		InsertionQueueEntry IQE = Para_m_IQ[ID].IQETopValidEntry();
		assert(IQE.eType() == InsertionQueueEntry::eTriCell);
		bool success;
		double TimeST = omp_get_wtime();
		success = splitTriangle(dynamic_cast<TriCell*>(IQE.pCCell()), true, ID);
		TimeST = omp_get_wtime() - TimeST;
		TTimeST += TimeST;
		// We might not have split yet in favor of splitting encroached
		// bdry faces.

		if (success) {
			++vertsAdded;
			Para_m_IQ[ID].vPopTopEntry();
#ifndef NDEBUG
			m_mesh->writeTempMesh();
#endif

		}

		// The following two lines are useful if the refinement is crashing,
		// because it gives you a sense of when to stop it to get a near-crash
		// picture of what's going on.
//		      if (vertsAdded % 5 == 0)
//		        logMessage(1, "%d  ", m_mesh->getNumVerts());
	}
	ThreadTime[ID][25] += TTimeSBE;
	ThreadTime[ID][26] += TTimeST;

	assert(
			(Para_m_IQ[ID].qIsQueueEmpty()
					&& Para_m_bdryFacesToSplit[ID].empty())
					|| vertsAdded >= maxVertsAdded);
//	Para_m_IQ[ID].clearQueue();
//	Para_m_bdryFacesToSplit[ID].clear();
	return vertsAdded;
}

Vert*
RefinementManager2D::InsertEntity(
		InsertionPointCalculator::ePointType pointType,
		GR_index_t& vertsAdded) {
	m_IPC->updateStats(pointType);
	//Else, insert the point in m_mesh.
	Vert* newVert = m_WI->insertPointInHull();
	newVert->setType(Vert::eInterior);
	newVert->setParentEntity(static_cast<GeometryEntity*>(NULL));
	// TODO Once a sizing field or equivalent has been added, it'll need to be
	// updated here.
	if (m_sizingField)
		m_sizingField->setLengthScale(newVert);

	// The priority queue is updated elsewhere.
	// We might not have split yet in favor of splitting encroached
	// bdry faces.
	++vertsAdded;
//
	m_mesh->sendEvents();
#ifndef NDEBUG
	m_mesh->writeTempMesh();
#endif
	return newVert;
}

/// each pass of parallel refinement consists of iterating amongst all ...
/// leaves, updating the insertionQueue from their TrisInQuad container meeting ...
/// the filter criteria, and then moving to the next pass with a finer filter limit.
GR_index_t RefinementManager2D::refineForQuality(const GR_index_t maxMeshVerts,
		bool Para_Flag, int NQLayer) {
	double ToTime = omp_get_wtime();
	if (Para_Flag) {
		//Mesh2D * M2D = dynamic_cast<Mesh2D*>(m_mesh);
		//writeVTKLegacy(*M2D, "outputBeforInsertion");
		// 1. A set of independent threads would be determined.
		// 2. threads will be allocated to the top leaves on the set.
		// 3. queueAllCells would update the insertionQueue and TrisInQuad of ...
		//    ThreadLeaves and other leaves respectively.
		// 4.
		IsInParallel = Para_Flag;
		logMessage(1, "In refinement; parallel flag is %s\n",
				IsInParallel ? "true" : "false");
		double QTime = omp_get_wtime();
		queueAllCells(Para_Flag, NQLayer);
		QTime = omp_get_wtime() - QTime;
		QueueingTime += QTime;
		FilterLevel = 3;
		Granularity = 1;
		Filter = 0;
		GR_index_t TotalvertsAdded = 0;
		GR_index_t RemainingVerts = 0;
		GR_index_t vertexAddedSoFar = 0;
		GR_index_t maxVertsAdded = maxMeshVerts - m_mesh->getNumVerts();
		if (maxMeshVerts < m_mesh->getNumVerts())
			maxVertsAdded = 0;
		Filter += Granularity;
		double ParallelTime = 0;
		SyncEntryTime += m_mesh->Mesh2D_sync_entry(ParallelTime, false);
		m_mesh->IsInParallel(true);
		double TimeThisPass = 0;
		omp_set_dynamic(0);
		omp_set_num_threads(NUM_PROCS);
#pragma omp parallel private(TimeThisPass) shared(RemainingVerts, TotalvertsAdded)
		{
#pragma omp single
			{
				int TotalvertsAddedOld = 0;
				int passes = 1;
				bool FlagNewLayer = true;
				do {
					TimeThisPass = omp_get_wtime();
//				StayInLoop = false;

					// The following loop will be parallelised using tasks.
					// Each thread grab a leaf from the back of IndependetSet, Then ...
					// it initialises a task to call "refinForQualityFilteredbyQuadtree"
					RemainingVerts = 0;
					if (!OkToGo)
						logMessage(2,
								"pass #%d: Refining %lu leaves with filter level of %d ... ",
								passes, QuadT->initial_quads.size(), Filter);
					while (!QuadT->IsWBLEmpty() || !IndependentSet.empty()) {
						double PaTime = omp_get_wtime();
						while (!IndependentSet.empty()) {
							std::set<Rect*>::reverse_iterator ritI =
									IndependentSet.rbegin();
							Rect * LeafI = *ritI;
							IndependentSet.erase(LeafI);
#pragma omp task firstprivate(LeafI) shared(RemainingVerts, TotalvertsAdded)
							{
								double thread_time = omp_get_wtime();
								int VertsAddedThisPass = 0;
								int ID = omp_get_thread_num();
								ThreadLeaves[ID] = LeafI;
								GR_index_t RVerts = 0;
								VertsAddedThisPass =
										refinForQualityFilteredbyQuadtree(ID,
												Filter, RVerts, maxMeshVerts);
#pragma omp atomic
								TotalvertsAdded += VertsAddedThisPass;

								thread_time = omp_get_wtime() - thread_time;
								ThreadTime[ID][100] += thread_time;
								ThreadInNum[ID][100] += VertsAddedThisPass;
							}
						}
#pragma omp taskwait
						PaTime = omp_get_wtime() - PaTime;
						JustParallelTime += (PaTime);

						double FISTime = omp_get_wtime();
						QuadT->FindIndSets(IndependentSet, FlagNewLayer);
						FISTime = omp_get_wtime() - FISTime;
						FIndSetsTime += FISTime;
					}

					OkToGo = false;
					double SendingTime = omp_get_wtime();
					m_mesh->sendEvents();
					SendingTime = omp_get_wtime() - SendingTime;
					SendEventsTime += SendingTime;
					if (OkToGo && RemainingVerts > 0) {
						IndependentSet.clear();
						QuadT->FillWBZ(false);
						continue;
					}
					if (Filter < FilterLevel) {
						Filter += Granularity;
						Granularity += 0;
					} else if (Filter <= 2 * FilterLevel) {
						Granularity += 1;
						Filter += Granularity;
					} else {
						Granularity += 0;
						Filter += Granularity;
					}

					TimeThisPass = omp_get_wtime() - TimeThisPass;
					vertexAddedSoFar = TotalvertsAdded - TotalvertsAddedOld;
					logMessage(2, "%d verticies inserted in %4.10f seconds\n",
							vertexAddedSoFar, TimeThisPass);
					TotalvertsAddedOld = TotalvertsAdded;

					double QCTime = omp_get_wtime();
					int OldQTSize = QuadT->initial_quads.size();
					passes++;
					QuadT->QuadTreeCreatorSecondPasses(passes);
					int NewQTSize = QuadT->initial_quads.size();
					QCTime = omp_get_wtime() - QCTime;
					SecondPassesTime += QCTime;
					if (NewQTSize == OldQTSize) {
						FlagNewLayer = false;
					} else {
						FlagNewLayer = true;
						QuadT->clearKTInLeaves();
					}
					IndependentSet.clear();
					QuadT->FillWBZ(FlagNewLayer);

				} while (TotalvertsAdded < maxVertsAdded
						&& vertexAddedSoFar > 0.005 * TotalvertsAdded);
			}				//********* End of Single ******
		} //**** End of Parallel ************
		m_mesh->IsInParallel(false);
		double SETime = omp_get_wtime();
		m_mesh->Mesh2D_sync_exit(false);
		SETime = omp_get_wtime() - SETime;
		SyncExitTime += SETime;
		ToTime = omp_get_wtime() - ToTime;
//		writeVTKLegacy(*M2D, "outputAfterInsertion");
		TotalRefineTime += ToTime;
		m_mesh->purgeAllEntities();
		m_mesh->IsInParallel(false);
		return TotalvertsAdded;
	} else // Reza: This else is simply to keep the serial version in place
	{
		ToTime = omp_get_wtime();
		assert(m_refineBdryForSize == false);
		logMessage(1, "Refining for quality\n");
		logMessage(1, "Apparently refining in serial\n");

		GR_index_t vertsAdded = 0;
		GR_index_t maxVertsAdded = maxMeshVerts - m_mesh->getNumVerts();
		// This next line handles the case of unsigned integer wraparound from
		// subtracting to get a negative.  Using GR_sindex_t might overflow instead
		// of underflowing.
		if (maxMeshVerts < m_mesh->getNumVerts())
			maxVertsAdded = 0;

		m_IQ.vCleanQueue();

		while (!m_IQ.qIsQueueEmpty() && vertsAdded < maxVertsAdded) {
			while (!m_bdryFacesToSplit.empty()) {
				InsertionQueueEntry IQE = m_bdryFacesToSplit.front();

				BdryEdgeBase* bdryEdge =
						dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace());
				if (!IQE.qStillInMesh()) {
					m_bdryFacesToSplit.pop_front();
					continue;
				}

				// This fails only if there's another edge we should split first.
				bool success = splitBdryEdge(bdryEdge);
				if (success) {
					++vertsAdded;
					m_bdryFacesToSplit.pop_front();
#ifndef NDEBUG
					m_mesh->writeTempMesh();
#endif
					if (vertsAdded >= maxVertsAdded)
						return vertsAdded;
				}

			}

			m_mesh->sendEvents();
			m_IQ.vCleanQueue();

			if (m_IQ.qIsQueueEmpty())
				break;
			// Get the first entity in the queue
			InsertionQueueEntry IQE = m_IQ.IQETopValidEntry();
			assert(IQE.eType() == InsertionQueueEntry::eTriCell);

			bool success = true;
			TriCell* cellToSplit = (dynamic_cast<TriCell*>(IQE.pCCell()));
			assert(cellToSplit);
			assert(m_bdryFacesToSplit.empty());
			if (cellToSplit->isSmallAngleAffected()) {
				m_IQ.vPopTopEntry();
				continue;
			}

			double insertionPoint[3];
			//  Calculate the Insertion Point
			InsertionPointCalculator::ePointType pointType =
					m_IPC->calcInsertionPoint(cellToSplit, insertionPoint);
			//Compute watson hull.
			m_WI->computeHull(insertionPoint, cellToSplit);
			//If point is encroaching, queue all encroached edges and return.
			if (m_WI->areBFacesEncroached()) {
				updateEncroachedBFaceQueue();
				success = false;
			}

			if (!m_WI->isHullValid()) {
				cellToSplit->calcCentroid(insertionPoint);
				m_WI->computeHull(insertionPoint, cellToSplit);
				if (m_WI->areBFacesEncroached()) {
					updateEncroachedBFaceQueue();
					success = false;
				}
			}
			// InsertEntity
			if (success) {
				Vert* newVert = InsertEntity(pointType, vertsAdded);
				assert(newVert->isValid());
			}
			// The following two lines are useful if the refinement is crashing,
			// because it gives you a sense of when to stop it to get a near-crash
			// picture of what's going on.
			//				      if (vertsAdded % 10000 == 0)
			//				        logMessage(1, "%d  ", m_mesh->getNumVerts());
		}

		assert(
				(m_IQ.qIsQueueEmpty() && m_bdryFacesToSplit.empty())
						|| vertsAdded >= maxVertsAdded);
		ToTime = omp_get_wtime() - ToTime;
		TotalRefineTime += ToTime;
		return vertsAdded;
	}
	assert(0);
	return 0;
}

double RefinementManager2D::FindAltitude(const Vert* pVA, const Vert* pVB,
		double circumcenter[], const Vert* pVC) const {
	/// Find the position of inserted point
	/// altitude

	//////////////////////////////////////////////
	//bisection method
	const double * const adA = pVA->getCoords();
	const double * const adB = pVB->getCoords();
	const double * const adC = pVC->getCoords();

	double midAB[2] = { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };
//	double h0=2.*sqrt(pow(adC[0]-midAB[0],2)+pow(adC[1]-midAB[1],2));
	double h0 = sqrt(dDIST_SQ_2D(circumcenter, midAB));
	double h1 = 0.;
//	printf("circumcenter: %6.4f\t%6.4f\n",circumcenter[0],circumcenter[1]);
	double t = 0.;
	double location[2] = { 0., 0. };
	double distA = 0.;
	double distC = 0.;
	double lsLoc = 0.;
	double answer = 1.e5;
	int b0, b1;
	double m_g = m_mesh->gradeparameter;
	t = (h0) / sqrt(dDIST_SQ_2D(circumcenter, midAB));
	location[0] = midAB[0] + t * (circumcenter[0] - midAB[0]);
	location[1] = midAB[1] + t * (circumcenter[1] - midAB[1]);
	distA = sqrt(dDIST_SQ_2D(location, adA));
	distC = sqrt(dDIST_SQ_2D(location, adC));
	lsLoc = std::min(pVC->getLengthScale() + distC / m_g,
			pVA->getLengthScale() + distA / m_g);
//	lsLoc = pVA->getLengthScale()+distA/m_g;
	lsLoc = std::min(lsLoc, pVB->getLengthScale() + distA / m_g);
	answer = (distA - 0.5 * (pVA->getLengthScale() + lsLoc));
	b0 = (answer > 0) ? 1 : -1;
	t = (h1) / sqrt(dDIST_SQ_2D(circumcenter, midAB));
	location[0] = midAB[0] + t * (circumcenter[0] - midAB[0]);
	location[1] = midAB[1] + t * (circumcenter[1] - midAB[1]);
	distA = sqrt(dDIST_SQ_2D(location, adA));
	distC = sqrt(dDIST_SQ_2D(location, adC));
	lsLoc = std::min(pVC->getLengthScale() + distC / m_g,
			pVA->getLengthScale() + distA / m_g);
//		lsLoc = pVA->getLengthScale()+distA/m_g;
	lsLoc = std::min(lsLoc, pVB->getLengthScale() + distA / m_g);
	answer = (distA - 0.5 * (pVA->getLengthScale() + lsLoc));
	b1 = (answer > 0) ? 1 : -1;
	if (b0 * b1 > 0) {
		return sqrt(dDIST_SQ_2D(circumcenter, midAB));
	}
	assert(b0 * b1 < 0);
	double hc = 0.5 * (h0 + h1);
	while (true) {
		hc = 0.5 * (h0 + h1);
		t = hc / sqrt(dDIST_SQ_2D(circumcenter, midAB));
		location[0] = midAB[0] + t * (circumcenter[0] - midAB[0]);
		location[1] = midAB[1] + t * (circumcenter[1] - midAB[1]);
		distA = sqrt(dDIST_SQ_2D(location, adA));
		distC = sqrt(dDIST_SQ_2D(location, adC));
		lsLoc = std::min(pVC->getLengthScale() + distC / m_g,
				pVA->getLengthScale() + distA / m_g);
//			lsLoc = pVA->getLengthScale()+distA/m_g;
		lsLoc = std::min(lsLoc, pVB->getLengthScale() + distA / m_g);
		answer = (distA - 0.5 * (pVA->getLengthScale() + lsLoc));
		if (fabs(answer) < 1e-6 * distA)
			break;
		else {
			int bc = (answer > 0) ? 1 : -1;
			if (bc * b0 < 0) {
				h1 = hc;
				b1 = bc;
			} else {
				h0 = hc;
				b0 = bc;
			}
		}

	}

	double alti = t * sqrt(dDIST_SQ_2D(circumcenter, midAB));
//
////	bool B0 = iFuzzyComp(start0,)
//	double alti;
	return alti;
}

void RefinementManager2D::InsertionCalculatorForAF(Face* face,
		TriCell* cellToSplit, double insertionPoint[]) const {
	cellToSplit->calcCircumcenter(insertionPoint);
	const double * const adA = face->getVert(0)->getCoords();
	const double * const adB = face->getVert(1)->getCoords();
	double adAB[2] = { 0.5 * (adA[0] + adB[0]), 0.5 * (adA[1] + adB[1]) };
	double grade = m_mesh->gradeparameter;
	//equilateral triangle
	double t1 = face->calcSize() / 2. / tan(30. / 180. * M_PI)
			/ sqrt(dDIST_SQ_2D(insertionPoint, adAB));
	//Find the length scale at circumcenter
	double hcenter = 1.e9;
	for (int ii = 0; ii < cellToSplit->getNumVerts(); ii++) {
		hcenter = std::min(hcenter,
				cellToSplit->calcCircumradius() / grade
						+ cellToSplit->getVert(ii)->getLengthScale());
	}
	//Size optimal triangle
	double alti1 = FindAltitude(face->getVert(0), face->getVert(1),
			insertionPoint, cellToSplit->getOppositeVert(face));
	double alti2 = FindAltitude(face->getVert(1), face->getVert(0),
			insertionPoint, cellToSplit->getOppositeVert(face));
	double alti = std::min(alti1, alti2);
	double t = alti / sqrt(dDIST_SQ_2D(insertionPoint, adAB));
//			double t = 1./2.*(face->getVert(0)->getLengthScale() + face->getVert(1)->getLengthScale())/2.*sqrt(3.)/sqrt(dDIST_SQ_2D(insertionPoint,adAB));
//							if(iFuzzyComp(face->calcSize(),(1. + 1./grade)*cellToSplit->getOppositeVert(face)->getLengthScale()) > - 1.){
//								t = t/(1. + 1./grade);
//							}
//							else
//								t = t*(1. + 1./grade);

	//Shape optimal triangle
	double t2 = sqrt(0.5 / (1 - cos(m_mesh->minangle)))
			* sqrt(dDIST_SQ_2D(adA,adB) / dDIST_SQ_2D(insertionPoint, adAB));
	//t2 supposed to be greater than t1 if minangle < 60
	if (t2 <= t1) {
		printf("circumcenter: %6.4f\t%6.4f\nmidpoint: %6.4f\t%6.4f\n",
				insertionPoint[0], insertionPoint[1], adAB[0], adAB[1]);
		const double * const adC =
				cellToSplit->getOppositeVert(face)->getCoords();
		printf("triangle: %6.4f\t%6.4f\n%6.4f\t%6.4f\n%6.4f\t%6.4f\n", adA[0],
				adA[1], adB[0], adB[1], adC[0], adC[1]);

	}
	assert(t2 > t1);
	double t3 = face->calcSize() * 0.5 * tan(m_mesh->minangle)
			/ sqrt(dDIST_SQ_2D(insertionPoint, adAB));
//	double t3 = face->calcSize()*0.5*tan(M_PI/4.)/sqrt(dDIST_SQ_2D(insertionPoint,adAB));
	if (m_mesh->gradeparameter > 100.)
		t = t1;
	t = iFuzzyComp(t, t1) > -1. ? t1 : t;
	t = iFuzzyComp(t, t2) > -1. ? t2 : t;
	t = iFuzzyComp(t, 1.) > -1. ? 1. : t;
	t = iFuzzyComp(t, t3) > -1. ? t : t3;

	insertionPoint[0] = adAB[0] + t * (insertionPoint[0] - adAB[0]);
	insertionPoint[1] = adAB[1] + t * (insertionPoint[1] - adAB[1]);
}

GR_index_t RefinementManager2D::refineForQualityonfaces(
		const GR_index_t maxMeshVerts, bool Para_Flag, int NQLayer) {
	(void) Para_Flag;
	(void) NQLayer;	//Not used, silence warnings.
	double ToTime = omp_get_wtime();
	// Reza: This else is simply to keep the serial version in place
	{
		ToTime = omp_get_wtime();
		assert(m_refineBdryForSize == false);
		logMessage(1, "Refining for quality\n");
		logMessage(1, "Apparently refining in serial\n");
		GR_index_t vertsAdded = 0;
		GR_index_t maxVertsAdded = maxMeshVerts - m_mesh->getNumVerts();
		// This next line handles the case of unsigned integer wraparound from
		// subtracting to get a negative.  Using GR_sindex_t might overflow instead
		// of underflowing.
		if (maxMeshVerts < m_mesh->getNumVerts())
			maxVertsAdded = 0;

		m_IQ.vCleanQueue();

		while (!m_IQ.qIsQueueEmpty() && vertsAdded < maxVertsAdded) {
			while (!m_bdryFacesToSplit.empty()) {
				InsertionQueueEntry IQE = m_bdryFacesToSplit.front();

				BdryEdgeBase* bdryEdge =
						(dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace()));
				if (!IQE.qStillInMesh()) {
					m_bdryFacesToSplit.pop_front();
					continue;
				}

				// This fails only if there's another edge we should split first.
				bool success = splitBdryEdge(bdryEdge);
				if (success) {
					++vertsAdded;
					m_bdryFacesToSplit.pop_front();
#ifndef NDEBUG
					m_mesh->writeTempMesh();
#endif
					if (vertsAdded >= maxVertsAdded)
						return vertsAdded;
				}

			}
			m_mesh->sendEvents();
			m_IQ.vCleanQueue();

			if (m_IQ.qIsQueueEmpty())
				break;
			// Get the first entity in the queue
			InsertionQueueEntry IQE = m_IQ.IQETopValidEntry();
			assert(IQE.eType() == InsertionQueueEntry::e2DFace);
			Face* face = m_IQ.IQETopValidEntry().pFFace();
			m_IQ.vPopTopEntry();
			m_IQ.setKeyBdryCondition(this->keybdryconditionNum);
			TriCell* cellToSplit = NULL;
			if (face->isDeleted()) {
				continue;
			}
			///Find the cell to split
			if (face->isBdryFace()) {
				if (face->getCell(0)->isBdryCell()) {
					if (m_mesh->DeterGoodTri(
							dynamic_cast<TriCell*>(face->getCell(1)), face)
							== false)
						cellToSplit =
								(dynamic_cast<TriCell*>(face->getCell(1)));
					else {
						continue;
					}
				} else if (face->getCell(1)->isBdryCell()) {
					if (m_mesh->DeterGoodTri(
							dynamic_cast<TriCell*>(face->getCell(0)), face)
							== false)
						cellToSplit =
								(dynamic_cast<TriCell*>(face->getCell(0)));
					else {
						continue;
					}
				}
			} else {
				TriCell* cell = (dynamic_cast<TriCell*>(face->getCell(0)));
				bool tag0 = m_mesh->DeterGoodTri(cell, face);
				cell = (dynamic_cast<TriCell*>(face->getCell(1)));
				bool tag1 = m_mesh->DeterGoodTri(cell, face);
				if (!(tag0 xor tag1)) {
//					m_mesh->sendEvents();
					continue;
				}
				assert(tag0 xor tag1);
				if (tag0 == false)
					cellToSplit = ((dynamic_cast<TriCell*>(face->getCell(0))));
				else
					cellToSplit = ((dynamic_cast<TriCell*>(face->getCell(1))));
			}
			if (cellToSplit->isSmallAngleAffected()) {
				continue;
			}
			assert(!cellToSplit->isSmallAngleAffected());
			///Cell To Split has been found.

			bool success = true;
			assert(cellToSplit);
			assert(m_bdryFacesToSplit.empty());
			double insertionPoint[3];
			//Calculate the insertion point
//			cellToSplit->calcCircumcenter(insertionPoint);
			InsertionCalculatorForAF(face, cellToSplit, insertionPoint);

			//To Check if need to merge two possible steiner points.
//			double possibleInsertion[3][3];
//			for (int i = 0; i < 3; i++){
//				bool lengthcheck = fabs(cellToSplit->getFace(i)->calcSize()-face->calcSize())/face->calcSize() < 1./20.;
//				lengthcheck = lengthcheck and (cellToSplit->getFace(i)->calcSize() == cellToSplit->calcShortestEdgeLength() or face->calcSize()== cellToSplit->calcShortestEdgeLength());
//				if (cellToSplit->getFace(i)!=face and lengthcheck){
//					cellToSplit->calcCircumcenter(possibleInsertion[i]);
//					InsertionCalculatorForAF(cellToSplit->getFace(i),cellToSplit,possibleInsertion[i]);
//				}
//				else{
//					for (int j = 0; j < 3 ; j++){
//						possibleInsertion[i][j] = INT_MAX;
//					}
//				}
//			}
//			//If they are close, merge.
//			int numclose = 0;
//			for (int i = 0; i < 3; i++){
//				if (possibleInsertion[i][0]!=INT_MAX){
//					if(dDIST_SQ_2D(possibleInsertion[i],insertionPoint) < 0.17*0.17*face->calcSize()*face->calcSize()){
//						numclose++;
//						for (int j = 0; j < 3 ; j++){
//							insertionPoint[j]+=possibleInsertion[i][j];
//						}
//					}
//				}
//			}
//			if (numclose > 0){
////				printf("numclose: %d\n",numclose);
//				for (int j = 0; j < 3 ; j++){
//					insertionPoint[j]/=(numclose+1);
//				}
//			}

			InsertionPointCalculator::ePointType pointType = m_IPC->eAF;

			//Compute watson hull.
			assert(cellToSplit->circumscribesPoint(insertionPoint));
			m_WI->computeHull(insertionPoint, cellToSplit);
			//If point is encroaching, queue all encroached edges and return.
			if (m_WI->areBFacesEncroached()) {
				updateEncroachedBFaceQueue();
				success = false;
			}

			if (!m_WI->isHullValid()) {
//				continue;
				cellToSplit->calcCentroid(insertionPoint);
//				cellToSplit->calcCircumcenter(insertionPoint);
				m_WI->computeHull(insertionPoint, cellToSplit);
				if (m_WI->areBFacesEncroached()) {
					updateEncroachedBFaceQueue();
					success = false;
				}
			}
			if (m_WI->isHullValid()) {
				if (m_WI->minLengthOfNewFaces(insertionPoint)
						< 0.7 * face->calcSize()) {
					//reject
					success = false;
				}
			}
			// InsertEntity

			if (success) {
				Vert* newVert = InsertEntity(pointType, vertsAdded);
				assert(newVert);
				m_mesh->sendEvents();

			}

		}

		assert(
				(m_IQ.qIsQueueEmpty() && m_bdryFacesToSplit.empty())
						|| vertsAdded >= maxVertsAdded);
		ToTime = omp_get_wtime() - ToTime;
		TotalRefineTime += ToTime;
		return vertsAdded;
	}
	assert(0);
	return 0;
}

GR_index_t RefinementManager2D::refineBoundaryForLength(
		const GR_index_t maxMeshVerts) {

	fixBoundaryEncroachment();
	m_refineBdryForSize = true;

	// new use of the flag, to make sure all new boundary faces are queued

	GR_index_t vertsAdded = 0;
	GR_index_t maxVertsAdded = maxMeshVerts - m_mesh->getNumVerts();
	// This next line handles the case of unsigned integer wraparound from
	// subtracting to get a negative.  Using GR_sindex_t might overflow instead
	// of underflowing.
	if (maxMeshVerts < m_mesh->getNumVerts())
		maxVertsAdded = 0;

	if (!m_sizingField)
		return 0;

	BdryEdgeBase* bdryEdge = NULL;

	//for(int i = m_mesh->iNumBdryFaces() - 1; i >= 0; i--) {
	GR_index_t i;
	for (i = 0; i < m_mesh->getNumTotalBdryFaces(); i++) {
		bdryEdge = dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		assert(bdryEdge);
		m_bdryFacesToSplit.push_back(InsertionQueueEntry(bdryEdge));

	}

	if (m_bdryFacesToSplit.empty())
		return 0;

	logMessage(1, "Refining boundary for length\n");

	while (!m_bdryFacesToSplit.empty()) {

		InsertionQueueEntry IQE = m_bdryFacesToSplit.front();
		bdryEdge = dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace());
		if (!IQE.qStillInMesh()) {
			m_bdryFacesToSplit.pop_front();
			continue;
		}
		//if the circumradius > average length scale,
		//or if the boundary edge is encroached.
		if (bdryEdge->calcLength()
				> (bdryEdge->getVert(0)->getLengthScale()
						+ bdryEdge->getVert(1)->getLengthScale())
				|| bdryEdge->isEncroached(m_mesh->getEncroachmentType())
						!= eClean) {
			if (splitBdryEdge(bdryEdge)) {
				vertsAdded++;
				m_bdryFacesToSplit.pop_front();
				m_mesh->sendEvents();
				if (vertsAdded >= maxVertsAdded)
					return vertsAdded;
			}
			//			assert(m_mesh->isValid());
		}

		else
			m_bdryFacesToSplit.pop_front();
	}
#ifndef NDEBUG
	//There should not be any encroached edges here.
	queueEncroachedBdryEdges();
	assert(m_bdryFacesToSplit.empty());
	m_bdryFacesToSplit.clear();
#endif
	//set this back to the status quo
	m_refineBdryForSize = false;

	return vertsAdded;
}

GR_index_t RefinementManager2D::refineBoundaryForLengthGeom(
		const GR_index_t maxMeshVerts) {
//	for(int i = 0; i < m_mesh->getNumVerts(); i++){
//		Vert * vert = m_mesh->getVert(i);
//		vert->printVertInfo();
//		printf("Length scale %f\n",vert->getLengthScale());
//	}

	fixBoundaryEncroachment();
	m_refineBdryForSize = true;
	GR_index_t vertsAdded = 0;
	GR_index_t maxVertsAdded = maxMeshVerts - m_mesh->getNumVerts();
	// This next line handles the case of unsigned integer wraparound from
	// subtracting to get a negative.  Using GR_sindex_t might overflow instead
	// of underflowing.
	if (maxMeshVerts < m_mesh->getNumVerts())
		maxVertsAdded = 0;

	if (!m_sizingField)
		return 0;

	BdryEdgeBase* bdryEdge = NULL;

	//for(int i = m_mesh->iNumBdryFaces() - 1; i >= 0; i--) {
	GR_index_t i;
	double minLengthScale = LARGE_DBL;
	for (i = 0; i < m_mesh->getNumTotalBdryFaces(); i++) {
		bdryEdge = dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		assert(bdryEdge);
		Vert * v0 = bdryEdge->getVert(0);
		Vert * v1 = bdryEdge->getVert(1);
//		v0->setLengthScale(fmin(v0->getLengthScale(), v1->getLengthScale()));
//		v1->setLengthScale(fmin(v0->getLengthScale(), v1->getLengthScale()));
		minLengthScale = fmin(v0->getLengthScale(), minLengthScale);
		minLengthScale = fmin(minLengthScale, v1->getLengthScale());

	}
	minLengthScale *= 1.0;
	for (i = 0; i < m_mesh->getNumTotalBdryFaces(); i++) {
		bdryEdge = dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		assert(bdryEdge);
		Vert * v0 = bdryEdge->getVert(0);
		Vert * v1 = bdryEdge->getVert(1);
//		v0->setLengthScale(fmin(v0->getLengthScale(), v1->getLengthScale()));
//		v1->setLengthScale(fmin(v0->getLengthScale(), v1->getLengthScale()));
		v0->setLengthScale(minLengthScale);
		v1->setLengthScale(minLengthScale);
//		v0->printVertInfo();
//		printf("Length scale %f\n",v0->getLengthScale());
//		v1->printVertInfo();
//		printf("Length scale %f\n",v1->getLengthScale());
	}
	for (i = 0; i < m_mesh->getNumTotalBdryFaces(); i++) {
		bdryEdge = dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		assert(bdryEdge);
		m_bdryFacesToSplit.push_back(InsertionQueueEntry(bdryEdge));

	}

	if (m_bdryFacesToSplit.empty())
		return 0;

	logMessage(1, "Refining boundary for length\n");

	while (!m_bdryFacesToSplit.empty()) {

		InsertionQueueEntry IQE = m_bdryFacesToSplit.front();
		bdryEdge = dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace());
		if (!IQE.qStillInMesh()) {
			m_bdryFacesToSplit.pop_front();
			continue;
		}
		//half the circumradius > average length scale,
		//or if the boundary edge is encroached.
		if (bdryEdge->calcLength()
				> bdryEdge->getVert(0)->getLengthScale()
						+ bdryEdge->getVert(1)->getLengthScale()
				|| bdryEdge->isEncroached(m_mesh->getEncroachmentType())
						!= eClean) {

			if (splitBdryEdge(bdryEdge)) {
				vertsAdded++;
				m_bdryFacesToSplit.pop_front();
				m_mesh->sendEvents();
				if (vertsAdded >= maxVertsAdded)
					return vertsAdded;
			}
			//			assert(m_mesh->isValid());
		}

		else
			m_bdryFacesToSplit.pop_front();
	}
//#ifndef NDEBUG
//	//There should not be any encroached edges here.
//	queueEncroachedBdryEdges();
//	assert(m_bdryFacesToSplit.empty());
//	m_bdryFacesToSplit.clear();
//#endif
	m_refineBdryForSize = false;
	return vertsAdded;
}

static void removeVertFromMesh2D(Mesh2D* const mesh, Vert* vertex) {
	logMessage(3, "Deleting vertex at %f %f\n", vertex->x(), vertex->y());
	assert(!vertex->isDeleted());
	assert(vertex->getVertType() == Vert::eInterior);

	Face* face = vertex->getHintFace();

	assert(face->getFaceLoc() == Face::eInterior);
	assert(face->getLeftCell()->getType() == Cell::eTriCell);

	int region = face->getLeftCell()->getRegion();
	int num_cells_before = mesh->getNumCells();
	int num_swaps;

#ifndef NDEBUG
	bool success = mesh->removeVert(vertex, num_swaps);
	assert(success);
#else
	mesh->removeVert(vertex, num_swaps);
#endif

	int num_cells_after = mesh->getNumCells();

	for (int i = num_cells_before; i < num_cells_after; i++) {

		Cell* cell = mesh->getCell(i);
		assert(!cell->isDeleted());
		assert(cell->getType() == Cell::eTriCell);

		cell->setRegion(region);
		assert(cell->getRegion() > 0);
	}
}

void RefinementManager2D::deleteVertsInEdgeBall(
		BdryEdgeBase* const edge_to_split) {
	double radius = 0.5 * edge_to_split->calcLength();
	CubitVector midpoint = edge_to_split->calcMidpoint();
	CubitVector vert_coord(0., 0., 0.);

	queue<Cell*> cells_to_test;
	set<Cell*> cells_tested;
	set<Vert*> verts_tested;
	set<Vert*> verts_to_delete;

	int i, j;

	for (i = 0; i < edge_to_split->getNumFaces(); i++)
		cells_to_test.push(
				edge_to_split->getFace(i)->getOppositeCell(edge_to_split));

	while (!cells_to_test.empty()) {
		Cell* cell = cells_to_test.front();
		cells_to_test.pop();

		assert(!cell->isDeleted());
		if (cell->getType() != Cell::eTriCell)
			continue;
		if (!cells_tested.insert(cell).second)
			continue;

		for (i = 0; i < 3; i++) {
			Vert* vertex = cell->getVert(i);

			if (vertex->getVertType() != Vert::eInterior)
				continue;

			if (verts_tested.insert(vertex).second) {
				vert_coord.set(vertex->x(), vertex->y(), 0.);

				if (iFuzzyComp((midpoint - vert_coord).length(), radius) != 1) {
					verts_to_delete.insert(vertex);

					for (j = 0; j < 3; j++)
						cells_to_test.push(
								cell->getFace(j)->getOppositeCell(cell));
				}
			}
		}
	}
	//No point in continuing if there is no vert to delete.
	if (verts_to_delete.empty())
		return;

	for (Vert *vert : verts_to_delete) {
		removeVertFromMesh2D(m_mesh, vert);
	}
}

#ifdef HAVE_MESQUITE
void AnisoEdgeRefinementManager2D::SetQueuesPriority(enum eInsertType Type) {

	switch (Type) {

	case eQualityEdge: //Edge Insertion based on quality

//	m_IQBF.setPriority(InsertionQueue::eBestImpr);
		m_IQBF.setPriority(InsertionQueue::eLongestEdge); /*In IQBF,give highest priority to longest edge(for edge insertion)...IQBF means Insertion queue boundary faces*/
		m_IQBF.setIsoQMetric(GRUMMP::TMOPQual::eShape);
		m_IQBF.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

		m_IQIF.setPriority(InsertionQueue::eBestImpr); /*In IQIF,give highest priority to edge that gives best improvement(for edge insertion)...IQIF means Insertion queue internal faces*/
		m_IQIF.setIsoQMetric(GRUMMP::TMOPQual::eShape);
		m_IQIF.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

		SetReceiveCreated(true);
		break;

	case eMetricLengthEdge: //Edge Insertion based on a max metric length

		m_IQBF.setPriority(InsertionQueue::eLongestEdge); /*In IQBF,give highest priority to longest edge (for edge insertion)*/
		m_IQBF.setIsoQMetric(GRUMMP::TMOPQual::eShape);
		m_IQBF.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

		m_IQIF.setPriority(InsertionQueue::eLongestEdge); /*In IQIF ,give highest priority to longest edge (for edge insertion)*/
		m_IQIF.setIsoQMetric(GRUMMP::TMOPQual::eShape);
		m_IQIF.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

		SetReceiveCreated(false);
		break;

	case eCoarsening: //Coarsening based on min metric length

		m_IQBF.setPriority(InsertionQueue::eShortestEdge); /*In IQBF while coarsening,give highest priority to shortest edge*/
		m_IQBF.setIsoQMetric(GRUMMP::TMOPQual::eShape);
		m_IQBF.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

//		  m_IQBF.setPriority(InsertionQueue::eInvalid);

		m_IQIF.setPriority(InsertionQueue::eShortestEdge); /*In IQIF while coarsening,give highest priority to shortest edge*/
		m_IQIF.setIsoQMetric(GRUMMP::TMOPQual::eShape);
		m_IQIF.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

		SetReceiveCreated(true);
		break;
	}
}
/*since edges are inserted so new faces and boundary faces are create and that is what follows next*/
void AnisoEdgeRefinementManager2D::receiveCreatedFaces(std::vector<Face*>& createdFaces) {
	logMessage(4,"AnisoRefinementManager2D received word about %zu newly created faces.\n",createdFaces.size());
	if (bReceiveCreated == false)
		return;

	std::vector<Face*>::iterator iter = createdFaces.begin(), iterEnd =
			createdFaces.end();
	for (; iter != iterEnd; ++iter) {
		Face* pF = dynamic_cast<Face*>(*iter);
		if (!pF->doFullCheck() || pF->getFaceLoc() != Face::eInterior)
			continue;
		addToPriorityQueue(pF);
	}
}

void AnisoEdgeRefinementManager2D::receiveCreatedBFaces(std::vector<BFace*>& createdBFaces) {
	logMessage(4,"AnisoRefinementManager2D received word about %zu newly created bfaces.\n",createdBFaces.size());
	if (bReceiveCreated == false)
		return;

	std::vector<BFace*>::iterator iter = createdBFaces.begin(), iterEnd =
			createdBFaces.end();
	for (; iter != iterEnd; ++iter) {
		BdryEdgeBase* bdryEdge = dynamic_cast<BdryEdgeBase*>(*iter);
		if (!bdryEdge->doFullCheck())
			continue;
		addToPriorityQueue(bdryEdge);
	}
}

bool AnisoEdgeRefinementManager2D::splitFace(Face* const FaceToSplit) {

	std::set<Vert*> setpV;
	Vert* pNVert;
	double insertionPoint[2];

	for (int i = 0; i < FaceToSplit->getNumVerts(); i++)
		setpV.insert(FaceToSplit->getVert(i));

	m_IPC->SetIsoQM(m_IQIF.getIsoQMetric());
	m_IPC->SetAnisoQM(m_IQIF.getAnisoQMetric());
	InsertionPointCalculator::ePointType pointType = m_IPC->calcInsertionPoint(
			FaceToSplit, insertionPoint);
	assert(pointType == InsertionPointCalculator::eAnisoEdge);

	pNVert = m_SI2D->insertPointOnFace(insertionPoint, FaceToSplit);

	InterpMetric(pNVert, setpV, m_Type);
	pNVert->setType(Vert::eInterior);



	return true;
}

bool AnisoEdgeRefinementManager2D::splitBdryEdge(BdryEdgeBase* bdryEdgeToSplit) {

	assert(bdryEdgeToSplit);
	assert(!bdryEdgeToSplit->isDeleted());

	std::set<Vert*> setpV;
	for (int i = 0; i < bdryEdgeToSplit->getNumVerts(); i++)
		setpV.insert(bdryEdgeToSplit->getVert(i));

	Face *pF = bdryEdgeToSplit->getFace(0);
	assert(pF->getFaceLoc() == Face::eBdryFace);
	double insertionPoint[2];
	m_IPC->SetIsoQM(m_IQBF.getIsoQMetric());
	m_IPC->SetAnisoQM(m_IQBF.getAnisoQMetric());
	InsertionPointCalculator::ePointType pointType = m_IPC->calcInsertionPoint(
			bdryEdgeToSplit, insertionPoint);
	assert(pointType == InsertionPointCalculator::eAnisoEdge);

	Vert* pNVert = m_SI2D->insertPointOnBdryFace(insertionPoint, pF);

	InterpMetric(pNVert, setpV, m_Type);
	pNVert->setType(Vert::eBdry);

	return true;
}

bool AnisoEdgeRefinementManager2D::addToPriorityQueue(Face* const pF) /*Add face to priority queue*/
{

	if (!pF->doFullCheck() || pF->getEntTopology() != iMesh_LINE_SEGMENT|| pF->getType() != Face::eEdgeFace || pF->isBdryFace()|| pF->getFaceLoc() == Face::ePseudoSurface)
		return false;
	return m_IQIF.qAddEntryAnisoEdge(InsertionQueueEntry(pF));
}

bool AnisoEdgeRefinementManager2D::addToPriorityQueue(BdryEdgeBase* const pBE) /*Add boundary face to priority queue*/
{

	if (pBE->isDeleted())
		return false;
	return m_IQBF.qAddEntryAnisoEdge(InsertionQueueEntry(pBE));
}

GR_index_t AnisoEdgeRefinementManager2D::queueAll() {

	GR_index_t FAdded = 0;
	GR_index_t BFAdded = 0;


	FAdded = queueFaces();
	BFAdded = queueBFaces();
	assert(FAdded + BFAdded<= (this->m_IQIF.iQueueLength()+ this->m_IQBF.iQueueLength()));
	return FAdded + BFAdded;
}

GR_index_t AnisoEdgeRefinementManager2D::queueFaces() {
	GR_index_t FacesAdded = 0;

	for (GR_index_t i = 0; i < m_mesh->getNumFaces(); i++) {
		Face *pF = m_mesh->getFace(i);
		if (!pF->doFullCheck() || pF->getFaceLoc() == Face::eBdryFace|| pF->getFaceLoc() == Face::ePseudoSurface)
			continue;
		if (addToPriorityQueue(pF))
			FacesAdded++;
	}
	assert(FacesAdded <= (this->m_IQIF.iQueueLength()));
	return FacesAdded;
}

GR_index_t AnisoEdgeRefinementManager2D::queueBFaces() {
	GR_index_t BFacesAdded = 0;
	for (GR_index_t i = 0; i < m_mesh->getNumBdryFaces(); i++) {
		BdryEdgeBase* bdryEdge = dynamic_cast<BdryEdgeBase*>(m_mesh->getBFace(i));
		assert(bdryEdge);
		if (!bdryEdge->doFullCheck())
			continue;
		if (addToPriorityQueue(bdryEdge)) {
			BFacesAdded++;
		}
	}
	assert(BFacesAdded <= (this->m_IQBF.iQueueLength()));
	return BFacesAdded;

}

GR_index_t AnisoEdgeRefinementManager2D::AnisoRefineMesh() {

	GR_index_t Added;
	GR_index_t Removed;
	bool bIterate;
	double dConvergePerc = 0.50;

	do {
		bIterate = false;
		Added = 0;
		Removed = 0;
		int iMeshSize = m_mesh->getNumVerts();

		Added += AnisoQualityEdgeInsertion();
		if (m_useGlobalSmoothing) {
			GlobalSmoothing();
		} else {
			LocalSmoothing(3);
		}

		m_pSM2D->clearQueue();


		m_pSM2D->swapAllFaces();

		m_pSM2D->clearQueue();

		Added += AnisoQualityEdgeInsertion(); /*during edge insertion, we add edges to get best improvement....this one based on quality*/

		Added += AnisoMetricLengthEdgeInsertion(); /*during edge insertion, we add edges to get best improvement...this one based on metric length*/

		Removed += AnisoMetricEdgeCoarsening(); /*while coarsening, we remove edges (remember those smaller than particular desired length*/

		m_pSM2D->clearQueue();
		m_pSM2D->swapAllFaces();
		m_pSM2D->clearQueue();

		if ((Removed > iMeshSize * dConvergePerc)|| (Added > iMeshSize * dConvergePerc))
			bIterate = true;

	} while (bIterate == true);

	AnisoQualityEdgeInsertion();
	return 0;
}

void AnisoEdgeRefinementManager2D::LocalSmoothing(int iIter) {

	m_pSMM->vSetTargetfromMetric();

	for (int i = 0; i < iIter; i++) {
		for (unsigned int j = 0; j < m_mesh->getNumVerts(); j++) {
			Vert *pVSmooth = m_mesh->getVert(j);
			if (!pVSmooth->isValid() || pVSmooth->isDeleted())
				continue;
			m_pSMM->smoothVert(pVSmooth);
		}
	}
}

void AnisoEdgeRefinementManager2D::GlobalSmoothing() {
	m_mesh->purgeAllEntities();
	GRUMMP::SmoothingManagerMesquite2D *pSMM2Da =new GRUMMP::SmoothingManagerMesquite2D(m_mesh);
	pSMM2Da->vSetTargetfromMetric();

	pSMM2Da->smoothAllVerts(1);
	delete pSMM2Da;
}

GR_index_t AnisoEdgeRefinementManager2D::AnisoQualityEdgeInsertion() {

	GR_index_t Add = 0;
	GR_index_t TotalAdd = 0;
	SetQueuesPriority(eQualityEdge);

	do {
		Add = 0;
		m_pSM2D->clearQueue();
		m_IQIF.emptyQueue();
		m_IQBF.emptyQueue();
		queueAll();

		Add += AnisoRefineAllQueued();
		assert(m_IQIF.qIsQueueEmpty());
		assert(m_IQBF.qIsQueueEmpty());

		{
			m_pSM2D->swapAllQueuedFacesRecursive();
		}
		continue;
		m_pSM2D->clearQueue();

		TotalAdd += Add;
		assert(m_mesh->isValid());
	} while (Add > 100);

	return TotalAdd;
}

GR_index_t AnisoEdgeRefinementManager2D::AnisoMetricLengthEdgeInsertion() {

	GR_index_t TotalAdd = 0;
	SetQueuesPriority(eMetricLengthEdge);
	m_pSM2D->clearQueue();
	m_IQIF.emptyQueue();
	m_IQBF.emptyQueue();

	queueAll();

	TotalAdd = AnisoRefineAllQueued();

	m_pSM2D->swapAllQueuedFacesRecursive();
	m_pSM2D->clearQueue();

	return TotalAdd;
}

GR_index_t AnisoEdgeRefinementManager2D::AnisoMetricEdgeCoarsening() {

	GR_index_t Removed = 0;
	GR_index_t TotalRemoved = 0;
	SetQueuesPriority(eCoarsening);
	do {
		Removed = 0;
		m_pSM2D->clearQueue();
		m_IQIF.emptyQueue();
		m_IQBF.emptyQueue();
		for (unsigned i = 0; i < m_mesh->getNumVerts(); i++) {

			m_mesh->getVert(i)->markToKeep();

		}

		queueAll();
		Removed = AnisoRefineAllQueued();
		TotalRemoved += Removed;
		assert(m_mesh->isValid());
	}

	while (Removed > 1500);

	return TotalRemoved;

}

GR_index_t AnisoEdgeRefinementManager2D::AnisoRefineAllQueued() {

	logMessage(4, "%d Bdry faces queued for refinement.\n",
			m_IQBF.iQueueLength());
	logMessage(4, "%d Int faces queued for refinement.\n",
			m_IQIF.iQueueLength());

	GR_index_t VertsAdded = 0;

	do {
		m_IQBF.vCleanQueue();
		while (!m_IQBF.qIsQueueEmpty()) {
			InsertionQueueEntry IQE = m_IQBF.IQETopValidEntry();
			assert(IQE.eType() == InsertionQueueEntry::eBdryEdge);
			BdryEdgeBase* bdryEdge = dynamic_cast<BdryEdgeBase*>(IQE.pBFBFace());

			if (!IQE.qStillInMesh()) {
				m_IQBF.vPopTopEntry();
				continue;
			}

			assert(bdryEdge->doFullCheck());
			bool success;
			if (m_IQBF.getPriority() == InsertionQueue::eShortestEdge) {
				assert(m_mesh->isValid());
				Face *pF = IQE.pBFBFace()->getFace(0);
				success = m_pAR->bRemoveVertexwithQual(pF);
			} else
				success = splitBdryEdge(bdryEdge);

			if (success)
				++VertsAdded;
			m_IQBF.vPopTopEntry();
			m_mesh->sendEvents();

		}

		m_IQIF.vCleanQueue();

		while (!m_IQIF.qIsQueueEmpty()) {

			InsertionQueueEntry IQE = m_IQIF.IQETopValidEntry();
			assert(IQE.eType() == InsertionQueueEntry::e2DFace);
			assert(m_mesh->isValid());

			assert(IQE.pFFace()->isValid());
			Face *pF = IQE.pFFace();
			if (!IQE.qStillInMesh() || !IQE.pFFace()->doFullCheck()|| pF->getFaceLoc()== Face::eBdryFace) {
				m_IQIF.vPopTopEntry();
				continue;
			}

			bool success = false;
			if (m_IQIF.getPriority() == InsertionQueue::eShortestEdge ) {
				assert(m_mesh->isValid());
				success = m_pAR->bRemoveVertexwithQual(pF);
			}
			else
				success = splitFace(dynamic_cast<Face*>(IQE.pFFace()));

			if (success == true)
				++VertsAdded;
			m_mesh->sendEvents();
			m_IQIF.vPopTopEntry();
		}

	}

	while (!m_IQBF.qIsQueueEmpty() || !m_IQIF.qIsQueueEmpty());
	return VertsAdded;
}

void AnisoEdgeRefinementManager2D::LinearInterpMetric(Vert* pV,std::set<Vert*> setpV) {

	std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();
	double dDist[setpV.size()];
	double dTDist = 0;
	double adNMetric[iNMetricValues];
	GR_index_t i = 0;

	for (; iter != iterEnd; ++iter) {
		Vert* psetV = (*iter);
		assert(psetV->isValid());
		dDist[i] = calcDistanceBetween(pV, psetV);
		dTDist += dDist[i];
		i++;
	}

	for (int ii = 0; ii < iNMetricValues; ii++)
		adNMetric[ii] = 0;

	iter = setpV.begin();
	iterEnd = setpV.end();

	i = 0;
	for (; iter != iterEnd; ++iter) {
		Vert* psetV = (*iter);
		for (int ii = 0; ii < iNMetricValues; ii++) {
			adNMetric[ii] += (dDist[i] / dTDist) * psetV->getMetric(ii);
		}
		i++;
	}

	pV->setMetric(adNMetric);

	return;
}

void AnisoEdgeRefinementManager2D::AverageMetric(Vert* pV,std::set<Vert*> setpV) {

	std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();
	double adNMetric[iNMetricValues];

	for (int ii = 0; ii < iNMetricValues; ii++)
		adNMetric[ii] = 0;

	for (; iter != iterEnd; ++iter) {
		Vert* psetV = (*iter);
		for (int ii = 0; ii < iNMetricValues; ii++)
			adNMetric[ii] += (1. / setpV.size()) * psetV->getMetric(ii);
	}
	pV->setMetric(adNMetric);
	return;
}

void AnisoEdgeRefinementManager2D::InterpMetric(Vert* pV, std::set<Vert*> setpV,enum eInterpType Type) {

	switch (Type) {
	case eDefault:
		AverageMetric(pV, setpV);
		break;
	case eAverage:
		AverageMetric(pV, setpV);
		break;
	case eLinear:
		LinearInterpMetric(pV, setpV);
		break;
	case eLogEuclidean:
		m_TQ->LogEuclideanInterpMetric(pV, setpV);
		break;
	}
	return;
}
#endif
}
