#include "iMesh_unitTest.hh"

void
vTestMeshWithSets(iMesh_Instance& instance)
{
  //    getNumOfType, getNumOfTopo, getAllVtxCoords,
  //    getVtxCoordIndex, getEntities, and getAdjEntities
  const int iNumType = iBase_ALL_TYPES - iBase_VERTEX;
  const int iNumTopo = iMesh_ALL_TOPOLOGIES - iMesh_POINT;
  iBase_EntitySetHandle apvSets[5];
  iBase_EntityHandle a2pvTypeEnts[iNumType][iNumType * 5];
  iBase_EntityIterator pvIter;
  iBase_EntityHandle a2pvTopoEnts[iNumTopo][iNumTopo * 5];
  int aiNumTypeEnts[iNumType] =
    { 0 };
  int aiNumTopoEnts[iNumTopo] =
    { 0 };
  int i, iGeoDim, qGotData = false, err;

  iMesh_getGeometricDimension(instance, &iGeoDim, &err);
  CHECK_ERR(err);

  // Get a different number of entities for each type and record how
  // many were retrieved --a2pvTypeEnts stores the entities,
  // aiNumTypeEnts stores how many for each type
  qGotData = false;

  for (i = 0; i < iNumType; i++) {
    if ((iGeoDim == 2) && (i > iBase_FACE))
      break;
    iMesh_initEntIter(instance, pvRootSet,
		      static_cast<iBase_EntityType>(i + iBase_VERTEX),
		      iMesh_ALL_TOPOLOGIES, 0, &pvIter, &err);
    CHECK_ERR(err);

    iBase_EntityHandle entity;
    iMesh_getNextEntIter(instance, pvIter, &entity, &qGotData, &err);
    CHECK_ERR(err);

    int j = 0;
    while (qGotData && aiNumTypeEnts[i] < 5 * (iNumType - i)) {
      a2pvTypeEnts[i][j] = entity;
      j++;
      aiNumTypeEnts[i] = j;
      iMesh_getNextEntIter(instance, pvIter, &entity, &qGotData, &err);
      CHECK_ERR(err);
    };
    iMesh_endEntIter(instance, pvIter, &err);
    CHECK_ERR(err);
  }

  // If there aren't the required number of vertices (20 to be exact) or
  // there aren't any of any of the other types of entities, bail out
  // because the subsequent tests would be useless and might not even
  // work.
  if ((aiNumTypeEnts[0] < 20)) {
    cout << "\nCould not get mesh data to work with... exiting...\n";
    TEST(data_preparation);
    return;
  }

  // Get a different number of entities for each topology and record how
  // many were retrieved -- a2pvTopoEnts stores the entities,
  // aiNumTopoEnts stores how many for each topology.
  qGotData = false;
  for (i = 0; i < iNumTopo; i++) {
    if ((iGeoDim == 2) && (i > iMesh_QUADRILATERAL))
      break;
    iMesh_initEntIter(instance, pvRootSet, iBase_ALL_TYPES,
		      static_cast<iMesh_EntityTopology>(i), 0, &pvIter, &err);
    CHECK_ERR(err);

    iBase_EntityHandle entity;
    iMesh_getNextEntIter(instance, pvIter, &entity, &qGotData, &err);
    CHECK_ERR(err);
    int j = 0;
    while (qGotData && aiNumTopoEnts[i] < 5 * (iNumTopo - i)) {
      a2pvTopoEnts[i][j] = entity;
      j++;
      aiNumTopoEnts[i] = j;
      iMesh_getNextEntIter(instance, pvIter, &entity, &qGotData, &err);
      CHECK_ERR(err);
    }
    iMesh_endEntIter(instance, pvIter, &err);
    CHECK_ERR(err);
  }

  // If there aren't the required number of points (55 to be exact) or
  // there aren't any of any of the other topologies of entities, bail
  // out because the subsequent tests would be useless and might not
  // even work
  if (aiNumTopoEnts[0] < 55) {
    cout << "\nCould not get mesh data to work with... exiting...\n";
    TEST(data_preparation);
    return;
  }

  // Create two lists and two sets, one for types and one for topologies each,
  // as well as an empty list
  for (i = 0; i < 5; i++) {
    iMesh_createEntSet(instance, i < 3, &apvSets[i], &err);
    CHECK_ERR(err);
  }

  // Add all of the type entities to list 1 and set 3
  for (i = 0; i < iNumType; i++) {
    for (int j = 0; j < aiNumTypeEnts[i]; j++) {
      iMesh_addEntToSet(instance, a2pvTypeEnts[i][j], apvSets[1], &err);
      CHECK_ERR(err);
      iMesh_addEntToSet(instance, a2pvTypeEnts[i][j], apvSets[3], &err);
      CHECK_ERR(err);
    }
  }

  // Double the number of vertices for the list
  for (int j = 0; j < aiNumTypeEnts[0]; j++)
    iMesh_addEntToSet(instance, a2pvTypeEnts[0][j], apvSets[1], &err);
  CHECK_ERR(err);

  // Add all of the topology entities to list 2 and set 4
  for (i = 0; i < iNumTopo; i++) {
    for (int j = 0; j < aiNumTopoEnts[i]; j++) {
      iMesh_addEntToSet(instance, a2pvTopoEnts[i][j], apvSets[2], &err);
      CHECK_ERR(err);
      iMesh_addEntToSet(instance, a2pvTopoEnts[i][j], apvSets[4], &err);
      CHECK_ERR(err);
    }
  }
  // Double the number of points for the list
  for (int j = 0; j < aiNumTopoEnts[0]; j++)
    iMesh_addEntToSet(instance, a2pvTopoEnts[0][j], apvSets[2], &err);
  CHECK_ERR(err);

  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsFromEntCounts(instance, pvRootSet, apvSets, aiNumTypeEnts,
			      aiNumTopoEnts);
  }

  cout << "Testing basic queries...                                  ";

  bool qNumOK0 = true, qNumOK1 = true, qNumOK2 = true, qNumOK3 = true, qNumOK4 =
      true;
  // Ensure that the correct number of each type/topology of entities is
  // returned

  // The vertices have to be handled separately because they were added
  // twice to the list
  qNumOK0 = qNumOK0 && qCheckNumOfType(instance, apvSets[0], iBase_VERTEX, 0);
  qNumOK1 = qNumOK1
      && qCheckNumOfType(instance, apvSets[1], iBase_VERTEX,
			 2 * aiNumTypeEnts[0]);
  qNumOK3 = qNumOK3
      && qCheckNumOfType(instance, apvSets[3], iBase_VERTEX, aiNumTypeEnts[0]);
  for (int iType = iBase_VERTEX + 1; iType < iBase_ALL_TYPES; iType++) {
    if ((iGeoDim == 2) && (iType > iBase_FACE))
      break;
    qNumOK0 = qNumOK0 && qCheckNumOfType(instance, apvSets[0], iType, 0);
    qNumOK1 = qNumOK1
	&& qCheckNumOfType(instance, apvSets[1], iType, aiNumTypeEnts[iType]);
    qNumOK3 = qNumOK3
	&& qCheckNumOfType(instance, apvSets[3], iType, aiNumTypeEnts[iType]);
  }

  // The points have to be handled separately because they were added
  // twice to the list
  qNumOK0 = qNumOK0 && qCheckNumOfTopo(instance, apvSets[0], iMesh_POINT, 0);
  qNumOK2 = qNumOK2
      && qCheckNumOfTopo(instance, apvSets[2], iMesh_POINT,
			 2 * aiNumTopoEnts[0]);
  qNumOK4 = qNumOK4
      && qCheckNumOfTopo(instance, apvSets[4], iMesh_POINT, aiNumTopoEnts[0]);
  for (int iTopo = iMesh_POINT + 1; iTopo < iMesh_ALL_TOPOLOGIES; iTopo++) {
    if ((iGeoDim == 2) && (i > iMesh_QUADRILATERAL))
      break;
    qNumOK0 = qNumOK0 && qCheckNumOfTopo(instance, apvSets[0], iTopo, 0);
    qNumOK2 = qNumOK2
	&& qCheckNumOfTopo(instance, apvSets[2], iTopo, aiNumTopoEnts[iTopo]);
    qNumOK4 = qNumOK4
	&& qCheckNumOfTopo(instance, apvSets[4], iTopo, aiNumTopoEnts[iTopo]);
  }
  TEST(qNumOK0 && qNumOK1 && qNumOK2 && qNumOK3 && qNumOK4);
  cout << "DONE\n";

  cout << "Testing entity retrieval...                               ";
  // Try to retrieve the entities for each type (sets 0, 1, 3)
  // Set #0 is empty

  for (int iType = iBase_VERTEX; iType < iBase_ALL_TYPES; iType++) {
    vCheckEntities(instance, apvSets[0], iType, iMesh_ALL_TOPOLOGIES,
		   a2pvTypeEnts[iType], 0);
    vCheckEntities(instance, apvSets[3], iType, iMesh_ALL_TOPOLOGIES,
		   a2pvTypeEnts[iType], aiNumTypeEnts[iType]);
    {
      iBase_EntityHandle *aEnts = new iBase_EntityHandle[2 * aiNumTypeEnts[0]];
      int iTotal = 0;
      for (; iTotal < aiNumTypeEnts[iType]; iTotal++) {
	aEnts[iTotal] = a2pvTypeEnts[iType][iTotal];
      }
      if (iType == 0) {
	// In this case, add the verts again.
	for (; iTotal < 2 * aiNumTypeEnts[0]; iTotal++) {
	  aEnts[iTotal] = a2pvTypeEnts[0][iTotal - aiNumTypeEnts[0]];
	}
      }
      vCheckEntities(instance, apvSets[1], iType, iMesh_ALL_TOPOLOGIES, aEnts,
		     iTotal);
      delete[] aEnts;
    }
  }

  // Try to retrieve the entities for each topo (sets 0, 2, 4)
  // Set #0 is empty

  for (int iTopo = iMesh_POINT; iTopo < iMesh_ALL_TOPOLOGIES; iTopo++) {
    vCheckEntities(instance, apvSets[0], iBase_ALL_TYPES, iTopo,
		   a2pvTopoEnts[iTopo], 0);
    vCheckEntities(instance, apvSets[4], iBase_ALL_TYPES, iTopo,
		   a2pvTopoEnts[iTopo], aiNumTopoEnts[iTopo]);
    {
      iBase_EntityHandle *aEnts = new iBase_EntityHandle[2 * aiNumTopoEnts[0]];
      int iTotal = 0;
      for (; iTotal < aiNumTopoEnts[iTopo]; iTotal++) {
	aEnts[iTotal] = a2pvTopoEnts[iTopo][iTotal];
      }
      if (iTopo == 0) {
	// In this case, add the verts again.
	for (; iTotal < 2 * aiNumTopoEnts[0]; iTotal++) {
	  aEnts[iTotal] = a2pvTopoEnts[0][iTotal - aiNumTopoEnts[0]];
	}
      }
      vCheckEntities(instance, apvSets[2], iBase_ALL_TYPES, iTopo, aEnts,
		     iTotal);
      delete[] aEnts;
    }
  }
  cout << "DONE\n";
  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsFromEntCounts(instance, pvRootSet, apvSets, aiNumTypeEnts,
			      aiNumTopoEnts);
  }

  cout << "Testing entity adjacency consistency...\n";
  bool qAdjCons, qInCons;

  for (int eEType = iBase_VERTEX; eEType < iBase_ALL_TYPES; eEType++) {
    for (int eETopo = iMesh_POINT; eETopo < iMesh_ALL_TOPOLOGIES; eETopo++) {
      if (a2qLegalTypeAndTopo[eEType][eETopo]) {
	if ((iGeoDim == 2)
	    && (eEType > iBase_FACE || eETopo > iMesh_QUADRILATERAL))
	  continue;
	int iNumEnt;
	iMesh_getNumOfTopo(instance, apvSets[2],
			   static_cast<iMesh_EntityTopology>(eETopo), &iNumEnt,
			   &err);
	CHECK_ERR(err);
	if (iNumEnt == 0)
	  continue;

	cout << "  Type: ";
	cout.width(9);
	cout << astrTypes[eEType] << "      with topology: ";
	cout.width(14);
	cout << astrTopologies[eETopo] << endl;

	// Now test for equivalence in adjacency operations.
	int eAdjType, eLastType = iBase_REGION;
	for (eAdjType = iBase_VERTEX; eAdjType <= eLastType; eAdjType++) {
	  // Skip any that are unsupported.
	  if (a2iAdjTable[eEType][eAdjType] == iBase_UNAVAILABLE)
	    continue;
	  qAdjCons = true;
	  qInCons = true;

	  iBase_EntityHandle * aAdjEnt = NULL, *aEnt = NULL;
	  int * aiOffset = NULL, *aiIndex = NULL;
	  int iNumAdj, iOffsetSize, iIndexSize, iEntSize;
	  int iAdjEntAlloc = 0;
	  int iOffsetAlloc = 0;
	  int iIndexAlloc = 0;
	  int iEntAlloc = 0;

	  // Get adjacencies using the entityset version of the index call.
	  iMesh_getAdjEntIndices(instance, apvSets[2], eEType, eETopo, eAdjType,
				 &aEnt, &iEntAlloc, &iEntSize, &aAdjEnt,
				 &iAdjEntAlloc, &iNumAdj, &aiIndex,
				 &iIndexAlloc, &iIndexSize, &aiOffset,
				 &iOffsetAlloc, &iOffsetSize, &err);
	  CHECK_ERR(err);

	  // Now check that the adjacencies match with the single entity call.
	  cout << "   Checking adjacency for entity type:\t";
	  cout.width(9);
	  cout << astrTypes[eAdjType] << "\t";
	  if (iNumEnt != iOffsetSize - 1) {
	    qAdjCons = false;
	  }
	  else {
	    for (int iEnt = 0; iEnt < iNumEnt && qAdjCons; iEnt++) {
	      int iNLocalAdj;
	      iBase_EntityHandle * aAdjEntHandles = NULL;
	      int iAdjEntHandlesAlloc = 0;
	      iMesh_getEntAdj(instance, aEnt[iEnt], eAdjType, &aAdjEntHandles,
			      &iAdjEntHandlesAlloc, &iNLocalAdj, &err);
	      CHECK_ERR(err);
	      // Number of adjacencies must match.
	      if (iNLocalAdj != (aiOffset[iEnt + 1] - aiOffset[iEnt])) {
		qAdjCons = false;
		continue;
	      }
	      int iOffset = aiOffset[iEnt];
	      for (int iAdj = 0; iAdj < iNLocalAdj; iAdj++) {
		iBase_EntityHandle global_handle = aAdjEnt[aiIndex[iOffset
		    + iAdj]];
		iBase_EntityHandle local_handle = aAdjEntHandles[iAdj];
		qAdjCons = qAdjCons && (global_handle == local_handle);
	      } // Loop over local adjacencies
	      UT_FREE(aAdjEntHandles);
	    } // Loop over entities
	    cout << "\tDONE\n";
	  }
	  UT_FREE(aAdjEnt);
	  UT_FREE(aEnt);
	  UT_FREE(aiOffset);
	  UT_FREE(aiIndex);
	} // Loop over adjacency types
	TEST(qAdjCons);
	TEST(qInCons);
      }
    }
  } // Done with double loop to verify consistency of adjacencies
  vDoSaveLoad(&instance, &pvRootSet);
  if (qSaveLoad) {
    vRecoverSetsFromEntCounts(instance, pvRootSet, apvSets, aiNumTypeEnts,
			      aiNumTopoEnts);
  }

  // getAllVtxCoords, getVtxCoordIndex
//   cout << "Checking consistency of vertex coordinate retrieval...    ";
//   int iNumAdjEnts, iNumVerts, iNumIndices, iNumTopologies,
//     iNumOffsets1, iNumOffsets2, iNumOffsets3,
//     iNumAllCoords, iNumCoords,
//     iNumInSet1, iNumInSet2;
//   bool qCoordsCons, qStorCons, qToposCons;
//   int SOTemp;

//   int SODefault;
//   iMesh_getDfltStorage(instance, &SODefault, &err);
//   CHECK_ERR(err);

//   qStorCons = true;
//   for ( int iSet = 0; iSet < 5 && qStorCons; iSet++ ) {
//     int * aiInSet1;
//     double *adAllCoords;
//     SOTemp = iBase_BLOCKED;
//     int iAllCoordsAlloc = 0;
//     int iInSet1Alloc = 0;
//     iMesh_getAllVtxCoords(instance,  apvSets[iSet],
// 			  &adAllCoords, &iAllCoordsAlloc, &iNumAllCoords,
// 			  &aiInSet1, &iInSet1Alloc, &iNumInSet1,
// 			  SOTemp, &err);
//     CHECK_ERR(err);
//     TEST ( iNumAllCoords == iNumInSet1*iGeoDim );
//     qStorCons = qStorCons && ( SOTemp == SODefault );

//     qInCons = qCoordsCons = ( iNumAllCoords == iNumInSet1*iGeoDim );
//     qToposCons = true;
//     for ( int eEType = iBase_VERTEX;
// 	  eEType <= iBase_ALL_TYPES; eEType++ ) {
//       for ( int eETopo = iMesh_POINT;
// 	    eETopo <= iMesh_ALL_TOPOLOGIES
// 	      && qInCons && qToposCons && qCoordsCons
// 	      ; eETopo++ ) {
// 	// Skip any illegal type-topology combinations
// 	if ( !a2qLegalTypeAndTopo[eEType][eETopo] )
// 	  continue;
// 	if ( (iGeoDim == 2) && (eEType == iBase_REGION ||
// 				(eETopo > iMesh_QUADRILATERAL &&
// 				 eETopo != iMesh_ALL_TOPOLOGIES)) )
// 	  continue;

// 	for ( int eAdjType = iBase_VERTEX;
// 	      eAdjType < iBase_ALL_TYPES; eAdjType++ ) {
// 	  // Skip any that are identities or unsupported.
// 	  if ( eEType == eAdjType ||
// 	       a2iAdjTable[eEType][eAdjType] == iBase_UNAVAILABLE)
// 	    continue;

// 	  iBase_EntityHandle * aAdjEnts;
// 	  int * aiOffsets1, *aiInSet2;
// 	  int iAdjEntsAlloc = 0;
// 	  int iOffsets1Alloc = 0;
// 	  int iInSet2Alloc = 0;
// 	  iMesh_getAdjEntities(instance,  apvSets[iSet], eEType,
// 			       eETopo, eAdjType,
// 			       &aAdjEnts, &iAdjEntsAlloc, &iNumAdjEnts,
// 			       &aiOffsets1, &iOffsets1Alloc, &iNumOffsets1,
// 			       &aiInSet2, &iInSet2Alloc, &iNumInSet2, &err);
// 	  CHECK_ERR(err);

// 	  iBase_EntityHandle *aVerts;
// 	  int * aiOffsets2;
// 	  if ( eAdjType == iBase_VERTEX ) {
// 	    aVerts = aAdjEnts;
// 	    iNumVerts = iNumAdjEnts;
// 	  }
// 	  else {
// 	    int iVertsAlloc = 0;
// 	    int iOffsets2Alloc = 0;
// 	    iMesh_getEntArrAdj(instance,  aAdjEnts, iNumAdjEnts, iBase_VERTEX,
// 			       &aVerts, &iVertsAlloc, &iNumVerts,
// 			       &aiOffsets2, &iOffsets2Alloc, &iNumOffsets2,
// 			       &err);
// 	    CHECK_ERR(err);
// 	  }
// 	  double * adCoords;
// 	  SOTemp = iBase_BLOCKED;
// 	  int iCoordsAlloc = 0;
// 	  iMesh_getVtxArrCoords(instance,  aVerts, iNumVerts, SOTemp,
// 				&adCoords, &iCoordsAlloc, &iNumCoords, &err);
// 	  CHECK_ERR(err);
// 	  qStorCons = qStorCons && ( SOTemp == SODefault );

// 	  int * aiOffsets3, *aiIndices;
// 	  int * aiTopologies;
// 	  int iOffsets3Alloc = 0;
// 	  int iIndicesAlloc = 0;
// 	  int iTopologiesAlloc = 0;
// 	  iMesh_getVtxCoordIndex(instance,  apvSets[iSet], eEType, eETopo,
// 				 eAdjType,
// 				 &aiOffsets3, &iOffsets3Alloc, &iNumOffsets3,
// 				 &aiIndices, &iIndicesAlloc, &iNumIndices,
// 				 &aiTopologies, &iTopologiesAlloc,
// 				 &iNumTopologies, &err);
// 	  CHECK_ERR(err);
// 	  if ( iNumTopologies == iNumAdjEnts ) {
// 	    for ( i = 0; i < iNumAdjEnts; i++ ) {
// 	      int iTopo;
// 	      iMesh_getEntTopo(instance, aAdjEnts[i], &iTopo, &err);
// 	      CHECK_ERR(err);
// 	      qToposCons = qToposCons && ( aiTopologies[(i)] == iTopo );
// 	    }
// 	  }
// 	  else {
// 	    qToposCons = false;
// 	  }
// 	  if ( (iNumOffsets3 == iNumAdjEnts+1 ) &&
// 	       (iNumVerts == iNumIndices) &&
// 	       (iNumCoords == iGeoDim*iNumVerts) ) {
// 	    if ( SOTemp == iBase_BLOCKED )
// 	      cout << "\nUsing blocked storage order...";
// 	    int iBegin = 0;
// 	    for ( i = 0; i < iNumAdjEnts && qCoordsCons; i++ ) {
// 	      if ( aiInSet2[i] == 1 ) {
// 		for ( int j = aiOffsets3[i]; j < aiOffsets3[i+1]; j++ ) {
// 		  int iCurrInd = iGeoDim*aiIndices[j];
// 		  iBegin = iGeoDim*j;
// 		  qCoordsCons = qCoordsCons
// 		    && ( adCoords[iBegin] == adAllCoords[iCurrInd] )
// 		    && ( adCoords[iBegin+1] == adAllCoords[iCurrInd+1] );
// 		  if ( iGeoDim == 3 )
// 		    qCoordsCons = qCoordsCons
// 		      && ( adCoords[iBegin+2] == adAllCoords[iCurrInd+2] );
// 		  int iIsCont;
// 		  iMesh_isEntContained(instance, apvSets[iSet], aVerts[j],
// 				       &iIsCont, &err);
// 		  CHECK_ERR(err);
// 		  qInCons = qInCons
// 		    && ( aiInSet1[iCurrInd/iGeoDim] == iIsCont );
// 		}
// 	      }
// 	      else {
// 		int iIsCont;
// 		iMesh_isEntContained(instance, apvSets[iSet], aAdjEnts[i],
// 				     &iIsCont, &err);
// 		CHECK_ERR(err);
// 		qInCons = qInCons && ( !iIsCont );
// 	      }
// 	    } // Done looping over adjacent ents
// 	  } // Done with if (sizes are okay) ...
// 	  else {
// 	    qCoordsCons = false;
// 	  }
// 	} // Loop over adjacent types
//       } // Loop over topos
//     } // Loop over types
//     TEST ( qCoordsCons );
//     TEST ( qInCons );
//     TEST ( qToposCons );
//   } // Loop over sets
//   TEST ( qStorCons );
//   cout << "DONE\n";

  if (qErrTests) {
    cout << "Testing error throwing...                                 ";
    // Error throwing tests simply try to create a situation where the interface
    // should, according to the spec, throw an error, and verify that:
    // (a) an error is thrown, and (b) the error is the one expected.
    // If an error is not thrown, the reversal of the action is attempted, if
    // it would interfere with further tests.

    // Invalid entity/entity set/tag handle tests usually involve passing in an
    // entity as an entity set or vice versa, etc., but sometimes try more
    // difficult tests such as passing in a pointer to character data as an
    // entity handle.

    int iDum;
    //  Getting number of type for invalid entity type
    iMesh_getNumOfType(instance, apvSets[0], iBase_ALL_TYPES + 1, &iDum, &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_TYPE);

    // Getting number of topology for invalid entity topology
    iMesh_getNumOfTopo(instance, apvSets[0], iMesh_ALL_TOPOLOGIES + 1, &iDum,
		       &err);
    CHECK_ERR2(err, iBase_INVALID_ENTITY_TOPOLOGY);

    cout << "DONE\n";
  }

  for (i = 0; i < 5; i++) {
    iMesh_destroyEntSet(instance, apvSets[i], &err);
    CHECK_ERR(err);
  }
}
