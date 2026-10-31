#include "iMesh_unitTest.hh"
#include <limits>

//
// Test that both 1 and 2D mesh's coordinates are 'handled'
// correctly in that nulls may be passed to set/get them,
// and arrays of them returned return correct number.
//
void
vTestLowerDimensionedCoordinates()
{
  bool supports1D = false;
  bool supports2D = false;
  for (int gdim = 1; gdim < 3; gdim++) {
    iMesh_Instance instance;
    int err;

    // Setup a mesh of the desired geometric dimension
    iMesh_newMesh("", &instance, &err, 0);
    iMesh_setGeometricDimension(instance, gdim, &err);
    if (err != iBase_SUCCESS) {
      cout << "Skipping geometric dimension " << gdim << " no support in impl."
	  << endl;
      iMesh_dtor(instance, &err);
      continue;
    }
    if (gdim == 1)
      supports1D = true;
    if (gdim == 2)
      supports2D = true;

    //
    // We test setting bogus coordinate data via the 
    // individual Vtx calls using two different bogus values.
    // If we wind up getting the same bogus values back from
    // the implementation both times, it means the implementation
    // actually did attempt to read and use the values that were
    // passed to it and that is considered an error.
    //
    bool createReadBadCoords[2] =
      { false, false };
    bool setReadBadCoords[2] =
      { false, false };
    for (int pass = 0; pass < 2; pass++) {
      // Create a vtx with specific values in 'extra' coords'
      double oneCoordA[3];
      double badVal =
	  (pass == 0) ?
	      std::numeric_limits<double>::quiet_NaN() :
	      std::numeric_limits<double>::infinity();
      for (int i = 0; i < gdim; i++)
	oneCoordA[i] = 0.0;
      for (int i = gdim; i < 3; i++)
	oneCoordA[i] = badVal;
      iBase_EntityHandle oneVert;
      iMesh_createVtx(instance, oneCoordA[0], oneCoordA[1], oneCoordA[2],
		      &oneVert, &err);
      CHECK_ERR(err);

      // Read back the values created
      double oneCoordB[3];
      iMesh_getVtxCoord(instance, oneVert, &oneCoordB[0], &oneCoordB[1],
			&oneCoordB[2], &err);
      CHECK_ERR(err);
      for (int i = gdim; i < 3; i++) {
	if (oneCoordB[i] == badVal)
	  createReadBadCoords[pass] = true;
      }

      // Ok now set all to zero and confirm they are all zero.
      iMesh_setVtxCoord(instance, oneVert, 0.0, 0.0, 0.0, &err);
      CHECK_ERR(err);
      iMesh_getVtxCoord(instance, oneVert, &oneCoordB[0], &oneCoordB[1],
			&oneCoordB[2], &err);
      CHECK_ERR(err);
      TEST(oneCoordB[0] == 0 && oneCoordB[1] == 0 && oneCoordB[2] == 0);

      // Now, use setVtx to set bad values
      iMesh_setVtxCoord(instance, oneVert, oneCoordA[0], oneCoordA[1],
			oneCoordA[2], &err);
      CHECK_ERR(err);

      // Read back values set
      iMesh_getVtxCoord(instance, oneVert, &oneCoordB[0], &oneCoordB[1],
			&oneCoordB[2], &err);
      CHECK_ERR(err);
      for (int i = gdim; i < 3; i++) {
	if (oneCoordB[i] == badVal)
	  setReadBadCoords[pass] = true;
      }
    }
    TEST(!(createReadBadCoords[0] && createReadBadCoords[1]));
    TEST(!(setReadBadCoords[0] && setReadBadCoords[1]));

    //
    // Now test VtxArr calls in similar manner.
    // 4 passes. 2 for Nan/infinity garbage times 2 for blocked/interleave
    //
    for (int pass = 0; pass < 4; pass++) {
      double badVal =
	  (pass % 2) ?
	      std::numeric_limits<double>::quiet_NaN() :
	      std::numeric_limits<double>::infinity();
      int storage_order = (pass / 2) ? iBase_BLOCKED : iBase_INTERLEAVED;
      int other_order = (pass / 2) ? iBase_INTERLEAVED : iBase_BLOCKED;

      // Always allocate as for 3D so we can test impl. to
      // see if it wound up visiting the 'extra' data
      double *tenCoordsA = new double[10 * 3];
      for (int i = 0; i < 10 * gdim; i++)
	tenCoordsA[i] = (double) (i / gdim);
      for (int i = 10 * gdim; i < 10 * 3; i++)
	tenCoordsA[i] = badVal;

      iBase_EntityHandle *tenVerts = 0;
      int tenVertsAlloc = 0, tenVertsSize = 0;
      iMesh_createVtxArr(instance, 10, storage_order, tenCoordsA, 10 * gdim,
			 &tenVerts, &tenVertsAlloc, &tenVertsSize, &err);
      CHECK_ERR(err);

      // Read back the values created in same storage order
      double *tenCoordsB = 0;
      int tenCoordsBAlloc = 0, tenCoordsBSize = 0;
      iMesh_getVtxArrCoords(instance, tenVerts, 10, storage_order, &tenCoordsB,
			    &tenCoordsBAlloc, &tenCoordsBSize, &err);
      CHECK_ERR(err);
      bool qCoordSizeOk = tenCoordsBSize == 10 * gdim;
      TEST(qCoordSizeOk);
      bool qCoordValsOk = true;
      // All values should be equal to what we created
      for (int i = 0; i < tenCoordsBSize; i++) {
	if (tenCoordsB[i] != tenCoordsA[i])
	  qCoordValsOk = false;
      }
      // 'bad' values in 'extra' portion is also an error
      for (int i = tenCoordsBSize; i < 10 * 3 && i < tenCoordsBAlloc; i++) {
	if (tenCoordsB[i] == badVal)
	  qCoordValsOk = false;
      }
      TEST(qCoordValsOk);

      // Read back the values created in other storage order
      free(tenCoordsB);
      tenCoordsB = 0;
      tenCoordsBAlloc = 0;
      tenCoordsBSize = 0;
      iMesh_getVtxArrCoords(instance, tenVerts, 10, other_order, &tenCoordsB,
			    &tenCoordsBAlloc, &tenCoordsBSize, &err);
      CHECK_ERR(err);
      qCoordSizeOk = tenCoordsBSize == 10 * gdim;
      TEST(qCoordSizeOk);
      qCoordValsOk = true;
      bool qCoordOrderOk = false;
      for (int i = 0; i < tenCoordsBSize; i++) {
	if (tenCoordsB[i] == badVal)
	  qCoordValsOk = false;
	if (storage_order == iBase_BLOCKED) {
	  if (tenCoordsB[i] != tenCoordsA[i / gdim + 10 * (i % gdim)])
	    qCoordValsOk = false;
	}
	else { // iBase_INTERLEAVED
	  if (tenCoordsB[i / gdim] != tenCoordsA[i])
	    qCoordValsOk = false;
	}
      }
      TEST(qCoordValsOk);

      // Ok, set them all to known, zero values, even in the 'extra'
      // values, if impl is handling them or not. Confirm they are all zero
      for (int i = 0; i < 10 * 3; i++)
	tenCoordsA[i] = 0.0;
      iMesh_setVtxArrCoords(instance, tenVerts, 10, storage_order, tenCoordsA,
			    10 * gdim, &err);
      CHECK_ERR(err);
      free(tenCoordsB);
      tenCoordsB = 0;
      tenCoordsBAlloc = 0;
      tenCoordsBSize = 0;
      iMesh_getVtxArrCoords(instance, tenVerts, 10, other_order, &tenCoordsB,
			    &tenCoordsBAlloc, &tenCoordsBSize, &err);
      CHECK_ERR(err);
      qCoordValsOk = true;
      for (int i = 0; i < tenCoordsBSize; i++) {
	if (tenCoordsB[i] != 0.0)
	  qCoordValsOk = false;
      }
      TEST(qCoordValsOk);

      // Now use setVtxArr to set (possibly) bad values
      for (int i = 0; i < 10 * gdim; i++)
	tenCoordsA[i] = (double) (i / gdim);
      for (int i = 10 * gdim; i < 10 * 3; i++)
	tenCoordsA[i] = badVal;
      // Now, use setVtx to set bad values
      iMesh_setVtxArrCoords(instance, tenVerts, 10, storage_order, tenCoordsA,
			    10 * gdim, &err);
      CHECK_ERR(err);

      // Read back values just set (same storage order)
      UT_FREE(tenCoordsB);
      tenCoordsBAlloc = 0;
      tenCoordsBSize = 0;
      iMesh_getVtxArrCoords(instance, tenVerts, 10, storage_order, &tenCoordsB,
			    &tenCoordsBAlloc, &tenCoordsBSize, &err);
      CHECK_ERR(err);
      qCoordSizeOk = tenCoordsBSize == 10 * gdim;
      TEST(qCoordSizeOk);
      qCoordValsOk = true;
      // All values should be equal to what we created
      for (int i = 0; i < tenCoordsBSize; i++) {
	if (tenCoordsB[i] != tenCoordsA[i])
	  qCoordValsOk = false;
      }
      // 'bad' values in 'extra' portion is also an error
      for (int i = tenCoordsBSize; i < 10 * 3 && i < tenCoordsBAlloc; i++) {
	if (tenCoordsB[i] == badVal)
	  qCoordValsOk = false;
      }
      TEST(qCoordValsOk);

      // Now, read back again using 'other' storage order
      free(tenCoordsB);
      tenCoordsB = 0;
      tenCoordsBAlloc = 0;
      tenCoordsBSize = 0;
      iMesh_getVtxArrCoords(instance, tenVerts, 10, other_order, &tenCoordsB,
			    &tenCoordsBAlloc, &tenCoordsBSize, &err);
      CHECK_ERR(err);
      qCoordSizeOk = tenCoordsBSize == 10 * gdim;
      TEST(qCoordSizeOk);
      qCoordValsOk = true;
      qCoordOrderOk = false;
      for (int i = 0; i < tenCoordsBSize; i++) {
	if (tenCoordsB[i] == badVal)
	  qCoordValsOk = false;
	if (storage_order == iBase_BLOCKED) {
	  if (tenCoordsB[i] != tenCoordsA[i / gdim + 10 * (i % gdim)])
	    qCoordValsOk = false;
	}
	else { // iBase_INTERLEAVED
	  if (tenCoordsB[i / gdim] != tenCoordsA[i])
	    qCoordValsOk = false;
	}
      }
      TEST(qCoordValsOk);
      delete[] tenCoordsA;
      UT_FREE(tenVerts);
      UT_FREE(tenCoordsB);
    }
    iMesh_dtor(instance, &err);
    CHECK_ERR(err);
  }
  TEST(supports2D);
}
