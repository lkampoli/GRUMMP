#include <memory>

#include "GR_BaryCentricLengthScale.h"
#include "GR_Geometry.h"
#include "GR_GeomComp.h"
#include "GR_GraphicsOut.h"
#include "GR_InsertionManager.h"
#include "GR_Length.h"
#include "GR_SurfaceInsertion.h"

#include "FacetEvalTool.hpp"
#include "FacetModifyEngine.hpp"
#include "FacetQueryEngine.hpp"
#include "FacetSurface.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitPoint.hpp"

namespace GRUMMP
{

  SurfaceInserter3D::SurfaceInserter3D(VolMesh * pVMesh,
				       std::shared_ptr<Length> len,
				       double grading, double scale,
				       const char strSurfParams[]) :
      SurfaceInserter(len, grading, scale, strSurfParams), m_mesh(pVMesh), m_surfMeshes(), m_currentRefVolume(
      NULL), m_isSurfaceClosed(false)
  {

    // watches moved verts to determine what has changed
    m_mesh->addObserver(this, Observable::vertMoved);
    m_mesh->allowBdryChanges();

    GR_index_t iLen = strlen(strSurfParams), iParam = 0;
    while (iParam < iLen) {
      char cType = strSurfParams[iParam++];
      switch (cType)
	{
	// swap decider type
	// doesn't actively do anything right now, based off of
	// 2D version.
	case 'd':
	  {
	    char sdType = strSurfParams[iParam++];
	    switch (sdType)
	      {
	      case 'd':
		break;
	      }
	    break;
	  }
	  // aniso
	case 'a':
	  {
	    m_isAnisotropic = true;
	    break;
	  }
	case 'o':
	  {
	    // this is experimental as offsetting surfaces
	    // is much harder than offsetting curves
	    // only works for surfaces that are "layers"
	    //-s o0x0y0z0.01n1
	    m_useOffsets = true;
	    char temp1, temp2, temp3, temp4;
	    int surfnum;
	    double x, y, z;
	    sscanf(&strSurfParams[iParam], "%d%c%lf%c%lf%c%lf%c%d", &surfnum,
		   &temp1, &x, &temp2, &y, &temp3, &z, &temp4,
		   &m_numberOfOffsets);
	    m_offsetCoordinate = std::pair<int, CubitVector>(
		surfnum, CubitVector(x, y, z));
	    if (temp1 != 'x' || temp2 != 'y' || temp3 != 'z' || temp4 != 'n') {
	      logMessage(1, "Warning, incorrect offset input string\n");
	    }

	  }
	  break;
	case 'p':
	  {
	    // -s p10s1s2s2s3...
	    m_useOffsets = true;
	    char * temp = const_cast<char*>(&strSurfParams[iParam]);
	    m_numberOfOffsets = strtoul(temp, &temp, 10);
	    iParam += temp - &strSurfParams[iParam];
	    char * pairs = strtok(temp, "s");
	    for (GR_index_t iP = 0; iP < m_numberOfOffsets; iP++) {
	      // collect pairs
	      int s1 = strtol(pairs, &temp, 10);
	      pairs = strtok(NULL, "s");
	      int s2 = strtol(pairs, &temp, 10);
	      pairs = strtok(NULL, "s");
	      m_offsetPairs.push_back(std::pair<int, int>(s1, s2));
	    }

	    m_offsetByNormals = true;
	  }
	  break;
	case 'q': // quasi-structured
	  {
	    m_isStructured = true;
	  }
	  break;
	}

    }

    SwapDec3D = new DelaunaySwapDecider3D(true, 90.);
    SwapMan3D = new SwapManager3D(SwapDec3D, m_mesh);
    SmoothMan3D = new OptMSSmoothingManager3D(m_mesh);

    // if the mesh has known subsegs, use them
    // otherwise infer them
    std::set<Subseg*> mesh_subsegs;
    m_mesh->getAllSubsegs(mesh_subsegs);
    RefMan3D =
	(m_mesh->isSubsegMapEmpty()) ?
	    new RefinementManager3D(m_mesh, m_IPC, NULL, eBall, m_isoLength,
				    OldPriorityCalculator::SIN_DIHED) :
	    new RefinementManager3D(m_mesh, m_IPC, &mesh_subsegs, eBall,
				    m_isoLength,
				    OldPriorityCalculator::SIN_DIHED);

    RefMan3D->setQueueForQuality();

    assert(m_mesh->subsegMapIsValid());

    initializeStats(m_mesh);

  }
  SurfaceInserter3D::~SurfaceInserter3D()
  {
    finalizeStats(*m_mesh);

    GeometryQueryTool::instance()->delete_geometry();

    for (GR_index_t i = 0; i < m_surfMeshes.size(); i++)
      delete m_surfMeshes[i];

    CGMApp::instance()->shutdown();

    m_mesh->removeObserver(this);

    delete SwapMan3D;
    delete SwapDec3D;
    delete SmoothMan3D;
    delete RefMan3D;
  }

  void
  SurfaceInserter3D::insertSurfaceInMesh(const char strSurface[],
					 const char strParams[])
  {

    char strFileName[FILE_NAME_LEN];
    bool qIsListOfSurfaces = !strstr(strSurface, ".stl");
    FILE * pFSurfaceList = NULL;
    GR_index_t numVolumes = 1;
    std::set<RefVolume*> surfacesToIgnore;

    // lets start reading in the geometry
    GeometryQueryTool *gqt = GeometryQueryTool::instance();
    FacetQueryEngine *fqe = FacetQueryEngine::instance();

    // takes in either a list of surfaces in a file
    if (qIsListOfSurfaces) {
      pFSurfaceList = fopen(strSurface, "r");
      fscanf(pFSurfaceList, "%u\n", &numVolumes);
    }

    logMessage(1, "Reading from %s\n", strSurface);

    for (GR_index_t iVolume = 0; iVolume < numVolumes; iVolume++) {

      // gets the surface filename
      if (qIsListOfSurfaces)
	fscanf(pFSurfaceList, "%" FILE_NAME_LEN_STR "s\n", strFileName);
      else
	strncpy(strFileName, strSurface, FILE_NAME_LEN);
      printf("strFileName: %s\n", strFileName);
      printf("Length: %d\n", FILE_NAME_LEN);

      CubitBoolean use_feature_angle = CUBIT_TRUE;
      double feature_angle = 135;
      double tolerance = 1.e-6;
      int interp_order = 4;
      CubitBoolean smooth_non_manifold = CUBIT_FALSE;
      CubitBoolean split_surfaces = CUBIT_FALSE;
      CubitBoolean stitch = CUBIT_FALSE;
      CubitBoolean improve = CUBIT_FALSE;
      DLIList<CubitQuadFacet*> quad_facet_list;
      DLIList<CubitFacet*> tri_facet_list;
      DLIList<Surface*> surface_list;

      //This is the call to CGM to build the geometry from
      //a facetted surface representation.

      fqe->import_facets(strFileName, use_feature_angle, feature_angle,
			 tolerance, interp_order, smooth_non_manifold,
			 split_surfaces, stitch, improve, quad_facet_list,
			 tri_facet_list, surface_list, STL_FILE);

    }
    numVolumes = gqt->num_bodies();
    logMessage(MSG_MANAGER, "A total of %d volumes to insert\n", numVolumes);
    if (pFSurfaceList)
      fclose(pFSurfaceList);

    // if we have offsets, lets prepare to ignore them in the sampling process
    // otherwise, we'll sample them normally
    for (GR_index_t iP = 0; iP < m_offsetPairs.size(); iP++) {
      DLIList<RefFace*> newSurfaces;
      logMessage(1, "Offsetting surface number %d to surface number %d\n",
		 m_offsetPairs[iP].first, m_offsetPairs[iP].second);
      assert(gqt->get_ref_volume(m_offsetPairs[iP].second + 1));
      gqt->get_ref_volume(m_offsetPairs[iP].second + 1)->ref_faces(newSurfaces);
      surfacesToIgnore.insert(
	  gqt->get_ref_volume(m_offsetPairs[iP].second + 1));
    }
    // this will loop over the whole algorithm

    // volume indexing begins with 1, not 0
    for (GR_index_t iRV = 0; iRV < numVolumes; iRV++) {

      DLIList<RefFace*> allSurfaces;
      std::map<Vert*, Vert *> vertMap;

      m_currentRefVolume = gqt->get_ref_volume(iRV + 1);
      assert(m_currentRefVolume);

      if (surfacesToIgnore.count(m_currentRefVolume) == 1)
	continue;

      // a lot of this stems from the original sampling code,
      // I simply follow along with Serge's work.
      m_currentRefVolume->ref_faces(allSurfaces);
      RefFace* surface;
      FacetSurface * facetSurface;

      for (int iSurface = 0; iSurface < allSurfaces.size(); ++iSurface) {
	surface = allSurfaces.get_and_step();
	facetSurface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());
	assert(facetSurface);
	facetSurface->get_eval_tool()->compare_tol(1.e-12);
      }

      // check the volume to determine if it is closed or not
      m_isSurfaceClosed = (iFuzzyComp(m_currentRefVolume->measure(), 0.0) == 1);
      // log this
      logMessage(1, "Inserting a %s surface.\n",
		 m_isSurfaceClosed ? "closed" : "open");
      writeVTKLegacyWithoutPurge(*m_mesh, "surfins_initmesh");

      m_mesh->disallowSwapRecursion();

      m_mesh->sendEvents();

      m_mesh->purgeAllEntities();

      // calls correct sampling function based on structured or not
      !m_isStructured ?
	  m_surfMeshes.size() * sampleSurface() :
	  m_surfMeshes.size() * sampleQuasiStructured();

      assert(m_mesh->subsegMapIsValid());

      numRemovedVerts += removeVertsNearSurface();

      assert(m_mesh->subsegMapIsValid());
      m_mesh->writeTempMesh();

      numInsertedVerts += insertPointsIntoMesh();
      m_mesh->writeTempMesh();

      // make sure subseg verts are updated correctly
      // this should eventually be passed to the mesh level
      m_mesh->sendEvents();
      assert(m_mesh->subsegMapIsValid());

      m_mesh->purgeAllEntities();
      SmoothMan3D->emptyQueue();

      assert(m_mesh->subsegMapIsValid());
      assert(m_mesh->isValid());

//		if (surface creates new region) {
      setRegions();
//		}

      GR_index_t vertsBeforeRefinement = m_mesh->getNumVerts();

      assert(m_mesh->isValid());
      refineAndSmooth(strParams);

      numRefinedVerts += m_mesh->getNumVerts() - vertsBeforeRefinement;

      m_mesh->printSubsegMap();
      m_mesh->sendEvents();
      m_mesh->writeTempMesh();
      assert(m_mesh->isValid());
      //		vertMap->clear();
      //		m_mesh->purgeAllEntities(vertMap);
      //		m_mesh->updateSubsegVertDataAfterPurge(*vertMap);
      m_mesh->purgeAllEntities();
    }
  }
  int
  SurfaceInserter3D::sampleSurface()
  {
    BaryLengthScale3D BLS(m_mesh);
    // we are going to use a set of FaceQueueEntries to define the surface
    // this is the raw output from a surfmeshbuilder, and thus closed surfaces
    // accommodate it
    FaceVec restrictedFaces;
    TetMeshBuilder* TMB = nullptr;
    std::set<Subseg*> subsegSet;

    logMessage(1, "Initializing Sample Locations\n");

    // Use the SurfMeshBuilder to do the sampling
    // This is a modified version that will sample by lengthscale
    SurfMeshBuilder * SMB = new SurfMeshBuilder(m_currentRefVolume,
						SurfMeshBuilder::FACET,
						0.45 * M_PI, m_mesh,
						CurveDiscretizer::LENGTHSCALE);

    SMB->sample_all_surfaces();
    // lets determine the new locations..
    if (m_isSurfaceClosed) {
      // then we need to use a tetMeshBuilder to ensure we have all the facets
      // for an open surface this is generally true, but not fully tested
      TMB = new TetMeshBuilder(SMB);
      TMB->buildMesh(&subsegSet);
      VolMesh * sampledMesh = TMB->get_mesh();

      // new length object used for refinement
      // this creates a surface mesh the right sizing based on the
      // background mesh. Previously was re-using this code

      std::shared_ptr<GRUMMP::Length> length3D = GRUMMP::Length::Create(m_mesh,
									1., 1.);

      length3D->setLengthScale(sampledMesh);
      {
	RefinementManager3D ref3D(sampledMesh, m_IPC, &subsegSet, eBall,
				  length3D);
	ref3D.refineMesh();
      }

      sampledMesh->getAllSubsegs(subsegSet);
      // create fake restricted face set to build surface mesh
      // of the tetrahedralization
      for (GR_index_t iBdryFace = 0; iBdryFace < sampledMesh->getNumBdryFaces();
	  iBdryFace++) {
	BFace * bdryFace = sampledMesh->getBFace(iBdryFace);
	if (bdryFace->isDeleted())
	  continue;
	restrictedFaces.push_back(
	    std::make_pair(
		new FaceQueueEntry(bdryFace->getFace(0)),
		dynamic_cast<RefFace*>(bdryFace->getTopologicalParent())));
      }

    }
    else {
      // just use a regular surfaceMesh
      SMB->refine_by_background_mesh(m_mesh);
      SMB->get_subsegs_from_bridges(subsegSet);
      SMB->get_restricted_faces(restrictedFaces);

    }

    // sampling is over, now build the surfmesh
    m_surfMeshes.push_back(new SurfMesh(restrictedFaces, subsegSet));
    writeVTKLegacyWithoutPurge(*m_surfMeshes.back(), "surfins_samplesurface");

    if (m_useOffsets)
      sampleOffsetSurfaces(restrictedFaces, subsegSet);

    // do a bit of cleanup first
    for (FaceVec::iterator it = restrictedFaces.begin();
	it != restrictedFaces.end(); ++it) {
      if (it->first)
	delete it->first;
    }

    restrictedFaces.clear();
    subsegSet.clear();

    if (TMB) {
      TMB->destroyInitMesh();
      delete TMB;
    }
    delete SMB;

    sampleForIntBdryIntersections();

    return m_surfMeshes[0]->getNumVerts();

  }
  void
  SurfaceInserter3D::sampleOffsetSurfaces(FaceVec & fqes,
					  std::set<Subseg*> & subsegSet)
  {

    GeometryQueryTool * gqt = GeometryQueryTool::instance();
    CubitVector offsetCoord = m_offsetCoordinate.second;

    for (GR_index_t i = 0; i < m_numberOfOffsets; i++) {

      if (!m_offsetPairs.empty()) {
	// makes sure we have unique vertices to project
	// this is quite expensive, given that m_surfMeshes[0] contains
	// a unique set of verts already in it
	// but we want a new set of vertices
	// this is tentatively easier, TODO revisit this

	std::set<Vert*> surfVerts;
	for (FaceVec::iterator it = fqes.begin(); it != fqes.end(); ++it) {
	  surfVerts.insert((*it).first->getVert(0));
	  surfVerts.insert((*it).first->getVert(1));
	  surfVerts.insert((*it).first->getVert(2));
	}
//			DLIList<RefFace*> oldSurfaces;
	DLIList<RefFace*> newSurfaces;
	// the +1 is due to the identifiers of CGM
//			gqt->get_ref_volume(m_offsetPairs[i].first+1)->ref_faces(oldSurfaces);
	gqt->get_ref_volume(m_offsetPairs[i].second + 1)->ref_faces(
	    newSurfaces);
//			assert(oldSurfaces.size() == newSurfaces.size());
	for (int iSurf = 0; iSurf < newSurfaces.size(); iSurf++) {

	  RefFace * newSurface = newSurfaces.get_and_step();

	  DLIList<RefVertex*> newRefVerts;
	  newSurface->ref_vertices(newRefVerts);

	  for (std::set<Vert*>::iterator it = surfVerts.begin();
	      it != surfVerts.end(); ++it) {
	    Vert * vert = *it;
	    CubitVector oldPoint(vert->getCoords());
	    CubitVector newPoint;

	    // if its close to a reference point, should match it up to the corresponding one
	    // find the closest one.

	    if (vert->getVertType() == Vert::eBdryApex) {
	      assert(vert->getParentTopology()->geometry_type() == 0);
	      newPoint =
		  newRefVerts[dynamic_cast<RefVertex*>(vert->getParentTopology())->id()
		      - 1]->coordinates();
	    }
	    else {
	      newSurface->find_closest_point_trimmed(oldPoint, newPoint);
	    }
	    vert->setCoords(3, newPoint);

	  }

	}
      }
      else {
	// this corresponds to where we are offsetting a surface
	// by just offsetting the sampling points and creating a new surface mesh
	// it basically loses the curvature info but this serves as a good
	// test for whether or not the rest of the algorithm works
	std::set<Vert*> surfVerts;
	for (FaceVec::iterator it = fqes.begin(); it != fqes.end(); ++it) {
	  surfVerts.insert((*it).first->getVert(0));
	  surfVerts.insert((*it).first->getVert(1));
	  surfVerts.insert((*it).first->getVert(2));
	}
	for (std::set<Vert*>::iterator it = surfVerts.begin();
	    it != surfVerts.end(); ++it) {
	  Vert * vert = *it;
	  CubitVector coords(vert->getCoords());
	  coords = coords + offsetCoord;
	  vert->setCoords(3, coords);
	}
      }

      m_surfMeshes.push_back(new SurfMesh(fqes, subsegSet));

#ifndef NDEBUG
      char offsetFileName[FILE_NAME_LEN];
      sprintf(offsetFileName, "surfins_offset%d", i);
      writeVTK(*m_surfMeshes[i + 1], offsetFileName, "");
      printf("writing m_surfMesh %d %p\n", i + 1, m_surfMeshes[i + 1]);
#endif

    }

  }

  Face*
  SurfaceInserter3D::attemptFaceRecovery(
      std::vector<Face*>& newFaces, const GR_index_t iC,
      Vert* const faceVerts[3], std::vector<Cell*>& facesUnableToRecover,
      Cell* const faceToInsert, GR_index_t& iFacesRecovered, Cell*& cellGuess)
  {
    m_mesh->recoverFace(faceVerts[0], faceVerts[1], faceVerts[2]);
    Face* newFace = findCommonFace(faceVerts[0], faceVerts[1], faceVerts[2],
    NULL,
				   true);
    // if we can't recover the face,
    // store it, and worry about it after.
    if (!newFace) {
      facesUnableToRecover.push_back(faceToInsert);
    }
    else {
      ++iFacesRecovered;
      newFace->lock();
      cellGuess = newFace->getLeftCell();
      assert(cellGuess->doFullCheck());
    }
    newFaces[iC] = newFace;
    return newFace;
  }

  int
  SurfaceInserter3D::insertPointsIntoMesh()
  {

    m_mesh->allowEdgeSwapping();
    m_mesh->allowBdryChanges();
    m_mesh->sendEvents();
    assert(m_mesh->subsegMapIsValid());

    // again, use a SwappingInserter in most places.
    SwapMan3D->clearQueue();

    SwappingInserter3D SwapIns3D(m_mesh, SwapMan3D);
    SwapIns3D.setForcedInsertion(true);

    // initial first guess, again replaceable by a search tree.
    Cell * cellGuess = m_mesh->getCell(0);
    GR_index_t iFacesRecovered = 0;
    GR_index_t iNumVertsInserted = 0;

    // store all the new faces we create inside a vector of vectors
    // so they can be used to create internal boundaries
    std::vector<std::vector<Face*>> newFaces(m_surfMeshes.size(),
					     std::vector<Face*>());
    // reserve the memory we are going to need
    for (GR_index_t iSurf = 0; iSurf < m_surfMeshes.size(); iSurf++) {
      newFaces[iSurf].reserve(m_surfMeshes[iSurf]->getNumCells());
    }

    for (GR_index_t iSurf = 0; iSurf < m_surfMeshes.size(); iSurf++) {

      std::map<Vert*, Vert*> insertedVertMap;
      std::map<Vert*, Vert*>::iterator vertMapIt;

      std::vector<Cell*> facesUnableToRecover;

      logMessage(1, "There are %d vertices to insert in surf mesh %u.\n",
		 m_surfMeshes[iSurf]->getNumVerts(), iSurf);
      logMessage(1, "There are %d faces to recover in surf mesh %u.\n",
		 m_surfMeshes[iSurf]->getNumCells(), iSurf);

      GR_index_t interval = (m_surfMeshes[iSurf]->getNumVerts() / 10);
      // these are for the recovery process
      Vert * faceVerts[3] =
	{ pVInvalidVert, pVInvalidVert, pVInvalidVert };

      for (GR_index_t iC = 0; iC < m_surfMeshes[iSurf]->getNumCells(); iC++) {
	Cell * faceToInsert = m_surfMeshes[iSurf]->getCell(iC);

	for (GR_index_t iVert = 0; iVert < 3; iVert++) {
	  Vert * vert = faceToInsert->getVert(iVert);

	  //if we haven't inserted already, then insert it
	  if (insertedVertMap.count(vert) == 0) {
	    Vert * newVert = pVInvalidVert;
	    double coords[3];
	    bool status; // used for findCell

	    BFace* nearestBFace = findNearestBoundaryCell(vert->getCoords(),
							  coords);
	    assert(nearestBFace);
	    double distance = dDIST3D(coords, vert->getCoords());
	    double lenScale = nearestBFace->getVert(0)->getLengthScale();
	    double distRatio = distance / lenScale;
	    assert(iFuzzyComp(distRatio, 0) != -1);

	    // TODO: This tolerance was pulled out of the air with no rational basis.
	    if (distRatio > 1.e-4) {
	      // Insert the point where it is.

	      // hopefully the cellguess is close, and this is quick
	      if (!cellGuess->doFullCheck())
		cellGuess = pCInvalidCell;

	      cellGuess = m_mesh->findCell(vert->getCoords(), cellGuess,
					   status);
	      newVert = SwapIns3D.insertPoint(vert->getCoords(), cellGuess);
	      assert(newVert);

	    }
	    else {
	      // Insert the point where it projects to the bdry

	      // find the nearest boundary face, check if its on a subseg, lots of global stuff thats probably
	      // unnecessary, but stems from the fact we don't store edges.
	      // Again, this falls under the category of stuff we shouldn't need to do once the mesh object
	      // properly acknowledges subsegments.
	      // projecting to the boundary is important for tolerances
	      // the sampled point may actually be just off due to the geometry
	      // and findNearestBoundaryCell accounts for that

	      Cell * nearestCell = nearestBFace->getFace(0)->getOppositeCell(
		  nearestBFace);
	      GR_index_t iV0, iV1;
	      Subseg * subseg = NULL;

	      // if we are on a subseg, check ahead of time.
	      // this is needed for manual point insertion
	      if (dynamic_cast<TetCell*>(nearestCell)->isPointOnEdge(coords,
								     iV0, iV1))
		subseg = m_mesh->areSubsegVertsConnected(
		    nearestCell->getVert(iV0), nearestCell->getVert(iV1));

	      newVert = SwapIns3D.insertPoint(coords, nearestCell);
	      assert(newVert);

	      if (subseg) {
		double splitParam;
		subseg->compute_refinement_split_data(coords, splitParam, true);
		m_mesh->createNewSubsegsBySplit(subseg, coords, splitParam,
						newVert);
	      }
	    }

	    assert(newVert);
	    if (++iNumVertsInserted % interval == 0)
	      logMessage(1, "%d Verts Inserted\n", iNumVertsInserted);

	    // need to be VERY careful here to not screw things up for future vertex removal
	    // this is because removal by contraction looks for an acceptable vertex to
	    // contract to, and uses vertex type as an indicator
	    switch (vert->getVertType())
	      {
	      case Vert::eBdryApex:
		newVert->setType(Vert::eBdryApex);
		break;
	      case Vert::eBdryCurve:
		newVert->setType(Vert::eBdryCurve);
		break;
	      case Vert::eBdry:
		newVert->setType(Vert::eBdryTwoSide);
		break;
		// this holds for interior or other types
	      default:
		newVert->setType(Vert::eBdry);
		break;
	      }

//	    m_isoLength->setLengthScale(newVert);

	    cellGuess =
		(newVert->getFace(0)->getLeftCell()->getType() == Cell::eTet) ?
		    newVert->getFace(0)->getLeftCell() :
		    newVert->getFace(0)->getRightCell(); // reset the guess to a new valid vert
	    assert(cellGuess->doFullCheck());

	    insertedVertMap.insert(std::make_pair(vert, newVert));
	    faceVerts[iVert] = newVert;
	    m_isoLength->setLengthScale(newVert);
	  }
	  else {
	    // find the vert we've already inserted.
	    assert(insertedVertMap.count(vert) == 1);
	    vertMapIt = insertedVertMap.find(vert);
	    faceVerts[iVert] = vertMapIt->second;
	  }

	}
	assert(cellGuess->doFullCheck());

	attemptFaceRecovery(newFaces[iSurf], iC, faceVerts,
			    facesUnableToRecover, faceToInsert, iFacesRecovered,
			    cellGuess);

	if (iFacesRecovered % interval == 0) {
	  logMessage(1, "%d Faces Recovered\n", iFacesRecovered);
	}
      } // Done trying to recover this face.

      m_mesh->writeTempMesh();
      // remove leftover verts once the rest are inserted
      removeLeftoverVerts();

      if (!facesUnableToRecover.empty()) {
	std::vector<Cell*> facesStillUnableToRecover;
	GR_index_t iFacesRecoveredPreviously = 0;
	do {
	  facesStillUnableToRecover.clear();
	  iFacesRecoveredPreviously = iFacesRecovered;
	  for (GR_index_t iC = 0; iC < facesUnableToRecover.size(); iC++) {
	    Cell *cell = facesUnableToRecover[iC];
	    for (int iV = 0; iV < 3; iV++) {
	      Vert * vert = cell->getVert(iV);
	      assert(insertedVertMap.count(vert) == 1);
	      vertMapIt = insertedVertMap.find(vert);
	      faceVerts[iV] = vertMapIt->second;
	    }
	    attemptFaceRecovery(newFaces[iSurf], iC, faceVerts,
				facesStillUnableToRecover, cell,
				iFacesRecovered, cellGuess);

	  }
	  logMessage(MSG_MANAGER, "Faces recovered this pass:  %d\n",
		     iFacesRecovered - iFacesRecoveredPreviously);
	  std::swap(facesStillUnableToRecover, facesUnableToRecover);
	}
	while (iFacesRecovered != iFacesRecoveredPreviously);

	for (GR_index_t iC = 0; iC < facesUnableToRecover.size(); iC++) {
	  Cell * cell = facesUnableToRecover[iC];

	  if (newFaces[iSurf][m_surfMeshes[iSurf]->getCellIndex(cell)] != NULL)
	    continue;
	  Face * newFace = recoverFaceBySwapping(
	      insertedVertMap[cell->getVert(0)],
	      insertedVertMap[cell->getVert(1)],
	      insertedVertMap[cell->getVert(2)]);
	  if (newFace) {
	    newFace->lock();
	    newFaces[iSurf][m_surfMeshes[iSurf]->getCellIndex(cell)] = newFace;
	  }
	  else {
	    facesStillUnableToRecover.push_back(cell);
	  }
	}
	std::swap(facesUnableToRecover, facesStillUnableToRecover);
	facesStillUnableToRecover.clear();
	// this is a recovery loop that works in pairs, and is fairly
	// specific to surface insertion, where I can swap an edge
	// to make a pair of faces recoverable.

	// Thus this relies on two faces to be recovered.
	// This allows for the swapping of the edge we can't recover

	// lets look for pairs of faces to try and just swap the edge
	if (!facesUnableToRecover.empty())
	  logMessage(1, "Unable to recover %lu faces in first pass\n",
		     facesUnableToRecover.size());

	for (GR_index_t iCA = 0; iCA < facesUnableToRecover.size(); iCA++) {
	  Cell * cellA = facesUnableToRecover[iCA];
	  Vert * vertsA[3];

	  // get the verts of cellA, store them for comparison.
	  for (GR_index_t iV = 0; iV < 3; iV++)
	    vertsA[iV] = insertedVertMap[cellA->getVert(iV)];
	  Face* checkFace = findCommonFace(vertsA[0], vertsA[1], vertsA[2],
	  NULL,
					   true);
	  if (checkFace) {
	    // Recovered this one already.  Skip to the next.
	    continue;
	  }

	  // look through the list afterwards to find if theres a face we match up with
	  for (GR_index_t iCB = iCA + 1; iCB < facesUnableToRecover.size();
	      iCB++) {
	    Cell * cellB = facesUnableToRecover[iCB];
	    Vert * vertsB[3];
	    for (GR_index_t iV = 0; iV < 3; iV++)
	      vertsB[iV] = insertedVertMap[cellB->getVert(iV)];
	    Face * sharedFace = NULL;

	    // if the two cells share a face, then we can attempt to swap.
	    for (GR_index_t iF = 0; iF < 3; iF++)
	      if (cellA->hasFace(cellB->getFace(iF)))
		sharedFace = cellB->getFace(iF);
	    // if they share a face, try and recover them both
	    if (sharedFace) {
	      std::pair<Face*, Face*> newFacePair(pFInvalidFace,
	      pFInvalidFace);
	      Vert * sharedVerts[] =
		{ insertedVertMap[sharedFace->getVert(0)],
		    insertedVertMap[sharedFace->getVert(1)] };
	      bool qRecovered = recoverFacePairByEdgeSwap(sharedFace,
							  sharedVerts, vertsA,
							  vertsB, newFacePair);
	      if (qRecovered) {

		newFaces[iSurf][m_surfMeshes[iSurf]->getCellIndex(cellA)] =
		    newFacePair.first;
		newFaces[iSurf][m_surfMeshes[iSurf]->getCellIndex(cellB)] =
		    newFacePair.second;
	      }
	      else {
		facesStillUnableToRecover.push_back(cellA);
	      }
	    }

	  }
	}
#ifndef NDEBUG
	writeVTK_CellSet(facesStillUnableToRecover,
			 "surfins_facesUnableToRecover.vtk");
	if (m_leftoverVertsToRemove.size() > 0)
	  logMessage(1, "Unable to remove %lu verts\n",
		     m_leftoverVertsToRemove.size());
	if (facesUnableToRecover.size() > 0)
	  logMessage(1, "Unable to recover %lu faces\n",
		     facesStillUnableToRecover.size());

	if (!facesStillUnableToRecover.empty() && iSurf == 0) {
	  writeVTKLegacyWithoutPurge(*m_mesh, "surfins_unable.vtk");
	  vFatalError("Unable to recover all faces", "surface insertion");
	}
#endif
      }
      assert(m_mesh->isValid());

      logMessage(1, "Inserted %d faces.\n", m_surfMeshes[iSurf]->getNumCells());

      // this is where all the boundary faces are created
      // they are created after locking and insertion as even internal
      // boundaries can cause problems in the recovery stage.
      for (GR_index_t iC = 0; iC < m_surfMeshes[iSurf]->getNumCells(); iC++) {

	assert(newFaces[iSurf][iC]);

	BFace * bdryFace = m_mesh->createBFace(newFaces[iSurf][iC]);
	bdryFace->setTopologicalParent(m_surfMeshes[iSurf]->getSurface(iC));

	assert(bdryFace->doFullCheck());
	assert(bdryFace->getType() == Cell::eIntTriBFace);
	/*
	 * This was intended to create prisms, but that requires
	 * createPrismCell to be modified to work on existing meshes
	 * and for offsets to handled
	 * TODO: If this feature is required

	 // lets do a test
	 assert(newFaces[iSurf-1][iC]);
	 Vert * prism[6];
	 prism[0] = newFaces[iSurf-1][iC]->getVert(0);
	 prism[1] = newFaces[iSurf-1][iC]->getVert(1);
	 prism[2] = newFaces[iSurf-1][iC]->getVert(2);
	 prism[3] = newFaces[iSurf][iC]->getVert(0);
	 prism[4] = newFaces[iSurf][iC]->getVert(1);
	 prism[5] = newFaces[iSurf][iC]->getVert(2);
	 printf("creating prism\n");
	 for(int iPris = 0; iPris < 6; iPris++){
	 prism[iPris]->printVertInfo();
	 }
	 bool qExist;
	 m_mesh->createPrismCell(qExist,prism[0],prism[1],prism[2],prism[3],prism[4],prism[5],10);
	 */

      }

      assert(m_mesh->subsegMapIsValid());

      // now set up the new subsegments
      std::set<Subseg*> surfMeshSubsegs;
      m_surfMeshes[iSurf]->getSubsegs(surfMeshSubsegs);

      for (std::set<Subseg*>::iterator it_sub = surfMeshSubsegs.begin();
	  it_sub != surfMeshSubsegs.end(); ++it_sub) {

	Subseg * subseg = *it_sub;
	Subseg * newSubseg;

	Vert* v0 = subseg->get_beg_vert();
	Vert* v1 = subseg->get_end_vert();

	vertMapIt = insertedVertMap.find(v0);
	v0 = (*vertMapIt).second;
	v0->setType(vertMapIt->first->getVertType());
	vertMapIt = insertedVertMap.find(v1);
	v1 = (*vertMapIt).second;
	v1->setType(vertMapIt->first->getVertType());

	if (m_mesh->areSubsegVertsConnected(v0, v1) == NULL) {

	  if (subseg->get_ref_edge()) {
	    newSubseg = new Subseg(v0, v1, subseg->get_beg_param(),
				   subseg->get_end_param(),
				   subseg->get_ref_edge());
	  }
	  else if (subseg->get_ref_face()) {
	    newSubseg = new Subseg(v0, v1, subseg->get_ref_face());
	  }
	  else {
	    newSubseg = new Subseg(v0, v1);
	  }
	  m_mesh->addSubseg(newSubseg);
	}
      }
      //		m_mesh->printSubsegMap();
      assert(m_mesh->subsegMapIsValid());
      removeLeftoverVerts();
      m_mesh->sendEvents();
    } // looped over all surfaces

    logMessage(1, "Unable to remove %lu verts\n",
	       m_leftoverVertsToRemove.size());

    m_mesh->evaluateQuality();
    removeVertsAfterInsertion();
    return iNumVertsInserted;
  }

  Face *
  SurfaceInserter3D::recoverFaceBySwapping(Vert * pV0, Vert * pV1, Vert * pV2)
  {

    Face * face = findCommonFace(pV0, pV1, pV2, NULL, true);

    // likely because it can't recover edges. So Swap
    // for quality around edges
    logMessage(1, "Swapping nearby faces ");
    // If for whatever reason, we can't get rid of an edge in our way,
    // try and swap every face connected to one of the vertices we are after

    GRUMMP::MinMaxDihedSwapDecider3D tempSwapDec3D(true);
    GRUMMP::SwapManager3D tempSwapMan3D(&tempSwapDec3D, m_mesh);
    GR_index_t iVert = 0;
    Vert * verts[] =
      { pV0, pV1, pV2 };
    GR_index_t iSwaps = 0;
    do {
      for (int iF = 0; iF < verts[iVert]->getNumFaces(); iF++) {
	Face * faceToSwap = verts[iVert]->getFace(iF);
	if (faceToSwap->isNotDeleted() && tempSwapMan3D.swapFace(faceToSwap)) {
	  m_mesh->recoverFace(pV0, pV1, pV2);
	  face = findCommonFace(pV0, pV1, pV2, NULL, true);
	  iSwaps++;
	}
	if (face)
	  break;
      }
    }
    while (!face && iVert++ < 2);
    logMessage(1, " ... Swapped %d times in desperation\n", iSwaps);
    return face;
  }
  bool
  SurfaceInserter3D::recoverFacePairByEdgeSwap(
      Face * sharedFace, Vert * sharedVerts[], Vert * vertsA[], Vert * vertsB[],
      std::pair<Face*, Face*> & facePair)
  {

    FaceSwapInfo2D FSI2D;
    FSI2D.calculateFaceInfo(sharedFace);
    if (FSI2D.isSwappable()) {
      // find new verts to form the cells
      // Swapping Manager isn't set up for surface meshes
      // so this takes its place.
      // the swap is not actually done
      // rather two difference faces are attempted to be recovered
      // if we cannot recover the edge, don't even bother
      // to avoid the awkward situation of one of the two faces recovered
      Vert * pV0 = sharedVerts[0];
      Vert * pV1 = sharedVerts[1];

      Vert * pV2 = NULL;
      Vert * pV3 = NULL;

      for (GR_index_t iV = 0; iV < 3; iV++)
	if (vertsA[iV] != pV0 && vertsA[iV] != pV1)
	  pV2 = vertsA[iV];

      for (GR_index_t iV = 0; iV < 3; iV++)
	if (vertsB[iV] != pV0 && vertsB[iV] != pV1)
	  pV3 = vertsB[iV];

      assert(pV2 && pV3);
      // if we can't even recover the potentially new edge,
      // don't even try this funny game
      if (m_mesh->recoverEdge(pV2, pV3, true) != 1)
	return false;

      // try and recover the new faces now
      m_mesh->recoverFace(pV0, pV2, pV3);
      Face * newFaceA = findCommonFace(pV0, pV2, pV3, NULL, true);

      m_mesh->recoverFace(pV1, pV2, pV3);
      Face * newFaceB = findCommonFace(pV1, pV2, pV3, NULL, true);
      if ((newFaceA && !newFaceB) || (!newFaceA && newFaceB)) {
	logMessage(
	    1,
	    "Warning: tried to recover a pair of faces by swapping edge only one face recovered!!\n");
	return false;

      }
      newFaceA->lock();
      newFaceB->lock();
      facePair.first = newFaceA;
      facePair.second = newFaceB;
      return true;
    }
    return false;
  }

  void
  SurfaceInserter3D::refineAndSmooth(const char strParams[])
  {

    GR_index_t iLen = strlen(strParams), iParam = 0;

    /// set everything up
    m_mesh->sendEvents();
    m_mesh->allowEdgeSwapping();
    m_mesh->allowSwapRecursion();
    m_mesh->evaluateQuality();
    assert(m_mesh->subsegMapIsValid());
    // this should be replaced with a way to just close queues
    SwapMan3D->clearQueue();

    //Identical to tetra input
    while (iParam < iLen) {
      char cType = strParams[iParam++];
      assert(cType);
      switch (cType)
	{
	case 'W':
	case 'w':  // Do some swapping
	  {
	    char cMeasure = strParams[iParam++];
	    GRUMMP::SwapDecider3D *pSD3D = NULL;
	    switch (cMeasure)
	      {
	      case 'I':
	      case 'i': // Insphere criterion
		logMessage(1, "Doing Delaunay face swapping...\n");
		m_mesh->setSwapType(eDelaunay);
		pSD3D = new GRUMMP::DelaunaySwapDecider3D(true);
		break;
	      case 'M':
	      case 'm': // Minmax dihedral criterion
		logMessage(
		    1,
		    "Doing face swapping to reduce maximum dihedral angle...\n");
		m_mesh->setSwapType(eMaxDihed);
		pSD3D = new GRUMMP::MinMaxDihedSwapDecider3D(true);
		break;
	      case 'S':
	      case 's': // Maxmin sine criterion
		logMessage(
		    1,
		    "Doing face swapping to increase minimum sine of dihedral angle...\n");
		m_mesh->setSwapType(eMinSine);
		pSD3D = new GRUMMP::MaxMinSineSwapDecider3D(true);
		break;
	      case 'U':
	      case 'u': // Uniform vertex degree
		logMessage(
		    1, "Doing face swapping to equalize vertex degree...\n");
		pSD3D = new GRUMMP::UniformDegreeSwapDecider3D(true);
		break;
	      default:
		vFatalError("Bad option for swapping criterion", "main()");
		break;
	      }
	    SwapMan3D->changeSwapDecider(pSD3D);
	    SwapMan3D->swapAllQueuedFaces(false);
	    SwapMan3D->clearQueue();

	    delete pSD3D;

	  }
	  break;
	case 'M':
	case 'm': // Do some smoothing
	  {
	    logMessage(1, "Doing Local Vertex Smoothing to improve quality.\n");
	    m_mesh->disallowBdryChanges(); // don't smooth boundaries
	    char cTech = strParams[iParam++];
	    int iFunc = atoi(&strParams[iParam++]);
	    SmoothMan3D->setSmoothingTechnique(cTech);
	    if (iFunc == 0)
	      iFunc = 10;
	    SmoothMan3D->setSmoothingGoal(iFunc);
	    SmoothMan3D->smoothAllQueuedVerts();
	    m_mesh->sendEvents();
	    m_mesh->allowBdryChanges(); // reset this
	  }
	  break;
	case 'r':
	case 'R':
	  {

	    writeVTKLegacyWithoutPurge(*m_mesh, "surfins_prerefine.vtk");

	    //			logMessage(1,"Broken %d cells\n",m_mesh->breakBadBdryCells(true));
	    //			m_mesh->disallowBdryChanges();
	    //			logMessage(1,"Repaired %d cells\n",m_mesh->repairBadCells(true));
	    //			m_mesh->allowBdryChanges();
	    RefMan3D->refineAllQueuedCellsForQuality(10.);

	    writeVTKLegacyWithoutPurge(*m_mesh, "surfins_postrefine.vtk");

	    m_mesh->evaluateQuality();
	  }
	  break;
	}
    }

    assert(m_mesh->subsegMapIsValid());
    m_mesh->disallowBdryChanges();

    assert(m_mesh->isValid());
    // this ensures cleanup is done properly, SwapManager likes to have a SwapDecider when cleaning up.
    SwapMan3D->changeSwapDecider(SwapDec3D);

  }
  int
  SurfaceInserter3D::removeVertsAfterInsertion()
  {

    std::set<Vert*> cantRemove;
    for (std::set<Vert*>::iterator it = m_vertsToRemoveAfterInsertion.begin();
	it != m_vertsToRemoveAfterInsertion.end();) {
      Vert * vertToRemove = *it;
      if (!vertToRemove->isValid())
	continue;
      if (removeVert(vertToRemove))
	m_vertsToRemoveAfterInsertion.erase(it++);
      else {
	cantRemove.insert(vertToRemove);
	++it;
      }
    }

    assert(m_mesh->subsegMapIsValid());

    logMessage(1, "Removed %lu Verts of %lu after surface recovery.\n",
	       m_vertsToRemoveAfterInsertion.size(),
	       cantRemove.size() + m_vertsToRemoveAfterInsertion.size());

    m_leftoverVertsToRemove.insert(cantRemove.begin(), cantRemove.end());
    return m_vertsToRemoveAfterInsertion.size() - cantRemove.size();
  }
  int
  SurfaceInserter3D::removeLeftoverVerts()
  {
    int numRemoved = 0;
    for (std::set<Vert*>::iterator it = m_leftoverVertsToRemove.begin();
	it != m_leftoverVertsToRemove.end();) {
      Vert * vertToRemove = *it;
      if (!vertToRemove->isValid())
	continue;
      if (removeVert(vertToRemove)) {
	numRemoved++;
	m_leftoverVertsToRemove.erase(it++);
      }
      else
	++it;
    }
    return numRemoved;
  }

  int
  SurfaceInserter3D::removeVertsNearSurface()
  {

#ifndef NDEBUG
    std::vector<CubitVector> pointList;
#endif

    BaryLengthScale3D BLS(m_mesh);

    std::set<Vert*> vertsToRemove;

    GR_index_t numTooClose = 0, numEncroach = 0;
    Cell * cellGuess = m_mesh->getCell(0);

    bool status;
    double adCircCent[3];
    double radius;
    for (GR_index_t iSurf = 0; iSurf < m_surfMeshes.size(); iSurf++) {
      std::set<Subseg*> subsegSet;
      m_surfMeshes[iSurf]->getSubsegs(subsegSet);
      // for each face, check if there is a vertex too close

      for (GR_index_t iC = 0; iC < m_surfMeshes[iSurf]->getNumCells(); iC++) {
	Cell * surfCell = m_surfMeshes[iSurf]->getCell(iC);
	calcCircumcenter3D(surfCell->getVert(0), surfCell->getVert(1),
			   surfCell->getVert(2), adCircCent);
	radius = dDIST3D(adCircCent, surfCell->getVert(0)->getCoords());
	logMessage(3, "Checking surf tri %3d.\n", iC);
	logMessage(3, "  Verts: (%8f %8f %8f), (%8f %8f %8f), (%8f %8f %8f)\n",
		   surfCell->getVert(0)->x(), surfCell->getVert(0)->y(),
		   surfCell->getVert(0)->z(), surfCell->getVert(1)->x(),
		   surfCell->getVert(1)->y(), surfCell->getVert(1)->z(),
		   surfCell->getVert(2)->x(), surfCell->getVert(2)->y(),
		   surfCell->getVert(2)->z());
	logMessage(3, "  Circumcent: (%8f %8f %8f)  Radius %8f\n",
		   adCircCent[0], adCircCent[1], adCircCent[2], radius);

	std::set<Vert*> neighVerts;
	// build the neighbourhood of verts connected to verts in a cell
	for (GR_index_t iV = 0; iV < 3; iV++) {
	  cellGuess = m_mesh->findCell(surfCell->getVert(iV)->getCoords(),
				       cellGuess, status);
	  for (GR_index_t iCV = 0; iCV < 4; iCV++) {
	    for (int iF = 0; iF < cellGuess->getVert(iCV)->getNumFaces();
		iF++) {
	      neighVerts.insert(
		  cellGuess->getVert(iCV)->getFace(iF)->getVert(0));
	      neighVerts.insert(
		  cellGuess->getVert(iCV)->getFace(iF)->getVert(1));
	      neighVerts.insert(
		  cellGuess->getVert(iCV)->getFace(iF)->getVert(2));
	    }
	  }
	}
	// for each neighboring vertex, if it hasn't already been checked
	// check if it encroaches a subseg, then check it
	// if deletion is requested, we know we have already seen that vertex
	for (std::set<Vert*>::iterator itv = neighVerts.begin();
	    itv != neighVerts.end(); ++itv) {
	  Vert * vertToCheck = *itv;
	  if (vertToCheck->isDeletionRequested()
	      || vertToCheck->getVertType() == Vert::eBdryApex)
	    continue;
	  logMessage(3, "  Checking vert @ (%8f %8f %8f)", vertToCheck->x(),
		     vertToCheck->y(), vertToCheck->z());

	  for (GR_index_t iFace = 0; iFace < 3; iFace++) {
	    if (surfCell->getFace(iFace)->getOppositeCell(surfCell)->isBdryCell()) {

	      Subseg * subseg = m_surfMeshes[iSurf]->areSubsegVertsConnected(
		  surfCell->getFace(iFace)->getVert(0),
		  surfCell->getFace(iFace)->getVert(1));

	      if (subseg
		  && SubsegManager::vert_inside_subseg_ball(vertToCheck,
							    subseg)) {
		//then we are encroached
		// markToDelete as an easy way to see if we've already checked this vertex
		vertToCheck->markToDelete();
		logMessage(3, " marked");
		break;
	      }
	    }
	  }
	  // if it encroaches a subseg, or a face
	  if (vertToCheck->isDeletionRequested()
	      || iFuzzyComp(dDIST3D(adCircCent, vertToCheck->getCoords()),
			    radius) < 1) {
	    vertsToRemove.insert(vertToCheck);
	    vertToCheck->markToDelete();
	    logMessage(3, " marked");
	    numEncroach++;

#ifndef NDEBUG
	    pointList.push_back(CubitVector(vertToCheck->getCoords()));
#endif

	  }
	  // Finally, check if its too close for removal later
	  else if ((BLS.getLengthScale(surfCell->getVert(0)) * 0.5
	      > dDIST3D(vertToCheck->getCoords(),
			surfCell->getVert(0)->getCoords()))
	      || (BLS.getLengthScale(surfCell->getVert(1)) * 0.5
		  > dDIST3D(vertToCheck->getCoords(),
			    surfCell->getVert(1)->getCoords()))
	      || (BLS.getLengthScale(surfCell->getVert(2)) * 0.5
		  > dDIST3D(vertToCheck->getCoords(),
			    surfCell->getVert(2)->getCoords()))) {
	    assert(!vertToCheck->isDeletionRequested());
	    std::pair<std::set<Vert*>::iterator, bool> result =
		m_vertsToRemoveAfterInsertion.insert(vertToCheck);
	    if (result.second) {
	      numTooClose++;
	      logMessage(3, " close");
#ifndef NDEBUG
//						pointList.push_back(
//								CubitVector(vertToCheck->getCoords()));
#endif
	    }
	  }
	  logMessage(3, "\n");
	} //neighVerts
      } // surfMesh Cells
    } //num surfs
    assert(numTooClose == m_vertsToRemoveAfterInsertion.size());

    std::set<Vert*>::iterator iter, iterEnd = vertsToRemove.end();
    for (iter = vertsToRemove.begin(); iter != iterEnd; ++iter) {
      Vert* tmp = *iter;
      size_t count = m_vertsToRemoveAfterInsertion.erase(tmp);
      assert(count == 1 || count == 0);
      numTooClose -= count;
    }

    assert(numTooClose == m_vertsToRemoveAfterInsertion.size());

    logMessage(3,
	       "Removing %u verts too close, %u encroach faces or subsegs.\n",
	       numTooClose, numEncroach);

#ifndef NDEBUG
    writeVTKPointSet(pointList, NULL, "surfins_removepts", "");
    pointList.clear();
#endif
    GR_index_t numVertsBefore = m_mesh->getNumVertsInUse();

    std::set<Vert*> vertsToTryAgain;
    GR_index_t numRemoved = 0;
    do {
      vertsToTryAgain.clear();
      numRemoved = 0;
      for (std::set<Vert*>::iterator it = vertsToRemove.begin();
	  it != vertsToRemove.end(); ++it) {
	bool removed = removeVert(*it);
	logMessage(2, "Vert @ (%10f %10f %10f) %s\n", (*it)->x(), (*it)->y(),
		   (*it)->z(), removed ? "removed" : "not removed");
	if (removed) {
	  numRemoved++;
	}
	else {
	  vertsToTryAgain.insert(*it);
	}
      }
      std::swap(vertsToTryAgain, vertsToRemove);
      SwapMan3D->swapAllQueuedFacesRecursive(); // Swap towards Delaunay
    }
    while (numRemoved != 0);
    m_leftoverVertsToRemove.insert(vertsToRemove.begin(), vertsToRemove.end());

    //	m_mesh->removeTaggedVerts(false,&m_leftoverVertsToRemove);

    GR_index_t numVertsAfter = m_mesh->getNumVertsInUse();
    // This next clause is already taken care of by elimination encroaching
    // verts from the too-close list.
//	for(std::set<Vert*>::iterator it = m_vertsToRemoveAfterInsertion.begin();
//			it != m_vertsToRemoveAfterInsertion.end(); ){
//		Vert * vert = *it;
//		std::set<Vert*>::iterator itTmp = it;
//		it++;
//		if(vert->isDeleted() || m_leftoverVertsToRemove.count(vert)) {
//			m_vertsToRemoveAfterInsertion.erase(itTmp);
//		}
//	}

    assert(m_mesh->subsegMapIsValid());
    logMessage(1, "\nMesh has %u verts before removal, %u verts after.\n",
	       numVertsBefore, numVertsAfter);
    logMessage(1, "Removing %lu Verts of %u\n",
	       numEncroach - m_leftoverVertsToRemove.size(), numEncroach);

    return numEncroach;

  }
// basically need to be way more careful to clean up subsegmap
// this should eventually be handled by the volMesh class.
  bool
  SurfaceInserter3D::removeVert(Vert * vert)
  {
    int iSwaps;
    if (vert->getVertType() != Vert::eBdryCurve) {
      return m_mesh->removeVert(vert, iSwaps);
    }
    else {
      std::pair<std::set<Subseg*>::const_iterator,
	  std::set<Subseg*>::const_iterator> subseg_its =
	  m_mesh->getSubsegRange(vert);

      bool qRemoved = m_mesh->removeVert(vert, iSwaps);

      if (qRemoved) {
	std::set<Subseg*>::const_iterator its = subseg_its.first;
	assert(its != subseg_its.second);
	Subseg * s0 = *its;
	Subseg * s1 = *(++its);
	assert(its != subseg_its.second);
	assert(++its == subseg_its.second);

	Vert * v0 =
	    (s0->get_beg_vert() == vert) ?
		s0->get_end_vert() : s0->get_beg_vert();
	Vert * v1 =
	    (s1->get_beg_vert() == vert) ?
		s1->get_end_vert() : s1->get_beg_vert();
	s0->mark_deleted();
	s1->mark_deleted();
	m_mesh->removeSubseg(s0);
	m_mesh->removeSubseg(s1);
	m_mesh->recoverEdge(v0, v1);
	m_mesh->addSubseg(new Subseg(v0, v1));
      }
      return qRemoved;
    }
  }

  void
  SurfaceInserter3D::sampleForIntBdryIntersections()
  {

    // theoretically, subsegs are in right place here.. thats fine, now to move bdry faces
    if (m_mesh->getNumIntBdryFaces() == 0 || m_useOffsets)
      return;

    std::set<std::pair<Vert *, Vert *>> subsegEdges;

    // get verts from subsegs
    // get initial set of subsegment verts,
    // make sure they are identified as bdry curve ones for this
    // In the subseg manager, when sampling, it will have identified
    // intersections with internal boundaries and labelled them as eBdryTwoSide

    // mark moved verts as "active"
    // this reuses that flag
    for (GR_index_t iV = 0; iV < m_surfMeshes[0]->getNumVerts(); iV++) {
      if (m_surfMeshes[0]->getVert(iV)->getVertType() == Vert::eBdryTwoSide) {
	m_surfMeshes[0]->getVert(iV)->markActive();
	// change the type to what they should be
	m_surfMeshes[0]->getVert(iV)->setType(Vert::eBdryApex);
      }
    }
    for (GR_index_t iC = 0; iC < m_surfMeshes[0]->getNumCells(); iC++) {
      Cell * cell = m_surfMeshes[0]->getCell(iC);
      RefFace * refFace = m_surfMeshes[0]->getSurface(iC);

      GR_index_t iMovedVerts = cell->getVert(0)->isActive()
	  + cell->getVert(1)->isActive() + cell->getVert(2)->isActive();
      // if two verts have been moved already, then don't check this cell.

      if (iMovedVerts < 2) {

	bool qEdge[3];
	qEdge[0] = qEdge[1] = qEdge[2] = false;
	CubitVector intersections[3];
	BFace * bfaces[3];
	for (GR_index_t iB = 0; iB < m_mesh->getNumIntBdryFaces(); iB++) {
	  bfaces[0] = bfaces[1] = bfaces[2] = m_mesh->getIntBFace(iB);
	  for (int iEdge = 0; iEdge < 3; iEdge++) {
	    Vert * vA = cell->getVert(iEdge % 3);
	    Vert * vB = cell->getVert((iEdge + 1) % 3);
	    if (vA->isActive() && vB->isActive())
	      continue;
	    qEdge[iEdge] = (qEdge[iEdge]
		|| computeRefFaceIntersectionsWithBFace(refFace,
							vA->getCoords(),
							vB->getCoords(),
							bfaces[iEdge],
							intersections[iEdge]));
	  }
	  if (qEdge[0] + qEdge[1] + qEdge[2] == 3)
	    break;
	}
	//			assert(qEdge[0] + qEdge[1] + qEdge[2] != 1);
	// find intersections between surface mesh and internal boundary faces.
	// basic idea -> move vertices
	// allow at most two vertices of a triangle to move, otherwise risk collapsing a triangle
	// keep track of moved vertices to prevent this.
	// given triangle can intersect at most in two places. There is one vertex that is connected to both
	// intersected edges. Call that the closest. Use this to guide the rest
	// project all intersections normal to surface, not to inersection point.
	for (int iEdge = 0; iEdge < 3; iEdge++) {
	  if (qEdge[iEdge]) {
	    Vert * vA = cell->getVert(iEdge % 3);
	    Vert * vB = cell->getVert((iEdge + 1) % 3);

	    if (vA->isActive() || vB->isActive())
	      continue;
	    double distA =
		(CubitVector(vA->getCoords()) - intersections[iEdge]).length();
	    double distB =
		(CubitVector(vB->getCoords()) - intersections[iEdge]).length();
	    if (distA < distB) {
	      vA->setCoords(3, intersections[iEdge]);
	      vA->markActive();
	      if (vA->getVertType() == Vert::eBdryCurve)
		vA->setType(Vert::eBdryApex);
	      if (vA->getVertType() == Vert::eInterior)
		vA->setType(Vert::eBdryCurve);
	    }
	    else if (distA > distB) {
	      vB->setCoords(3, intersections[iEdge]);
	      vB->markActive();
	      if (vB->getVertType() == Vert::eBdryCurve)
		vB->setType(Vert::eBdryApex);
	      if (vB->getVertType() == Vert::eInterior)
		vB->setType(Vert::eBdryCurve);
	    }
	  }
	}
      }
      for (GR_index_t iFace = 0; iFace < 3; iFace++) {
	Face * face = cell->getFace(iFace);
	Vert * vA = face->getVert(0);
	Vert * vB = face->getVert(1);
	if (vA->isActive() && vB->isActive()) {
	  // check the opposite cell to see if its been collapsed by moving
	  // we really don't care about a collapsed triangle
	  // all we care about is that we get the subsegments right, as we'll
	  // re-sample based on this discretization.

	  Vert * oppVert = face->getOppositeCell(cell)->getOppositeVert(face);

	  // If a triangle is collapsed, ignore the longest edge of it
	  // and the subsegments are fine.
	  if (oppVert->isActive()
	      && iFuzzyComp(face->calcSize(), dDist(vA, oppVert)) >= 0
	      && iFuzzyComp(face->calcSize(), dDist(vB, oppVert)) >= 0)
	    continue;

	  if (subsegEdges.insert(std::pair<Vert *, Vert *>(vA, vB)).second
	      && subsegEdges.insert(std::pair<Vert *, Vert *>(vB, vA)).second)
	    m_surfMeshes[0]->addSubseg(new Subseg(vA, vB, refFace));

	}
      }
    }

#ifndef NDEBUG
    m_surfMeshes[0]->writeSubsegMap();
    writeVTK(*m_surfMeshes[0], "surfins_movedsamplesurface", "");
#endif

    if (subsegEdges.empty())
      return;

    FaceVec restrictedFaces;
    std::set<Subseg*> subsegSet;

    logMessage(1, "Found Internal Boundary Intersections\n");

    // now we need to re-sample the mesh, with the intersections forming a new subsegment
    assert(m_surfMeshes[0]->isValid());
    m_surfMeshes[0]->purgeAllEntities();
    SurfMeshBuilder SMB(m_surfMeshes[0], m_currentRefVolume,
			SurfMeshBuilder::FACET);

    delete m_surfMeshes[0];
    m_surfMeshes.clear();

    if (m_isSurfaceClosed) {
      TetMeshBuilder TMB(&SMB);
      TMB.buildMesh(&subsegSet);
      VolMesh * volMesh = TMB.get_mesh();
      volMesh->purgeAllEntities();
      for (GR_index_t i = 0; i < volMesh->getNumBdryFaces(); i++) {
	BFace * bdryFace = volMesh->getBFace(i);
	restrictedFaces.push_back(
	    std::make_pair(
		new FaceQueueEntry(bdryFace->getFace(0)),
		dynamic_cast<RefFace*>(bdryFace->getTopologicalParent())));
      }
      m_surfMeshes.push_back(new SurfMesh(restrictedFaces, subsegSet));

    }
    else {

      SMB.refine_by_background_mesh(m_mesh);
      SMB.get_subsegs_from_bridges(subsegSet);
      SMB.get_restricted_faces(restrictedFaces);
      m_surfMeshes.push_back(new SurfMesh(restrictedFaces, subsegSet));

    }

    for (FaceVec::iterator it = restrictedFaces.begin();
	it != restrictedFaces.end(); ++it) {
      if (it->first)
	delete it->first;
    }

#ifndef NDEBUG
    writeVTK(*m_surfMeshes[0], "surfins_newsamplesurface", "");
    assert(m_mesh->isValid());
#endif

  }

  BFace*
  SurfaceInserter3D::findNearestBoundaryCell(const double * coords,
					     double * projectedCoords)
  {

    //	Need to find the boundaryFace that is closest;
    //  Assumes there is one and the point is close.
    // TODO: Tree this
    double dist = LARGE_DBL, new_dist = 0.;
    BFace * nearestBFace = pBFInvalidBFace;
    // find the closest one by projection
    for (GR_index_t iF = 0; iF < m_mesh->getNumTotalBdryFaces(); iF++) {
      BFace * bdryFace = m_mesh->getBFace(iF);
      if (bdryFace->isDeleted())
	continue;

      bdryFace->findClosestPointOnBFace(coords, projectedCoords);
      new_dist = dDIST3D(projectedCoords, coords);
      if (new_dist < dist) {
	dist = new_dist;
	nearestBFace = bdryFace;
      }
    }
    for (GR_index_t i = 0; i < 3; i++)
      projectedCoords[i] = coords[i];

    nearestBFace->getFace(0)->projectOntoFace(projectedCoords);
    assert(nearestBFace->doFullCheck());
    return nearestBFace;

  }
  void
  SurfaceInserter3D::setRegions()
  {
    // find the max region in the mesh
    // For each region, identify exactly one cell that is in that region.
    Cell* pCRegionReps[iMaxRegion];
    for (int ii = 0; ii < iMaxRegion; ++ii) {
      pCRegionReps[ii] = nullptr;
    }

    int iRegion = 0;
    for (GR_index_t iC = 0; iC < m_mesh->getNumCells(); iC++) {
      Cell *pC = m_mesh->getCell(iC);
      int iThisReg = pC->getRegion();
      iRegion = max(iRegion, iThisReg);
      if (pCRegionReps[iThisReg] == nullptr) {
	pCRegionReps[iThisReg] = pC;
      }
      pC->setRegion(iInvalidRegion);
    }
    // Floodfill all existing regions.
    for (int iReg = 1; iReg <= iRegion; iReg++) {
      logMessage(1, "Labeling existing region %d\n", iReg);
      setRegionFromCell(pCRegionReps[iReg], iReg);
    }

    iRegion++;
    // check all internal boundaries, recoloring.
    for (GR_index_t iB = 0; iB < m_mesh->getNumIntBdryFaces(); iB++) {
      BFace * bFace = m_mesh->getIntBFace(iB);
      bFace->getBdryCondition();
      Cell *cell0 = bFace->getFace(0)->getOppositeCell(bFace);
      Cell *cell1 = bFace->getFace(1)->getOppositeCell(bFace);
      if (cell0->getRegion() == iInvalidRegion) {
	logMessage(1, "Relabelling new region as region %d\n", iRegion);
	setRegionFromCell(cell0, iRegion++);
      }
      if (cell1->getRegion() == iInvalidRegion) {
	logMessage(1, "Relabelling new region as region %d\n", iRegion);
	setRegionFromCell(cell1, iRegion++);
      }
      if (iRegion == iMaxRegion)
	break;

    }
  }

  int
  SurfaceInserter3D::sampleQuasiStructured()
  {

    logMessage(1, "Using quasi-structured sampling\n");
    // now what...
    // lets try a facet based approach. Iterating is stupid.
    DLIList<RefFace*> all_surfaces;
    m_currentRefVolume->ref_faces(all_surfaces);
    int num_surfaces = all_surfaces.size();
    RefFace* surface;
    FacetSurface * facet_surface;
    // lets just do this for every edge to see what happens..
    std::vector<CubitVector> interiorPoints;
    std::vector<CubitVector> curvePoints;

    std::set<CubitVector, CubitVectorCompare> uniquePoints;
    std::set<std::pair<Vert*, Vert*> > edges;
    int numSamplingPoints;

    CubitVector align(1.0, 0.0, 0.0);
    double tol = 1.e-6;
    for (int i = 0; i < num_surfaces; ++i) {
      DLIList<RefEdge*> rEdges;
      DLIList<RefVertex*> rVerts;
      surface = all_surfaces.get_and_step();
      facet_surface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());

      surface->ref_edges(rEdges);
      surface->ref_vertices(rVerts);
      CubitVector edge0, edge1, edge2;
      CubitVector point0, point1, point2;

      for (GR_index_t iF = 0; iF < m_mesh->getNumFaces(); iF++) {
	DLIList<CubitVector*> intersects;
	Face * face = m_mesh->getFace(iF);
	bool qE0 =
	    edges.insert(
		std::make_pair<Vert*, Vert*>(face->getVert(1),
					     face->getVert(2))).second;
	bool qE1 =
	    edges.insert(
		std::make_pair<Vert*, Vert*>(face->getVert(0),
					     face->getVert(2))).second;
	bool qE2 =
	    edges.insert(
		std::make_pair<Vert*, Vert*>(face->getVert(0),
					     face->getVert(1))).second;

	point0.set(face->getVert(0)->getCoords());
	point1.set(face->getVert(1)->getCoords());
	point2.set(face->getVert(2)->getCoords());

	edge0 = ~(point1 - point2);
	edge1 = ~(point0 - point2);
	edge2 = ~(point0 - point1);

	bool q0 = (face->getVert(1)->getVertType() != Vert::eInterior
	    && face->getVert(2)->getVertType() != Vert::eInterior);
	bool q1 = (face->getVert(0)->getVertType() != Vert::eInterior
	    && face->getVert(2)->getVertType() != Vert::eInterior);
	bool q2 = (face->getVert(0)->getVertType() != Vert::eInterior
	    && face->getVert(1)->getVertType() != Vert::eInterior);

	double dot0 = fabs(edge0 % align);
	double dot1 = fabs(edge1 % align);
	double dot2 = fabs(edge2 % align);

	//			double normal[3];
	//			face->calcNormal(normal);
	//			NORMALIZE3D(normal);
	//			printf("normal %f %f %f %d\n",iF,normal[0],normal[1],normal[2]);
	//			if(fabs(normal[0]) < 0.999999 && fabs(normal[1]) < 0.999999 && fabs(normal[2]) < 0.999999) continue;
	// find edge that is tangential to either y or z
	if (dot0 > (1. - tol) && qE0) {
	  //			double edgelength[] = {(point1-point2).length(),(point0-point2).length(),(point1-point0).length()};
	  //			if(edgelength[0] < edgelength[1] && edgelength[0] < edgelength[2]){
	  facet_surface->get_eval_tool()->get_intersections(point1, point2,
							    intersects, true);
	  //				printf("intersects 0 size %d\n",intersects.size());
	}
	else if (dot1 > (1. - tol) && qE1) {
	  //			if(edgelength[1] < edgelength[2] && edgelength[1] < edgelength[0]){
	  facet_surface->get_eval_tool()->get_intersections(point0, point2,
							    intersects, true);
	  //				printf("intersects 1 size %d\n",intersects.size());
	}
	else if (dot2 > (1. - tol) && qE2) {
	  //			if(edgelength[2] < edgelength[1] && edgelength[2] < edgelength[0]){
	  facet_surface->get_eval_tool()->get_intersections(point0, point1,
							    intersects, true);
	  //				printf("intersects 2 size %d\n",intersects.size());
	}
	for (int iInt = 0; iInt < intersects.size(); iInt++) {
	  CubitVector * cv = intersects.get_and_step();
	  //				CubitVector newpoint;
	  //				surface->find_closest_point_trimmed(CubitVector(*cv),newpoint);
	  CubitVector newpoint(*cv);
	  if ((q0 && qE0) || (q1 && qE1) || (q2 && qE2))
	    curvePoints.push_back(newpoint);
	  else
	    interiorPoints.push_back(newpoint);

	}

      }
      // lets get unique ones
      uniquePoints.insert(interiorPoints.begin(), interiorPoints.end());
      interiorPoints.assign(uniquePoints.begin(), uniquePoints.end());
      uniquePoints.clear();
      uniquePoints.insert(curvePoints.begin(), curvePoints.end());
      curvePoints.assign(uniquePoints.begin(), uniquePoints.end());
      uniquePoints.clear();
      writeVTKPointSet(interiorPoints, NULL, "tempInteriors", "");
      writeVTKPointSet(curvePoints, NULL, "tempBdry", "");

      printf("Num sampling points %lu\n",
	     interiorPoints.size() + curvePoints.size());
      numSamplingPoints = interiorPoints.size();
    }
    // need a new surface mesh.

    SurfMeshBuilder * smb = new SurfMeshBuilder(m_currentRefVolume,
						SurfMeshBuilder::FACET,
						0.99 * M_PI, m_mesh,
						CurveDiscretizer::QUASISTRUCT);

    //	smb->sample_all_surfaces();
    // this needs to be taken care of
    // TODO
    smb->insert_interior_points_into_sampled_surf(interiorPoints);
    //	smb->insert_boundary_points_into_sampled_surf(curvePoints);
    smb->output_stitched_surface("tempsurface");

    FaceVec fqes;
    std::set<Subseg*> subsegSet;

    smb->get_subsegs_from_bridges(subsegSet);
    smb->get_restricted_faces(fqes);
    m_surfMeshes[0] = new SurfMesh(fqes, subsegSet);

    writeVTK(*m_surfMeshes[0], "tempsurf");
    return numSamplingPoints;
  }
//void SurfaceInserter3D::offsetSurfByNormals(){
//	DLIList<RefFace*> all_surfaces;
//	m_refVolume->ref_faces(all_surfaces);
//	int num_surfaces = all_surfaces.size();
//	RefFace* surface;
//	FacetSurface * facet_surface;
//
//	CubitVector align(1.0,0.0,0.0);
//	printf("OFFSETTING \n");
//	for(int i = 0; i < num_surfaces; ++i) {
//		DLIList<RefEdge*> rEdges;
//		DLIList<RefVertex*> rVerts;
//		surface = all_surfaces.get_and_step();
//		facet_surface = dynamic_cast<FacetSurface*>(surface->get_surface_ptr());
//		DLIList<CubitPoint*> facetPoints;
//		DLIList<CubitFacet*> facets;
//		facet_surface->get_my_facets(facets,facetPoints);
//		for(int iP = 0; iP < facetPoints.size(); iP++){
//			CubitPoint * point = facetPoints.get_and_step();
//			point->set(CubitVector(2.0*point->x(),2.0*point->y(),2.0*point->z()));
//		}
//	}
//	GraphicsOut go_surf("CGMsurf.vtk");
//	GraphicsOut go_curv("CGMcurv.vtk");
//	go_surf.output_faces();
//	go_curv.output_edges();
//}
// end of class
}// end of namespace
