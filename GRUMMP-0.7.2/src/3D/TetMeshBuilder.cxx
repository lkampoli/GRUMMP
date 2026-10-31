#include "GR_TetMeshBuilder.h"
#include "GR_Constrain.h"
#include "GR_Geometry.h"
#include "GR_Mesh.h"
#include "GR_SurfaceSampler.h"
#include "GR_Subseg.h"
#include "GR_SurfMeshBuilder.h"
#include "GR_SwapDecider.h"
#include "GR_SwapManager.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"
#include "GR_VolMesh.h"
#include "GR_BaseQueueEntry.h"
#include "GR_AdaptPred.h"
#include "GR_RefinementManager3D.h"
#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "FacetEvalTool.hpp"
#include "FacetSurface.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"

#include <algorithm>
#include <deque>
#include <iterator>
#include <map>
#include <set>
#include <vector>
#include <utility>

using std::deque;
using std::map;
using std::multimap;
using std::pair;
using std::set;
using std::vector;

//The following three are private:
TetMeshBuilder::TetMeshBuilder() :
    m_mesh_init(NULL), m_surf_mesh(NULL), m_pDel3D(NULL), m_pSM3D(NULL), m_pSI3D(
    NULL), m_WI(NULL), m_sample_vert_dict(), m_subseg_vert_dict()
{
  assert(0);
}

TetMeshBuilder::TetMeshBuilder(const TetMeshBuilder&) :
    m_mesh_init(NULL), m_surf_mesh(NULL), m_pDel3D(NULL), m_pSM3D(NULL), m_pSI3D(
    NULL), m_WI(NULL), m_sample_vert_dict(), m_subseg_vert_dict()
{
  assert(0);
}

TetMeshBuilder&
TetMeshBuilder::operator=(const TetMeshBuilder&)
{
  assert(0);
  return *this;
}

//Only public constructor
TetMeshBuilder::TetMeshBuilder(SurfMeshBuilder* surf_mesh) :
    m_mesh_init(NULL), m_surf_mesh(surf_mesh), m_pDel3D(NULL), m_pSM3D(
    NULL), m_pSI3D(NULL), m_WI(NULL), m_sample_vert_dict(), m_subseg_vert_dict()
{

  init_mesh(dynamic_cast<RefVolume*>(surf_mesh->getEntity()));

  m_pDel3D = new GRUMMP::DelaunaySwapDecider3D(true);
  m_pSM3D = new GRUMMP::SwapManager3D(m_pDel3D, m_mesh_init);
  m_pSI3D = new GRUMMP::SwappingInserter3D(m_mesh_init, m_pSM3D);
  m_pSI3D->setForcedInsertion(true);
  // uses new class of watson inserter
  m_WI = new GRUMMP::WatsonInserter(m_mesh_init);

}

//Destructor
TetMeshBuilder::~TetMeshBuilder()
{
  if (m_WI)
    delete m_WI;
  if (m_pSI3D)
    delete m_pSI3D;
  if (m_pSM3D)
    delete m_pSM3D;
  if (m_pDel3D)
    delete m_pDel3D;
  // Here's the problem: this is sometimes also used externally.  Shared
  // pointers are looking better all the time...
//	if (m_mesh_init) delete m_mesh_init;
}

void
TetMeshBuilder::destroyInitMesh()
{
  delete m_mesh_init;
}

bool
TetMeshBuilder::buildMesh(set<Subseg*>* subsegs)
{

  insert_sample_verts();

  bool success = recover_boundary();

  assert(m_mesh_init->isValid());
  logMessage(2, "Passed validity check.\n");

  //m_surf_mesh->output_stitched_surface("surf.mesh");

#ifndef NDEBUG
  int bef_num_verts = m_mesh_init->getNumVerts();
  m_mesh_init->purgeAllEntities();
  int aft_num_verts = m_mesh_init->getNumVerts();
  assert(bef_num_verts == aft_num_verts);
  assert(m_mesh_init->isValid());
  logMessage(2, "Passed validity check.\n");
#endif

  success = clean_exterior();
  if (!success) {
    logMessage(1, "Unable to recover closed surface\n");
    return false;
  }
  //After the boundary is recovered, swap the interior
  //to obtain constrained Delaunay tets. We are not sure
  //this is guaranteed to work in theory, but in practice,
  //it never failed.

#ifndef NDEBUG
  //If the number of verts does not change during swap, we are safe.
  GR_index_t dbg_num_verts = m_mesh_init->getNumVerts();
#endif

  GRUMMP::DelaunaySwapDecider3D SDDel(false);
  GRUMMP::SwapManager3D DelSwap(&SDDel, m_mesh_init);
  DelSwap.swapAllFaces();

#ifndef NDEBUG
  assert(dbg_num_verts == m_mesh_init->getNumVerts());
#endif

#ifndef NDEBUG

  //I want to make sure that all the pointers to vertices
  //remain constant through the mesh purge. If the purge is called, then
  //the pointers to the vertices in m_subseg_vert_dict will
  //remain valid i.e. they will still point to the correct object.

  //	std::map<Vert*, Vert*> vert_map;
  //	mesh->purgeAllEntities(&vert_map);
  //	std::map<Vert*, Vert*>::const_iterator
  //	it = vert_map.begin(), it_end = vert_map.end();
  //	for( ; it != it_end; ++it) assert(it->first == it->second);

#else

  //	m_mesh_init->purgeAllEntities();

#endif

  assert(m_mesh_init->isValid());
  assert(mesh_is_const_delaunay());

  //If we want to get the subsegments, then output them.
  if (subsegs)
    get_subseg_copy(subsegs);
  if (!success)
    return false;
  return true;
}
bool
TetMeshBuilder::mesh_is_const_delaunay() const
{

  //still have to figure out how to reliably detect
  //unswappable configs so that they do not affect refinement.
  return true;

  //Check that the mesh is truely constrained Delaunay at this point.
  //For every cell, check the neighboring cells to see if the vertex
  //opposing the face they share is inside the cell's circumsphere.

  assert(m_mesh_init->isValid());

  Cell* dbg_cell1, *dbg_cell2;
  Face* dbg_face;
  Vert* dbg_vert;

  for (int i = m_mesh_init->getNumCells() - 1; i >= 0; --i) {

    dbg_cell1 = m_mesh_init->getCell(i);
    if (dbg_cell1->isDeleted())
      continue;
    assert(dbg_cell1->getType() == Cell::eTet);

    for (int j = 0; j < 4; ++j) {

      dbg_face = dbg_cell1->getFace(j);
      assert(!dbg_face->isDeleted());
      assert(dbg_face->hasCell(dbg_cell1));
      dbg_cell2 =
	  dbg_face->getLeftCell() == dbg_cell1 ?
	      dbg_face->getRightCell() : dbg_face->getLeftCell();
      assert(dbg_face->hasCell(dbg_cell2));
      assert(dbg_cell1 != dbg_cell2);
      if (dbg_cell2->getType() != Cell::eTet)
	continue;
      dbg_vert = dbg_cell2->getOppositeVert(dbg_face);

      //Check if the vertex is in the circumsphere.
      TetCell* tet = dynamic_cast<TetCell*>(dbg_cell1);
      assert(tet->isValid() && !tet->isDeleted());
      if (tet->circumscribesPointAdaptive(dbg_vert)) {

	//Because of tie-break procedure in the swapping routines, it is possible that
	//if the vertex is close enough to the circumsphere, the predicate used
	//in TetCell does not provide the expected result, bail out if this is the case.
	if (isInsphere(dbg_cell1->getVert(0), dbg_cell1->getVert(0),
		       dbg_cell1->getVert(0), dbg_cell1->getVert(0), dbg_vert)
	    == 0)
	  continue;

	Vert* pVVertA, *pVVertB, *pVVertC, *pVVertD, *pVVertE, *pVPivot0,
	    *pVPivot1, *pVOther;
	TetCell* apTCTets[4];
	int iNTets;
	printf(
	    "Face category = %d\n",
	    dynamic_cast<TriFace*>(dbg_face)->categorizeFace(pVVertA, pVVertB,
							     pVVertC, pVVertD,
							     pVVertE, apTCTets,
							     iNTets, pVPivot0,
							     pVPivot1,
							     pVOther));

	vFatalError("Could not build a constrained Delaunay Tetrahedralization",
		    "TetMeshBuilder::mesh_is_const_delaunay()");
	abort();

      }

    }

  }

  logMessage(1, "The tetrahedralization is constrained Delaunay\n");
  return true;

}

VolMesh*
TetMeshBuilder::get_mesh() const
{
  return m_mesh_init;
}

void
TetMeshBuilder::insert_sample_verts()
{
  qAdaptPred = true;
  assert(m_mesh_init);

  //Stick all the sampling vertices into the mesh.

  typedef pair<VertexDict::iterator, bool> InsertResult;

  { //We begin with subseg vertices

    vector<Subseg*> subsegs;
    m_surf_mesh->get_subsegs(subsegs);
    vector<Subseg*>::iterator it_sub = subsegs.begin(), it_sub_end =
	subsegs.end();
    for (; it_sub != it_sub_end; ++it_sub) {

      Subseg* subseg = *it_sub;
      Vert *sub_beg_vert = subseg->get_beg_vert(), *sub_end_vert =
	  subseg->get_end_vert();

      InsertResult ins_beg = m_subseg_vert_dict.insert(
	  std::make_pair(sub_beg_vert, static_cast<Vert*>(NULL)));
      InsertResult ins_end = m_subseg_vert_dict.insert(
	  std::make_pair(sub_end_vert, static_cast<Vert*>(NULL)));

      Vert *new_beg_vert = ins_beg.first->second, *new_end_vert =
	  ins_end.first->second;

      if (ins_beg.second) {

	if (ins_end.second) {
	  assert(!new_beg_vert && !new_end_vert);
	  m_WI->computeHull(sub_beg_vert->getCoords(),
			    getSeedCellFromVerts(sub_beg_vert));
	  new_beg_vert = ins_beg.first->second = m_WI->insertPointInHull(
	      m_mesh_init->createVert(sub_beg_vert));
	  m_WI->computeHull(sub_end_vert->getCoords(),
			    getSeedCellFromVerts(sub_end_vert, new_beg_vert));
	  new_end_vert = ins_end.first->second = m_WI->insertPointInHull(
	      m_mesh_init->createVert(sub_end_vert));
	}
	else {
	  assert(!new_beg_vert && new_end_vert);
	  m_WI->computeHull(sub_beg_vert->getCoords(),
			    getSeedCellFromVerts(sub_beg_vert, new_end_vert));
	  new_beg_vert = ins_beg.first->second = m_WI->insertPointInHull(
	      m_mesh_init->createVert(sub_beg_vert));
	}

      }

      else {

	if (ins_end.second) {
	  assert(new_beg_vert && !new_end_vert);
	  m_WI->computeHull(sub_end_vert->getCoords(),
			    getSeedCellFromVerts(sub_end_vert, new_beg_vert));
	  new_end_vert = ins_end.first->second = m_WI->insertPointInHull(
	      m_mesh_init->createVert(sub_end_vert));
	}
	else {
	  assert(new_beg_vert && new_end_vert);
	}

      }

      assert(new_beg_vert && new_end_vert); //makes sure that vertices exist for the subseg.

      vector<SubsegBridge*> subseg_bridges;
      subseg->get_children(subseg_bridges);
      vector<SubsegBridge*>::iterator it_bridge = subseg_bridges.begin(),
	  it_bridge_end = subseg_bridges.end();

      for (; it_bridge != it_bridge_end; ++it_bridge) {

	SubsegBridge* bridge = *it_bridge;
	Vert* bridge_beg_vert = bridge->get_vert(0);
	Vert* bridge_end_vert = bridge->get_vert(1);

	assert(
	    iFuzzyComp( dDIST3D(sub_beg_vert->getCoords(), bridge_beg_vert->getCoords() ), 0.) == 0);
	assert(
	    iFuzzyComp( dDIST3D(sub_end_vert->getCoords(), bridge_end_vert->getCoords() ), 0.) == 0);

	m_sample_vert_dict.insert(
	    std::make_pair(bridge_beg_vert, new_beg_vert));
	m_sample_vert_dict.insert(
	    std::make_pair(bridge_end_vert, new_end_vert));

      }

    }

  }

  { //Proceed with face samplers, being careful not to insert the same vertex twice.

    typedef vector<SurfaceSampler*> SamplerVec;
    typedef vector<FaceQueueEntry*> SurfTriVec;

    SamplerVec face_samplers;
    m_surf_mesh->get_face_samplers(face_samplers);
    SamplerVec::iterator it_sampler = face_samplers.begin(), it_sampler_end =
	face_samplers.end();

    for (; it_sampler != it_sampler_end; ++it_sampler) {

      SurfaceSampler* sampler = *it_sampler;

      SurfTriVec surf_tri;
      sampler->get_restricted_faces(surf_tri);
      SurfTriVec::iterator it_tri = surf_tri.begin(), it_tri_end =
	  surf_tri.end();

      for (; it_tri != it_tri_end; ++it_tri) {

	Face * face = (*it_tri)->getFace();
	delete *it_tri;
	assert(face->getNumVerts() == 3);
	InsertResult ins_res[3];

	for (int i = 0; i < 3; ++i) {
	  ins_res[i] = m_sample_vert_dict.insert(
	      std::make_pair(face->getVert(i), static_cast<Vert*>(NULL)));
	}

	if (!ins_res[0].second && !ins_res[1].second && !ins_res[2].second)
	  continue;

	Vert* new_vert[] =
	  { ins_res[0].first->second, ins_res[1].first->second,
	      ins_res[2].first->second };
	// This is modified code from Serge's work
	// In each case, the vertes are inserted based on newly created vertices,
	// To account for the fact that we no longer do memory the same way
	// the sample vert dict erases the entry corresponding to the face-queued vert and replaces it
	// with a new entry that simply points to itself
	// This is a bit redundant, but the resulting vert dict is correct

	if (!new_vert[0]) {
	  if (!new_vert[1]) {
	    if (!new_vert[2]) {
	      assert(!new_vert[0] && !new_vert[1] && !new_vert[2]);
	      m_WI->computeHull(face->getVert(2)->getCoords(),
				getSeedCellFromVerts(face->getVert(2)));
	      new_vert[2] = ins_res[2].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(2)));
	      m_WI->computeHull(
		  face->getVert(1)->getCoords(),
		  getSeedCellFromVerts(face->getVert(1), new_vert[2]));
	      new_vert[1] = ins_res[1].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(1)));
	      m_WI->computeHull(
		  face->getVert(0)->getCoords(),
		  getSeedCellFromVerts(face->getVert(0), new_vert[2],
				       new_vert[1]));
	      new_vert[0] = ins_res[0].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(0)));
	      m_sample_vert_dict.erase(face->getVert(0));
	      m_sample_vert_dict[new_vert[0]] = new_vert[0];
	      m_sample_vert_dict.erase(face->getVert(1));
	      m_sample_vert_dict[new_vert[1]] = new_vert[1];
	      m_sample_vert_dict.erase(face->getVert(2));
	      m_sample_vert_dict[new_vert[2]] = new_vert[2];
	    }
	    else {
	      assert(!new_vert[0] && !new_vert[1] && new_vert[2]);
	      m_WI->computeHull(
		  face->getVert(1)->getCoords(),
		  getSeedCellFromVerts(face->getVert(1), new_vert[2]));
	      new_vert[1] = ins_res[1].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(1)));
	      m_WI->computeHull(
		  face->getVert(0)->getCoords(),
		  getSeedCellFromVerts(face->getVert(0), new_vert[2],
				       new_vert[1]));
	      new_vert[0] = ins_res[0].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(0)));
	      m_sample_vert_dict.erase(face->getVert(0));
	      m_sample_vert_dict[new_vert[0]] = new_vert[0];
	      m_sample_vert_dict.erase(face->getVert(1));
	      m_sample_vert_dict[new_vert[1]] = new_vert[1];
	    }
	  }
	  else {
	    if (!new_vert[2]) {
	      assert(!new_vert[0] && new_vert[1] && !new_vert[2]);
	      m_WI->computeHull(
		  face->getVert(2)->getCoords(),
		  getSeedCellFromVerts(face->getVert(2), new_vert[1]));
	      new_vert[2] = ins_res[2].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(2)));
	      m_WI->computeHull(
		  face->getVert(0)->getCoords(),
		  getSeedCellFromVerts(face->getVert(0), new_vert[1],
				       new_vert[2]));
	      new_vert[0] = ins_res[0].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(0)));
	      m_sample_vert_dict.erase(face->getVert(0));
	      m_sample_vert_dict[new_vert[0]] = new_vert[0];
	      m_sample_vert_dict.erase(face->getVert(2));
	      m_sample_vert_dict[new_vert[2]] = new_vert[2];
	    }
	    else {
	      assert(!new_vert[0] && new_vert[1] && new_vert[2]);
	      m_WI->computeHull(
		  face->getVert(0)->getCoords(),
		  getSeedCellFromVerts(face->getVert(0), new_vert[1],
				       new_vert[2]));
	      new_vert[0] = ins_res[0].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(0)));
	      m_sample_vert_dict.erase(face->getVert(0));
	      m_sample_vert_dict[new_vert[0]] = new_vert[0];
	    }
	  }
	}

	else {
	  if (!new_vert[1]) {
	    if (!new_vert[2]) {
	      assert(new_vert[0] && !new_vert[1] && !new_vert[2]);
	      m_WI->computeHull(
		  face->getVert(2)->getCoords(),
		  getSeedCellFromVerts(face->getVert(2), new_vert[0]));
	      new_vert[2] = ins_res[2].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(2)));
	      m_WI->computeHull(
		  face->getVert(1)->getCoords(),
		  getSeedCellFromVerts(face->getVert(1), new_vert[2],
				       new_vert[0]));
	      new_vert[1] = ins_res[1].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(1)));

	      m_sample_vert_dict.erase(face->getVert(1));
	      m_sample_vert_dict[new_vert[1]] = new_vert[1];
	      m_sample_vert_dict.erase(face->getVert(2));
	      m_sample_vert_dict[new_vert[2]] = new_vert[2];
	    }
	    else {
	      assert(new_vert[0] && !new_vert[1] && new_vert[2]);
	      m_WI->computeHull(
		  face->getVert(1)->getCoords(),
		  getSeedCellFromVerts(face->getVert(1), new_vert[2],
				       new_vert[0]));
	      new_vert[1] = ins_res[1].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(1)));

	      m_sample_vert_dict.erase(face->getVert(1));
	      m_sample_vert_dict[new_vert[1]] = new_vert[1];

	    }
	  }
	  else {
	    if (!new_vert[2]) {
	      assert(new_vert[0] && new_vert[1] && !new_vert[2]);
	      m_WI->computeHull(
		  face->getVert(2)->getCoords(),
		  getSeedCellFromVerts(face->getVert(2), new_vert[1],
				       new_vert[0]));
	      new_vert[2] = ins_res[2].first->second = m_WI->insertPointInHull(
		  m_mesh_init->createVert(face->getVert(2)));

	      m_sample_vert_dict.erase(face->getVert(2));
	      m_sample_vert_dict[new_vert[2]] = new_vert[2];
	    }
	    else {
	      assert(new_vert[0] && new_vert[1] && new_vert[2]);
	    }
	  }
	}

	assert(new_vert[0] && new_vert[1] && new_vert[2]);
      }
      surf_tri.clear();
    }
  }

  qAdaptPred = false;
  typedef vector<SurfaceSampler*> SamplerVec;
  SamplerVec face_samplers;
  m_surf_mesh->get_face_samplers(face_samplers);
  SamplerVec::iterator it_sampler = face_samplers.begin(), it_sampler_end =
      face_samplers.end();

  // Here, all restricted faces are correct in the sampler, but the vertices don't agree
  // since they have been mapped.  This simply fixes the vertices.
  for (; it_sampler != it_sampler_end; ++it_sampler) {
    SurfaceSampler* sampler = *it_sampler;
    sampler->reset_restricted_faces();
  }

}

void
TetMeshBuilder::init_mesh(RefVolume * const refVolume)
{

  assert(m_surf_mesh);

  //Obtain all sampling vertices, compute a bounding box and build an initial mesh

  CubitBox box = refVolume->bounding_box();

  CubitVector mini(box.minimum()), maxi(box.maximum()), mean(0., 0., 0.), diff(
      0., 0., 0.);

  vector<Vert*> subseg_verts, sample_verts;
  m_surf_mesh->get_subseg_verts(subseg_verts);
  m_surf_mesh->get_sample_verts(sample_verts);

  //   printf("number of subseg verts = %d\n", static_cast<int>(subseg_verts.size()) );

  std::for_each(subseg_verts.begin(), subseg_verts.end(),
		ComputeBoundingBox(&mini, &maxi));
  std::for_each(sample_verts.begin(), sample_verts.end(),
		ComputeBoundingBox(&mini, &maxi));

  mean.set(0.5 * (mini + maxi));
  diff.set(maxi - mini);

  // Make the region roughly cubical
  double max_aspect = 0.1;
  if (diff.x() < max_aspect * diff.y())
    diff.x(max_aspect * diff.y());
  if (diff.x() < max_aspect * diff.z())
    diff.x(max_aspect * diff.z());
  if (diff.y() < max_aspect * diff.x())
    diff.y(max_aspect * diff.x());
  if (diff.y() < max_aspect * diff.z())
    diff.y(max_aspect * diff.z());
  if (diff.z() < max_aspect * diff.x())
    diff.z(max_aspect * diff.x());
  if (diff.z() < max_aspect * diff.y())
    diff.z(max_aspect * diff.y());

  mini.set(mean - diff);
  maxi.set(mean + diff);

  //   printf("mini: %lf %lf %lf\n", mini.x(), mini.y(), mini.z());
  //   printf("maxi: %lf %lf %lf\n", maxi.x(), maxi.y(), maxi.z());

  m_mesh_init = new VolMesh(CubitBox(mini, maxi));

  assert(m_mesh_init->isSimplicial());

}

bool
TetMeshBuilder::recover_boundary()
{

  logMessage(1, "Recovering surface mesh...\n");

  set<RecoverEdge> subcurves, surface_edges;
  set<RecoverFace> subfacets;
  set<RecoverFace>::iterator it_facet;

  RecoverFace my_face;
  vector<RecoverFace> faces_not_recovered;

  if (is_bdry_constrained(subcurves, surface_edges, subfacets))
    return true;

  while (!subcurves.empty() || !surface_edges.empty() || !subfacets.empty()) {

    if (!subcurves.empty() || !surface_edges.empty()) {
      if (recover_all_edges(subcurves, surface_edges))
	return true;
      if (is_bdry_constrained(subcurves, surface_edges, subfacets))
	return true;
      assert(subcurves.empty() && surface_edges.empty());
    }
    else if (subfacets.empty()) {
      assert(subcurves.empty() && surface_edges.empty());
      return true;
    }

    for (it_facet = subfacets.begin(); it_facet != subfacets.end();
	++it_facet) {
      my_face = *it_facet;
      if (recover_face_flip(my_face) == -1)
	faces_not_recovered.push_back(my_face);
    }

    //If we get here, face flipping alone has failed to recover the boundary
    //so let's try inserting for get the missing faces back.

    while (!faces_not_recovered.empty()) {
      my_face = faces_not_recovered.back();
      recover_face_split(my_face);
      faces_not_recovered.pop_back();
    }

    if (is_bdry_constrained(subcurves, surface_edges, subfacets))
      return true;

  }

  return true;

}

bool
TetMeshBuilder::recover_all_edges(set<RecoverEdge>& bdry_edges_to_recover,
				  set<RecoverEdge>& face_edges_to_recover)
{

  if (bdry_edges_to_recover.empty() && face_edges_to_recover.empty())
    return false;

  typedef set<RecoverEdge> EdgeSet;
  typedef set<RecoverFace> FaceSet;
  typedef EdgeSet::iterator ItRE;

  EdgeSet edges(bdry_edges_to_recover);
  FaceSet faces; //this is a dummy set.

  bool allow_insertion = false;
  int num_passes, num_fruitless_passes, num_inserts, num_changes;
  int edge_recovered;

  m_mesh_init->disallowSwapRecursion();

  if (edges.empty())
    edges.insert(face_edges_to_recover.begin(), face_edges_to_recover.end());

  do {

    num_passes = 0;
    num_fruitless_passes = 0;

    do {

      assert(!edges.empty());

      ++num_passes;
      ++num_fruitless_passes;
      num_changes = 0;
      num_inserts = 0;

      ItRE it_edge = edges.begin(), it_edge_end = edges.end(), it_edge_tmp;

      for (; it_edge != it_edge_end;) {

	if (allow_insertion) {
	  edge_recovered = 0;
	  // 	  edge_recovered = recover_edge_split(*it_edge);
	  // 	  ++num_inserts;
	}
	else
	  edge_recovered = recover_edge_flip(*it_edge);

	if (edge_recovered == 1) {
	  num_fruitless_passes = 0;
	  it_edge_tmp = it_edge++;
	  edges.erase(it_edge_tmp);
	  ++num_changes;
	}
	else {
	  if (edge_recovered == -1)
	    ++num_changes;
	  ++it_edge;
	}

      } // Done checking non-constrained segments

      if (edges.empty()) {

	if (is_bdry_constrained(bdry_edges_to_recover, face_edges_to_recover,
				faces))
	  return true;

	edges.clear();
	if (bdry_edges_to_recover.empty())
	  edges.insert(face_edges_to_recover.begin(),
		       face_edges_to_recover.end());
	else
	  edges.insert(bdry_edges_to_recover.begin(),
		       bdry_edges_to_recover.end());

	// If swapping fixed all the edges, don't turn off insertion,
	// or you could loop forever trying to do surface recovery.
	if (allow_insertion && num_inserts != 0) {
	  assert(num_fruitless_passes == 0);
	  num_passes = 0;
	  allow_insertion = false;
	}

      }

    }
    while (num_fruitless_passes < 2 && num_changes && (!edges.empty())
	&& num_passes < 10);

    if (edges.empty())
      return false;
    if (!m_mesh_init->areBdryChangesAllowed()) {
      assert(0);
    }

    if (allow_insertion) {

      // Since allowing insertion wasn't enough, force a split of
      // all edges not yet constrained.  Enable swap recursion in
      // the volume mesh while splitting.

      ItRE it_edge = edges.begin(), it_edge_end = edges.end(), it_edge_tmp;

      m_mesh_init->allowSwapRecursion();

      for (; it_edge != it_edge_end; ++it_edge) {

	recover_edge_split(*it_edge);
      }

      m_mesh_init->disallowSwapRecursion();

      allow_insertion = false;

      if (is_bdry_constrained(bdry_edges_to_recover, face_edges_to_recover,
			      faces))
	return true;

      edges.clear();
      if (bdry_edges_to_recover.empty())
	edges.insert(face_edges_to_recover.begin(),
		     face_edges_to_recover.end());
      else
	edges.insert(bdry_edges_to_recover.begin(),
		     bdry_edges_to_recover.end());

    }

    else {
      allow_insertion = true;
    }

  }
  while (!edges.empty() && num_passes < 20);

  if (!edges.empty())
    vFatalError("Failed to recover edges\n", "Edge recovery");

  return false;

}

int
TetMeshBuilder::recover_edge_flip(const RecoverEdge& edge_to_recover)
{

  Vert* v0 = edge_to_recover.get_v0();
  Vert* v1 = edge_to_recover.get_v1();

  int num_passes = 0, change, num_skips, ret_value;
  Face* faces_skipped[20];

  do {

    ++num_passes;
    logMessage(4, "Pass %d to remove edge from v0 = %p to v1 = %p\n",
	       num_passes, v0, v1);

    bool hit_edge = false;
    int ream_result = m_mesh_init->clearPipe(v0, v1, hit_edge, faces_skipped,
					     num_skips);

    if (ream_result == 1)
      return 1; // Pipe is unclogged
    if (ream_result == 0) {
      change = 1;
      ret_value = -1;
    } //change made
    else {
      assert(ream_result == -1);
      ret_value = -1;
    }

    if (hit_edge) { //try from the other end
      hit_edge = false;
      num_skips = 0;

      //pasto here???
      ream_result = m_mesh_init->clearPipe(v0, v1, hit_edge, faces_skipped,
					   num_skips);

      if (ream_result == 1)
	return 1; // Pipe is unclogged; not bloody likely
      if (ream_result == 0) {
	change = 1;
	ret_value = -1;
      } //change made
      else {
	assert(ream_result == -1);
	ret_value = -1;
      }
    }

    // If an edge was hit, try to remove it
    if (hit_edge)
      logMessage(3, "Unresolved edge hit!  Edge not recovered.\n");

  }
  while (change && num_passes < 5);

  return ret_value;

}

int
TetMeshBuilder::recover_edge_split(const RecoverEdge& edge_to_recover)
{

  static int global_num_splits = 0;

  //First, split the subseg in m_surf_mesh. One (or possibly many) vertices
  //will be inserted. These will also correspond to sample verts inserted in
  //corresponding samplers. We have to make sure that all these vertices are
  //properly put in our dictionaries.

  VertVec new_verts;

  Subseg* subseg_to_split = edge_to_recover.get_subseg();
  TriFace* face_to_split;
  RefFace* insert_surf;

  if (!subseg_to_split) { //Case of an in-surface edge

    face_to_split = edge_to_recover.get_face();
    insert_surf = edge_to_recover.get_surface();
    assert(face_to_split && insert_surf);
    //		face_to_split->printFaceInfo();
    //		printf("sample vert dict sizeA %d\n",m_sample_vert_dict.size());
    if (!face_to_split->isRestrictedDelaunay())
      return 0;
    if (face_to_split->isDeleted())
      return 0;

    global_num_splits += m_surf_mesh->split_face_and_sample(
	new FaceQueueEntry(face_to_split), insert_surf, &new_verts);

  }

  else { //Case of a subsegment or subcurve.

    assert(subseg_to_split);

    if (subseg_to_split->is_deleted())
      return 0;

    global_num_splits += m_surf_mesh->split_subseg_and_sample(subseg_to_split,
							      &new_verts);

  }

  insert_new_verts(new_verts);

  printf(
      "Needed to split edge to recover, total number of edge insertions = %d\n",
      global_num_splits);
  //m_surf_mesh->reset_restricted_faces();
  return 1;

}

int
TetMeshBuilder::recover_face_flip(const RecoverFace& face_to_recover)
{

  //@@ Find all faces and edges violating the given surface face.
  Vert* vert0 = face_to_recover.get_v0();
  Vert* vert1 = face_to_recover.get_v1();
  Vert* vert2 = face_to_recover.get_v2();
  assert(vert0 && vert1 && vert2);

  // First, make sure that all the edges of the surface face are
  // present.  If they aren't, this means that they are interior edges
  // to a surface facet (and so weren't recovered earlier when doing
  // bdry subsegments); and that there is a vertex encroaching on them
  // (which makes it so that the Delaunay tetrahedralization doesn't
  // match the constrained DT).  If the edge is missing, recover it
  // without point insertion (this should always be possible, but note
  // that there's still a check to be sure that it works).
  int edge_a_ok = recover_edge_flip(RecoverEdge(vert0, vert1));
  int edge_b_ok = recover_edge_flip(RecoverEdge(vert1, vert2));
  int edge_c_ok = recover_edge_flip(RecoverEdge(vert2, vert0));
  if (edge_a_ok != 1 || edge_b_ok != 1 || edge_c_ok != 1) {
    printf("cannot recover edge: %d %d %d\n", edge_a_ok, edge_b_ok, edge_c_ok);
    return recover_face_split(face_to_recover);
    //     printf("cannot recover edge: %d %d %d\n", edge_a_ok, edge_b_ok, edge_c_ok);
    //     return -1;
    //     vFatalError("Failed to recover a non-subsegment bdry edge.",
    // 		"3D constrained Delaunay tetrahedralization.");

  }

  set<BadEdge> bad_edges;
  set<Face*> bad_faces;
  m_mesh_init->identifyBadEdgesAndFaces(vert0, vert1, vert2, bad_edges,
					bad_faces);

  //@@ For each face in the list, swap it if possible.
  // Swap only if the face is of type T32 and all of its equatorial
  // points are on the same side of or in the plane of the surface face.
  bool change, global_change = false;

  do {

    change = false;

    deque<Face*> bad_face_queue;
    std::copy(bad_faces.begin(), bad_faces.end(),
	      std::back_inserter(bad_face_queue));

    while (!bad_face_queue.empty()) {

      Face* face = bad_face_queue.front();
      bad_face_queue.pop_front();
      if (face->isDeleted())
	continue;

      TriFace* triface = dynamic_cast<TriFace*>(face);
      assert(triface);

      Vert *vertA, *vertB, *vertC, *vertD, *vertE, *vertJ0, *vertJ1, *vertJ2;

      int num_tets;
      TetCell* tet_cells[4];

      eFaceCat face_cat = triface->categorizeFace(vertA, vertB, vertC, vertD,
						  vertE, tet_cells, num_tets,
						  vertJ0, vertJ1, vertJ2);

      if (face_cat == eT32) {

	assert(num_tets == 3);

	Face *face_delete[3];
	face_delete[0] = findCommonFace(tet_cells[0], tet_cells[1]);
	face_delete[1] = findCommonFace(tet_cells[1], tet_cells[2]);
	face_delete[2] = findCommonFace(tet_cells[2], tet_cells[0]);
	assert(face_delete[0]->isValid());
	assert(face_delete[1]->isValid());
	assert(face_delete[2]->isValid());

	int orient_D = checkOrient3D(vert0, vert1, vert2, vertD);
	int orient_E = checkOrient3D(vert0, vert1, vert2, vertE);
	if (orient_D * orient_E != -1)
	  continue;

	assert(
	    triface == face_delete[0] || triface == face_delete[1]
		|| triface == face_delete[2]);

	change = global_change = true;

#ifndef NDEBUG
	assert(
	    m_mesh_init->reconfTet32_deprecated(tet_cells, triface, vertA,
						vertB, vertC, vertD, vertE)
		== 1);
#else
	m_mesh_init->reconfTet32_deprecated(tet_cells, triface, vertA, vertB, vertC, vertD, vertE);
#endif

	// Check each of the three faces to be sure that either
	// they've been deleted or all of their verts lie on the
	// same side of the target face (including ties).

	for (int i = 0; i < 3; i++) {

	  Face *face_test = face_delete[i];

	  if (face_test->isDeleted())
	    bad_faces.erase(face_test);
	  else {
	    Vert *vert_test[] =
	      { face_test->getVert(0), face_test->getVert(1),
		  face_test->getVert(2) };
	    int orient[] =
	      { checkOrient3D(vert0, vert1, vert2, vert_test[0]), checkOrient3D(
		  vert0, vert1, vert2, vert_test[1]), checkOrient3D(
		  vert0, vert1, vert2, vert_test[2]) };

	    if ((orient[0] * orient[1] != -1) && (orient[1] * orient[2] != -1)
		&& (orient[2] * orient[0] != -1))
	      bad_faces.erase(face_test);
	  }

	}

      } // T32 face

    } // End loop over all bad faces

    m_mesh_init->identifyBadEdgesAndFaces(vert0, vert1, vert2, bad_edges,
					  bad_faces);

    // Do this until no more progress can be made by simple swapping or
    // until all offending edges are removed

  }
  while (change && !bad_faces.empty());

  if (bad_faces.empty())
    return 1;
  else
    return -1;

}

int
TetMeshBuilder::recover_face_split(const RecoverFace& face_to_recover)
{

  static int global_num_splits = 0;

  RefFace* insert_surf = face_to_recover.get_surface();
  TriFace* face_to_split = face_to_recover.get_face();
  assert(insert_surf && face_to_split);

  if (face_to_split->isDeleted())
    return 0;

  VertVec new_verts;
  global_num_splits += m_surf_mesh->split_face_and_sample(
      new FaceQueueEntry(face_to_split), insert_surf, &new_verts);

  insert_new_verts(new_verts);

  logMessage(
      1,
      "Needed to split face to recover, total number of face insertions = %d\n",
      global_num_splits);

  return 1;

}

void
TetMeshBuilder::insert_new_verts(const VertVec& new_verts)
{

  assert(m_mesh_init->isSimplicial());

  typedef multimap<Vert*, Vert*> VertMap;
  VertMap new_sub_verts;

  Vert *new_vert, *sub_vert, *sam_vert;

  VertVec::const_iterator itv = new_verts.begin(), itv_end = new_verts.end();

  for (; itv != itv_end; ++itv) {

    sub_vert = itv->first;
    sam_vert = itv->second;
    assert(sam_vert);

    if (sub_vert) {
      new_sub_verts.insert(std::make_pair(sub_vert, sam_vert));

    }
    else {
      new_vert = m_mesh_init->createVert(sam_vert);

      // The restricted faces have changed to the new verts, this is a hack for now
      m_sample_vert_dict.insert(std::make_pair(new_vert, new_vert));
      insert_vert_and_swap(new_vert);

      new_vert->setType(Vert::ePseudoSurface);

      assert(dDIST3D(sam_vert->getCoords(), new_vert->getCoords()) < 1.e-12);
      assert(sam_vert->getParentTopology() == new_vert->getParentTopology());
      assert(sam_vert->getVertType() == new_vert->getVertType());
      assert(new_vert->getVertType() == Vert::ePseudoSurface);
    }

  }

  VertMap::iterator itm, itm_beg, itm_end;

  while (!new_sub_verts.empty()) {

    itm_beg = new_sub_verts.begin();
    itm_end = new_sub_verts.upper_bound(itm_beg->first);

    sub_vert = itm_beg->first;
    new_vert = m_mesh_init->createVert(sub_vert);

    m_subseg_vert_dict.insert(std::make_pair(sub_vert, new_vert));
    insert_vert_and_swap(new_vert);

    new_vert->setType(Vert::eBdryCurve);

    for (itm = itm_beg; itm != itm_end; itm++) {

      sam_vert = itm->second;
      assert(sub_vert && sam_vert);
      m_sample_vert_dict.insert(std::make_pair(sam_vert, new_vert));

      assert(dDIST3D(sub_vert->getCoords(), sam_vert->getCoords()) < 1.e-12);
      assert(sub_vert->getVertType() == sam_vert->getVertType());
      assert(sub_vert->getParentTopology() == sam_vert->getParentTopology());

    }

    assert(dDIST3D(sub_vert->getCoords(), new_vert->getCoords()) < 1.e-12);
    assert(sub_vert->getParentTopology() == new_vert->getParentTopology());
    assert(sub_vert->getVertType() == new_vert->getVertType());
    assert(new_vert->getVertType() == Vert::eBdryCurve);
    //
    new_sub_verts.erase(itm_beg, itm_end);

  }
  typedef vector<SurfaceSampler*> SamplerVec;
  SamplerVec face_samplers;
  m_surf_mesh->get_face_samplers(face_samplers);
  SamplerVec::iterator it_sampler = face_samplers.begin(), it_sampler_end =
      face_samplers.end();
  // Here, all restricted faces are correct in the sampler, but the vertices don't agree
  // because of the swapping. this simply fixes the vertices.
  for (; it_sampler != it_sampler_end; ++it_sampler) {
    SurfaceSampler* sampler = *it_sampler;
    sampler->reset_restricted_faces();
  }

}

bool
TetMeshBuilder::insert_vert_and_swap(Vert* const vertex)
{

  assert(m_mesh_init->getSwapType() == eDelaunay);

  //vertex is created and initialized, but has not been inserted in the mesh yet.
  //Well, now is the time to insert it. With all the swaps, our mesh is probably
  //no longer Delaunay, so we cannot use our favorite Watson insertion scheme.
  //We will insert and swap (locally) to obtain a (locally) Delaunay mesh.

  assert(vertex->getParentTopology());

  double new_coords[] =
    { vertex->x(), vertex->y(), vertex->z() };

  int cell_idx = 0;
  Cell *cell = m_mesh_init->getCell(cell_idx);

  while (!cell->doFullCheck())
    cell = m_mesh_init->getCell(++cell_idx);

  // Any live tet will do, since the domain is concave.
  Vert* newVert = m_pSI3D->insertPoint(new_coords, cell, vertex);
  assert(newVert == vertex);
  return newVert->isValid();

}

bool
TetMeshBuilder::is_bdry_constrained(set<RecoverEdge>& bdry_edges_to_recover,
				    set<RecoverEdge>& face_edges_to_recover,
				    set<RecoverFace>& faces_to_recover) const
{

  //Clear all containers
  bdry_edges_to_recover.clear();
  face_edges_to_recover.clear();
  faces_to_recover.clear();

  //Fill these sets based on current mesh configuration:
  std::set<RecoverEdge> tet_mesh_edges, surf_mesh_bdry_edges,
      surf_mesh_face_edges;
  std::set<RecoverFace> tet_mesh_faces, surf_mesh_faces;

  int num_cells = m_mesh_init->getNumCells();

  //Start by filling the two tetrahedral mesh containers:
  for (int i = 0; i < num_cells; i++) {

    Cell* cell = m_mesh_init->getCell(i);
    if (cell->isDeleted())
      continue;

    assert(cell->isValid());
    assert(cell->getNumFaces() == 4);
    assert(cell->getNumVerts() == 4);

    tet_mesh_edges.insert(RecoverEdge(cell->getVert(0), cell->getVert(1)));
    tet_mesh_edges.insert(RecoverEdge(cell->getVert(0), cell->getVert(2)));
    tet_mesh_edges.insert(RecoverEdge(cell->getVert(0), cell->getVert(3)));
    tet_mesh_edges.insert(RecoverEdge(cell->getVert(1), cell->getVert(2)));
    tet_mesh_edges.insert(RecoverEdge(cell->getVert(1), cell->getVert(3)));
    tet_mesh_edges.insert(RecoverEdge(cell->getVert(2), cell->getVert(3)));

    for (int j = 0; j < 4; j++) {

      Face* face = cell->getFace(j);
      assert(face->isValid());
      assert(face->getNumVerts() == 3);

      tet_mesh_faces.insert(
	  RecoverFace(face->getVert(0), face->getVert(1), face->getVert(2)));

    }

  }

  //Now do the same with the surface mesh:
  typedef vector<Subseg*> SubsegVector;
  SubsegVector subsegs;
//	typedef vector< pair<FaceQueueEntry*, RefFace*> > FaceVec;
  FaceVec rest_faces;

  m_surf_mesh->get_subsegs(subsegs);
  m_surf_mesh->get_restricted_faces(rest_faces);

  SubsegVector::iterator its = subsegs.begin(), its_end = subsegs.end();
  FaceVec::iterator itf = rest_faces.begin(), itf_end = rest_faces.end();

  for (; its != its_end; ++its) {

    Subseg* sub = *its;
    assert(!sub->is_deleted());
    Vert *v0 = sub->get_beg_vert(), *v1 = sub->get_end_vert();

    assert(
	m_subseg_vert_dict.count(v0) == 1 && m_subseg_vert_dict.count(v1) == 1);

    RecoverEdge recover_edge(m_subseg_vert_dict.find(v0)->second,
			     m_subseg_vert_dict.find(v1)->second);
    //Keeping track of which subseg is represented by recover_edge for subsequent splits.
    recover_edge.set_subseg(sub);

    surf_mesh_bdry_edges.insert(recover_edge);

  }

  for (; itf != itf_end; ++itf) {

    FaceQueueEntry* fqe = itf->first;
    RefFace* ref_face = itf->second;

    assert(!fqe->isDeleted());

    assert(fqe->getNumVerts() == 3);

    Vert *v0 = fqe->getVert(0), *v1 = fqe->getVert(1), *v2 = fqe->getVert(2);

    assert(
	m_sample_vert_dict.count(v0) == 1 && m_sample_vert_dict.count(v1) == 1
	    && m_sample_vert_dict.count(v2) == 1);

#ifndef NDEBUG
    if (v0->getVertType() == Vert::eBdryApex)
      assert(dynamic_cast<RefVertex*>(v0->getParentTopology()));
    else if (v0->getVertType() == Vert::eBdryCurve)
      assert(dynamic_cast<RefEdge*>(v0->getParentTopology()));
    else if (v0->getVertType() == Vert::ePseudoSurface)
      assert(dynamic_cast<RefFace*>(v0->getParentTopology()));
    else
      assert(0);

    if (v1->getVertType() == Vert::eBdryApex)
      assert(dynamic_cast<RefVertex*>(v1->getParentTopology()));
    else if (v1->getVertType() == Vert::eBdryCurve)
      assert(dynamic_cast<RefEdge*>(v1->getParentTopology()));
    else if (v1->getVertType() == Vert::ePseudoSurface)
      assert(dynamic_cast<RefFace*>(v1->getParentTopology()));
    else
      assert(0);

    if (v2->getVertType() == Vert::eBdryApex)
      assert(dynamic_cast<RefVertex*>(v2->getParentTopology()));
    else if (v2->getVertType() == Vert::eBdryCurve)
      assert(dynamic_cast<RefEdge*>(v2->getParentTopology()));
    else if (v2->getVertType() == Vert::ePseudoSurface)
      assert(dynamic_cast<RefFace*>(v2->getParentTopology()));
    else
      assert(0);
#endif

    TriFace* triface = dynamic_cast<TriFace*>(fqe->getFace());
    assert(triface);

    RecoverEdge edge1(m_sample_vert_dict.find(v0)->second,
		      m_sample_vert_dict.find(v1)->second);
    RecoverEdge edge2(m_sample_vert_dict.find(v0)->second,
		      m_sample_vert_dict.find(v2)->second);
    RecoverEdge edge3(m_sample_vert_dict.find(v1)->second,
		      m_sample_vert_dict.find(v2)->second);

    if (surf_mesh_bdry_edges.find(edge1) == surf_mesh_bdry_edges.end()) {
      edge1.set_face(triface);
      edge1.set_surface(ref_face);
      surf_mesh_face_edges.insert(edge1);
    }
    if (surf_mesh_bdry_edges.find(edge2) == surf_mesh_bdry_edges.end()) {
      edge2.set_face(triface);
      edge2.set_surface(ref_face);
      surf_mesh_face_edges.insert(edge2);
    }
    if (surf_mesh_bdry_edges.find(edge3) == surf_mesh_bdry_edges.end()) {
      edge3.set_face(triface);
      edge3.set_surface(ref_face);
      surf_mesh_face_edges.insert(edge3);
    }

    RecoverFace face_to_recover(m_sample_vert_dict.find(v0)->second,
				m_sample_vert_dict.find(v1)->second,
				m_sample_vert_dict.find(v2)->second);

    face_to_recover.set_surface(ref_face);
    face_to_recover.set_face(triface);

    surf_mesh_faces.insert(face_to_recover);
    delete fqe;
  }

  std::set_difference(
      surf_mesh_bdry_edges.begin(), surf_mesh_bdry_edges.end(),
      tet_mesh_edges.begin(), tet_mesh_edges.end(),
      std::inserter(bdry_edges_to_recover, bdry_edges_to_recover.begin()));

  std::set_difference(
      surf_mesh_face_edges.begin(), surf_mesh_face_edges.end(),
      tet_mesh_edges.begin(), tet_mesh_edges.end(),
      std::inserter(face_edges_to_recover, face_edges_to_recover.begin()));

  std::set_difference(
      surf_mesh_faces.begin(), surf_mesh_faces.end(), tet_mesh_faces.begin(),
      tet_mesh_faces.end(),
      std::inserter(faces_to_recover, faces_to_recover.begin()));

  logMessage(1, "Number of bdry edges = %lu\n", surf_mesh_bdry_edges.size());
  logMessage(1, "Number of face edges = %lu\n", surf_mesh_face_edges.size());

  //Issue a report on the results
  if (bdry_edges_to_recover.empty() && face_edges_to_recover.empty()
      && faces_to_recover.empty()) {
    logMessage(1, "All edges and faces are properly constrained.\n");
    return true;
  }
  else {
    logMessage(
	1, "Missing %lu input subcurves, %lu/%lu total edges, %lu/%lu faces.\n",
	bdry_edges_to_recover.size(),
	bdry_edges_to_recover.size() + face_edges_to_recover.size(),
	surf_mesh_bdry_edges.size() + surf_mesh_face_edges.size(),
	faces_to_recover.size(), rest_faces.size());
    return false;
  }

}

bool
TetMeshBuilder::clean_exterior()
{

  bool success = true;
  //Must start with a purge because new vertices might
  //have been inserted and faces flipped during boundary
  //recovery. This leads to deleted faces that might
  //correspond to a restricted face (basically, there
  //are more than one copy of the face, one is valid, the
  //others are deleted).

  m_mesh_init->purgeAllEntities();
  assert(m_mesh_init->isValid());
  logMessage(2, "Passed validity check.\n");

  //Mark all the cells with an invalid region tag.
  for (int i = m_mesh_init->getNumCells() - 1; i >= 0; i--) {
    m_mesh_init->getCell(i)->setRegion(iInvalidRegion);
#ifndef NDEBUG
    if (!m_mesh_init->getCell(i)->isDeleted())
      assert(m_mesh_init->getCell(i)->doFullCheck());
#endif
  }

  //The first 8 vertices are the ones connected to the
  //outside bounding box. They need to be deleted.
  for (int i = 0; i < 8; i++) {
    m_mesh_init->getVert(i)->markAsDeleted();
    m_mesh_init->deleteVert(m_mesh_init->getVert(i));
  }

  //Mark the outside boundary faces for deletion.
  assert(m_mesh_init->getNumBdryFaces() == 12);
  for (int i = m_mesh_init->getNumBdryFaces() - 1; i >= 0; i--) {
    BFace *bface = m_mesh_init->getBFace(i);
    bface->getFace()->setFaceLoc(Face::eBdryFace);
    bface->setRegion(iInvalidRegion);
    m_mesh_init->deleteBFace(bface);
  }

  //Build the set of restricted faces. These faces, which should
  //appear in the initial tetrahedralization, will become the boundary
  //faces of the mesh.

  map<FaceQueueEntry*, RefFace*, FaceQueueEntryComp> face_surface_map;

  {

    FaceQueueEntry* fqe;
    map<RecoverFace, RefFace*> restrict_faces;
    map<RecoverFace, RefFace*>::iterator itrf;

//		typedef vector<pair<FaceQueueEntry*, RefFace*>> FaceVec;
    FaceVec rest_faces;
    m_surf_mesh->get_restricted_faces(rest_faces);
    FaceVec::iterator itf = rest_faces.begin(), itf_end = rest_faces.end();

    for (; itf != itf_end; ++itf) {

      fqe = itf->first;
      assert(!fqe->isDeleted());
      assert(fqe->getNumVerts() == 3);

      Vert *v0 = fqe->getVert(0), *v1 = fqe->getVert(1), *v2 = fqe->getVert(2);
      assert(
	  m_sample_vert_dict.count(v0) == 1 && m_sample_vert_dict.count(v1) == 1
	      && m_sample_vert_dict.count(v2) == 1);

      restrict_faces.insert(
	  std::make_pair(
	      RecoverFace(m_sample_vert_dict.find(v0)->second,
			  m_sample_vert_dict.find(v1)->second,
			  m_sample_vert_dict.find(v2)->second),
	      itf->second));
      delete fqe;
    }
    rest_faces.clear();
    //#ifndef NDEBUG
    GR_index_t num_restrict = 0;
    //#endif
    Face * old_face;
    for (int i = m_mesh_init->getNumFaces() - 1; i >= 0; i--) {

      old_face = m_mesh_init->getFace(i);

      old_face->unsetRestrictedDelaunay();

      RecoverFace rec_face(old_face->getVert(0), old_face->getVert(1),
			   old_face->getVert(2));
      itrf = restrict_faces.find(rec_face);

      if (itrf != restrict_faces.end()) {

	old_face->setRestrictedDelaunay();
	set_surf_compare_tol(itrf->second, 1.e-14);
	face_surface_map.insert(
	    std::make_pair(new FaceQueueEntry(old_face), itrf->second));
	//#ifndef NDEBUG
	++num_restrict;
	//#endif
      }

    }

    assert(num_restrict == restrict_faces.size());
    assert(restrict_faces.size() == face_surface_map.size());

  }

  //Find cells to delete by walking. Do not cross boundary
  //faces and restricted faces. Each visited set will have
  //the cells from an entire region.

  {

    bool found_exterior = false;
    int region_tag = 0, this_region_tag;

    Cell* cell, *next_cell;
    Face* face;

    set<Cell*> all_cells, visited;
    set<Cell*>::iterator itc, itc_end;
    deque<Cell*> to_visit;

#ifndef NDEBUG
    GR_index_t num_valid_cells = 0, num_visited_cells = 0;
#endif

    for (int i = m_mesh_init->getNumCells() - 1; i >= 0; i--) {
      cell = m_mesh_init->getCell(i);
      assert(cell->getRegion() == iInvalidRegion);
      if (cell->isDeleted())
	continue;
      all_cells.insert(cell);
    }

#ifndef NDEBUG
    num_valid_cells = all_cells.size();
#endif

    while (!all_cells.empty()) {

      cell = *(all_cells.begin());
      to_visit.push_back(cell);

      while (!to_visit.empty()) {

	cell = to_visit.front();
	to_visit.pop_front();

	assert(!cell->isDeleted());
	assert(cell->getType() == Cell::eTet);
	assert(cell->getNumFaces() == 4);

	if (!visited.insert(cell).second)
	  continue;
	all_cells.erase(cell);

	for (int i = 0; i < 4; ++i) {

	  face = cell->getFace(i);
	  assert(!face->isDeleted());

	  if (face->isRestrictedDelaunay())
	    continue; //do not cross restrict delaunay faces.
	  next_cell =
	      cell == face->getRightCell() ?
		  face->getLeftCell() : face->getRightCell();

	  if (next_cell->isValid() && next_cell->getType() == Cell::eTet)
	    to_visit.push_back(next_cell);

	}

      }

      if (!found_exterior) {
	BFace* bface = m_mesh_init->getBFace(0);
	assert(bface->isDeleted());
	assert(bface->getNumFaces() == 1);
	cell =
	    bface->getFace(0)->getRightCell() == bface ?
		bface->getFace(0)->getLeftCell() :
		bface->getFace(0)->getRightCell();
	if (visited.count(cell) == 1) {
	  found_exterior = true;
	  this_region_tag = iOutsideRegion;
	}
	else {
	  this_region_tag = ++region_tag;
	}
      }
      else
	this_region_tag = ++region_tag;

      itc_end = visited.end();
      for (itc = visited.begin(); itc != itc_end; ++itc) {
	assert(this_region_tag != iMaxRegion);
	cell = *itc;
	cell->setRegion(this_region_tag);
      }

#ifndef NDEBUG
      num_visited_cells += visited.size();
#endif

      visited.clear();

    }

    if (region_tag == 0) {
      logMessage(1, "Bounding surface is not watertight\n");
      success = false;
      return success;
    }

    assert(num_valid_cells == num_visited_cells);
    assert(found_exterior);

  }

  //Mark faces for deletion. All the faces between two deleted
  //tets must be deleted. Then add the boundary faces where
  //they are needed.

  {

    Face* face;
    Cell* lcell, *rcell, *cell;

    // This next loop appears to be completely unnecessary.
    //     //First pass over the faces. Identify and mark
    //     //those that must be deleted.

    //     for(int i = mesh->iNumFaces() - 1; i >= 0; --i) {

    //       face = mesh->pFFace(i);
    //       if(face->qDeleted()) continue;
    //       lcell = face->pCCellLeft();
    //       rcell = face->pCCellRight();

    //       if(lcell->eType() == Cell::eTet && rcell->eType() == Cell::eTet) {
    // 	if(lcell->iRegion() == iOutsideRegion &&
    // 	   rcell->iRegion() == iOutsideRegion) face->vMarkDeleted();
    //       }
    //       else face->vMarkDeleted();

    //     }

    //Second pass over the faces. Those marked restricted
    //will become boundary faces (either interior or exterior).
    //As new faces get created in the loop, make sure we do not
    //over the original number of faces.

    int num_faces = m_mesh_init->getNumFaces();
#ifndef NDEBUG
    GR_index_t num_restrict_faces = 0;
#endif

    for (int i = 0; i < num_faces; ++i) {

      face = m_mesh_init->getFace(i);

      if (face->isDeleted())
	continue;

      //If the face is tagged restricted,
      //create the boundary face for it.
      if (face->isRestrictedDelaunay()) {

	assert(!face->isDeleted());

	face->unsetRestrictedDelaunay();
#ifndef NDEBUG
	FaceQueueEntry * fqe_debug = new FaceQueueEntry(face);
	assert(face_surface_map.count(fqe_debug) == 1);
	delete fqe_debug;
	++num_restrict_faces;
#endif
	lcell = face->getLeftCell();
	rcell = face->getRightCell();
	// It's possible that one of the cells has been deleted.  That's
	// okay; that side is the exterior.
	assert(rcell->isValid() || lcell->isValid());
	assert(
	    (!rcell->isValid() || rcell->getType() == Cell::eTet)
		&& (!lcell->isValid() || lcell->getType() == Cell::eTet));

	//Internal boundary.  Need to duplicate the face and add an
	//internal bdry face between the two.
	if (lcell->isValid() && (lcell->getRegion() != iOutsideRegion)
	    && rcell->isValid() && (rcell->getRegion() != iOutsideRegion)) {

	  //Create a new internal tri boundary face

	  // FIX ME  CFOG, 19 Oct, 2010
	  // This is a bit of a hack:  create a fake internal BFace to
	  // attach the patch to, then clone this face and stick an internal
	  // bdry face in between.  Really should be able to create a bface,
	  // then classify it onto a piece of the geometry; or pass the
	  // geometry info to the creation routine directly instead of
	  // attached to a bdry face.

	  BFace *bfaceTmp = m_mesh_init->getNewBFace(3, true);
	  bfaceTmp->setTopologicalParent(
	      face_surface_map.find(new FaceQueueEntry(face))->second);

	  BFace *bface = m_mesh_init->createBFace(face, bfaceTmp);
	  m_mesh_init->deleteBFace(bfaceTmp);

	  //Set the vertex type correctly for surface verts
	  //(should now be ePseudoSurface, we want eBdryTwoSide)
	  assert(bface->getNumVerts() == 3);
	  for (int ii = 0; ii < 3; ++ii) {
	    assert(
		bface->getVert(ii)->getVertType() == Vert::eBdryApex
		    || bface->getVert(ii)->getVertType() == Vert::eBdryCurve
		    || bface->getVert(ii)->getVertType() == Vert::eBdryTwoSide
		    || //this case if already changed
		    bface->getVert(ii)->getVertType() == Vert::ePseudoSurface);
	    if (bface->getVert(ii)->getVertType() == Vert::ePseudoSurface)
	      bface->getVert(ii)->setType(Vert::eBdryTwoSide);
	  }

	  assert(lcell->doFullCheck() == 1);
	  assert(rcell->doFullCheck() == 1);
	  assert(bface->doFullCheck() == 1);
	  assert(face->doFullCheck());
	}

	//A regular boundary face.
	else {
	  bool qLeftIsOutside =
	      rcell->isValid() && rcell->getRegion() != iOutsideRegion;
	  assert(
	      qLeftIsOutside || (lcell->isValid() && lcell->getRegion() != iOutsideRegion));

	  cell = (qLeftIsOutside) ? lcell : rcell;
	  if (cell)
	    m_mesh_init->deleteCell(cell);

	  BFace* bface = m_mesh_init->createBFace(face);
	  FaceQueueEntry * fqe = new FaceQueueEntry(face);
	  bface->setTopologicalParent(face_surface_map.find(fqe)->second);
	  delete fqe;
	  assert(bface->getType() == Cell::eTriBFace);

	  //Set the vertex types correctly for surface verts
	  //(should now be ePseudoSurface, we want eBdry)
	  assert(bface->getNumVerts() == 3);
	  for (int ii = 0; ii < 3; ++ii) {
	    assert(
		bface->getVert(ii)->getVertType() == Vert::eBdryApex
		    || bface->getVert(ii)->getVertType() == Vert::eBdryCurve
		    || bface->getVert(ii)->getVertType() == Vert::eBdry || //this case if already changed
		    bface->getVert(ii)->getVertType() == Vert::ePseudoSurface);
	    if (bface->getVert(ii)->getVertType() == Vert::ePseudoSurface)
	      bface->getVert(ii)->setType(Vert::eBdry);
	  }

	}

      }

    }

#ifndef NDEBUG
    //Making sure I created the same number of boundary faces
    //as there were restricted faces.
    GR_index_t count = 0;
    for (int i = m_mesh_init->getNumTotalBdryFaces() - 1; i >= 0; --i) {
      BFace* bface = m_mesh_init->getBFace(i);
      if (!bface->isDeleted()) {
	++count;
	assert(bface->getTopologicalParent());
      }
    }
    assert(count == num_restrict_faces);
#endif

  }

  //Some post-verification to make sure everything is correct.

  for (int i = m_mesh_init->getNumCells() - 1; i >= 0; --i) {
    Cell* cell = m_mesh_init->getCell(i);
    if (cell->getRegion() == iOutsideRegion)
      m_mesh_init->deleteCell(cell);
    /*			if(cell->isDeleted()) continue;
     assert(cell->getType() == Cell::eTet);
     assert(cell->doFullCheck() == 1);
     for(int j = 0; j < 4; ++j) assert(!cell->getFace(j)->isDeleted());*/
  }
#ifndef NDEBUG
  for (int i = m_mesh_init->getNumFaces() - 1; i >= 0; --i) {

    Face* face = m_mesh_init->getFace(i);
    if (face->isDeleted())
      continue;

    face->unsetRestrictedDelaunay();
    assert(!face->getLeftCell()->isDeleted());
    assert(!face->getRightCell()->isDeleted());

    if (face->getNumCells() == 1) {
      assert(
	  face->getFaceLoc() == Face::eBdryFace
	      || face->getFaceLoc() == Face::eBdryTwoSide);
      if (face->getLeftCell()->getType() == Cell::eTet) {
	assert(
	    face->getRightCell()->getType() == Cell::eTriBFace
		|| face->getRightCell()->getType() == Cell::eIntTriBFace);
      }
      else {
	assert(face->getRightCell()->getType() == Cell::eTet);
	assert(
	    face->getLeftCell()->getType() == Cell::eTriBFace
		|| face->getLeftCell()->getType() == Cell::eIntTriBFace);
      }
    }
    else {
      assert(face->getNumCells() == 2);
      assert(face->getFaceLoc() == Face::eInterior);
      assert(
	  face->getLeftCell()->getType() == Cell::eTet
	      && face->getRightCell()->getType() == Cell::eTet);
    }

    assert(face->doFullCheck() == 1);

  }

#endif
  for (map<FaceQueueEntry*, RefFace*, FaceQueueEntryComp>::iterator it =
      face_surface_map.begin(); it != face_surface_map.end(); ++it)
    delete it->first;
  face_surface_map.clear();
  //Finally, let's purge the deleted entities

  map<Vert*, Vert*> vert_map;
  m_mesh_init->purgeAllEntities(&vert_map);
  assert(m_mesh_init->isValid());

  update_vert_dict(vert_map);
  return success;
}

void
TetMeshBuilder::update_vert_dict(const std::map<Vert*, Vert*>& verts_old_to_new)
{

  Vert *old_vert;
  VertexDict::iterator it, it_end;

  for (int i = 0; i < 2; ++i) {

    switch (i)
      {
      case 0:
	it = m_subseg_vert_dict.begin(), it_end = m_subseg_vert_dict.end();
	break;
      case 1:
	it = m_sample_vert_dict.begin(), it_end = m_sample_vert_dict.end();
	break;
      default:
	assert(0);
	break;
      }

    for (; it != it_end; ++it) {

      old_vert = it->second;
      assert(verts_old_to_new.count(old_vert) == 1);

      it->second = verts_old_to_new.find(old_vert)->second;
      assert(!it->second->isDeleted());

    }

  }

}

void
TetMeshBuilder::get_subseg_copy(std::set<Subseg*>* subseg_copy) const
{

  if (!subseg_copy)
    return;

  assert(m_surf_mesh);

  std::vector<Subseg*> subsegs;
  m_surf_mesh->get_subsegs(subsegs);

  Subseg* my_subseg;
  Vert *v0, *v1;

  std::vector<Subseg*>::iterator it = subsegs.begin(), it_end = subsegs.end();

  for (; it != it_end; ++it) {

    my_subseg = *it;

    //The vertex of my_subseg are not in the tet mesh, but those in SubsegManager.
    //The correspondance is stored in m_subseg_vert_dict.
    //Find the replacements.
    v0 = my_subseg->get_beg_vert();
    v1 = my_subseg->get_end_vert();

    assert(v0 != v1);
    assert(
	m_subseg_vert_dict.count(v0) == 1 && m_subseg_vert_dict.count(v1) == 1);

    v0 = m_subseg_vert_dict.find(v0)->second;
    v1 = m_subseg_vert_dict.find(v1)->second;

    //Now, create the subseg copy and store it in the set.
    //We do not care about the SubsegBridges... they are
    //no longer useful for boundary-recovered tet meshes.
    subseg_copy->insert(
	new Subseg(v0, v1, my_subseg->get_beg_param(),
		   my_subseg->get_end_param(), my_subseg->get_ref_edge()));

  }

}

void
TetMeshBuilder::ComputeBoundingBox::operator()(const Vert* const vertex)
{

  m_mini->x(std::min(m_mini->x(), vertex->x()));
  m_mini->y(std::min(m_mini->y(), vertex->y()));
  m_mini->z(std::min(m_mini->z(), vertex->z()));

  m_maxi->x(std::max(m_maxi->x(), vertex->x()));
  m_maxi->y(std::max(m_maxi->y(), vertex->y()));
  m_maxi->z(std::max(m_maxi->z(), vertex->z()));

}

CubitVector
TetMeshBuilder::surf_normal_at_point(RefFace* const surface,
				     const CubitVector& point,
				     const double tolerance) const
{

  Surface* surf = surface->get_surface_ptr();
  CubitVector dummy, normal;

  FacetSurface* facet_surf = dynamic_cast<FacetSurface*>(surf);
  if (facet_surf) {
    facet_surf->get_eval_tool()->compare_tol(tolerance);
    facet_surf->closest_point(point, &dummy, &normal);
  }
  else {
    vFatalError("Only support for facet modeler is implemented at the moment",
		"TetMeshBuilder::surf_normal_at_point()");
  }

  return normal;

}

void
TetMeshBuilder::set_surf_compare_tol(RefFace* const surface,
				     const double tolerance) const
{

  Surface* surf = surface->get_surface_ptr();
  FacetSurface* facet_surf = dynamic_cast<FacetSurface*>(surf);

  if (facet_surf)
    facet_surf->get_eval_tool()->compare_tol(tolerance);
  else
    vFatalError("Only support for facet modeler is implemented at the moment",
		"TetMeshBuilder::surf_normal_at_point()");

}
Cell *
TetMeshBuilder::getSeedCellFromVerts(const Vert* const vertex,
				     const Vert* const vert_guess1,
				     const Vert* const vert_guess2,
				     const Vert* const vert_guess3)
{

  assert(vertex->isValid());

  set<const Cell*> guesses;
  set<Cell*> tmp_cells;
  set<Vert*> tmp_verts;

  if (vert_guess1) {
    findNeighborhoodInfo(vert_guess1, tmp_cells, tmp_verts);
    std::copy(tmp_cells.begin(), tmp_cells.end(),
	      std::inserter(guesses, guesses.begin()));
  }
  if (vert_guess2) {
    findNeighborhoodInfo(vert_guess2, tmp_cells, tmp_verts);
    std::copy(tmp_cells.begin(), tmp_cells.end(),
	      std::inserter(guesses, guesses.begin()));
  }
  if (vert_guess3) {
    findNeighborhoodInfo(vert_guess3, tmp_cells, tmp_verts);
    std::copy(tmp_cells.begin(), tmp_cells.end(),
	      std::inserter(guesses, guesses.begin()));
  }
  assert(vertex->getSpaceDimen() == 3);

  set<const Cell*>::iterator it = guesses.begin(), it_end = guesses.end();

  Cell* seed = NULL;

  for (; it != it_end; ++it) {

    //Not too elegant I admit:
    Cell* cell = const_cast<Cell*>(*it);
    assert(cell && !cell->isDeleted());

    if (dynamic_cast<TetCell*>(cell)->circumscribesPointAdaptive(vertex)) {
      seed = cell;
      break;
    }

  }

  if (!seed) {
    int num_cells = m_mesh_init->getNumCells();
    bool cis = false;
    for (int i = 0; i < num_cells; i++) {

      seed = m_mesh_init->getCell(i);
      if (!seed || seed->isDeleted())
	continue;
      cis = dynamic_cast<TetCell*>(seed)->circumscribesPointAdaptive(vertex);
      if (cis)
	break;
    }
    if (!cis)
      vFatalError(
	  "The point you are trying to insert is outside the meshed domain.",
	  "TetMeshBuilder::find_seed_cell_naively()");
  }
  assert(seed);
  return seed;

}
void
TetMeshBuilder::print_markings() const
{

  double sum = 0.;
  double region_sum[iOutsideRegion + 1];
  std::fill(region_sum, region_sum + iOutsideRegion + 1, 0.);

  for (GR_index_t i = 0; i < m_mesh_init->getNumCells(); i++) {

    Cell* cell = m_mesh_init->getCell(i);
    if (cell->isDeleted())
      continue;

    double size = cell->calcSize();
    int region = cell->getRegion();
    region_sum[region] += size;

    if (!(region == iOutsideRegion || region == iInvalidRegion))
      sum += size;

    logMessage(
	1,
	"Cell %3u. Verts %3u, %3u, %3u, %3u.  Size: %10.4g (%10.4g) Region %3d (%10.4g).\n",
	i, m_mesh_init->getVertIndex(cell->getVert(0)),
	m_mesh_init->getVertIndex(cell->getVert(1)),
	m_mesh_init->getVertIndex(cell->getVert(2)),
	m_mesh_init->getVertIndex(cell->getVert(3)), size, sum, region,
	region_sum[region]);

  }

}

