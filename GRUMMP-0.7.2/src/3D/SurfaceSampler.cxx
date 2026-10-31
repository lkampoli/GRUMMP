#include "GR_SurfaceSampler.h"
#include "GR_Cell.h"
#include "GR_Face.h"
#include "GR_FacetSurfIntersect.h"
#include "GR_Geometry.h"
#include "GR_Mesh.h"
#include "GR_SamplingQueue.h"
#include "GR_Subseg.h"
#include "GR_SurfIntersect.h"
#include "GR_SurfMeshBuilder.h"
#include "GR_Vec.h"
#include "GR_Vertex.h"
#include "GR_VolMesh.h"

#include "CubitBox.hpp"
#include "CubitVector.hpp"
#include "FacetEvalTool.hpp"
#include "FacetSurface.hpp"
#include "GeometryDefines.h"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVertex.hpp"

#include <map>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

using std::map;
using std::multimap;
using std::multiset;
using std::pair;
using std::queue;
using std::set;
using std::vector;

SurfaceSampler::SurfaceSampler(RefFace* const surface,
			       const SurfMeshBuilder::ModelerType modeler_type) :
    m_modeler_type(modeler_type), m_mesh(NULL), m_WI(NULL), m_surface(surface), m_surf_intersect(
    NULL), m_sample_verts(), m_restricted_faces(), m_restrict_initialized(
	false), m_new_cells(), m_temporary_fqes()
{

  assert(modeler_type != SurfMeshBuilder::UNKNOWN);
  assert(modeler_type == SurfMeshBuilder::FACET); //temporary
  //Initialize the mesh used by the sampler.
  initialize_mesh();

  if (modeler_type == SurfMeshBuilder::FACET) {
    m_surf_intersect = new FacetSurfIntersect(surface);
  }
  else { //temporary... I hope.
    vFatalError("Only facet based modeler supported for now.",
		"SurfaceSampler::SurfaceSampler()");
  }

  m_WI = new GRUMMP::WatsonInserter(m_mesh);
  m_mesh->addObserver(this, GRUMMP::Observable::cellCreated);
}

SurfaceSampler::~SurfaceSampler()
{
  std::set<FaceQueueEntry*>::iterator it = m_temporary_fqes.begin(), it_end =
      m_temporary_fqes.end();
  for (; it != it_end; ++it)
    if (m_restricted_faces.count(*it) == 0)
      delete *it;
  RestrictedFaceDict::iterator itf = m_restricted_faces.begin(), itf_end =
      m_restricted_faces.end();

  for (; itf != itf_end; ++itf) {
    if (itf->first)
      delete itf->first;
    if (itf->second)
      delete itf->second;
  }
  SampleVertDict::iterator itv = m_sample_verts.begin(), itv_end =
      m_sample_verts.end();

  for (; itv != itv_end; ++itv)
    if (itv->second)
      delete itv->second;

  m_restricted_faces.clear();
  m_sample_verts.clear();
  //Delete what has been instantiated here.
  if (m_mesh) {
    m_mesh->removeObserver(this);
    delete m_mesh;
  }

  if (m_surf_intersect)
    delete m_surf_intersect;
  if (m_WI)
    delete m_WI;
}

void
SurfaceSampler::initialize_mesh()
{

  //Computes the bouding box of the surface.
  //Builds a tetrahedral inside a rectangular prism
  //roughly the size of this box.

  assert(!m_mesh);

  CubitBox box = m_surface->bounding_box();

  CubitVector box_dim = box.maximum() - box.minimum();
  assert(
      iFuzzyComp(box_dim.x(), 0.) != -1 && iFuzzyComp(box_dim.y(), 0.) != -1
	  && iFuzzyComp(box_dim.z(), 0.) != -1);

  double large_dim = MAX3(box_dim.x(), box_dim.y(), box_dim.z());

  CubitVector tolerance(large_dim, large_dim, large_dim);

  CubitVector mini = box.minimum() - tolerance;
  CubitVector maxi = box.maximum() + tolerance;
  m_mesh = new VolMesh(CubitBox(mini, maxi));
}

//This is an initialization function. As such, it first deletes all
//existing restricted triangulation data (if any exists). It populates
//the two vectors with the restricted faces and the sample vertices
//currently in the sampler.

void
SurfaceSampler::init_restricted_delaunay(vector<VertWithTag>& sample_verts,
					 vector<RestrictFace>& restricted_faces)
{

  assert(m_mesh);
  assert(!m_restrict_initialized);
  assert(sample_verts.empty());
  assert(restricted_faces.empty());

  //Remove all existing restricted data.
  for (int i = m_mesh->getNumFaces() - 1; i >= 0; i--)
    m_mesh->getFace(i)->unsetRestrictedDelaunay();

  RestrictedFaceDict::iterator itf = m_restricted_faces.begin(), itf_end =
      m_restricted_faces.end();

  for (; itf != itf_end; ++itf) {
    if (itf->first)
      delete itf->first;
    if (itf->second)
      delete itf->second;
  }
  m_restricted_faces.clear();

  //(re)Compute restricted triangulation.
  Face* face;
#ifndef NDEBUG
  int num_bdry_faces = 0;
#endif
  output_restricted_delaunay("rdel_during.vtk");
  logMessage(MSG_MANAGER, "Number of verts: %ul\n", m_mesh->getNumVerts());

  for (int i = m_mesh->getNumFaces() - 1; i >= 0; i--) {

    face = m_mesh->getFace(i);
    if (face->isDeleted())
      continue;
    if (face->getFaceLoc() == Face::eBdryFace) {
#ifndef NDEBUG
      assert(
	  face->getVert(0)->getVertType() == Vert::eBBox
	      && face->getVert(1)->getVertType() == Vert::eBBox
	      && face->getVert(2)->getVertType() == Vert::eBBox);
      ++num_bdry_faces;
#endif
      continue;
    }

    VoronoiEdge edge(face->getLeftCell(), face->getRightCell());
    compute_restricted_data(edge, restricted_faces);
  }

  //We should have twelve boundary faces (TWELVE ONLY)
  //on the boundary of the bounding box.
  assert(num_bdry_faces == 12);

  //The restricted triangulation is now initialized.
  m_restrict_initialized = true;

  //If no restricted faces were found, and if there are no subsegments
  //i.e. a smooth surface without initial restricted faces, then need to insert
  //arbitrary points on the surface to initializing sampling algo.
  if (restricted_faces.empty() && m_sample_verts.empty()) {
    seed_surface(restricted_faces);
  }

  sample_verts.reserve(m_sample_verts.size());
  std::copy(m_sample_verts.begin(), m_sample_verts.end(),
	    std::back_inserter(sample_verts));

}

//void SurfaceSampler::
//split_box_attached(vector<RestrictFace>& restrict_faces) {
//
//	//I don't think this is necessary.
//	//Plus this is more of a guideline implementation, I think
//	//pieces are still missing, so it might not even work.
//	assert(0);
//
//	FaceQueueEntry* face;
//	bool  split_occured = false;
//	vector<RestrictFace> new_restrict_faces;
//	vector<RestrictFace>::const_iterator it, it_end;
//
//	do {
//
//		it     = restrict_faces.begin();
//		it_end = restrict_faces.end();
//
//		for( ; it != it_end; ++it) {
//
//			face = it->first;
//			if(face->isDeleted()) continue;
//
//			if( face->getVert(0)->getVertType() == Vert::eBBox ||
//					face->getVert(1)->getVertType() == Vert::eBBox ||
//					face->getVert(2)->getVertType() == Vert::eBBox ) {
//
//				//Face is connected to the bounding box, needs to be split.
//				split_occured = true;
//				new_restrict_faces.clear();
//
//				assert(!it->second.empty()); //it->second is the set of intersections.
//				TriFace* tri_face = dynamic_cast<TriFace*>(face->getFace()); assert(tri_face);
//				SurfIntersect::Intersection far_intersect;
//				std::for_each(it->second.begin(), it->second.end(),
//						SurfIntersect::FarthestIntersect(tri_face->calcCircumcenter(), &far_intersect));
//
//				CubitVector insert_coord;
//				far_intersect.get_triangle()->
//						project_to_surface( far_intersect.get_intersect_coord(), &insert_coord );
//
//				Vert dummy_vert;
//				dummy_vert.setCoords(3, insert_coord);
//				insert_sample_point(&dummy_vert, face->getVert(0), NULL, NULL,
//						&new_restrict_faces);
//
//			}
//
//		}
//
//		std::copy( new_restrict_faces.begin(), new_restrict_faces.end(),
//				std::back_inserter(restrict_faces) );
//
//	} while(!new_restrict_faces.empty());
//
//	if(!split_occured) return;
//
//	assert(new_restrict_faces.empty());
//
//	it = restrict_faces.begin();
//	it_end = restrict_faces.end();
//
//	for( ; it != it_end; ++it)
//		if(!it->first->isDeleted())
//			new_restrict_faces.push_back( std::make_pair(it->first, it->second) );
//
//	restrict_faces.swap(new_restrict_faces);
//
//}

//To accelerate intersection search, we only intersect the Voronoi edge with
//the facetted model (there should always be a facetted model available in CGM,
//independently of what modeler is used). Once this is done, we eliminate double
//(triple, quadruple, ...) entries (this will occur if an edge intersects a
//triangle on its boundary).

//The function returns all the faces whose Voronoi edge intersects the surface
//once or more, accompanied by at set of these intersections.

void
SurfaceSampler::compute_restricted_data(
    const VoronoiEdge& voronoi_edge, vector<RestrictFace>& new_restrict_faces)
{

  assert(m_surf_intersect);
  assert(!voronoi_edge.get_face()->isDeleted());
  assert(!voronoi_edge.get_face()->getLeftCell()->isDeleted());
  assert(!voronoi_edge.get_face()->getRightCell()->isDeleted());

  //Get face associated with Voronoi edge.
  Face* face = voronoi_edge.get_face();
  RestrictedFaceTag* face_tag = NULL;
  assert(face->getNumVerts() == 3);

  //Finds intersection with Voronoi edge and surface.
  SurfIntersect::IntersectSet intersections;
  m_surf_intersect->intersect_with_voronoi_edge(voronoi_edge, intersections);

  //Sets tags and intersection info depending on face's current status.
  if (face->isRestrictedDelaunay()) {

    //Face was restricted Delaunay.
    FaceQueueEntry * tmp_fqe = new FaceQueueEntry(face);
    RestrictedFaceDict::iterator it_rest_face = m_restricted_faces.find(
	tmp_fqe);

    if (it_rest_face == m_restricted_faces.end()) {
      face_tag = new RestrictedFaceTag(tmp_fqe);
    }
    else {
      face_tag = it_rest_face->second;
      delete tmp_fqe;
    }
    if (!intersections.empty()) { //Face remains restricted Delaunay.
      face_tag->reset_tag();
      new_restrict_faces.push_back(
	  std::make_pair(face_tag->get_fqe(), intersections));
    }

    else { //Face no longer restricted Delaunay.
      for (int i = 0; i < 3; i++) {

	if (m_sample_verts.count(face->getVert(i)) == 1) {
	  assert(m_sample_verts.count(face->getVert(i)) == 1);
	  m_sample_verts.find(face->getVert(i))->second->remove_ring_face(
	      face_tag->get_fqe());
	}
      }
      if (face_tag) {
	delete face_tag;
      }
      if (it_rest_face != m_restricted_faces.end()) {
	m_restricted_faces.erase(it_rest_face);

      }
      face->unsetRestrictedDelaunay();
      return;
    }

  }
  else { //Face was not restricted Delaunay.

    // m_restricted_faces.erase(face); // May be a no-op...

    if (!intersections.empty()) { //Face becomes restricted Delaunay.
      FaceQueueEntry * tmp_fqe = new FaceQueueEntry(face);
      face->setRestrictedDelaunay();
      face_tag = new RestrictedFaceTag(tmp_fqe);
      m_restricted_faces.insert(std::make_pair(tmp_fqe, face_tag));
      new_restrict_faces.push_back(std::make_pair(tmp_fqe, intersections));

      for (int i = 0; i < 3; i++) { //Face does not change status.
	SampleVertDict::iterator it_vert;
	it_vert = m_sample_verts.find(face->getVert(i));
	if (it_vert != m_sample_verts.end())
	  it_vert->second->add_ring_face(tmp_fqe);
      }
      m_temporary_fqes.insert(tmp_fqe);
    }

    else {
      return;
    }

  }

  //Set intersection data.
  //It is possible that many intersections were found for the Voronoi edge.
  //If that is the case, use the farthest from face's circumcenter and set
  //accordingly in the face tag.
  assert(face_tag);
  assert(face_tag->get_face() == face);

  TriFace* tri_face = dynamic_cast<TriFace*>(face);
  assert(tri_face);

  SurfIntersect::Intersection far_intersect;
  std::for_each(
      intersections.begin(),
      intersections.end(),
      SurfIntersect::FarthestIntersect(tri_face->calcCircumcenter(),
				       &far_intersect));

  face_tag->set_intersect(far_intersect);

}

void
SurfaceSampler::seed_surface(vector<RestrictFace>& restrict_faces)
{

  assert(m_surface);
  assert(m_surf_intersect);

  SurfTri* triangle;
  Vert insert_vert, *vertex;
  SampleVertTag* vertex_tag;
  CubitVector insert_coord, insert_normal, triangle_centroid;

  SurfIntersect::SurfTriSet surf_tris;
  m_surf_intersect->get_triangles(surf_tris);

  SurfIntersect::SurfTriSet::iterator it_tri = surf_tris.begin(), it_tri_end =
      surf_tris.end();

  SampleVertDict::iterator itv, itv_end;

  vector<RestrictFace> tmp_vec;

  while (it_tri != it_tri_end) {

    tmp_vec.clear();

    triangle = *it_tri;
    triangle_centroid = triangle->get_centroid();
    triangle->project_to_surface(triangle_centroid, &insert_coord,
				 &insert_normal);

    insert_vert.setCoords(3, insert_coord);
    insert_vert.setParentTopology(m_surface);
    insert_vert.setType(Vert::ePseudoSurface);

    vertex = insert_sample_point(&insert_vert, NULL, NULL, NULL, &tmp_vec);
    vertex_tag = get_vert_tag(vertex);
    vertex_tag->set_normal(insert_normal);
    vertex_tag->set_curvature(triangle->curvature_at_coord(triangle_centroid));

    std::copy(tmp_vec.begin(), tmp_vec.end(),
	      std::back_inserter(restrict_faces));

    itv = m_sample_verts.begin();
    itv_end = m_sample_verts.end();
    for (; itv != itv_end; ++itv)
      if (itv->second->num_ring_faces() == 0) {
	++it_tri;
	break;
      }

    if (itv != itv_end)
      continue;
    if (!restrict_faces.empty())
      break;

  }

  assert(!restrict_faces.empty());

}

Vert*
SurfaceSampler::insert_sample_point(const Vert* const vertex,
				    const Vert* const guess1,
				    const Vert* const guess2,
				    const Vert* const guess3,
				    vector<RestrictFace>* new_restrict_faces)
{

  set<const Vert*> vert_guesses;
  if (guess1)
    vert_guesses.insert(guess1);
  if (guess2)
    vert_guesses.insert(guess2);
  if (guess3)
    vert_guesses.insert(guess3);

  return insert_sample_point(vertex, vert_guesses, new_restrict_faces);

}

Vert*
SurfaceSampler::insert_sample_point(const Vert* const vertex,
				    const Face* const face, int attemptNum,
				    vector<RestrictFace>* new_restrict_faces)
{

  Vert* new_vert = NULL;

  set<Cell*> guesses;
  guesses.insert(face->getLeftCell());
  guesses.insert(face->getRightCell());

  double oldLoc[3];
  CubitVector newPoint;
  // compare vertex to face points;

  switch (attemptNum)
    {
    case 0:
      {
	//		calcCircumcenter3D(face->getVert(0),face->getVert(1),face->getVert(2),oldLoc);
	//		double radius = dDIST3D(oldLoc,face->getVert(0)->getCoords());
	//		for(int iV = 0; iV < 3; iV++){
	//			if(dDIST3D(face->getVert(iV)->getCoords(),vertex->getCoords())/radius < 0.1) //don't insert
	//				return NULL;
	//		}
	new_vert = insert_sample_point(vertex, guesses, true,
				       new_restrict_faces);
	// hacky attempt to fix things up...
	// try other point locations if the other one doesn't work.
	// circumcenter, then centroid, then midpoint of longest edge

      }
      break;
    case 1:
      {

	calcCircumcenter3D(face->getVert(0), face->getVert(1), face->getVert(2),
			   oldLoc);
	m_surface->find_closest_point_trimmed(CubitVector(oldLoc), newPoint);

	const_cast<Vert*>(vertex)->setCoords(3, newPoint);

	double radius = dDIST3D(oldLoc, face->getVert(0)->getCoords());
	for (int iV = 0; iV < 3; iV++) {
	  if (dDIST3D(face->getVert(iV)->getCoords(), vertex->getCoords())
	      / radius < 0.1) //don't insert
	    return NULL;
	}
	new_vert = insert_sample_point(vertex, guesses, true,
				       new_restrict_faces);
      }
      break;

    case 2:
      {
	face->calcCentroid(oldLoc);
	//	calcCircumcenter3D(face->getVert(0),face->getVert(1),face->getVert(2),oldLoc);
	m_surface->find_closest_point_trimmed(CubitVector(oldLoc), newPoint);
	const_cast<Vert*>(vertex)->setCoords(3, newPoint);

	double radius = dDIST3D(oldLoc, face->getVert(0)->getCoords());
	for (int iV = 0; iV < 3; iV++) {
	  if (dDIST3D(face->getVert(iV)->getCoords(), vertex->getCoords())
	      / radius < 0.1) //don't insert
	    return NULL;
	}
	new_vert = insert_sample_point(vertex, guesses, true,
				       new_restrict_faces);
      }
      break;
    case 3:
      {
	double length = 0, newlength;
	int iMax = 0;
	for (int i = 0; i < 3; i++) {
	  newlength = dDist(face->getVert(i), face->getVert((i + 1) % 3));
	  if (newlength > length) {
	    iMax = i;
	    length = newlength;
	  }
	}
	oldLoc[0] = 0.5
	    * (face->getVert(iMax)->x() + face->getVert((iMax + 1) % 3)->x());
	oldLoc[1] = 0.5
	    * (face->getVert(iMax)->y() + face->getVert((iMax + 1) % 3)->y());
	oldLoc[2] = 0.5
	    * (face->getVert(iMax)->z() + face->getVert((iMax + 1) % 3)->z());
	m_surface->find_closest_point_trimmed(CubitVector(oldLoc), newPoint);
	const_cast<Vert*>(vertex)->setCoords(3, newPoint);

	double radius = dDIST3D(oldLoc, face->getVert(0)->getCoords());
	for (int iV = 0; iV < 3; iV++) {
	  if (dDIST3D(face->getVert(iV)->getCoords(), vertex->getCoords())
	      / radius < 0.1) //don't insert
	    return NULL;
	}

	new_vert = insert_sample_point(vertex, guesses, true,
				       new_restrict_faces);
      }
      break;
    case 4:
      {
	std::set<const Vert*> vert_guesses;
	vert_guesses.insert(face->getVert(0));
	vert_guesses.insert(face->getVert(1));
	vert_guesses.insert(face->getVert(2));

	new_vert = insert_sample_point(vertex, vert_guesses,
				       new_restrict_faces);
      }
      break;
    }
  if (new_vert)
    return new_vert;
  //	else if(attemptNum < 5) insert_sample_point(vertex,face,attemptNum++,new_restrict_faces);
  else
    return pVInvalidVert;
}

Vert*
SurfaceSampler::insert_sample_point(const Vert* const vertex,
				    const set<const Vert*>& vert_guesses,
				    vector<RestrictFace>* new_restrict_faces)
{

  set<Cell*> seed_guesses;
  set<Cell*> cells_tmp;
  set<Vert*> verts_tmp;

  set<const Vert*>::iterator vert = vert_guesses.begin(), vert_end =
      vert_guesses.end();

  for (; vert != vert_end; ++vert) {

    Vert* guess_vert = const_cast<Vert*>(*vert);
    assert(m_sample_verts.find(guess_vert) != m_sample_verts.end());
    cells_tmp.clear();
    verts_tmp.clear();
    findNeighborhoodInfo(guess_vert, cells_tmp, verts_tmp);
    std::copy(cells_tmp.begin(), cells_tmp.end(),
	      std::inserter(seed_guesses, seed_guesses.begin()));

  }
  return insert_sample_point(vertex, seed_guesses, false, new_restrict_faces);

}

Vert*
SurfaceSampler::insert_sample_point(const Vert* const vertex,
				    const set<Cell*>& seed_guesses,
				    const bool bail_on_bad_guesses,
				    vector<RestrictFace>* new_restrict_faces)
{

  qAdaptPred = true;
  Vert* new_vert = NULL;

  //Find if there is a valid Watson seed among the guesses.
  Cell* seed = NULL;

  set<Cell*>::reverse_iterator cell = seed_guesses.rbegin(), cell_end =
      seed_guesses.rend();

  for (; cell != cell_end; ++cell) {
    Cell* guess = *cell;
    assert(guess);
    assert(!guess->isDeleted());
    // this is hacky cause of the const everything.
    if (dynamic_cast<TetCell*>(guess)->circumscribesPoint(
	vertex->getCoords())) {
      seed = guess;
      break;
    }
  }

  //Can choose to bail if no valid seed is found in the set.
  //In this case, return a NULL pointer.
  if (bail_on_bad_guesses && !seed)
    return static_cast<Vert*>(NULL);

  //If a seed was found, use it for insertion,
  //If not, we will have to find a seed naively.

  if (seed) {
    assert(seed->getType() == Cell::eTet);
    assert(
	dynamic_cast<TetCell*>(seed)->circumscribesPoint(vertex->getCoords()));

    //		new_vert = m_mesh->insert_vert_in_mesh(vertex, seed, &new_cells);
  }
  else {

    //		if(!seed_guesses.empty()){
    int num_cells = m_mesh->getNumCells();
    bool cis = false;
    for (int i = 0; i < num_cells; i++) {

      seed = m_mesh->getCell(i);
      if (!seed || seed->isDeleted())
	continue;
      cis = dynamic_cast<TetCell*>(seed)->circumscribesPoint(
	  vertex->getCoords());
      if (cis)
	break;
    }

    //		} else
    //			seed = *seed_guesses.begin();
    //		bool status;
    //		seed = m_mesh->findCell(vertex->getCoords(),seed,status);
    if (!cis)
      if (dynamic_cast<TetCell*>(seed)->circumscribesPoint(vertex->getCoords()))
	vFatalError(
	    "The point you are trying to insert is outside the meshed domain.",
	    "SurfaceSampler::find_seed_cell()");

  }
  new_vert = m_mesh->createVert(vertex);
  m_WI->computeHull(new_vert->getCoords(), seed);
  m_WI->insertPointInHull(new_vert);
  m_mesh->sendEvents();
  qAdaptPred = false;
  assert(new_vert->isValid());

  if (m_restrict_initialized) {

    RestrictedFaceDict::iterator itf = m_restricted_faces.begin(), itf_end =
	m_restricted_faces.end();

    while (itf != itf_end) {
      FaceQueueEntry * fqe = itf->first;
      if (fqe->isDeleted()) {
	RestrictedFaceTag* face_tag = itf->second;
	for (int i = 0; i < 3; i++) {
	  if (m_sample_verts.count(fqe->getVert(i)) == 1) {
	    m_sample_verts.find(fqe->getVert(i))->second->remove_ring_face(fqe);
	  }
	}
	if (face_tag) {
	  delete face_tag;
	}
	RestrictedFaceDict::iterator toErase = itf++;
	m_restricted_faces.erase(toErase);
      }
      else {
	++itf;
      }

    }

  }

  assert(new_vert);
  assert(m_sample_verts.count(new_vert) == 0);
  //Create the sample vertex tag and insert it into the dictionary.
  SampleVertTag* vertex_tag = new SampleVertTag(new_vert);
  m_sample_verts.insert(std::make_pair(new_vert, vertex_tag));

  //Update the Voronoi diag and obtain new edges.
  vector<VoronoiEdge> new_voronoi_edges;
  build_voronoi_edges(m_new_cells, new_voronoi_edges);

  //Update restricted triangulation data based
  //on the new Voronoi edges. Gets new restricted faces.

  vector<RestrictFace> rest_faces;

  if (m_restrict_initialized) {
    //print_surface_sampler();
    vector<VoronoiEdge>::iterator edge = new_voronoi_edges.begin(), edge_end =
	new_voronoi_edges.end();
    int iv = 0;
    for (; edge != edge_end; ++edge) {
      compute_restricted_data(*edge, rest_faces);
      iv++;
    }
  }

  if (new_restrict_faces) {
    assert(new_restrict_faces->empty());
    new_restrict_faces->swap(rest_faces);
  }

  return new_vert;

}
Vert*
SurfaceSampler::insert_arbitrary_point(const Vert* const vertex, Cell * guess,
				       vector<RestrictFace>* new_restrict_faces)
{
  bool status;
  if (!guess) {
    for (GR_index_t iC = 0; iC < m_mesh->getNumCells(); iC++) {
      guess = m_mesh->getCell(iC);
      if (guess->isValid() && !guess->isDeleted())
	break;
    }
  }
  guess = m_mesh->findCell(vertex->getCoords(), guess, status);
  set<Cell*> cell_guesses;
  cell_guesses.insert(guess);
  Vert * newVert = insert_sample_point(vertex, cell_guesses,
				       new_restrict_faces);
  assert(newVert->isValid());
  guess = newVert->getFace(0)->getCell(0);
  return newVert;

}
void
SurfaceSampler::get_sample_verts(vector<Vert*>& sample_verts) const
{

  //Note that the first 8 vertices are always
  //the corner of the bounding box. We do not want them.

  Vert* vertex;

  for (int i = m_mesh->getNumVerts() - 1; i >= 8; i--) {

    vertex = m_mesh->getVert(i);
    if (vertex->isDeleted())
      continue;

    assert(m_sample_verts.count(vertex) == 1);
    sample_verts.push_back(vertex);

  }

}

// Note. new face queue entries need to be deleted manually, they are not self cleaning
void
SurfaceSampler::get_restricted_faces(
    vector<FaceQueueEntry*>& restricted_faces) const
{
  Face* face;

  for (int i = m_mesh->getNumFaces() - 1; i >= 0; i--) {

    face = m_mesh->getFace(i);

    if (face->isDeleted())
      continue;

    if (face->isRestrictedDelaunay()) {
      restricted_faces.push_back(new FaceQueueEntry(face));
    }
  }

}
void
SurfaceSampler::reset_restricted_faces() const
{
  RestrictedFaceDict::const_iterator itf = m_restricted_faces.begin();
  RestrictedFaceDict::const_iterator itf_end = m_restricted_faces.end();

  for (; itf != itf_end; ++itf) {

    if (itf->first->getFace()->isRestrictedDelaunay()) {
      itf->first->reset();
      assert(!itf->first->isDeleted());

    }
    else {
      delete itf->first;
    }
  }

}
SampleVertTag*
SurfaceSampler::get_vert_tag(Vert* const vertex) const
{

  assert(m_sample_verts.count(vertex) == 1);
  return m_sample_verts.find(vertex)->second;

}

RestrictedFaceTag*
SurfaceSampler::get_face_tag(FaceQueueEntry* fqe) const
{

#ifndef NDEBUG
  assert(fqe->getFace()->isRestrictedDelaunay());
  assert(m_restricted_faces.count(fqe) == 1);
#endif

  return m_restricted_faces.find(fqe)->second;

}

CubitVector
SurfaceSampler::normal_at(const Vert* const vertex,
			  const double tolerance) const
{

  assert(m_surface);
  assert(m_sample_verts.count(const_cast<Vert*>(vertex)) == 1);

  Surface* surf = m_surface->get_surface_ptr();
  CubitVector normal, dummy;
  double old_tol;
  FacetSurface* facet_surf;

  switch (m_modeler_type)
    {

    case SurfMeshBuilder::FACET:
      facet_surf = dynamic_cast<FacetSurface*>(surf);
      assert(facet_surf);
      old_tol = facet_surf->get_eval_tool()->compare_tol();
      facet_surf->get_eval_tool()->compare_tol(tolerance);
      facet_surf->closest_point(CubitVector(vertex->getCoords()), &dummy,
				&normal);
      facet_surf->get_eval_tool()->compare_tol(old_tol);

      break;

    case SurfMeshBuilder::UNKNOWN:
      vFatalError("Modeler type not set properly.",
		  "SurfaceSampler::normal_at()");
      break;

    case SurfMeshBuilder::ACIS:
    case SurfMeshBuilder::OCC:
    default:
      vFatalError("Only support for facet modeler is currently implemented.",
		  "SurfaceSampler::normal_at()");
      break;

    }

  return normal;

}

void
SurfaceSampler::project_to_surface(const CubitVector& point,
				   CubitVector& point_on_surf,
				   CubitVector& normal_at_point,
				   const double tolerance) const
{

  Surface* surf = m_surface->get_surface_ptr();

  FacetSurface* facet_surf;

  switch (m_modeler_type)
    {

    case SurfMeshBuilder::FACET:
      facet_surf = dynamic_cast<FacetSurface*>(surf);
      assert(facet_surf);
      facet_surf->get_eval_tool()->compare_tol(tolerance);
      facet_surf->closest_point(point, &point_on_surf, &normal_at_point);
      break;

    case SurfMeshBuilder::UNKNOWN:
      vFatalError("Modeler type not set properly.",
		  "SurfaceSampler::normal_at()");
      break;

    case SurfMeshBuilder::ACIS:
    case SurfMeshBuilder::OCC:
    default:
      vFatalError("Only support for facet modeler is currently implemented.",
		  "SurfaceSampler::project_to_surface()");
      break;

    }

}

//void SurfaceSampler::
//restricted_neighborhood(const Vert* const vertex,
//		vector<Face*>* const neigh_faces,
//		vector<Vert*>* const neigh_verts) const {
//
//	//computes the 1-star of vertex in the restricted triangulation.
//	//Not really used anymore since the star is stored in the vertex tag.
//
//	assert(neigh_faces);
//	neigh_faces->clear();
//
//	int i;
//
//	Face* face = vertex->getFace(0);
//	assert(face->isValid() && face->hasVert(vertex));
//
//	Cell* cell = face->getLeftCell();
//	if(cell->isBdryCell()) cell = face->getRightCell();
//	assert(!cell->isBdryCell());
//	assert(cell->getNumFaces() == 4);
//
//	queue< std::pair<Face*, Cell*> > face_to_traverse;
//	set<Face*>                       face_traversed;
//
//	for(i = 0; i < 4; i++) {
//		face = cell->getFace(i);
//		if(face->hasVert(vertex) && face_traversed.insert(face).second)
//			face_to_traverse.push( std::make_pair(face, face->getOppositeCell(cell)) );
//	}
//
//	while(!face_to_traverse.empty()) {
//
//		face = face_to_traverse.front().first;
//		cell = face_to_traverse.front().second;
//		assert(face->hasVert(vertex));
//		assert(cell->hasFace(face));
//		assert(!cell->isBdryCell());
//		assert(cell->getNumFaces() == 4);
//
//		face_to_traverse.pop();
//
//		if(face->isRestrictedDelaunay())
//			neigh_faces->push_back(face);
//
//		for(i = 0; i < 4; i++) {
//			face = cell->getFace(i);
//			if(face->hasVert(vertex) && face_traversed.insert(face).second)
//				face_to_traverse.push( std::make_pair(face, face->getOppositeCell(cell)) );
//		}
//
//	}
//
//	if(neigh_verts) {
//
//		Vert* vert;
//		set<Vert*> dummy;
//
//		vector<Face*>::iterator
//		it     = neigh_faces->begin(),
//		it_end = neigh_faces->end();
//
//		for( ; it != it_end; ++it) {
//
//			face = *it;
//			assert(face->getNumVerts() == 3);
//
//			for(i = 0; i < 3; i++) {
//				vert = face->getVert(i);
//				if(vert != vertex && dummy.insert(vert).second)
//					neigh_verts->push_back(vert);
//			}
//
//		}
//
//	}
//
//}

void
SurfaceSampler::build_voronoi_edges(const vector<Cell*>& new_cells,
				    vector<VoronoiEdge>& new_voronoi_edges)
{

  int i;
  Face* face;
  Cell *cell, *cell1, *cell2;
  CubitVector dummy;

  set<Face*> faces_visited;

  typedef map<Cell*, CubitVector> VisitedCells;
  VisitedCells cells_visited;

  typedef std::pair<VisitedCells::iterator, bool> InsertData;
  InsertData ID1, ID2;

  vector<Cell*>::const_iterator it_cell = new_cells.begin(), it_cell_end =
      new_cells.end();

  for (; it_cell != it_cell_end; ++it_cell) {

    cell = *it_cell;
    assert(!cell->isDeleted());

    for (i = cell->getNumFaces() - 1; i >= 0; i--) {

      face = cell->getFace(i);

      if (faces_visited.insert(face).second) {

	cell1 = face->getLeftCell();
	cell2 = face->getRightCell();

	ID1 = cells_visited.insert(std::make_pair(cell1, dummy));
	ID2 = cells_visited.insert(std::make_pair(cell2, dummy));

	if (ID1.second)
	  VoronoiEdge::voronoi_vert_coord(cell1, ID1.first->second);
	if (ID2.second)
	  VoronoiEdge::voronoi_vert_coord(cell2, ID2.first->second);
	//	if((ID1.first->second - ID2.first->second).length() > 1e-10)
	new_voronoi_edges.push_back(
	    VoronoiEdge(face, ID1.first->second, ID2.first->second));

      }

    }

  }

}

void
SurfaceSampler::attach_bridge(Vert* const vertex, SubsegBridge* const bridge)
{

  assert(m_sample_verts.count(vertex) == 1);
  m_sample_verts.find(vertex)->second->add_bridge(bridge);

}

void
SurfaceSampler::remove_bridge(Vert* const vertex, SubsegBridge* const bridge)
{

  assert(m_sample_verts.count(vertex) == 1);
  m_sample_verts.find(vertex)->second->remove_bridge(bridge);

}

void
SurfaceSampler::clean_vert_tags()
{

  //Removes all tags whose Vert is deleted.
  //Also removes the deleted faces from the 1-star stored in the tag.

  Vert* vertex;
  SampleVertTag* vertex_tag;
  SampleVertDict dummy_map;

  SampleVertDict::iterator itv = m_sample_verts.begin(), itv_end =
      m_sample_verts.end();

  for (; itv != itv_end; ++itv) {

    vertex = itv->first;
    vertex_tag = itv->second;

    if (vertex->isDeleted()) {
      if (vertex_tag)
	delete vertex_tag;
    }
    else {
      vertex_tag->clean_tag();
      dummy_map.insert(std::make_pair(vertex, vertex_tag));
    }

  }

  m_sample_verts.swap(dummy_map);

#ifndef NDEBUG

  Face* face;

  itv = m_sample_verts.begin();
  itv_end = m_sample_verts.end();

  for (; itv != itv_end; ++itv) {

    vertex = itv->first;
    vertex_tag = itv->second;
    assert(vertex == vertex_tag->get_vertex());

    FaceQueueEntrySet::iterator itf = vertex_tag->get_ring_faces().begin(),
	itf_end = vertex_tag->get_ring_faces().end();

    for (; itf != itf_end; ++itf) {
      FaceQueueEntry * fqe = *itf;
      face = fqe->getFace();
      if (!face->hasVert(vertex) || !face->isRestrictedDelaunay()
	  || m_restricted_faces.count(fqe) == 0) {
	vertex_tag->remove_ring_face(fqe);
	delete fqe;

      }
      //
      //			assert(face && face->isValid());
      //			assert(!face->isDeleted());
      //			assert(face->hasVert(vertex));
      //			assert(face->isRestrictedDelaunay());
      //			assert(m_restricted_faces.count(face) == 1);
    }

  }

#endif

}

void
SurfaceSampler::clean_face_tags()
{

  //Removes all tags whose Face is deleted.

  FaceQueueEntry* fqe;
  RestrictedFaceTag* face_tag;
  RestrictedFaceDict dummy_map;

  RestrictedFaceDict::iterator it = m_restricted_faces.begin(), it_end =
      m_restricted_faces.end();

  for (; it != it_end; ++it) {

    fqe = it->first;
    face_tag = it->second;

    if (fqe->isDeleted()) {
      if (face_tag)
	delete face_tag;
      delete fqe;
    }
    else {
      face_tag->clean_tag();
      dummy_map.insert(std::make_pair(fqe, face_tag));
    }

  }

  m_restricted_faces.swap(dummy_map);

}

void
SurfaceSampler::purge_sampler()
{
  clean_vert_tags();
  clean_face_tags();
}
//
//void SurfaceSampler::
//purge_sampler(VertMap& vert_map, FaceMap& face_map,
//		CellMap& cell_map) {
////
////	//Removes deleted entities from the dictionaries
////	//as well as deleted entities from tags.
//	clean_vert_tags();
//	clean_face_tags();
////
////#ifndef NDEBUG
////
////	//Making sure tag cleaning went well.
////
////	int num_faces = m_mesh->get_mesh()->getNumFaces(), num_undeleted = 0;
////
////	for(int i = 0; i < num_faces; i++) {
////
////		Face* face = m_mesh->get_mesh()->getFace(i);
////		assert(face->isValid());
////
////		if( face->isDeleted() ) {
////			assert( m_restricted_faces.count(face) == 0 );
////			continue;
////		}
////
////		++num_undeleted;
////		Vert* vert0 = face->getVert(0), *vert1 = face->getVert(1), *vert2 = face->getVert(2);
////		SampleVertTag *tag0, *tag1, *tag2;
////
////		if( vert0->getVertType() == Vert::eBBox ||
////				vert1->getVertType() == Vert::eBBox ||
////				vert2->getVertType() == Vert::eBBox ) {
////			assert( !face->isRestrictedDelaunay() );
////			assert( m_restricted_faces.count(face) == 0 );
////			continue;
////		}
////		else {
////			tag0 = get_vert_tag(vert0);
////			tag1 = get_vert_tag(vert1);
////			tag2 = get_vert_tag(vert2);
////		}
////
////		if(face->isRestrictedDelaunay()) {
////			assert( m_restricted_faces.count(face) == 1 );
////			assert( tag0->has_ring_face(face) );
////			assert( tag1->has_ring_face(face) );
////			assert( tag2->has_ring_face(face) );
////		}
////		else {
////			assert( m_restricted_faces.count(face) == 0 );
////			assert( !tag0->has_ring_face(face) );
////			assert( !tag1->has_ring_face(face) );
////			assert( !tag2->has_ring_face(face) );
////		}
////	}
////
////#endif
////
////	//Make sure the containers are empty then,
////	//PURGE deleted entities in the Tet mesh, mapping old pointers to new pointers.
////	vert_map.clear(); face_map.clear(); cell_map.clear();
////	m_mesh->get_mesh()->purgeAllEntities(&vert_map, &face_map, &cell_map);
////
////	assert( m_mesh->get_mesh()->getNumVerts() == vert_map.size() );
////	assert( m_mesh->get_mesh()->getNumFaces() == face_map.size() );
////	assert( m_mesh->get_mesh()->getNumCells() == cell_map.size() );
////
//
////
////#ifndef NDEBUG
////
////	for(itv = vert_map.begin(); itv != itv_end; ++itv) {
////
////		Vert* this_vert = itv->first;
////		int vert_count = m_sample_verts.count(this_vert), vert_index;
////
////		switch(vert_count) {
////		case 0:
////			//The first 8 vertices should be the corner of the bounding
////			//box and as such, do not appear in the vertex dictionary.
////			assert(this_vert->getVertType() == Vert::eBBox);
////			vert_index = m_mesh->get_mesh()->getVertIndex(this_vert);
////			assert( vert_index >= 0 && vert_index < 8 );
////			break;
////		case 1:
////			break;
////		default:
////			assert(0); break;
////		}
////
////	}
////
////#endif
////
//	//Based on the maps, UPDATE THE DATA STRUCTURES.
//
//	//Update the sample vert dictionary and the associated tags.
//
////	Vert* old_vert, *new_vert;
////	Face* old_face, *new_face;
////	SampleVertTag*     vertex_tag;
////	RestrictedFaceTag* face_tag;
////
////	SampleVertDict::iterator     itv, itvdict;
////	RestrictedFaceDict::iterator itfdict;
////
//////	set<Face*> new_ring_faces;
////	set<Face*>::iterator itface, itface_end;
////	set<SubsegBridge*>::iterator itbridge, itbridge_end;
//
//	//Start by updating the vertices and their tags:
//	//the ring faces in the tag, then the bridges.
//
////	for(itv = m_sample_verts.begin(); itv != m_sample_verts.end(); ++itv) {
////
////		old_vert = itv->first;
////		if(old_vert->getVertType() == Vert::eBBox) continue;
////
////		new_vert = itv->second;
////
////		itvdict = m_sample_verts.find(old_vert);
////		assert(itvdict != m_sample_verts.end());
////		assert(itvdict->first == old_vert);
////
////		vertex_tag = itvdict->second;
////		assert(vertex_tag->get_vertex() == old_vert);
////		vertex_tag->set_vertex(new_vert);
////
////		//Replace the ring faces with their new post-purge address.
////		itface     = vertex_tag->get_ring_faces().begin();
////		itface_end = vertex_tag->get_ring_faces().end();
////		new_ring_faces.clear();
////		printf("vt %d %d\n",vertex_tag->get_ring_faces().size(), new_ring_faces.size());
////		for( ; itface != itface_end; ++itface) {
////
////			itf = face_map.find( *itface );
////			assert( itf != face_map.end() );
////			new_face = itf->second;
////			assert( new_face->hasVert(new_vert) );
////			assert( new_face->isRestrictedDelaunay() );
////			new_ring_faces.insert(itf->second);
////
////		}
////
////		assert( vertex_tag->get_ring_faces().size() == new_ring_faces.size() );
////		vertex_tag->replace_all_ring_faces(new_ring_faces);
////
////		//If the vertex address did not change, nothing more to do.
////		if( old_vert == new_vert ) continue;
////
////		//Update vertex address in the subseg bridges.
////		itbridge     = vertex_tag->get_bridges().begin();
////		itbridge_end = vertex_tag->get_bridges().end();
////
////		for( ; itbridge != itbridge_end; ++itbridge)
////			(*itbridge)->replace_vert(old_vert, new_vert);
////
////		//Finally, replace the entry in the vertex dictionnary.
////		m_sample_verts.erase(itvdict);
////		m_sample_verts.insert( std::make_pair(new_vert, vertex_tag) );
////
////	}
//
//	//Now, update the face dictonary and tags.
//
////		old_face = itf->first;
////		new_face = itf->second;
////		if(old_face == new_face) continue;
////
////		itfdict = m_restricted_faces.find(old_face);
////
////		//Should only happen if face is NOT restricted Delaunay
////		if(itfdict == m_restricted_faces.end()) {
////			assert(!new_face->isRestrictedDelaunay());
////			continue;
////		}
////
////		assert(new_face->isRestrictedDelaunay());
////		face_tag = itfdict->second;
////		face_tag->set_face(new_face);
//
//		//We only have to replace the face tag.
////		m_restricted_faces.erase(itfdict);
////		m_restricted_faces.insert( std::make_pair(new_face, face_tag) );
//
////	}
//
//}

void
SurfaceSampler::output_restricted_delaunay(const char* const filename) const
{

  FILE* out_file = fopen(filename, "w");
  if (!out_file)
    vFatalError("Cannot open output file for writing",
		"RefFaceSampler::output_restricted_delaunay");

  typedef vector<FaceQueueEntry*> FacePrint;
  typedef vector<Vert*> VertPrint;
  typedef map<Vert*, int> VertIndex;

  FacePrint faces_to_print;
  VertPrint verts_to_print;
  VertIndex vert_indices;
  FacePrint::iterator itf, itf_end;
  VertPrint::iterator itv, itv_end;

  //Compute the faces intersected by Voronoi edges. Also obtain vertices to print.

  int counter = 1;

  get_restricted_faces(faces_to_print);
  itf = faces_to_print.begin();
  itf_end = faces_to_print.end();

  for (; itf != itf_end; ++itf) {

    Face* face = (*itf)->getFace();

    for (int i = 0; i < face->getNumVerts(); i++) {

      Vert* vert = face->getVert(i);

      if (vert_indices.insert(std::make_pair(vert, counter)).second) {
	verts_to_print.push_back(vert);
	++counter;
      }

    }

  }

  fprintf(out_file, "# vtk DataFile Version 1.0\n");
  fprintf(out_file, "Surface sample data\n");
  fprintf(out_file, "ASCII\n");
  fprintf(out_file, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(out_file, "POINTS %d float\n",
	  static_cast<int>(verts_to_print.size()));

  itv = verts_to_print.begin();
  itv_end = verts_to_print.end();

  for (; itv != itv_end; ++itv) {

    Vert* vert = *itv;

    fprintf(out_file, "%lf %lf %lf\n", vert->x(), vert->y(), vert->z());

  }

  itf = faces_to_print.begin();

  int numTris = static_cast<int>(itf_end - itf);
  fprintf(out_file, "CELLS %d %d\n", numTris, 4 * numTris);

  for (; itf != itf_end; ++itf) {
    Face* face = (*itf)->getFace();
    fprintf(out_file, "3 %d %d %d\n",
	    vert_indices.find(face->getVert(0))->second - 1,
	    vert_indices.find(face->getVert(1))->second - 1,
	    vert_indices.find(face->getVert(2))->second - 1);
  }

  fprintf(out_file, "CELL_TYPES %d\n", numTris);
  for (int ii = 0; ii < numTris; ii++) {
    fprintf(out_file, "5\n");
  }

  //The following code outputs normals at every vertex.
  //For this to work, the normal must be properly set in
  //the SampleVertTag, something that is not done anymore.
  //Code is therefore commented out for now.

  //   fprintf(out_file, "Normals\n");
  //   fprintf(out_file, "%d\n", static_cast<int>(verts_to_print.size()));

  //   itv     = verts_to_print.begin();
  //   itv_end = verts_to_print.end();

  //   for( ; itv != itv_end; ++itv) {

  //     Vert* vert = *itv;
  //     assert(m_sample_verts.count(vert) == 1);

  //     SampleVertTag* vert_tag = (m_sample_verts.find(vert))->second;
  //     CubitVector* normal = vert_tag->get_normal();

  //     fprintf(out_file, "%lf %lf %lf\n",
  // 	    normal->x(), normal->y(), normal->z());

  //   }

  //   fprintf(out_file, "NormalAtVertices\n");
  //   fprintf(out_file, "%d\n", static_cast<int>(verts_to_print.size()));

  //   itv     = verts_to_print.begin();
  //   itv_end = verts_to_print.end();
  //   counter = 0;

  //   for( ; itv != itv_end; ++itv) {
  //     ++counter;
  //     fprintf(out_file, "%d %d\n", counter, counter);
  //   }
  fclose(out_file);

}

void
SurfaceSampler::output_faces(const char* const filename,
			     const set<Face*>& faces_to_print) const
{

  vector<Face*> faces;
  std::copy(faces_to_print.begin(), faces_to_print.end(),
	    std::back_inserter(faces));

  output_faces(filename, faces);

}

void
SurfaceSampler::output_faces(const char* const filename,
			     const vector<Face*>& faces_to_print) const
{

  FILE* out_file = fopen(filename, "w");
  if (!out_file)
    vFatalError("Cannot open output file for writing",
		"RefFaceSampler::output_restricted_delaunay");

  typedef vector<Face*> FacePrint;
  typedef vector<Vert*> VertPrint;
  typedef map<Vert*, int> VertIndex;

  VertPrint verts_to_print;
  VertIndex vert_indices;
  FacePrint::const_iterator itf, itf_end;
  VertPrint::const_iterator itv, itv_end;

  //Compute the faces intersected by Voronoi edges. Also obtain vertices to print.

  int counter = 1;

  itf = faces_to_print.begin();
  itf_end = faces_to_print.end();

  for (; itf != itf_end; ++itf) {

    Face* face = *itf;

    for (int i = 0; i < face->getNumVerts(); i++) {

      Vert* vert = face->getVert(i);

      if (vert_indices.insert(std::make_pair(vert, counter)).second) {
	verts_to_print.push_back(vert);
	++counter;
      }

    }

  }

  fprintf(out_file, "MeshVersionFormatted 1\n");
  fprintf(out_file, "Dimension 3\n");

  fprintf(out_file, "Vertices\n");
  fprintf(out_file, "%d\n", static_cast<int>(verts_to_print.size()));

  itv = verts_to_print.begin();
  itv_end = verts_to_print.end();

  for (; itv != itv_end; ++itv) {

    Vert* vert = *itv;

    fprintf(out_file, "%lf %lf %lf 2\n", vert->x(), vert->y(), vert->z());

  }

  itf = faces_to_print.begin();

  fprintf(out_file, "Triangles\n");
  fprintf(out_file, "%d\n", static_cast<int>(itf_end - itf));

  for (; itf != itf_end; ++itf) {

    Face* face = *itf;

    fprintf(out_file, "%d %d %d 2\n",
	    vert_indices.find(face->getVert(0))->second,
	    vert_indices.find(face->getVert(1))->second,
	    vert_indices.find(face->getVert(2))->second);

  }

  fprintf(out_file, "End\n");

  fclose(out_file);

}

void
SurfaceSampler::receiveCreatedCells(std::vector<Cell*>& createdCells)
{
  m_new_cells.clear();
  logMessage(3, "SurfaceSampler received word about %zu newly created cells.\n",
	     createdCells.size());
  std::vector<Cell*>::iterator iter = createdCells.begin(), iterEnd =
      createdCells.end();
  for (; iter != iterEnd; ++iter) {
    m_new_cells.push_back(*iter);
  }
}
void
SurfaceSampler::print_surface_sampler()
{
  printf("SS --------- surface area: %f\n", m_surface->area());
  SampleVertDict::iterator itv = m_sample_verts.begin(), itv_end =
      m_sample_verts.end();
  //		printf("Sample Verts: %d\n",m_sample_verts.size());
  for (; itv != itv_end; ++itv) {
    itv->first->printVertInfo();
    const CubitVector* cv = itv->second->get_normal();
    printf("normal %f %f %f curv %f\n", (cv->x()), (cv->y()), (cv->z()),
	   itv->second->get_curvature());
    std::set<FaceQueueEntry*, FaceQueueEntryComp> face_set =
	itv->second->get_ring_faces();
    //			printf("Num Faces %d\n",face_set.size());
    std::set<SubsegBridge*> ss_set = itv->second->get_bridges();
    for (std::set<FaceQueueEntry*, FaceQueueEntryComp>::iterator it =
	face_set.begin(); it != face_set.end(); ++it) {
      (*it)->getFace()->printFaceInfo();
    }
  }
  //	printf("Restricted Faces: %d\n",m_restricted_faces.size());
  //	RestrictedFaceDict::iterator itf = m_restricted_faces.begin(),
  //			itf_end = m_restricted_faces.end();
  //	//
  //	for( ; itf != itf_end; ++itf){
  //		printf("First %p %p %d %d %d\n",itf->first,itf->first->getFace(),
  //				itf->first->isDeleted(),itf->first->getFace()->isDeleted(),itf->first->getFace()->isRestrictedDelaunay());
  //		itf->first->getFace()->printFaceInfo();
  //		//		printf("%p ",itf->first->getStoredVert(0));
  //		//		itf->first->getStoredVert(0)->printVertInfo();
  //		//		printf("%p ",itf->first->getStoredVert(1));
  //		//		itf->first->getStoredVert(1)->printVertInfo();
  //		//		printf("%p ",itf->first->getStoredVert(2));
  //		//		itf->first->getStoredVert(2)->printVertInfo();
  //		//
  //		//		printf("Second \n");
  //		//		const Face * face =itf->second->get_face();
  //		//		const_cast<Face*>(face)->printFaceInfo();
  //		//			    	CubitVector cv = itf->second->get_face_normal();
  //		//			 tw   	printf("normal %f %f %f\n",(cv.x()),(cv.y()),(cv.z()));
  //		//	}
  //	}
}

////////////

VoronoiEdge::VoronoiEdge(Face* const face, Cell* const cell1, Cell* const cell2) :
    m_face(face), m_beg(NULL), m_end(NULL)
{

  assert(m_face == findCommonFace(cell1, cell2));
  voronoi_vert_coord(cell1, m_beg);
  voronoi_vert_coord(cell2, m_end);

}

VoronoiEdge::VoronoiEdge(Cell* const cell1, Cell* const cell2) :
    m_face(findCommonFace(cell1, cell2)), m_beg(), m_end()
{

  assert(m_face);
  voronoi_vert_coord(cell1, m_beg);
  voronoi_vert_coord(cell2, m_end);

}

void
VoronoiEdge::voronoi_vert_coord(Cell* const cell, CubitVector& coord)
{

  double circum_coord[3];

  switch (cell->getType())
    {

    case Cell::eTet:
      {
	TetCell* tet_cell = dynamic_cast<TetCell*>(cell);
	assert(tet_cell);
	tet_cell->calcCircumcenter(circum_coord);
	break;
      }
    case Cell::eTriBFace:
      {
	TriBFace* tri_face = dynamic_cast<TriBFace*>(cell);
	assert(tri_face);
	tri_face->calcCircumcenter(circum_coord);
	break;
      }
    default:
      vFatalError(
	  "Unexpected cell type\n",
	  "VoronoiEdge::VoronoiEdge(const Cell* const, const Cell* const)");
      break;
    }

  coord.set(circum_coord);

}
