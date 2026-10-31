//#include "GR_TetMeshInitializer.h"
//#include "GR_AdaptPred.h"
//#include "GR_Geometry.h"
//#include "GR_Cell.h"
//#include "GR_Vertex.h"
//#include "GR_VolMesh.h"
//#include "CubitBox.hpp"
//#include "CubitVector.hpp"
//#include "DLIList.hpp"
//#include "KDDTree.hpp"
//
//#include <algorithm>
//#include <functional>
//#include <set>
//#include <vector>
//
//using std::set;
//using std::vector;
//
//TetMeshInitializer::
//TetMeshInitializer()
//  : m_mesh(NULL), m_using_tree(false), m_cell_tree(NULL) { assert(0); }
//
//TetMeshInitializer::
//TetMeshInitializer(const TetMeshInitializer&)
//  : m_mesh(NULL), m_using_tree(false), m_cell_tree(NULL) { assert(0); }
//
//TetMeshInitializer& TetMeshInitializer::
//operator=(const TetMeshInitializer&) { assert(0); return *this; }
//
//TetMeshInitializer::
//TetMeshInitializer(const CubitBox& box,
//		   const bool using_tree)
//  : m_mesh(new VolMesh(box)), m_using_tree(using_tree),
//    m_cell_tree(m_using_tree ? new KDDTree<TetCell*>() : NULL) {
//
//  for(GR_index_t i = 0; m_using_tree && i < m_mesh->getNumCells(); i++) {
//    TetCell* tet_cell = dynamic_cast<TetCell*>(m_mesh->getCell(i));
//    if(tet_cell && !tet_cell->isDeleted()) m_cell_tree->add(tet_cell);
//  }
//
//}
//
//TetMeshInitializer::
//~TetMeshInitializer() {
//
//  if(m_cell_tree) delete m_cell_tree;
//  if(m_mesh)      delete m_mesh;
//
//}
//
//Vert* TetMeshInitializer::
//insert_vert_in_mesh(const Vert* const vertex,
//		    std::vector<Cell*>* new_cells) {
//
//  assert(vertex->getSpaceDimen() == 3);
//
//  Cell* seed_cell = NULL;
//
//  if(m_using_tree) {
//    assert(m_cell_tree);
//    seed_cell = find_seed_cell_tree(vertex);
//  }
//  else {
//    seed_cell = find_seed_cell_naively(vertex);
//  }
//
//  assert(seed_cell);
//
//  return insert_vert_in_mesh(vertex, seed_cell, new_cells);
//
//}
//
//Vert* TetMeshInitializer::
//insert_vert_in_mesh(const Vert* const vertex,
//		    const Vert* const vert_guess1,
//		    const Vert* const vert_guess2,
//		    const Vert* const vert_guess3,
//		    vector<Cell*>* new_cells)
//{
//
//  assert(vertex->isValid());
//
//  set<const Cell*> guesses; set<Cell*>tmp_cells; set<Vert*> tmp_verts;
//
//  if(vert_guess1) {
//    findNeighborhoodInfo(vert_guess1, tmp_cells, tmp_verts);
//    std::copy(tmp_cells.begin(), tmp_cells.end(),
//	      std::inserter(guesses, guesses.begin()));
//  }
//  if(vert_guess2) {
//    findNeighborhoodInfo(vert_guess2, tmp_cells, tmp_verts);
//    std::copy(tmp_cells.begin(), tmp_cells.end(),
//	      std::inserter(guesses, guesses.begin()));
//  }
//  if(vert_guess3) {
//    findNeighborhoodInfo(vert_guess3, tmp_cells, tmp_verts);
//    std::copy(tmp_cells.begin(), tmp_cells.end(),
//	      std::inserter(guesses, guesses.begin()));
//  }
//
//  return insert_vert_in_mesh(vertex, guesses, new_cells);
//
//}
//
//Vert* TetMeshInitializer::
//insert_vert_in_mesh(const Vert* const vertex,
//		    const set<const Cell*>& seed_guesses,
//		    vector<Cell*>* new_cells)
//{
//
//  assert(vertex->getSpaceDimen() == 3);
//
//  set<const Cell*>::iterator it     = seed_guesses.begin(),
//                             it_end = seed_guesses.end();
//
//  Cell* seed = NULL;
//
//  for( ; it != it_end; ++it) {
//
//    //Not too elegant I admit:
//    Cell* cell = const_cast<Cell*>(*it);
//    assert(cell && !cell->isDeleted());
//
//    if(cell_is_seed(vertex, cell)) {
//      seed = cell;
//      break;
//    }
//
//  }
//
//  if(!seed)
//    seed = find_seed_cell_naively(vertex);
//
//  assert(seed);
//
//  return insert_vert_in_mesh(vertex, seed, new_cells);
//
//}
//
//Vert* TetMeshInitializer::
//insert_vert_in_mesh(const Vert* const vertex,
//		    const Cell* known_seed,
//		    vector<Cell*>* new_cells)
//{
//
//  assert(cell_is_seed(vertex, known_seed));
//  Vert* new_vertex = m_mesh->createVert(vertex);
//  //Compute Watson Data.
//  VolMesh::WatsonData watson_data;
//  m_mesh->computeWatsonData_deprecated(new_vertex, known_seed,
//			      watson_data, false, false);
//
//#ifndef NDEBUG
//
//  assert(!watson_data.cells_to_remove.empty());
//  assert(watson_data.encroached_bfaces.empty());
//
//  if(watson_data.cells_to_remove.size() == 1)
//    assert(watson_data.faces_to_remove.empty());
//  else
//    assert(!watson_data.faces_to_remove.empty());
//
//#endif
//
//  if(m_using_tree) {
//
//    assert(m_cell_tree);
//    //Remove cells to delete from m_cell_tree.
//    std::for_each( watson_data.cells_to_remove.begin(), watson_data.cells_to_remove.end(),
//		   std::bind1st(std::mem_fun(&TetMeshInitializer::remove_from_cell_tree), this) );
//    //Proceed with insertion.
//
//    vector<Cell*> local_new_cells;
//    qAdaptPred = true;
//    m_mesh->insertWatsonInterior_deprecated(new_vertex,
//				  watson_data.cells_to_remove,
//				  watson_data.faces_to_remove,
//				  watson_data.hull_faces,
//				  &local_new_cells);
//    qAdaptPred = false;
//
//    //Add newly created cells to m_cell_tree.
//    std::for_each( local_new_cells.begin(), local_new_cells.end(),
//		   std::bind1st(std::mem_fun(&TetMeshInitializer::add_to_cell_tree), this) );
//    if(new_cells) new_cells->swap(local_new_cells);
//
//  }
//
//  else {
//
//    qAdaptPred = true;
//
//    m_mesh->insertWatsonInterior_deprecated(new_vertex,
//				  watson_data.cells_to_remove,
//				  watson_data.faces_to_remove,
//				  watson_data.hull_faces,
//				  new_cells);
//    qAdaptPred = false;
//
//  }
//  return new_vertex;
//
//}
//
//bool TetMeshInitializer::
//cell_is_seed(const Vert* const vertex,
//	     const Cell* const cell) {
//
//  assert(cell->getType() == Cell::eTet);
//  return ( isInsphere(cell->getVert(0)->getCoords(),
//		     cell->getVert(1)->getCoords(),
//		     cell->getVert(2)->getCoords(),
//		     cell->getVert(3)->getCoords(),
//		     vertex->getCoords()) == 1 );
//
//}
//
//void TetMeshInitializer::
//purge_mesh() {
//
//  assert(m_mesh);
//  assert(m_cell_tree);
//
//  m_mesh->purgeAllEntities();
//
//  //After purging the mesh, we need to rebuild the search tree entirely.
//
//  if(m_using_tree) {
//
//    assert(m_cell_tree);
//    delete m_cell_tree;
//
//    m_cell_tree = new KDDTree<TetCell*>();
//
//    for(GR_index_t i = 0; i < m_mesh->getNumCells(); i++) {
//      TetCell* tet_cell = dynamic_cast<TetCell*>(m_mesh->getCell(i));
//      if(tet_cell && !tet_cell->isDeleted()) m_cell_tree->add(tet_cell);
//    }
//
//  }
//
//}
//
//void TetMeshInitializer::
//add_to_cell_tree(Cell* const cell) {
//
//  TetCell* tet_cell = dynamic_cast<TetCell*>(cell);
//  assert(tet_cell && !tet_cell->isDeleted());
//
//  m_cell_tree->add(tet_cell);
//
//}
//
//void TetMeshInitializer::
//remove_from_cell_tree(Cell* const cell) {
//
//  TetCell* tet_cell = dynamic_cast<TetCell*>(cell);
//  assert(tet_cell);
//
//  m_cell_tree->remove(tet_cell);
//
//}
//
//Cell* TetMeshInitializer::
//find_seed_cell_tree(const Vert* const vertex) const {
//
//  //printf("finding seed with tree\n");
//
//  assert(m_using_tree && m_cell_tree);
//
//  DLIList<TetCell*> seed_candidates;
//  CubitVector point_loc(vertex->x(), vertex->y(), vertex->z());
//
//  m_cell_tree->find(CubitBox(point_loc, point_loc), seed_candidates);
//
//  int size = seed_candidates.size();
//
//  for(int i = 0; i < size; i++) {
//
//    TetCell* cell = seed_candidates.next(i);
//    if(cell->isDeleted()) continue;
//
//    if(cell_is_seed(vertex, cell)) return cell;
//
//  }
//
//  return find_seed_cell_naively(vertex);
//
//}
//
//Cell* TetMeshInitializer::
//find_seed_cell_naively(const Vert* const vertex) const {
//
//  int num_cells = m_mesh->getNumCells();
//
//  for(int i = 0; i < num_cells; i++) {
//
//    Cell* cell = m_mesh->getCell(i);
//    if(!cell || cell->isDeleted()) continue;
//    if(cell_is_seed(vertex, cell)) return cell;
//
//  }
//
//#ifndef NDEBUG
//
//  logMessage(1,"Tried to insert a vertex at: ");
//	 vertex->printVertInfo();
//
//  assert(m_mesh->isValid());
//
//  for(int i = 0; i < num_cells; i++) {
//
//    Cell* cell = m_mesh->getCell(i);
//    if(!cell || cell->isDeleted()) continue;
//
//    int in_tet = isInTet(vertex, cell);
//
//    if(in_tet == -1) continue;
//
//    logMessage(1,"cell = %p: ", cell);
//    if(in_tet == 0) logMessage(1,"INVALID TET\n");
//    if(in_tet == 1) logMessage(1,"COINCIDES WITH VERT\n");
//    if(in_tet == 2 || in_tet == 3) logMessage(1,"COINCIDES WITH FACE\n");
//    if(in_tet == 4) logMessage(1,"IN TET\n");
//
//  }
//
//#endif
//
//  vFatalError("The point you are trying to insert is outside the meshed domain.",
//	      "TetMeshInitializer::find_seed_cell()");
//  return pCInvalidCell; // This exit is never taken.
//}
//
//void TetMeshInitializer::
//find_cells_in_box(const CubitBox& box,
//		  std::vector<TetCell*>& cells_in_box) const {
//
//  DLIList<TetCell*> query_return;
//  m_cell_tree->find(box, query_return);
//
//  int size = query_return.size();
//
//  for(int i = 0; i < size; i++)
//    cells_in_box.push_back(query_return.next(i));
//
//}
