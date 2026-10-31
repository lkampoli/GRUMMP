#ifndef GR_TETMESHINITIALIZER_H
#define GR_TETMESHINITIALIZER_H 1
//
//#include "GR_config.h"
//#include <vector>
//#include <set>
//
//class Cell;
//class CubitBox;
//class CubitVector;
//class TetCell;
//class Vert;
//class VolMesh;
//template <class X> class KDDTree;
//
////This class is used to initialize and insert vertices into
////a tetrahedral mesh instantiated somewhere else.
//
//class TetMeshInitializer {
//
//  //A pointer to the volume mesh we wish to initialize.
// public:
//  VolMesh* m_mesh;
// private:
//  //Flag that lets one know if a tree exist to help find seed cells.
//  bool m_using_tree;
//
//  //A spatial indexing tree to accelerate seed cell computation.
//  KDDTree<TetCell*>* m_cell_tree;
//
// public:
//
//  TetMeshInitializer(const CubitBox& box,
//		     const bool using_tree = false);
//
//  ~TetMeshInitializer();
//
//  //Inserts a new_vertex, which is a copy of vertex, in m_mesh.
//  //The new_vertex is located at vertex's coords.
//  //It preserves the vertex's "properties" (flags, parent topology, etc.).
//  Vert* insert_vert_in_mesh(const Vert* const vertex,
//			    std::vector<Cell*>* new_cells = NULL);
//
//  Vert* insert_vert_in_mesh(const Vert* const vertex,
//			    const Vert* const vert_guess1,
//			    const Vert* const vert_guess2 = NULL,
//			    const Vert* const vert_guess3 = NULL,
//			    std::vector<Cell*>* new_cells = NULL);
//
//  Vert* insert_vert_in_mesh(const Vert* const vertex,
//			    const std::set<const Cell*>& seed_guesses,
//			    std::vector<Cell*>* new_cells = NULL);
//
//  Vert* insert_vert_in_mesh(const Vert* const vertex,
//			    const Cell* known_seed,
//			    std::vector<Cell*>* new_cells = NULL);
//
//  //Checks if vertex is inside cell's circumsphere.
//  static bool cell_is_seed(const Vert* const vertex,
//			   const Cell* const cell);
//
//  //Returns the pointer to the volume mesh.
//  VolMesh* get_mesh() const { return m_mesh; }
//
//  //Purges deleted entities from the mesh and rebuilds the tree.
//  void purge_mesh();
//
// private:
//
//  //These are made private.
//  TetMeshInitializer();
//  TetMeshInitializer(const TetMeshInitializer&);
//  TetMeshInitializer& operator=(const TetMeshInitializer&);
//
//  //Adds a cell to the cell tree.
//  void add_to_cell_tree(Cell* const cell);
//
//  //Removes a cell from the cell tree.
//  void remove_from_cell_tree(Cell* const cell);
//
//  //Find a TetCell to be used as seed in Watson insertion using the spatial indexing tree.
//  Cell* find_seed_cell_tree(const Vert* const vertex) const;
//
//  //Find a TetCell to be used as seed in Watson insertion by looking at every cell.
//  Cell* find_seed_cell_naively(const Vert* const vertex) const;
//
//  //Find all the cells having their bounding box intersected by box;
//  void find_cells_in_box(const CubitBox& box,
//			 std::vector<TetCell*>& cells_in_box) const;
//
//};

#endif
