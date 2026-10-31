#ifndef GRTS_Workset_hh
#define GRTS_Workset_hh

#include <set>
#include <list>

#include "GR_Classes.h"
#include "GR_Mesh.h"
#include "GR_iMesh_Classes.h"

namespace ITAPS_GRUMMP
{
  /**\brief A class for iterating over entities in blocks.
   *
   * As expected, this class implements iMesh array iterators, which are
   * an example of the Iterator pattern.
   */
  class Workset : public iBase_EntityArrIterator_Private {
    Workset(const Workset&);
    Workset&
    operator=(const Workset&);
  protected:
    /// The size of data block to be returned.
    int iSize;
    /// Spatial dimension of the mesh that we're iterating over.
    int iDim;
    /**\brief An enumeration of data type that we can iterate over.
     *
     * In addition to treating tris and quads differently in 2D than in
     * 3D, this enumeration also contains entries for all faces, all
     * cells, and all entities.  eImpossible is used for things like
     * iteration over regions in a 2D mesh.  eUndefined is the default,
     * but this value is always overwritten in the constructor (assuming
     * the constructor succeeds).
     */
    enum eWSDataType {
      eVerts,
      eEdgeFaces,
      eTriFaces,
      eQuadFaces,
      ePolygonFaces,
      eAllFaces,
      eTriCells,
      eQuadCells,
      ePolygonCells,
      eTetCells,
      ePyrCells,
      ePrismCells,
      eHexCells,
      eSeptaCells,
      ePolyhedra,
      eAllCells2D,
      eAllCells3D,
      eEverything,
      eImpossible,
      eUndefined
    };
    /// Data type that we're iterating over.
    enum eWSDataType eDT;
    /// For iteration over all entities, type we're currently iterating over.
    enum eWSDataType eCurrentDT;
  public:
    /// Constructor based on type, topology, block size, and mesh dimension.
    Workset(const int eEType, const int eETopo, const int iSzIn,
	    const int iDim);

    /// Destructor; doesn't need to clean up anything.

    virtual
    ~Workset()
    {
    }

    /// Retrieve the next chunk of data from the iterator.
    virtual bool
    qGetNextBlock(iBase_EntityHandle** entity_handles,
		  int* entity_handles_allocated, int* entity_handles_size)
		  = 0;

    /// Reset iterator to the beginning.
    virtual bool
    qReset() = 0;
  };

/// A class for block iteration of RefImpl meshes.
  class MeshWorkset : public Workset {
    MeshWorkset(const MeshWorkset&);
    MeshWorkset&
    operator=(const MeshWorkset&);
    // State information for worksets in meshes.
    /// The mesh being iterated over.
    Mesh *pM;
    /**\brief The index of the next entity.
     *
     * The iterator also knows the type being iterated over.  Both
     * pieces of information are required.
     */
    int iNext;
    /**\brief The index of the last entity.
     *
     * The iterator also knows the type being iterated over.  Both
     * pieces of information are required.
     */
    int iTotal;
  public:
    /// Construct a block iterator for this mesh.
    MeshWorkset(Mesh * const pMIn, const int eEType, const int eETopo,
		const int iSzIn);

    // Inherited from Workset
    virtual bool
    qGetNextBlock(iBase_EntityHandle** entity_handles,
		  int* entity_handles_allocated, int* entity_handles_size);
    virtual bool
    qReset();
  };

/// A class for block iteration over RefImpl (unordered) sets.
  class EntitySetWorkset : public Workset {
    EntitySetWorkset(const EntitySetWorkset&);
    EntitySetWorkset&
    operator=(const EntitySetWorkset&);
    // State information for worksets in meshes.
    EntitySet *pES;

    // These iterators over STL sets are all of different type, so there
    // are a bunch of them instead of just an integer or two.

    /// Iterator over vertices in the set.
    std::set<VertRec>::iterator itVR;

    /// Iterator over edges in the set.
    std::set<EdgeFace*>::iterator itEdgeF;

    /// Iterator over triangular faces in the set (3D mesh).
    std::set<TriFace*>::iterator itTriF;

    /// Iterator over quadrilateral faces in the set (3D mesh).
    std::set<QuadFace*>::iterator itQuadF;

    /// Iterator over triangular cells in the set (2D mesh).
    std::set<TriCell*>::iterator itTriC;

    /// Iterator over quadrilateral cells in the set (3D mesh).
    std::set<QuadCell*>::iterator itQuadC;

    /// Iterator over tetrahedra in the set.
    std::set<TetCell*>::iterator itTetC;

    /// Iterator over pyramids in the set.
    std::set<PyrCell*>::iterator itPyrC;

    /// Iterator over prisms in the set.
    std::set<PrismCell*>::iterator itPrismC;

    /// Iterator over hexahedra in the set.
    std::set<HexCell*>::iterator itHexC;

    /// Last iterator value for vertices in the set.
    std::set<VertRec>::iterator itVREnd;

    /// Last iterator value for edges in the set.
    std::set<EdgeFace*>::iterator itEdgeFEnd;

    /// Last iterator value for triangular faces in the set (3D mesh).
    std::set<TriFace*>::iterator itTriFEnd;

    /// Last iterator value for quadrilateral faces in the set (3D mesh).
    std::set<QuadFace*>::iterator itQuadFEnd;

    /// Last iterator value for triangular cells in the set (2D mesh).
    std::set<TriCell*>::iterator itTriCEnd;

    /// Last iterator value for quadrilateral cells in the set (3D mesh).
    std::set<QuadCell*>::iterator itQuadCEnd;

    /// Last iterator value for tetrahedra in the set.
    std::set<TetCell*>::iterator itTetCEnd;

    /// Last iterator value for pyramids in the set.
    std::set<PyrCell*>::iterator itPyrCEnd;

    /// Last iterator value for prisms in the set.
    std::set<PrismCell*>::iterator itPrismCEnd;

    /// Last iterator value for hexahedra in the set.
    std::set<HexCell*>::iterator itHexCEnd;

  public:
    /// Construct a block iterator for the given set.
    EntitySetWorkset(EntitySet * const pESIn, const int eEType,
		     const int eETopo, const int iSzIn);

    /// Destructor.

    ~EntitySetWorkset();

    // Inherited from Workset
    bool
    qReset();
    bool
    qGetNextBlock(iBase_EntityHandle** entity_handles,
		  int* entity_handles_allocated, int* entity_handles_size);

    /// Updates the iterator when an entity is deleted from the set.
    void
    vDelete(Entity* pEToDelete);
  };

/// A class for block iteration over RefImpl (ordered) sets.
  class EntityListWorkset : public Workset {
    EntityListWorkset(const EntityListWorkset&);
    EntityListWorkset&
    operator=(const EntityListWorkset&);
    // State information for worksets in meshes.
    /// The EntityList being iterated over.
    EntityList *pEL;
    /// The number of returnable entities encountered in the current iteration.
    int iFound;
    /// Total number of returnable entities expected in the current iteration.
    int iTotal;
    /**\brief A flag to tell qGetNextBlock to increment the STL iterator
     * before getting the next entity.
     */
    bool qContinue;
    /// The STL iterator used internally to move through the data.
    std::list<Entity*>::iterator iter;
    /// The STL iterator used internally to mark the end of the data.
    std::list<Entity*>::iterator iterEnd;

    /// Update a type-based iterator for entity deletion from the set.
    void
    vAdjustIterForDeletionByType(std::list<Entity*>::iterator itToDelete,
				 const int iType);
    /// Update a topology-based iterator for entity deletion from the set.
    void
    vAdjustIterForDeletionByTopo(std::list<Entity*>::iterator itToDelete,
				 const int iTopo);
    /// Update a topology-based iterator for entity addition to the set.
    void
    vAddEntsByTopo(iBase_EntityHandle** entity_handles, int iTopo, int& iCurr);
    /// Update a type-based iterator for entity addition to the set.
    void
    vAddEntsByType(iBase_EntityHandle** entity_handles, int iType, int& iCurr);
  public:
    /// Constructs an ordered set iterator with given type, topo, size, and set.
    EntityListWorkset(EntityList * const pESIn, const int eEType,
		      const int eETopo, const int iSzIn);
    /// Destructor

    ~EntityListWorkset();

    // Inherited from Workset
    virtual bool
    qGetNextBlock(iBase_EntityHandle** entity_handles,
		  int* entity_handles_allocated, int* entity_handles_size);
    virtual bool
    qReset();

    /// EntityList is a friend so that vDelete and vAdded can be private.
    friend class EntityList;
  private:
    /// Delete the entity that the given iterator points to.
    void
    vDelete(std::list<Entity*>::iterator itToDelete);

    /// Add an entity to the iterator.
    void
    vAdded(Entity* pEAdded);
  };
}
#endif
