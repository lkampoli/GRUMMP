#ifndef GRTS_EntIter_hh
#define GRTS_EntIter_hh

#include <set>
#include <list>
#include <assert.h>

#include "GR_Classes.h"
#include "GR_Mesh.h"
#include "GR_iMesh_Classes.h"

namespace ITAPS_GRUMMP
{
  /**\brief A class for iterating over entities one at a time.
   *
   * As expected, this class implements iMesh single-entity iterators,
   * which are an example of the Iterator pattern.  Note that it is not,
   * semantically or syntactically, the same as STL iterators.  The most
   * notable differences are the way in which the end of iteration is
   * indicated and the way in which the entity associated with an
   * iterator is returned.
   */
  class EntIter : public iBase_EntityIterator_Private {
  protected:
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
    enum eWSDataType eExtraDT;
    /// The entity the iterator currently points to.
    Entity *pE;
  private:
    /// Default constructor disabled.
    EntIter();
    /// Copy constructor disabled.
    EntIter(const EntIter&);
    /// Operator= disabled.
    const EntIter&
    operator=(const EntIter&);
  public:
    /// Constructs an iterator for the given type, topology, and mesh dimension.
    EntIter(const int eEType, const int eETopo, const int iDim);
    /// Destructor; doesn't need to clean up anything.
    virtual
    ~EntIter()
    {
    }
    /// Increment the iterator.
    virtual bool
    qIncrement() = 0;
    /// Return the data pointed to by the iterator.
    iBase_EntityHandle
    pvData() const
    {
      return pE;
    }
    /// Reset the iterator to the first entity of the given type and topology.
    virtual bool
    qReset() = 0;
  };

/// A class for single-entity iteration over RefImpl meshes.
  class MeshIter : public EntIter {
    MeshIter(const MeshIter&);
    MeshIter&
    operator=(const MeshIter&);
    // State information for worksets in meshes.
    /// The mesh we're iterating over.
    const Mesh *pM;
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
    int iLast;
  public:
    /// Construct an iterator over a RefImpl mesh.
    MeshIter(Mesh * const pMIn, const int eEType, const int eETopo);

    // Inherited from EntIter
    virtual bool
    qIncrement();
    virtual bool
    qReset();
  };

/// A class for single-entity iteration over RefImpl (unordered) sets.
  class EntitySetIter : public EntIter {
    EntitySetIter(const EntitySetIter&);
    EntitySetIter&
    operator=(const EntitySetIter&);
    // State information for worksets in entity sets.
    /// The set being iterated over.
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

  public:
    /// Constructs an unordered set iterator for a given type, topo, and set.
    EntitySetIter(EntitySet * const pESIn, const int eEType, const int eETopo);

    /// Destructor.
    ~EntitySetIter();

    // Inherited from EntIter
    virtual bool
    qIncrement();
    virtual bool
    qReset();

    /// Updates the iterator when an entity is deleted from the set.
    void
    vDelete(Entity* pEToDelete);
  };

/// A class for single-entity iteration over RefImpl (ordered) sets.
  class EntityListIter : public EntIter {
    EntityListIter(const EntityListIter&);
    EntityListIter&
    operator=(const EntityListIter&);
    // State information for worksets in entity lists.
    /// The EntityList being iterated over.
    EntityList *pEL;
    /// The number of returnable entities encountered in the current iteration.
    int iFound;
    /// Total number of returnable entities expected in the current iteration.
    int iTotal;
    /// The STL iterator used internally to move through the data.
    std::list<Entity*>::iterator iter;

    /// Update a type-based iterator for entity deletion from the set.
    void
    vAdjustIterForDeletionByType(std::list<Entity*>::iterator itToDelete,
				 const int iType);
    /// Update a topology-based iterator for entity deletion from the set.
    void
    vAdjustIterForDeletionByTopo(std::list<Entity*>::iterator itToDelete,
				 const int iTopo);
  public:
    /// Constructs an ordered set iterator for a given type, topo, and set.
    EntityListIter(EntityList * const pELIn, const int eEType,
		   const int eETopo);

    /// Destructor.
    ~EntityListIter();

    // Inherited from EntIter.
    virtual bool
    qIncrement();
    virtual bool
    qReset();

    /// Delete the entity that the given iterator points to.
    void
    vDelete(std::list<Entity*>::iterator itToDelete);

    /// Add an entity to the iterator.
    void
    vAdded(Entity* pEAdded);
  };
}
#endif
