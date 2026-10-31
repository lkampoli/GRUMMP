#ifndef GRTS_EntitySet_hh
#define GRTS_EntitySet_hh

#include <vector>
#include <set>
#include <list>
#include <assert.h>
#include "GR_Classes.h"
#include "GR_iMesh_Classes.h"

namespace ITAPS_GRUMMP
{
  /**\brief A class for vertices that are in a set implicitly.
   *
   * The iMesh spec says that downward adjacency calls on a set other
   * than the root set that request vertices must return all vertices on
   * the closure of entities in the set, even if those vertices haven't
   * been explicitly added to the set.  This data structure exists to
   * provide a way to store this data.  Each vertex is reference
   * counted, so that if all the higher-dimension entities that have the
   * vert on their closure are removed from the set, the vert can be
   * removed, too.
   */
  class VertRec {
  private:
    Vert *pV;
    int iC;
  public:
    /// Constructor
    VertRec(Vert * const pVIn = pVInvalidVert) :
	pV(pVIn), iC(1)
    {
    }
    /// Return the vertex associated with this record.
    operator Vert*() const
    {
      return pV;
    }

    /// Return the vertex associated with this record.
    Vert*
    pVVert() const
    {
      return pV;
    }

    /// Increase the reference count for this vertex
    int
    iIncrementCount()
    {
      return ++iC;
    }

    /// Decrease the reference count for this vertex
    int
    iDecrementCount()
    {
      return --iC;
    }

    /// Return the reference count for this vertex
    int
    iCount() const
    {
      return iC;
    }

    /// Reset the reference count for this vertex to a given value
    void
    vSetCount(const int i)
    {
      iC = i;
    }

    /// Comparison between vertex records.
    bool
    operator<(const VertRec& VR) const
    {
      return (pV < VR.pV);
    }
  };

  /**\brief A base class for all RefImpl entity sets.
   *
   * In addition to providing an interface to the specialized entity set
   * classes, this class directly handles set membership and set
   * heirarchies.
   */
  class EntitySetBase : public ::iBase_EntitySetHandle_Private {
  protected:
    std::set<EntitySetBase*> spESParents, spESChildren;
    std::set<EntitySetBase*> spESMembers, spESOwners;
    int iMeshTopoDim;
    bool qVertNumsOK;

    void
    vRebuildAuxiliaryVertList();
  private:
    const EntitySetBase&
    operator=(const EntitySetBase&)
    {
      assert(0);
      return (*this);
    }
  public:
    /// Construct for a given topological dimension
    EntitySetBase(int iDimIn) :
	spESParents(), spESChildren(), spESMembers(), spESOwners(), iMeshTopoDim(
	    iDimIn), qVertNumsOK(false)
    {
      assert(iDimIn == 2 || iDimIn == 3);
    }

    /// Copy constructor
    EntitySetBase(const EntitySetBase& ESB) :
	iBase_EntitySetHandle_Private(), spESParents(), spESChildren(), spESMembers(), spESOwners(), iMeshTopoDim(
	    ESB.iMeshTopoDim), qVertNumsOK(ESB.qVertNumsOK)
    {
      std::set<EntitySetBase*>::iterator iter = ESB.spESMembers.begin(),
	  iterEnd = ESB.spESMembers.end();
      for (; iter != iterEnd; iter++) {
	(*iter)->spESOwners.insert(this);
	spESMembers.insert(*iter);
      }
    }

    /// Destructor
    virtual
    ~EntitySetBase();

    /// Return the topological dimension of the mesh this set belongs to.
    int
    iMeshTopologicalDimension() const
    {
      return iMeshTopoDim;
    }

    /// Implements iMesh_getNumEntSets
    int
    iGetNumEntSets(const int iNumHops) const;

    /// Implements iMesh_getNumChldn
    int
    iGetNumChildren(const int iNumHops) const;

    /// Implements iMesh_getNumPrnts
    int
    iGetNumParents(const int iNumHops) const;

    /// Checks for descendants at any depth, to prevent circularity.
    bool
    qIsMyDescendant(EntitySetBase* pESB) const;

    /// Implements iMesh_isEntSetContained
    bool
    qIsMyMember(EntitySetBase* const pESB) const
    {
      return spESMembers.find(pESB) != spESMembers.end();
    }

    /// Implements iMesh_isChildOf
    bool
    qIsMyChild(EntitySetBase* const pESB) const
    {
      return spESChildren.find(pESB) != spESChildren.end();
    }

    /// Implements one side of iMesh_addPrntChld
    void
    vAddParent(EntitySetBase* const pESB)
    {
      spESParents.insert(pESB);
    }

    /// Implements one side of iMesh_addPrntChld
    void
    vAddChild(EntitySetBase* const pESB)
    {
      spESChildren.insert(pESB);
    }

    /// Implements one side of iMesh_rmvPrntChld
    void
    vRemoveParent(EntitySetBase* const pESB)
    {
      spESParents.erase(pESB);
    }

    /// Implements one side of iMesh_rmvPrntChld
    void
    vRemoveChild(EntitySetBase* const pESB)
    {
      spESChildren.erase(pESB);
    }

    /// Implements one side of iMesh_rmvEntSet
    void
    vRemoveMember(EntitySetBase* const pESB)
    {
      spESMembers.erase(pESB);
    }

    /// Implements one side of iMesh_rmvEntSet
    void
    vRemoveOwner(EntitySetBase* const pESB)
    {
      spESOwners.erase(pESB);
    }

    /// Internal call for retrieving set members (recursively)
    void
    vGetEntitySets(std::set<EntitySetBase*>& spESRes, const int iNumHops) const;

    /// Implements iMesh_getEntSets
    void
    vGetEntitySets(iBase_EntitySetHandle** ent_set_handles,
		   int* ent_set_handles_allocated, int* iNumSets,
		   const int iNumHops) const;

    /// Implements iMesh_addEntSet
    void
    vAddEntitySet(void* ent_set_handle);

    /// Internal helper function for iMesh_addEntSet
    void
    vAddEntitySet(EntitySetBase* pESB);

    /// Implements iMesh_rmvEntSet
    void
    vRemoveEntitySet(void* ent_set_handle);

    /// Implements iMesh_getChldn
    void
    vGetChildren(const int iHops, iBase_EntitySetHandle** ent_set_handles,
		 int* ent_set_handles_allocated, int* iNumChildren) const;

    /// Implements iMesh_getPrnts
    void
    vGetParents(const int iHops, iBase_EntitySetHandle** ent_set_handles,
		int* ent_set_handles_allocated, int* iNumParents) const;

  private:
    /// Internal function to help with iMesh_getChldn
    void
    vGetChildren(std::set<EntitySetBase*>& spESResult, const int iHops) const;

    /// Internal function to help with iMesh_getPrnts
    void
    vGetParents(std::set<EntitySetBase*>& spESResult, const int iHops) const;

  public:
    /// Implements iMesh_getEntities with non-root set argument
    virtual void
    vGetEntities(const int entity_type, const int entity_topology,
		 iBase_EntityHandle** entity_handles,
		 int *entity_handles_allocated,
		 int *entity_handles_size) const = 0;

    /// Implements iMesh_addEntArrToSet
    virtual void
    vAddEntities(const iBase_EntityHandle entity_handles[],
		 const int iNumEnt) = 0;

    /// Implements iMesh_addEntToSet
    virtual void
    vAddEntity(void* const entity_handle)
    = 0;

    /// Implements iMesh_rmvEntArrFromSet
    virtual void
    vRemoveEntities(const iBase_EntityHandle entity_handles[],
		    const int iNumEnt) = 0;

    /// Implements iMesh_rmvEntFromSet
    void
    vRemoveEntity(void* const entity_handle)

    {
      Entity* pE = static_cast<Entity*>(entity_handle);
      vRemoveEntity(pE);
    }

    /// Internal helper function for iMesh_rmvEntFromSet
    virtual void
    vRemoveEntity(Entity* const pE) = 0;

    /// Implements iMesh_isEntContained
    virtual bool
    qIsEntityInSet(void* const entity_handle) const = 0;

    /// Returns number of vertices in set explicitly
    virtual int
    iNumVerts() const = 0;

    /// Returns number of vertices in set implicitly
    virtual int
    iNumAuxVerts() const = 0;

    /// Returns number of edges in set
    virtual int
    iNumEdges() const = 0;

    /// Returns number of triangular faces in set
    virtual int
    iNumTriFaces() const = 0;

    /// Returns number of quadrilateral faces in set
    virtual int
    iNumQuadFaces() const = 0;

    /// Returns number of triangular cells in set
    virtual int
    iNumTriCells() const = 0;

    /// Returns number of quadrilateral cells in set
    virtual int
    iNumQuadCells() const = 0;

    /// Returns number of tetrahedral cells in set
    virtual int
    iNumTetCells() const = 0;

    /// Returns number of pyramidal cells in set
    virtual int
    iNumPyrCells() const = 0;

    /// Returns number of prismatic cells in set
    virtual int
    iNumPrismCells() const = 0;

    /// Returns number of hexahedral cells in set
    virtual int
    iNumHexCells() const = 0;

    /// Returns true if the set is an ordered set.
    virtual bool
    qIsOrdered() const = 0;

    /// Implements iMesh_intersect
    friend void
    vIntersect(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	       EntitySetBase*& pESBRes);
  };

  class EntitySet;
  class EntitySetIter;
  class EntitySetWorkset;
  class EntityList;
  class EntityListIter;
  class EntityListWorkset;

  /**\brief Specialization of entity sets for the unique, unordered case.
   */
  class EntitySet : public EntitySetBase {
  private:
    std::set<VertRec> sVR;
    std::set<VertRec> sVRAux;
    std::set<EdgeFace*> spEdgeF;
    std::set<TriFace*> spTriF;
    std::set<QuadFace*> spQuadF;
    std::set<TriCell*> spTriC;
    std::set<QuadCell*> spQuadC;
    std::set<TetCell*> spTetC;
    std::set<PyrCell*> spPyrC;
    std::set<PrismCell*> spPrismC;
    std::set<HexCell*> spHexC;
    std::set<EntitySetIter*> spIters;
    std::set<EntitySetWorkset*> spWorksets;

    void
    vAddVertex(Vert* const pV, const bool qFromEntity);
    void
    vRemoveVertex(Vert* const pV, const bool qFromEntity);

    void
    vRebuildAuxiliaryVertList();

    const EntitySet&
    operator=(const EntitySet&)
    {
      assert(0);
      return (*this);
    }
  public:
    /// Default constructor
    EntitySet(const int iDim = 2);

    /// Construct a set from a list; needed for some of the set boolean ops.
    EntitySet(const EntityList& EL);

    /// Copy constructor
    EntitySet(const EntitySet& ES);

    /// Destructor
    ~EntitySet()
    {
    }

    virtual void
    vGetEntities(const int entity_type, const int entity_topology,
		 iBase_EntityHandle** entity_handles,
		 int *entity_handles_allocated, int *entity_handles_size) const;

  private:
    /// Internal helper function to facillitate implementation of iMesh_getEnts
    void
    vGetEntities(const int entity_type, const int entity_topology,
		 std::set<Entity*>& spE) const;
  public:
    virtual void
    vAddEntities(const iBase_EntityHandle entity_handles[], const int iNumEnt);
    virtual void
    vAddEntity(void* const entity_handle);
    virtual void
    vRemoveEntity(Entity* const pE);
    virtual void
    vRemoveEntities(const iBase_EntityHandle entity_handles[],
		    const int iNumEnt);
    virtual bool
    qIsEntityInSet(void* const entity_handle) const;

    virtual bool
    qIsOrdered() const
    {
      return false;
    }

    virtual int
    iNumVerts() const
    {
      return sVR.size();
    }
    virtual int
    iNumAuxVerts() const
    {
      return sVRAux.size();
    }
    virtual int
    iNumEdges() const
    {
      return spEdgeF.size();
    }
    virtual int
    iNumTriFaces() const
    {
      return spTriF.size();
    }
    virtual int
    iNumQuadFaces() const
    {
      return spQuadF.size();
    }
    virtual int
    iNumTriCells() const
    {
      return spTriC.size();
    }
    virtual int
    iNumQuadCells() const
    {
      return spQuadC.size();
    }
    virtual int
    iNumTetCells() const
    {
      return spTetC.size();
    }
    virtual int
    iNumPyrCells() const
    {
      return spPyrC.size();
    }
    virtual int
    iNumPrismCells() const
    {
      return spPrismC.size();
    }
    virtual int
    iNumHexCells() const
    {
      return spHexC.size();
    }

    /// Implements one case of iMesh_unite
    friend void
    vUnite(const EntitySet& ES1, const EntitySet& ES2, EntitySet& ESRes);

    /// Implements one case of iMesh_subtract
    friend void
    vSubtract(const EntitySet& ES1, const EntitySet& ES2, EntitySet& ESRes);

    /// Implements one case of iMesh_intersect
    friend void
    vIntersect(const EntitySet& ES1, const EntitySet& ES2, EntitySet& ESRes);

    friend class EntitySetIter;
    friend class EntitySetWorkset;
  };

  /**\brief Specialization of entity sets for the non-unique, ordered case.
   */
  class EntityList : public EntitySetBase {
  private:
    std::list<Entity*> lpEntities;
    std::set<VertRec> sVRAux;
    std::set<EntityListIter*> spIters;
    std::set<EntityListWorkset*> spWorksets;

    int iVerts, iEdges, iTriFaces, iQuadFaces, iTriCells, iQuadCells, iTetCells,
	iPyrCells, iPrismCells, iHexCells;

    void
    vAddVertex(Vert* const pV);
    void
    vRemoveVertex(Vert* const pV);

    void
    vRebuildAuxiliaryVertList();

    const EntityList&
    operator=(const EntityList&)
    {
      assert(0);
      return (*this);
    }
  public:
    EntityList(const int iDim = 2);
    EntityList(const EntityList& EL);
    ~EntityList()
    {
    }

    virtual void
    vGetEntities(const int entity_type, const int entity_topology,
		 iBase_EntityHandle** entity_handles,
		 int *entity_handles_allocated, int *entity_handles_size) const;
  private:
    /// Internal helper function to facillitate implementation of iMesh_getEnts
    void
    vGetEntities(const int entity_type, const int entity_topology,
		 std::vector<Entity*>& vpE) const;
  public:
    virtual void
    vAddEntities(const iBase_EntityHandle entity_handles[], const int iNumEnt);
    virtual void
    vAddEntity(void* const entity_handle);
    virtual void
    vRemoveEntity(Entity* pE);
    virtual void
    vRemoveEntities(const iBase_EntityHandle entity_handles[],
		    const int iNumEnt);
    virtual bool
    qIsEntityInSet(void* entity_handle) const;

    virtual bool
    qIsOrdered() const
    {
      return true;
    }

    virtual int
    iNumVerts() const
    {
      return iVerts;
    }
    virtual int
    iNumAuxVerts() const
    {
      return sVRAux.size();
    }
    virtual int
    iNumEdges() const
    {
      return iEdges;
    }
    virtual int
    iNumTriFaces() const
    {
      return iTriFaces;
    }
    virtual int
    iNumQuadFaces() const
    {
      return iQuadFaces;
    }
    virtual int
    iNumTriCells() const
    {
      return iTriCells;
    }
    virtual int
    iNumQuadCells() const
    {
      return iQuadCells;
    }
    virtual int
    iNumTetCells() const
    {
      return iTetCells;
    }
    virtual int
    iNumPyrCells() const
    {
      return iPyrCells;
    }
    virtual int
    iNumPrismCells() const
    {
      return iPrismCells;
    }
    virtual int
    iNumHexCells() const
    {
      return iHexCells;
    }

    /// Implements one case of iMesh_unite
    friend void
    vUnite(const EntityList& EL1, const EntityList& EL2, EntityList& ELRes);

    /// Implements one case of iMesh_subtract
    friend void
    vSubtract(const EntityList& EL1, const EntityList& EL2, EntityList& ELRes);

    /// Implements one case of iMesh_intersect
    friend void
    vIntersect(const EntityList& EL1, const EntityList& EL2, EntityList& ELRes);

    // Makes construction of sets from lists -much- easier.
    friend
    EntitySet::EntitySet(const EntityList& EL);

    friend class EntityListIter;
    friend class EntityListWorkset;
  };

/// Implements iMesh_unite through type-specific calls
  void
  vUnite(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	 EntitySetBase*& pESBRes);

/// Implements iMesh_subtract through type-specific calls
  void
  vSubtract(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	    EntitySetBase*& pESBRes);

/// Implements iMesh_intersect through type-specific calls
  void
  vIntersect(const EntitySetBase* pESB1, const EntitySetBase* pESB2,
	     EntitySetBase*& pESBRes);
}

#endif
