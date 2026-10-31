#ifndef included_iMesh_GRUMMP_Base_hh
#define included_iMesh_GRUMMP_Base_hh

#include <assert.h>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include "GR_iMesh_Classes.h"
#include "GR_Mesh.h"
#include "GR_Mesh2D.h"
#include "GR_VolMesh.h"
// #include "iMesh_GRUMMP_EntitySet.hh"
// #include "iMesh_GRUMMP_Workset.hh"
// #include "iMesh_GRUMMP_EntIter.hh"
// #include "iMesh_GRUMMP_Tag.hh"
#include "iMesh_GRUMMP_misc.hh"

using namespace ITAPS_GRUMMP;

namespace ITAPS_GRUMMP
{
  /**\brief A base class for the GRUMMP iMesh classes.
   *
   * In addition to providing a base class for the 2D and volume mesh
   * iMesh_GRUMMP classes, this class maintains information about the
   * sets, iterators, and tags that are outstanding for the mesh, so
   * that these can be updated as appropriate when entities are deleted,
   * etc.
   */
  class iMesh_Base : public iBase_EntitySetHandle_Private {
  private:
    iMesh_Base(iMesh_Base&);
    iMesh_Base&
    operator=(const iMesh_Base&);
  protected:
    /// A pointer to the GRUMMP Mesh2D, Surface, or VolMesh object that stores
    /// the mesh database.
    Mesh *pM;
    bool externalData;
    int m_err;
  public:
    /// Return a pointer to the Mesh2D/Surface/VolMesh object to allow calls to it.
    Mesh *
    pMesh()
    {
      return pM;
    }

    /// An enumerator for the type of mesh.
    enum eMesh_Type {
      ePlane = 2, eSurface = -3, eVolume = 3, eDefault = -1
    };

    /// A factory that returns on demand an iMesh_Plane/Surface/Volume object. 
    static iMesh_Base*
    Create(const enum eMesh_Type eMTIn);

    /// Return most recent error code
    int
    getError() const
    {
      return m_err;
    }

    /// Set error code (called from C wrapper functions that aren't
    /// member functions of iMesh_Base or its descendants.
    void
    setError(const int err)
    {
      m_err = err;
    }
  protected:
    /// What type is this mesh?
    enum eMesh_Type eMT;

    /// A container for entity sets belonging to this mesh.
    std::set<EntitySetBase*> spEntSets;

    /// A container for block iterators belonging to this mesh.
    std::set<Workset*> spSetWorksets;

    /// A container for single entity iterators belonging to this mesh.
    std::set<EntIter*> spSetIters;

  public:
    /// Remove an entity from all sets and remove all its tags.
    void
    vRemoveAllReferencesTo(Entity* const pE);

    // Default constructor
    iMesh_Base(const enum eMesh_Type);

    // Virtual destructor (required)
    virtual
    ~iMesh_Base()
    {
      vDestroySetsTagsIters();
      if (pM && !externalData) {
	delete pM;
	pM = NULL;
      }
    }

  public:
    /// Recompute the entity counts; sometimes required to ensure correctness.
    void
    updateEntityCounts();

    /// Return number of vertices actually in use.
    unsigned int
    getNumVertices() const
    {
      return pM->getNumVertsInUse();
    }

    /// Return number of edges actually in use.
    unsigned int
    getNumLineSegments() const
    {
      return pM->getNumEdgesInUse();
    }

    /// Return number of triangles actually in use.
    unsigned int
    getNumTriangles() const
    {
      return pM->getNumTrisInUse();
    }

    /// Return number of quadrilaterals actually in use.
    unsigned int
    getNumQuadrilaterals() const
    {
      return pM->getNumQuadsInUse();
    }

    /// Return number of non-standard polygons actually in use.
    unsigned int
    getNumOtherPolygons() const
    {
      return 0;
    }

    /// Return number of tetrahedra actually in use.
    unsigned int
    getNumTetrahedra() const
    {
      return pM->getNumTetsInUse();
    }

    /// Return number of pyramids actually in use.
    unsigned int
    getNumPyramids() const
    {
      return pM->getNumPyrsInUse();
    }

    /// Return number of prisms actually in use.
    unsigned int
    getNumPrisms() const
    {
      return pM->getNumPrismsInUse();
    }

    /// Return number of hexahedra actually in use.
    unsigned int
    getNumHexahedra() const
    {
      return pM->getNumHexesInUse();
    }

    /// Return number of septahedra actually in use.
    unsigned int
    getNumSeptahedra() const
    {
      return 0;
    }

    /// Return number of non-standard polyhedra actually in use.
    unsigned int
    getNumOtherPolyhedra() const
    {
      return 0;
    }

  public:
    /// Does the work for iMesh_setGeometricDimension
    void
    setGeometricDimension(const int dim_in);

    /// Does the work for iMesh_createVert
    Vert*
    createVert(const double x, const double y, const double z)
    {
      if (eMT == ePlane) {
	return (pM->createVert(x, y, 0));
      }
      else {
	return (pM->createVert(x, y, z));
      }
    }

    /// Clean up all outstanding sets, tags, and iterators.
    void
    vDestroySetsTagsIters();

    /// Does the work for iMesh_destroyEntSet
    void
    destroyEntSet(/*in*/iBase_EntitySetHandle entity_set);

    /// Does the work for iMesh_isList
    bool
    isList(/*in*/iBase_EntitySetHandle entity_set);

    /// Does the work for iMesh_getNumEntSets
    int
    getNumEntSets(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/int num_hops);

    /// Does the work for iMesh_getEntSets
    void
    getEntSets(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/int num_hops,
	       /*inout*/iBase_EntitySetHandle** contained_entset_handles,
	       /*inout*/int* contained_entset_handles_allocated,
	       /*out*/int* contained_entset_handles_size);

    /// Does the work for iMesh_addEntToSet
    void
    addEntToSet(/*in*/iBase_EntityHandle entity_handle,
    /*inout*/iBase_EntitySetHandle entity_set);

    /// Does the work for iMesh_rmvEntFromSet
    void
    rmvEntFromSet(/*in*/iBase_EntityHandle entity_handle,
    /*inout*/iBase_EntitySetHandle entity_set);

    /// Does the work for iMesh_addEntArrToSet
    void
    addEntArrToSet(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/const int entity_handles_size,
		   /*inout*/iBase_EntitySetHandle entity_set);

    /// Does the work for iMesh_rmvEntArrFromSet
    void
    rmvEntArrFromSet(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/const int entity_handles_size,
		     /*inout*/iBase_EntitySetHandle entity_set);

    /// Does the work for iMesh_addEntSet
    void
    addEntSet(/*in*/iBase_EntitySetHandle entity_set_to_add,
    /*inout*/iBase_EntitySetHandle entity_set_handle);

    /// Does the work for iMesh_rmvEntSet
    void
    rmvEntSet(/*in*/iBase_EntitySetHandle entity_set_to_remove,
    /*inout*/iBase_EntitySetHandle entity_set_handle);

    /// Does the work for iMesh_isEntContained
    bool
    isEntContained(/*in*/iBase_EntitySetHandle containing_entity_set,
    /*in*/iBase_EntityHandle entity_handle);

    /// Does the work for iMesh_isEntSetContained
    bool
    isEntSetContained(/*in*/iBase_EntitySetHandle containing_entity_set,
    /*in*/iBase_EntitySetHandle contained_entity_set);

    /// Does the work for iMesh_addPrntChld
    void
    addPrntChld(/*inout*/iBase_EntitySetHandle parent_entity_set,
    /*inout*/iBase_EntitySetHandle child_entity_set);

    /// Does the work for iMesh_rmvPrntChld
    void
    rmvPrntChld(/*inout*/iBase_EntitySetHandle parent_entity_set,
    /*inout*/iBase_EntitySetHandle child_entity_set);

    /// Does the work for iMesh_isChildOf
    bool
    isChildOf(/*in*/iBase_EntitySetHandle parent_entity_set,
    /*in*/iBase_EntitySetHandle child_entity_set);

    /// Does the work for iMesh_getNumChld
    int
    getNumChld(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/int num_hops);

    /// Does the work for iMesh_getNumPrnt
    int
    getNumPrnt(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/int num_hops);

    /// Does the work for iMesh_getChldn
    void
    getChldn(/*in*/iBase_EntitySetHandle from_entity_set,
    /*in*/int num_hops,
	     /*inout*/iBase_EntitySetHandle** entity_set_handles,
	     /*inout*/int* entity_set_handles_allocated,
	     /*out*/int* entity_set_handles_size);

    /// Does the work for iMesh_getPrnts
    void
    getPrnts(/*in*/iBase_EntitySetHandle from_entity_set,
    /*in*/int num_hops,
	     /*inout*/iBase_EntitySetHandle** entity_set_handles,
	     /*inout*/int* entity_set_handles_allocated,
	     /*out*/int* entity_set_handles_size);

    /// Does the work for iMesh_subtract
    void
    subtract(/*in*/iBase_EntitySetHandle entity_set_1,
    /*in*/iBase_EntitySetHandle entity_set_2,
	     /*out*/iBase_EntitySetHandle* result_entity_set);

    /// Does the work for iMesh_intersect
    void
    intersect(/*in*/iBase_EntitySetHandle entity_set_1,
    /*in*/iBase_EntitySetHandle entity_set_2,
	      /*out*/iBase_EntitySetHandle* result_entity_set);

    /// Does the work for iMesh_unite
    void
    unite(/*in*/iBase_EntitySetHandle entity_set_1,
    /*in*/iBase_EntitySetHandle entity_set_2,
	  /*out*/iBase_EntitySetHandle* result_entity_set);

    /// Creates an entity set with a copy of the whole database in it.
    void
    vCopyRootSet(/*out*/EntitySetBase*& pESBResult,
    /*in*/const EntitySetBase* const pESBCompare = NULL);

    /// Returns an array with all the entity handles in it.
    void
    vCopyRootSet(/*out*/iBase_EntityHandle** entity_handles,
		 int* entity_handles_allocated, int* entity_handles_size);

    // TAGS TAGS TAGS TAGS TAGS TAGS TAGS TAGS TAGS TAGS TAGS TAGS
  private:
    /// Container for tag data.
    std::set<TagBase*> spTags;
  public:
    /**\brief Write all tag data.
     *
     * The args provide data to (a) find the index of a tagged entity
     * and (b) find the index of an entity tag value.
     */
    void
    writeTags(const char tagFileName[], const iBase_EntityHandle allHandles[],
	      const int allHandles_size,
	      const iBase_EntitySetHandle allSetHandles[],
	      const int allSetHandles_size);

    /// This function is used to check whether a tag handle is valid.
    bool
    isTag(void* theTag)
    {
      return (spTags.find(static_cast<TagBase*>(theTag)) != spTags.end());
    }

    /// Does the work for iMesh_createTag
    void
    createTag(/*in*/const char* tag_name,
    /*in*/const int tag_name_len,
	      /*in*/const int iTagSize,
	      /*in*/const int type,
	      /*out*/iBase_TagHandle* tagHandle);

    /// Does the work for iMesh_destroyTag
    void
    destroyTag(/*in*/iBase_TagHandle tagHandle, /*in*/bool forced);

    /// Does the work for iMesh_getTagName
    void
    getTagName(/*in*/iBase_TagHandle tagHandle, char* name, int name_len);

    /// Does the work for iMesh_getTagSizeValues
    int
    getTagSizeValues(/*in*/iBase_TagHandle tagHandle);

    /// Does the work for iMesh_getTagSizeBytes
    int
    getTagSizeBytes(/*in*/iBase_TagHandle tagHandle);

    /// Does the work for iMesh_getTagHandle
    iBase_TagHandle
    getTagHandle(/*in*/const char *name, const int name_len);

    /// Does the work for iMesh_getTagType
    int
    getTagType(/*in*/iBase_TagHandle tagHandle);

    /// Does the work for iMesh_getData
    void
    getData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle,
	    /*inout*/char** tag_value,
	    /*inout*/int* tag_value_allocated,
	    /*out*/int* tag_size);

    /// Does the work for iMesh_getIntData
    int
    getIntData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getDblData
    double
    getDblData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getEHData
    iBase_EntityHandle
    getEHData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getESHData
    iBase_EntitySetHandle
    getESHData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_setData
    void
    setData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle,
	    /*in*/const char* tag_value, /*in*/int tag_size);

    /// Does the work for iMesh_setIntData
    void
    setIntData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle,
	       /*in*/int tag_value);

    /// Does the work for iMesh_setDblData
    void
    setDblData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle,
	       /*in*/double tag_value);

    /// Does the work for iMesh_setEHData
    void
    setEHData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle,
	      /*in*/iBase_EntityHandle tag_value);

    /// Does the work for iMesh_setESHData
    void
    setESHData(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle,
	       /*in*/iBase_EntitySetHandle tag_value);

    /// Does the work for iMesh_getAllTags
    void
    getAllTags(/*in*/iBase_EntityHandle entity_handle,
    /*inout*/iBase_TagHandle** tag_handles,
	       /*inout*/int* tag_handles_allocated,
	       /*out*/int* tag_handles_size);

    /// Does the work for iMesh_rmvTag
    void
    rmvTag(/*in*/iBase_EntityHandle entity_handle,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getArrData
    void
    getArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
	       /*in*/iBase_TagHandle tag_handle,
	       /*inout*/char** tag_value,
	       /*inout*/int* tag_value_allocated,
	       /*out*/int* tag_value_size);

    /// Does the work for iMesh_getIntArrData
    void
    getIntArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		  /*in*/iBase_TagHandle tag_handle,
		  /*inout*/int** tag_value,
		  /*inout*/int* tag_value_allocated,
		  /*out*/int* tag_value_size);

    /// Does the work for iMesh_getDblArrData
    void
    getDblArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		  /*in*/iBase_TagHandle tag_handle,
		  /*inout*/double** tag_value,
		  /*inout*/int* tag_value_allocated,
		  /*out*/int* tag_value_size);

    /// Does the work for iMesh_getEHArrData
    void
    getEHArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		 /*in*/iBase_TagHandle tag_handle,
		 /*inout*/iBase_EntityHandle** tag_value,
		 /*inout*/int* tag_value_allocated,
		 /*out*/int* tag_value_size);

    /// Does the work for iMesh_getESHArrData
    void
    getESHArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		  /*in*/iBase_TagHandle tag_handle,
		  /*inout*/iBase_EntitySetHandle** tag_value,
		  /*inout*/int* tag_value_allocated,
		  /*out*/int* tag_value_size);

    /// Does the work for iMesh_setArrData
    void
    setArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
	       /*in*/iBase_TagHandle tag_handle,
	       /*in*/const char tag_values[],
	       /*in*/int tag_values_size);

    /// Does the work for iMesh_setIntArrData
    void
    setIntArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		  /*in*/iBase_TagHandle tag_handle,
		  /*in*/const int tag_values[],
		  /*in*/int tag_values_size);

    /// Does the work for iMesh_setDblArrData
    void
    setDblArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		  /*in*/iBase_TagHandle tag_handle,
		  /*in*/const double tag_values[],
		  /*in*/int tag_values_size);

    /// Does the work for iMesh_setEHArrData
    void
    setEHArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		 /*in*/iBase_TagHandle tag_handle,
		 /*in*/const iBase_EntityHandle tag_values[],
		 /*in*/int tag_values_size);

    /// Does the work for iMesh_setESHArrData
    void
    setESHArrData(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
		  /*in*/iBase_TagHandle tag_handle,
		  /*in*/const iBase_EntitySetHandle tag_values[],
		  /*in*/int tag_values_size);

    /// Does the work for iMesh_rmvArrTag
    void
    rmvArrTag(/*in*/const iBase_EntityHandle entity_handles[],
    /*in*/int entity_handles_size,
	      /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_setEntSetData
    void
    setEntSetData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle,
		  /*in*/const char* tag_value,
		  /*in*/const int tag_size);

    /// Does the work for iMesh_setEntSetIntData
    void
    setEntSetIntData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle,
		     /*in*/int tag_value);

    /// Does the work for iMesh_setEntSetDblData
    void
    setEntSetDblData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle,
		     /*in*/double tag_value);

    /// Does the work for iMesh_setEntSetEHData
    void
    setEntSetEHData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle,
		    /*in*/iBase_EntityHandle tag_value);

    /// Does the work for iMesh_setEntSetESHData
    void
    setEntSetESHData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle,
		     /*in*/iBase_EntitySetHandle tag_value);

    /// Does the work for iMesh_getEntSetData
    void
    getEntSetData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle,
		  /*inout*/char** tag_value,
		  /*inout*/int* tag_value_allocated,
		  /*out*/int* tag_size);

    /// Does the work for iMesh_getEntSetIntData
    int
    getEntSetIntData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getEntSetDblData
    double
    getEntSetDblData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getEntSetEHData
    iBase_EntityHandle
    getEntSetEHData(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getEntSetESHData
    iBase_EntitySetHandle
    getEntSetESHData(
    /*in*/iBase_EntitySetHandle entity_set,
		     /*in*/iBase_TagHandle tag_handle);

    /// Does the work for iMesh_getAllEntSetTags
    void
    getAllEntSetTags(/*in*/iBase_EntitySetHandle entity_set,
    /*inout*/iBase_TagHandle** tag_handles,
		     /*inout*/int* tag_handles_allocated,
		     /*out*/int* tag_handles_size);

    /// Does the work for iMesh_rmvEntSetTag
    void
    rmvEntSetTag(/*in*/iBase_EntitySetHandle entity_set,
    /*in*/iBase_TagHandle tag_handle);

  private:
    /// Helper function for finding tags by name.
    bool
    qFindTag(const char name[], const int name_len, TagBase*& ptbResult);

  public:
    // These pure virtual functions must be over-ridden by the volume,
    // plane, and surface classes.  Obviously.
    //
    // These are the only functions that differ between the three,
    // somewhat remarkably...

  public:
    /// Does the work for iMesh_load
    virtual void
    load(const char *name, const char *options, int* err) = 0;

    /// Does the work for iMesh_save
    virtual void
    save(const char *name, const char *options, int* err) = 0;

    /// Does the work for iMesh_getGeometricDimension
    virtual int
    getGeometricDimension() = 0;

    /// Does the work for iMesh_getAdjTable
    virtual void
    getAdjTable(int* adjacency_table) = 0;

    /// Does the work for iMesh_setAdjTable
    virtual void
    setAdjTable(int* adjacency_table, int* err) = 0;

    /// Does the work for iMesh_getEntities
    virtual void
    getEntities(const int entity_type, const int entity_topology,
		iBase_EntityHandle** entity_handles,
		int* entity_handles_allocated, int* entity_handles_size,
		int* err) = 0;

    /// Does the work for iMesh_getEntAdj
    virtual void
    getEntAdj(const iBase_EntityHandle entity_handle,
	      const int entity_type_requested,
	      iBase_EntityHandle** adj_entity_handles,
	      int* adj_entity_handles_allocated, int* adj_entity_handles_size,
	      int *err) = 0;

    /// Does the work for iMesh_getEntArrAdj
    virtual void
    getEntArrAdj(const iBase_EntityHandle* entity_handles,
		 const int entity_handles_size, const int entity_type_requested,
		 iBase_EntityHandle** adj_entity_handles,
		 int* adj_entity_handles_allocated,
		 int* adj_entity_handles_size, int** offset,
		 int* offset_allocated, int* offset_size, int *err) = 0;

    /// Does the work for iMesh_getEnt2ndAdj
    virtual void
    getEnt2ndAdj(const iBase_EntityHandle entity_handle,
		 const int bridge_entity_type, const int requested_entity_type,
		 iBase_EntityHandle** adjacent_entities,
		 int* adjacent_entities_allocated, int* adjacent_entities_size,
		 int* err) = 0;

    /// Does the work for iMesh_iBase_EntityHandle
    void
    getEntArr2ndAdj(iBase_EntityHandle const* entity_handles,
		    const int entity_handles_size, const int bridge_entity_type,
		    const int requested_entity_type,
		    iBase_EntityHandle** adj_entity_handles,
		    int* adj_entity_handles_allocated,
		    int* adj_entity_handles_size, int** offset,
		    int* offset_allocated, int* offset_size, int* err);

    /// Does the work for iMesh_createEnt
    virtual void
    createEnt(const int new_entity_topology,
	      const iBase_EntityHandle* lower_order_entity_handles,
	      const int lower_order_entity_handles_size,
	      iBase_EntityHandle* new_entity_handle, int* status, int *err) = 0;

    /// Does the work for iMesh_createEntArr
    virtual void
    createEntArr(const int new_entity_topology,
		 const iBase_EntityHandle* lower_order_entity_handles,
		 const int lower_order_entity_handles_size,
		 iBase_EntityHandle** new_entity_handles,
		 int* new_entity_handles_allocated,
		 int* new_entity_handles_size, int** status,
		 int* status_allocated, int* status_size, int *err) = 0;

    /// Does the work for iMesh_deleteEnt
    virtual void
    deleteEnt(iBase_EntityHandle entity_handle, int *err) = 0;
    // deleteEntArr calls deleteEnt, so it isn't needed here.

    /// Create the internal mesh data structure
    virtual void
    vCreateMesh()
    {
      assert(0);
    }
  public:
    /// True iff the given set exists for this mesh.
    bool
    qHasSet(EntitySetBase* ent_set)
    {
      return (spEntSets.find(ent_set) != spEntSets.end());
    }

    /// True iff the given array iterator exists for this mesh.
    bool
    qHasWorkset(Workset* ws)
    {
      return (spSetWorksets.find(ws) != spSetWorksets.end());
    }

    /// True iff the given single entity iterator exists for this mesh.
    bool
    qHasIter(EntIter* iter)
    {
      return (spSetIters.find(iter) != spSetIters.end());
    }

    /// Add an array iterator to the collection for this mesh.
    void
    addWorkset(Workset* ws)
    {
      spSetWorksets.insert(ws);
    }

    /// Add a single entity iterator to the collection for the mesh.
    void
    addEntIter(EntIter* iter)
    {
      spSetIters.insert(iter);
    }

    /// Add an entity set to the collection for the mesh.
    void
    addEntSet(EntitySetBase* set)
    {
      spEntSets.insert(set);
    }

    /// Remove an array iterator from the collection for this mesh.
    void
    removeWorkset(Workset* ws, int* err)
    {
      try {
	spSetWorksets.erase(ws);
      }
      catch (iBase_ErrorType& iET_) {
	*err = iET_;
	setError(*err);
	return;
      }
    }

    /// Remove a single entity iterator from the collection for the mesh.
    void
    removeEntIter(EntIter* iter, int* err)
    {
      try {
	spSetIters.erase(iter);
      }
      catch (iBase_ErrorType& iET_) {
	*err = iET_;
	setError(*err);
	return;
      }
    }

    /// Remove an entity set from the collection for the mesh.
    void
    removeEntSet(EntitySetBase* set)
    {
      spEntSets.erase(set);
    }
  };
// end iMesh_GRUMMP class

  /**\brief This class specializes iMesh_Base for 2D meshes.
   *
   * Functions that are inherited from the base class aren't documented
   * again here.
   */
  class iMesh_Plane : public iMesh_Base {
  public:
    /// Constructor called from iMesh_newMesh
    iMesh_Plane(const eMesh_Type eMT_in = ePlane) :
	iMesh_Base(eMT_in)
    {
    }

    /// Constructor that adds an iMesh wrapper around a Mesh2D object
    iMesh_Plane(Mesh2D *grummpMesh2D);

    /// Destructor
    ~iMesh_Plane()
    {
    }

    // Inherited from iMesh_Base
    virtual void
    vCreateMesh();
    virtual void
    load(const char *name, const char *options, int* err);
    virtual void
    save(const char *name, const char *options, int* err);

    virtual int
    getGeometricDimension()
    {
      return 2;
    }
    virtual void
    setAdjTable(int* adjacency_table, int* err);
    virtual void
    getAdjTable(int* adjacency_table);
    virtual void
    getEntities(const int entity_type, const int entity_topology,
		iBase_EntityHandle** entity_handles,
		int* entity_handles_allocated, int* entity_handles_size,
		int* err);
    virtual void
    getEntAdj(const iBase_EntityHandle entity_handle,
	      const int entity_type_requested,
	      iBase_EntityHandle** adj_entity_handles,
	      int* adj_entity_handles_allocated, int* adj_entity_handles_size,
	      int *err);
    virtual void
    getEntArrAdj(const iBase_EntityHandle* entity_handles,
		 const int entity_handles_size, const int entity_type_requested,
		 iBase_EntityHandle** adj_entity_handles,
		 int* adj_entity_handles_allocated,
		 int* adj_entity_handles_size, int** offset,
		 int* offset_allocated, int* offset_size, int *err);
    virtual void
    getEnt2ndAdj(const iBase_EntityHandle entity_handle,
		 const int bridge_entity_type, const int requested_entity_type,
		 iBase_EntityHandle** adjacent_entities,
		 int* adjacent_entities_allocated, int* adjacent_entities_size,
		 int* err);
    virtual void
    createEnt(const int new_entity_topology,
	      const iBase_EntityHandle* lower_order_entity_handles,
	      const int lower_order_entity_handles_size,
	      iBase_EntityHandle* new_entity_handle, int* status, int *err);
    virtual void
    createEntArr(const int new_entity_topology,
		 const iBase_EntityHandle* lower_order_entity_handles,
		 const int lower_order_entity_handles_size,
		 iBase_EntityHandle** new_entity_handles,
		 int* new_entity_handles_allocated,
		 int* new_entity_handles_size, int** status,
		 int* status_allocated, int* status_size, int *err);
    virtual void
    deleteEnt(iBase_EntityHandle entity_handle, int *err);
  };

  /**\brief This class specializes iMesh_Plane for surface meshes.
   *
   * This is done this way because nearly all of the specialized
   * functions for the 2D class can be reused.
   */
  class iMesh_Surface : public iMesh_Plane {
  public:
    iMesh_Surface() :
	iMesh_Plane(eSurface)
    {
    }
    /// Constructor that adds an iMesh wrapper around a SurfMesh object
    iMesh_Surface(SurfMesh *grummpSurfMesh);

    ~iMesh_Surface()
    {
    }

    virtual void
    vCreateMesh();
    virtual void
    load(const char *name, const char *options, int* err);
    virtual void
    save(const char *name, const char *options, int* err);
    virtual int
    getGeometricDimension()
    {
      return 3;
    }
  };

  /**\brief This class specializes iMesh_Base for 3D meshes.
   *
   * Functions that are inherited from the base class aren't documented
   * again here.
   */
  class iMesh_Volume : public iMesh_Base {
  public:
    /// Constructor
    iMesh_Volume() :
	iMesh_Base(eVolume)
    {
    }
    /// Constructor that adds an iMesh wrapper around a VolMesh object
    iMesh_Volume(VolMesh *grummpVolMesh);
    /// Destructor
    ~iMesh_Volume()
    {
    }

    // Inherited from iMesh_Base
    virtual void
    vCreateMesh();
    virtual void
    load(const char *name, const char *options, int* err);
    virtual void
    save(const char *name, const char *options, int* err);

    virtual int
    getGeometricDimension()
    {
      return 3;
    }
    virtual void
    setAdjTable(int* adjacency_table, int* err);
    virtual void
    getAdjTable(int* adjacency_table);
    virtual void
    getEntities(const int entity_type, const int entity_topology,
		iBase_EntityHandle** entity_handles,
		int* entity_handles_allocated, int* entity_handles_size,
		int* err);
    virtual void
    getEntAdj(const iBase_EntityHandle entity_handle,
	      const int entity_type_requested,
	      iBase_EntityHandle** adj_entity_handles,
	      int* adj_entity_handles_allocated, int* adj_entity_handles_size,
	      int *err);
    virtual void
    getEntArrAdj(const iBase_EntityHandle* entity_handles,
		 const int entity_handles_size, const int entity_type_requested,
		 iBase_EntityHandle** adj_entity_handles,
		 int* adj_entity_handles_allocated,
		 int* adj_entity_handles_size, int** offset,
		 int* offset_allocated, int* offset_size, int *err);
    virtual void
    getEnt2ndAdj(const iBase_EntityHandle entity_handle,
		 const int bridge_entity_type, const int requested_entity_type,
		 iBase_EntityHandle** adjacent_entities,
		 int* adjacent_entities_allocated, int* adjacent_entities_size,
		 int* err);
    virtual void
    createEnt(const int new_entity_topology,
	      const iBase_EntityHandle* lower_order_entity_handles,
	      const int lower_order_entity_handles_size,
	      iBase_EntityHandle* new_entity_handle, int* status, int *err);
    virtual void
    createEntArr(const int new_entity_topology,
		 const iBase_EntityHandle* lower_order_entity_handles,
		 const int lower_order_entity_handles_size,
		 iBase_EntityHandle** new_entity_handles,
		 int* new_entity_handles_allocated,
		 int* new_entity_handles_size, int** status,
		 int* status_allocated, int* status_size, int *err);
    virtual void
    deleteEnt(iBase_EntityHandle entity_handle, int *err);
  };

}

/**\brief iMesh_Instance_Private is a wrapper class around a
 * iMesh_Base*.
 *
 * At some point, I tried to remove this extra wrapper, but it didn't
 * turn out to be as easy as I'd have thought...
 *
 * The class has C binding so that a pointer to it can be typedef'ed to
 * iMesh_Instance in the iMesh.h header file.
 */

extern "C"
{
  class iMesh_Instance_Private : public iBase_Taggable {
  private:
    iMesh_Instance_Private(iMesh_Instance_Private&);
    iMesh_Instance_Private&
    operator=(const iMesh_Instance_Private&);
    iMesh_Base::eMesh_Type m_eMT;
  public:
    /// The underlying data that runs all this.
    iMesh_Base *pGB;

    /// Constructor, which initializes the iMesh_Base data.
    iMesh_Instance_Private(const enum iMesh_Base::eMesh_Type eMT) :
	iBase_Taggable(), m_eMT(eMT), pGB(NULL)
    {
    }

    iMesh_Instance_Private(iMesh_Base * const pIMB) :
	iBase_Taggable(), m_eMT(iMesh_Base::eDefault), pGB(pIMB)
    {
      Mesh::eMeshType GMT = pIMB->pMesh()->getType();
      switch (GMT)
	{
	case Mesh::eMesh2D:
	  m_eMT = iMesh_Base::ePlane;
	  break;
	case Mesh::eSurfMesh:
	  m_eMT = iMesh_Base::eSurface;
	  break;
	case Mesh::eVolMesh:
	  m_eMT = iMesh_Base::eVolume;
	  break;
	}
    }

    /// Destructor
    ~iMesh_Instance_Private()
    {
      if (pGB)
	delete pGB;
    }

    /// Create an iMesh_Base, using the factory.
    void
    vInitiMesh(const enum iMesh_Base::eMesh_Type eMT = iMesh_Base::eDefault)
    {
      if (eMT == m_eMT && pGB)
	return;

      if (eMT == iMesh_Base::eDefault) {
	// Do nothing
      }
      else {
	// Can't reliably change at this point, apparently; try anyway.
	if (pGB)
	  delete pGB;
	m_eMT = eMT;
      }
      pGB = iMesh_Base::Create(m_eMT);
      pGB->vCreateMesh();
    }
    int
    getError() const
    {
      return pGB->getError();
    }
  };

} // end extern "C"

#endif
