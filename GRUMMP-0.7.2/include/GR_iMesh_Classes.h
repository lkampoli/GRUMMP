#ifndef GR_iMesh_Classes_h
#define GR_iMesh_Classes_h

#include "iBase.h"
#include "iMesh.h"

/**\brief This file contains forward declarations of all iMesh classes.
 * Only those that are actually declared here as base classes are
 * documented in this file.
 */

namespace ITAPS_GRUMMP
{
/// A base class for all things that iMesh can apply tags to.
  class iBase_Taggable {
  public:
    iBase_Taggable()
    {
    }
    virtual
    ~iBase_Taggable()
    {
    }
  };

  class iMesh_Base;
  class iMesh_Plane;
  class iMesh_Surface;
  class iMesh_Volume;

  class EntIter;
  class MeshIter;
  class EntitySetIter;
  class EntityListIter;

  class Workset;
  class MeshWorkset;
  class EntitySetWorkset;
  class EntityListWorkset;

  class VertRec;

  class EntitySetBase;
  class EntitySet;
  class EntityList;

  class TagBase;
}
;

/// Forward declaration of a class for tags.
class iBase_TagHandle_Private {
public:
  virtual
  ~iBase_TagHandle_Private()
  {
  }
};

/// Forward declaration of a class for tags.
class iBase_EntityHandle_Private : public ITAPS_GRUMMP::iBase_Taggable {
public:
  virtual
  ~iBase_EntityHandle_Private()
  {
  }
};

/// A base class for all the RefImpl entity sets.
class iBase_EntitySetHandle_Private : public ITAPS_GRUMMP::iBase_Taggable {
public:
  iBase_EntitySetHandle_Private()
  {
  }
  virtual
  ~iBase_EntitySetHandle_Private()
  {
  }
};

class iBase_EntityArrIterator_Private {
public:
  virtual
  ~iBase_EntityArrIterator_Private()
  {
  }
};

class iBase_EntityIterator_Private {
public:
  virtual
  ~iBase_EntityIterator_Private()
  {
  }
};

class iMesh_Instance_Private;

#endif
