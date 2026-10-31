#ifndef GR_Entity_h
#define GR_Entity_h

// The following class is defined primarily so that TSTT entities will
// work smoothly.  In particular, with a base Entity class that all
// GRUMMP mesh pieces derive from, in the TSTT interface (a) many entity
// queries can be handled trivially by polymorphism (perhaps easier than
// in a class that has internal state, interestingly) and (b) more
// importantly, when a function requires the return of a collection of
// entity handles, I can return voided pointers to GRUMMP objects, and
// never have to address management of memory for newly-minted entities.

#include "GR_config.h"
#include "GR_Classes.h"

#ifdef ITAPS
#include "GR_iMesh_Classes.h"
typedef iBase_EntityHandle_Private GRUMMP_Entity;
class Entity : public iBase_EntityHandle_Private {
#else

// If we're not building for ITAPS, then there's no iBase_EntityHandle_Private struct around to inherit from.
  typedef Entity GRUMMP_Entity;
  class Entity {

#endif
public:
  typedef enum {
    eMin, eInterior = eMin, eGhosted, ePartBdry, eForeign, eUnknown, eMax
  } Entity_Status_t;
private:
  Entity_Status_t m_stat :7;
  bool m_owned :1;
public:
  Entity() :
      m_stat(eInterior), m_owned(true)
  {
  }
  virtual
  ~Entity()
  {
  }
  /// Ensure a sane initial state for face data
  virtual void
  resetAllData() = 0;
  virtual const Vert *
  getVert(const int iV) const = 0;
  virtual const Face *
  getFace(const int iF) const = 0;
  virtual const Cell *
  getCell(const int iC) const = 0;
  Vert *
  getVert(const int iV)
  {
    return const_cast<Vert*>(static_cast<const Entity*>(this)->getVert(iV));
  }
  Face *
  getFace(const int iF)
  {
    return const_cast<Face*>(static_cast<const Entity*>(this)->getFace(iF));
  }
  Cell *
  getCell(const int iC)
  {
    return const_cast<Cell*>(static_cast<const Entity*>(this)->getCell(iC));
  }

  virtual void
  getAllVertHandles(GRUMMP_Entity* aHandle[]) const = 0;
  virtual void
  getAllFaceHandles(GRUMMP_Entity* aHandle[]) const = 0;
  virtual void
  getAllCellHandles(GRUMMP_Entity* aHandle[]) const = 0;

  virtual int
  getNumVerts() const = 0;
  virtual int
  getNumCells() const = 0;
  virtual int
  getNumFaces() const = 0;
  virtual int
  getEntType() const = 0;
  virtual int
  getEntTopology() const = 0;
  virtual void
  markAsDeleted() = 0;
  virtual bool
  isDeleted() const = 0;

  bool
  isOwned() const
  {
    return m_owned;
  }
  bool
  setOwnership(const bool newVal)
  {
    bool oldVal = m_owned;
    m_owned = newVal;
    return oldVal;
  }

  Entity_Status_t
  getEntStatus() const
  {
    return m_stat;
  }
  Entity_Status_t
  setEntStatus(const Entity_Status_t newVal)
  {
    Entity_Status_t oldVal = m_stat;
    m_stat = newVal;
    return oldVal;
  }
};

#endif // GR_Entity_h
