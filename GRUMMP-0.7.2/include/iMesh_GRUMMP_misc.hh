#ifndef iMesh_GRUMMP_misc_hh
#define iMesh_GRUMMP_misc_hh
#include "stdlib.h"
#include "iBase.h"

#include <algorithm>

// The do-while(0) business below is so that these macros can be used as
// statements, even in un-{}'ed if clauses or loops.

/// Make sure ITAPS arrays have the right size, or die trying.
template<class T>
  void
  vEnsureArraySize(T** array, int* alloc_size, const int req_size)

  {
    if (*alloc_size == 0) {
      if (req_size > 0) {
	*array = reinterpret_cast<T*>(calloc(req_size, sizeof(T)));
	*alloc_size = req_size;
      }
    }
    else if (*alloc_size < req_size) {
      throw(iBase_BAD_ARRAY_SIZE);
    }
  }

// Some of the internal functions throw exceptions; this macro
// converts those exceptions into iMesh error return values.
#define TRY(a)					\
  do {						\
    try {					\
      a;					\
    } catch (iBase_ErrorType& iET_) {		\
      *err = iET_;				\
      instance->pGB->setError(*err);		\
      return;					\
    }						\
  } while (0)

// This is a variant on the TRY macro, specialized for the array size
// check above.
#define TRY_ARRAY_SIZE(array_, alloc_size_, req_size_, type_)	\
  do {									\
    try {								\
      vEnsureArraySize<type_>(array_, alloc_size_, req_size_);		\
    } catch (iBase_ErrorType& ET_) {					\
      *err = ET_;							\
      setError(*err);							\
      return;								\
    }									\
    catch (...) {							\
      *err = iBase_FAILURE;						\
      setError(*err);							\
      return;								\
    }									\
  } while (0)

// This is a variant on the TRY macro, specialized for the array size
// check above.
#define TRY_ARRAY_SIZE2(array_, alloc_size_, req_size_, type_)	\
  do {									\
    try {								\
      vEnsureArraySize<type_>(array_, alloc_size_, req_size_);		\
    } catch (iBase_ErrorType& ET_) {					\
      *err = ET_;							\
      instance->pGB->setError(*err);					\
      return;								\
    }									\
    catch (...) {							\
      *err = iBase_FAILURE;						\
      instance->pGB->setError(*err);					\
      return;								\
    }									\
  } while (0)

#define GR_rootSet (reinterpret_cast<iBase_EntitySetHandle>(INT_MIN + 1))

// #ifdef NDEBUG

// // For production code, there are a bunch of checks that get compiled
// // out; a decent optimizer will make this all go away.
// #define CHECK_ENTITY(pMG_, entity_, err_) do {} while (0)
// #define CHECK_ENTITY_SET(pMG_, entity_set_, err_) do {} while (0)
// #define CHECK_TAG(pMG_, ptbTheTag_, err_) do {} while (0)
// #define CHECK_ENT_ITER(pMG_, pEI_, err_) do {} while (0)
// #define CHECK_ENTARR_ITER(pMG_, pWs_, err_) do {} while (0)
// #define CHECK_CYCLIC(pESB_parent_, pESB_child_, err_) do {} while (0)
// #define CHECK_FOR_NULL_TERMINATION(string_, len_, err_)	do {} while (0)

// #else

// Check whether an entity really does exist in this iMesh database.
#define CHECK_ENTITY(pMG_, entity_, err_)				\
  do {									\
    if ( entity_ == NULL || !pMG_->pGB || !pMG_->pGB->pMesh() ||	\
	 pMG_->pGB->pMesh()->getNumVerts() == 0 ||			\
	 !pMG_->pGB->pMesh()->isValidEntHandle(entity_) ) {		\
      *err_ = iBase_INVALID_ENTITY_HANDLE;				\
      pMG_->pGB->setError(*err);					\
      return;								\
    }									\
  } while (0)

// Check whether an entity set really does exist in this iMesh database.
#define CHECK_ENTITY_SET(pMG_, entity_set_, err_)			\
  do {									\
    void *pvSet = reinterpret_cast<void*>(entity_set_);		\
    if ( (pvSet != GR_rootSet) &&					\
	 !pMG_->pGB->qHasSet						\
	 (dynamic_cast<ITAPS_GRUMMP::EntitySetBase*>(entity_set_))) {	\
      *err_ = iBase_INVALID_ENTITYSET_HANDLE;				\
      pMG_->pGB->setError(*err);					\
      return;								\
    }									\
  } while (0)

// Check whether a tag really does exist in this iMesh database.
#define CHECK_TAG(pMG_, ptbTheTag_, err_)				\
  do {									\
    if ( !pMG_->pGB->isTag(ptbTheTag_) ) {				\
      *err_ = iBase_INVALID_TAG_HANDLE;					\
      pMG_->pGB->setError(*err);					\
      return;								\
    }									\
} while (0)

// Check whether an iterator really does exist in this iMesh database.
#define CHECK_ENT_ITER(pMG_, pEI_, err_)				\
  do {									\
    if ( !pMG_->pGB->qHasIter(pEI_)) {					\
      *err_ = iBase_INVALID_ITERATOR_HANDLE;				\
      pMG_->pGB->setError(*err);					\
      return;								\
    }									\
  } while (0)

// Check whether an array iterator really does exist in this iMesh database.
#define CHECK_ENTARR_ITER(pMG_, pWs_, err_)				\
  do {									\
    if ( !pMG_->pGB->qHasWorkset(pWs_))	{				\
      *err_ = iBase_INVALID_ITERATOR_HANDLE;				\
      pMG_->pGB->setError(*err);					\
      return;								\
    }									\
  } while (0)

// Check whether a string is null terminated, as the C spec requires.
#define CHECK_FOR_NULL_TERMINATION(string_, len_, err_)			\
  do {									\
    if (string_[len_] != (char)(0)) {					\
      *err_ = iBase_FAILURE;						\
      pMG_->pGB->setError(*err);					\
      return;								\
    }									\
  } while (0)

// #endif

// Confirm that the provided mesh type and topology are compatible.
#define CHECK_TYPE_AND_TOPO(type_, topo_, err_)				\
  do {									\
    if ( type_ < iBase_VERTEX || type_ > iBase_ALL_TYPES ) {		\
      *err_ = iBase_INVALID_ENTITY_TYPE;				\
      instance->pGB->setError(*err);					\
      return;								\
    }									\
    else if ( topo_ < iMesh_POINT || topo_ > iMesh_ALL_TOPOLOGIES ) {	\
      *err_ = iBase_INVALID_ENTITY_TOPOLOGY;				\
      instance->pGB->setError(*err);					\
      return;								\
    }									\
    else {								\
      switch (type_) {							\
      case iBase_VERTEX:						\
	if (topo_ != iMesh_POINT && topo_ != iMesh_ALL_TOPOLOGIES) {	\
	  *err_ = iBase_BAD_TYPE_AND_TOPO;				\
	  instance->pGB->setError(*err);				\
	  return;							\
	}								\
	break;								\
      case iBase_EDGE:							\
	if (topo_ != iMesh_LINE_SEGMENT &&				\
	    topo_ != iMesh_ALL_TOPOLOGIES) {				\
	  *err_ = iBase_BAD_TYPE_AND_TOPO;				\
	  instance->pGB->setError(*err);				\
	  return;							\
	}								\
	break;								\
      case iBase_FACE:							\
	if (topo_ != iMesh_TRIANGLE &&					\
	    topo_ != iMesh_QUADRILATERAL &&				\
	    topo_ != iMesh_POLYGON &&					\
	    topo_ != iMesh_ALL_TOPOLOGIES) {				\
	  *err_ = iBase_BAD_TYPE_AND_TOPO;				\
	  instance->pGB->setError(*err);				\
	  return;							\
	}								\
	break;								\
      case iBase_REGION:						\
	if (topo_ != iMesh_TETRAHEDRON &&				\
	    topo_ != iMesh_PYRAMID &&					\
	    topo_ != iMesh_PRISM &&					\
	    topo_ != iMesh_HEXAHEDRON &&				\
	    topo_ != iMesh_SEPTAHEDRON &&				\
	    topo_ != iMesh_POLYHEDRON &&				\
	    topo_ != iMesh_ALL_TOPOLOGIES) {				\
	  *err_ = iBase_BAD_TYPE_AND_TOPO;				\
	  instance->pGB->setError(*err);				\
	  return;							\
	}								\
	break;								\
      case iBase_ALL_TYPES:						\
	break;								\
      }									\
    }									\
  } while (0)

// Ensure that a string will be null terminated, even if it wasn't before.
#define NULL_TERMINATE(a_in_, len_, a_)		\
  do {						\
    snprintf(a_, 1023, "%s", a_in_);		\
    a_[len_] = 0;				\
    int i__ = len_-1;				\
    while (i__ >= 0 && a_[i__] == ' ') {	\
      a_[i__--] = 0;				\
    };						\
  } while(false);

// Find the location of data in the array [start, end); it is an error
// for data -not- to be present.
template<typename T>
  int
  findArrayIndex(const T& data, const T* const start, const T* const end)
  {
    T const *loc = std::find(start, end, data);
    int index = loc - start;
    assert(index >= 0 && index < end - start);
    return index;
  }

#endif
