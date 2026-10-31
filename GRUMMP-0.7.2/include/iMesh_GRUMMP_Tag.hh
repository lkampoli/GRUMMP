// File:          "iMesh_GRUMMP_Tag.hh"
// Modified by:    Ben Kerby
// First Modified: May 2, 2005

// The tag class enables entities to be tagged with a variety of kinds of data.
// It utilizes a hash map that maps from entity handles to value the
// corresponding value.  Each entity handle is passed in as a void* and is
// assumed to be a pointer of some kind, but is not checked for any kind of
// validity.  The value for each entity is unique (i.e. if another value is
// passed in with the same entity handle, it will be overwritten).

// The data incoming can be any number of values of a standard type.  The
// number of values is passed in as an argument to the map constructor.  The
// map is responsible for copying any incoming data into its own internal
// array.

// Collision resolution is handled through chaining, so each item in the map
// array is the head of a linked list and that is why they are all initialized
// to NULL.

#ifndef iMesh_GRUMMP_Tag_hh
#define iMesh_GRUMMP_Tag_hh

#include <set>
#include <map>
#include <string>
#include <algorithm>
#include "GR_iMesh_Classes.h"
#include <string.h>

namespace ITAPS_GRUMMP
{

/// Max average number of entries per bin in the hash map.
  static const float MAX_ALPHA = 3.0;

/// Min average number of entries per bin in the hash map.
  static const float MIN_ALPHA = 0.5;

/// Min size of the initial hash map (can be overridden at construction time).
  static const unsigned int MIN_MAP_SIZE = 256;

/// Ratio by which the size of the hash map grows or shrinks.
  static const unsigned int GROWTH_FACTOR = 2;

  /**\class HashMap implements an STL map using hash tables.
   *
   * This hash map is used in RefImpl to store tags with faster
   * retrieval than a std::map.  std::hash_map is officially an SGI
   * extension, so we rolled our own.
   *
   * Any given key is hashed to an integer, which is an index into an
   * array of linked lists.  The links in these lists each store the
   * data (of class Value) for one key, and a pointer to the next entry
   * in the linked list.  Expected average retrieval time is
   * O(MAX_ALPHA) + the cost of the hash function.  This compares very
   * favorably to O(log n) for std::map, where n is the total number of
   * entities with this tag.
   *
   * A number of constants for changing the behaviour of the map (to
   * optimize it, if necessary).  The map has an initial size of
   * MIN_MAP_SIZE, unless specified othewise as a constructor parameter.
   * When the number of items in the array exceeds MAX_ALPHA times the
   * current size of the array, the array is increased in size by the
   * GROWTH_FACTOR.  Similarly, when the  number of items goes below
   * MIN_ALPHA times the size of the array, the size of the array is
   * decreased by the GROWTH_FACTOR.  It will not, however, go below the
   * MIN_MAP_SIZE. Because of the way the hash function works, both
   * MIN_MAP_SIZE and GROWTH_FACTOR should be powers of 2.
   */
  template<class Key, class Value>
    class HashMap {
      HashMap(const HashMap&);
      HashMap&
      operator=(const HashMap&);
      /**\brief The number of class Value items that are stored in a given tag.
       *
       * iNumValues must not change during the lifetime of the map
       * or very bad things will happen.
       */
      const int iNumValues;

      /**\brief Size of the array into which hashed data is stored.
       *
       * This can change if the hash table gets too full.
       */
      unsigned int uiArrSize;

      /// Number of keys in the hash table (objects with this tag).
      unsigned int uiNumKeys;

      /// The basic unit of each linked list, a link
      struct Link {
	Key linkKey;
	Link* plNextLink;
	// Has an unsized array at the end
	Value aLinkValue[];
      };

      /// Array in which tag data is actually stored.
      Link** aplMapArray;

    public:
      /// Create a hash map with a given data size and initial hash table size.
      HashMap(const int iValuesIn, unsigned int uiInitSize = MIN_MAP_SIZE) :
	  iNumValues(iValuesIn), uiArrSize(uiInitSize), uiNumKeys(0), aplMapArray(
	  NULL)
      {
	// Create a new array of heads of linked lists
	aplMapArray = new Link*[uiInitSize];
	// Ensure that all of the linked lists are marked as empty
	std::fill(aplMapArray, aplMapArray + uiArrSize,
		  static_cast<Link*>(NULL));
      }

      /// Destroy the hash map.
      ~HashMap()
      {
	Link *plCursor, *plPrev;
	// For each entry in the array...
	for (unsigned int i = 0; i < uiArrSize; i++)
	  // ... traverse the linked list...
	  for (plCursor = aplMapArray[i]; plCursor != NULL;) {
	    plPrev = plCursor;
	    plCursor = plCursor->plNextLink;
	    // ... deleting as you go
	    deleteLink(plPrev);
	  }
	delete[] aplMapArray;
      }

      /// Returns true iff the hash map is empty.
      bool
      qEmpty() const
      {
	// Pre:  NONE
	// Post: If the map is empty, true is returned; otherwise false
	//       is returned
	return (uiNumKeys == 0);
      }

      /// Returns true iff this key is already in the map.
      bool
      qFind(const Key& keyIn) const
      // Pre:  NONE
      // Post: If a key-value pair with the same key has been added
      //       to the map, returns true, else returns false
      {
	unsigned int uiIndex = uiHash(keyIn);
	// Check the linked list at the index obtained
	Link* plCursor = aplMapArray[uiIndex];
	while (plCursor != NULL) {
	  if (plCursor->linkKey == keyIn)
	    return true;
	  plCursor = plCursor->plNextLink;
	}
	return false;
      }

      /// Adds this key to the map, or overwrites pre-existing data.
      void
      vInsert(const Key& keyIn, const Value aValueIn[])
      // Pre:  NONE
      // Post: Inserts the key-value pair given into the map.  If a
      //       pair with the same key already exists, overwrites the
      //       existing value with the one desired
      {
	// Get the hash index
	unsigned int uiIndex = uiHash(keyIn);
	// A couple of special cases to deal with first:
	//    - there are no items at that index
	if (aplMapArray[uiIndex] == NULL)
	  aplMapArray[uiIndex] = newLink(iNumValues, keyIn, aValueIn);
	//    - the first item has the same key
	else if (aplMapArray[uiIndex]->linkKey == keyIn) {
	  for (int i = 0; i < iNumValues; i++)
	    aplMapArray[uiIndex]->aLinkValue[i] = aValueIn[i];
	  // Return to prevent incrementing the number of keys
	  return;
	}
	else {
	  // Traverse the linked list
	  Link *plCursor = aplMapArray[uiIndex];
	  // Until we reach the last link
	  while (plCursor->plNextLink != NULL) {
	    // If the next link has the same key...
	    if (plCursor->plNextLink->linkKey == keyIn) {
	      plCursor = plCursor->plNextLink;
	      // ... overwrite the values...
	      for (int i = 0; i < iNumValues; i++)
		plCursor->aLinkValue[i] = aValueIn[i];
	      // ... and quit to prevent incrementing the number of keys
	      return;
	    }
	    // Move on to the next link
	    plCursor = plCursor->plNextLink;
	  }
	  // If we get here, plCursor has the last link in the chain, so
	  // create a new link and slap it on the end
	  plCursor->plNextLink = newLink(iNumValues, keyIn, aValueIn);
	}
	// Increment the number of keys and check to see if growth is necessary
	uiNumKeys++;
	if (uiNumKeys > MAX_ALPHA * uiArrSize)
	  vGrow();
      }

      /// Retrieves data for the given key.
      bool
      qGet(const Key& keyIn, Value aValueIn[]) const
      // Pre:  aValueIn is an array with a size of at least iNumValues
      // Post: If a pair with the same key is stored in the map, aValueIn
      //       is written with the value stored therein and true is returned.
      //       If the key is not found, false is returned (and aValueIn is
      //       left untouched).
      {
	// Get the index of the linked list to look in
	unsigned int uiIndex = uiHash(keyIn);
	// Traverse the linked list at that index
	Link* plCursor = aplMapArray[uiIndex];
	while (plCursor != NULL) {
	  // If the desired key is found...
	  if (plCursor->linkKey == keyIn) {
	    // ... copy the value(s) and return true
	    for (int i = 0; i < iNumValues; i++)
	      aValueIn[i] = plCursor->aLinkValue[i];
	    return true;
	  }
	  plCursor = plCursor->plNextLink;
	}
	// If we get here the key was not found, so return false
	return false;
      }

      /// Removes all entries with the given data.
      void
      vEraseByValue(const Value valueIn[/*iNumValues*/])
      // Pre:  NONE
      // Post: If one or more pairs with the given value exists in the map,
      //       they are removed.
      {
	// Have to traverse the entire map, because it could be anywhere.
	Link *plCursor, *plPrev;
	// For each entry in the array...
	for (unsigned int i = 0; i < uiArrSize; i++) {
	  // ... traverse the linked list...
	  for (plCursor = aplMapArray[i]; plCursor != NULL;) {
	    plPrev = plCursor;
	    plCursor = plCursor->plNextLink;
	    bool qMatch = true;
	    for (int ii = 0; ii < iNumValues && qMatch; ii++) {
	      qMatch = (valueIn[ii] == plCursor->aLinkValue[ii]);
	    }
	    if (qMatch) {
	      deleteLink(plCursor);
	      plCursor = plPrev->plNextLink;
	    }
	  }
	}
      }

      /// Removes the entry with the given key, if it exists.
      bool
      qErase(const Key& keyIn)
      // Pre:  NONE
      // Post: If a pair with the given key exists in the map, it is removed
      //       and true is returned; otherwise, false is returned
      {
	// Get the index to start at
	unsigned int uiIndex = uiHash(keyIn);
	// Two special cases to deal with:
	//    - the linked list is empty
	if (aplMapArray[uiIndex] == NULL)
	  return false;
	//    - the first item in the linked list is the one to delete
	if (aplMapArray[uiIndex]->linkKey == keyIn) {
	  Link* temp = aplMapArray[uiIndex]->plNextLink;
	  deleteLink(aplMapArray[uiIndex]);
	  aplMapArray[uiIndex] = temp;
	  uiNumKeys--;
	  if (uiNumKeys < MIN_ALPHA * uiArrSize && uiArrSize > MIN_MAP_SIZE)
	    vShrink();
	  return true;
	}
	// If we get here, we know that there is a first item in the list and it is
	// not the one we're looking for
	Link *plCursor, *plPrev;
	plCursor = aplMapArray[uiIndex]->plNextLink;
	plPrev = aplMapArray[uiIndex];
	while (plCursor != NULL) {
	  // If we find the key...
	  if (plCursor->linkKey == keyIn) {
	    // ... reroute the pointers around that link...
	    plPrev->plNextLink = plCursor->plNextLink;
	    // ... and deallocate the memory
	    deleteLink(plCursor);
	    uiNumKeys--;
	    // Shrink if necessary, but not below the MIN_MAP_SIZE
	    if (uiNumKeys < MIN_ALPHA * uiArrSize && uiArrSize > MIN_MAP_SIZE)
	      vShrink();
	    return true;
	  }
	  // Continue the traversal (key was not found at this link)
	  plPrev = plCursor;
	  plCursor = plCursor->plNextLink;
	}
	// If we get here, the key was not found, so return false
	return false;
      }

    private:

      /**\brief Create new link, and optionally add to the head of a linked list.
       *
       * Creates a new Link with the speicifed size, key, value, and pointer
       * to the next link (optional).  'malloc' is used instead of the 'new'
       * operator in order to reserve space for the unsized array that is
       * part of the Link struct.
       */
      Link*
      newLink(int iSize, const Key& keyIn, const Value aValueIn[], Link* next =
      NULL)
      {
	void* newMem = malloc(sizeof(Link) + iSize * sizeof(Value));
	if (newMem == NULL)
	  throw std::bad_alloc();
	Link* pLNew = reinterpret_cast<Link*>(newMem);
	pLNew->linkKey = keyIn;
	pLNew->plNextLink = next;
	for (int i = 0; i < iSize; i++)
	  pLNew->aLinkValue[i] = aValueIn[i];
	return pLNew;
      }

      /// Delete a link (already disconnected from the linked list).
      inline void
      deleteLink(Link* theTarget)
      // Must use 'free' to delete links because the memory was allocated
      // using 'malloc'
      {
	free(theTarget);
      }

      /// Essentially stirs the bits for an unsigned integer uiTarget
      inline void
      vMix(unsigned int& uiTarget) const
      {
	uiTarget += uiTarget << 13;
	uiTarget ^= uiTarget >> 7;
	uiTarget += uiTarget << 3;
	uiTarget ^= uiTarget >> 17;
	uiTarget += uiTarget << 5;
      }

      /**\brief Hash a Key into the range [0,uiArrSize).
       *
       * Does not require tables of random integers.  This works for both
       * 32- and 64-bit pointers.  No serious testing has been done to
       * confirm that this hash function produces a uniform hash of
       * representative pointer values.
       */
      unsigned int
      uiHash(const Key& keyIn) const
      {
	unsigned int iInput = 0;
	if (sizeof(Key) <= sizeof(int)) {
	  memcpy(&iInput, &keyIn, sizeof(int));
	}
	else {
	  assert(sizeof(Key) <= 2 * sizeof(int));
	  int iHi, iLo;
	  memcpy(&iHi, &keyIn, sizeof(int));
	  memcpy(&iLo, (reinterpret_cast<const char*>(&keyIn)) + sizeof(int),
		 sizeof(int));
	  iInput = iHi + iLo;
	}
	// Split the data into individual bytes
	const unsigned int aiData[] =
	  { iInput / (256 * 65536), (iInput / 65536) % 256, (iInput / 256)
	      % 256, iInput % 256 };
	// A large prime to mess things up
	const unsigned int p = 16777619;

	// Some random sort of bits to start with
	unsigned int uiHashed = 0x811C9DC5;

	// For each byte, XOR it with the hashed value and mix it up
	for (int i = 0; i < 4; i++) {
	  uiHashed = (uiHashed ^ aiData[i]) * p;
	  vMix(uiHashed);
	}

	// Truncate it to the size required
	return uiHashed % uiArrSize;
      }

      /// A helper function, used to move the links in the chain to their
      /// new location for the grow and shrink functions (and only
      /// those). 
      void
      vAdd(Link* pLink)
      {
	unsigned int uiIndex = uiHash(pLink->linkKey);
	pLink->plNextLink = NULL;
	if (aplMapArray[uiIndex] == NULL) {
	  aplMapArray[uiIndex] = pLink;
	}
	else {
	  Link* plCursor = aplMapArray[uiIndex];
	  while (plCursor->plNextLink != NULL)
	    plCursor = plCursor->plNextLink;
	  plCursor->plNextLink = pLink;
	}
      }

      /**\brief Increases the size of the array, rehashing all the elements.
       * If the memory allocation of the new array fails, it attempts to
       * recover the original array.
       */
      void
      vGrow()
      {
	Link** aplOldArray = aplMapArray;
	unsigned int uiOldSize = uiArrSize;
	uiArrSize *= GROWTH_FACTOR;
	try {
	  aplMapArray = new Link*[uiArrSize];
	}
	catch (std::bad_alloc& b) {
	  aplMapArray = aplOldArray;
	  uiArrSize = uiOldSize;
	  throw(b);
	}
	std::fill(aplMapArray, aplMapArray + uiArrSize,
		  static_cast<Link*>(NULL));
	Link *plOldCursor, *plTemp;
	for (unsigned int uiOldIndex = 0; uiOldIndex < uiOldSize;
	    uiOldIndex++) {
	  plOldCursor = aplOldArray[uiOldIndex];
	  while (plOldCursor != NULL) {
	    plTemp = plOldCursor->plNextLink;
	    vAdd(plOldCursor);
	    plOldCursor = plTemp;
	  }
	}
	delete[] aplOldArray;
      }

      /**\brief Decreases the size of the array, rehashing all the elements.
       *
       * If the memory allocation of the new array fails, it attempts
       * to recover the original array.
       */
      void
      vShrink()
      {
	Link** aplOldArray = aplMapArray;
	unsigned int uiOldSize = uiArrSize;
	uiArrSize /= GROWTH_FACTOR;
	try {
	  aplMapArray = new Link*[uiArrSize];
	}
	catch (std::bad_alloc& b) {
	  aplMapArray = aplOldArray;
	  uiArrSize = uiOldSize;
	  throw(b);
	}
	std::fill(aplMapArray, aplMapArray + uiArrSize,
		  static_cast<Link*>(NULL));
	Link *plOldCursor, *plTemp;
	for (unsigned int uiOldIndex = 0; uiOldIndex < uiOldSize;
	    uiOldIndex++) {
	  plOldCursor = aplOldArray[uiOldIndex];
	  while (plOldCursor != NULL) {
	    plTemp = plOldCursor->plNextLink;
	    vAdd(plOldCursor);
	    plOldCursor = plTemp;
	  }
	}
	delete[] aplOldArray;
      }
    };

  /**\brief A base class for all iMesh tags.
   *
   * The TagBase class allows pointers to all types of tags to be stored
   * in the same set and allows non-data-specific methods to be called
   * without knowing what kind of data the tag is templatized for.
   */
  class TagBase : public iBase_TagHandle_Private {
    /// The iBase tag name
    std::string strName;
    /// The data type of this tag.
    int tType;
  public:
    /// Constructor.
    TagBase(std::string strNameIn, const int tTypeIn) :
	strName(strNameIn), tType(tTypeIn)
    {
    }

    /// Destructor.
    virtual
    ~TagBase()
    {
    }

    /// Comparator, so that tags can be stored in a set.
    inline friend int
    operator<(const TagBase& TA, const TagBase&TB)
    {
      return TA.strName < TB.strName;
    }

    /// Returns the tag name as a C/Fortran string.
    void
    vGetTagName(char* name, int name_len) const
    {
      snprintf(name, name_len, "%s", strName.data());
    }

    /// Returns the tag type.
    int
    iGetTagType() const
    {
      return tType;
    }

    /// Returns the number of objects that have this tag applied to them.
    virtual int
    iGetSize() const = 0;

    /// Returns true iff there are no objects with this tag.
    virtual bool
    qIsEmpty() const = 0;

    /// Returns true iff this object is tagged with this tag.
    virtual bool
    qIsTagged(iBase_Taggable* pHandle) const = 0;

    /// Untag entity and return true if it was tagged; otherwise return false.
    virtual bool
    qUntagEntity(iBase_Taggable* pHandle) = 0;

    /**\brief Untag an array of entities
     *
     * Return value is the last entity that was -not- tagged, or the
     * array size if all given entities were tagged.
     */
    virtual int
    iUntagEntityArr(const iBase_EntityHandle pHandle[], int iArrSize) = 0;
  };

/// Class Tag is a template class implementing all iBase tag types.
  template<class T>
    class Tag : public TagBase {
      /**\brief The number of values for each key.
       *
       * The values are of type T.  The key is of type iBase_Taggable*.
       */
      const int iSize;

      /// A HashMap to store the tag data.
      HashMap<iBase_Taggable*, T> mHandleVal;
    public:
      /// Create a tag of given name, data size, and a data type.
      Tag(std::string strNameIn, int sizeIn, const int tTypeIn) :
	  TagBase(strNameIn, tTypeIn), iSize(sizeIn), mHandleVal(sizeIn)
      {
      }

      /// The destructor does nothing explicitly.
      ~Tag()
      {
      }

      // Inherited from TagBase
      int
      iGetSize() const
      {
	return iSize;
      }
      bool
      qIsEmpty() const
      {
	return mHandleVal.qEmpty();
      }

      /**\brief Retrieve tag values for a given object.
       *
       * Called indirectly by iMesh_getData, iMesh_get[Int|Dbl|EH]Data and
       * the entity set variants.
       */
      bool
      getTagValue(iBase_Taggable* pHandle, T value[]) const
      // PRE:  'pHandle' points to an entity or entity set for which the
      //       tag value is desired
      // POST: If the entity or entity set has already been tagged (is
      //       contained in a tag value set), the value it was tagged
      //       with is stored in 'value', otherwise'value' is set to
      //       the default zero for its type
      {
	return mHandleVal.qGet(pHandle, value);
      }

      /**\brief Retrieve tag values for an array of objects.
       *
       * Called indirectly by iMesh_getArrData, iMesh_get[Int|Dbl|EH]ArrData and
       * the entity set variants.
       */
      int
      getTagValueArr(const iBase_EntityHandle pHandle[], const int iArrSize,
		     T value[]) const
		     // PRE:  'pHandle' points to an entity or entity set for which the
		     //       tag value is desired
		     // POST: If the entity or entity set has already been tagged (is
		     //       contained in a tag value set), the value it was tagged
		     //       with is stored in 'value', otherwise 'value' is set to
		     //       the default zero for its type
      {
	for (int i = 0; i < iArrSize; i++)
	  if (!mHandleVal.qGet(pHandle[i], value + i * iSize))
	    return i;
	return iArrSize;
      }

      /**\brief Set tag values for a given object.
       *
       * Called indirectly by iMesh_setData, iMesh_set[Int|Dbl|EH]Data and
       * the entity set variants.
       */
      void
      setTagValue(iBase_Taggable* pHandle, const T *value)
      // PRE:  'pHandle' points to an entity or entity set to be tagged,
      //       and 'value' contains the desired tag value
      // POST: If the enitity or entity set  has already been tagged (is
      //       present in some value set),  it is untagged (removed from
      //       its current tag value set); either way it is then tagged
      //       with 'value' (added to the desired tag value set)
      {
	mHandleVal.vInsert(pHandle, value);
      }

      /**\brief Set tag values for an array of objects.
       *
       * Called indirectly by iMesh_setArrData, iMesh_set[Int|Dbl|EH]ArrData and
       * the entity set variants.
       */
      void
      setTagValueArr(const iBase_EntityHandle pHandle[], const int iArrSize,
		     const T value[])
      // PRE:  'pHandle' is an array of pointers to  entities or entity sets
      //       to be tagged, 'iArrSize' contains the number of elements in
      //       'pHandle' and 'value', and 'value' is a flattened array of desired tag
      //       values for the corresponding entities/entity sets in 'pHandle'
      // POST: If the entity/entity set pointed to by 'pHandle[i]' has been
      //       tagged (is already contained in a tag value set), it is untagged
      //       (removed from that set); either way it is then tagged with the
      //       value in 'value[i]' (added to the set for'value[i]')
      {
	for (int i = 0; i < iArrSize; i++) {
	  mHandleVal.vInsert(pHandle[i], value + i * iGetSize());
	}
      }

      /**\brief Remove all tags with the given value.
       *
       * In practice, used only to remove entity handle tags when an
       * entity is deleted.
       */
      void
      vRemoveTagValue(const T value[/*iSize*/])
      {
	mHandleVal.vEraseByValue(value);
      }

      /// Return true iff the handle has a value for this tag.
      bool
      qIsTagged(iBase_Taggable* pHandle) const
      // PRE:  'pHandle' points to an entity or entity set
      // POST: If the entity or entity set has been tagged (is contained
      //       in the tag value set), true is returned; otherwise false
      //       is returned
      {
	return mHandleVal.qFind(pHandle);
      }

      /// Remove the tag from this entity.
      bool
      qUntagEntity(iBase_Taggable* pHandle)
      // PRE:  'pHandle' points to a valid entity or entity set
      // POST: If the entity or entity set has been tagged (is in a tag value
      //       set), it is untagged (removed from that set) and 'true' is
      //       returned; if the entity/entity set has not been tagged (is not
      //       in any tag value sets), 'false' is returned
      {
	return mHandleVal.qErase(pHandle);
      }

      /// Remove the tag from these entities.
      int
      iUntagEntityArr(const iBase_EntityHandle pHandle[], int iArrSize)
      // PRE:  'pHandle' is an array of pointers to valid entities or entity
      //       sets, and 'iArrSize' contains the number of elements in 'pHandle'
      //       and the array of bools 'found'
      // POST: For each of the entities/entity sets, if it is already tagged
      //       (present in a tag value set), it is untagged (removed from that
      //       set) and 'found' at that index is set to true; otherwise 'found'
      //       at that index is set to false
      {
	int iErr = iArrSize;
	for (int i = 0; i < iArrSize; i++)
	  if (!mHandleVal.qErase(pHandle[i]))
	    iErr = i;
	return iErr;
      }
    };
}

#endif
