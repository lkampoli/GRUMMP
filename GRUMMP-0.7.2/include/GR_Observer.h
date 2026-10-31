#ifndef GR_Observer_h
#define GR_Observer_h

#include "GR_Observable.h"

namespace GRUMMP
{
  class Observer {
  public:
    Observer()
    {
    }
    virtual
    ~Observer()
    {
    }
    virtual void
    detachFromObservable()
    {
    }

  public:
    // By default, all these receiving functions do nothing.  Individual
    // Observer derived classes will have to override some of these to
    // do useful things.
    virtual void
    receiveDeletedVerts(std::vector<Vert*>& /*deletedVerts*/)
    {
    }
    virtual void
    receiveDeletedFaces(std::vector<Face*>& /*deletedFaces*/)
    {
    }
    virtual void
    receiveDeletedCells(std::vector<Cell*>& /*deletedCells*/)
    {
    }
    virtual void
    receiveDeletedBFaces(std::vector<BFace*>& /*deletedBFaces*/)
    {
    }

    virtual void
    receiveCreatedVerts(std::vector<Vert*>& /*createdVerts*/)
    {
    }
    virtual void
    receiveCreatedFaces(std::vector<Face*>& /*createdFaces*/)
    {
    }
    virtual void
    receiveCreatedCells(std::vector<Cell*>& /*createdCells*/)
    {
    }
    virtual void
    receiveCreatedBFaces(std::vector<BFace*>& /*createdBFaces*/)
    {
    }

    virtual void
    receiveMovedVerts(std::vector<Vert*>& /*movedVerts*/)
    {
    }
  };
}

#endif
