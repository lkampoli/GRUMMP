#ifndef GR_COARSENMANAGER_H
#define GR_COARSENMANAGER_H 1

#include "GR_BaseQueueEntry.h"
#include "GR_BFace.h"
#include "GR_Face.h"
#include "GR_SwapManager.h"
#include "GR_SwapDecider.h"
#include "GR_Subseg.h"
#include "GR_TreeConstructor.h"

#include <deque>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <functional>
#include <tr1/unordered_map>

class Length_Old;

namespace GRUMMP
{

  /*
   * This functionality factored out of mesh member functions.
   */
  class CoarsenManager : public Observer {

    //The mesh to refine.
    Mesh* m_mesh;

    //The length scale calculator.
    std::shared_ptr<Length> m_length;

    // For unidirectional coarsening
    bool m_qFillBdry;
    int m_intSkip;

    //Make these private (for now).
    CoarsenManager();
    CoarsenManager(const CoarsenManager&);
    CoarsenManager&
    operator=(const CoarsenManager&);

    int
    markChain(Vert * const pVSurf, Vert * const pVFirst,
	      const VertConnect aVC[], const double dMaxDist,
	      const int iSkipNum);
    void
    deactivateChain(Vert * const pVSurf, Vert * const pVFirst,
		    const VertConnect aVC[], const double dMaxDist) const;
    void
    addNormalVerts(const VertConnect aVC[],
		   const Vert::VertType VTVertTypeToStartFrom);
    void
    addToCoarseVerts(const VertConnect aVC[]) const;
    void
    selectCoarseVerts(VertConnect aVC[]);
    void
    deactivateAll() const;
    int
    activate(const Vert::VertType VT) const;
    void
    initConflictGraph(VertConnect aVC[]) const;
    void
    updateConflictGraph(VertConnect aVC[]) const;
    void
    removeTaggedVerts(bool smooth = true, std::set<Vert*> *cantRemove =
    NULL);

  public:

    //The constructor...
    CoarsenManager(Mesh* const mesh, std::shared_ptr<Length> length) :
	m_mesh(mesh), m_length(length), m_qFillBdry(true), m_intSkip(3)
    {
      m_mesh->addObserver(
	  this,
	  Observable::cellCreated | Observable::cellDeleted
	      | Observable::bfaceCreated); // might need to change this enum list, not sure what it does
    }

    //... and destructor.
    ~CoarsenManager()
    {
      if (m_mesh)
	m_mesh->removeObserver(this);
    }

    void
    coarsen(const double dLengthRatio = 2);
    void
    coarsenToLengthScale(const bool qHeal = false);

    int
    getInteriorSkip() const
    {
      return m_intSkip;
    }

    void
    setInteriorSkip(int intSkip)
    {
      m_intSkip = intSkip;
    }

    bool
    fillBoundary() const
    {
      return m_qFillBdry;
    }

    void
    setFillBoundary(bool qFillBdry)
    {
      m_qFillBdry = qFillBdry;
    }
  };

}
#endif
