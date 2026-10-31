/*
 * GR_InterfaceTracker.h
 *
 *  Created on: Apr 12, 2017
 *      Author: cfog
 */

#ifndef APPS_SHOCK_TRACK_GR_INTERFACETRACKER_H_
#define APPS_SHOCK_TRACK_GR_INTERFACETRACKER_H_

#include "GR_Observer.h"
#include "GR_InsertionManager.h"

namespace GRUMMP_Apps
{

  class InterfaceTracker : public GRUMMP::Observer {
    InterfaceTracker(const InterfaceTracker&); // Null
    InterfaceTracker&
    operator=(const InterfaceTracker&); // Null
  public:
    InterfaceTracker();
    ~InterfaceTracker();
    /// Read a volume mesh, in any format GRUMMP reads.
    void
    readVolumeMesh(const char volumeFilename[]);
    /// Read an STL file and add it to the mesh as an internal surface.
    void
    addInterfaceToMesh(const char surfaceFilename[]);
    /// Split surface open (duplicating the two sides)
    void
    splitInterface();
    /// Output the volume mesh to VTK and Tetgen formats
    void
    writeVolumeMesh();
    /// Output the number of faces on a shock surface
    GR_index_t
    outputNumberOfShockFaces();
    /// Update the mesh connectivity when the shock has moved.
    void
    updateMesh();

  private:
    VolMesh* m_compMesh;
    SurfMesh* m_surface;
    std::shared_ptr<GRUMMP::Length> m_length;
    GRUMMP::SwappingInserter3D* m_pSI3D;
    std::map<Vert*, std::pair<Vert*, Vert*> > m_surfToVolMap;
    int m_numVars;
    double *m_solnData;
    static constexpr double m_minDist = 0.4, m_maxDist = 0.6;

    bool
    isCloseToShock(Vert* const pV, const ADT& BBoxTree, BFace* bdryFaces[],
		   const double ratio);
    bool
    makeVertexPhantom(Vert* const pV, int& newSwaps);
    bool
    makePhantomReal(Vert* const pV);
    void
    improveMesh();
    void
    swapNearShocks(GRUMMP::SwapManager3D& SM3D);
  };

  void
  readAndManipulateSurface(const char fileName[]);
} // end namespace

#endif /* APPS_SHOCK_TRACK_GR_INTERFACETRACKER_H_ */
