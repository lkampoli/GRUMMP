#ifndef GR_LENGTHANISO_H
#define GR_LENGTHANISO_H

#include "GR_config.h"
#include <utility>
#include <queue>
#include <set>
#include "CubitBox.hpp"
#include "CubitVector.hpp"

class Cell;
class Vert;
class Mesh;
class ADT;

class LengthAniso {

public:

  enum eLengthType {
    eFunc, eNeigh, eMesh
  };

protected:

  //Struct that cleans up queries to Mesh::vNeighborhood.
  struct Neighborhood {
    Vert* my_vert;
    std::set<Vert*> verts;
    std::set<Cell*> cells;
    Neighborhood() :
	my_vert(NULL), verts(), cells()
    {
    }
    void
    clear()
    {
      my_vert = NULL, verts.clear();
      cells.clear();
    }
  private:
//    Neighborhood(const Neighborhood& neigh);
    Neighborhood&
    operator=(const Neighborhood& neigh);
  };

  //What type of length scale are we using?
  eLengthType m_lengthType;

  //Refinement and grading constants.
  double m_refine, m_grade, m_minLength;

  void
  (*m_lengthFunc)(const double adLoc[], double adRetVal[]);

  //Interpolate length metrics in this mesh.
  Mesh* m_lengthMesh;

  //Tree of mesh cells for efficient navigation.
  ADT* m_BBoxTree;

  //To queue vertices requiring a grading check.
  std::queue<Vert*> m_lengthCheckQueue;

  //A neighborhood object declared here to save time on initilization.
  Neighborhood m_neigh;

  LengthAniso(); // Never defined
  LengthAniso&
  operator=(const LengthAniso&); // Never defined

  void
  initLengthFunc(Vert* const vertex);

  void
  initLengthNeigh(Vert* const vertex);

  void
  initLengthMesh(Vert* const vertex);

  void
  setLengthScaleFunc(Vert* const vertex);

  void
  setLengthScaleNeigh(Vert* const vertex);

  void
  setLengthScaleMesh(Vert* const vertex);

  void
  setGradedLength();

  void
  getNeighborhood(Vert* const vertex);

  void
  evaluateLengthFromFunc(Vert* const vertex, double* adMetric) const;

  virtual void
  evaluateLengthFromNeighbors(Vert* const vertex, double* adMetric,
			      const std::set<Vert*>& neigh_verts) const = 0;

  virtual void
  evaluateLengthFromMesh(Vert* const vertex, double* adMetric) const = 0;

public:

  //Public constructor.
  //length_type: the type of length scale used:  eAutomatic, eTree.
  //refine: the refinement constant,
  //grade: the grading constant
  LengthAniso(const eLengthType& lengthType, double refine, double grading,
	      double min_length = -LARGE_DBL);

  //Destructor
  virtual
  ~LengthAniso();

  //Initializes the length metric of vertices in mesh.
  //CAUTION: AUTOMATIC type can only have boundary verts in the mesh.
  void
  initLength(const Mesh* const mesh);

  //Set length metric for all vertices in the mesh.
  void
  setLengthScale(const Mesh* const mesh);

  //Computes the length metric for a given vertex
  void
  setLengthScale(Vert* const vertex);

  //Provides a function to calculate metric from vertex coordinates.
  void
  provideFunc(void
  (*func)(const double adLoc[], double adRetVal[]));

  //Provides a background mesh for length scale interpolation.
  virtual void
  provideMesh(Mesh* pM, char* strLenFileName = NULL) = 0;
};

#endif
