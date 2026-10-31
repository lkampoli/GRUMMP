#ifndef GR_LENGTH_H
#define GR_LENGTH_H

#include "GR_config.h"

#include <memory>
#include <queue>
#include <set>
#include <utility>

#include "CubitBox.hpp"
#include "CubitVector.hpp"

#include "GR_Classes.h"

namespace GRUMMP
{

  class Length {

  private:

    //Struct that cleans up queries to Mesh::vNeighborhood.
    struct Neighborhood {
    private:
//    Neighborhood(const Neighborhood& neigh);
      Neighborhood&
      operator=(const Neighborhood& neigh);
    public:
      Vert* my_vert;
      std::set<Vert*> verts;
      std::set<Cell*> cells;
      std::set<BFace*> bfaces;
    public:
      Neighborhood() :
	  my_vert(NULL), verts(), cells(), bfaces()
      {
      }
      void
      clear()
      {
	my_vert = NULL, verts.clear();
	cells.clear();
	bfaces.clear();
      }
      void
      getNeighborhood(Vert* const vertex);
    };
  protected:
    static const double defaultMinLength;

    //Refinement and grading constants.
    double m_refine, m_grade, m_minLength;

    // Interpolate length scales in this mesh.  This needs to be an
    // actual copy of the mesh passed in.
    std::unique_ptr<Mesh> m_lengthMesh;

    //Tree of mesh cells for efficient navigation.
    ADT* m_BBoxTree;

    // Should length scale be interpolated?  Very useful for adaptation with
    // large ratio change in length scale.
    bool m_interpolateLS;

  private:
    Length&
    operator=(const Length&); // Never defined

  protected:
    //To queue vertices requiring a grading check.
    std::set<Vert*> m_checkLength;
    void
    setGradedLength();

    //Initializes the length scale of vertices in mesh.
    void
    initLength();

  private:
    double
    evaluateLengthFromNeighbors(Vert* const vertex,
				const std::set<Vert*>& neigh_verts) const;

  protected:
    Length(const double& refine, const double& grading,
	   const double& minLength = defaultMinLength);

    virtual double
    distToNearestBdryEntity(const Vert* const vert,
			    const std::set<Cell*>& cells) const = 0;
  public:
    // Factory for creating descendents of Length.
    static std::shared_ptr<Length>
    Create(const Mesh* const & pM, const double& refine, const double& grading,
	   const double& minLength = defaultMinLength);

    //Destructor
    virtual
    ~Length();

    //Computes the length scale for a given vertex
    //(neighbor vertices' length scale must be properly initialized).
    void
    setLengthScale(Vert* const vertex);

    //Computes the length scale for all verts in a mesh
    void
    setLengthScale(Mesh* const mesh);

    // Computes the length scale for a given location that is -not- a vertex
    virtual double
    queryLengthScale(const double coords[]) const = 0;

    //Multiplies all length scales by dLengthRatio.
    void
    scale(double dLengthRatio);

    //Adds an entry to the length tree.
    virtual void
    addToLengthMesh(double length, double x_loc, double y_loc,
		    double z_loc = 0) = 0;

    // Set and read whether to interpolate length scale within cells.  Should
    // only be true for adaptation with large ratio of length scales before
    // and after.
    bool
    isInterpolatedLS() const;
    void
    setInterpolateLS(bool interpolateLS);

    // A temporary function for testing.  Will become obsolete. TODO: Get rid of it.
    Mesh*
    getMesh() const
    {
      return m_lengthMesh.get();
    }
  };

  class Length2D : public Length {
  protected:
    virtual double
    distToNearestBdryEntity(const Vert* const vert,
			    const std::set<Cell*>& cells) const;
  public:
    Length2D(const Mesh* const & pM, const double& refine,
	     const double& grading, const double& minLength = defaultMinLength);
    ~Length2D()
    {
    }

    double
    queryLengthScale(const double coords[]) const;

    void
    addToLengthMesh(double length, double x_loc, double y_loc,
		    double z_loc = 0);
  };

  class Length3D : public Length {
  protected:
    virtual double
    distToNearestBdryEntity(const Vert* const vert,
			    const std::set<Cell*>& cells) const;
  public:
    Length3D(const Mesh* const & pM, const double& refine,
	     const double& grading, const double& minLength = defaultMinLength);
    ~Length3D()
    {
    }

    // Computes the length scale for a given location that is -not- a vertex
    virtual double
    queryLengthScale(const double coords[]) const;
    void
    addToLengthMesh(double length, double x_loc, double y_loc,
		    double z_loc = 0);
  };
} // namespace GRUMMP

#endif
