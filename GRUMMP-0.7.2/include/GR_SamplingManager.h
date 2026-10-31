#ifndef GR_SAMPLINGMANAGER_H
#define GR_SAMPLINGMANAGER_H 1

#include "GR_config.h"
#include "CubitVector.hpp"
#include <map>

class RefFace;
class RefFaceSampler;
class SubsegManager;
class Vert;

//This class acts as a layer between the RefFaceSamplers
//and the TetMeshBuilder.

class SamplingManager {

  //Maps geometry definition to its sampler.
  typedef std::map<RefFace*, RefFaceSampler*> SamplerMap;

  //The map itself.
  SamplerMap m_sampler_map;

public:

  //Struct defined to return insertion data.
  //Might be extended if need be.
  struct NextInsertData {
    CubitVector insert_coord;
    RefFace* insert_face;
  };

  //Constructor.
  SamplingManager();

  //Destructor.
  ~SamplingManager();

  //Returns the location of the next vertex to insert.
  void
  get_next_insert_data(NextInsertData& next_insert_data) const;

  //Creates a copy of vertex and inserts it into the FaceSampler of 
  //ref_face. Returns the address of the vertex copy.
  Vert*
  insert_vertex(const Vert* const vertex, const RefFace* const ref_face);

  //Returns the number of samplers (which is the size of m_sampler_map)
  int
  num_samplers() const;

  void
  output_stitched_surface(const char* const filename) const;

  //Outputs the restricted Delaunay of all face samplers in m_sampler_map to file.
  void
  output_all_restricted_triangulations() const;

private:

  SamplingManager(const SamplingManager&);
  SamplingManager&
  operator=(const SamplingManager&);

  //Initializes a FaceSampler for every RefFace in the geometric model.
  void
  init_face_samplers();

  //Deletes all the FaceSamplers and clears m_sampler_map contents.
  void
  destroy_face_samplers();

};

#endif
