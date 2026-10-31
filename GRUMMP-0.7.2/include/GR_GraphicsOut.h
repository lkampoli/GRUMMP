#ifndef GR_GRAPHICS_OUT_H
#define GR_GRAPHICS_OUT_H

#include "GR_config.h"
#include <vector>
#include <string>

class GraphicsOut {

  //Mostly used for debugging. It outputs
  //CGM facetted points / curves / surfaces to 
  //a medit formatted file.

  std::string output_filename;
  struct triplet {
    int ind1, ind2, ind3;
  };

public:

  GraphicsOut() :
      output_filename()
  {
  }
  GraphicsOut(const std::string& file) :
      output_filename(file)
  {
  }
  GraphicsOut(const char* file) :
      output_filename(file)
  {
  }
  ~GraphicsOut()
  {
  }

  //Outputs all the components
  void
  output_edges() const;
  void
  output_faces() const;

private:

  //Make the copy constructor and operator= private. There should
  //not be any use for them anyway.

  GraphicsOut(const GraphicsOut&); // Never defined
  GraphicsOut&
  operator=(const GraphicsOut&); // Never defined

};

#endif
