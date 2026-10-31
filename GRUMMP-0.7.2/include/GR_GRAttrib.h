#ifndef GR_GRATTRIB_H
#define GR_GRATTRIB_H

#include "GR_config.h"
#include "GR_misc.h"
#include "CubitDefines.h"
#include "CubitString.hpp"

class CubitSimpleAttrib;
class CubitString;

class GRAttrib {
  GRAttrib(GRAttrib&);
  GRAttrib&
  operator=(GRAttrib&);

public:

  //Constructor, Destructor.
  GRAttrib(CubitSimpleAttrib* csa);
  ~GRAttrib();

  //Returns a new CubitSimpleAttrib (look out for leaks)
  CubitSimpleAttrib*
  get_csa() const;

  //Check if two attribs are equal.
  bool
  equals(CubitSimpleAttrib* csa) const;

  //Write to a file at the current file offset
  CubitStatus
  save(FILE* save_file) const;

  //Read from file starting at current file offset
  //Returns a new GRAttrib object. Careful for leaks.
  static GRAttrib*
  restore(FILE* restore_file, unsigned int endian);

  CubitString
  name() const
  {
    return string_array[0];
  }

  CubitString* string_array;
  double* double_array;
  int* int_array;

  int n_strings;
  int n_doubles;
  int n_integers;

  GRAttrib* list_next;

private:

  //Private constructor for use by restore(FILE*)
  GRAttrib(int string_count, CubitString* strings, int double_count,
	   double* doubles, int int_count, int* integers);

};

#endif
