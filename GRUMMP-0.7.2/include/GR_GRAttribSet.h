#ifndef GR_GRATTRIBSET_H
#define GR_GRATTRIBSET_H

#include "GR_config.h"
#include "GR_misc.h"
#include "DLIList.hpp"

//Forward declarations
class CubitSimpleAttrib;
class CubitString;
class GRAttrib;

class GRAttribSet {
  GRAttribSet&
  operator=(const GRAttribSet&);
  GRAttribSet(const GRAttribSet&);
public:

  GRAttribSet() :
      list_head(NULL)
  {
  }

  ~GRAttribSet()
  {
    remove_all_attributes();
  }

  void
  append_attribute(CubitSimpleAttrib* csa);

  void
  remove_attribute(CubitSimpleAttrib* csa);

  void
  remove_all_attributes();

  CubitStatus
  get_attributes(DLIList<CubitSimpleAttrib*>& attrib_list) const;

  CubitStatus
  get_attributes(const CubitString& name,
		 DLIList<CubitSimpleAttrib*>& attrib_list) const;

  CubitStatus
  save_attributes(FILE* f) const;

  CubitStatus
  restore_attributes(FILE* f, unsigned int endian);

  int
  attribute_count() const;

private:

  GRAttrib* list_head;

};

#endif
