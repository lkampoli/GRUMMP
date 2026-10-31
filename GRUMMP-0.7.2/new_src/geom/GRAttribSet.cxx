#include "GR_GRAttribSet.h"
#include "GR_GRAttrib.h"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"

void
GRAttribSet::append_attribute(CubitSimpleAttrib* csa)
{

  GRAttrib* new_attrib = new GRAttrib(csa);
  new_attrib->list_next = list_head;
  list_head = new_attrib;

}

void
GRAttribSet::remove_attribute(CubitSimpleAttrib* csa)
{

  if (!list_head) {
    return;
  }

  GRAttrib* attrib = 0;
  if (list_head->equals(csa)) {
    attrib = list_head;
    list_head = list_head->list_next;
    delete attrib;
    return;
  }

  for (GRAttrib* prev = list_head; prev->list_next; prev = prev->list_next) {
    if (prev->list_next->equals(csa)) {
      attrib = prev->list_next;
      prev->list_next = attrib->list_next;
      delete attrib;
      return;
    }
  }

}

void
GRAttribSet::remove_all_attributes()
{

  while (list_head) {
    GRAttrib* dead = list_head;
    list_head = dead->list_next;
    delete dead;
  }

}

CubitStatus
GRAttribSet::get_attributes(DLIList<CubitSimpleAttrib*>& attrib_list) const
{

  for (GRAttrib* attrib = list_head; attrib; attrib = attrib->list_next) {
    attrib_list.append(attrib->get_csa());
  }
  return CUBIT_SUCCESS;

}

CubitStatus
GRAttribSet::get_attributes(const CubitString& name,
			    DLIList<CubitSimpleAttrib*>& attrib_list) const
{

  for (GRAttrib* attrib = list_head; attrib; attrib = attrib->list_next) {
    if (attrib->name() == name) {
      attrib_list.append(attrib->get_csa());
    }
  }
  return CUBIT_SUCCESS;

}

CubitStatus
GRAttribSet::save_attributes(FILE* f) const
{

  GRAttrib* curr_attrib;
  CubitStatus status = CUBIT_SUCCESS;

  //Save number of attribs
  unsigned int size = attribute_count();
  NCubitFile::CIOWrapper wrapper(f);
  wrapper.Write(&size, 1);

  //Save each attrib
  for (curr_attrib = list_head; curr_attrib;
      curr_attrib = curr_attrib->list_next) {
    if (!curr_attrib->save(f)) {
      status = CUBIT_FAILURE;
    }
  }

  return status;

}

CubitStatus
GRAttribSet::restore_attributes(FILE* f, unsigned endian)
{

  GRAttrib* curr_attrib;

  //Read number of attribs
  unsigned int size;
  NCubitFile::CIOWrapper wrapper(endian, f);
  wrapper.Read(&size, 1);

  for (unsigned i = 0; i < size; i++) {
    curr_attrib = GRAttrib::restore(f, endian);

    if (!curr_attrib) { // file corrupt?  don't try to read any more
      return CUBIT_FAILURE;
    }

    curr_attrib->list_next = list_head;
    list_head = curr_attrib;
  }

  return CUBIT_SUCCESS;

}

int
GRAttribSet::attribute_count() const
{
  int count = 0;
  for (GRAttrib* attrib = list_head; attrib; attrib = attrib->list_next) {
    count++;
  }
  return count;
}

