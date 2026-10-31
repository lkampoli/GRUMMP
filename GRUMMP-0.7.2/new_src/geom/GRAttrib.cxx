#include <string.h>
#include "GR_GRAttrib.h"
#include "CubitString.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitFileIOWrapper.hpp"
#include "DLIList.hpp"

GRAttrib::GRAttrib(CubitSimpleAttrib* csa) :
    string_array(NULL), double_array(NULL), int_array(NULL), n_strings(0), n_doubles(
	0), n_integers(0), list_next(0)
{
  //Constructor - copy info from CubitSimpleAttrib

  int i;

  //Reset all lists
  csa->string_data_list()->reset();
  csa->int_data_list()->reset();
  csa->double_data_list()->reset();

  //Save counts
  n_strings = csa->string_data_list()->size();
  n_doubles = csa->double_data_list()->size();
  n_integers = csa->int_data_list()->size();

  //Allocate arrays, but don't try to allocate zero-length arrays
  string_array = n_strings ? new CubitString[n_strings] : NULL;
  double_array = n_doubles ? new double[n_doubles] : NULL;
  int_array = n_integers ? new int[n_integers] : NULL;

  //Copy data into arrays
  for (i = 0; i < n_strings; i++) {
    string_array[i] = *csa->string_data_list()->next(i);
  }
  for (i = 0; i < n_integers; i++) {
    int_array[i] = *csa->int_data_list()->next(i);
  }
  for (i = 0; i < n_doubles; i++) {
    double_array[i] = *csa->double_data_list()->next(i);
  }

}

//Private constructor for use by restore(FILE*)
GRAttrib::GRAttrib(int string_count, CubitString strings[], int double_count,
		   double doubles[], int int_count, int integers[]) :
    string_array(strings), double_array(doubles), int_array(integers), n_strings(
	string_count), n_doubles(double_count), n_integers(int_count), list_next(
	0)
{
}

// Destructor - Free arrays
GRAttrib::~GRAttrib()
{

  //Deleteing NULL pointers is okay.
  delete[] int_array;
  delete[] double_array;
  delete[] string_array;

}

//Copy this into a new CubitSimpleAttrib
CubitSimpleAttrib*
GRAttrib::get_csa() const
{

  // Set initial list size
  DLIList<CubitString*> string_list(n_strings);
  DLIList<int> int_list(n_integers);
  DLIList<double> double_list(n_doubles);

  // Don't need to 'new' objects in DLIList because
  // CSA will make copies.  Just put addresses in list.
  int i;
  for (i = 0; i < n_strings; i++) {
    string_list.append(&string_array[i]);
  }
  for (i = 0; i < n_integers; i++) {
    int_list.append(int_array[i]);
  }
  for (i = 0; i < n_doubles; i++) {
    double_list.append(double_array[i]);
  }

  return new CubitSimpleAttrib(&string_list, &double_list, &int_list);

}

//Compare to a CubitSimpleAttrib
bool
GRAttrib::equals(CubitSimpleAttrib* csa) const
{

  //Compare counts
  if (csa->int_data_list()->size() != n_integers
      || csa->double_data_list()->size() != n_doubles
      || csa->string_data_list()->size() != n_strings) {
    return false;
  }

  //Compare strings first.
  int i;
  csa->string_data_list()->reset();
  for (i = 0; i < n_strings; i++) {
    if (string_array[i] != *csa->string_data_list()->next(i)) {
      return false;
    }
  }

  //Compare integers
  csa->int_data_list()->reset();
  for (i = 0; i < n_integers; i++) {
    if (int_array[i] != *csa->int_data_list()->next(i)) {
      return false;
    }
  }

  //Compare doubles with some fuzz.
  csa->double_data_list()->reset();
  for (i = 0; i < n_doubles; i++) {
    if (fabs(double_array[i] - *csa->double_data_list()->next(i)) > 1.e-200) {
      return false;
    }
  }

  //Everything is equal, return true.
  return true;

}

//Write to a file at the current file offset
CubitStatus
GRAttrib::save(FILE *save_file) const
{

  if (save_file == NULL) {
    PRINT_ERROR("Problem saving MBG attributes: null FILE ptr\n");
    return CUBIT_FAILURE;
  }

  NCubitFile::CIOWrapper wrapper(save_file);

  //Write a version number for the attribute data
  unsigned int attrib_version = 1;
  wrapper.Write(&attrib_version, 1);

  //Write the number of strings, number of doubles, and number of integers
  int counts[3] = {n_strings, n_doubles, n_integers};
  wrapper.Write(reinterpret_cast<unsigned int*>(counts), 3);

  //Write the string data
  int i;
  for (i = 0; i < n_strings; i++) {
    wrapper.Write(string_array[i].c_str());
  }
  //Write the doubles
  wrapper.Write(double_array, n_doubles);
  //Write the integers
  wrapper.Write(reinterpret_cast<unsigned int*>(int_array), n_integers);

  return CUBIT_SUCCESS;

}

//Read from file starting at current file offset
GRAttrib*
GRAttrib::restore(FILE *restore_file, unsigned int endian)
{

  if (restore_file == NULL) {
    return NULL;
  }

  NCubitFile::CIOWrapper wrapper(endian, restore_file);

  //Read a version number for the attribute data
  unsigned int version;
  wrapper.Read(&version, 1);

  //Only version 1 exists.
  if (version != 1) {
    PRINT_ERROR("Wrong GRAttrib version : %u\n", version );
    return NULL;
  }

  //Read the number of strings, number of doubles, and number of integers
  int counts[3];
  wrapper.Read(reinterpret_cast<unsigned int*>(counts), 3);
  int num_strings = counts[0];
  int num_doubles = counts[1];
  int num_ints = counts[2];

  //Allocate arrays, but don't try to allocate zero-length array
  CubitString* strings = num_strings ? new CubitString[num_strings] : NULL;
  double *doubles = num_doubles ? new double[num_doubles] : NULL;
  int *ints = num_ints ? new int[num_ints] : NULL;

  //Read the string data
  int i;
  for (i = 0; i < num_strings; i++) {
    strings[i] = CubitString(wrapper.Read());
  }

  //Write the doubles
  wrapper.Read(doubles, num_doubles);

  //Write the integers
  wrapper.Read(reinterpret_cast<unsigned int*>(ints), num_ints);

  return new GRAttrib(num_strings, strings, num_doubles, doubles, num_ints,
		      ints);

}

