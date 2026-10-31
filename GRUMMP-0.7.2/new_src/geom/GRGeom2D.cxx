#include "CubitVector.hpp"  //CGM vector class.
#include "GR_GRGeom2D.h"    //This class header file
#include "GR_GRPoint.h"     //GRUMMP point class.
#include "GR_GRCurve.h"     //GRUMMP curve class.
#include "GR_GRCoEdge.h"    //GRUMMP co-edge class.
#include "GR_GRLoop.h"      //GRUMMP loop class.
#include "DLIList.hpp"
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>

GRGeom2D::GRGeom2D() :
    geom_initialized(false), n_points(0), n_curves(0), n_constructed_curves(0), point_list(), curve_list(), coedge_list(), loop_list(), v_points(), v_curves()
{
  assert(iInvalidRegion == 0);
}

GRGeom2D::GRGeom2D(const char* const bdry_file) :
    geom_initialized(false), n_points(0), n_curves(0), n_constructed_curves(0), point_list(), curve_list(), coedge_list(), loop_list(), v_points(), v_curves()
{
  assert(iInvalidRegion == 0);

  char bdry_file_name[FILE_NAME_LEN];
  makeFileName(bdry_file_name, "%s.bdry", bdry_file,
		  "GRGeom2D::GRGeom2D(const string)");

  FILE* f = open_bdry_file(bdry_file_name);
  read_num_entities(f);
  read_points(f);
  read_curves(f);
  fclose(f);

  point_list.clear();
  curve_list.clear();
  coedge_list.clear();
  loop_list.clear();

  build_geometry();

}

GRGeom2D::~GRGeom2D()
{

  std::list<GRPoint*>::iterator itp;
  for (itp = point_list.begin(); itp != point_list.end(); itp++) {
    if (*itp)
      delete *itp;
    *itp = NULL;
  }

  std::list<GRCurve*>::iterator itc;
  for (itc = curve_list.begin(); itc != curve_list.end(); itc++) {
    //if((*itc)->get_geometry()) (*itc)->delete_geometry();
    if (*itc)
      delete *itc;
    *itc = NULL;
  }

  std::list<GRCoEdge*>::iterator itce;
  for (itce = coedge_list.begin(); itce != coedge_list.end(); itce++) {
    if (*itce)
      delete *itce;
    *itce = NULL;
  }

  std::list<GRLoop*>::iterator itl;
  for (itl = loop_list.begin(); itl != loop_list.end(); itl++) {
    if (*itl)
      delete *itl;
    *itl = NULL;
  }

  point_list.clear();
  curve_list.clear();
  coedge_list.clear();
  loop_list.clear();
  v_points.clear();
  v_curves.clear();
}

GRPoint*
GRGeom2D::build_point(double x_coord, double y_coord)
{

  GRPoint* point = new GRPoint(x_coord, y_coord, 0.);
  point_list.push_back(point);
  return point;

}

GRCurve*
GRGeom2D::build_line(GRPoint* const point1, GRPoint* const point2,
		     int region_left, int region_right, int bound_left,
		     int bound_right)
{

  GRLine* line = new GRLine(point1->coordinates(), point2->coordinates());
  GRCurve* curve = new GRCurve(line, point1, point2);
  point1->add_curve(curve);
  point2->add_curve(curve);
  curve_list.push_back(curve);

  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);
  string_list.append(&reg_left);
  double_list.append(0.);
  int_list.append(region_left);
  string_list.append(&bdry_left);
  double_list.append(0.);
  int_list.append(bound_left);
  string_list.append(&reg_right);
  double_list.append(0.);
  int_list.append(region_right);
  string_list.append(&bdry_right);
  double_list.append(0.);
  int_list.append(bound_right);

  CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
  curve->append_simple_attribute_virt(&csa);

  return curve;

}

void
GRGeom2D::build_geometry()
{

  GR_index_t i;
  std::map<GR_index_t, GRPoint*> created_points;
  std::map<GR_index_t, std::vector<GRCoEdge*> > region_curves;

  for (i = 0; i < v_curves.size(); i++)
    build_curve(i, created_points, region_curves);

  if (n_constructed_curves != n_curves) {
    logMessage(
	0,
	"ERROR: The number of constructed curves is not the expected number.\n");
    logMessage(0, "Constructed: %u, expected: %u\n", n_constructed_curves,
	       n_curves);
    abort();
  }

  std::map<GR_index_t, std::vector<GRCoEdge*> >::iterator it;
  for (it = region_curves.begin(); it != region_curves.end(); it++) {
    std::vector<GRLoop*> created_loops;
    build_loop(it->second, created_loops);
    verify_region_conditions(it->first, created_loops);
  }

  set_initialized();
  //verification();

}

FILE*
GRGeom2D::open_bdry_file(const char* const bdry_file) const
{

  FILE* file = fopen(bdry_file, "r");
  if (file == NULL) {
    logMessage(0, "ERROR: Unable to open 2D boundary file: %s\n", bdry_file);
    abort();
  }
  return file;
}

void
GRGeom2D::read_num_entities(FILE* const f)
{

  char buffer[256];
  vGetLineOrAbort(buffer, 256, f);
  if (2 != sscanf(buffer, "%u %u", &n_points, &n_curves)) {
    logMessage(0, "ERROR: Unable to read the number of points or curves.\n");
  }

  if (n_points <= 0) {
    vFatalError("Invalid number of boundary points",
		"Geom2DCreator::read_num_entities(FILE*)");
  }
  if (n_curves <= 0) {
    vFatalError("Invalid number of boundary curves",
		"Geom2DCreator::read_num_entities(FILE*)");
  }

  logMessage(0, "Number of points: %u \n", n_points);
  logMessage(0, "Number of edges: %u \n", n_curves);

}

void
GRGeom2D::read_points(FILE* const f)
{
  GR_index_t i;
  char buffer[256];
  v_points.clear();

  for (i = 0; i < n_points; i++) {
    vGetLineOrAbort(buffer, 256, f);

    double x = -LARGE_DBL, y = -LARGE_DBL, z = 0.;
    if (2 != sscanf(buffer, "%lf%lf", &x, &y)) {
      logMessage(0, "Some coordinates in boundary file are invalid\n");
      logMessage(0, "Invalid line reads as: %s\n", buffer);
      abort();
    }
    //		x = floor(1e8*x)/1e8;
    //		y = floor(1e8*y)/1e8;
    //		printf("x %.16f y %.16f\n",x,y);
    CubitVector cv(x, y, z);

    v_points.push_back(cv);

  }

  if (v_points.size() != n_points) {
    logMessage(
	0,
	"The number of specified points is different from the expected number.\n");
    logMessage(0, "Expected number: %u, Specified number: %zu\n", n_points,
	       v_points.size());
    abort();
  }

}

void
GRGeom2D::read_curves(FILE* const f)
{

  int max_size = 262144;

  char curve_input[max_size], dum[max_size];
  v_curves.clear();

  //Reading line by line allowing line wrap
  while (!feof(f)) {

    get_line(curve_input, max_size, f);

    if (strlen(curve_input) == 0) {
      assert(feof(f));
      continue;
    }

    char curve_type[32];
    sscanf(curve_input, "%32s", curve_type);
    if (!valid_patch_type(curve_type)) {
      logMessage(0, "ERROR: %s is not a valid curve type.\n", curve_type);
      abort();
    }

    fpos_t pos;
    fgetpos(f, &pos);

    while (1) {

      //Need to take a peek at the next line.
      if (!feof(f))
	get_line(dum, max_size, f);
      else
	break;

      sscanf(dum, "%32s", curve_type);

      if (!valid_patch_type(curve_type)) {
	strcat(curve_input, " ");
	strcat(curve_input, dum);
	fgetpos(f, &pos);
      }
      else {
	fsetpos(f, &pos);
	break;
      }

    }

    //Curve data all in one string, now tokenizing it.
    curve_info ci;
    ci.curve_type = CubitString(strtok(curve_input, " "));
    ci.l_cond = CubitString(strtok(NULL, " "));
    ci.left = atoi(strtok(NULL, " "));
    ci.r_cond = CubitString(strtok(NULL, " "));
    ci.right = atoi(strtok(NULL, " "));

    char* tok = strtok(NULL, " ");
    while (tok != NULL) {
      ci.curve_data.push_back(atof(tok));
      tok = strtok(NULL, " ");
    }
    v_curves.push_back(ci);
  }
}

void
GRGeom2D::get_line(char buffer[], const int size, FILE *f) const
{

  int i;
  do {
    do {
      if (NULL == fgets(buffer, size, f)) {
	strcpy(buffer, "");
	return;
      }
    }
    while (buffer[0] == '!' || buffer[0] == '#' || buffer[0] == '\n'
	|| buffer[0] == '\000'); //Skip comment lines
    for (i = (strlen(buffer)) - 1; i >= 0 && isspace(buffer[i]); i--)
      buffer[i] = '\000'; //Replace whitespaces by NULL characters
  }
  while (strlen(buffer) == 0); //Empty line, get another one

}

bool
GRGeom2D::valid_patch_type(const char* const str) const
{

  if (!strcmp("vertex", str))
    return true;
  if (!strcmp("polyline", str))
    return true;
  if (!strcmp("arc", str))
    return true;
  if (!strcmp("longarc", str))
    return true;
  if (!strcmp("circle", str))
    return true;
  if (!strcmp("bezier", str))
    return true;
  if (!strcmp("spline", str))
    return true;
  return false;

}

void
GRGeom2D::build_curve(
    const GR_index_t curve_id, std::map<GR_index_t, GRPoint*>& created_points,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves)
{

  curve_info* ci = &(v_curves[curve_id]);

  //Verifying validity of the input: Checking if region and boundary
  //specifications are correct.
  CubitString r("r"), b("b");

  if (ci->l_cond == r) {
    if (ci->left <= 0 || ci->left >= 127) {
      logMessage(
	  0,
	  "ERROR: Invalid region index on left side; must be between 1 and 126\n");
      logMessage(0, "Region index specified on left side: %d\n", ci->left);
      abort();
    }
  }
  else if (ci->l_cond == b) {
    if (ci->left <= 0 || ci->left >= 127) {
      logMessage(
	  0,
	  "ERROR: Invalid boundary index on left side; must be between 1 and 126\n");
      logMessage(0, "Boundary index specified on left side: %d\n", ci->left);
      abort();
    }
  }
  else {
    logMessage(
	0, "ERROR: Invalid condition specified on left side; must be r or b\n");
    logMessage(0, "Condition specified on left side: %s\n", ci->l_cond.c_str());
    abort();
  }

  if (ci->r_cond == r) {
    if (ci->right <= 0 || ci->right >= 127) {
      logMessage(
	  0,
	  "ERROR: Invalid region index on right side; must be between 1 and 126\n");
      logMessage(0, "Region index specified on right side: %d\n", ci->right);
      abort();
    }
  }
  else if (ci->r_cond == b) {
    if (ci->right <= 0 || ci->right >= 127) {
      logMessage(
	  0,
	  "ERROR: Invalid boundary index on right side; must be between 1 and 126\n");
      logMessage(0, "Boundary index specified on right side: %d\n", ci->right);
      abort();
    }
  }
  else {
    logMessage(
	0,
	"ERROR: Invalid condition specified on right side; must be r or b\n");
    logMessage(0, "Condition specified on right side: %s\n",
	       ci->r_cond.c_str());
    abort();
  }

  //   if(ci->l_cond == r && ci->r_cond == r) {
  //     if(ci->left == ci->right) {
  //       vMessage(0, "ERROR: Same region specified on both sides of the curve.\n");
  //       abort();
  //     }
  //   }

  if (ci->l_cond == b && ci->r_cond == b) {
    logMessage(
	0, "ERROR: Boundary condition specified on both sides of the curve.\n");
    logMessage(0, "This feature is not yet supported by GRUMMP.\n");
    abort();
  }

  if (!strcmp("vertex", ci->curve_type.c_str()))
    build_vertex(ci, created_points, region_curves);
  else if (!strcmp("polyline", ci->curve_type.c_str()))
    build_polyline(ci, created_points, region_curves);
  else if (!strcmp("arc", ci->curve_type.c_str()))
    build_arc(ci, created_points, region_curves);
  else if (!strcmp("longarc", ci->curve_type.c_str()))
    build_arc(ci, created_points, region_curves, false);
  else if (!strcmp("circle", ci->curve_type.c_str()))
    build_circle(ci, created_points, region_curves);
  else if (!strcmp("bezier", ci->curve_type.c_str()))
    build_bezier(ci, created_points, region_curves);
  else if (!strcmp("spline", ci->curve_type.c_str()))
    build_interp(ci, created_points, region_curves);
  else if (!strcmp("interp", ci->curve_type.c_str()))
    build_interp(ci, created_points, region_curves);
  else {
    logMessage(0, "ERROR: %s is an unknown curve type\n",
	       ci->curve_type.c_str());
    logMessage(
	0,
	"Allowed curve types are: vertex, polyline, arc, longarc, circle, bezier and spline");
    abort();
  }

}

void
GRGeom2D::build_vertex(const curve_info* const ci,
		       std::map<GR_index_t, GRPoint*>& created_points,
		       std::map<GR_index_t, std::vector<GRCoEdge*> >&)
{

  if (strcmp(ci->l_cond.c_str(), "r") != 0
      || strcmp(ci->r_cond.c_str(), "r") != 0) {
    logMessage(
	0,
	"ERROR: A vertex curve must have a region specifier on both sides.\n");
    logMessage(0, "Vertex creation failed.\n");
    abort();
  }

  if (ci->left != ci->right) {
    logMessage(
	0,
	"ERROR: A vertex curve must have the same region specifier on both sides.\n");
    logMessage(0, "Vertex creation failed.\n");
    abort();
  }

  if (ci->curve_data.size() != 1) {
    logMessage(0, "ERROR: More than one point given to define the vertex.\n");
    logMessage(0, "Vertex creation failed.\n");
    abort();
  }

  GR_index_t pt = static_cast<GR_index_t>(round(ci->curve_data.front()));

  if (pt > n_points - 1) {
    logMessage(0, "ERROR: Point specified to define the vertex is invalid.\n");
    logMessage(0, "Vertex creation failed.\n");
    abort();
  }

  std::map<GR_index_t, GRPoint*>::iterator it = created_points.find(pt);
  if (it == created_points.end()) {
    GRPoint* point = new GRPoint(v_points[pt]);
    point_list.push_back(point);
    n_constructed_curves++;
    created_points.insert(std::make_pair(pt, point));
  }
  else {
    logMessage(
	0,
	"ERROR: Trying to create a vertex with an already existing point.\n");
    logMessage(0, "Vertex creation failed.\n");
    abort();
  }

}

void
GRGeom2D::build_polyline(
    const curve_info* const ci, std::map<GR_index_t, GRPoint*>& created_points,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves)
{

  if (ci->curve_data.size() < 3) {
    logMessage(0, "ERROR: Insufficient information to create a polyline.\n");
    logMessage(0, "Polyline creation failed.\n");
    abort();
  }

  GR_index_t n_pts = static_cast<GR_index_t>(ci->curve_data.front());

  if (n_pts < 2) {
    logMessage(
	0, "ERROR: At least two points are required to create a polyline.\n");
    logMessage(0, "Polyline creation failed.\n");
    abort();
  }

  if (ci->curve_data.size() != n_pts + 1) {
    logMessage(
	0,
	"ERROR: The number of points provided does not correspond to the expected number.\n");
    logMessage(0, "Excepted number: %u. Provided number: %zu.\n", n_pts,
	       ci->curve_data.size() - 1);
    logMessage(0, "Polyline creation failed.\n");
    abort();
  }

  for (GR_index_t i = 1; i < ci->curve_data.size() - 1; i++) {

    GR_index_t pt[2];
    GRPoint* p[2];

    pt[0] = static_cast<GR_index_t>(ci->curve_data[i]);
    pt[1] = static_cast<GR_index_t>(ci->curve_data[i + 1]);

    if (pt[0] == pt[1]) {
      logMessage(0, "ERROR: Polyline start and end point are the same.\n");
      logMessage(0, "Polyline creation failed.\n");
      abort();
    }

    if (pt[0] > n_points - 1 || pt[1] > n_points - 1) {
      logMessage(0, "ERROR: One or both of the polyline points are invalid.\n");
      logMessage(0, "Polyline creation failed.\n");
      abort();
    }

    for (int j = 0; j < 2; j++) {
      std::map<GR_index_t, GRPoint*>::iterator it = created_points.find(pt[j]);
      if (it == created_points.end()) {
	GRPoint* point = new GRPoint(v_points[pt[j]]);
	point_list.push_back(point);
	created_points.insert(std::make_pair(pt[j], point));
	p[j] = point;
      }
      else
	p[j] = it->second;
    }

    //Building the line.
    GRLine* line = new GRLine(v_points[pt[0]], v_points[pt[1]]);
    GRCurve* curve = new GRCurve(line, p[0], p[1]);
    p[0]->add_curve(curve);
    p[1]->add_curve(curve);
    curve_list.push_back(curve);

    //Building co-edge(s) and add attributes for left and right side region/bdry.

    DLIList<CubitString*> string_list;
    DLIList<double> double_list;
    DLIList<int> int_list;

    CubitString cond_string("CURVE_CONDITION");
    CubitString reg_left("REG_LEFT");
    CubitString bdry_left("BDRY_LEFT");
    CubitString reg_right("REG_RIGHT");
    CubitString bdry_right("BDRY_RIGHT");

    string_list.append(&cond_string);
    double_list.append(0.);
    int_list.append(0);

    if (strcmp(ci->l_cond.c_str(), "r") == 0) {
      build_coedge(ci->left, curve, region_curves);
      string_list.append(&reg_left);
      double_list.append(0.);
      int_list.append(ci->left);
      string_list.append(&bdry_left);
      double_list.append(0.);
      int_list.append(0);
    }
    else if (strcmp(ci->l_cond.c_str(), "b") == 0) {
      string_list.append(&reg_left);
      double_list.append(0.);
      int_list.append(0);
      string_list.append(&bdry_left);
      double_list.append(0.);
      int_list.append(ci->left);
    }
    else
      assert(0);

    if (strcmp(ci->r_cond.c_str(), "r") == 0) {
      build_coedge(ci->right, curve, region_curves);
      string_list.append(&reg_right);
      double_list.append(0.);
      int_list.append(ci->right);
      string_list.append(&bdry_right);
      double_list.append(0.);
      int_list.append(0);
    }
    else if (strcmp(ci->r_cond.c_str(), "b") == 0) {
      string_list.append(&reg_right);
      double_list.append(0.);
      int_list.append(0);
      string_list.append(&bdry_right);
      double_list.append(0.);
      int_list.append(ci->right);
    }
    else
      assert(0);

    n_constructed_curves++;

    CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
    curve->append_simple_attribute_virt(&csa);

  }

}

void
GRGeom2D::build_arc(
    const curve_info* const ci, std::map<GR_index_t, GRPoint*>& created_points,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves,
    const bool short_arc)
{

  int i;

  if (ci->curve_data.size() != 3) {
    logMessage(0,
	       "ERROR: Information provided to create an arc is incorrect.\n");
    logMessage(0, "Arc creation failed.\n");
    abort();
  }

  double radius = ci->curve_data[0];

  if (radius <= 0.) {
    logMessage(0, "ERROR: Radius provided is %lf. It must be greater than 0.\n",
	       radius);
    logMessage(0, "Arc creation failed.\n");
    abort();
  }

  GR_index_t pt[2];
  GRPoint* p[2];

  for (i = 0; i < 2; i++) {

    pt[i] = static_cast<GR_index_t>(ci->curve_data[i + 1]);

    if (pt[i] > n_points - 1) {
      logMessage(0, "ERROR: One of the specified point index is invalid.\n");
      logMessage(0, "Arc creation failed.\n");
      abort();
    }

    std::map<GR_index_t, GRPoint*>::iterator it = created_points.find(pt[i]);
    if (it == created_points.end()) {
      GRPoint* point = new GRPoint(v_points[pt[i]]);
      point_list.push_back(point);
      created_points.insert(std::make_pair(pt[i], point));
      p[i] = point;
    }
    else
      p[i] = it->second;

  }

  if (pt[0] == pt[1]) {
    logMessage(
	0,
	"ERROR: Arc start and end point are the same. Use circle instead.\n");
    logMessage(0, "Arc creation failed.\n");
    abort();
  }

  //Building the arc
  GRArc* arc = new GRArc(v_points[pt[0]], v_points[pt[1]], radius, short_arc);
  GRCurve* curve = new GRCurve(arc, p[0], p[1]);
  p[0]->add_curve(curve);
  p[1]->add_curve(curve);
  curve_list.push_back(curve);

  //Building co-edge(s) and add attributes for left and right side region/bdry.
  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);

  if (strcmp(ci->l_cond.c_str(), "r") == 0) {
    build_coedge(ci->left, curve, region_curves);
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(ci->left);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->l_cond.c_str(), "b") == 0) {
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(ci->left);
  }
  else
    assert(0);

  if (strcmp(ci->r_cond.c_str(), "r") == 0) {
    build_coedge(ci->right, curve, region_curves);
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(ci->right);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->r_cond.c_str(), "b") == 0) {
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(ci->right);
  }
  else
    assert(0);

  n_constructed_curves++;

  CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
  curve->append_simple_attribute_virt(&csa);

}

void
GRGeom2D::build_circle(
    const curve_info* const ci, std::map<GR_index_t, GRPoint*>&,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves)
{

  if (ci->curve_data.size() != 2) {
    logMessage(
	0, "ERROR: Information provided to create a circle is incorrect.\n");
    logMessage(0, "Circle creation failed.\n");
    abort();
  }

  double radius = ci->curve_data[0];

  if (radius <= 0.) {
    logMessage(0, "ERROR: Radius provided is %lf. It must be greater than 0.\n",
	       radius);
    logMessage(0, "Circle creation failed.\n");
    abort();
  }

  GR_index_t center_point = static_cast<GR_index_t>(ci->curve_data[1]);
  if (center_point > n_points - 1) {
    logMessage(0, "ERROR: The center point index is invalid.\n");
    logMessage(0, "Circle creation failed.\n");
    abort();
  }

  //The circle is actually made of two 180 degrees arcs

  CubitVector center_coord = v_points[center_point];
  CubitVector coord1 = center_coord;
  CubitVector coord2 = center_coord;
  coord1.x(coord1.x() + radius);
  coord2.x(coord2.x() - radius);

  GRPoint* point1 = new GRPoint(coord1);
  GRPoint* point2 = new GRPoint(coord2);
  point_list.push_back(point1);
  point_list.push_back(point2);

  GRArc* arc1 = new GRArc(coord1, coord2, radius, true);
  GRArc* arc2 = new GRArc(coord2, coord1, radius, true);
  GRCurve* curve1 = new GRCurve(arc1, point1, point2);
  GRCurve* curve2 = new GRCurve(arc2, point2, point1);

  point1->add_curve(curve1);
  point1->add_curve(curve2);
  point2->add_curve(curve1);
  point2->add_curve(curve2);
  curve_list.push_back(curve1);
  curve_list.push_back(curve2);

  //Building co-edge(s) and add attributes for left and right side region/bdry.
  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);

  if (strcmp(ci->l_cond.c_str(), "r") == 0) {
    build_coedge(ci->left, curve1, region_curves);
    build_coedge(ci->left, curve2, region_curves);
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(ci->left);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->l_cond.c_str(), "b") == 0) {
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(ci->left);
  }
  else
    assert(0);

  if (strcmp(ci->r_cond.c_str(), "r") == 0) {
    build_coedge(ci->right, curve1, region_curves);
    build_coedge(ci->right, curve2, region_curves);
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(ci->right);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->r_cond.c_str(), "b") == 0) {
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(ci->right);
  }
  else
    assert(0);

  n_constructed_curves++;

  CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
  curve1->append_simple_attribute_virt(&csa);
  curve2->append_simple_attribute_virt(&csa);

}

void
GRGeom2D::build_bezier(
    const curve_info* const ci, std::map<GR_index_t, GRPoint*>& created_points,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves)
{

  int i;

  if (ci->curve_data.size() != 4) {
    logMessage(
	0,
	"ERROR: Information provided to create a bezier curve is incorrect.\n");
    logMessage(0, "Bezier creation failed.\n");
    abort();
  }

  GR_index_t pt[4];
  GRPoint* p[2];

  for (i = 0; i < 4; i++) {
    pt[i] = static_cast<int>(ci->curve_data[i]);

    if (pt[i] > n_points - 1) {
      logMessage(0, "ERROR: One of the specified point index is invalid.\n");
      logMessage(0, "Bezier creation failed.\n");
      abort();
    }

    if (i == 0 || i == 1) {
      std::map<GR_index_t, GRPoint*>::iterator it = created_points.find(pt[i]);
      if (it == created_points.end()) {
	GRPoint* point = new GRPoint(v_points[pt[i]]);
	point_list.push_back(point);
	created_points.insert(std::make_pair(pt[i], point));
	p[i] = point;
      }
      else
	p[i] = it->second;
    }
  }

  CubitVector input_points[4];

  //   printf("point 1: %lf %lf\n", v_points[pt[0]].x(), v_points[pt[0]].y());
  //   printf("point 2: %lf %lf\n", v_points[pt[1]].x(), v_points[pt[1]].y());
  //   printf("point 3: %lf %lf\n", v_points[pt[2]].x(), v_points[pt[2]].y());
  //   printf("point 4: %lf %lf\n", v_points[pt[3]].x(), v_points[pt[3]].y());

  input_points[0] = v_points[pt[0]];
  input_points[1] = v_points[pt[2]];
  input_points[2] = v_points[pt[3]];
  input_points[3] = v_points[pt[1]];

  //Building the bezier curve
  GRBezier* bez = new GRBezier(input_points);
  GRCurve* curve = new GRCurve(bez, p[0], p[1]);
  p[0]->add_curve(curve);
  p[1]->add_curve(curve);
  curve_list.push_back(curve);

  //Building co-edge(s) and add attributes for left and right side region/bdry.
  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);

  if (strcmp(ci->l_cond.c_str(), "r") == 0) {
    build_coedge(ci->left, curve, region_curves);
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(ci->left);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->l_cond.c_str(), "b") == 0) {
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(ci->left);
  }
  else
    assert(0);

  if (strcmp(ci->r_cond.c_str(), "r") == 0) {
    build_coedge(ci->right, curve, region_curves);
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(ci->right);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->r_cond.c_str(), "b") == 0) {
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(ci->right);
  }
  else
    assert(0);

  n_constructed_curves++;

  CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
  curve->append_simple_attribute_virt(&csa);

}

void
GRGeom2D::build_interp(
    const curve_info* const ci, std::map<GR_index_t, GRPoint*>& created_points,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves)
{
  GR_index_t num_pts = static_cast<GR_index_t>(ci->curve_data[0]);

  if (num_pts < 3) {
    logMessage(
	0, "ERROR: An interpolated spline requires at least three points.\n");
    logMessage(0, "Interpolated spline creation failed\n");
  }

  if (ci->curve_data.size() != num_pts + 1) {
    logMessage(
	0,
	"ERROR: The number of points provided does not correspond to the expected number.\n");
    logMessage(0, "Excepted number: %u. Provided number: %zu.\n", num_pts,
	       (ci->curve_data.size()) - 1);
    logMessage(0, "Interpolaed spline creation failed.\n");
    abort();
  }

  std::vector<int> control_points;
  GR_index_t i;
  for (i = 1; i < num_pts + 1; i++) {
    GR_index_t this_point = static_cast<GR_index_t>(ci->curve_data[i]);
    if (this_point > n_points - 1) {
      logMessage(0, "ERROR: One of the specified point index is invalid.\n");
      logMessage(0, "Interpolated spline creation failed.\n");
      abort();
    }

    control_points.push_back(this_point);

  }

  int beg = control_points.front();
  int end = control_points.back();

  GRPoint* beg_pt;
  GRPoint* end_pt;

  std::map<GR_index_t, GRPoint*>::iterator it;
  it = created_points.find(beg);
  if (it == created_points.end()) {
    beg_pt = new GRPoint(v_points[beg]);
    point_list.push_back(beg_pt);
    created_points.insert(std::make_pair(beg, beg_pt));
  }
  else
    beg_pt = it->second;

  it = created_points.find(end);
  if (it == created_points.end()) {
    end_pt = new GRPoint(v_points[end]);
    point_list.push_back(end_pt);
    created_points.insert(std::make_pair(end, end_pt));
  }
  else
    end_pt = it->second;

  CubitVector* coords = new CubitVector[control_points.size()];
  unsigned int j;
  for (j = 0; j < control_points.size(); j++) {
    int index = control_points[j];
    coords[j] = v_points[index];
  }

  GRInterp* interp = new GRInterp(coords, control_points.size());
  GRCurve* curve = new GRCurve(interp, beg_pt, end_pt);
  beg_pt->add_curve(curve);
  end_pt->add_curve(curve);
  curve_list.push_back(curve);

  delete[] coords;

  //Building co-edge(s) and add attributes for left and right side region/bdry.
  DLIList<CubitString*> string_list;
  DLIList<double> double_list;
  DLIList<int> int_list;

  CubitString cond_string("CURVE_CONDITION");
  CubitString reg_left("REG_LEFT");
  CubitString bdry_left("BDRY_LEFT");
  CubitString reg_right("REG_RIGHT");
  CubitString bdry_right("BDRY_RIGHT");

  string_list.append(&cond_string);
  double_list.append(0.);
  int_list.append(0);

  if (strcmp(ci->l_cond.c_str(), "r") == 0) {
    build_coedge(ci->left, curve, region_curves);
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(ci->left);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->l_cond.c_str(), "b") == 0) {
    string_list.append(&reg_left);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_left);
    double_list.append(0.);
    int_list.append(ci->left);
  }
  else
    assert(0);

  if (strcmp(ci->r_cond.c_str(), "r") == 0) {
    build_coedge(ci->right, curve, region_curves);
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(ci->right);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(0);
  }
  else if (strcmp(ci->r_cond.c_str(), "b") == 0) {
    string_list.append(&reg_right);
    double_list.append(0.);
    int_list.append(0);
    string_list.append(&bdry_right);
    double_list.append(0.);
    int_list.append(ci->right);
  }
  else
    assert(0);

  n_constructed_curves++;

  CubitSimpleAttrib csa(&string_list, &double_list, &int_list);
  curve->append_simple_attribute_virt(&csa);

}

void
GRGeom2D::build_coedge(
    int region, GRCurve* const curve,
    std::map<GR_index_t, std::vector<GRCoEdge*> >& region_curves)
{

  GRCoEdge* ce = new GRCoEdge(curve, static_cast<LoopSM*>(NULL), CUBIT_UNKNOWN);
  coedge_list.push_back(ce);
  curve->add_coedge(ce);

  std::map<GR_index_t, std::vector<GRCoEdge*> >::iterator it =
      region_curves.find(region);

  if (it == region_curves.end()) {
    std::vector<GRCoEdge*> ce_vector;
    ce_vector.push_back(ce);
    region_curves.insert(std::make_pair(region, ce_vector));
  }
  else
    it->second.push_back(ce);

}

void
GRGeom2D::build_loop(std::vector<GRCoEdge*>& coedge_vector,
		     std::vector<GRLoop*>& created_loops)
{

  //The input vector contains all the co-edges attached to surface to
  //which we are building loops. Also gives orientation to the co-edges

  unsigned int i;

  //co-edges already added to a loop are true.
  bool added[coedge_vector.size()];
  std::fill(added, &added[coedge_vector.size()], false);
  created_loops.clear();

  while (1) {

    //The created loop
    GRLoop* loop = new GRLoop();

    //List of coedges for this loop;
    DLIList<CoEdgeSM*> ce;

    //Finding the first co-edge not already in a loop. Set the sense
    //to forward and set subsequent ce-edge senses based on this one.
    bool found = false;
    GRPoint* start_pt;
    GRPoint* end_pt;

    for (i = 0; i < coedge_vector.size() && !found; i++) {
      if (!added[i]) {
	found = true;
	added[i] = true;
	coedge_vector[i]->set_sense(CUBIT_FORWARD);
	coedge_vector[i]->add_loop(loop);

	GRCurve* c = dynamic_cast<GRCurve*>(coedge_vector[i]->get_curve());
	start_pt = dynamic_cast<GRPoint*>(c->get_start_point());
	end_pt = dynamic_cast<GRPoint*>(c->get_end_point());
	ce.insert(coedge_vector[i]);
      }
    }

    //All co-edges have been used - all loops are formed: return.
    if (!found) {
      delete loop;
      return;
    }

    //Self-closed curves
    if (start_pt == end_pt) {
      if (ce.size() != 1)
	vFatalError("The co-edge list should only contain one element",
		    "GRGeom2D::build_loops(std::vector<GRCoEdge*>&)");
      loop->add_coedges(ce);
      loop_list.push_back(loop);
      created_loops.push_back(loop);
      continue;
    }

    //Completing loops from our seed co-edge.
    GRPoint* p1 = NULL;
    GRPoint* p2 = NULL;
    while (p1 != start_pt && p2 != start_pt) {

      found = false;
      for (i = 0; i < coedge_vector.size(); i++) {

	if (added[i])
	  continue;

	GRCurve* c = dynamic_cast<GRCurve*>(coedge_vector[i]->get_curve());
	p1 = dynamic_cast<GRPoint*>(c->get_start_point());
	p2 = dynamic_cast<GRPoint*>(c->get_end_point());

	//This is a special case, a self-closed curve connected to
	//another loop. Construct a new loop independently.
	if (p1 == p2) {
	  added[i] = true;
	  DLIList<CoEdgeSM*> ce_list;
	  ce_list.insert(coedge_vector[i]);
	  GRLoop* closed_loop = new GRLoop(ce_list);
	  coedge_vector[i]->set_sense(CUBIT_FORWARD);
	  coedge_vector[i]->add_loop(closed_loop);
	  continue;
	}

	if (p1 == end_pt) {
	  found = true;
	  added[i] = true;
	  end_pt = p2;
	  coedge_vector[i]->set_sense(CUBIT_FORWARD);
	  coedge_vector[i]->add_loop(loop);
	  ce.insert(coedge_vector[i]);
	}
	else if (p2 == end_pt) {
	  found = true;
	  added[i] = true;
	  end_pt = p1;
	  coedge_vector[i]->set_sense(CUBIT_REVERSED);
	  coedge_vector[i]->add_loop(loop);
	  ce.insert(coedge_vector[i]);
	}

	if (found)
	  break;

      }

      if (!found) {
	logMessage(0,
		   "ERROR: Could not built a closed loops with the curves.\n");
	logMessage(0, "Geometry creation failed.\n");
	abort();
      }

    }

    loop->add_coedges(ce);
    loop_list.push_back(loop);
    created_loops.push_back(loop);

  }

}

std::pair<GRCurve *, GRCurve *>
GRGeom2D::offset_curve_linear(GR_index_t curveNumber,
			      const CubitVector & offset)
{

  if (curveNumber >= num_curves())
    return std::pair<GRCurve *, GRCurve *>(static_cast<GRCurve*>(NULL),
					   static_cast<GRCurve*>(NULL));
  std::list<GRCurve*>::iterator itc = curve_list.begin();
  GR_index_t iC = 0;
  GRCurve * curveToOffset = NULL;
  while (iC <= curveNumber) {
    if (iC == curveNumber) {
      curveToOffset = (*itc);
      break;
    }
    iC++, itc++;

  }
  assert(curveToOffset);
  // now we need to create new points and offset

  GRCurveGeom * cgToOffset = curveToOffset->get_curve_geom();
  // need to offset points
  GR_index_t n_pts = cgToOffset->get_num_points();

  curve_info * old_ci = &(v_curves[curveNumber]);
  // create new curve info
  curve_info ci;
  ci.curve_type = old_ci->curve_type;
  ci.l_cond = old_ci->l_cond;
  ci.left = old_ci->left;
  ci.r_cond = old_ci->r_cond;
  ci.right = old_ci->right;
  ci.curve_data.push_back(n_pts);

  for (GR_index_t iP = 0; iP < n_pts; iP++) {
    CubitVector point = cgToOffset->get_coord(iP);
    v_points.push_back(point + offset);
    ci.curve_data.push_back(n_points + iP);

  }
  n_points += n_pts;
  v_curves.push_back(ci);
  assert(v_points.size() == n_points);
  std::map<GR_index_t, GRPoint*> created_points;
  std::map<GR_index_t, std::vector<GRCoEdge*> > region_curves;
  build_curve(v_curves.size() - 1, created_points, region_curves);

  return std::pair<GRCurve *, GRCurve *>(curveToOffset, curve_list.back());
}
std::map<GRCurve *, GRCurve *>
GRGeom2D::offset_curves(std::list<GRCurve*> & curvesToOffset, double dist)
{
  std::map<GRCurve *, GRCurve *> curveMap;
  std::list<GRCurve*>::iterator itc, itc_next, itc_prev;
  std::list<GRCurve*> newCurves;
  // assume we are going to offset ALL points.
  //	assert(n_points == point_list.size());

  //	for(GR_index_t iPoint = 0; iPoint < n_points; iPoint++ ){
  //		pointNormalCount[iPoint] = 0;
  //	}
  //	printf("num points %d %d\n",point_list.size(),v_points.size());
  std::map<GRPoint*, CubitVector> pointNormals;
  for (itc = curvesToOffset.begin(); itc != curvesToOffset.end(); itc++) {
    GRCurve * curveToOffset = *itc;
    //		printf("curve %p ",curveToOffset);

    DLIList<GRPoint*> local_point_list;
    // assert(curveToOffset->get_curve_geom()->get_curve_type()
    //        == GRCurveGeom::LINE);
    curveToOffset->get_points(local_point_list);
    CubitVector normal;
    //		printf("normal %p %f %f %f\n",curveToOffset,normal.x(),normal.y(),normal.z());
    for (int i = 0; i < local_point_list.size(); ++i) {
      GRPoint * local_point = local_point_list.get_and_step();
      //			printf("local pt %p ",local_point);
      curveToOffset->get_curve_geom()->unit_normal(0.0, normal);
      std::map<GRPoint*, CubitVector>::iterator pnit = pointNormals.find(
	  local_point);
      if (pnit == pointNormals.end()) {
	pointNormals.insert(
	    std::pair<GRPoint*, CubitVector>(
		local_point, CubitVector(normal.x(), normal.y(), 3.0)));

      }
      else {
	// cheating way, store angle between them
	double dot = pnit->second % normal;
	//				assert(dot > -1e-10);
	pnit->second = 0.5 * (pnit->second + normal);
	pnit->second.set(pnit->second.x(), pnit->second.y(), dot);
	//				printf("dot %f %f %f %f\n",dot,pnit->second.x(),pnit->second.y(),pnit->second.z());
      }
    }
    //		printf("\n");

  }
  //	printf("point no %d\n",pointNormals.size());
  // mapping by integers, this is going to get tricky
  std::map<GRPoint*, int> pointMap;
  // new points created
  for (std::map<GRPoint*, CubitVector>::iterator it = pointNormals.begin();
      it != pointNormals.end(); it++) {
    //		printf("Curve %p %f\n", it->first,it->second.z());
    CubitVector normal = it->second;
    CubitVector point = CubitVector(it->first->X(), it->first->Y(), 0.0);

    if (normal.z() < 2.0) {
      double projdist = dist / cos(0.5 * acos(normal.z()));
      if (isnanl(projdist))
	projdist = dist;
      //						printf("proj dist %f %f %f\n", dist,projdist,normal.z());
      normal.set(normal.x(), normal.y(), 0.0);
      normal.normalize();
      //							CubitVector tempPoint = point+projdist*normal;
      //							printf("temp %f %f | %f %f \n",point.x(),point.y(), tempPoint.x(),tempPoint.y());
      v_points.push_back(point + projdist * normal);

    }
    else {
      normal.set(normal.x(), normal.y(), 0.0);
      //						CubitVector tempPoint = point+dist*normal;

      //							printf("temp %f %f | %f %f \n",point.x(),point.y(), tempPoint.x(),tempPoint.y());

      v_points.push_back(point + dist * normal);
    }
    pointMap.insert(std::pair<GRPoint*, int>(it->first, v_points.size() - 1));
    n_points++;
  }

  int curveNumber = 0;
  // need to create new points
  std::map<GR_index_t, GRPoint*> created_points;
  std::map<GR_index_t, std::vector<GRCoEdge*> > region_curves;
  for (itc = curvesToOffset.begin(); itc != curvesToOffset.end();
      itc++, curveNumber++) {
    GRCurve * curveToOffset = *itc;
    curve_info * old_ci = &(v_curves[curveNumber]);
    // create new curve info
    curve_info ci;
    ci.curve_type = old_ci->curve_type;
    ci.l_cond = old_ci->l_cond;
    ci.left = old_ci->left;
    ci.r_cond = old_ci->r_cond;
    ci.right = old_ci->right;
    ci.curve_data.push_back(2);
    DLIList<GRPoint*> local_point_list;
    curveToOffset->get_points(local_point_list);

    for (int i = 0; i < local_point_list.size(); ++i) {
      GRPoint * local_point = local_point_list.get_and_step();
      assert(pointMap.find(local_point) != pointMap.end());
      //			printf("int %d ",pointMap.find(local_point)->second);
      ci.curve_data.push_back(pointMap.find(local_point)->second);
    }
    //		printf("\n");
    v_curves.push_back(ci);

    build_curve(v_curves.size() - 1, created_points, region_curves);
    curveMap.insert(
	std::pair<GRCurve *, GRCurve *>(curveToOffset, curve_list.back()));

  }
  return curveMap;
}

void
GRGeom2D::verify_region_conditions(int region,
				   std::vector<GRLoop*>& created_loops)
{

  //At this point, we must make sure that all the co-edges of
  //the loop have the same region on one side.

  if (!loop_list.size()) {
    logMessage(0, "ERROR: No loops were formed.\n");
    logMessage(0, "Must at least have one loop to create a mesh.\n");
    logMessage(0, "Geometry creation failed.\n");
    abort();
  }

  std::vector<GRLoop*>::iterator itl;
  int largest_loop = 0;

  for (itl = created_loops.begin(); itl != created_loops.end(); itl++) {

    int condition_side = 0;

    DLIList<GRCoEdge*> ce_list;
    (*itl)->get_coedges(ce_list);
    ce_list.reset();

    GRCoEdge* ce = ce_list.next(0);
    GRCurve* c = dynamic_cast<GRCurve*>(ce->get_curve());
    assert(c);

    int lr = c->region_left(), rr = c->region_right();
    largest_loop = std::max(largest_loop, ce_list.size());

    switch (ce->sense())
      {
      case CUBIT_FORWARD:
	{
	  if (lr == region && rr == region) {
	    condition_side = 3;
	    break;
	  }
	  if (rr == region) {
	    condition_side = 2;
	    break;
	  }
	  if (lr == region) {
	    condition_side = 1;
	    break;
	  }
	  assert(0);
	  break;
	}
      case CUBIT_REVERSED:
	{
	  if (lr == region && rr == region) {
	    condition_side = 3;
	    break;
	  }
	  if (lr == region) {
	    condition_side = 2;
	    break;
	  }
	  if (rr == region) {
	    condition_side = 1;
	    break;
	  }
	  assert(0);
	  break;
	}
      default:
	{
	  assert(0);
	  break;
	}
      }

    for (int i = 1; i < ce_list.size(); i++) {

      ce = ce_list.next(i);
      c = dynamic_cast<GRCurve*>(ce->get_curve());
      assert(c);

      lr = c->region_left();
      rr = c->region_right();

      switch (condition_side)
	{

	case 1:
	  {
	    switch (ce->sense())
	      {
	      case CUBIT_FORWARD:
		{
		  if (lr != region) {
		    logMessage(
			0,
			"ERROR: Region %d has an incorrect curve assignment.\n",
			region);
		    logMessage(0, "Geometry creation failed.\n");
		    abort();
		  }
		  break;
		}
	      case CUBIT_REVERSED:
		{
		  if (rr != region) {
		    logMessage(
			0,
			"ERROR: Region %d has an incorrect curve assignment.\n",
			region);
		    logMessage(0, "Geometry creation failed.\n");
		    abort();
		  }
		  break;
		}
	      default:
		{
		  assert(0);
		  break;
		}
	      }
	    break;
	  }

	case 2:
	  {
	    switch (ce->sense())
	      {
	      case CUBIT_FORWARD:
		{
		  if (rr != region) {
		    logMessage(
			0,
			"ERROR: Region %d has an incorrect curve assignment.\n",
			region);
		    logMessage(0, "Geometry creation failed.\n");
		    abort();
		  }
		  break;
		}
	      case CUBIT_REVERSED:
		{
		  if (lr != region) {
		    logMessage(
			0,
			"ERROR: Region %d has an incorrect curve assignment.\n",
			region);
		    logMessage(0, "Geometry creation failed.\n");
		    abort();
		  }
		  break;
		}
	      default:
		{
		  assert(0);
		  break;
		}
	      }
	    break;
	  }

	case 3:
	  {
	    if (lr != region || rr != region) {
	      logMessage(
		  0,
		  "ERROR: This curve should have a region on both sides - OR - \n");
	      logMessage(
		  0, "a curve has a region on both sides and should not.\n");
	      logMessage(0, "Geometry creation failed.\n");
	      abort();
	    }
	    break;
	  }

	default:
	  {
	    assert(0);
	    break;
	  }
	}

    }

  }

  //   if(largest_loop < 2) {
  //     vMessage(0, "ERROR: No loops of region %d can define a closed region.\n", region);
  //     vMessage(0, "Geometry creation failed.\n");
  //     abort();
  //   }

}

void
GRGeom2D::writeGeometryVTK(const char strBaseFileName[])
{
  GR_index_t num_inter_points = 101;
  FILE *pFOutFile = NULL;

  std::list<GRCurve*> curves;
  get_curves(curves);
  GR_index_t numCurves = curves.size();
  char strFileName[1024];
  sprintf(strFileName, "%s.vtk", strBaseFileName);
  printf("Writing Geometry File %s\n", strFileName);

  pFOutFile = fopen(strFileName, "w");
  //	pFOutFileText = fopen(strFileNameText,"w");
  if (NULL == pFOutFile)
    vFatalError("Couldn't open Vtk output file for writing", "2d mesh output");
  // Write the VTK header details
  //-------------------------------------------------------
  fprintf(pFOutFile, "# vtk DataFile Version 1.0");
  fprintf(pFOutFile, "\nOVMesh Tri example");
  fprintf(pFOutFile, "\nASCII\n");
  fprintf(pFOutFile, "\nDATASET UNSTRUCTURED_GRID");
  fprintf(pFOutFile, "\nPOINTS %u float", num_inter_points * numCurves);
  double minLength = LARGE_DBL;
  //-------------------------------------
  // write 2d functional data
  //-------------------------------------
  std::list<GRCurve*>::iterator itc;
  for (itc = curves.begin(); itc != curves.end(); itc++) {
    GRCurve * newCurve = (*itc);
    GRCurveGeom * cg = newCurve->get_curve_geom();
    minLength = min(cg->arc_length(), minLength);
    double dp = cg->max_param() / (num_inter_points - 1);
    CubitVector vec;
    for (GR_index_t i = 0; i < num_inter_points; ++i) {
      cg->coord_at_param(dp * i, vec);
      fprintf(pFOutFile, "\n%16.8f %16.8f %16.8f", vec.x(), vec.y(), 0.0);
      //			fprintf(pFOutFileText, "%16.8f %16.8f\n",vec.x(),vec.y());
    }
  }
  //	fclose(pFOutFileText);
  // Treat the cubic as a single cell, type 4 (connected line?)

  fprintf(pFOutFile, "\n\nCELLS %u %u\n", numCurves,
	  (1 + num_inter_points) * numCurves);
  for (GR_index_t ic = 0; ic < numCurves; ic++) {

    fprintf(pFOutFile, "%u ", num_inter_points);
    for (GR_index_t i = ic * num_inter_points; i < num_inter_points * (ic + 1);
	++i) {
      fprintf(pFOutFile, " %10d", i);
    }
    fprintf(pFOutFile, "\n");
  }
  fprintf(pFOutFile, "\nCELL_TYPES %u\n", numCurves);
  for (GR_index_t ic = 0; ic < numCurves; ic++) {

    fprintf(pFOutFile, "4\n");
  }
  fprintf(pFOutFile, "\n");
  fclose(pFOutFile);
  printf("Minimum Boundary Input Length %.16f\n", minLength);
}

