#include "GR_GraphicsOut.h"
#include "GeometryQueryTool.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "GMem.hpp"
#include "CubitVector.hpp"
#include "DLIList.hpp"

#include <utility>
#include <map>
#include <iterator>

void
GraphicsOut::output_edges() const
{

  assert(output_filename.size() != 0);

  GeometryQueryTool* gqt = GeometryQueryTool::instance();

  DLIList<RefEdge*> edge_list;
  gqt->ref_edges(edge_list);

  int num_edges = edge_list.size();
  edge_list.reset();

  std::vector<CubitVector> all_points;
  std::vector<std::pair<int, int> > all_edges;

  for (int i = 0; i < num_edges; i++) {
    RefEdge* edge = edge_list.next(i);
    GMem gmem;
    edge->get_graphics(gmem);

    int num_points = gmem.pointListCount;
    int num_facets = gmem.fListCount;

    GPoint* point_array = gmem.point_list();
    int* facet_array = gmem.facet_list();

    std::map<int, int> local_to_global;
    std::map<int, int>::iterator it;

    for (int j = 0; j < num_facets;) {
      int num_pts = facet_array[j++];
      assert(num_pts == 2);
      int local[2] =
	{ -1 };
      int global[2] =
	{ -1 };
      local[0] = facet_array[j++];
      local[1] = facet_array[j++];
      assert(local[0] > -1 && local[0] < num_points);
      assert(local[1] > -1 && local[1] < num_points);

      for (int k = 0; k < 2; k++) {
	it = local_to_global.find(local[k]);
	if (it == local_to_global.end()) {
	  CubitVector coord(point_array[local[k]].x, point_array[local[k]].y,
			    point_array[local[k]].z);
	  std::back_inserter(all_points) = coord;
	  global[k] = static_cast<int>(all_points.size()) - 1;
	  local_to_global.insert(std::make_pair(local[k], global[k]));
	}
	else {
	  global[k] = it->second;
	}
      }

      std::back_inserter(all_edges) = std::make_pair(global[0] + 1,
						     global[1] + 1);
    }

  }

  FILE* out_file = fopen(output_filename.c_str(), "w");
  fprintf(out_file, "# vtk DataFile Version 1.0\n");
  fprintf(out_file, "Geom edge data\n");
  fprintf(out_file, "ASCII\n");
  fprintf(out_file, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(out_file, "POINTS %d float\n", static_cast<int>(all_points.size()));
  for (unsigned int i = 0; i < all_points.size(); i++) {
    fprintf(out_file, "%lf %lf %lf\n", all_points[i].x(), all_points[i].y(),
	    all_points[i].z());
  }

  unsigned int nEdges = all_edges.size();
  fprintf(out_file, "Cells %d %d\n", nEdges, nEdges * 3);
  for (unsigned int i = 0; i < nEdges; i++) {
    assert(all_edges[i].first > 0);
    assert(all_edges[i].second > 0);
    fprintf(out_file, "2 %d %d\n", all_edges[i].first - 1,
	    all_edges[i].second - 1);
  }

  fprintf(out_file, "CELL_TYPES %d\n", nEdges);
  for (unsigned int i = 0; i < nEdges; i++) {
    fprintf(out_file, "3\n");
  }
  fclose(out_file);

}

void
GraphicsOut::output_faces() const
{

  assert(output_filename.size() != 0);

  GeometryQueryTool* gqt = GeometryQueryTool::instance();

  DLIList<RefFace*> face_list;
  gqt->ref_faces(face_list);

  int num_faces = face_list.size();
  face_list.reset();

  std::vector<CubitVector> all_points;
  std::vector<triplet> all_faces;

  for (int i = 0; i < num_faces; i++) {
    RefFace* face = face_list.next(i);
    GMem gmem;
    face->get_graphics(gmem);

    int num_points = gmem.pointListCount;
    int num_facets = gmem.fListCount;

    GPoint* point_array = gmem.point_list();
    int* facet_array = gmem.facet_list();

    std::map<int, int> local_to_global;
    std::map<int, int>::iterator it;

    for (int j = 0; j < num_facets;) {
      int num_pts = facet_array[j++];
      assert(num_pts == 3);
      int local[3] =
	{ -1 };
      int global[3] =
	{ -1 };
      local[0] = facet_array[j++];
      local[1] = facet_array[j++];
      local[2] = facet_array[j++];
      assert(local[0] > -1 && local[0] < num_points);
      assert(local[1] > -1 && local[1] < num_points);
      assert(local[2] > -1 && local[2] < num_points);

      for (int k = 0; k < 3; k++) {
	it = local_to_global.find(local[k]);
	if (it == local_to_global.end()) {
	  CubitVector coord(point_array[local[k]].x, point_array[local[k]].y,
			    point_array[local[k]].z);
	  std::back_inserter(all_points) = coord;
	  global[k] = static_cast<int>(all_points.size()) - 1;
	  local_to_global.insert(std::make_pair(local[k], global[k]));
	}
	else {
	  global[k] = it->second;
	}
      }

      triplet t;
      t.ind1 = global[0] + 1;
      t.ind2 = global[1] + 1;
      t.ind3 = global[2] + 1;
      std::back_inserter(all_faces) = t;
    }

  }

  FILE* out_file = fopen(output_filename.c_str(), "w");
  fprintf(out_file, "# vtk DataFile Version 1.0\n");
  fprintf(out_file, "Geom tri data\n");
  fprintf(out_file, "ASCII\n");
  fprintf(out_file, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(out_file, "POINTS %d float\n", static_cast<int>(all_points.size()));
  for (unsigned int i = 0; i < all_points.size(); i++) {
    fprintf(out_file, "%lf %lf %lf\n", all_points[i].x(), all_points[i].y(),
	    all_points[i].z());
  }

  unsigned int nTri = all_faces.size();
  fprintf(out_file, "Cells %d %d\n", nTri, nTri * 4);
  for (unsigned int i = 0; i < nTri; i++) {
    assert(all_faces[i].ind1 > 0);
    assert(all_faces[i].ind2 > 0);
    assert(all_faces[i].ind3 > 0);
    fprintf(out_file, "3 %d %d %d\n", all_faces[i].ind1 - 1,
	    all_faces[i].ind2 - 1, all_faces[i].ind3 - 1);
  }

  fprintf(out_file, "CELL_TYPES %d\n", nTri);
  for (unsigned int i = 0; i < nTri; i++) {
    fprintf(out_file, "5\n");
  }
  fclose(out_file);

}
