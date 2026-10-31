#include "GR_RefinementManager3D.h"

#include "GR_config.h"
#include "GR_BaseQueueEntry.h"
#include "GR_BFace.h"
#include "GR_BFace3D.h"
#include "GR_Cell.h"
#include "GR_Geometry.h"
#include "GR_InsertionManager.h"
#include "GR_SmoothingManager.h"
#include "GR_Mesh.h"

#include "GR_Subseg.h"
#include "GR_SubsegManager.h"
#include "GR_Vertex.h"
#include "GR_VolMesh.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

using std::deque;
using std::map;
using std::multimap;
using std::pair;
using std::queue;
using std::set;
using std::vector;

OldPriorityCalculator::~OldPriorityCalculator()
{
}

void
OldPriorityCalculator::set_length_calculator(
    std::shared_ptr<GRUMMP::Length> const length)
{

  assert(length);
  m_length = length;

}

void
OldPriorityCalculator::compute_priority(const Cell* const cell,
					double& length_prio,
					double& shape_prio) const
{

  //Start by computing the size priority.
  length_prio = compute_length_priority(cell);

  //If length_prio == 0, then we will queue based on
  //shape quality (i.e. size is correct). Otherwise, cell is either
  //too big (length_prio > 0.), or too small (this should only
  //happens around small angles, length_prio < 0.).

  if (iFuzzyComp(length_prio, 0.) != 0) {
    shape_prio = -LARGE_DBL;
    return;
  }
  shape_prio = compute_shape_priority(cell);

}
double
OldPriorityCalculator::compute_shape_priority(const Cell * const cell) const
{
  switch (m_type)
    {

    case EDGE_RADIUS:

      //The inverse of the classical circumradius-to-shortest edge ratio.
      //Use the inverse so that the lowest value will be the worst quality (not the inverse).
      //Leaves slivers behind, but comes with a termination guarantee.

      //For a regular tetrahedron, the ratio is 1. / sqrt(3./8.),
      //which is the maximum achievable value (and best possible quality).
      //The value returned by the function is scaled by this number
      //so that it is between 0 and 1.

      return edge_radius_ratio(cell);

    case SIN_DIHED:

      //The minimum sine of a tet's six dihedreal angles.
      //Leaves behind needle tets.

      //Obviously, returns values between 0 and 1, so no need to normalize.
      //Best possible quality is sin(arccos(1./3.)).

      return min_sin_dihed(cell);

    case BIAS_SIN_DIHED:

      //A biased version of the minimum sine of a tet's dihedral angles neasure.
      //Exactly like the minimum sine measure, but if a dihedral angle is obtuse,
      //multiply its sine by 0.7, therefore penalizing large dihedral angles
      //without affecting small diheds.

      //Values are between 0 and 1.

      return bias_min_sin_dihed(cell);

    case VOLUME_LENGTH:

      //The volume-length measure: (6*sqrt(2)) * Volume / (length(rms))^3
      //Proposed by Parthasarathy, Graichen and Hathaway.
      //Values are between 0 and 1.

      return volume_length(cell);

    case NORMALIZED_SHAPE_RATIO:

      //This is basically the inradius to circumradius ratio
      //with a normalization factor.
      //Returns 3 * inradius / circumradius. Scaled this way, the ratio
      //is always between 0 and 1 ( == 1 for a regular tet).
      //Proposed by Cavendish, Field and Frey...
      //or is it David A. Field, A generic Delaunay triangulation algorithm
      //for finite element meshes. (1991)

      return normalized_shape_ratio(cell);

    case MIN_SOLID:

      //The minimum solid angle of a tetrahedron.
      //A regular tet has a min solid angle of 3 * acos(1./3.) - M_PI
      //and the value returned is scaled with this ideal value so that the
      //return value is between 0 an 1.

      return min_solid(cell);

    case MEAN_RATIO:

      //The mean ratio as defined by Liu and Joe.
      //Value is between 0 and 1.

      return mean_ratio(cell);

    case VOLUME_AREA_RATIO:

      //A measure that relates the volume of the tet
      //to the area of its faces. Proposed by
      //Cougny, Shephard and Georges, Explicit node point
      //smoothing within the octree mesh generator (1990)

      return volume_area_ratio(cell);

    case INRAD_ALTITUDE:

      //Ratio of the inradius to the longest altitude. Might
      //be a bit costly to compute.

      return inrad_altitude_ratio(cell);

    case INRAD_LENGTH:

      //Ratio of the inradius to the longest edge.

      return inrad_longest_edge_ratio(cell);

    case SIZE_ONLY:

    	// No shape quality calculated.
    	return 0;

    default:
      vFatalError("Invalid priority type",
		  "PriorityCalculator::compute_priority");
      return 0;

    }

}

double
OldPriorityCalculator::compute_length_priority(const Cell* const cell) const
{
  static const double unqueueablePriority = -LARGE_DBL;
  if (!m_length)
    return 0.;

  assert(cell && cell->isValid() && !cell->isDeleted());
  if (cell->getType() != Cell::eTet)
    return unqueueablePriority;
  assert(cell->getNumVerts() == 4);

  double current_size = (2. / sqrt(3.)) * cell->calcCircumradius();
  double target_size = 0.;

  const Vert* vert;
  double vert_ls;

  for (int ii = 0; ii < 4; ++ii) {

    vert = cell->getVert(ii);
    //If vert is at the apex of a small angle, bail out right now.
    if (vert->isSmallAngleVert()) {
      assert(vert->getVertType() == Vert::eBdryApex);
      return unqueueablePriority;
    }

    vert_ls = vert->getLengthScale();
    assert(vert_ls > 0.);

    target_size += vert_ls;

  }

  //divide it by the number of vertices.
  target_size *= 0.25;

  if (current_size > target_size)
    return current_size / target_size;

  else if (current_size < target_size / 15.)
    //This hack to prevent cascading refinement
    //near small angle in the input.
    return unqueueablePriority;

  else
    //Size is correct, returns 0. so we do not queue solely based on size.
    return 0.;

}

double
OldPriorityCalculator::edge_radius_ratio(const Cell* const cell)
{

  //Scale the ratio by sqrt(3/8) so that it is between 0 and 1.

  assert(cell->getType() == Cell::eTet);

  return (sqrt(0.375) * cell->calcShortestEdgeLength()
      / cell->calcCircumradius());

}

double
OldPriorityCalculator::min_sin_dihed(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);
  double all_diheds[6], mini_sin_dihed = LARGE_DBL;
  int i, num_diheds;

  //Get the dihedral angles in radians.
  cell->calcAllDihed(all_diheds, &num_diheds, false);

  for (i = 0; i < num_diheds; ++i)
    mini_sin_dihed = std::min(mini_sin_dihed, sin(all_diheds[i]));

  return mini_sin_dihed;

}

double
OldPriorityCalculator::bias_min_sin_dihed(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);
  double all_diheds[6], mini_sin_dihed = LARGE_DBL;
  int i, num_diheds;

  //Get the dihedral angles in radians.
  cell->calcAllDihed(all_diheds, &num_diheds, false);

  for (i = 0; i < num_diheds; ++i) {
    if (all_diheds[i] > 90.)
      mini_sin_dihed = std::min(mini_sin_dihed, 0.7 * sin(all_diheds[i]));
    else
      mini_sin_dihed = std::min(mini_sin_dihed, sin(all_diheds[i]));
  }

  return mini_sin_dihed;

}

double
OldPriorityCalculator::volume_length(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);

  double edge0[] =
  adDIFF3D(cell->getVert(0)->getCoords(),
      cell->getVert(1)->getCoords()), edge1[] =
  adDIFF3D(cell->getVert(0)->getCoords(),
      cell->getVert(2)->getCoords()), edge2[] =
  adDIFF3D(cell->getVert(0)->getCoords(),
      cell->getVert(3)->getCoords()), edge3[] =
  adDIFF3D(cell->getVert(1)->getCoords(),
      cell->getVert(2)->getCoords()), edge4[] =
  adDIFF3D(cell->getVert(1)->getCoords(),
      cell->getVert(3)->getCoords()), edge5[] =
  adDIFF3D(cell->getVert(2)->getCoords(),
      cell->getVert(3)->getCoords());

  double length = sqrt(
      (1. / 6.)
	  * (dMAG3D_SQ(edge0) + dMAG3D_SQ(edge1) + dMAG3D_SQ(edge2)
	      + dMAG3D_SQ(edge3) + dMAG3D_SQ(edge4) + dMAG3D_SQ(edge5)));

  return ((6. * sqrt(2.) * cell->calcSize()) / (length * length * length));

}

double
OldPriorityCalculator::min_solid(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);
  double all_solids[4], mini_solid = LARGE_DBL;
  int i, num_solids;

  //Get the solid angles in steradians.
  cell->calcAllSolid(all_solids, &num_solids, false);

  for (i = 0; i < num_solids; ++i)
    mini_solid = std::min(mini_solid, all_solids[i]);

  //Scale by 3 * GR_acos(1./3.) - M_PI, the min solid for a regular tet.
  double scale = 3 * GR_acos(1. / 3.) - M_PI;

  return (mini_solid / scale);

}

double
OldPriorityCalculator::normalized_shape_ratio(const Cell* const cell)
{

  //Returns 3 * inradius / circumradius. Scaled this way, the ratio
  //is always between 0 and 1 ( == 1 for a regular tet).
  //Proposed by Cavendish, Field and Frey.

  assert(cell->getType() == Cell::eTet);
  return 3. * cell->calcInradius() / cell->calcCircumradius();

}

double
OldPriorityCalculator::mean_ratio(const Cell* const cell)
{

  //This measure is proposed in A. Liu and B. Joe (1994)

  assert(cell->getType() == Cell::eTet);
  double volume = cell->calcSize();
  double sum_length_squared = 0.;
  const Vert *v0, *v1;

  for (int i = 0; i < 6; ++i) {

    switch (i)
      {
      case 0:
	v0 = cell->getVert(0);
	v1 = cell->getVert(1);
	break;
      case 1:
	v0 = cell->getVert(0);
	v1 = cell->getVert(2);
	break;
      case 2:
	v0 = cell->getVert(0);
	v1 = cell->getVert(3);
	break;
      case 3:
	v0 = cell->getVert(1);
	v1 = cell->getVert(2);
	break;
      case 4:
	v0 = cell->getVert(1);
	v1 = cell->getVert(3);
	break;
      case 5:
	v0 = cell->getVert(2);
	v1 = cell->getVert(3);
	break;
      }

    sum_length_squared += dDIST_SQ_3D(v0->getCoords(), v1->getCoords());

  }

  return (12. * (pow(3. * volume, 2. / 3.)) / sum_length_squared);

}

double
OldPriorityCalculator::volume_area_ratio(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);

  double volume = cell->calcSize();
  double sum_area_sq = 0., area;

  const Face* face;

  for (int i = 0; i < 4; ++i) {
    face = cell->getFace(i);
    assert(face->isValid() && !face->isDeleted());
    area = face->calcSize();
    sum_area_sq += area * area;
  }

  return (2187. * (volume * volume * volume * volume)
      / (sum_area_sq * sum_area_sq * sum_area_sq));

}

double
OldPriorityCalculator::inrad_altitude_ratio(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);

  double minArea = LARGE_DBL;

  for (int i = 0; i < 4; ++i) {

    const Face* face = cell->getFace(i);
    assert(face->isValid() && !face->isDeleted());
    double area = face->calcSize();
    minArea = std::min(area, minArea);
  }

  double volume = cell->calcSize();
  double altitude = volume / (3 * minArea);
  return (4. * cell->calcInradius() / altitude);

}

double
OldPriorityCalculator::inrad_longest_edge_ratio(const Cell* const cell)
{

  assert(cell->getType() == Cell::eTet);

  double long_edge = -LARGE_DBL;
  const Vert *v0, *v1;

  for (int i = 0; i < 6; ++i) {

    switch (i)
      {
      case 0:
	v0 = cell->getVert(0);
	v1 = cell->getVert(1);
	break;
      case 1:
	v0 = cell->getVert(0);
	v1 = cell->getVert(2);
	break;
      case 2:
	v0 = cell->getVert(0);
	v1 = cell->getVert(3);
	break;
      case 3:
	v0 = cell->getVert(1);
	v1 = cell->getVert(2);
	break;
      case 4:
	v0 = cell->getVert(1);
	v1 = cell->getVert(3);
	break;
      case 5:
	v0 = cell->getVert(2);
	v1 = cell->getVert(3);
	break;
      }

    long_edge = std::max(long_edge,
			 dDIST_SQ_3D(v0->getCoords(), v1->getCoords()));

  }

  return (2. * sqrt(6.) * cell->calcInradius() / sqrt(long_edge));

}

RefinementQueue::~RefinementQueue()
{

  if (m_calculator)
    delete m_calculator;
  empty_queue();

}

Cell*
RefinementQueue::get_top_valid_cell()
{

  if (empty())
    return static_cast<Cell*>(NULL);

  while (top()->isDeleted()) {
    if (!pop_top_entry())
      return static_cast<Cell*>(NULL);
  }

  assert(!top()->get_cell()->isDeleted());
  return top()->get_cell();

}

bool
RefinementQueue::pop_top_entry()
{

  if (!empty()) {
    CellQueueEntry* dummy = top();
    pop();
    if (dummy)
      delete dummy;
    return true;
  }

  return false;

}
void
RefinementQueue::add_quality_entry(Cell* const cell)
{
  if (!cell->isValid() || cell->isDeleted())
    return;
  double shape_prio = m_calculator->compute_shape_priority(cell);
  //If cell is within the acceptable shape priority range,
  //we do not want to queue it.
  if (shape_prio < m_min_shape_priority || shape_prio > m_max_shape_priority)
    return;

#ifndef NDEBUG
  //I added this just to make sure the minimal length priority
  //is always bigger than the maximum shape priority. (see a few
  //lines above where I add 1.e4 to the length_prio when queueing).
  //This is not absolutely necessary, just more like a sanity check.
  if (comp.max_on_top())
    assert(shape_prio < 1.e4);
#endif

  push(new CellQueueEntry(cell, shape_prio));
}
void
RefinementQueue::add_entry(Cell* const cell)
{

  if (!cell->isValid() || cell->isDeleted())
    return;

  double length_prio, shape_prio;
  m_calculator->compute_priority(cell, length_prio, shape_prio);

  //We do something different depending on the value of length_prio:
  // - If negative: Never queue the cell (probably in a small angle corner).
  // - If zero: Cell is small enough w/r to sizing field.
  // - If positive: Cell is too big, need to be split regardless of shape.
  int length_comp = iFuzzyComp(length_prio, 0.);

  if (length_comp == 1) {
    if (comp.max_on_top())
      push(new CellQueueEntry(cell, 1.e4 + length_prio));
    else
      push(new CellQueueEntry(cell, -length_prio));
    return;
  }
  else if (length_comp == 0)
    assert(iFuzzyComp(shape_prio, 0.) != -1);
  else
    return; //if the length_prio is negative, we do not want to queue this cell.
  //If cell is within the acceptable shape priority range,
  //we do not want to queue it.
  if (shape_prio <= m_min_shape_priority || shape_prio >= m_max_shape_priority)
    return;

#ifndef NDEBUG
  //I added this just to make sure the minimal length priority
  //is always bigger than the maximum shape priority. (see a few
  //lines above where I add 1.e4 to the length_prio when queueing).
  //This is not absolutely necessary, just more like a sanity check.
  if (comp.max_on_top())
    assert(shape_prio < 1.e4);
#endif

  push(new CellQueueEntry(cell, shape_prio));

}
void
RefinementQueue::clean_queue()
{

  //Removes the deleted entities from the priority queue.

  if (empty())
    return;

  CellQueueEntry* entry;
  EntryVec valid, deleted;
  EntryVec::iterator it = c.begin(), it_end = c.end();

  GR_index_t num_total = size();
  valid.reserve(num_total);
  deleted.reserve(num_total);

  for (; it != it_end; ++it) {
    entry = *it;

    if (entry->isDeleted())
      deleted.push_back(entry);
    else
      valid.push_back(entry);
  }

#ifndef NDEBUG
  GR_index_t num_valid = valid.size(), num_deleted = deleted.size();
  assert(num_total == num_valid + num_deleted);
#endif

  //Put the valid entries in the queue container and sort.
  c.swap(valid);

  //	std::make_heap(c.begin(), c.end(), CellQueueEntry::PrioOrder());
  //	print_queue();
  //Delete the deleted entries.
  for (EntryVec::iterator it2 = deleted.begin(); it2 != deleted.end(); ++it2) {
    if ((*it2))
      delete (*it2);
  }
#ifndef NDEBUG
  assert(size() == num_valid);
#endif

  //print_queue();
}

void
RefinementQueue::print_queue()
{

  //Prints the queue, but must empty it while doing so.
  //Calling this function aborts the program run.

  CellQueueEntry* entry;
  logMessage(2, "Min Queue Prio %f Max %f\n", m_min_shape_priority,
	     m_max_shape_priority);
  do {
    entry = top();
    pop_top_entry();
    logMessage(2, "cell = %p, priority = %.16e, deleted = %d\n",
	       entry->get_cell(), entry->get_priority(), entry->isDeleted());
    for (int iV = 0; iV < entry->get_cell()->getNumVerts(); iV++) {
      entry->get_cell()->getVert(iV)->printVertInfo();
    }

  }
  while (!empty());

  vFatalError("Printing the priority queue emptied it, cannot do anything else",
	      "RefinementQueue::print_queue()");

}

////////////////////////////////////////////////////////////////
///// Declaration of TetMeshRefiner functions starts here //////
////////////////////////////////////////////////////////////////
namespace GRUMMP
{

  RefinementManager3D::RefinementManager3D(
      VolMesh* const mesh, InsertionPointCalculator * const IPC,
      const std::set<Subseg*>* const subsegs, const eEncroachType encroach_type,
      std::shared_ptr<GRUMMP::Length> length3D,
      const OldPriorityCalculator::PrioType prio_type) :
      m_mesh(mesh), m_IPC(IPC), m_Del3D(NULL), m_DelSwap(NULL), m_SI3D(NULL), m_isoLength(
	  length3D), m_refineForQuality(false), m_subseg_queue(), m_bdryFacesToSplit(), m_cell_queue(
      NULL), m_refineBdryForSize(false), m_lockedCells(), m_lockedBdryFaces(), m_lockedSubsegs()
  {

    assert(m_mesh);
    assert(m_mesh->isValid());
    m_mesh->addObserver(
	this,
	GRUMMP::Observable::cellCreated | GRUMMP::Observable::bfaceCreated
	    | GRUMMP::Observable::cellDeleted);
    // Initialize the WatsonInserter
    m_WI3D = new GRUMMP::WatsonInserter3D(mesh);
    m_Del3D = new GRUMMP::DelaunaySwapDecider3D(true);
    m_DelSwap = new GRUMMP::SwapManager3D(m_Del3D, m_mesh);
    m_SI3D = new GRUMMP::SwappingInserter3D(m_mesh, m_DelSwap);
    m_SI3D->setForcedInsertion(true);

    //Only support ball encroachment for now. In fact, I am not sure I am going
    //to implement lens encroachment. It complicates things enormously
    //(deleting the vertices inside a subsegment's ball might lead to a mesh
    //that is not constrained Delaunay causing an algorithm's failure). In practice,
    //I believe the same quality bounds can be reached, maybe with a few more insertions.
    assert(encroach_type == eBall);

    //Initializing some of the mesh parameter.
    m_mesh->setSwapType(eDelaunay);
    m_mesh->setEncroachmentType(encroach_type);

    //Initialize the subsegment map.
    m_mesh->initializeSubsegs(subsegs);

    //Create the priority queue object.
    //Set the ordering as well as the bounds based on the priority type.
    bool max_prio_on_top;
    double min_shape_prio, max_shape_prio;

    switch (prio_type)
      {

      case OldPriorityCalculator::EDGE_RADIUS:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 3. * sqrt(2.) / 8.;
	break;
      case OldPriorityCalculator::SIN_DIHED:
      case OldPriorityCalculator::BIAS_SIN_DIHED:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = sin((M_PI / 180.) * 14.);
	break;
      case OldPriorityCalculator::MIN_SOLID:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.2;
	break;
      case OldPriorityCalculator::VOLUME_LENGTH:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.35;
	break;
      case OldPriorityCalculator::NORMALIZED_SHAPE_RATIO:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.3;
	break;
      case OldPriorityCalculator::MEAN_RATIO:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.45;
	break;
      case OldPriorityCalculator::VOLUME_AREA_RATIO:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.02;
	break;
      case OldPriorityCalculator::INRAD_ALTITUDE:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.6;
	break;
      case OldPriorityCalculator::INRAD_LENGTH:
	max_prio_on_top = false;
	min_shape_prio = 0.;
	max_shape_prio = 0.35;
	break;
      case OldPriorityCalculator::SIZE_ONLY:
    	  max_prio_on_top = false;
    	  min_shape_prio = 0;
    	  max_shape_prio = 0;
	break;
      default:
	assert(0);
	break;

      }

    //Create an empty queue.
    m_cell_queue = new RefinementQueue(prio_type, max_prio_on_top,
				       min_shape_prio, max_shape_prio);
    //Pass the length scale pointer to the refinement queue.
    m_cell_queue->set_length_calculator(m_isoLength);

    //************ Parallel variables *************
    NumCreatedCells = 0;
    OctT = new OctTree(m_mesh);
    IsInParallel = false;
    lck = new omp_lock_t[NUM_PROCS];
    CheckLock = new bool[NUM_PROCS];
    OKTOGO = false;
    for (int ID = 0; ID < NUM_PROCS; ID++) {
      ThreadLeaves[ID] = NULL;
      Para_m_WI3D[ID] = new GRUMMP::WatsonInserter3D(mesh);
      Para_m_cell_queue[ID] = new RefinementQueue(prio_type, max_prio_on_top,
						  min_shape_prio,
						  max_shape_prio);
      Para_m_SI3D[ID] = new GRUMMP::SwappingInserter3D(m_mesh, m_DelSwap);
      Para_m_SI3D[ID]->setForcedInsertion(true);
//            Para_m_length[ID] = new Length3D(*m_length);
//
//            Para_m_length[ID]->set_subseg_map(m_mesh->getSubsegMap());
//            Para_m_length[ID]->initLength(m_mesh);
//            Para_m_cell_queue[ID]->set_length_calculator(Para_m_length[ID]);
      Para_m_cell_queue[ID]->set_length_calculator(m_isoLength);
      omp_init_lock(&lck[ID]);
      CheckLock[ID] = false;
      ThreadTime[ID][99] = 0;
      ThreadInNum[ID][99] = 0;
    }
    TotalRefineTime = 0;
    JustParallelTime = 0;
    SendEventsTime = 0;
    SendCCTime = 0;
    SendCBTime = 0;
    CumWRTime = 0;
    SecondPassesTime = 0;
    FIndSetsTime = 0;
    //*********************************************

    //We are done. The queues will be filled when calling the refine() function.

  }

  RefinementManager3D::~RefinementManager3D()
  {

    //Deleted m_length and m_cell_queue, that's it...
    //	if(m_length) delete m_length;  // This is just a copy; don't delete it.
    if (m_cell_queue)
      delete m_cell_queue;
    if (m_SI3D)
      delete m_SI3D;
    if (m_DelSwap)
      delete m_DelSwap;
    if (m_Del3D)
      delete m_Del3D;
    if (m_WI3D)
      delete m_WI3D;
    if (m_mesh)
      m_mesh->removeObserver(this);

    delete OctT;
    for (int i = 0; i < NUM_PROCS; i++) {
      omp_destroy_lock(&lck[i]);
      delete Para_m_WI3D[i];
      delete Para_m_SI3D[i];
      delete Para_m_cell_queue[i];
//            delete Para_m_length[i];
    }
    delete[] lck;
    delete[] CheckLock;
  }

//This function initializes the subsegment map. It acts differently depending
//on which geometry type is used. For GRUMMP's native geometry, no subsegments
//are given, therefore they must be infered from the geometry and the mesh that
//is passed to the TetMeshRefiner. Unfortunately, because of that, some topological
//configurations, such as scratch curves or dangling edges cannot be supported. A
//subsegment can only be recognized along the boundary of two surfaces.

//For other geometry types, subsegments are built during the surface sampling phase
//and should be explicitely available. This is the easy case.

  void
  RefinementManager3D::initializeSubsegs(
      const set<Subseg*>* const given_subsegs)
  {

    //This has been moved to a VolMesh member function - this version should never be called.
    assert(0);

    logMessage(MSG_MANAGER, "Initialize subsegs\n");
    assert(m_mesh);
    if (!m_mesh->isSubsegMapEmpty())
      return;
    assert(m_mesh->isSubsegMapEmpty());

    if (given_subsegs) {
      //The case where subsegments are given. This is easy...
      //Note that only the pointer gets added to m_subseg_map (not a copy of the object).
      //The memory gets de-allocated only when the subsegment gets deleted from m_mesh.
      //Be careful with that... (should maybe create a local copy if the subsegs
      //in the map must be used somewhere else, in another object).

      Subseg* subseg;
      set<Subseg*>::const_iterator it = given_subsegs->begin(), it_end =
	  given_subsegs->end();

      for (; it != it_end; ++it) {

	subseg = *it;

#ifndef NDEBUG
	assert(subseg->get_beg_vert()->isValid());
	assert(subseg->get_end_vert()->isValid());
	assert(subseg->get_beg_vert() != subseg->get_end_vert());
	assert(!subseg->get_beg_vert()->isDeleted());
	assert(!subseg->get_end_vert()->isDeleted());
	m_mesh->hasVert(subseg->get_beg_vert());
	m_mesh->hasVert(subseg->get_end_vert());
#endif

	m_mesh->addSubseg(subseg);

      }

    }

    else {

      //Subsegments are not given. We must infer them from
      //the mesh and the underlying geometry.

      //This case should only happen with GRUMMP's native
      //geometry interface (and the following code is
      //designed as such). The code does not store geometric
      //entity pointers for vertices (makes it impossible to
      //identify orphaned subsegments inside a boundary surface).

      assert(!given_subsegs);

      //To be a subsegment, an edge must have both its vertices
      //either on a curve, or at an apex point.

      //Go through all the mesh boundary faces and identify all
      //the edges satisfying this condition. Store them in a set
      //as TetRefinerEdge's, this way each edge is stored
      //once and only once.
      BFace* bface;
      Vert* vert[2];

      TetRefinerEdgeBFace edge;
      set<TetRefinerEdgeBFace> all_edges;
      set<TetRefinerEdgeBFace>::iterator it, it_end;

      GR_index_t i, j, num_bfaces = m_mesh->getNumTotalBdryFaces();
      bool qBdryPatchesInferred = m_mesh->isBdryFromMesh();

#ifndef NDEBUG
      int num_bdry_verts = 0;
      for (i = 0; i < m_mesh->getNumVerts(); ++i) {
	Vert* my_vert = m_mesh->getVert(i);
	if (my_vert->isDeleted())
	  continue;
	if (my_vert->getVertType() == Vert::eBdryApex
	    || my_vert->getVertType() == Vert::eBdryCurve)
	  ++num_bdry_verts;
      }
#endif

      for (i = 0; i < num_bfaces; ++i) {

	bface = m_mesh->getBFace(i);
	if (bface->isDeleted())
	  continue;

	assert(
	    bface->getType() == Cell::eTriBFace
		|| bface->getType() == Cell::eIntTriBFace);
	assert(bface->getNumVerts() == 3);

	for (j = 0; j < 3; ++j) {

	  switch (j)
	    {
	    case 0:
	      vert[0] = bface->getVert(0);
	      vert[1] = bface->getVert(1);
	      break;
	    case 1:
	      vert[0] = bface->getVert(1);
	      vert[1] = bface->getVert(2);
	      break;
	    case 2:
	      vert[0] = bface->getVert(2);
	      vert[1] = bface->getVert(0);
	      break;
	    default:
	      assert(0);
	      break;
	    }

	  if ((vert[0]->getVertType() == Vert::eBdryApex
	      || vert[0]->getVertType() == Vert::eBdryCurve)
	      && (vert[1]->getVertType() == Vert::eBdryApex
		  || vert[1]->getVertType() == Vert::eBdryCurve)) {

	    edge = TetRefinerEdgeBFace(vert[0], vert[1]);
	    it = all_edges.find(edge);

	    if (it == all_edges.end()) { //edge not in the set.
	      edge.add_bface(bface);
	      all_edges.insert(edge);
	    }
	    else { //edge already in set.
	      edge = *it;
	      assert(!edge.get_bfaces().empty());
	      edge.add_bface(bface);
	      all_edges.erase(it);
	      all_edges.insert(edge);
	    }

	  }

	  else
	    continue; //definitely not a subseg, move on!

	}

      }

      //Now we have a set of unique edges that are potentially
      //subsegments. Go through this set and identify which edges
      //are truly subsegments.

      Subseg* subseg;
      it = all_edges.begin();
      it_end = all_edges.end();

      for (; it != it_end; ++it) {

	subseg = NULL;
	vert[0] = it->vert(0);
	vert[1] = it->vert(1);

	//This piece of code should only be used when using GRUMMP
	//without an external geometry library (i.e. CGM). As such,
	//there should not be geometry entity stored with the vertices.
	assert(!vert[0]->getParentTopology());
	assert(!vert[1]->getParentTopology());

	switch (it->get_bfaces().size())
	  {

	  case 0: //These two cases should never happen with GRUMMP's native geometry,
	  case 1: //otherwise something is wrong.
	    assert(0);
	    break;

	  case 2:
	    {
	      BFace* bface0 = *((it->get_bfaces().begin()));
	      BFace* bface1 = *(++(it->get_bfaces().begin()));
	      assert(bface0 && bface1);
	      assert(bface0 != bface1);
	      assert(!bface0->isDeleted() && !bface1->isDeleted());
	      assert(bface0->getPatchPointer());
	      assert(bface1->getPatchPointer());
	      double dDot = 0.;
	      double norm0[3], norm1[3];
	      bface0->calcUnitNormal(norm0);
	      bface1->calcUnitNormal(norm1);
	      dDot = dDOT3D(norm0, norm1);

	      //				logMessage(2,"qq\n");
	      //				bface0->getFace(0)->printFaceInfo();
	      //				bface1->getFace(0)->printFaceInfo();
	      if ((bface0->getPatchPointer() != bface1->getPatchPointer()
		  && !qBdryPatchesInferred)
		  || (fabs(dDot) < 0.99999999 && qBdryPatchesInferred))
		subseg = new Subseg(vert[0], vert[1]);
	      break;
	    }

	  default: //If there are more than 2 distinct bfaces, we definitely have a subseg.
	    subseg = new Subseg(vert[0], vert[1]);
	    break;

	  }

	if (subseg)
	  m_mesh->addSubseg(subseg);

      }

      //Must make sure that every vertex on a curve or at an apex
      //has an entry in the subseg map.
      // not necessarily true;
      if (qBdryPatchesInferred) {
	assert(num_bdry_verts == m_mesh->subsegMapSize());
      }
    } //End case where subsegments are not given.

    //Assert that the entries in the subseg map are consistent.
    assert(m_mesh->subsegMapIsValid());
    m_mesh->printSubsegMap();
  }

  bool
  RefinementManager3D::clearTopSubsegQueue()
  {

    if (!IsInParallel) // Serial Section
    {
      //Clears deleted subsegs from the top of the queue.
      //Returns true if there is still something in the queue.
      //Otherwise, returns false.

      while (!m_subseg_queue.empty()) {
	if (m_subseg_queue.front()->is_deleted()
	    || m_lockedSubsegs.count(m_subseg_queue.front()) == 1)
	  m_subseg_queue.pop_front();
	else
	  return true;

      }

      //	m_subseg_map.release_deleted_subsegs();
      return false;
    }
    else // Parallel Section
    {
//		assert(omp_in_parallel() != 0);
      int ID = omp_get_thread_num();
      while (!Para_m_subseg_queue[ID].empty()) {
	if (Para_m_subseg_queue[ID].front()->is_deleted()
	    || m_lockedSubsegs.count(m_subseg_queue.front()) == 1)
	  Para_m_subseg_queue[ID].pop_front();
	else
	  return true;
      }

      //	m_subseg_map.release_deleted_subsegs();
      return false;
    }
  }

  bool
  RefinementManager3D::clearTopBdryFaceQueue()
  {

    if (!IsInParallel) // Serial
    {
      //Clears deleted subfaces from the top of the queue.
      //Also removes subfaces that are considered too small with
      //respect to the length scale.

      BFace* bface;
      InsertionQueueEntry IQE;
      while (!m_bdryFacesToSplit.empty()) {

	IQE = m_bdryFacesToSplit.front();
	if (!IQE.qStillInMesh() || m_lockedBdryFaces.count(IQE) == 1
	    || !isBdryFaceEncroached(IQE.pBFBFace())) {
	  m_bdryFacesToSplit.pop_front();
	}
	else {
	  //To prevent infinitely recursive insertion, don't split
	  //BFace's that are too small relative to the local length
	  //scale (at this point, 15 times smaller is counted as "too
	  //small"). In this situation, the face is removed from the queue.
	  bface = IQE.pBFBFace();
	  double circ_cent[3];
	  bface->calcCircumcenter(circ_cent);
	  double circ_rad = dDIST3D(bface->getVert(0)->getCoords(), circ_cent),
	      size_now = (2. / sqrt(2.)) * circ_rad, size_limit = 0.;
	  for (int ii = 0; ii < 3; ++ii) {
	    size_limit += bface->getVert(ii)->getLengthScale();

	  }

	  if (size_now < size_limit / 45.) { //45 == 15 * 3 ; 3 = num verts; 15 = size factor
	    m_bdryFacesToSplit.pop_front();
	  }
	  else
	    //if we get here, then the subface can safely be split.
	    return true;
	}

      }

      return false;
    }
    else // Parallel
    {
//		assert(omp_in_parallel() != 0);
      int ID = omp_get_thread_num();
      //Clears deleted subfaces from the top of the queue.
      //Also removes subfaces that are considered too small with
      //respect to the length scale.

      BFace* bface;
      InsertionQueueEntry IQE;
      while (!Para_m_bdryFacesToSplit[ID].empty()) {

	IQE = Para_m_bdryFacesToSplit[ID].front();
	if (!IQE.qStillInMesh() || m_lockedBdryFaces.count(IQE) == 1) {
	  Para_m_bdryFacesToSplit[ID].pop_front();
	}
	else {
	  //To prevent infinitely recursive insertion, don't split
	  //BFace's that are too small relative to the local length
	  //scale (at this point, 15 times smaller is counted as "too
	  //small"). In this situation, the face is removed from the queue.
	  bface = IQE.pBFBFace();
	  double circ_cent[3];
	  bface->calcCircumcenter(circ_cent);
	  double circ_rad = dDIST3D(bface->getVert(0)->getCoords(), circ_cent),
	      size_now = (2. / sqrt(2.)) * circ_rad, size_limit = 0.;
	  for (int ii = 0; ii < 3; ++ii) {
	    size_limit += bface->getVert(ii)->getLengthScale();

	  }

	  if (size_now < size_limit / 45.) { //45 == 15 * 3 ; 3 = num verts; 15 = size factor
	    Para_m_bdryFacesToSplit[ID].pop_front();
	  }
	  else
	    //if we get here, then the subface can safely be split.
	    return true;
	}

      }

      return false;

    }
  }

  int
  RefinementManager3D::refineSubsegsForLength()
  {

    assert(m_isoLength);
    assert(m_subseg_queue.empty());

    Subseg* subseg;
    double length, vert0_ls, vert1_ls;
    std::pair<Subseg*, Subseg*> new_subsegs;

    set<Subseg*> all_subsegs;
    m_mesh->getAllSubsegs(all_subsegs);

    //Start by stuffing all the subsegments into m_subseg_queue.
    //Will pop those deleted and those short enough in the loop.
    std::copy(all_subsegs.begin(), all_subsegs.end(),
	      std::back_inserter(m_subseg_queue));

    int num_inserts = 0;

    //The refinement loop. Split subsegments until they are short enough.

    do {

      if (clearTopSubsegQueue()) {

	subseg = m_subseg_queue.front();
	assert(!subseg->is_deleted());

	length = subseg->get_length();
	vert0_ls = subseg->get_beg_vert()->getLengthScale();
	vert1_ls = subseg->get_end_vert()->getLengthScale();

	if (length > vert0_ls + vert1_ls) {

	  if (splitSubseg(subseg, &new_subsegs)) {
	    m_subseg_queue.push_back(new_subsegs.first);
	    m_subseg_queue.push_back(new_subsegs.second);
	    ++num_inserts;
	  }
	  else
	    m_subseg_queue.pop_front();

	}
	else
	  m_subseg_queue.pop_front();

      }

    }
    while (!m_subseg_queue.empty());

    return num_inserts;

  }

  int
  RefinementManager3D::refineBoundaryForLength()
  {

    m_refineBdryForSize = true;
    // Start by refining subsegs for length, then move onto boundary faces..
    // easier to have this as one operation.
    logMessage(1, "Refining subsegments to match the sizing field... ");
    int num_inserts = refineSubsegsForLength();
    logMessage(1, "done - inserted %d vertices.\n", num_inserts);
    logMessage(1, "Refining subfaces to match the sizing field... ");
    // now refine the rest of the boundary faces, refining subsegs based on
    // encroachment
    assert(m_isoLength);
    //These should be empty in the grand scheme of things,
    //although the function should still work if they are not...
    assert(m_subseg_queue.empty());
    while (!m_bdryFacesToSplit.empty()) {
      m_bdryFacesToSplit.pop_front();
    }
    assert(m_bdryFacesToSplit.empty());

    Subseg* subseg;
    num_inserts = 0;

    //Start by fixing subsegment's encroachement.
    //Note that subsegs getting encroached by the splits
    //will be queued automagically in the split_subseg function.

    queueEncroachedSubsegs();

    while (!m_subseg_queue.empty()) {

      subseg = m_subseg_queue.front();
      m_subseg_queue.pop_front();

      if (subseg->is_deleted())
	continue;

      num_inserts += splitSubseg(subseg);

    }
    //Empty the subface queue. It currently contains subfaces that
    //were found to be encroached while fixing subsegment's encroachment.
    //We do not want to fix subface's encroachment here, just split the
    //subfaces that are too big.
    while (!m_bdryFacesToSplit.empty()) {
      m_bdryFacesToSplit.pop_front();
    }
    assert(m_bdryFacesToSplit.empty());
    assert(m_subseg_queue.empty()); //we just fixed encroachment.
    //Refine subfaces for length.
    BFace* bface;
    bool qSplit;

    double center[3], radius, current_size, target_size, sqrt_of_two = sqrt(2.),
	one_third = 1. / 3.;

    //Put all the valid subfaces in the queue.
    for (GR_index_t ii = 0; ii < m_mesh->getNumTotalBdryFaces(); ++ii) {
      bface = m_mesh->getBFace(ii);
      if (bface->isDeleted())
	continue;

      m_bdryFacesToSplit.push_back(InsertionQueueEntry(bface));
    }
    //Refinement loop.
    while (true) {
      //		logMessage(1,"Refinement Queue: %lu subsegs %lu bdryfaces %d cells\n",
      //						m_subseg_queue.size(), m_bdryFacesToSplit.size(),
      //						m_cell_queue->queue_size());
      //Fix subsegment encroachment.
      if (!m_subseg_queue.empty()) {
	subseg = m_subseg_queue.front();
	m_subseg_queue.pop_front();
	if (subseg->is_deleted())
	  continue;
	num_inserts += splitSubseg(subseg);

      }

      else if (!m_bdryFacesToSplit.empty()) {
	InsertionQueueEntry IQE = m_bdryFacesToSplit.front();
	if (!IQE.qStillInMesh()) {
	  m_bdryFacesToSplit.pop_front();
	  continue;
	}
	bface = IQE.pBFBFace();

	//Start by computing the current size and the target size.
	assert(bface->getNumVerts() == 3);

	bface->calcCircumcenter(center);
	radius = dDIST3D(center, bface->getVert(0)->getCoords());
	current_size = sqrt_of_two * radius;
	target_size = one_third
	    * (bface->getVert(0)->getLengthScale()
		+ bface->getVert(1)->getLengthScale()
		+ bface->getVert(2)->getLengthScale());
	//If the subface is too big, then split it and queue the newly created subfaces.
	if (current_size > target_size) {
	  qSplit = splitBdryFace(bface);
	  if (qSplit) {
	    num_inserts++;
	    //A split actually occured (no encroachment on subsegs).
	    //Pop the subface and queue the newly created subfaces for further size check.
	  }
	  m_bdryFacesToSplit.pop_front();

	}

	//The subface is small enough, pop it and forget about it for now.
	else {
	  m_bdryFacesToSplit.pop_front();

	}
      }

      else
	break;

    }

    //	for(GR_index_t ii = 0; ii < m_mesh->getNumTotalBdryFaces(); ++ii) {
    //		bface = m_mesh->getBFace(ii);
    //		if(bface->isDeleted()) continue;
    //		bface->calcCircumcenter(center);
    //				radius = dDIST3D(center, bface->getVert(0)->getCoords());
    //				current_size = sqrt_of_two * radius;
    //				target_size  = one_third * (bface->getVert(0)->getLengthScale() +
    //						bface->getVert(1)->getLengthScale() + bface->getVert(2)->getLengthScale());
    //				printf("curr %f targ %f\n",current_size,target_size);
    //	}
    //The queues should definitely be empty at this point.
    //Note that the some boundary faces can still be encroached,
    //as we made no attempt to fix their encroachment.
    assert(m_subseg_queue.empty());
    assert(m_bdryFacesToSplit.empty());
    m_refineBdryForSize = false;
    logMessage(1, "done - inserted %d vertices.\n", num_inserts);

    return num_inserts;

  }

  int
  RefinementManager3D::runRefinementForEncroachedSubsegs()
  {
    int num_inserts = 0;

    Subseg* subseg;

    while (true) {
      if (clearTopSubsegQueue()) {
	subseg = m_subseg_queue.front();
	assert(!subseg->is_deleted());

	num_inserts += splitSubseg(subseg);

	assert(subseg->is_deleted());

      }
      else
	break;
    }
    return num_inserts;
  }
  int
  RefinementManager3D::runRefinementForEncroachedBdryFaces()
  {
    int num_inserts = 0;

    bool qSplit;

    InsertionQueueEntry BFaceIQE;
    InsertionQueueEntry lastBFaceIQE;
//        CellQueueEntry* lastCellEntry = NULL;

    while (true) {
      if (clearTopBdryFaceQueue()) {

	BFaceIQE = m_bdryFacesToSplit.front();
	assert(BFaceIQE.qStillInMesh());

	if (BFaceIQE == lastBFaceIQE && lastBFaceIQE.qStillInMesh()) {
	  m_bdryFacesToSplit.pop_front();
	  continue;
	}

	qSplit = splitBdryFace(BFaceIQE.pBFBFace());

	if (qSplit) {
	  ++num_inserts;
//                    lastCellEntry = NULL;
	}
	else {
	  lastBFaceIQE = BFaceIQE;
	}

      }
      else
	break;
    }
    return num_inserts;
  }
  int
  RefinementManager3D::runRefinementLoopForEachLeaf(int ID, int FLevel,
						    GR_index_t& RVerts)
  {
    bool initial_clean = false;
    int num_inserts = 0;

    bool qSplit;
    SplitType split_type;
    Subseg* subseg;

    InsertionQueueEntry BFaceIQE;
    InsertionQueueEntry lastBFaceIQE;
    CellQueueEntry* lastCellEntry = NULL;
    CellQueueEntry* cellEntry;
//        queueAllCellsThreadWise(ID, Filter);

    RVerts = 0;        //Para_m_cell_queue[ID]->queue_size();

//        if (!OKTOGO) logMessage(
//                3,
//                "Refinement Queue in leaf number %d : %lu subsegs %lu bdryfaces %d cells at level %d\n",
//                ThreadLeaves[ID]->Identifier, Para_m_subseg_queue[ID].size(),
//                Para_m_bdryFacesToSplit[ID].size(),
//                Para_m_cell_queue[ID]->queue_size(), FLevel);
    queueAllCellsThreadWise(ID, FLevel);

    while (true) {

      RVerts += Para_m_cell_queue[ID]->queue_size();
      if (clearTopSubsegQueue()) {
	subseg = Para_m_subseg_queue[ID].front();
	assert(!subseg->is_deleted());

	num_inserts += splitSubsegInParallel(subseg);
	assert(subseg->is_deleted());
	lastCellEntry = NULL;

      }

      else if (clearTopBdryFaceQueue()) {

	BFaceIQE = Para_m_bdryFacesToSplit[ID].front();
	assert(BFaceIQE.qStillInMesh());

	if (BFaceIQE == lastBFaceIQE && lastBFaceIQE.qStillInMesh()) {
	  Para_m_bdryFacesToSplit[ID].pop_front();
	  continue;
	}
	//			BFaceIQE.pBFBFace()->printCellInfo();

	qSplit = splitBdryFaceInParallel(BFaceIQE.pBFBFace());

	if (qSplit) {
	  ++num_inserts;
	  lastCellEntry = NULL;
	}
	else {
	  lastBFaceIQE = BFaceIQE;
	}

      }

      else if (!Para_m_cell_queue[ID]->empty()) {
	if (!initial_clean) {
	  initial_clean = true;
	  Para_m_cell_queue[ID]->clean_queue();
	}
	if (Para_m_cell_queue[ID]->empty())
	  break;

	cellEntry = Para_m_cell_queue[ID]->top();
	if (cellEntry->isDeleted() || cellEntry->is_equal(lastCellEntry)
	    || !cellEntry->get_cell()->doFullCheck()) {
	  Para_m_cell_queue[ID]->pop();
	  delete cellEntry;
	  lastCellEntry = NULL;
	  continue;
	}
	lastCellEntry = cellEntry;
	split_type = splitTetInParallel(cellEntry->get_cell());

	if (split_type == CELL) {
	  lastCellEntry = NULL;
	  ++num_inserts;
	}
	else if (split_type == LOCKED) {
	  Para_m_bdryFacesToSplit[ID].clear();
	  Para_m_subseg_queue[ID].clear();
	  Para_m_cell_queue[ID]->pop();
	  delete cellEntry;
	  lastCellEntry = NULL;
	  continue;
	}
	else {
	  assert(split_type == NONE);
	  lastCellEntry = cellEntry;
	}
      }
      else
	break; //all queues are empty, stop here!
    }

    logMessage(3, "inserts %d verticies at leaf %d\n", num_inserts,
	       ThreadLeaves[ID]->Identifier);
    assert(Para_m_cell_queue[ID]->empty());
    assert(Para_m_bdryFacesToSplit[ID].empty());
    assert(Para_m_subseg_queue[ID].empty());

    return num_inserts;
  }
  int
  RefinementManager3D::runRefinementLoop()
  {

    if (!IsInParallel) {
      bool initial_clean = false;
      int num_inserts = 0;

      bool qSplit;
      SplitType split_type;
      Subseg* subseg;

      InsertionQueueEntry BFaceIQE;
      InsertionQueueEntry lastBFaceIQE;
      CellQueueEntry* lastCellEntry = NULL;
      CellQueueEntry* cellEntry;

      while (true) {
	logMessage(MSG_DEBUG,
		   "Refinement Queue: %lu subsegs %lu bdryfaces %d cells\n",
		   m_subseg_queue.size(), m_bdryFacesToSplit.size(),
		   m_cell_queue->queue_size());
	if (m_mesh->getNumVerts() % 5000 == 0) {
	  logMessage(MSG_SERVICE,
		     "Refinement Queue: %lu subsegs %lu bdryfaces %d cells\n",
		     m_subseg_queue.size(), m_bdryFacesToSplit.size(),
		     m_cell_queue->queue_size());
	}

	if (clearTopSubsegQueue()) {
	  logMessage(MSG_DEBUG, "Splitting subseg\n");

	  subseg = m_subseg_queue.front();
	  assert(!subseg->is_deleted());

	  num_inserts += splitSubseg(subseg);
	  assert(subseg->is_deleted());
	  lastCellEntry = NULL;

	}

	else if (clearTopBdryFaceQueue()) {
	  logMessage(MSG_DEBUG, "Splitting bdry face\n");
	  BFaceIQE = m_bdryFacesToSplit.front();
	  assert(BFaceIQE.qStillInMesh());

	  if (BFaceIQE == lastBFaceIQE && lastBFaceIQE.qStillInMesh()) {
	    m_bdryFacesToSplit.pop_front();
	    continue;
	  }
	  //			BFaceIQE.pBFBFace()->printCellInfo();

	  qSplit = splitBdryFace(BFaceIQE.pBFBFace());

	  if (qSplit) {
	    ++num_inserts;
	    lastCellEntry = NULL;
	  }
	  else {
	    lastBFaceIQE = BFaceIQE;
	  }

	}

	else if (!m_cell_queue->empty()) {
	  if (!initial_clean) {
	    initial_clean = true;
	    m_cell_queue->clean_queue();
	  }

	  if (m_cell_queue->empty())
	    break;

	  cellEntry = m_cell_queue->top();
	  while (cellEntry->isDeleted() || cellEntry->is_equal(lastCellEntry)) {
	    m_cell_queue->pop();
	    delete cellEntry;
	    lastCellEntry = NULL;
	    if (m_cell_queue->empty()) {
	      cellEntry = NULL;
	      break;
	    }
	    cellEntry = m_cell_queue->top();
	    //continue;
	  }
	  if (cellEntry) {
	    lastCellEntry = cellEntry;
	    split_type = splitTet(cellEntry->get_cell());

	    if (split_type == CELL) {
	      lastCellEntry = NULL;
	      ++num_inserts;
	    }
	    else if (split_type == LOCKED) {
	      m_bdryFacesToSplit.clear();
	      m_subseg_queue.clear();
	      m_cell_queue->pop();
	      delete cellEntry;
	      lastCellEntry = NULL;
	      continue;
	    }
	    else {
	      assert(split_type == NONE);
	      lastCellEntry = cellEntry;
	    }
	  }
	}
	else
	  break; //all queues are empty, stop here!
      }

      return num_inserts;
    }
    else // Parallel section
    {
      int num_inserts = 0;
      FilterLevel = 3;
      Granularity = 1;
      Filter = 0;
      GR_index_t RemainingVerts = 0;
      Filter += Granularity;
      double SETime = 0;
      m_mesh->IsInParallel(true);
      m_mesh->Mesh3D_sync_entry(SETime, false);
      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared (num_inserts,RemainingVerts)
      {
#pragma omp single
	{
	  int num_insertsOld = 0;
	  int num_inserts_so_far = 0;
	  int passes = 1;
	  bool FlagNewLayer = true;
	  do {
	    double TimeThisPass = omp_get_wtime();
	    RemainingVerts = 0;
	    if (!OKTOGO)
	      logMessage(
		  3,
		  "pass #%d: Refining %lu leaves with filter level of %d ... ",
		  passes, OctT->OctLeafList.size(), Filter);
	    while (!OctT->IsWBLEmpty() || !IndependentSet.empty()) {
	      double PaTime = omp_get_wtime();
	      while (!IndependentSet.empty()) {
		std::set<Rect3D*>::reverse_iterator ritI =
		    IndependentSet.rbegin();
		Rect3D * LeafI = *ritI;
		IndependentSet.erase(*ritI);
#pragma omp task firstprivate(LeafI) shared(RemainingVerts)
		{
		  double thread_time = omp_get_wtime();
		  int num_inserts_para = 0;
		  int ID = omp_get_thread_num();
		  ThreadLeaves[ID] = LeafI;
		  GR_index_t RVerts = 0;

		  num_inserts_para = runRefinementLoopForEachLeaf(ID, Filter,
								  RVerts);

#pragma omp atomic
		  num_inserts += num_inserts_para;
#pragma omp atomic
		  RemainingVerts += RVerts;
		  thread_time = omp_get_wtime() - thread_time;
		  ThreadTime[ID][100] += thread_time;
		  ThreadInNum[ID][100] += num_inserts_para;
		} // End of Task
	      } // End of IndependentSet
#pragma omp taskwait
	      PaTime = omp_get_wtime() - PaTime;
	      JustParallelTime += (PaTime);

	      double FISTime = omp_get_wtime();
	      OctT->FindIndSets(IndependentSet, FlagNewLayer);
	      FISTime = omp_get_wtime() - FISTime;
	      FIndSetsTime += FISTime;
	    } // End of one Refining pass

	    OKTOGO = false;
	    double SendingTime = omp_get_wtime();
	    ThreadID = omp_get_thread_num();
	    m_mesh->sendEvents();
	    SendingTime = omp_get_wtime() - SendingTime;
	    SendEventsTime += SendingTime;
	    if (OKTOGO && RemainingVerts > 0) {
	      IndependentSet.clear();
	      OctT->FillWBZ();
	      continue;
	    }
	    if (Filter < FilterLevel) {
	      Filter += Granularity;
	      Granularity += 0;
	    }
	    else if (Filter <= 2 * FilterLevel) {
	      Granularity += 1;
	      Filter += Granularity;
	    }
	    else {
	      Granularity += 0;
	      Filter += Granularity;
	    }

	    TimeThisPass = omp_get_wtime() - TimeThisPass;
	    num_inserts_so_far = num_inserts - num_insertsOld;
	    logMessage(3, "%d verticies inserted in %4.10f seconds,\n",
		       num_inserts_so_far, TimeThisPass);
	    num_insertsOld = num_inserts;

	    double QCTime = omp_get_wtime();
	    int OldQTSize = OctT->OctLeafList.size();
	    passes++;
	    OctT->OctTreeCreatorSecondPasses(passes);
	    int NewQTSize = OctT->OctLeafList.size();
	    QCTime = omp_get_wtime() - QCTime;
	    SecondPassesTime += QCTime;

	    if (NewQTSize == OldQTSize)
	      FlagNewLayer = false;
	    else {
	      FlagNewLayer = true;
	      OctT->clearKTInLeaves();
	    }
	    IndependentSet.clear();
	    OctT->FillWBZ(true);
	  }
	  while (RemainingVerts > 0 && num_inserts_so_far > 0.001 * num_inserts);
	} // End of Single
      } // End of parallel
      m_mesh->IsInParallel(false);
      m_mesh->Mesh3D_sync_Exit(false);
      return num_inserts;
    }
    assert(0);
  }
  void
  RefinementManager3D::refineMesh(bool Para_Flag, int NOL)
  {
    IsInParallel = Para_Flag;
    if (!IsInParallel) {

      //The queues should be empty as we have
      //not done anything yet to the mesh.

      assert(m_cell_queue);
      assert(m_subseg_queue.empty());
      assert(m_bdryFacesToSplit.empty());

      queueAllCells();

      CellQueueEntry* entry;

      int num_inserts = 0;

      logMessage(1, "Running interior pre-refinement... ");

//		m_mesh->setWriteTempMesh(true);
      while (!m_cell_queue->empty()) {

	entry = m_cell_queue->top();
	if (entry->isDeleted()) {
	  m_cell_queue->pop();
	  delete entry;
	  continue;
	}
	SplitType split_type = splitTet(entry->get_cell());
	if (split_type == CELL) {
	  ++num_inserts;
	  if (num_inserts % 500 == 0) {
	    logMessage(MSG_SERVICE, "inserted %d vertices.\n", num_inserts);
//					exit(1);
	  }
	}
	else {
	  m_cell_queue->pop();
	  delete entry;

	}

      }
      logMessage(1, "inserted %d vertices.\n", num_inserts);
      //Next, we fix encroachment and refine the interior for quality
      //and size. Note however than until the boundaries have not been
      //refined for size, the interior will not achieve the desired
      //size targets because we average a tet's four vertices to determine
      //its size.
      logMessage(
	  1,
	  "Fixing initial encroachment and running first refinement stage... ");
      //
      while (!m_bdryFacesToSplit.empty()) {
	m_bdryFacesToSplit.pop_front();
      }
      m_subseg_queue.clear();
      assert(m_cell_queue->empty());
      assert(m_bdryFacesToSplit.empty());

      queueEncroachedSubsegs();
      queueEncroachedBdryFaces();
      queueAllCells();

      num_inserts = runRefinementLoop();

      logMessage(1, "inserted %d vertices.\n", num_inserts);

      //Now that the boundary is not encroached and the mesh is of high quality,
      //we refine the boundary until it achieves the desired length scale.
      assert(m_subseg_queue.empty());
      assert(m_bdryFacesToSplit.empty());
      assert(m_cell_queue->empty());

      num_inserts = refineBoundaryForLength();

      num_inserts = 0;

      //Finally, refine the mesh for size and quality while keeping the
      //the boundary entities unencroached.
      //	std::map<Vert* , Vert* > *vertMap = new std::map<Vert* , Vert *>();
      //	m_mesh->purgeAllEntities(vertMap);
      //	m_mesh->updateSubsegVertDataAfterPurge(*vertMap);
      //	delete vertMap;
      //	vertMap->clear();
      //
      //		m_mesh->purgeAllEntities();

      assert(m_subseg_queue.empty());
      assert(m_bdryFacesToSplit.empty());
//		setQueueForQuality();
      queueAllCells();
      //	m_cell_queue->print_queue();

      logMessage(1, "Final refinement to achieve size and quality bounds... ");
      double ToTime = omp_get_wtime();
      num_inserts = runRefinementLoop();
      ToTime = omp_get_wtime() - ToTime;
      TotalRefineTime += ToTime;
      logMessage(1, "inserted %d vertices.\n", num_inserts);
      assert(m_mesh->isValid());

      //	std::map<Vert* , Vert* > *vertMap = new std::map<Vert* , Vert *>();
      //	m_mesh->purgeAllEntities(vertMap);
      //	m_mesh->updateSubsegVertDataAfterPurge(*vertMap);
      //	vertMap->clear();
      //	delete vertMap;
      //	m_mesh->purgeAllEntities();

      //	print_quality_data();
      //	m_mesh->breakBadBdryCells(true);
    }
    else {
      refineForQuality(Para_Flag, NOL);
    }
  }

  void
  RefinementManager3D::refineForQuality(bool Para_Flag, int NOL)
  {

    IsInParallel = false;
    //The queues should be empty as we have
    //not done anything yet to the mesh.

    assert(m_cell_queue);
    assert(m_subseg_queue.empty());
    assert(m_bdryFacesToSplit.empty());

    queueAllCells(false);

    CellQueueEntry* entry;
    int num_inserts = 0;

    logMessage(1, "Running interior pre-refinement... ");

    while (!m_cell_queue->empty()) {

      entry = m_cell_queue->top();
      if (entry->isDeleted()) {
	m_cell_queue->pop();
	delete entry;
	continue;
      }
      SplitType split_type = splitTet(entry->get_cell());
      if (split_type == CELL) {
	++num_inserts;
      }
      else {
	m_cell_queue->pop();
	delete entry;

      }
      assert(m_mesh->isValid());

    }
    logMessage(1, "inserted %d vertices.\n", num_inserts);

    //Next, we fix encroachment and refine the interior for quality
    //and size. Note however than until the boundaries have not been
    //refined for size, the interior will not achieve the desired
    //size targets because we average a tet's four vertices to determine
    //its size.

    //	m_mesh->breakBadBdryCells(true);

    logMessage(
	1,
	"Fixing initial encroachment and running first refinement stage... ");
    //
    while (!m_bdryFacesToSplit.empty()) {
      m_bdryFacesToSplit.pop_front();
    }
    m_subseg_queue.clear();
//		assert(m_newBdryFaces.empty());
    assert(m_cell_queue->empty());
    assert(m_bdryFacesToSplit.empty());

    queueEncroachedSubsegs();
    queueEncroachedBdryFaces();
    while (true) {
      if (clearTopSubsegQueue())
	runRefinementForEncroachedSubsegs();
      else if (clearTopBdryFaceQueue())
	runRefinementForEncroachedBdryFaces();
      else
	break;
//			if (RFFES == 0 && RFFEB == 0) break;
    }
    IsInParallel = false;
    queueAllCells(false);
    num_inserts = runRefinementLoop();

    logMessage(1, "inserted %d vertices.\n", num_inserts);
    //Now that the boundary is not encroached and the mesh is of high quality,
    //we refine the boundary until it achieves the desired length scale.
    assert(m_subseg_queue.empty());
    assert(m_bdryFacesToSplit.empty());
    assert(m_cell_queue->empty());

    num_inserts = refineBoundaryForLength();

    num_inserts = 0;

    //Finally, refine the mesh for size and quality while keeping the
    //the boundary entities unencroached.
    m_mesh->purgeAllEntities();

    assert(m_subseg_queue.empty());
    assert(m_bdryFacesToSplit.empty());
    IsInParallel = Para_Flag;
    queueAllCells(Para_Flag, NOL);
    //	m_cell_queue->print_queue();

    logMessage(1, "Final refinement to achieve quality bound... ");
    double ToTime = omp_get_wtime();
    num_inserts = runRefinementLoop();
    ToTime = omp_get_wtime() - ToTime;
    TotalRefineTime += ToTime;
    logMessage(1, "inserted %d vertices.\n", num_inserts);

    m_mesh->purgeAllEntities();
  }

  void
  RefinementManager3D::refineAllQueuedCellsForQuality(double angle)
  {
    assert(m_cell_queue);
    assert(m_subseg_queue.empty());
    m_bdryFacesToSplit.clear();
    assert(m_bdryFacesToSplit.empty());

    //	m_cell_queue->empty_queue();
    assert(m_mesh);

    m_cell_queue->set_max_shape_priority(sin((M_PI / 180.) * angle));
    setQueueForQuality();
    //	m_DelSwap->clearQueue();
    //	m_cell_queue->empty_queue();
    //	queueAllCells();
    //	assert(m_mesh->isValid());
    //	std::vector<Cell*> newCells;
    //	CellQueueEntry* entry;
    //	int num_inserts = 0;
    //	logMessage(1,"Running interior pre-refinement on %d cells ... ",m_cell_queue->queue_size());
    //
    //	while(!m_cell_queue->empty()) {
    //
    //		entry = m_cell_queue->top();
    //		if(entry->isDeleted()) {
    //			m_cell_queue->pop();
    //			delete entry;
    //			continue;
    //		}
    //		SplitType split_type = splitTet(entry->get_cell(),false);
    //		if(split_type == CELL){
    //			++num_inserts;
    //		} else{
    //			newCells.push_back(entry->get_cell());
    //			m_cell_queue->pop();
    //			delete entry;
    //		}
    //
    //	}
    //	logMessage(1,"inserted %d vertices.\n", num_inserts);
    //	while(!m_bdryFacesToSplit.empty()){
    //		m_bdryFacesToSplit.pop_front();
    //	}
    //	m_subseg_queue.clear();
    //	m_bdryFacesToSplit.clear();
    //	m_cell_queue->empty_queue();

    queueEncroachedBdryFaces();
    queueEncroachedSubsegs();
    queueAllCells();
    //	m_lockedBdryFaces.clear();
    //	m_lockedCells.clear();
    //	m_lockedSubsegs.clear();
    //	m_newBdryFaces.clear();

    GR_index_t num_inserts = 0;

    bool qSplit, initial_clean = false;
    SplitType split_type;
    Subseg * subseg;
    BFace * bface;
    InsertionQueueEntry BFaceIQE;
    InsertionQueueEntry lastBFaceIQE;
    CellQueueEntry * entry, *last_cell = NULL;
    logMessage(1, "Refining for minimum dihedral ...\n");
    logMessage(4, "%lu subsegs %lu bfaces %d cells %lu %lu %lu \n",
	       m_subseg_queue.size(), m_bdryFacesToSplit.size(),
	       m_cell_queue->queue_size(), m_lockedSubsegs.size(),
	       m_lockedBdryFaces.size(), m_lockedCells.size());
    while (true) {

      logMessage(4, "%lu subsegs %lu bfaces %d cells %lu %lu %lu \n",
		 m_subseg_queue.size(), m_bdryFacesToSplit.size(),
		 m_cell_queue->queue_size(), m_lockedSubsegs.size(),
		 m_lockedBdryFaces.size(), m_lockedCells.size());

      if (clearTopSubsegQueue()) {
	subseg = m_subseg_queue.front();
	assert(!subseg->is_deleted());

	if (splitSubseg(subseg)) {
	  num_inserts++;
	}
      }

      else if (clearTopBdryFaceQueue()) {
	BFaceIQE = m_bdryFacesToSplit.front();
	bface = BFaceIQE.pBFBFace();
	assert(BFaceIQE.qStillInMesh());

	if (BFaceIQE == lastBFaceIQE && lastBFaceIQE.qStillInMesh()) {

	  m_bdryFacesToSplit.pop_front();
	  //				last_cell = NULL;
	  continue;
	}

	qSplit = splitBdryFace(bface);

	if (qSplit) {
	  ++num_inserts;

	}
	else {
	  //				last_cell = entry;
	  lastBFaceIQE = BFaceIQE;
	  m_bdryFacesToSplit.pop_front();
	}

      }

      else if (!m_cell_queue->empty()) {
	if (!initial_clean) {
	  initial_clean = true;
	  m_cell_queue->clean_queue();
	}

	if (m_cell_queue->empty())
	  break;
	entry = m_cell_queue->top();
	if (entry->isDeleted()
	    || (last_cell != NULL && entry->is_equal(last_cell))
	    || m_lockedCells.count(InsertionQueueEntry(entry->get_cell()))
		== 1) {
	  m_cell_queue->pop();
	  delete entry;
	  last_cell = NULL;
	  continue;
	}
	last_cell = entry;
	split_type = splitTet(entry->get_cell());

	if (split_type == CELL) {
	  last_cell = NULL;
	  ++num_inserts;
	}
	else if (split_type == LOCKED) {
	  m_bdryFacesToSplit.clear();
	  m_subseg_queue.clear();
	  m_cell_queue->pop();
	  delete entry;
	  //				last_cell = NULL;
	  continue;
	}
	else {
	  assert(split_type == NONE);
	  last_cell = entry;
	}
      }
      else
	break; //all queues are empty, stop here!
    }

  }

  void
  RefinementManager3D::printQualityData() const
  {

    Cell* cell;
    int num_cells = m_mesh->getNumCells();

    int num_real_cells = 0;
    double dummy;

    int num_diheds, total_num_diheds = 0;
    double min_dihed = LARGE_DBL, max_dihed = -LARGE_DBL, sum_dihed = 0,
	all_diheds[6];

    int num_solids, total_num_solids = 0;
    double min_solid = LARGE_DBL, max_solid = -LARGE_DBL, sum_solid = 0,
	all_solids[4];

    double min_vol_length_ratio = LARGE_DBL, sum_vol_length_ratio = 0.;
    double min_radius_ratio = LARGE_DBL, sum_radius_ratio = 0.;

    double radius_edge = LARGE_DBL, max_radius_edge = -LARGE_DBL,
	sum_radius_edge = 0.;

    for (int i = 0; i < num_cells; ++i) {

      cell = m_mesh->getCell(i);
      if (cell->isDeleted())
	continue;

      ++num_real_cells;

      { //dihedral angles
	cell->calcAllDihed(all_diheds, &num_diheds);
	for (int j = 0; j < num_diheds; ++j) {
	  ++total_num_diheds;
	  min_dihed = std::min(min_dihed, all_diheds[j]);
	  max_dihed = std::max(max_dihed, all_diheds[j]);
	  sum_dihed += all_diheds[j];
	}
      }

      { //solid angles
	cell->calcAllSolid(all_solids, &num_solids);
	for (int j = 0; j < num_solids; ++j) {
	  ++total_num_solids;
	  min_solid = std::min(min_solid, all_solids[j]);
	  max_solid = std::max(max_solid, all_solids[j]);
	  sum_solid += all_solids[j];
	}
      }

      { //volume length
	dummy = OldPriorityCalculator::volume_length(cell);
	min_vol_length_ratio = std::min(min_vol_length_ratio, dummy);
	sum_vol_length_ratio += dummy;
      }

      { //inradius:circumradius.
	dummy = OldPriorityCalculator::normalized_shape_ratio(cell);
	min_radius_ratio = std::min(min_radius_ratio, dummy);
	sum_radius_ratio += dummy;
      }

      { //circumradius-to-shortest edge ratio.
	dummy = OldPriorityCalculator::edge_radius_ratio(cell);
	radius_edge = sqrt(0.375) / dummy;
	max_radius_edge = std::max(max_radius_edge, radius_edge);
	sum_radius_edge += radius_edge;
      }

    }

    logMessage(2, "\n--- Quality statistics ---\n\n");

    logMessage(2, "Number of cells in the mesh = %d\n\n", num_real_cells);

    logMessage(2, "***Dihedral angles***\n");
    logMessage(2, "min = %lf, max = %lf, average = %lf\n", min_dihed, max_dihed,
	       sum_dihed / total_num_diheds);

    logMessage(2, "***Solid angles***\n");
    logMessage(2, "min = %lf, max = %lf, average = %lf\n", min_solid, max_solid,
	       sum_solid / total_num_solids);

    logMessage(2, "***Volume-length ration\n");
    logMessage(2, "min = %e, average = %e\n", min_vol_length_ratio,
	       sum_vol_length_ratio / num_real_cells);

    logMessage(2, "***Aspect ratio***\n");
    logMessage(2, "min = %e, average = %e\n", min_radius_ratio,
	       sum_radius_ratio / num_real_cells);

    logMessage(2, "***Circumradius-to-shortest edge***\n");
    logMessage(2, "max = %e, average = %e\n", max_radius_edge,
	       sum_radius_edge / num_real_cells);

  }

  bool
  RefinementManager3D::splitSubseg(Subseg* const subseg_to_split,
				   pair<Subseg*, Subseg*>* new_subsegs)
  {

    assert(subseg_to_split && subseg_to_split->is_not_deleted());

    // Compute insertion coordinates:

    double insert_coord[3], insert_param;
    subseg_to_split->compute_refinement_split_data(insert_coord, insert_param);

    // Compute the Watson cavity:

    //Find the seed cells (start point for cavity computation).
    vector<Cell*> seed_cells;
    vector<BFace*> seed_bfaces;
    // for now use this one,
    if (!findSeedCellsForSubseg(insert_coord, subseg_to_split, seed_cells,
				seed_bfaces)) {
      subseg_to_split->mark_deleted();
      m_mesh->removeSubseg(subseg_to_split);
      return false;
    }

    assert(!seed_bfaces.empty());

    m_WI3D->computeSubsegHull(insert_coord, seed_cells, seed_bfaces);

    std::vector<InsertionQueueEntry> encroachedBdryFaces;

    if (m_WI3D->areBFacesEncroached()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  m_WI3D->getEncroachedBFaces();

      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	encroachedBdryFaces.push_back(InsertionQueueEntry(*iter));
      }
    }
    //At this point, we know that the vertex will be inserted and the subsegment split.
    //Create the new vertex and the two new subsegments and proceed with insertion.
    Vert* new_vert;

    //Insert new_vert in m_mesh, updating connectivity table for
    //all new entities. Obtain the newly created cells and bdry faces
    //through observing
    if (!m_WI3D->isHullValid()) {
      logMessage(2, "Invalid Subseg Hull -> Using Swapping Inserter\n");
      new_vert = m_SI3D->insertPointOnSubseg(insert_coord,
					     subseg_to_split->get_beg_vert(),
					     subseg_to_split->get_end_vert());
      if (new_vert == pVInvalidVert) {
	//			logMessage(1,"SwappingInserter failed in splitting subseg\n");
	//			logMessage(1,"unable to insert %f %f %f into the subseg \n",insert_coord[0],insert_coord[1],insert_coord[2]);
	//					subseg_to_split->get_beg_vert()->printVertInfo();
	//					subseg_to_split->get_end_vert()->printVertInfo();
	//					logMessage(1,"see temphull.vtk,tempcavity.vtk,tempmesh.vtk\n");
	//			m_WI3D->writeHullVTK("temphull.vtk");
	//			m_WI3D->writeCavityVTK("tempcavity.vtk");
	//			writeVTK_LegacyWithoutPurge(*m_mesh,"tempmesh");
	//			assert(0);
	//			printf("locking subseg %p\n",subseg_to_split);
	//			subseg_to_split->get_beg_vert()->printVertInfo();
	//			subseg_to_split->get_end_vert()->printVertInfo();
	m_lockedSubsegs.insert(subseg_to_split);
	return false;
      }
    }
    else {

      new_vert = m_WI3D->insertPointInHull();
    }
    Subseg *new_subseg0 = new Subseg(), *new_subseg1 = new Subseg();
    subseg_to_split->refinement_split(insert_coord, insert_param, new_vert,
				      new_subseg0, new_subseg1);

    new_vert->setType(Vert::eBdryCurve);
    assert(new_vert->getVertType() == Vert::eBdryCurve);

    // 6 - Update the subsegment map.
    assert(new_subseg0->is_not_deleted());
    assert(new_subseg1->is_not_deleted());
    assert(subseg_to_split->is_deleted());

    m_mesh->addSubseg(new_subseg0);
    m_mesh->addSubseg(new_subseg1);
    m_mesh->removeSubseg(subseg_to_split);

    // Compute and set the new vertex length scale.

    assert(m_isoLength);
    //m_length->set_curv_vert_length(new_vert, new_subseg0, new_subseg1);
    m_isoLength->setLengthScale(new_vert);

    // Update the queues.
    //Queue the encroached subsegments at the FRONT of m_subseg_queue.

    if (m_WI3D->areSubsegsEncroached()) {
      std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator> range =
	  m_WI3D->getEncroachedSubsegs();
      for (std::set<Subseg*>::iterator iter = range.first; iter != range.second;
	  ++iter)
	if ((*iter)->is_not_deleted() && (*iter) != subseg_to_split)
	  m_subseg_queue.push_back(*iter);
    }

    for (GR_index_t iEncBdryFace = 0; iEncBdryFace < encroachedBdryFaces.size();
	iEncBdryFace++) {
      if (encroachedBdryFaces[iEncBdryFace].qStillInMesh())
	m_bdryFacesToSplit.push_back(encroachedBdryFaces[iEncBdryFace]);
    }

    m_mesh->sendEvents();

    //Do not forget to check the new subsegs and subfaces for encroachment.
    if (isSubsegEncroached(new_subseg0))
      m_subseg_queue.push_front(new_subseg0);
    if (isSubsegEncroached(new_subseg1))
      m_subseg_queue.push_front(new_subseg1);

    //If we want to return the new subsegments:
    if (new_subsegs)
      *new_subsegs = std::make_pair(new_subseg0, new_subseg1);

    m_mesh->writeTempMesh();

    return true;

  }
  bool
  RefinementManager3D::splitSubsegInParallel(
      Subseg* const subseg_to_split, pair<Subseg*, Subseg*>* new_subsegs)
  {

//	assert(omp_in_parallel() != 0);
    assert(subseg_to_split && subseg_to_split->is_not_deleted());
    int ID = omp_get_thread_num();

    // Compute insertion coordinates:

    double insert_coord[3], insert_param;
    subseg_to_split->compute_refinement_split_data(insert_coord, insert_param);

    // Compute the Watson cavity:

    //Find the seed cells (start point for cavity computation).
    vector<Cell*> seed_cells;
    vector<BFace*> seed_bfaces;
    // for now use this one,
    if (!findSeedCellsForSubseg(insert_coord, subseg_to_split, seed_cells,
				seed_bfaces)) {
      subseg_to_split->mark_deleted();
#pragma omp critical(removeSubseg)
      {
	m_mesh->removeSubseg(subseg_to_split);
      }
      return false;
    }

    assert(!seed_bfaces.empty());

    Para_m_WI3D[ID]->computeSubsegHull(insert_coord, seed_cells, seed_bfaces);

    std::vector<InsertionQueueEntry> encroachedBdryFaces;

    if (Para_m_WI3D[ID]->areBFacesEncroached()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  Para_m_WI3D[ID]->getEncroachedBFaces();

      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	encroachedBdryFaces.push_back(InsertionQueueEntry(*iter));
      }
    }
    //At this point, we know that the vertex will be inserted and the subsegment split.
    //Create the new vertex and the two new subsegments and proceed with insertion.
    Vert* new_vert;

    //Insert new_vert in m_mesh, updating connectivity table for
    //all new entities. Obtain the newly created cells and bdry faces
    //through observing
    if (!Para_m_WI3D[ID]->isHullValid()) {
      logMessage(2, "Invalid Subseg Hull -> Using Swapping Inserter\n");
//		while (!omp_test_lock(&lck[ID])){}
//		CheckLock[ID] = true;
      new_vert = Para_m_SI3D[ID]->insertPointOnSubseg(
	  insert_coord, subseg_to_split->get_beg_vert(),
	  subseg_to_split->get_end_vert());
//		CheckLock[ID] = false;
//		omp_unset_lock(&lck[ID]);
      if (new_vert == pVInvalidVert) {
	//			logMessage(1,"SwappingInserter failed in splitting subseg\n");
	//			logMessage(1,"unable to insert %f %f %f into the subseg \n",insert_coord[0],insert_coord[1],insert_coord[2]);
	//					subseg_to_split->get_beg_vert()->printVertInfo();
	//					subseg_to_split->get_end_vert()->printVertInfo();
	//					logMessage(1,"see temphull.vtk,tempcavity.vtk,tempmesh.vtk\n");
	//			m_WI3D->writeHullVTK("temphull.vtk");
	//			m_WI3D->writeCavityVTK("tempcavity.vtk");
	//			writeVTK_LegacyWithoutPurge(*m_mesh,"tempmesh");
	//			assert(0);
	//			printf("locking subseg %p\n",subseg_to_split);
	//			subseg_to_split->get_beg_vert()->printVertInfo();
	//			subseg_to_split->get_end_vert()->printVertInfo();
	Para_m_lockedSubsegs[ID].insert(subseg_to_split);
	return false;
      }
    }
    else {

      new_vert = Para_m_WI3D[ID]->insertPointInHull();

    }
    Subseg *new_subseg0 = new Subseg(), *new_subseg1 = new Subseg();

    subseg_to_split->refinement_split(insert_coord, insert_param, new_vert,
				      new_subseg0, new_subseg1);

    new_vert->setType(Vert::eBdryCurve);
    assert(new_vert->getVertType() == Vert::eBdryCurve);

    // 6 - Update the subsegment map.
    assert(new_subseg0->is_not_deleted());
    assert(new_subseg1->is_not_deleted());
    assert(subseg_to_split->is_deleted());
#pragma omp critical(removeSubseg)
    {
      m_mesh->addSubseg(new_subseg0);
      m_mesh->addSubseg(new_subseg1);
      m_mesh->removeSubseg(subseg_to_split);
      // Compute and set the new vertex length scale.

//		assert(m_length);
      //m_length->set_curv_vert_length(new_vert, new_subseg0, new_subseg1);
      m_isoLength->setLengthScale(new_vert);
    }
//	Para_m_length[ID]->set_curv_vert_length(new_vert, new_subseg0, new_subseg1);

    // Update the queues.
    //Queue the encroached subsegments at the FRONT of m_subseg_queue.

    if (Para_m_WI3D[ID]->areSubsegsEncroached()) {
      std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator> range =
	  Para_m_WI3D[ID]->getEncroachedSubsegs();
      for (std::set<Subseg*>::iterator iter = range.first; iter != range.second;
	  ++iter)
	if ((*iter)->is_not_deleted() && (*iter) != subseg_to_split)
	  Para_m_subseg_queue[ID].push_back(*iter);
    }

    for (GR_index_t iEncBdryFace = 0; iEncBdryFace < encroachedBdryFaces.size();
	iEncBdryFace++) {
      if (encroachedBdryFaces[iEncBdryFace].qStillInMesh())
	Para_m_bdryFacesToSplit[ID].push_back(
	    encroachedBdryFaces[iEncBdryFace]);
    }

    //Do not forget to check the new subsegs and subfaces for encroachment.
    if (isSubsegEncroached(new_subseg0))
      Para_m_subseg_queue[ID].push_front(new_subseg0);
    if (isSubsegEncroached(new_subseg1))
      Para_m_subseg_queue[ID].push_front(new_subseg1);

//	checkNewBdryFacesForEncroachment();

    //If we want to return the new subsegments:
    if (new_subsegs)
      *new_subsegs = std::make_pair(new_subseg0, new_subseg1);

#ifndef NDEBUG
    m_mesh->writeTempMesh();
#endif

    return true;

  }

  bool
  RefinementManager3D::splitBdryFace(BFace* bfaceToSplit, bool splitAtEdge)
  {

    //	printf("splitting bdry face %p %d\n",bfaceToSplit,m_lockedBdryFaces.count(bfaceToSplit));
    //	bfaceToSplit->printCellInfo();
    assert(bfaceToSplit->doFullCheck());
    assert(m_mesh);
    assert(
	bfaceToSplit && bfaceToSplit->isValid() && !bfaceToSplit->isDeleted());
    assert(
	bfaceToSplit->getType() == Cell::eTriBFace
	    || bfaceToSplit->getType() == Cell::eIntTriBFace);
    // Compute insertion coordinates (circumcenter of the TriBFace):
    double insert_coord[3];
    if (splitAtEdge == false)
      bfaceToSplit->calcSplitPoint(insert_coord);
    else {
      //		bfaceToSplit->calcCentroid(insert_coord);
      double length = 0, newlength;
      int iMax = 0;
      for (int i = 0; i < 3; i++) {
	newlength = dDist(bfaceToSplit->getVert(i),
			  bfaceToSplit->getVert((i + 1) % 3));
	if (newlength > length) {
	  iMax = i;
	  length = newlength;
	}

      }
      insert_coord[0] = 0.5
	  * (bfaceToSplit->getVert(iMax)->x()
	      + bfaceToSplit->getVert((iMax + 1) % 3)->x());
      insert_coord[1] = 0.5
	  * (bfaceToSplit->getVert(iMax)->y()
	      + bfaceToSplit->getVert((iMax + 1) % 3)->y());
      insert_coord[2] = 0.5
	  * (bfaceToSplit->getVert(iMax)->z()
	      + bfaceToSplit->getVert((iMax + 1) % 3)->z());
      if (bfaceToSplit->getTopologicalParent()) {
	RefFace* refFace =
	    dynamic_cast<RefFace*>(bfaceToSplit->getTopologicalParent());
	CubitVector point(insert_coord), closest_point;
	refFace->find_closest_point_trimmed(point, closest_point);
	closest_point.get_xyz(insert_coord);
      }
    }

    //Compute the cavity itself. March out from the seeds finding
    //cells containing new_vert in their circumsphere. Never cross boundary faces.
    // -> Find the cells to remove from the (constrained Delaunay) mesh.
    // -> Find the boundary faces hit during the march (but not crossed).
    // -> Find the subsegments encroached by new_vert.

    m_WI3D->compute3DHull(insert_coord, bfaceToSplit, bfaceToSplit);

    if (m_WI3D->getNumCellsToRemove() == 1) {
      Cell * temp_cell = bfaceToSplit->getFace(0)->getOppositeCell(
	  bfaceToSplit);
      if (isInTet(insert_coord, temp_cell) == -1) {
	bfaceToSplit->calcCentroid(insert_coord);
	m_WI3D->compute3DHull(insert_coord, bfaceToSplit, bfaceToSplit);

      }
    }

    //If an encroached subsegment was found, we cannot go ahead
    //with the subface split. Queue the first offending subseg and exit.
    if (m_WI3D->areSubsegsEncroached()) {
      std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator> range =
	  m_WI3D->getEncroachedSubsegs();
      for (std::set<Subseg*>::iterator iter = range.first; iter != range.second;
	  ++iter)
	if ((*iter)->is_not_deleted())
	  m_subseg_queue.push_back(*iter);
      //		m_subseg_queue.push_back( *(range.first) );
      return false;
    }
    //		printf("inscoord %f %f %f %p\n",insert_coord[0],insert_coord[1],insert_coord[2],bfaceToSplit);

    std::vector<InsertionQueueEntry> encroachedBdryFaces;
    if (m_WI3D->areBFacesEncroached()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  m_WI3D->getEncroachedBFaces();

      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	//			m_bdryFacesToSplit.push_back(InsertionQueueEntry(*iter));
	encroachedBdryFaces.push_back(InsertionQueueEntry(*iter));
      }
    }

    //At this point we know that the vertex will be inserted.
    //Create it and insert it in m_mesh data structure.
    Vert * newVert;

    //	printf("inscoord %f %f %f %p\n",insert_coord[0],insert_coord[1],insert_coord[2],bfaceToSplit);
    //	bfaceToSplit->printCellInfo();
    //	logMessage(2,"splitting bdry face\n");
    //	bfaceToSplit->getFace(0)->printFaceInfo();

    //Insert new_vert in m_mesh, updating connectivity table for
    //all new entities. Obtain the newly created cells and bdry faces.
    // Insertion.
    if (!m_WI3D->isHullValid()) {
      // Don't do this!!
      if (splitAtEdge == false) {
	return splitBdryFace(bfaceToSplit, true);
      }

      //		new_vert = m_mesh->createVert(insert_coord);
      //had issues with points lying on an edge, so this is hacky but called so rarely..
      Face * face = bfaceToSplit->getFace(0);
      TetCell * cell = dynamic_cast<TetCell*>(face->getOppositeCell(
	  bfaceToSplit));

      GR_index_t iV0, iV1;

      if (cell->isPointOnEdge(insert_coord, iV0, iV1)
	  && face->hasVert(cell->getVert(iV0))
	  && face->hasVert(cell->getVert(iV1))) {
	newVert = m_SI3D->insertPointOnEdge(insert_coord, cell->getVert(iV0),
					    cell->getVert(iV1));
      }
      else {
	newVert = m_SI3D->insertPointOnFace(insert_coord, face);
      }
      if (newVert == pVInvalidVert) {

	//			logMessage(1,"SwappingInserter failed in bface\n");
	//			logMessage(1,"unable to insert %f %f %f into face in tempface.vtk\n",insert_coord[0],insert_coord[1],insert_coord[2]);
	//			logMessage(1,"see temphull.vtk,tempcavity.vtk,tempmesh.vtk\n");
	//
	//			writeVTK_SingleFace(face,"tempface.vtk");
	//			m_WI3D->writeHullVTK("temphull.vtk");
	//			m_WI3D->writeCavityVTK("tempcavity.vtk");
	//
	//			writeVTK_LegacyWithoutPurge(*m_mesh,"tempmesh");
	//			assert(0);
	//			printf("locking bdry face %p\n",bfaceToSplit);
	//			bfaceToSplit->printCellInfo();
	m_lockedBdryFaces.insert(InsertionQueueEntry(bfaceToSplit));
	for (int iF = 0; iF < bfaceToSplit->getNumFaces(); iF++)
	  m_lockedCells.insert(
	      InsertionQueueEntry(
		  bfaceToSplit->getFace(iF)->getOppositeCell(bfaceToSplit)));
	return false;
      }

    }
    else {

      newVert = m_mesh->createVert(insert_coord);
      m_WI3D->insertPointInHull(newVert);
    }
    //
    for (GR_index_t iEncBdryFace = 0; iEncBdryFace < encroachedBdryFaces.size();
	iEncBdryFace++) {
      if (encroachedBdryFaces[iEncBdryFace].qStillInMesh())
	m_bdryFacesToSplit.push_back(encroachedBdryFaces[iEncBdryFace]);
    }

    m_mesh->sendEvents();

    // Properly set the new vert type
    if (bfaceToSplit->getType() == Cell::eTriBFace) {
      newVert->setType(Vert::eBdry);
    }
    else {
      newVert->setType(Vert::eBdryTwoSide);
    }

    assert(
	newVert->getVertType() == Vert::eBdry
	    || newVert->getVertType() == Vert::eBdryTwoSide);

    // Compute new vertex length scale.

    assert(m_isoLength);
    //m_length->set_surf_vert_length(newVert);
    m_isoLength->setLengthScale(newVert);

    return true;

  }

  bool
  RefinementManager3D::splitBdryFaceInParallel(BFace* bfaceToSplit,
					       bool splitAtCentroid)
  {

    int ID = omp_get_thread_num();
    //	printf("splitting bdry face %p %d\n",bfaceToSplit,m_lockedBdryFaces.count(bfaceToSplit));
    //	bfaceToSplit->printCellInfo();
    assert(bfaceToSplit->doFullCheck());
    assert(m_mesh);
    assert(
	bfaceToSplit && bfaceToSplit->isValid() && !bfaceToSplit->isDeleted());
    assert(
	bfaceToSplit->getType() == Cell::eTriBFace
	    || bfaceToSplit->getType() == Cell::eIntTriBFace);
    // Compute insertion coordinates (circumcenter of the TriBFace):
    double insert_coord[3];
    if (splitAtCentroid == false)
      bfaceToSplit->calcSplitPoint(insert_coord);
    else {
      //		bfaceToSplit->calcCentroid(insert_coord);
      double length = 0, newlength;
      int iMax = 0;
      for (int i = 0; i < 3; i++) {
	newlength = dDist(bfaceToSplit->getVert(i),
			  bfaceToSplit->getVert((i + 1) % 3));
	if (newlength > length) {
	  iMax = i;
	  length = newlength;
	}

      }
      insert_coord[0] = 0.5
	  * (bfaceToSplit->getVert(iMax)->x()
	      + bfaceToSplit->getVert((iMax + 1) % 3)->x());
      insert_coord[1] = 0.5
	  * (bfaceToSplit->getVert(iMax)->y()
	      + bfaceToSplit->getVert((iMax + 1) % 3)->y());
      insert_coord[2] = 0.5
	  * (bfaceToSplit->getVert(iMax)->z()
	      + bfaceToSplit->getVert((iMax + 1) % 3)->z());
      if (bfaceToSplit->getTopologicalParent()) {
	RefFace* refFace =
	    dynamic_cast<RefFace*>(bfaceToSplit->getTopologicalParent());
	CubitVector point(insert_coord), closest_point;
	refFace->find_closest_point_trimmed(point, closest_point);
	closest_point.get_xyz(insert_coord);
      }
    }

    //Compute the cavity itself. March out from the seeds finding
    //cells containing new_vert in their circumsphere. Never cross boundary faces.
    // -> Find the cells to remove from the (constrained Delaunay) mesh.
    // -> Find the boundary faces hit during the march (but not crossed).
    // -> Find the subsegments encroached by new_vert.

    Para_m_WI3D[ID]->compute3DHull(insert_coord, bfaceToSplit, bfaceToSplit);

    if (Para_m_WI3D[ID]->getNumCellsToRemove() == 1) {
      Cell * temp_cell = bfaceToSplit->getFace(0)->getOppositeCell(
	  bfaceToSplit);
      if (isInTet(insert_coord, temp_cell) == -1) {
	bfaceToSplit->calcCentroid(insert_coord);
	Para_m_WI3D[ID]->compute3DHull(insert_coord, bfaceToSplit,
				       bfaceToSplit);

      }
    }

    //If an encroached subsegment was found, we cannot go ahead
    //with the subface split. Queue the first offending subseg and exit.
    if (Para_m_WI3D[ID]->areSubsegsEncroached()) {
      std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator> range =
	  Para_m_WI3D[ID]->getEncroachedSubsegs();
      for (std::set<Subseg*>::iterator iter = range.first; iter != range.second;
	  ++iter)
	if ((*iter)->is_not_deleted())
	  Para_m_subseg_queue[ID].push_back(*iter);
      //		m_subseg_queue.push_back( *(range.first) );
      return false;
    }

    std::vector<InsertionQueueEntry> encroachedBdryFaces;
    if (Para_m_WI3D[ID]->areBFacesEncroached()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  Para_m_WI3D[ID]->getEncroachedBFaces();

      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	//			m_bdryFacesToSplit.push_back(InsertionQueueEntry(*iter));
	encroachedBdryFaces.push_back(InsertionQueueEntry(*iter));
      }

    }

    //At this point we know that the vertex will be inserted.
    //Create it and insert it in m_mesh data structure.
    Vert* new_vert;

    //
    //	logMessage(2,"splitting bdry face\n");
    //	bfaceToSplit->getFace(0)->printFaceInfo();

    //Insert new_vert in m_mesh, updating connectivity table for
    //all new entities. Obtain the newly created cells and bdry faces.
    // Insertion.

    if (!Para_m_WI3D[ID]->isHullValid()) {
      if (splitAtCentroid == false) {
	return splitBdryFaceInParallel(bfaceToSplit, true);
      }
      //		new_vert = m_mesh->createVert(insert_coord);
      //had issues with points lying on an edge, so this is hacky but called so rarely..
      Face * face = bfaceToSplit->getFace(0);
      TetCell * cell = dynamic_cast<TetCell*>(face->getOppositeCell(
	  bfaceToSplit));
      double adBary[4];

      cell->calcBarycentricCoords(insert_coord, adBary);
      int i, iLowBary = 0;
      for (i = 0; i <= 3; i++)
	if (adBary[i] < 1e-10 && adBary[i] > -1e-10) {
	  adBary[i] = 0;
	  iLowBary++;
	}
      if (iLowBary == 2) {

	Vert * pV0, *pV1;
	if (adBary[0] < 1e-10) {
	  if (adBary[1] < 1e-10) {
	    pV0 = cell->getVert(2);
	    pV1 = cell->getVert(3);
	  }
	  else if (adBary[2] < 1e-10) {
	    pV0 = cell->getVert(1);
	    pV1 = cell->getVert(3);
	  }
	  else {
	    pV0 = cell->getVert(2);
	    pV1 = cell->getVert(1);
	  }
	}
	else if (adBary[1] < 1e-10) {
	  if (adBary[2] < 1e-10) {
	    pV0 = cell->getVert(0);
	    pV1 = cell->getVert(3);
	  }
	  else {
	    pV0 = cell->getVert(2);
	    pV1 = cell->getVert(0);
	  }
	}
	else {
	  pV0 = cell->getVert(0);
	  pV1 = cell->getVert(1);
	}
	assert(pV0->isValid());
	assert(pV1->isValid());

	if (face->hasVert(pV0) && face->hasVert(pV1)) {
	  new_vert = Para_m_SI3D[ID]->insertPointOnEdge(insert_coord, pV0, pV1);
	}
	else {
	  new_vert = Para_m_SI3D[ID]->insertPointOnFace(insert_coord, face);
	}

      }
      else {
	new_vert = Para_m_SI3D[ID]->insertPointOnFace(insert_coord, face);
      }
      if (new_vert == pVInvalidVert) {

	//			logMessage(1,"SwappingInserter failed in bface\n");
	//			logMessage(1,"unable to insert %f %f %f into face in tempface.vtk\n",insert_coord[0],insert_coord[1],insert_coord[2]);
	//			logMessage(1,"see temphull.vtk,tempcavity.vtk,tempmesh.vtk\n");
	//
	//			writeVTK_SingleFace(face,"tempface.vtk");
	//			m_WI3D->writeHullVTK("temphull.vtk");
	//			m_WI3D->writeCavityVTK("tempcavity.vtk");
	//
	//			writeVTK_LegacyWithoutPurge(*m_mesh,"tempmesh");
	//			assert(0);
	//			printf("locking bdry face %p\n",bfaceToSplit);
	//			bfaceToSplit->printCellInfo();
	Para_m_lockedBdryFaces[ID].insert(InsertionQueueEntry(bfaceToSplit));
	for (int iF = 0; iF < bfaceToSplit->getNumFaces(); iF++)
	  Para_m_lockedCells[ID].insert(
	      InsertionQueueEntry(
		  bfaceToSplit->getFace(iF)->getOppositeCell(bfaceToSplit)));
	return false;
      }

    }
    else {
      new_vert = m_mesh->createVert(insert_coord);
      Para_m_WI3D[ID]->insertPointInHull(new_vert);

    }

    for (GR_index_t iEncBdryFace = 0; iEncBdryFace < encroachedBdryFaces.size();
	iEncBdryFace++) {
      if (encroachedBdryFaces[iEncBdryFace].qStillInMesh())
	Para_m_bdryFacesToSplit[ID].push_back(
	    encroachedBdryFaces[iEncBdryFace]);
    }

    // Properly set the new vert type
    if (bfaceToSplit->getNumFaces() == 1) {
      new_vert->setType(Vert::eBdry);
    }
    else {
      new_vert->setType(Vert::eBdryTwoSide);
    }

    assert(
	new_vert->getVertType() == Vert::eBdry
	    || new_vert->getVertType() == Vert::eBdryTwoSide);

//	checkNewBdryFacesForEncroachment();
#pragma omp critical(removeSubseg)
    {
      //m_length->set_surf_vert_length(new_vert);
      m_isoLength->setLengthScale(new_vert);
    }

#ifndef NDEBUG
    m_mesh->writeTempMesh();
#endif

    return true;

  }
  RefinementManager3D::SplitType
  RefinementManager3D::splitTet(Cell* const cellToSplit, bool splitAtCentroid)
  {
    //printf("splitting tet %p\n",cellToSplit);
    assert(m_mesh);
    assert(cellToSplit && cellToSplit->isValid() && !cellToSplit->isDeleted());

    //Compute the insertion location:
    double insert_coord[3];
    TetCell* tet_cell = dynamic_cast<TetCell*>(cellToSplit);
    assert(tet_cell->isValid());

    //The circumcenter of a flat tet (most probably made of 4 (close to) cocircular vertices)
    //is hard to compute. In fact, a perfectly flat tet has an infinite number of circumspheres,
    //therefore an infinite number of circumcenters. In these conditions, we insert the new
    //vertex at the centroid of the flat tet.
    if (checkOrient3D(cellToSplit->getVert(1), cellToSplit->getVert(2),
		      cellToSplit->getVert(3), cellToSplit->getVert(0)) == 0
	|| splitAtCentroid) {
      tet_cell->calcCentroid(insert_coord);
    }
    else {
      m_IPC->calcInsertionPoint(tet_cell, insert_coord);
      //		tet_cell->calcCircumcenter(insert_coord);
    }
    //Do not create a vertex in the mesh right away since we do not know
    //at this point if the vertex will actually get inserted.

    //The circumcenter of the tet must necessarily be inside the circumsphere.

    //Compute the Watson cavity:
    m_WI3D->compute3DHull(insert_coord, cellToSplit);

    //The checkOrient3D used above to detect flat tets is not sensitive enough
    //to catch all "degenerate" tets (in the sense that they cannot correctly
    //compute a correct insertion coord). The following hack takes advantage
    //of the fact that when only one cell was found to contain the new vertex
    //in its circumcenter, then the vertex must be strictly inside the cell.
    //If not, then the insertion coord must be incorrect.
    //
    if (m_WI3D->getNumCellsToRemove() == 1 && !m_WI3D->areSubsegsEncroached()) {
      if (isInTet(insert_coord, tet_cell) != 4) {
	tet_cell->calcCentroid(insert_coord);
	m_WI3D->compute3DHull(insert_coord, cellToSplit);
      }
    }

    //Must now check the boundary faces
    //hit during the walk for encroachment.
    // TODO: This is overkill, I'm confident of it. Don't know what else to do though
    if (m_WI3D->areBFacesHit()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  m_WI3D->getHitBFaces();
      BFace * bface;
      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	bface = *iter;
	if (bface->isPointEncroaching(insert_coord,
				      m_mesh->getEncroachmentType(), true)) {
	  m_bdryFacesToSplit.push_front(InsertionQueueEntry(bface));
	  //				return NONE;
	}

	//If the boundary face we hit is still encroached,
	//inserting the point will likely cause an orientation problem
	//in the mesh (should fail an assertion in correct_cavity).
	//This situation occurs when encroachment could not be fixed due
	//to a small angle and we had to bail out. Just refuse to split
	//the tet and move on.
	//This is also caught when we "pre-refine" the interior before
	//starting the real refinement.
	else if (bface->isEncroached(m_mesh->getEncroachmentType())) {
	  return LOCKED;
	}
      }
    }
    //If an encroached subsegment was found, queue it and exit.
    if (m_WI3D->areSubsegsEncroached()) {
      std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator> range =
	  m_WI3D->getEncroachedSubsegs();
      for (std::set<Subseg*>::iterator iter = range.first; iter != range.second;
	  ++iter)
	if ((*iter)->is_not_deleted())
	  m_subseg_queue.push_back(*iter);
      //		m_subseg_queue.push_back( *(range.first) );
      return NONE;
    }

    if (m_WI3D->areBFacesEncroached()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  m_WI3D->getEncroachedBFaces();

      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	m_bdryFacesToSplit.push_back(InsertionQueueEntry(*iter));
      }
      return NONE;
    }

    //If we got to this point, then the tetrahedron can safely be split.
    //First, compute the faces to remove and faces of hull from the set
    //of cells to remove.

    //Create the new vertex in the mesh data structures and
    //insert it using watson insertion.
    // clear out the queues so we know what cells are newly created
    Vert* new_vert;
    if (!m_WI3D->isHullValid()) {
      logMessage(2,
		 "Splitting Bdryface Hull Invalid, trying Swapping Inserter\n");
      if (isInTet(insert_coord, cellToSplit) > 0)
	new_vert = m_SI3D->insertPoint(insert_coord, cellToSplit);
      else
	new_vert = pVInvalidVert;
      if (new_vert == pVInvalidVert) {

	if (splitAtCentroid == false) {
	  return splitTet(cellToSplit, true);
	}
	//			logMessage(1,"SwappingInserter failed in cell\n");
	//			logMessage(1,"unable to insert %f %f %f into cell in tempcell.vtk\n",insert_coord[0],insert_coord[1],insert_coord[2]);
	//			writeVTK_SingleCell(cell,"tempcell.vtk");
	//			m_WI3D->writeHullVTK("temphull.vtk");
	//			m_WI3D->writeCavityVTK("tempcavity.vtk");
	//			writeVTK_LegacyWithoutPurge(*m_mesh,"tempmesh");
	//			logMessage(1,"see temphull.vtk,tempcavity.vtk,tempmesh.vtk\n");
	//
	//			assert(0);
	m_lockedCells.insert(InsertionQueueEntry(cellToSplit));
	return LOCKED;
      }
      if (new_vert->getVertType() == Vert::eInterior) {
	m_isoLength->setLengthScale(new_vert);
	//m_length->set_volu_vert_length(new_vert);
      }
      else {
	m_isoLength->setLengthScale(new_vert);
	//m_length->set_surf_vert_length(new_vert);
      }
    }
    else {
      new_vert = m_mesh->createVert(insert_coord);
      m_WI3D->insertPointInHull(new_vert);
      //m_length->set_volu_vert_length(new_vert);
      m_isoLength->setLengthScale(new_vert);
    }

    m_mesh->sendEvents();

    m_mesh->writeTempMesh();

    return CELL;
  }

  RefinementManager3D::SplitType
  RefinementManager3D::splitTetInParallel(Cell* const cellToSplit,
					  bool splitAtCentroid)
  {
    //printf("splitting tet %p\n",cellToSplit);
//	assert(omp_in_parallel() != 0);
    int ID = omp_get_thread_num();
    assert(m_mesh);
    assert(cellToSplit && cellToSplit->isValid() && !cellToSplit->isDeleted());

    //Compute the insertion location:
    double insert_coord[3];
    TetCell* tet_cell = dynamic_cast<TetCell*>(cellToSplit);
    assert(tet_cell->isValid());

    //The circumcenter of a flat tet (most probably made of 4 (close to) cocircular vertices)
    //is hard to compute. In fact, a perfectly flat tet has an infinite number of circumspheres,
    //therefore an infinite number of circumcenters. In these conditions, we insert the new
    //vertex at the centroid of the flat tet.
    if (checkOrient3D(cellToSplit->getVert(1), cellToSplit->getVert(2),
		      cellToSplit->getVert(3), cellToSplit->getVert(0)) == 0
	|| splitAtCentroid) {
      tet_cell->calcCentroid(insert_coord);
    }
    else {
      tet_cell->calcCircumcenter(insert_coord);
    }
    //Do not create a vertex in the mesh right away since we do not know
    //at this point if the vertex will actually get inserted.

    //The circumcenter of the tet must necessarily be inside the circumsphere.

    //Compute the Watson cavity:
    Para_m_WI3D[ID]->compute3DHull(insert_coord, cellToSplit);

    //The checkOrient3D used above to detect flat tets is not sensitive enough
    //to catch all "degenerate" tets (in the sense that they cannot correctly
    //compute a correct insertion coord). The following hack takes advantage
    //of the fact that when only one cell was found to contain the new vertex
    //in its circumcenter, then the vertex must be strictly inside the cell.
    //If not, then the insertion coord must be incorrect.
    //
    if (Para_m_WI3D[ID]->getNumCellsToRemove() == 1
	&& !Para_m_WI3D[ID]->areSubsegsEncroached()) {
      if (isInTet(insert_coord, tet_cell) != 4) {
	tet_cell->calcCentroid(insert_coord);
	Para_m_WI3D[ID]->compute3DHull(insert_coord, cellToSplit);
      }
    }

    //Must now check the boundary faces
    //hit during the walk for encroachment.
    // TODO: This is overkill, I'm confident of it. Don't know what else to do though
    if (Para_m_WI3D[ID]->areBFacesHit()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  Para_m_WI3D[ID]->getHitBFaces();
      BFace * bface;
      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	bface = *iter;
	if (bface->isPointEncroaching(insert_coord,
				      m_mesh->getEncroachmentType(), true)) {
	  Para_m_bdryFacesToSplit[ID].push_front(InsertionQueueEntry(bface));
	  //				return NONE;
	}

	//If the boundary face we hit is still encroached,
	//inserting the point will likely cause an orientation problem
	//in the mesh (should fail an assertion in correct_cavity).
	//This situation occurs when encroachment could not be fixed due
	//to a small angle and we had to bail out. Just refuse to split
	//the tet and move on.
	//This is also caught when we "pre-refine" the interior before
	//starting the real refinement.
	else if (bface->isEncroached(m_mesh->getEncroachmentType())) {
	  return LOCKED;
	}
      }
    }
    //If an encroached subsegment was found, queue it and exit.
    if (Para_m_WI3D[ID]->areSubsegsEncroached()) {
      std::pair<std::set<Subseg*>::iterator, std::set<Subseg*>::iterator> range =
	  Para_m_WI3D[ID]->getEncroachedSubsegs();
      for (std::set<Subseg*>::iterator iter = range.first; iter != range.second;
	  ++iter)
	if ((*iter)->is_not_deleted())
	  Para_m_subseg_queue[ID].push_back(*iter);
      //		m_subseg_queue.push_back( *(range.first) );
      return NONE;
    }

    if (Para_m_WI3D[ID]->areBFacesEncroached()) {
      std::pair<std::set<BFace*>::iterator, std::set<BFace*>::iterator> range =
	  Para_m_WI3D[ID]->getEncroachedBFaces();

      for (std::set<BFace*>::iterator iter = range.first; iter != range.second;
	  ++iter) {
	Para_m_bdryFacesToSplit[ID].push_back(InsertionQueueEntry(*iter));
      }
      return NONE;
    }

    //If we got to this point, then the tetrahedron can safely be split.
    //First, compute the faces to remove and faces of hull from the set
    //of cells to remove.

    //Create the new vertex in the mesh data structures and
    //insert it using watson insertion.
    // clear out the queues so we know what cells are newly created
    Vert* new_vert;
    if (!Para_m_WI3D[ID]->isHullValid()) {
      logMessage(2,
		 "Splitting Bdryface Hull Invalid, trying Swapping Inserter\n");
      if (isInTet(insert_coord, cellToSplit) > 0) {

	new_vert = Para_m_SI3D[ID]->insertPoint(insert_coord, cellToSplit);

      }
      else
	new_vert = pVInvalidVert;
      if (new_vert == pVInvalidVert) {

	if (splitAtCentroid == false) {
	  return splitTetInParallel(cellToSplit, true);
	}
	//			logMessage(1,"SwappingInserter failed in cell\n");
	//			logMessage(1,"unable to insert %f %f %f into cell in tempcell.vtk\n",insert_coord[0],insert_coord[1],insert_coord[2]);
	//			writeVTK_SingleCell(cell,"tempcell.vtk");
	//			m_WI3D->writeHullVTK("temphull.vtk");
	//			m_WI3D->writeCavityVTK("tempcavity.vtk");
	//			writeVTK_LegacyWithoutPurge(*m_mesh,"tempmesh");
	//			logMessage(1,"see temphull.vtk,tempcavity.vtk,tempmesh.vtk\n");
	//
	//			assert(0);
	Para_m_lockedCells[ID].insert(InsertionQueueEntry(cellToSplit));
	return LOCKED;
      }
#pragma omp critical(removeSubseg)
      {
	if (new_vert->getVertType() == Vert::eInterior) {
	  //m_length->set_volu_vert_length(new_vert);
	  m_isoLength->setLengthScale(new_vert);
	}
	else {
	  //m_length->set_surf_vert_length(new_vert);
	  m_isoLength->setLengthScale(new_vert);
	}
      }
    }
    else {
      new_vert = m_mesh->createVert(insert_coord);
      Para_m_WI3D[ID]->insertPointInHull(new_vert);
#pragma omp critical(removeSubseg)
      {
	//m_length->set_volu_vert_length(new_vert);
	m_isoLength->setLengthScale(new_vert);
      }
//		Para_m_length[ID]->set_volu_vert_length(new_vert);
    }

#ifndef NDEBUG
    m_mesh->writeTempMesh();
#endif

    return CELL;
  }
  void
  RefinementManager3D::queueEncroachedBdryFaces()
  {
    assert(m_mesh);
    if (m_mesh->areBdryChangesAllowed()) {
      //Go through all boundary faces and check if they are encroached.
      int num_bfaces = m_mesh->getNumTotalBdryFaces();
      for (int i = 0; i < num_bfaces; ++i) {
	BFace * bface = m_mesh->getBFace(i);
	if (bface->isDeleted())
	  continue;
	if (isBdryFaceEncroached(bface))
	  m_bdryFacesToSplit.push_back(InsertionQueueEntry(bface));

      }
      logMessage(3, "\nQueued %lu bdryfaces \n", m_bdryFacesToSplit.size());
    }
    else {
      logMessage(MSG_SERVICE, "\nBdry changes aren't allowed; encroached bdry faces not queued.\n");
    }
  }

  bool
  RefinementManager3D::isBdryFaceEncroached(const BFace* const bface) const
  {
    assert(m_mesh);
    assert(bface && !bface->isDeleted());
    assert(
	m_mesh->getEncroachmentType() == eBall
	    || m_mesh->getEncroachmentType() == eLens);

    return (bface->isEncroached(m_mesh->getEncroachmentType()) != eClean);

  }

  void
  RefinementManager3D::checkNewBdryFacesForEncroachment()
  {

    if (!IsInParallel) {
      assert(0); // Dan has erased this function!
    }
    else {
      int ID = omp_get_thread_num();
      //Do not forget to check the new subfaces for encroachment.
      for (vector<BFace*>::iterator it = Para_m_newBdryFaces[ID].begin();
	  it != Para_m_newBdryFaces[ID].end(); ++it) {
	BFace* bface = *it;
	if (bface->isDeleted())
	  continue;
	if (isBdryFaceEncroached(bface))
	  Para_m_bdryFacesToSplit[ID].push_back(InsertionQueueEntry(bface));

      }
      Para_m_newBdryFaces[ID].clear();
    }
  }
  void
  RefinementManager3D::queueEncroachedSubsegs()
  {

    assert(m_mesh);
    // Theoretically possible to have no subsegs on the subseg map
    //	assert( !m_subseg_map.empty() );

    //Go through all the subsegments and check cells attached to
    //them for a vertex that is inside the subsegment's diametral sphere.

    if (m_mesh->areBdryChangesAllowed()) {
      Subseg* subseg;
      set<Subseg*> all_subsegs;
      set<Subseg*>::iterator it, it_end;

      m_mesh->getAllSubsegs(all_subsegs);
      it = all_subsegs.begin();
      it_end = all_subsegs.end();

      for (; it != it_end; ++it) {
	subseg = *it;
	if (subseg->is_deleted())
	continue;
	if (isSubsegEncroached(subseg)) {
	  m_subseg_queue.push_back(subseg);
	}
      }
      logMessage(3, "\nQueued %lu subsegs\n", m_subseg_queue.size());
    }
    else {
      logMessage(MSG_SERVICE, "\nBdry can't be changed; no subsegs queued.\n");
    }
 
  }

  void
  RefinementManager3D::queueAllCells(bool OMP_FLAG, int NQLayers)
  {

    if (!OMP_FLAG) // serial section
    {

      assert(m_mesh);
      assert(m_cell_queue);

      Cell* cell;
      int num_cells = m_mesh->getNumCells();

      for (int i = 0; i < num_cells; ++i) {
	cell = m_mesh->getCell(i);
	if (cell->isDeleted() || cell->getType() != Cell::eTet)
	  continue;
	assert(cell->getType() == Cell::eTet);
	if (m_refineForQuality)
	  m_cell_queue->add_quality_entry(cell);
	else
	  m_cell_queue->add_entry(cell);
      }
      logMessage(3, "\nQueued %u tets\n", m_cell_queue->queue_size());

    }
    else // Parallel Section
    {
      assert(m_mesh);
      double MaxCircR = 0;
      double XCirB[2] =
	{ 0 };
      double YCirB[2] =
	{ 0 };
      double ZCirB[2] =
	{ 0 };
      double CirPos[3] =
	{ 0 };

      Cell* cell;
      int num_cells = m_mesh->getNumCells();
      for (int i = 0; i < num_cells; i++) {
	cell = m_mesh->getCell(i);
	if (cell->isDeleted())
	  continue;
	assert(cell->getType() == Cell::eTet);
	TetCell* tet = dynamic_cast<TetCell*>(cell);
	if (tet->calcCircumradius() > MaxCircR)
	  MaxCircR = tet->calcCircumradius();
	tet->calcCircumcenter(CirPos);
	if (CirPos[0] >= XCirB[1])
	  XCirB[1] = CirPos[0];
	if (CirPos[0] <= XCirB[0])
	  XCirB[0] = CirPos[0];
	if (CirPos[1] >= YCirB[1])
	  YCirB[1] = CirPos[1];
	if (CirPos[1] <= YCirB[0])
	  YCirB[0] = CirPos[1];
	if (CirPos[2] >= ZCirB[1])
	  ZCirB[1] = CirPos[2];
	if (CirPos[2] <= ZCirB[0])
	  ZCirB[0] = CirPos[2];
      }

      OctT->UpdateMaxCirR(MaxCircR);
      OctT->CalcBoundingBox();
      OctT->UpdateBBound(XCirB, YCirB, ZCirB);
      OctT->Initial_OctTree_Creator(true, NQLayers);

      for (int i = 0; i < num_cells; ++i) {
	cell = m_mesh->getCell(i);
	TetCell* tet = dynamic_cast<TetCell*>(cell);
	if (cell->isDeleted())
	  continue;
	assert(cell->getType() == Cell::eTet);
	OctT->AddTetToRect(tet);
      }

      OctT->FillWBZ(true);
      OctT->FindIndSets(IndependentSet, true);
      OctT->writeVTKOcttree();
    }
  }
  GR_index_t
  RefinementManager3D::queueAllCellsThreadWise(int ID, int Fil)
  {
    std::list<TetCell*> BL;
    GR_index_t NumAdded = 0;
    logMessage(4, "Level of refinement is: %d\n", Fil);
    if (ThreadLeaves[ID]->TetsInOct.empty()
	&& ThreadLeaves[ID]->m_newBdryFaces.empty())
      return 0;
    std::list<TetCell*>::iterator itS = ThreadLeaves[ID]->TetsInOct.begin();
    for (; itS != ThreadLeaves[ID]->TetsInOct.end(); itS++) {
      TetCell* pTC = *itS;
      if (!pTC->doFullCheck())
	continue;
      bool checkQ = OctT->QContainsCell(pTC, ThreadLeaves[ID]);
      if (checkQ) {
//			int ff = OctT->FilterTet(pTC,ThreadLeaves[ID]);
//			if (ff <= Filter)
//			{
	NumAdded++;
	if (m_refineForQuality)
	  Para_m_cell_queue[ID]->add_quality_entry(pTC);
	else
	  Para_m_cell_queue[ID]->add_entry(pTC);
//			}
//			else
//				BL.push_back(pTC);
      }
      else
	continue;
    }
    ThreadLeaves[ID]->TetsInOct.clear();
    ThreadLeaves[ID]->TetsInOct.swap(BL);

    std::vector<BFace*>::iterator ItBS =
	ThreadLeaves[ID]->m_newBdryFaces.begin(), ItBE =
	ThreadLeaves[ID]->m_newBdryFaces.end();
    for (; ItBS != ItBE; ItBS++) {
      BFace* BF = *ItBS;
      if (!BF->doFullCheck())
	continue;
      Face * pF = BF->getFace();
      TetCell* pTC = dynamic_cast<TetCell*>(pF->getOppositeCell(BF));
      if (!pTC->doFullCheck())
	continue;
      if (BF->isEncroached(m_mesh->getEncroachmentType()) != eClean)
	Para_m_bdryFacesToSplit[ID].push_back(dynamic_cast<BdryEdgeBase*>(BF));
    }
    ThreadLeaves[ID]->m_newBdryFaces.clear();
    return NumAdded;
  }

  bool
  RefinementManager3D::isSubsegEncroached(const Subseg* const subseg) const
  {

    //This is designed to work with constrained Delaunay tetrahedra only!

    assert(subseg && subseg->is_not_deleted());

    Vert *vert0 = subseg->get_beg_vert(), *vert1 = subseg->get_end_vert(),
	*this_vert;
    set<Cell*> neigh_cells;
    set<Vert*> neigh_verts, test_verts;
    findNeighborhoodInfo(vert0, neigh_cells, neigh_verts);
    assert(neigh_verts.count(vert1) == 1);

    int i;
    Cell* cell;
    set<Cell*>::iterator itc = neigh_cells.begin(), itc_end = neigh_cells.end();

    for (; itc != itc_end; ++itc) {

      cell = *itc;
      if (!cell->hasVert(vert1))
	continue;
      assert(cell->getType() == Cell::eTet && cell->getNumVerts() == 4);
      assert(cell->hasVert(vert0) && cell->hasVert(vert1));

      for (i = 0; i < 4; ++i) {
	this_vert = cell->getVert(i);
	if (this_vert == vert0 || this_vert == vert1)
	  continue;
	if (test_verts.insert(this_vert).second
	    && SubsegManager::vert_inside_subseg_ball(this_vert, subseg))
	  return true;
      }

    }

    return false;

  }

  bool
  RefinementManager3D::findSeedCellsForSubseg(double * const newLoc,
					      Subseg* const subseg,
					      std::vector<Cell*>& seed_cells,
					      std::vector<BFace*>& seed_bfaces)
  {

    Vert *vert0 = subseg->get_beg_vert(), *vert1 = subseg->get_end_vert();

    assert(vert0 && vert1);
    assert(vert0->isValid() && !vert0->isDeleted());
    assert(vert1->isValid() && !vert1->isDeleted());
    assert(vert0 != vert1);

    //Find all tets sharing vert0 and vert1
    set<Cell*> neigh_cells;
    set<Vert*> neigh_verts;
    findNeighborhoodInfo(vert0, neigh_cells, neigh_verts);

    // this is hacky, but included to prevent the case where we've swapped/eliminated a
    // subseg without knowing it, and this is good check
    // subseg is gone
    if (neigh_verts.count(vert1) == 0)
      return false;

    Cell* cell;

    set<Cell*>::iterator it = neigh_cells.begin(), it_end = neigh_cells.end();
    for (; it != it_end; ++it) {
      cell = *it;
      assert(cell->hasVert(vert0));
      if (cell->hasVert(vert1)) {
	TetCell* tet = dynamic_cast<TetCell*>(cell);
	assert(tet->isValid());
	if (tet->circumscribesPoint(newLoc))
	  seed_cells.push_back(cell);
	Cell* next_cell;
	for (int i = 0; i < 4; i++) {
	  BFace* bface;
	  next_cell = tet->getFace(i)->getOppositeCell(cell);
	  if (next_cell->isValid()) {
	    switch (next_cell->getType())
	      {

	      case Cell::eTriBFace:
	      case Cell::eIntTriBFace:

		bface = dynamic_cast<BFace*>(next_cell);
		assert(bface);
		if (bface->hasVert(vert0) && bface->hasVert(vert1))
		  seed_bfaces.push_back(bface);

		break;

	      default:
		break;
	      }
	  }
	}
      }
    }
    return true;
  }

  void
  RefinementManager3D::receiveCreatedCells(std::vector<Cell*>& createdCells)
  {
    double CCTime = omp_get_wtime();
    if (!IsInParallel) {

      logMessage(
	  MSG_DEBUG,
	  "RefinementManager3D received word about %zu newly created cells.\n",
	  createdCells.size());
      std::vector<Cell*>::iterator iter = createdCells.begin(), iterEnd =
	  createdCells.end();
      NumCreatedCells += createdCells.size();
      for (; iter != iterEnd; ++iter) {
	if ((*iter)->isValid() && !(*iter)->isDeleted()) {
	  if (m_refineForQuality)
	    m_cell_queue->add_quality_entry(*iter);
	  else
	    m_cell_queue->add_entry(*iter);
	}
      }
    }
    else {
      logMessage(
	  MSG_DEBUG,
	  "RefinementManager3D received word about %zu newly created cells.\n",
	  createdCells.size());
//            omp_lock_t *OLK = new omp_lock_t[OctT->OctLeafList.size()];

      omp_set_nested(1);
      omp_set_dynamic(0);
      omp_set_num_threads(NUM_PROCS);
#pragma omp parallel shared(createdCells)
      {
//#pragma omp for schedule(guided)
//            	for (GR_index_t i = 0 ; i<OctT->OctLeafList.size() ; i++)
//            		omp_init_lock(&(OLK[i]));
#pragma omp for schedule(guided)
	for (GR_index_t iC = 0; iC < createdCells.size(); iC++) {
	  Cell* CellT = createdCells[iC];
	  if (!CellT->doFullCheck())
	    continue;
	  TetCell* tet = dynamic_cast<TetCell*>(CellT);
	  Rect3D* LeafB = OctT->WhichRect(tet);
//                    omp_set_lock(&(OLK[LeafB->Identifier]));
//                    LeafB->UpdateCellLists(tet);
//                    omp_unset_lock(&(OLK[LeafB->Identifier]));

	  LeafB->Push_Back_Cell(tet);

	}
      }
//    		for (GR_index_t i =0 ; i<OctT->OctLeafList.size() ; i++)
//    			omp_destroy_lock(&(OLK[i]));
//    		delete OLK;
      Remaining_Cells.clear();
    }
    CCTime = omp_get_wtime() - CCTime;
    SendCCTime += CCTime;
  }

  void
  RefinementManager3D::receiveCreatedBFaces(std::vector<BFace*>& createdBFaces)
  {
    double BCTime = omp_get_wtime();
    if (!IsInParallel) {
      logMessage(
	  3,
	  "RefinementManager3D received word about %zu newly created bfaces.\n",
	  createdBFaces.size());
      std::vector<BFace*>::iterator iter = createdBFaces.begin(), iterEnd =
	  createdBFaces.end();
      for (; iter != iterEnd; ++iter) {
	BFace * bface = *iter;
	if (m_refineBdryForSize || isBdryFaceEncroached(bface))
	  m_bdryFacesToSplit.push_back(InsertionQueueEntry(bface));
      }
    }
    else {
      logMessage(
	  3,
	  "RefinementManager3D received word about %zu newly created bfaces.\n",
	  createdBFaces.size());
      std::vector<BFace*>::iterator itBRS = createdBFaces.begin(), itBRE =
	  createdBFaces.end();
      for (; itBRS != itBRE; itBRS++) {
	BFace * BFR = *itBRS;
	if (!BFR->doFullCheck())
	  continue;
	Face * pFR = BFR->getFace();
	Cell * Test = pFR->getOppositeCell(BFR);
	TetCell * pTCR = dynamic_cast<TetCell*>(Test);
	if (!pTCR->doFullCheck())
	  continue;
	Rect3D* LeafBR = OctT->WhichRect(pTCR);
	if (isBdryFaceEncroached(BFR))
	  LeafBR->m_newBdryFaces.push_back(BFR);
      }
    }
    BCTime = omp_get_wtime() - BCTime;
    SendCBTime += BCTime;
  }
//void RefinementManager3D::clearQueue(){
//	while(!m_subseg_queue.empty()){
//		m_subseg_queue.pop_front();
//	}
//	while(!m_bdryFacesToSplit.empty()){
//		m_bdryFacesToSplit.pop_front();
//	}
//	while(!m_cell_queue->empty()){
//		m_cell_queue->pop_top_entry();
//	}
//assert(m_cell_queue->empty());
//assert(m_newBdryFaces.empty());
//assert(m_subseg_queue.empty());
//}

#ifdef HAVE_MESQUITE
  GR_index_t
  AnisoEdgeRefinementManager3D::AnisoRefineMesh()
  {

    GR_index_t Added = 0;
    GR_index_t Removed = 0;

//	Added += AnisoQualityEdgeInsertion();

    GlobalSmoothing();
//	LocalSmoothing(3);

    SetReceiveCreated(false);
    m_pSM3D_2->clearQueue();
    m_pSM3D_2->swapAllFaces();
    m_pSM3D_2->clearQueue();

    Added += AnisoQualityEdgeInsertion();

    Added += AnisoMetricLengthEdgeInsertion();

    Removed += AnisoMetricEdgeCoarsening();

    SetReceiveCreated(false);
    m_pSM3D_2->clearQueue();
    m_pSM3D_2->swapAllFaces();
    m_pSM3D_2->clearQueue();

    AnisoQualityEdgeInsertion();

    return 0;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::AnisoQualityEdgeInsertion()
  {

    GR_index_t Add = 0;
    GR_index_t TotalAdd = 0;
    SetQueuesPriority(eQualityEdge);

    m_pSM3D->clearQueue();
    m_pSM3D_2->clearQueue();
    m_IQIE.emptyQueue();
    m_IQBE.emptyQueue();
    queueAll();

    SetReceiveCreated(true);
    Add += AnisoRefineAllQueued();
    SetReceiveCreated(false);

    m_pSM3D->swapAllQueuedFaces();
    m_pSM3D->clearQueue();
    m_pSM3D_2->swapAllQueuedFacesRecursive();
    m_pSM3D_2->clearQueue();

    TotalAdd += Add;
    assert(m_VMesh->isValid());

    return TotalAdd;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::AnisoMetricLengthEdgeInsertion()
  {

    GR_index_t TotalAdd = 0;
    SetQueuesPriority(eMetricLengthEdge);
    m_pSM3D->clearQueue();
    m_IQIE.emptyQueue();
    m_IQBE.emptyQueue();
    queueAll();
    SetReceiveCreated(false);

    TotalAdd = AnisoRefineAllQueued();
    m_pSM3D_2->swapAllQueuedFacesRecursive();
    m_pSM3D_2->clearQueue();
    m_VMesh->isValid();

    return TotalAdd;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::AnisoMetricEdgeCoarsening()
  {

    GR_index_t iRemoved = 0;
    GR_index_t iRemovedTotal = 0;
    SetReceiveCreated(true);
    SetQueuesPriority(eCoarsening);

    do {
      iRemoved = 0;
      m_IQIE.emptyQueue();
      m_IQBE.emptyQueue();

      for (unsigned i = 0; i < m_VMesh->getNumVerts(); i++)
	m_VMesh->getVert(i)->markToKeep();

      queueAll();
      iRemoved += AnisoRefineAllQueued();
      iRemovedTotal += iRemoved;

    }
    while (iRemoved > 100);

    return iRemovedTotal;
  }

  void
  AnisoEdgeRefinementManager3D::GlobalSmoothing()
  {

    GRUMMP::SmoothingManagerMesquite3D *pSMM3D =
	new GRUMMP::SmoothingManagerMesquite3D(m_VMesh, true);
    pSMM3D->vSetTargetfromMetric();
    pSMM3D->smoothAllVerts(1);
    delete pSMM3D;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::AnisoRefineAllQueued()
  {

    logMessage(1, "%d Bdry faces queued for refinement.\n",
	       m_IQBE.iQueueLength());
    logMessage(1, "%d Int faces queued for refinement.\n",
	       m_IQIE.iQueueLength());

    m_IQBE.vCleanQueue();
    m_IQIE.vCleanQueue();

    GR_index_t Verts;
    GR_index_t Totalverts = 0;

    do {
      Verts = 0;

      m_IQBE.vCleanQueue();

      while (!m_IQBE.qIsQueueEmpty()) {

	InsertionQueueEntry IQE = m_IQBE.IQETopValidEntry();
	assert(IQE.eType() == InsertionQueueEntry::eSubSeg);
	Vert* pV0 = IQE.pVVert(0);
	Vert* pV1 = IQE.pVVert(1);

	if (!IQE.qStillInMesh()) {
	  m_IQBE.vPopTopEntry();
	  continue;
	}

	bool success;
	if (m_IQBE.getPriority() == InsertionQueue::eShortestEdge) {
	  assert(m_VMesh->isValid());
	  success = m_pAR->bRemoveVertexwithQual(pV0, pV1);
	}
	else
	  success = splitBdryEdge(pV0, pV1);

	if (success)
	  ++Verts;

	m_IQBE.vPopTopEntry();
	m_VMesh->sendEvents();
      }

      m_IQIE.vCleanQueue();

      while (!m_IQIE.qIsQueueEmpty()) {

	InsertionQueueEntry IQE = m_IQIE.IQETopValidEntry();
	assert(IQE.eType() == InsertionQueueEntry::eSubSeg);
	Vert* pV0 = IQE.pVVert(0);
	Vert* pV1 = IQE.pVVert(1);

	if (!IQE.qStillInMesh()) {
	  m_IQIE.vPopTopEntry();
	  continue;
	}

	bool success;
	if (m_IQIE.getPriority() == InsertionQueue::eShortestEdge) {
	  assert(m_VMesh->isValid());
	  success = m_pAR->bRemoveVertexwithQual(pV0, pV1);
	}
	else
	  success = splitEdge(pV0, pV1);

	if (success)
	  ++Verts;

	m_IQIE.vPopTopEntry();
	m_VMesh->sendEvents();
      }

      Totalverts += Verts;

    }
    while (!m_IQBE.qIsQueueEmpty() || !m_IQIE.qIsQueueEmpty());
    return Totalverts;
  }

  void
  AnisoEdgeRefinementManager3D::SetQueuesPriority(enum eInsertType Type)
  {

    switch (Type)
      {

      case eQualityEdge: //Edge Insertion based on quality

	m_IQBE.setPriority(InsertionQueue::eBestImpr);
	m_IQBE.setIsoQMetric(GRUMMP::TMOPQual::eShape);
	m_IQBE.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

	m_IQIE.setPriority(InsertionQueue::eBestImpr);
	m_IQIE.setIsoQMetric(GRUMMP::TMOPQual::eShape);
	m_IQIE.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

	break;

      case eMetricLengthEdge: ////Edge Insertion based on a max metric length

	m_IQBE.setPriority(InsertionQueue::eLongestEdge);
	m_IQBE.setIsoQMetric(GRUMMP::TMOPQual::eShape);
	m_IQBE.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

	m_IQIE.setPriority(InsertionQueue::eLongestEdge);
	m_IQIE.setIsoQMetric(GRUMMP::TMOPQual::eShape);
	m_IQIE.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

	break;

      case eCoarsening: ////Coarsening based on min metric length

	m_IQBE.setPriority(InsertionQueue::eShortestEdge);
	m_IQBE.setIsoQMetric(GRUMMP::TMOPQual::eShape);
	m_IQBE.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

	m_IQIE.setPriority(InsertionQueue::eShortestEdge);
	m_IQIE.setIsoQMetric(GRUMMP::TMOPQual::eShape);
	m_IQIE.setAnisoQMetric(GRUMMP::TMOPQual::eShapeOrient);

	break;
      }
  }

  void
  AnisoEdgeRefinementManager3D::receiveCreatedFaces(
      std::vector<Face*>& createdFaces)
  {

    if (bReceiveCreated == true) {
      logMessage(
	  3,
	  "AnisoRefinementManager2D received word about %zu newly created bfaces.\n",
	  createdFaces.size());

      std::vector<Face*>::iterator iter = createdFaces.begin(), iterEnd =
	  createdFaces.end();
      for (; iter != iterEnd; ++iter) {
	Face* pF = dynamic_cast<Face*>(*iter);
	if (pF->isDeleted() || pF->getFaceLoc() != Face::eInterior)
	  continue;
	if (!pF->doFullCheck())
	  continue;

	Vert *pV0 = pF->getVert(0);
	Vert *pV1 = pF->getVert(1);
	Vert *pV2 = pF->getVert(2);
	addToPriorityQueueEdge(pV0, pV1);
	addToPriorityQueueEdge(pV0, pV2);
	addToPriorityQueueEdge(pV1, pV2);
      }
    }
  }

  void
  AnisoEdgeRefinementManager3D::receiveCreatedBFaces(
      std::vector<BFace*>& createdBFaces)
  {

    if (bReceiveCreated == true) {
      logMessage(
	  3,
	  "AnisoRefinementManager2D received word about %zu newly created bfaces.\n",
	  createdBFaces.size());

      std::vector<BFace*>::iterator iter = createdBFaces.begin(), iterEnd =
	  createdBFaces.end();
      for (; iter != iterEnd; ++iter) {

	BFace* pBFace = *iter;
	if (pBFace->isDeleted())
	  continue;
	if (!pBFace->doFullCheck())
	  continue;

	Vert *pV0 = pBFace->getVert(0);
	Vert *pV1 = pBFace->getVert(1);
	Vert *pV2 = pBFace->getVert(2);
	addToPriorityQueueBdryEdge(pV0, pV1);
	addToPriorityQueueBdryEdge(pV0, pV2);
	addToPriorityQueueBdryEdge(pV1, pV2);
      }
    }
  }

  bool
  AnisoEdgeRefinementManager3D::splitEdge(Vert * const pV0, Vert* const pV1)
  {

    double insertionPoint[3];
    std::set<Vert*> setpV;

    setpV.insert(pV0);
    setpV.insert(pV1);

    InsertionPointCalculator::ePointType pointType = m_IPC->calcInsertionPoint(
	pV0, pV1, insertionPoint);
    assert(pointType == InsertionPointCalculator::eAnisoEdge);

    Vert* pNVert = m_SI3D->insertPointOnEdge(insertionPoint, pV0, pV1);

    InterpMetric(pNVert, setpV, m_Type);
    pNVert->setType(Vert::eInterior);

    return true;
  }

  bool
  AnisoEdgeRefinementManager3D::splitBdryEdge(Vert * const pV0, Vert* const pV1)
  {

    double insertionPoint[3];
    std::set<Vert*> setpV;

    setpV.insert(pV0);
    setpV.insert(pV1);

    InsertionPointCalculator::ePointType pointType = m_IPC->calcInsertionPoint(
	pV0, pV1, insertionPoint);
    assert(pointType == InsertionPointCalculator::eAnisoEdge);

    Vert* pNVert = m_SI3D->insertPointOnEdge(insertionPoint, pV0, pV1);
    InterpMetric(pNVert, setpV, m_Type);

    if (pV0->getVertType() == Vert::eBdryApex
	&& pV1->getVertType() == Vert::eBdryApex) {
      pNVert->setType(Vert::eBdryApex);
    }
    else if (pV0->getVertType() == Vert::eBdryCurve
	&& pV1->getVertType() == Vert::eBdryCurve) {
      pNVert->setType(Vert::eBdryCurve);
    }
    else
      pNVert->setType(Vert::eBdry);

    return true;
  }

  bool
  AnisoEdgeRefinementManager3D::addToPriorityQueueBdryEdge(
      Vert* const pV0, Vert* const pV1/*Edge*/)
  {

    if (pV0->isDeleted() || pV1->isDeleted())
      return false;

    std::set<Cell*> spCJunk;
    std::set<Vert*> spVNeigh;
    std::set<Vert*>::iterator iterV, iterVEnd;
    findNeighborhoodInfo(pV0, spCJunk, spVNeigh);

    iterV = spVNeigh.begin();
    iterVEnd = spVNeigh.end();

    bool bEdge = false;

    for (; iterV != iterVEnd; ++iterV) {
      Vert* psetV = (*iterV);
      if (psetV == pV1)
	bEdge = true;
    }

    if (bEdge == false)
      return false;

    if ((pV0->getVertType() == Vert::eBdry
	|| pV0->getVertType() == Vert::eBdryApex
	|| pV0->getVertType() == Vert::eBdryCurve)
	&& (pV1->getVertType() == Vert::eBdry
	    || pV1->getVertType() == Vert::eBdryApex
	    || pV1->getVertType() == Vert::eBdryCurve)) {
      return m_IQBE.qAddEntryAnisoEdge(InsertionQueueEntry(pV0, pV1));
    }
    return false;
  }

  bool
  AnisoEdgeRefinementManager3D::addToPriorityQueueEdge(Vert* const pV0,
						       Vert* const pV1/*Edge*/)
  {

    if (pV0->isDeleted() || pV1->isDeleted())
      return false;
    std::set<Cell*> spCJunk;
    std::set<Vert*> spVNeigh;
    std::set<Vert*>::iterator iterV, iterVEnd;
    findNeighborhoodInfo(pV0, spCJunk, spVNeigh);

    iterV = spVNeigh.begin();
    iterVEnd = spVNeigh.end();

    bool bEdge = false;

    for (; iterV != iterVEnd; ++iterV) {
      Vert* psetV = (*iterV);
      if (psetV == pV1)
	bEdge = true;
    }

    if (bEdge == false)
      return false;

    if ((pV0->getVertType() == Vert::eInterior
	|| pV1->getVertType() == Vert::eInterior)) {
      return m_IQIE.qAddEntryAnisoEdge(InsertionQueueEntry(pV0, pV1));
    }
    return false;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::queueAll()
  {

    GR_index_t EAdded = 0;
    GR_index_t BEAdded = 0;

    EAdded = queueEdges();
    BEAdded = queueBEdges();

    assert(
	EAdded + BEAdded
	    <= (this->m_IQIE.iQueueLength() + this->m_IQBE.iQueueLength()));
    return EAdded + BEAdded;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::queueBEdges()
  {
    GR_index_t BEdgesAdded = 0;
    for (GR_index_t i = 0; i < m_VMesh->getNumBdryFaces(); i++) {
      BFace* pBFace = m_VMesh->getBFace(i);

      assert(pBFace);
      if (pBFace->isDeleted())
	continue;

      Vert *pV0 = pBFace->getVert(0);
      Vert *pV1 = pBFace->getVert(1);
      Vert *pV2 = pBFace->getVert(2);

      if (addToPriorityQueueBdryEdge(pV0, pV1))
	BEdgesAdded++;

      if (addToPriorityQueueBdryEdge(pV1, pV2))
	BEdgesAdded++;

      if (addToPriorityQueueBdryEdge(pV2, pV0))
	BEdgesAdded++;
    }
    return BEdgesAdded;
  }

  GR_index_t
  AnisoEdgeRefinementManager3D::queueEdges()
  {
    GR_index_t EdgesAdded = 0;
    for (GR_index_t i = 0; i < m_VMesh->getNumFaces(); i++) {
      Face* pFace = m_VMesh->getFace(i);

      if (pFace->isDeleted())
	continue;

      Vert *pV0 = pFace->getVert(0);
      Vert *pV1 = pFace->getVert(1);
      Vert *pV2 = pFace->getVert(2);

      if (addToPriorityQueueEdge(pV0, pV1))
	EdgesAdded++;

      if (addToPriorityQueueEdge(pV1, pV2))
	EdgesAdded++;

      if (addToPriorityQueueEdge(pV2, pV0))
	EdgesAdded++;
    }
    return EdgesAdded;
  }

  void
  AnisoEdgeRefinementManager3D::LinearInterpMetric(Vert* pV,
						   std::set<Vert*> setpV)
  {

    std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();
    double dDist[setpV.size()];
    double dTDist = 0;
    double adNMetric[iNMetricValues];
    GR_index_t i = 0;

    for (; iter != iterEnd; ++iter) {
      Vert* psetV = (*iter);
      assert(psetV->isValid());
      dDist[i] = calcDistanceBetween(pV, psetV);
      dTDist += dDist[i];
      i++;
    }

    for (int ii = 0; ii < iNMetricValues; ii++)
      adNMetric[ii] = 0;

    iter = setpV.begin();
    iterEnd = setpV.end();

    i = 0;
    for (; iter != iterEnd; ++iter) {
      Vert* psetV = (*iter);
      for (int ii = 0; ii < iNMetricValues; ii++) {
	adNMetric[ii] += (dDist[i] / dTDist) * psetV->getMetric(ii);
      }
      i++;
    }

    pV->setMetric(adNMetric);

    return;
  }

  void
  AnisoEdgeRefinementManager3D::AverageMetric(Vert* pV, std::set<Vert*> setpV)
  {

    std::set<Vert*>::iterator iter = setpV.begin(), iterEnd = setpV.end();
    double adNMetric[iNMetricValues];

    for (int ii = 0; ii < iNMetricValues; ii++)
      adNMetric[ii] = 0;

    for (; iter != iterEnd; ++iter) {
      Vert* psetV = (*iter);
      for (int ii = 0; ii < iNMetricValues; ii++)
	adNMetric[ii] += (1. / setpV.size()) * psetV->getMetric(ii);
    }

    pV->setMetric(adNMetric);
    return;
  }

  void
  AnisoEdgeRefinementManager3D::InterpMetric(Vert* pV, std::set<Vert*> setpV,
					     enum eInterpType Type)
  {

    switch (Type)
      {
      case eDefault:
	AverageMetric(pV, setpV);
	break;
      case eAverage:
	AverageMetric(pV, setpV);
	break;
      case eLinear:
	LinearInterpMetric(pV, setpV);
	break;
      case eLogEuclidean:
	m_TQ->LogEuclideanInterpMetric(pV, setpV);
	break;
      }
    return;
  }
#endif
}
