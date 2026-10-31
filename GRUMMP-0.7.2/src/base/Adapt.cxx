#include "GR_CoarsenManager.h"
#include "GR_InsertionQueue.h"
#include "GR_InsertionManager.h"
#include "GR_SmoothingManager.h"
#include "GR_Mesh.h"

void
Mesh::adaptToLengthScale(const eEncroachType eET,
			 std::shared_ptr<GRUMMP::Length> const sizing_field)
{
  logMessage(1, "Making mesh Delaunay\n");
  bool qResult = makeDelaunay();
  if (!qResult) {
    vWarning("Couldn't make mesh fully Delaunay before coarsening!\n");
    vWarning("This MIGHT still turn out okay...\n");
  }

  logMessage(1, "Coarsening\n");
  GRUMMP::CoarsenManager CM(this, sizing_field);
  CM.coarsenToLengthScale();

  logMessage(1, "Ensuring that mesh is Delaunay\n");
  qResult = makeDelaunay();
  if (!qResult) {
    vFatalError(
	"Couldn't make mesh Delaunay!  Insertion will fail; bailing out!\n",
	"Mesh::vAdaptToLengthScale");
  }

  logMessage(1, "Refining\n");
  setEncroachmentType(eET);
  if (getType() == Mesh::eMesh2D) {
    Mesh2D *pM2D = dynamic_cast<Mesh2D*>(this);
    assert(pM2D);
    GRUMMP::InsertionPointCalculator *pIPC = new GRUMMP::CircumcenterIPC();
    GRUMMP::RefinementManager2D RM2D(pM2D, pIPC, sizing_field);
    RM2D.setQualityTarget(2.0 * sin(M_PI / 6) / sqrt(3));
    RM2D.setRefineForSize(true);
    RM2D.refineMesh();
  }
  else {
    // Unlikely to work, even in 3D.  That is, I think the use of
    // InsertionQueue to drive refinement is actually dead now.

    // Set up insertion queue and insert.
//    InsertionQueue IQ(this);
//      IQ.vQualityRefine();
  }
  purgeAllEntities();

  GRUMMP::OptMSSmoothingManager* OMS = GRUMMP::OptMSSmoothingManager::Create(
      this);
  if (OMS) {
    OMS->smoothAllVerts(2);
    delete OMS;
  }
  assert(isValid());
}
