#include "iMesh.h"
#include "iMesh_GRUMMP.hh"

#include <vector>
#include "GR_Entity.h"
#include "GR_Geometry.h"
#include "GR_SurfMesh.h"
#include "iMesh_GRUMMP.hh"
#include "iMesh_GRUMMP_misc.hh"
#include "iMesh_GRUMMP_EntIter.hh"
#include "iMesh_GRUMMP_Workset.hh"

// Only three functions are different for surface meshes versus 2D meshes....

using namespace ITAPS_GRUMMP;

iMesh_Surface::iMesh_Surface(SurfMesh *pSM) :
    iMesh_Plane(eSurface)
{
  pM = pSM;
  externalData = true;
}

void
iMesh_Surface::vCreateMesh()
{
  SurfMesh *pSM = new SurfMesh();
  pM = pSM;
}

void
iMesh_Surface::load(const char *name, const char *, int* err)
{
  m_err = *err = iBase_SUCCESS;
  SurfMesh *pSM = new SurfMesh(name);
  pM = pSM;
}

void
iMesh_Surface::save(const char *name, const char *, int* err)
{
  m_err = *err = iBase_SUCCESS;
  SurfMesh *pSM = dynamic_cast<SurfMesh*>(pM);
  writeFile_Surface(*pSM, name);
}

