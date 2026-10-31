/*
 * Length.cxx
 *
 *  Created on: Jun 13, 2017
 *      Author: cfog
 */

#include <memory>

#include "GR_ADT.h"
#include "GR_Length.h"
#include "GR_Mesh.h"

namespace GRUMMP
{
  const double Length::defaultMinLength = 1.e-12;

  Length::Length(const double& refine, const double& grading,
		 const double& minLength) :
      m_refine(refine), m_grade(grading), m_minLength(minLength), m_lengthMesh(), m_BBoxTree(
	  nullptr), m_interpolateLS(false), m_checkLength()
  {
  }

  Length::~Length()
  {
    delete m_BBoxTree;
  }

  std::shared_ptr<Length>
  Length::Create(const Mesh* const & pM, const double& refine,
		 const double& grading, const double& minLength)
  {
    std::shared_ptr<Length> len;
    if (pM->getVert(0)->getSpaceDimen() == 2) {
      std::shared_ptr<Length> temp(
	  new Length2D(pM, refine, grading, minLength));
      len = temp;
    }
    else {
      std::shared_ptr<Length> temp(
	  new Length3D(pM, refine, grading, minLength));
      len = temp;
    }
    // Length is already initialized and copied to the source mesh, so don't
    // do it again.
    return len;
  }

void Length::initLength() {
	// Find the RMS length of all edges incident on each vertex.
	assert(m_lengthMesh);
	GR_index_t numVerts = m_lengthMesh->getNumVerts();
	Neighborhood neigh;
	for (GR_index_t iV = 0; iV < numVerts; iV++) {
		Vert *pV = m_lengthMesh->getVert(iV);
		neigh.getNeighborhood(pV);
		double len = LARGE_DBL;
		for (Vert* pVCand : neigh.verts) {
			double temp = dDIST3D(pV->getCoords(), pVCand->getCoords());
			len = std::min(len, temp);
		}
		if (pV->isBdryVert()) {
			double len2 = distToNearestBdryEntity(pV, neigh.cells);
			len = std::min(len, len2);
		}
		double size = 0;
		for (Cell* pCCand : neigh.cells) {
			size += pCCand->calcSize();
		}
		len = std::max(len, m_minLength);
		pV->setLengthScale(len / m_refine);
		m_checkLength.insert(pV);
	}
	setGradedLength();
}

  void
  Length::setGradedLength()
  {
    Length::Neighborhood neigh;

    logMessage(MSG_DEBUG, "Grading length scales; g = %f\n", m_grade);

    while (!m_checkLength.empty()) {
      std::set<Vert*>::iterator iter = m_checkLength.begin();
      Vert *vertex = *iter;
      m_checkLength.erase(iter);

      neigh.getNeighborhood(vertex);

      double length_scale = vertex->getLengthScale();
      double Lipschitz_length = LARGE_DBL;

      std::set<Vert*>::iterator it = neigh.verts.begin(), it_end =
	  neigh.verts.end();
      for (; it != it_end; ++it) {
	Vert *neigh_vertex = *it;
	double neigh_ls = neigh_vertex->getLengthScale();
	double neigh_dist = calcDistanceBetween(vertex, neigh_vertex);

	assert(iFuzzyComp(neigh_ls, 0.) == 1);
	assert(iFuzzyComp(neigh_dist, 0.) == 1);

	Lipschitz_length = std::min(Lipschitz_length,
				    neigh_ls + (neigh_dist / m_grade));
      }

      assert(iFuzzyComp(m_minLength, 0.) == 1);
      assert(iFuzzyComp(Lipschitz_length, 0.) == 1);

      double new_length_scale = std::max(
	  m_minLength, std::min(Lipschitz_length, length_scale));

      if (new_length_scale < length_scale) {
	// If we're here, then we need to update the length scale and
	// consider checking the neighbors, too.
	vertex->setLengthScale(new_length_scale);

	for (it = neigh.verts.begin(); it != it_end; ++it) {
	  Vert* neigh_vertex = *it;
	  double neigh_ls = neigh_vertex->getLengthScale();
	  if (neigh_ls > new_length_scale)
	    m_checkLength.insert(neigh_vertex);
	}
      }
    }
  }

  bool
  Length::isInterpolatedLS() const
  {
    return m_interpolateLS;
  }

  void
  Length::setInterpolateLS(bool interpolateLs)
  {
    m_interpolateLS = interpolateLs;
  }

  void
  Length::Neighborhood::getNeighborhood(Vert* const vertex)
  {
    my_vert = vertex;

    if (vertex->isBdryVert()) {
      bool bdry_vert = true;
      findNeighborhoodInfo(vertex, cells, verts, &bfaces, &bdry_vert);
    }
    else {
      findNeighborhoodInfo(vertex, cells, verts);
    }
  }

  void
  Length::setLengthScale(Vert* const vertex)
  {
    vertex->setLengthScale(queryLengthScale(vertex->getCoords()));
  }

  void
  Length::setLengthScale(Mesh* const mesh)
  {
    for (GR_index_t iV = 0; iV < mesh->getNumVerts(); iV++) {
      Vert * pV = mesh->getVert(iV);
      setLengthScale(pV);
    }
  }

  void
  Length::scale(double dLengthRatio)
  {
	m_refine *= dLengthRatio;
	m_grade *= dLengthRatio;
    const Mesh* const mesh = m_lengthMesh.get();
    for (GR_index_t iV = 0; iV < mesh->getNumVerts(); iV++) {
      Vert *pV = mesh->getVert(iV);
      double len = std::max(m_minLength, pV->getLengthScale() * dLengthRatio);
      pV->setLengthScale(len);
    }
  }

} // namespace GRUMMP

