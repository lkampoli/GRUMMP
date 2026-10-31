/*
 * PriorityCalculator.cxx
 *
 *  Created on: 2013-06-11
 *      Author: cfog
 */

#include "GR_Cell.h"
#include "GR_CellPriorityQueue.h"
#include "GR_PriorityCalculator.h"
#include "GR_QualMeasure.h"
#include "GR_Vertex.h"

namespace GRUMMP
{
  double
  IsotropicPrioCalc2D::calcPriority(const SimplexCell* const cell) const
  {
    double cellRadius = cell->calcCircumradius();
    double size = M_SQRT2 * cellRadius;
    double correctSize = 0;
    for (int i = 0; i < 3; i++) {
      Vert const *vert = cell->getVert(i);
      correctSize += vert->getLengthScale();
    }
    correctSize /= 3;

    double qual = 0;
    if (size * 10 < correctSize || cell->isSmallAngleAffected()) {
      qual = CellPriorityQueue::dNeverQueuedPriority;
    }
    else if (size > correctSize) {
      qual = -size / correctSize;
    }
    qual += m_qual->eval(cell->getVert(0), cell->getVert(1), cell->getVert(2));

    return qual;
  }
}
