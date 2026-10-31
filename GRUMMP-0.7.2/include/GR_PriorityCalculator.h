/*
 * GR_PriorityCalculator.h
 *
 *  Created on: 2013-06-11
 *      Author: cfog
 */

#ifndef GR_PRIORITYCALCULATOR_H_
#define GR_PRIORITYCALCULATOR_H_

#include <memory>

#include "GR_QualMeasure.h"

namespace GRUMMP
{
  class PriorityCalculator {
  public:
    PriorityCalculator()
    {
    }
    virtual
    ~PriorityCalculator()
    {
    }

    virtual double
    calcPriority(const SimplexCell* const cell) const = 0;
  };

  class IsotropicPrioCalc2D : public PriorityCalculator {
  public:
    IsotropicPrioCalc2D(std::shared_ptr<GRUMMP::QualMeasure> qual,
			std::shared_ptr<Length2D> len) :
	m_qual(qual), m_len(len)
    {
    }
    ~IsotropicPrioCalc2D()
    {
    }
    virtual double
    calcPriority(const SimplexCell* const cell) const;
  private:
    std::shared_ptr<GRUMMP::QualMeasure> m_qual;
    std::shared_ptr<Length2D> m_len;
  };
}

#endif /* GR_PRIORITYCALCULATOR_H_ */
