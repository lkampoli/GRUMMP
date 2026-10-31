#ifndef GR_Util_h
#define GR_Util_h

#include <math.h>

#include <assert.h>

#define VALIDATE_INPUT(expr) assert(expr);

namespace GRUMMP
{
  static double epsilon = 1.e-12;
  inline bool
  fuzzyGreater(const double a, const double b)
  {
    double sum = fabs(a) + fabs(b);
    double diff = a - b;
    if (sum < epsilon)
      return false;
    else
      return (diff / sum > epsilon);
  }
}

#endif
