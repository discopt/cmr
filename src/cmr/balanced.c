#include <cmr/balanced.h>

#include <cmr/camion.h>

#include "matrix_internal.h"
#include "camion_internal.h"

#include <stdlib.h>
#include <assert.h>
#include <time.h>

CMR_ERROR CMRtestBalancedEnumeration(CMR* cmr, CMR_CHRMAT* matrix, bool* pisBalanced, CMR_SUBMAT** psubmatrix,
  double timeLimit)
{
  assert(cmr);
  assert(matrix);

  // clock_t totalClock = clock();
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  printf("DEBUG: CMRtestBalancedEnumeration called!\n");

  return CMR_OKAY;
}

CMR_ERROR CMRtestBalancedGraph(CMR* cmr, CMR_CHRMAT* matrix, bool* pisBalanced, CMR_SUBMAT** psubmatrix,
  double timeLimit)
{
  assert(cmr);
  assert(matrix);

  // clock_t totalClock = clock();
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  printf("DEBUG: CMRtestBalancedGraph called!\n");

  return CMR_OKAY;
}
CMR_ERROR CMRtestBalanced(CMR* cmr, CMR_CHRMAT* matrix, bool* pisBalanced, CMR_SUBMAT** psubmatrix,
  double timeLimit)
{
  assert(cmr);
  assert(matrix);

  // clock_t totalClock = clock();
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  printf("DEBUG: CMRtestBalanced called!\n");

  return CMR_OKAY;
}
