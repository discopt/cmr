#include <cmr/tu.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "camion_internal.h"
#include "regular_internal.h"

#include <stdlib.h>
#include <assert.h>

#include "interface.h"

CMR_ERROR CMRtestTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_DEC** pdec,
  CMR_SUBMAT** psubmatrix, bool checkPlanarity, bool completeTree)
{
  assert(cmr);
  assert(matrix);

  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  CMR_CALL( CMRinterfaceTU(cmr, matrix, pisTotallyUnimodular, pdec, psubmatrix) );
  
  CMR_MINOR* minor = NULL;
//   CMR_CALL( CMRtestRegular(cmr, matrix, true, pisTotallyUnimodular, pdec, psubmatrix ? &minor : NULL, checkPlanarity,
//     completeTree) );
  if (minor)
  {
    assert(minor->numPivots == 0);
    *psubmatrix = minor->remainingSubmatrix;
    minor->remainingSubmatrix = NULL;
  }
  
  return CMR_OKAY;
}
