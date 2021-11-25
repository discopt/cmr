#include <cmr/tu.h>

#include <cmr/camion.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "camion_internal.h"
#include "regular_internal.h"

#include <stdlib.h>
#include <assert.h>

CMR_ERROR CMRtuInitParameters(CMR_TU_PARAMETERS* params)
{
  assert(params);

  CMRregularInitParameters(&params->regular);

  return CMR_OKAY;
}

CMR_ERROR CMRtestTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_DEC** pdec,
  CMR_SUBMAT** psubmatrix, CMR_TU_PARAMETERS* params)
{
  assert(cmr);
  assert(matrix);

  CMR_TU_PARAMETERS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRtuInitParameters(&defaultParams) );
    params = &defaultParams;
  }

  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  CMR_CALL( CMRtestCamionSigned(cmr, matrix, pisTotallyUnimodular, psubmatrix) );
  if (!*pisTotallyUnimodular)
    return CMR_OKAY;

  CMR_MINOR* minor = NULL;
  // TODO: run regularity check with ternary = true.
  CMR_CALL( CMRtestRegular(cmr, matrix, false, pisTotallyUnimodular, pdec, psubmatrix ? &minor : NULL,
    &params->regular, NULL) );
  if (minor)
  {
    assert(minor->numPivots == 0);
    *psubmatrix = minor->remainingSubmatrix;
    minor->remainingSubmatrix = NULL;
  }

  return CMR_OKAY;
}
