#include <cmr/tu.h>

#include <cmr/camion.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "camion_internal.h"
#include "regular_internal.h"
#include "hereditary_property.h"

#include <stdlib.h>
#include <assert.h>
#include <time.h>

CMR_ERROR CMRparamsTotalUnimodularityInit(CMR_TU_PARAMETERS* params)
{
  assert(params);

  CMR_CALL( CMRparamsRegularInit(&params->regular) );

  return CMR_OKAY;
}

CMR_ERROR CMRstatsTotalUnimodularityInit(CMR_TU_STATISTICS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  CMR_CALL( CMRstatsCamionInit(&stats->camion) );
  CMR_CALL( CMRstatsRegularInit(&stats->regular) );

  return CMR_OKAY;
}

CMR_ERROR CMRstatsTotalUnimodularityPrint(FILE* stream, CMR_TU_STATISTICS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Total unimodularity recognition:\n");
    prefix = "  ";
  }

  char subPrefix[256];
  snprintf(subPrefix, 256, "%scamion ", prefix);
  CMR_CALL( CMRstatsCamionPrint(stream, &stats->camion, subPrefix) );
  snprintf(subPrefix, 256, "%sregularity ", prefix);
  CMR_CALL( CMRstatsRegularPrint(stream, &stats->regular, subPrefix) );

  fprintf(stream, "%stotal: %ld in %f seconds\n", prefix, stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}

static
CMR_ERROR tuTest(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Some matrix to be tested for total unimodularity. */
  void* data,                 /**< Additional data (must be \c NULL). */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \p matrix is totally unimodular. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a proper non-totally unimodular submatrix of \p matrix. */
  double timeLimit                /**< Time limit to impose. */
)
{
  assert(cmr);
  assert(matrix);
  assert(pisTotallyUnimodular);
  assert(!psubmatrix || !*psubmatrix);

  CMR_TU_STATISTICS* stats = (CMR_TU_STATISTICS*) data;

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "tuTest called for a %dx%d matrix\n", matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  *pisTotallyUnimodular = true;
  clock_t time = clock();

  CMR_CALL( CMRtestCamionSigned(cmr, matrix, pisTotallyUnimodular, NULL, stats ? &stats->camion : NULL,
    timeLimit) );

  if (*pisTotallyUnimodular)
  {
    CMR_REGULAR_PARAMETERS params;
    CMR_CALL( CMRparamsRegularInit(&params) );
    double remainingTime = timeLimit - ((clock() - time) * 1.0 / CLOCKS_PER_SEC);
    CMR_CALL( CMRtestRegular(cmr, matrix, false, pisTotallyUnimodular, NULL, NULL, &params,
      stats ? &stats->regular : NULL, remainingTime) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_DEC** pdec,
  CMR_SUBMAT** psubmatrix, CMR_TU_PARAMETERS* params, CMR_TU_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );

  CMR_TU_PARAMETERS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRparamsTotalUnimodularityInit(&defaultParams) );
    params = &defaultParams;
  }

  clock_t totalClock = clock();
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  CMR_CALL( CMRtestCamionSigned(cmr, matrix, pisTotallyUnimodular, psubmatrix, stats ? &stats->camion : NULL,
    timeLimit) );

  if (!*pisTotallyUnimodular)
  {
    if (stats)
    {
      stats->totalCount++;
      stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
    }
    return CMR_OKAY;
  }

  // TODO: run regularity check with ternary = true.

  double remainingTime = timeLimit - ((clock() - totalClock) * 1.0 / CLOCKS_PER_SEC);
  CMR_CALL( CMRtestRegular(cmr, matrix, false, pisTotallyUnimodular, pdec, NULL, &params->regular,
    stats ? &stats->regular : NULL, remainingTime) );

  if (!*pisTotallyUnimodular && psubmatrix)
  {
    assert(!*psubmatrix);
    remainingTime = timeLimit - (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRtestHereditaryPropertySimple(cmr, matrix, tuTest, stats, psubmatrix, remainingTime) );
  }

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}
