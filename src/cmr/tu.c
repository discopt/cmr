#include <cmr/tu.h>

#include <cmr/camion.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "camion_internal.h"
#include "regular_internal.h"

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

  fprintf(stream, "%stotal: %ld in %f seconds\n", prefix, stats->totalCount,
    stats->totalTime);

  return CMR_OKAY;
}

CMR_ERROR CMRtestTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_DEC** pdec,
  CMR_SUBMAT** psubmatrix, CMR_TU_PARAMETERS* params, CMR_TU_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);

  CMR_TU_PARAMETERS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRparamsTotalUnimodularityInit(&defaultParams) );
    params = &defaultParams;
  }

  clock_t totalClock = 0;
  if (stats)
    totalClock = clock();
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  CMR_CALL( CMRtestCamionSigned(cmr, matrix, pisTotallyUnimodular, psubmatrix, stats ? &stats->camion : NULL) );
  if (!*pisTotallyUnimodular)
    return CMR_OKAY;

  CMR_MINOR* minor = NULL;
  // TODO: run regularity check with ternary = true.
  CMR_CALL( CMRtestRegular(cmr, matrix, false, pisTotallyUnimodular, pdec, psubmatrix ? &minor : NULL,
    &params->regular, stats ? &stats->regular : NULL) );

  if (minor)
  {
    assert(minor->numPivots == 0);
    *psubmatrix = minor->remainingSubmatrix;
    minor->remainingSubmatrix = NULL;
  }

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}
