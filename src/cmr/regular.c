// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "env_internal.h"
#include "seymour_internal.h"

CMR_ERROR CMRregularParamsInit(CMR_REGULAR_PARAMS* params)
{
  assert(params);

  CMR_CALL( CMRseymourParamsInit(&(params->seymour)) );

  return CMR_OKAY;
}

CMR_ERROR CMRregularStatsInit(CMR_REGULAR_STATS* stats)
{
  assert(stats);

  CMR_CALL( CMRseymourStatsInit(&(stats->seymour)) );

  return CMR_OKAY;
}

CMR_ERROR CMRregularStatsPrint(FILE* stream, CMR_REGULAR_STATS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  CMR_CALL( CMRseymourStatsPrint(stream, &(stats->seymour), prefix) );

  return CMR_OKAY;
}


CMR_ERROR CMRregularTest(CMR* cmr, CMR_CHRMAT* matrix, bool *pisRegular, CMR_SEYMOUR_NODE** proot,
  CMR_MINOR** pminor, CMR_REGULAR_PARAMS* params, CMR_REGULAR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);

  CMR_REGULAR_PARAMS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRregularParamsInit(&defaultParams) );
    params = &defaultParams;
  }

  CMR_SUBMAT* submatrix = NULL;
  if (!CMRchrmatIsBinary(cmr, matrix, pminor ? &submatrix : NULL))
  {
    if (pisRegular)
      *pisRegular = false;
    if (pminor)
      CMR_CALL( CMRminorCreate(cmr, pminor, 0, submatrix, CMR_MINOR_TYPE_ENTRY) );
    return CMR_OKAY;
  }

  CMR_SEYMOUR_NODE* root = NULL;
  CMR_CALL( CMRseymourDecompose(cmr, matrix, false, &root, &(params->seymour), stats ? &(stats->seymour) : NULL,
    timeLimit) );
  int8_t regularity = CMRseymourRegularity(root);
  if (regularity)
    *pisRegular = regularity > 0;
  if (proot)
    *proot = root;
  else
    CMR_CALL( CMRseymourRelease(cmr, &root) );

  return CMR_OKAY;
}

CMR_ERROR CMRregularCompleteDecomposition(CMR* cmr, CMR_SEYMOUR_NODE* dec, CMR_REGULAR_PARAMS* params,
  CMR_REGULAR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(dec);
  assert(timeLimit > 0);

  CMR_REGULAR_PARAMS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRregularParamsInit(&defaultParams) );
    params = &defaultParams;
  }

  CMR_CALL( CMRregularityCompleteDecomposition(cmr, dec, &(params->seymour), stats ? &(stats->seymour) : NULL,
    timeLimit) );

  return CMR_OKAY;
}

CMR_ERROR CMRregularRefineDecomposition(CMR* cmr, size_t numNodes, CMR_SEYMOUR_NODE** nodes, CMR_REGULAR_PARAMS* params,
  CMR_REGULAR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(timeLimit > 0);

  if (numNodes == 0)
    return CMR_OKAY;
  assert(nodes);

  CMR_REGULAR_PARAMS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRregularParamsInit(&defaultParams) );
    params = &defaultParams;
  }

  CMR_CALL( CMRregularityRefineDecomposition(cmr, numNodes, nodes, &(params->seymour), stats ? &(stats->seymour) : NULL,
    timeLimit) );

  return CMR_OKAY;
}

