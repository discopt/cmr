// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "env_internal.h"
#include "matroid_internal.h"
#include "regularity_internal.h"

CMR_ERROR CMRregularParamsInit(CMR_REGULAR_PARAMS* params)
{
  assert(params);

  params->directGraphicness = true;
  params->seriesParallel = true;
  params->planarityCheck = false;
  params->treeFlags = CMR_REGULAR_TREE_FLAGS_DEFAULT;
  params->threeSumPivotChildren = false;
  params->threeSumStrategy = CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS /* TODO: Later no pivots. */
    | CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_WIDE | CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_MIXED
    | CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_WIDE | CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_MIXED;
  params->graphs = CMR_DEC_CONSTRUCT_NONE;

  return CMR_OKAY;
}

CMR_ERROR CMRregularStatsInit(CMR_REGULAR_STATS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  CMR_CALL( CMRspStatsInit(&stats->seriesParallel) );
  CMR_CALL( CMRgraphicStatsInit(&stats->graphic) );
  CMR_CALL( CMRnetworkStatsInit(&stats->network) );
  CMR_CALL( CMRcamionStatsInit(&stats->camion) );
  stats->sequenceExtensionCount = 0;
  stats->sequenceExtensionTime = 0.0;  
  stats->sequenceGraphicCount = 0;
  stats->sequenceGraphicTime = 0.0;
  stats->enumerationCount = 0;
  stats->enumerationTime = 0.0;
  stats->enumerationCandidatesCount = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRregularStatsPrint(FILE* stream, CMR_REGULAR_STATS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Regular matrix recognition:\n");
    prefix = "  ";
  }

  char subPrefix[256];
  snprintf(subPrefix, 256, "%sseries-parallel ", prefix);
  CMR_CALL( CMRspStatsPrint(stream, &stats->seriesParallel, subPrefix) );
  snprintf(subPrefix, 256, "%s(co)graphic ", prefix);
  CMR_CALL( CMRgraphicStatsPrint(stream, &stats->graphic, subPrefix) );
  snprintf(subPrefix, 256, "%s(co)network ", prefix);
  CMR_CALL( CMRnetworkStatsPrint(stream, &stats->network, subPrefix) );
  snprintf(subPrefix, 256, "%scamion ", prefix);
  CMR_CALL( CMRcamionStatsPrint(stream, &stats->camion, subPrefix) );

  fprintf(stream, "%ssequence extensions: %lu in %f seconds\n", prefix, (unsigned long)stats->sequenceExtensionCount,
    stats->sequenceExtensionTime);
  fprintf(stream, "%ssequence (co)graphic: %lu in %f seconds\n", prefix, (unsigned long)stats->sequenceGraphicCount,
    stats->sequenceGraphicTime);
  fprintf(stream, "%senumeration: %lu in %f seconds\n", prefix, (unsigned long)stats->enumerationCount,
    stats->enumerationTime);
  fprintf(stream, "%s3-separation candidates: %lu in %f seconds\n", prefix,
    (unsigned long)stats->enumerationCandidatesCount, stats->enumerationTime);
  fprintf(stream, "%stotal: %lu in %f seconds\n", prefix, (unsigned long)stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}


CMR_ERROR CMRregularTest(CMR* cmr, CMR_CHRMAT* matrix, bool *pisRegular, CMR_SEYMOUR_NODE** pdec,
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
      CMR_CALL( CMRminorCreate(cmr, pminor, 0, submatrix) );
    return CMR_OKAY;
  }

  CMR_CALL( CMRregularityTest(cmr, matrix, false, pisRegular, pdec, pminor, params, stats, timeLimit) );

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

  CMR_CALL( CMRregularityCompleteDecomposition(cmr, dec, params, stats, timeLimit) );

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

  CMR_CALL( CMRregularityRefineDecomposition(cmr, numNodes, nodes, params, stats, timeLimit) );

  return CMR_OKAY;
}

