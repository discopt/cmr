// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "env_internal.h"
#include "dec_internal.h"
#include "regular_internal.h"

CMR_ERROR CMRparamsRegularInit(CMR_REGULAR_PARAMETERS* params)
{
  assert(params);

  params->directGraphicness = true;
  params->seriesParallel = true;
  params->planarityCheck = false;
  params->completeTree = false;
  params->matrices = CMR_DEC_CONSTRUCT_NONE;
  params->transposes = CMR_DEC_CONSTRUCT_NONE;
  params->graphs = CMR_DEC_CONSTRUCT_NONE;

  return CMR_OKAY;
}

CMR_ERROR CMRstatsRegularInit(CMR_REGULAR_STATISTICS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  CMR_CALL( CMRstatsSeriesParallelInit(&stats->seriesParallel) );
  CMR_CALL( CMRstatsGraphicInit(&stats->graphic) );
  CMR_CALL( CMRstatsNetworkInit(&stats->network) );
  stats->sequenceExtensionCount = 0;
  stats->sequenceExtensionTime = 0.0;  
  stats->sequenceGraphicCount = 0;
  stats->sequenceGraphicTime = 0.0;
  stats->enumerationCount = 0;
  stats->enumerationTime = 0.0;
  stats->enumerationCandidatesCount = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRstatsRegularPrint(FILE* stream, CMR_REGULAR_STATISTICS* stats, const char* prefix)
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
  CMR_CALL( CMRstatsSeriesParallelPrint(stream, &stats->seriesParallel, subPrefix) );
  snprintf(subPrefix, 256, "%s(co)graphic ", prefix);
  CMR_CALL( CMRstatsGraphicPrint(stream, &stats->graphic, subPrefix) );
  snprintf(subPrefix, 256, "%s(co)network ", prefix);
  CMR_CALL( CMRstatsNetworkPrint(stream, &stats->network, subPrefix) );

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

/**
 * \brief Tests a 2-connected binary or ternary matrix for regularity.
 */

static
CMR_ERROR testRegularTwoConnected(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC* dec,                   /**< Decomposition node. */
  bool ternary,                   /**< Whether signs matter. */
  bool *pisRegular,               /**< Pointer for storing whether \p matrix is regular. */
  CMR_MINOR** pminor,             /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMETERS* params, /**< Parameters for the computation. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

static
CMR_ERROR testRegularThreeConnectedWithSequence(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC* dec,                   /**< Decomposition node. */
  bool ternary,                   /**< Whether signs matter. */
  bool *pisRegular,               /**< Pointer for storing whether \p matrix is regular. */
  CMR_MINOR** pminor,             /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMETERS* params, /**< Parameters for the computation. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
)
{
  assert(cmr);
  assert(dec);
  assert(dec->matrix);
  assert(dec->nestedMinorsMatrix);
  assert(dec->nestedMinorsRowsOriginal);
  assert(dec->nestedMinorsColumnsOriginal);
  assert(dec->nestedMinorsSequenceNumRows);
  assert(dec->nestedMinorsSequenceNumColumns);
  assert(dec->nestedMinorsLength > 0);

  clock_t time = clock();

  CMRdbgMsg(6, "Testing binary %dx%d 3-connected matrix with given nested sequence of 3-connected minors for regularity.\n",
    dec->matrix->numRows, dec->matrix->numColumns);

#if defined(CMR_DEBUG)
  CMR_CALL( CMRdecPrintSequenceNested3ConnectedMinors(cmr, dec, stdout) );
#endif /* CMR_DEBUG */

  CMRconsistencyAssert( CMRdecConsistency(dec, false) );

  CMR_CHRMAT* nestedMinorsTranspose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, dec->nestedMinorsMatrix, &nestedMinorsTranspose) );

  /* Test sequence for graphicness. */
  double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  size_t lastGraphicMinor = 0;
  CMR_GRAPH* graph = NULL;
  CMR_ELEMENT* graphEdgeLabels = NULL;
  CMR_CALL( CMRregularSequenceGraphic(cmr, dec->nestedMinorsMatrix, nestedMinorsTranspose,
    dec->nestedMinorsRowsOriginal, dec->nestedMinorsColumnsOriginal, dec->nestedMinorsLength,
    dec->nestedMinorsSequenceNumRows, dec->nestedMinorsSequenceNumColumns, &lastGraphicMinor, &graph,
    &graphEdgeLabels, stats, remainingTime) );

  if (graph)
  {
    assert(graphEdgeLabels);

    CMRdbgMsg(8, "Matrix is graphic.\n");
    dec->type = CMR_DEC_GRAPHIC;
    dec->graph = graph;
    dec->graphForest = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &dec->graphForest, dec->matrix->numRows) );
    dec->graphCoforest = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &dec->graphCoforest, dec->matrix->numColumns) );

    assert(CMRgraphNumEdges(graph) == dec->matrix->numRows + dec->matrix->numColumns);
    for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
      iter = CMRgraphEdgesNext(graph, iter))
    {
      CMR_GRAPH_EDGE edge = CMRgraphEdgesEdge(graph, iter);
      CMR_ELEMENT element = graphEdgeLabels[edge];
      if (CMRelementIsRow(element))
        element = dec->nestedMinorsRowsOriginal[CMRelementToRowIndex(element)];
      else
        element = dec->nestedMinorsColumnsOriginal[CMRelementToColumnIndex(element)];
      if (CMRelementIsRow(element))
        dec->graphForest[CMRelementToRowIndex(element)] = edge;
      else
        dec->graphCoforest[CMRelementToColumnIndex(element)] = edge;
    }

    CMR_CALL( CMRfreeBlockArray(cmr, &graphEdgeLabels) );
  }
  else
  {
    dec->flags &= ~CMR_DEC_IS_GRAPHIC;
  }

  size_t lastCographicMinor = 0;
  if (!dec->graph || params->planarityCheck)
  {
    /* Test sequence for cographicness. */
    CMR_GRAPH* cograph = NULL;
    CMR_ELEMENT* cographEdgeLabels = NULL;
    remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRregularSequenceGraphic(cmr, nestedMinorsTranspose, dec->nestedMinorsMatrix,
      dec->nestedMinorsColumnsOriginal, dec->nestedMinorsRowsOriginal, dec->nestedMinorsLength,
      dec->nestedMinorsSequenceNumColumns, dec->nestedMinorsSequenceNumRows, &lastCographicMinor, &cograph,
      &cographEdgeLabels, stats, remainingTime) );

    if (cograph)
    {
      assert(cographEdgeLabels);

      CMRdbgMsg(8, "Matrix is graphic.\n");
      dec->type = (dec->type == CMR_DEC_GRAPHIC) ? CMR_DEC_PLANAR : CMR_DEC_COGRAPHIC;
      dec->cograph = cograph;
      dec->cographForest = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &dec->cographForest, dec->matrix->numColumns) );
      dec->cographCoforest = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &dec->cographCoforest, dec->matrix->numRows) );

      assert(CMRgraphNumEdges(cograph) == dec->matrix->numRows + dec->matrix->numColumns);
      for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(cograph); CMRgraphEdgesValid(cograph, iter);
        iter = CMRgraphEdgesNext(cograph, iter))
      {
        CMR_GRAPH_EDGE edge = CMRgraphEdgesEdge(cograph, iter);
        CMR_ELEMENT element = cographEdgeLabels[edge];
        if (CMRelementIsRow(element))
          element = dec->nestedMinorsColumnsOriginal[CMRelementToRowIndex(element)];
        else
          element = dec->nestedMinorsRowsOriginal[CMRelementToColumnIndex(element)];
        if (CMRelementIsRow(element))
          dec->cographCoforest[CMRelementToRowIndex(element)] = edge;
        else
          dec->cographForest[CMRelementToColumnIndex(element)] = edge;
      }

      CMR_CALL( CMRfreeBlockArray(cmr, &cographEdgeLabels) );
    }
    else
    {
      dec->flags &= ~CMR_DEC_IS_COGRAPHIC;
    }
  }

  if (!dec->graph && !dec->cograph)
  {
    CMRdbgMsg(8, "Checking for R10.\n");
    bool isR10;
    CMR_CALL( CMRregularThreeConnectedIsR10(cmr, dec, &isR10) );

    if (!isR10)
    {
      CMRdbgMsg(8, "Neither graphic (until %ld) nor cographic (until %ld), nor R_10! First non-graphic and non-cographic minor is %ld.\n",
        lastGraphicMinor, lastCographicMinor,
        lastGraphicMinor > lastCographicMinor ? (lastGraphicMinor+1) : (lastCographicMinor + 1) );

      remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
      CMR_CALL( CMRregularSearchThreeSeparation(cmr, dec, nestedMinorsTranspose, ternary,
        lastGraphicMinor > lastCographicMinor ? (lastGraphicMinor+1) : (lastCographicMinor + 1), NULL, params, stats,
        remainingTime) );

      if (dec->type == CMR_DEC_IRREGULAR)
      {
        CMRdbgMsg(8, "Minor determined to be irregular.\n");
        *pisRegular = false;
      }
      else
      {
        CMRdbgMsg(8, " Encountered a 3-separation.\n");
        assert(dec->numChildren == 2);

#if defined(CMR_DEBUG)
        CMR_CALL( CMRchrmatPrintDense(cmr, dec->children[0]->matrix, stdout, '0', true) );
        CMR_CALL( CMRchrmatPrintDense(cmr, dec->children[1]->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

        remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
        CMR_CALL( testRegularTwoConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params, stats,
          remainingTime) );
        
        if (params->completeTree || *pisRegular)
        {
          remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
          CMR_CALL( testRegularTwoConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params, stats,
            remainingTime) );
        }
      }
    }
  }

  CMR_CALL( CMRchrmatFree(cmr, &nestedMinorsTranspose) );

  return CMR_OKAY;
}

static
CMR_ERROR testRegularTwoConnected(CMR* cmr, CMR_DEC* dec, bool ternary, bool *pisRegular, CMR_MINOR** pminor,
  CMR_REGULAR_PARAMETERS* params, CMR_REGULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(dec);
  assert(dec->matrix);

  CMRdbgMsg(2, "Testing binary %dx%d 2-connected matrix for regularity.\n", dec->matrix->numRows,
    dec->matrix->numColumns);

  clock_t time = clock();
  CMR_SUBMAT* submatrix = NULL;

  if (params->directGraphicness || dec->matrix->numRows <= 3 || dec->matrix->numColumns <= 3)
  {
    /* We run the almost-linear time algorithm. Otherwise, graphicness is checked later for the 3-connected components. */

    if (params->directGraphicness || params->planarityCheck || dec->matrix->numRows > 3)
    {
      CMRdbgMsg(4, "Checking for graphicness...");
      bool isGraphic;
      CMR_CALL( CMRregularTestGraphic(cmr, &dec->matrix, &dec->transpose, ternary, &isGraphic, &dec->graph,
        &dec->graphForest, &dec->graphCoforest, &dec->graphArcsReversed, &submatrix, stats, timeLimit) );
      if (isGraphic)
      {
        CMRdbgMsg(0, " graphic.\n");
        dec->type = CMR_DEC_GRAPHIC;
        if (!params->planarityCheck)
          return CMR_OKAY;
      }
      CMRdbgMsg(0, " NOT graphic.\n");
    }

    CMRdbgMsg(4, "Checking for cographicness...");
    bool isCographic;
    double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRregularTestGraphic(cmr, &dec->transpose, &dec->matrix, ternary, &isCographic, &dec->cograph,
      &dec->cographForest, &dec->cographCoforest, &dec->cographArcsReversed, &submatrix, stats, remainingTime) );
    if (isCographic)
    {
      CMRdbgMsg(0, " cographic.\n");
      dec->type = (dec->type == CMR_DEC_GRAPHIC) ? CMR_DEC_PLANAR : CMR_DEC_COGRAPHIC;
      return CMR_OKAY;
    }
    CMRdbgMsg(0, " NOT cographic.\n");

    if (submatrix)
      CMR_CALL( CMRsubmatTranspose(submatrix) );
  }

  if (dec->nestedMinorsMatrix)
  {
    CMRdbgMsg(4, "A partial sequence of nested minors is already known.\n");

#if defined(CMR_DEBUG)
    CMR_CALL( CMRdecPrintSequenceNested3ConnectedMinors(cmr, dec, stdout) );
#endif /* CMR_DEBUG */

    double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRregularExtendNestedMinorSequence(cmr, dec, ternary, &submatrix, params, stats, remainingTime) );
    
    /* Handling of the resulting sequence or 2-separation is done at the end. */
  }
  else
  {
    CMRdbgMsg(4, "Splitting off series-parallel elements...");
    double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRregularDecomposeSeriesParallel(cmr, &dec, ternary, &submatrix, params, stats, remainingTime) );

    if (dec->type == CMR_DEC_IRREGULAR)
    {
      CMRdbgMsg(0, " NOT regular.\n");
      return CMR_OKAY;
    }
    else if (dec->type == CMR_DEC_SERIES_PARALLEL)
    {
      CMRdbgMsg(0, " series-parallel.\n");
      return CMR_OKAY;
    }

    if (dec->type == CMR_DEC_TWO_SUM)
    {
      CMRdbgMsg(0, " Encountered a 2-separation.\n");
      assert(dec->numChildren == 2);
      remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
      CMR_CALL( testRegularTwoConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params, stats,
        remainingTime) );

      if (params->completeTree || *pisRegular)
      {
        remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
        CMR_CALL( testRegularTwoConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params, stats,
          remainingTime) );
      }

      return CMR_OKAY;
    }

    CMR_SUBMAT* wheelSubmatrix = submatrix;
    CMRdbgMsg(0, " Found a W_%d minor.\n", wheelSubmatrix->numRows);
    submatrix = NULL;

    /* No 2-sum found, so we have a wheel submatrix. */

    remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRregularConstructNestedMinorSequence(cmr, dec, ternary, wheelSubmatrix, &submatrix, params, stats,
      remainingTime) );
    CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  }

  if (dec->type == CMR_DEC_IRREGULAR)
  {
    CMRdbgMsg(0, " NOT regular.\n");
    return CMR_OKAY;
  }

  if (dec->type == CMR_DEC_TWO_SUM)
  {
    CMRdbgMsg(0, " Encountered a 2-separation.\n");
    assert(dec->numChildren == 2);

    double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( testRegularTwoConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params, stats,
      remainingTime) );

    if (params->completeTree || *pisRegular)
    {
      remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
      CMR_CALL( testRegularTwoConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params, stats,
        remainingTime) );
    }

    return CMR_OKAY;
  }

  double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  CMR_CALL( testRegularThreeConnectedWithSequence(cmr, dec, ternary, pisRegular, pminor, params, stats,
    remainingTime) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestRegular(CMR* cmr, CMR_CHRMAT* matrix, bool ternary, bool *pisRegular, CMR_DEC** pdec,
  CMR_MINOR** pminor, CMR_REGULAR_PARAMETERS* params, CMR_REGULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(params);

  CMRdbgMsg(0, "Testing %s %dx%d matrix for regularity.\n", ternary ? "ternary" : "binary", matrix->numRows,
    matrix->numColumns);
#if defined(CMR_DEBUG)
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
#endif /* CMR_DEBUG */

  clock_t time = clock();
  if (stats)
    stats->totalCount++;

  CMR_DEC* dec = NULL;
  CMR_CALL( CMRdecCreate(cmr, NULL, matrix->numRows, NULL, matrix->numColumns, NULL, &dec) );
  dec->matrix = matrix;
  assert(dec);

  CMR_CALL( CMRregularDecomposeOneSum(cmr, dec) );

  CMRdbgMsg(2, "1-sum decomposition yields %d components.\n", dec->numChildren == 0 ? 1 : dec->numChildren);
#if defined(CMR_DEBUG)
  CMR_CALL( CMRdecPrint(cmr, dec, stdout, 2, true, true, true) );
#endif /* CMR_DEBUG */

  bool isRegular = true;
  double remainingTime;
  if (dec->numChildren)
  {
    for (size_t c = 0; c < dec->numChildren; ++c)
    {
      if (isRegular || params->completeTree)
      {
        bool childIsRegular = true;
        remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
        CMR_CALL( testRegularTwoConnected(cmr, dec->children[c], ternary, &childIsRegular, pminor, params, stats,
          remainingTime) );
        if (pminor && *pminor)
          CMR_CALL( CMRdecTranslateMinorToParent(dec->children[c], *pminor) );

        isRegular = isRegular && childIsRegular;
      }
    }
  }
  else
  {
    remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( testRegularTwoConnected(cmr, dec, ternary, &isRegular, pminor, params, stats, remainingTime) );
  }

  if (pisRegular)
    *pisRegular = isRegular;

  CMR_CALL( CMRdecComputeRegularity(dec) );

  /* If the root must have a matrix, we eventually copy the one passed by the user. */
  if (pdec && (params->matrices == CMR_DEC_CONSTRUCT_ALL
    || (params->matrices == CMR_DEC_CONSTRUCT_LEAVES)) && dec->numChildren == 0)
  {
    dec->matrix = NULL;
    CMR_CALL( CMRchrmatCopy(cmr, matrix, &dec->matrix) );
  }
  else
    dec->matrix = NULL;

  /* Either store or free the decomposition. */
  if (pdec)
    *pdec = dec;
  else
    CMR_CALL( CMRdecFree(cmr, &dec) );

  if (stats)
  {
    stats->totalTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestBinaryRegular(CMR* cmr, CMR_CHRMAT* matrix, bool *pisRegular, CMR_DEC** pdec, CMR_MINOR** pminor,
  CMR_REGULAR_PARAMETERS* params, CMR_REGULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(pisRegular);

  CMR_REGULAR_PARAMETERS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRparamsRegularInit(&defaultParams) );
    params = &defaultParams;
  }

  CMR_SUBMAT* submatrix = NULL;
  if (!CMRchrmatIsBinary(cmr, matrix, pminor ? &submatrix : NULL))
  {
    *pisRegular = false;
    if (pminor)
      CMR_CALL( CMRminorCreate(cmr, pminor, 0, submatrix) );
    return CMR_OKAY;
  }

  CMR_CALL( CMRtestRegular(cmr, matrix, false, pisRegular, pdec, pminor, params, stats, timeLimit) );

  return CMR_OKAY;
}

