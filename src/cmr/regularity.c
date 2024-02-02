// #define CMR_DEBUG /** Uncomment to debug this file. */

#include "env_internal.h"
#include "regularity_internal.h"

#include <time.h>


CMR_ERROR CMRregularityTaskCreateRoot(CMR* cmr, CMR_MATROID_DEC* dec, DecompositionTask** ptask,
  CMR_REGULAR_PARAMS* params, CMR_REGULAR_STATS* stats, clock_t startClock, double timeLimit)
{
  assert(cmr);
  assert(dec);
  assert(ptask);
  assert(params);

  CMR_CALL( CMRallocBlock(cmr, ptask) );
  DecompositionTask* task = *ptask;

  task->dec = dec;
  task->next = NULL;

  task->params = params;
  task->stats = stats;
  task->startClock = startClock;
  task->timeLimit = timeLimit;

  return CMR_OKAY;
}

CMR_ERROR CMRregularityTaskFree(CMR* cmr, DecompositionTask** ptask)
{
  assert(cmr);
  assert(ptask);

  if (*ptask == NULL)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlock(cmr, ptask) );

  return CMR_OKAY;
}

// /**
//  * \brief Tests a 2-connected binary or ternary matrix for regularity.
//  */
//
// static
// CMR_ERROR testRegularTwoConnected(
//   CMR* cmr,                   /**< \ref CMR environment. */
//   CMR_MATROID_DEC* dec,       /**< Decomposition node. */
//   bool ternary,               /**< Whether signs matter. */
//   bool *pisRegular,           /**< Pointer for storing whether \p matrix is regular. */
//   CMR_MINOR** pminor,         /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
//   CMR_REGULAR_PARAMS* params, /**< Parameters for the computation. */
//   CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
//   double timeLimit            /**< Time limit to impose. */
// );

// static
// CMR_ERROR testRegularThreeConnectedWithSequence(
//   CMR* cmr,                   /**< \ref CMR environment. */
//   CMR_MATROID_DEC* dec,       /**< Decomposition node. */
//   bool ternary,               /**< Whether signs matter. */
//   bool *pisRegular,           /**< Pointer for storing whether \p matrix is regular. */
//   CMR_MINOR** pminor,         /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
//   CMR_REGULAR_PARAMS* params, /**< Parameters for the computation. */
//   CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
//   double timeLimit            /**< Time limit to impose. */
// )
// {
//   assert(cmr);
//   assert(dec);
//   assert(dec->matrix);
// //   assert(dec->nestedMinorsMatrix);
// //   assert(dec->nestedMinorsRowsOriginal);
// //   assert(dec->nestedMinorsColumnsOriginal);
// //   assert(dec->nestedMinorsSequenceNumRows);
// //   assert(dec->nestedMinorsSequenceNumColumns);
// //   assert(dec->nestedMinorsLength > 0);
//
//   clock_t time = clock();
//
//   CMRdbgMsg(6, "Testing binary %dx%d 3-connected matrix with given nested sequence of 3-connected minors for regularity.\n",
//     dec->matrix->numRows, dec->matrix->numColumns);
//
// #if defined(CMR_DEBUG)
// //   CMR_CALL( CMRdecPrintSequenceNested3ConnectedMinors(cmr, dec, stdout) );
// #endif /* CMR_DEBUG */
//
//   CMRconsistencyAssert( CMRmatroiddecConsistency(dec, false) );
//
//   CMR_CHRMAT* nestedMinorsTranspose = NULL;
// //   CMR_CALL( CMRchrmatTranspose(cmr, dec->nestedMinorsMatrix, &nestedMinorsTranspose) );
//
//   /* Test sequence for graphicness. */
//   double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//   size_t lastGraphicMinor = 0;
//   CMR_GRAPH* graph = NULL;
//   CMR_ELEMENT* graphEdgeLabels = NULL;
// //   CMR_CALL( CMRregularSequenceGraphic(cmr, dec->nestedMinorsMatrix, nestedMinorsTranspose,
// //     dec->nestedMinorsLength, dec->nestedMinorsSequenceNumRows, dec->nestedMinorsSequenceNumColumns,
// //     &lastGraphicMinor, &graph, &graphEdgeLabels, stats, remainingTime) );
//
//   if (graph)
//   {
//     assert(graphEdgeLabels);
//
//     CMRdbgMsg(8, "Matrix is graphic.\n");
// //     dec->type = CMR_DEC_GRAPHIC;
//     dec->graph = graph;
//     dec->graphForest = NULL;
//     CMR_CALL( CMRallocBlockArray(cmr, &dec->graphForest, dec->matrix->numRows) );
//     dec->graphCoforest = NULL;
//     CMR_CALL( CMRallocBlockArray(cmr, &dec->graphCoforest, dec->matrix->numColumns) );
//
//     assert(CMRgraphNumEdges(graph) == dec->matrix->numRows + dec->matrix->numColumns);
//     for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
//       iter = CMRgraphEdgesNext(graph, iter))
//     {
//       CMR_GRAPH_EDGE edge = CMRgraphEdgesEdge(graph, iter);
//       CMR_ELEMENT element = graphEdgeLabels[edge];
// //       if (CMRelementIsRow(element))
// //         element = dec->nestedMinorsRowsOriginal[CMRelementToRowIndex(element)];
// //       else
// //         element = dec->nestedMinorsColumnsOriginal[CMRelementToColumnIndex(element)];
//       if (CMRelementIsRow(element))
//         dec->graphForest[CMRelementToRowIndex(element)] = edge;
//       else
//         dec->graphCoforest[CMRelementToColumnIndex(element)] = edge;
//     }
//
//     CMR_CALL( CMRfreeBlockArray(cmr, &graphEdgeLabels) );
//   }
//   else
//   {
// //     dec->flags &= ~CMR_DEC_IS_GRAPHIC;
//   }
//
//   size_t lastCographicMinor = 0;
//   if (!dec->graph || params->planarityCheck)
//   {
//     /* Test sequence for cographicness. */
//     CMR_GRAPH* cograph = NULL;
//     CMR_ELEMENT* cographEdgeLabels = NULL;
//     remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
// //     CMR_CALL( CMRregularSequenceGraphic(cmr, nestedMinorsTranspose, dec->nestedMinorsMatrix,
// //       dec->nestedMinorsLength, dec->nestedMinorsSequenceNumColumns, dec->nestedMinorsSequenceNumRows,
// //       &lastCographicMinor, &cograph, &cographEdgeLabels, stats, remainingTime) );
//
//     if (cograph)
//     {
//       assert(cographEdgeLabels);
//
//       CMRdbgMsg(8, "Matrix is graphic.\n");
// //       dec->type = (dec->type == CMR_DEC_GRAPHIC) ? CMR_DEC_PLANAR : CMR_DEC_COGRAPHIC;
//       dec->cograph = cograph;
//       dec->cographForest = NULL;
//       CMR_CALL( CMRallocBlockArray(cmr, &dec->cographForest, dec->matrix->numColumns) );
//       dec->cographCoforest = NULL;
//       CMR_CALL( CMRallocBlockArray(cmr, &dec->cographCoforest, dec->matrix->numRows) );
//
//       assert(CMRgraphNumEdges(cograph) == dec->matrix->numRows + dec->matrix->numColumns);
//       for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(cograph); CMRgraphEdgesValid(cograph, iter);
//         iter = CMRgraphEdgesNext(cograph, iter))
//       {
//         CMR_GRAPH_EDGE edge = CMRgraphEdgesEdge(cograph, iter);
//         CMR_ELEMENT element = cographEdgeLabels[edge];
// //         if (CMRelementIsRow(element))
// //           element = dec->nestedMinorsColumnsOriginal[CMRelementToRowIndex(element)];
// //         else
// //           element = dec->nestedMinorsRowsOriginal[CMRelementToColumnIndex(element)];
//         if (CMRelementIsRow(element))
//           dec->cographCoforest[CMRelementToRowIndex(element)] = edge;
//         else
//           dec->cographForest[CMRelementToColumnIndex(element)] = edge;
//       }
//
//       CMR_CALL( CMRfreeBlockArray(cmr, &cographEdgeLabels) );
//     }
//     else
//     {
// //       dec->flags &= ~CMR_DEC_IS_COGRAPHIC;
//     }
//   }
//
//   if (!dec->graph && !dec->cograph)
//   {
//     CMRdbgMsg(8, "Checking for R10.\n");
//     bool isR10;
//     CMR_CALL( CMRregularThreeConnectedIsR10(cmr, dec, &isR10) );
//
//     if (!isR10)
//     {
//       CMRdbgMsg(8, "Neither graphic (until %ld) nor cographic (until %ld), nor R_10! First non-graphic and non-cographic minor is %ld.\n",
//         lastGraphicMinor, lastCographicMinor,
//         lastGraphicMinor > lastCographicMinor ? (lastGraphicMinor+1) : (lastCographicMinor + 1) );
//
//       remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//       CMR_CALL( CMRregularSearchThreeSeparation(cmr, dec, nestedMinorsTranspose, ternary,
//         lastGraphicMinor > lastCographicMinor ? (lastGraphicMinor+1) : (lastCographicMinor + 1), NULL, stats,
//         remainingTime) );
//
//       if (dec->type == CMR_MATROID_DEC_TYPE_IRREGULAR)
//       {
//         CMRdbgMsg(8, "Minor determined to be irregular.\n");
//         *pisRegular = false;
//       }
//       else
//       {
//         CMRdbgMsg(8, " Encountered a 3-separation.\n");
//         assert(dec->numChildren == 2);
//
// #if defined(CMR_DEBUG)
//         CMR_CALL( CMRchrmatPrintDense(cmr, dec->children[0]->matrix, stdout, '0', true) );
//         CMR_CALL( CMRchrmatPrintDense(cmr, dec->children[1]->matrix, stdout, '0', true) );
// #endif /* CMR_DEBUG */
//
//         remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//         CMR_CALL( testRegularTwoConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params, stats,
//           remainingTime) );
//
//         if (params->completeTree || *pisRegular)
//         {
//           remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//           CMR_CALL( testRegularTwoConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params, stats,
//             remainingTime) );
//         }
//       }
//     }
//   }
//
//   CMR_CALL( CMRchrmatFree(cmr, &nestedMinorsTranspose) );
//
//   return CMR_OKAY;
// }
//
// static
// CMR_ERROR testRegularTwoConnected(CMR* cmr, CMR_MATROID_DEC* dec, bool ternary, bool *pisRegular, CMR_MINOR** pminor,
//   CMR_REGULAR_PARAMS* params, CMR_REGULAR_STATS* stats, double timeLimit)
// {
//   assert(cmr);
//   assert(dec);
//   assert(dec->matrix);
//
//   CMRdbgMsg(2, "Testing binary %dx%d 2-connected matrix for regularity.\n", dec->matrix->numRows,
//     dec->matrix->numColumns);
//
//   clock_t time = clock();
//   CMR_SUBMAT* submatrix = NULL;
//
//   if (params->directGraphicness || dec->matrix->numRows <= 3 || dec->matrix->numColumns <= 3)
//   {
//     /* We run the almost-linear time algorithm. Otherwise, graphicness is checked later for the 3-connected components. */
//
//     if (params->directGraphicness || params->planarityCheck || dec->matrix->numRows > 3)
//     {
//       CMRdbgMsg(4, "Checking for graphicness...");
//       bool isGraphic;
//       CMR_CALL( CMRregularTestGraphic(cmr, &dec->matrix, &dec->transpose, ternary, &isGraphic, &dec->graph,
//         &dec->graphForest, &dec->graphCoforest, &dec->graphArcsReversed, &submatrix, stats, timeLimit) );
//       if (isGraphic)
//       {
//         CMRdbgMsg(0, " graphic.\n");
// //         dec->type = CMR_MATROID_DEC_TYPE_GRAPHIC;
//         if (!params->planarityCheck)
//           return CMR_OKAY;
//       }
//       CMRdbgMsg(0, " NOT graphic.\n");
//     }
//
//     CMRdbgMsg(4, "Checking for cographicness...");
//     bool isCographic;
//     double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//     CMR_CALL( CMRregularTestGraphic(cmr, &dec->transpose, &dec->matrix, ternary, &isCographic, &dec->cograph,
//       &dec->cographForest, &dec->cographCoforest, &dec->cographArcsReversed, &submatrix, stats, remainingTime) );
//     if (isCographic)
//     {
//       CMRdbgMsg(0, " cographic.\n");
// //       dec->type = (dec->type == CMR_DEC_GRAPHIC) ? CMR_DEC_PLANAR : CMR_DEC_COGRAPHIC;
//       return CMR_OKAY;
//     }
//     CMRdbgMsg(0, " NOT cographic.\n");
//
//     if (submatrix)
//       CMR_CALL( CMRsubmatTranspose(submatrix) );
//   }
//
//   if (dec->nestedMinorsMatrix)
//   {
//     CMRdbgMsg(4, "A partial sequence of nested minors is already known.\n");
//
// #if defined(CMR_DEBUG)
//     CMR_CALL( CMRdecPrintSequenceNested3ConnectedMinors(cmr, dec, stdout) );
// #endif /* CMR_DEBUG */
//
//     double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//     CMR_CALL( CMRregularExtendNestedMinorSequence(cmr, dec, ternary, &submatrix, stats, remainingTime) );
//
//     /* Handling of the resulting sequence or 2-separation is done at the end. */
//   }
//   else
//   {
//     CMRdbgMsg(4, "Splitting off series-parallel elements...");
//     double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//     CMR_CALL( CMRregularDecomposeSeriesParallel(cmr, &dec, ternary, &submatrix, params, stats, remainingTime) );
//
//     if (dec->type == CMR_MATROID_DEC_TYPE_IRREGULAR)
//     {
//       CMRdbgMsg(0, " NOT regular.\n");
//       return CMR_OKAY;
//     }
//     else if (dec->type == CMR_MATROID_DEC_TYPE_SERIES_PARALLEL)
//     {
//       CMRdbgMsg(0, " series-parallel.\n");
//       return CMR_OKAY;
//     }
//
//     if (dec->type == CMR_MATROID_DEC_TYPE_TWO_SUM)
//     {
//       CMRdbgMsg(0, " Encountered a 2-separation.\n");
//       assert(dec->numChildren == 2);
//       remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//       CMR_CALL( testRegularTwoConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params, stats,
//         remainingTime) );
//
//       if (params->completeTree || *pisRegular)
//       {
//         remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//         CMR_CALL( testRegularTwoConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params, stats,
//           remainingTime) );
//       }
//
//       return CMR_OKAY;
//     }
//
//     CMR_SUBMAT* wheelSubmatrix = submatrix;
//     CMRdbgMsg(0, " Found a W_%d minor.\n", wheelSubmatrix->numRows);
//     submatrix = NULL;
//
//     /* No 2-sum found, so we have a wheel submatrix. */
//
//     remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//     CMR_CALL( CMRregularConstructNestedMinorSequence(cmr, dec, ternary, wheelSubmatrix, &submatrix, stats,
//       remainingTime) );
//     CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
//   }
//
//   if (dec->type == CMR_MATROID_DEC_TYPE_IRREGULAR)
//   {
//     CMRdbgMsg(0, " NOT regular.\n");
//     return CMR_OKAY;
//   }
//
//   if (dec->type == CMR_MATROID_DEC_TYPE_TWO_SUM)
//   {
//     CMRdbgMsg(0, " Encountered a 2-separation.\n");
//     assert(dec->numChildren == 2);
//
//     double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//     CMR_CALL( testRegularTwoConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params, stats,
//       remainingTime) );
//
//     if (params->completeTree || *pisRegular)
//     {
//       remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//       CMR_CALL( testRegularTwoConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params, stats,
//         remainingTime) );
//     }
//
//     return CMR_OKAY;
//   }
//
//   double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
//   CMR_CALL( testRegularThreeConnectedWithSequence(cmr, dec, ternary, pisRegular, pminor, params, stats,
//     remainingTime) );
//
//   return CMR_OKAY;
// }

/**
 * \brief Runs a task for processing the associated decomposition node.
 */

static
CMR_ERROR CMRregularityTaskRun(
  CMR* cmr,                         /**< \ref CMR environment. */
  DecompositionTask* task,          /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionTask** punprocessed  /**< Pointer to head of list of unprocessed tasks. */
)
{
  assert(cmr);
  assert(task);
  assert(punprocessed);

  CMRdbgMsg(2, "Processing %p.\n", task);

  if (!task->dec->testedTwoConnected)
  {
    CMRdbgMsg(4, "Searching for 1-separations.\n");
    CMR_CALL( CMRregularitySearchOneSum(cmr, task, punprocessed) );
  }
  else if (!task->dec->graphicness
    && (task->params->directGraphicness || task->dec->matrix->numRows <= 3 || task->dec->matrix->numColumns <= 3))
  {
    CMRdbgMsg(4, "Testing directly for %s.\n", task->dec->isTernary ? "being network" : "graphicness");
    CMR_CALL( CMRregularityTestGraphicness(cmr, task, punprocessed) );
  }
  else if (!task->dec->cographicness
    && (task->params->directGraphicness || task->dec->matrix->numRows <= 3 || task->dec->matrix->numColumns <= 3))
  {
    CMRdbgMsg(4, "Testing directly for %s.\n", task->dec->isTernary ? "being conetwork" : "cographicness");
    CMR_CALL( CMRregularityTestCographicness(cmr, task, punprocessed) );
  }
  else if (!task->dec->testedR10)
  {
    CMRdbgMsg(4, "Testing for being R_10.\n");
    CMR_CALL( CMRregularityTestR10(cmr, task, punprocessed) );
  }
  else if (!task->dec->testedSeriesParallel)
  {
    CMRdbgMsg(4, "Testing for series-parallel reductions.\n");
    CMR_CALL( CMRregularityDecomposeSeriesParallel(cmr, task, punprocessed) );
  }
  else if (task->dec->denseMatrix)
  {
    CMRdbgMsg(4, "Attempting to construct a sequence of nested minors.\n");
    CMR_CALL( CMRregularityExtendNestedMinorSequence(cmr, task, punprocessed) );
  }
  else if (task->dec->nestedMinorsMatrix && (task->dec->nestedMinorsLastGraphic == SIZE_MAX))
  {
    CMRdbgMsg(4, "Testing along the sequence for %s.\n", task->dec->isTernary ? "being network" : "graphicness");
    CMR_CALL( CMRregularityNestedMinorSequenceGraphicness(cmr, task, punprocessed) );
  }
  else if (task->dec->nestedMinorsMatrix && (task->dec->nestedMinorsLastCographic == SIZE_MAX))
  {
    CMRdbgMsg(4, "Testing for %s along the sequence.\n", task->dec->isTernary ? "being conetwork" : "cographicness");
    CMR_CALL( CMRregularityNestedMinorSequenceCographicness(cmr, task, punprocessed) );
  }
  else
  {
    CMRdbgMsg(4, "Searching for 3-separations along the sequence.\n");
    CMR_CALL( CMRregularityNestedMinorSequenceSearchThreeSeparation(cmr, task, punprocessed) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRregularityTest(CMR* cmr, CMR_CHRMAT* matrix, bool ternary, bool *pisRegular, CMR_MATROID_DEC** pdec,
  CMR_MINOR** pminor, CMR_REGULAR_PARAMS* params, CMR_REGULAR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(params);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "Testing %s %dx%d matrix for regularity.\n", ternary ? "ternary" : "binary", matrix->numRows,
    matrix->numColumns);
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
#endif /* CMR_DEBUG */

  clock_t time = clock();
  if (stats)
    stats->totalCount++;

  CMR_MATROID_DEC* root = NULL;
  CMR_CALL( CMRmatroiddecCreateMatrixRoot(cmr, &root, ternary, matrix) );
  assert(root);

  DecompositionTask* headUnprocessedTasks = NULL;
  CMR_CALL( CMRregularityTaskCreateRoot(cmr, root, &headUnprocessedTasks, params, stats, time, timeLimit) );
  while (headUnprocessedTasks)
  {
    DecompositionTask* task = headUnprocessedTasks;
    headUnprocessedTasks = task->next;

    CMR_CALL( CMRregularityTaskRun(cmr, task, &headUnprocessedTasks) );
  }

  CMR_CALL( CMRmatroiddecSetAttributes(root) );
  assert(root->regularity != 0);
  if (pisRegular)
    *pisRegular = root->regularity > 0;

  /* Either store or free the decomposition. */
  if (pdec)
    *pdec = root;
  else
    CMR_CALL( CMRmatroiddecFree(cmr, &root) );

  if (stats)
    stats->totalTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;

  return CMR_OKAY;
}


