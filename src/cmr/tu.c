// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/tu.h>

#include "matrix_internal.h"
#include "block_decomposition.h"
#include "camion_internal.h"
#include "regularity_internal.h"
#include "hereditary_property.h"

#include <stdlib.h>
#include <assert.h>
#include <time.h>

CMR_ERROR CMRtuParamsInit(CMR_TU_PARAMS* params)
{
  assert(params);

  params->algorithm = CMR_TU_ALGORITHM_DECOMPOSITION;
  params->directCamion = false;
  CMR_CALL( CMRregularParamsInit(&params->regular) );

  return CMR_OKAY;
}

CMR_ERROR CMRtuStatsInit(CMR_TU_STATS* stats)
{
  assert(stats);

  return CMRregularStatsInit(stats);
}

CMR_ERROR CMRtuStatsPrint(FILE* stream, CMR_TU_STATS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  return CMRregularStatsPrint(stream, stats, prefix);
}

static
CMR_ERROR tuTest(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Some matrix to be tested for total unimodularity. */
  void* data,                 /**< Additional data (must be \c NULL). */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \p matrix is totally unimodular. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a proper non-totally unimodular submatrix of \p matrix. */
  double timeLimit            /**< Time limit to impose. */
)
{
  assert(cmr);
  assert(matrix);
  assert(pisTotallyUnimodular);
  assert(!psubmatrix || !*psubmatrix);

  CMR_UNUSED(psubmatrix); /* TODO: Make use of submatrices. */

  CMR_TU_STATS* stats = (CMR_TU_STATS*) data;

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "tuTest called for a %dx%d matrix\n", matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
  fflush(stdout);
#endif /* CMR_DEBUG */

  *pisTotallyUnimodular = true;
  clock_t time = clock();

  CMR_TU_PARAMS params; /* TODO: We should supply some params?! */
  CMR_CALL( CMRtuParamsInit(&params) );
  double remainingTime = timeLimit - ((clock() - time) * 1.0 / CLOCKS_PER_SEC);
  CMR_CALL( CMRtuTest(cmr, matrix, pisTotallyUnimodular, NULL, NULL, &params,
    stats ? stats : NULL, remainingTime) );

  return CMR_OKAY;
}


/**
 * \brief Recursively assigns +1 or -1 to each row that is part of the subset to test Ghouila-Houri.
 *
 * \return whether a feasible assignment was found.
 */

static
bool testPartitionSearch(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  int8_t* selection,          /**< Array with selection. */
  size_t current,             /**< Index to decide for selection. */
  int* columnSum              /**< Array for computing column sums. */
)
{
  assert(cmr);
  assert(matrix);
  assert(selection);

  while (current < matrix->numRows && selection[current] == 0)
    ++current;

  if (current < matrix->numRows)
  {
    /* Recurse by keeping current row a +1. */
    bool found = testPartitionSearch(cmr, matrix, selection, current + 1, columnSum);
    if (found)
      return true;

    /* Recurse by making current row a -1. */
    size_t first = matrix->rowSlice[current];
    size_t beyond = matrix->rowSlice[current + 1];

    selection[current] = -1;
    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] -= 2 * matrix->entryValues[i];

    found = testPartitionSearch(cmr, matrix, selection, current + 1, columnSum);

    selection[current] = +1;
    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] += 2 * matrix->entryValues[i];

    return found;
  }
  else
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
    {
      int sum = columnSum[column];
      if (sum < -1 || sum > +1)
        return false;
    }
    return true;
  }
}

/**
 * \brief Recursively selects a row subset and tests Ghouila-Houri for each.
 *
 * \return 1 if totally unimodular, 0 if not, and -1 if time limit was reached.
 */

static
int testPartitionSubset(
  CMR* cmr,           /**< \ref CMR environment */
  CMR_CHRMAT* matrix, /**< Matrix \f$ M \f$. */
  int8_t* selection,  /**< Array with selection. */
  size_t current,     /**< Index to decide for selection. */
  int* columnSum,     /**< Array for computing column sums. */
  clock_t startClock, /**< Clock for start for computation. */
  double timeLimit    /**< Time limit for computation. */
)
{
  assert(cmr);
  assert(matrix);
  assert(selection);

  if (current < matrix->numRows)
  {
    /* Recurse by not selecting a column. */
    selection[current] = 0;
    int result = testPartitionSubset(cmr, matrix, selection, current + 1, columnSum, startClock, timeLimit);
    if (result <= 0)
      return result;

    /* Recurse by selecting a column unless we need to abort. */
    selection[current] = 1;
    size_t first = matrix->rowSlice[current];
    size_t beyond = matrix->rowSlice[current + 1];
    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] += matrix->entryValues[i];

    result = testPartitionSubset(cmr, matrix, selection, current + 1, columnSum, startClock, timeLimit);

    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] -= matrix->entryValues[i];

    return result;
  }
  else
  {
    if (((clock() * 1.0 / CLOCKS_PER_SEC) - startClock) > timeLimit)
      return -1;

#if defined(CMR_DEBUG)
    size_t count = 0;
    for (size_t row = 0; row < matrix->numRows; ++row)
    {
      if (selection[row])
        ++count;
    }
    CMRdbgMsg(4, "-> Selected a row subset of cardinality %zu.\n", count);
#endif /* CMR_DEBUG */

    bool foundPartition = testPartitionSearch(cmr, matrix, selection, 0, columnSum);
    return foundPartition ? 1 : 0;
  }
}

/**
 * \brief Partition test based on Ghouila-Houri.
 */

static
CMR_ERROR testPartition(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  double timeLimit            /**< Time limit to impose. */
)
{
  assert(cmr);
  assert(matrix);
  assert(pisTotallyUnimodular);

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "testPartition called for a %zux%zu matrix\n", matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
  fflush(stdout);
#endif /* CMR_DEBUG */

  CMR_ERROR error = CMR_OKAY;

  /* Consider transpose if this has fewer rows. */
  if (matrix->numRows > matrix->numColumns)
  {
    CMR_CHRMAT* transpose = NULL;
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
    CMR_CALL( testPartition(cmr, transpose, pisTotallyUnimodular, timeLimit) );
    CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    return CMR_OKAY;
  }

  int8_t* selection = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &selection, matrix->numRows) );
  int* columnSum = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnSum, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnSum[column] = 0;

  clock_t startClock = clock();
  int result = testPartitionSubset(cmr, matrix, selection, 0, columnSum, startClock, timeLimit);
  if (result < 0)
    error = CMR_ERROR_TIMEOUT;
  else
    *pisTotallyUnimodular = (result > 0);

  CMR_CALL( CMRfreeStackArray(cmr, &columnSum) );
  CMR_CALL( CMRfreeStackArray(cmr, &selection) );

  return error;
}

CMR_ERROR CMRtuTest(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_MATROID_DEC** pdec,
  CMR_SUBMAT** psubmatrix, CMR_TU_PARAMS* params, CMR_TU_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );

  CMR_TU_PARAMS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRtuParamsInit(&defaultParams) );
    params = &defaultParams;
  }

  clock_t totalClock = clock();

  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  double remainingTime = timeLimit - ((clock() - totalClock) * 1.0 / CLOCKS_PER_SEC);

  CMRdbgMsg(0, "CMRtuTest called with algorithm = %d.\n", params->algorithm);

  if (params->algorithm == CMR_TU_ALGORITHM_DECOMPOSITION)
  {
    if (params->directCamion)
    {
      CMR_CALL( CMRcamionTestSigns(cmr, matrix, pisTotallyUnimodular, psubmatrix, stats ? &stats->camion : NULL,
        remainingTime) );

      if (!*pisTotallyUnimodular)
      {
        if (stats)
        {
          stats->totalCount++;
          stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
        }
        return CMR_OKAY;
      }
    }

    double remainingTime = timeLimit - ((clock() - totalClock) * 1.0 / CLOCKS_PER_SEC);

    CMR_CALL( CMRregularityTest(cmr, matrix, !params->directCamion, pisTotallyUnimodular, pdec, NULL, &params->regular,
      stats, remainingTime) );

    if (!*pisTotallyUnimodular && psubmatrix)
    {
      assert(!*psubmatrix);
      remainingTime = timeLimit - (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
      CMR_CALL( CMRtestHereditaryPropertySimple(cmr, matrix, tuTest, stats, psubmatrix, remainingTime) );
    }
  }
  else if (params->algorithm == CMR_TU_ALGORITHM_SUBMATRIX)
  {
    assert(!"Not implemented");
  }
  else if (params->algorithm == CMR_TU_ALGORITHM_PARTITION)
  {
    CMR_CALL( testPartition(cmr, matrix, pisTotallyUnimodular, remainingTime) );
  }
  else
  {
    return CMR_ERROR_INVALID;
  }

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}
