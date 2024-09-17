// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/tu.h>

#include "matrix_internal.h"
#include "block_decomposition.h"
#include "camion_internal.h"
#include "hereditary_property.h"
#include "seymour_internal.h"

#include <stdlib.h>
#include <assert.h>
#include <time.h>

CMR_ERROR CMRtuParamsInit(CMR_TU_PARAMS* params)
{
  assert(params);

  params->algorithm = CMR_TU_ALGORITHM_DECOMPOSITION;
  params->ternary = true;
  params->camionFirst = true;
  params->naiveSubmatrix = false;
  CMR_CALL( CMRseymourParamsInit(&params->seymour) );

  return CMR_OKAY;
}

CMR_ERROR CMRtuStatsInit(CMR_TU_STATS* stats)
{
  assert(stats);

  CMR_CALL( CMRseymourStatsInit(&stats->seymour) );
  CMR_CALL( CMRcamionStatsInit(&stats->camion) );

  stats->enumerationRowSubsets = 0;
  stats->enumerationColumnSubsets = 0;
  stats->enumerationTime = 0.0;

  stats->partitionRowSubsets = 0;
  stats->partitionColumnSubsets = 0;
  stats->partitionTime = 0.0;

  return CMR_OKAY;
}

CMR_ERROR CMRtuStatsPrint(FILE* stream, CMR_TU_STATS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "TU matrix recognition:\n");
    prefix = "  ";
  }

  char subPrefix[256];
  snprintf(subPrefix, 256, "%sseymour ", prefix);
  CMR_CALL( CMRseymourStatsPrint(stream, &stats->seymour, subPrefix) );
  snprintf(subPrefix, 256, "%scamion ", prefix);
  CMR_CALL( CMRcamionStatsPrint(stream, &stats->camion, subPrefix) );

  fprintf(stream, "%seulerian enumeration row subsets: %lu\n", prefix, (unsigned long)stats->enumerationRowSubsets);
  fprintf(stream, "%seulerian enumeration column subsets: %lu\n", prefix, (unsigned long)stats->enumerationColumnSubsets);
  fprintf(stream, "%seulerian enumeration time: %f\n", prefix, stats->enumerationTime);

  fprintf(stream, "%spartition row subsets: %lu\n", prefix, (unsigned long)stats->partitionRowSubsets);
  fprintf(stream, "%spartition column subsets: %lu\n", prefix, (unsigned long)stats->partitionColumnSubsets);
  fprintf(stream, "%spartition time: %f\n", prefix, stats->partitionTime);

  return CMR_OKAY;
}

static
CMR_ERROR tuDecomposition(
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
  CMRdbgMsg(0, "tuDecomposition called for a %dx%d matrix\n", matrix->numRows, matrix->numColumns);
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
 * \brief Data for enumeration.
 */

typedef struct
{
  CMR* cmr;                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix;           /**< Matrix \f$ M \f$. */
  bool* pisTotallyUnimodular;   /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  CMR_SUBMAT** psubmatrix;      /**< Pointer for storing a minimal nonbalanced submatrix (may be \c NULL). */
  CMR_TU_STATS* stats;          /**< Statistics for the computation (may be \c NULL). */
  double timeLimit;             /**< Time limit to impose. */
  bool isTransposed;            /**< Whether we're dealing with the transposed matrix. */
  clock_t startClock;           /**< Clock for when we started. */
  size_t cardinality;           /**< Current cardinality of row/column subsets. */
  size_t* subsetRows;           /**< Array for the enumerated row subset. */
  size_t* usableColumns;        /**< Array of columns usable for enumeration. */
  size_t numUsableColumns;      /**< Length of usableColumns. */
  size_t* subsetUsable;         /**< Array for the enumerated subset from usable columns. */
  size_t* columnsNumNonzeros;   /**< Array with the number of nonzeros per column. */
  size_t* rowsNumNonzeros;      /**< Array with the number of nonzeros per row. */
  int sumEntries;               /**< Sum of the entries of the selected submatrix. */
} CMR_TU_ENUMERATION;


/**
 * \brief Recursive enumeration of column subsets and subsequent testing of Eulerian submatrix.
 */

static
CMR_ERROR tuEulerianColumns(
  CMR_TU_ENUMERATION* enumeration,  /**< Enumeration information. */
  size_t numColumns                 /**< Number of already selected columns. */
)
{
  assert(enumeration);

  if (numColumns < enumeration->cardinality)
  {
    /* Recursion step: pick one column and enumerate over the other columns. */

    size_t firstUsable = (numColumns == 0) ? 0 : enumeration->subsetUsable[numColumns - 1] + 1;
    size_t beyondUsable = enumeration->numUsableColumns - enumeration->cardinality + numColumns + 1;
    for (size_t usable = firstUsable; usable < beyondUsable; ++usable)
    {
      size_t column = enumeration->usableColumns[usable];
      assert(enumeration->columnsNumNonzeros[column] % 2 == 0);

      CMRdbgMsg(12 + enumeration->cardinality + numColumns,
        "Selecting column %zu as usable %zu of %zu.\n", column, usable, enumeration->numUsableColumns);
      enumeration->subsetUsable[numColumns] = usable;

      /* Increment row nonzero counters. */
      for (size_t i = 0; i < enumeration->cardinality; ++i)
      {
        size_t row = enumeration->subsetRows[i];
        size_t entry;
        CMRchrmatFindEntry(enumeration->matrix, row, column, &entry);
        if (entry != SIZE_MAX)
        {
          enumeration->sumEntries += enumeration->matrix->entryValues[entry];
          enumeration->rowsNumNonzeros[row]++;
        }
      }

      /* Recurse. */
      CMR_CALL( tuEulerianColumns(enumeration, numColumns + 1) );
      if (!*enumeration->pisTotallyUnimodular)
        return CMR_OKAY;

      /* Decrement row nonzero counters. */
      for (size_t i = 0; i < enumeration->cardinality; ++i)
      {
        size_t row = enumeration->subsetRows[i];
        size_t entry;
        CMRchrmatFindEntry(enumeration->matrix, row, column, &entry);
        if (entry != SIZE_MAX)
        {
          enumeration->sumEntries -= enumeration->matrix->entryValues[entry];
          enumeration->rowsNumNonzeros[row]--;
        }
      }
    }
  }
  else
  {
    if (enumeration->stats)
    {
      if (enumeration->isTransposed)
        enumeration->stats->enumerationRowSubsets++;
      else
        enumeration->stats->enumerationColumnSubsets++;
    }

    CMRdbgMsg(13 + 2 * enumeration->cardinality, "Sum over nonzeros is %d.\n", enumeration->sumEntries);

    if (enumeration->sumEntries % 4 != 0)
    {
      /* Now check if we're Eulerian. */
      bool isEulerian = true;
      for (size_t i = 0; i < enumeration->cardinality; ++i)
      {
        size_t row = enumeration->subsetRows[i];
        if (enumeration->rowsNumNonzeros[row] % 2 != 0)
        {
          isEulerian = false;
          break;
        }
      }

      CMRdbgMsg(14 + 2 * enumeration->cardinality, "Eulerian = %d.\n", isEulerian);

      if (isEulerian)
      {
        *(enumeration->pisTotallyUnimodular) = false;
        if (enumeration->psubmatrix)
        {
          CMR_SUBMAT* submatrix = NULL;
          CMR_CALL( CMRsubmatCreate(enumeration->cmr, enumeration->cardinality, enumeration->cardinality, &submatrix) );
          for (size_t i = 0; i < enumeration->cardinality; ++i)
          {
            submatrix->rows[i] = enumeration->subsetRows[i];
            submatrix->columns[i] = enumeration->usableColumns[enumeration->subsetUsable[i]];
          }

          *(enumeration->psubmatrix) = submatrix;
        }
      }
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Recursive enumeration of row subsets and subsequent testing of total unimodularity.
 *
 * \returns \c true if we can stop (due to time limit or a found violating submatrix).
 */

static
bool tuEulerianRows(
  CMR_TU_ENUMERATION* enumeration,  /**< Enumeration information. */
  size_t numRows                    /**< Number of already selected rows. */
)
{
  assert(enumeration);

  if (numRows < enumeration->cardinality)
  {
    /* Recursion step: pick one row and enumerate over the other rows. */

    size_t firstRow = (numRows == 0) ? 0 : enumeration->subsetRows[numRows - 1] + 1;
    size_t beyondRow = enumeration->matrix->numRows - enumeration->cardinality + numRows + 1;
    for (size_t row = firstRow; row < beyondRow; ++row)
    {
      CMRdbgMsg(10 + numRows, "Selecting row %zu as row #%zu.\n", row, numRows);
      enumeration->subsetRows[numRows] = row;

      /* Iterate over chosen row and increment column nonzero counters. */
      size_t first = enumeration->matrix->rowSlice[row];
      size_t beyond = enumeration->matrix->rowSlice[row + 1];
      for (size_t e = first; e < beyond; ++e)
        enumeration->columnsNumNonzeros[enumeration->matrix->entryColumns[e]]++;

      /* Recurse. */
      CMR_CALL( tuEulerianRows(enumeration, numRows + 1) );
      if (!*(enumeration->pisTotallyUnimodular) || enumeration->timeLimit <= 0)
        return CMR_OKAY;

      /* Decrement column nonzero counters again. */
      for (size_t e = first; e < beyond; ++e)
        enumeration->columnsNumNonzeros[enumeration->matrix->entryColumns[e]]--;
    }
  }
  else
  {
    if (enumeration->stats)
    {
      if (enumeration->isTransposed)
        enumeration->stats->enumerationColumnSubsets++;
      else
        enumeration->stats->enumerationRowSubsets++;
    }

    double remainingTime = enumeration->timeLimit - (clock() - enumeration->startClock) * 1.0 / CLOCKS_PER_SEC;
    if (remainingTime <= 0)
    {
      enumeration->timeLimit = 0;
      return CMR_OKAY;
    }

    assert(enumeration->sumEntries == 0);
    for (size_t i = 0; i < enumeration->cardinality; ++i)
      assert( enumeration->rowsNumNonzeros[enumeration->subsetRows[i]] == 0 );

    enumeration->numUsableColumns = 0;
    for (size_t column = 0; column < enumeration->matrix->numColumns; ++column)
    {
      if (enumeration->columnsNumNonzeros[column] % 2 == 0)
        enumeration->usableColumns[enumeration->numUsableColumns++] = column;
    }

    CMRdbgMsg(10 + numRows, "Number of columns with an even number of nonzeros: %zu.\n", enumeration->numUsableColumns);

    /* Skip enumeraton of columns if the number of columns with an even number of nonzeros is too low. */
    if (enumeration->numUsableColumns >= enumeration->cardinality)
    {
      CMR_CALL( tuEulerianColumns(enumeration, 0) );
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Submatrix test.
 */

static
CMR_ERROR tuEulerian(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool isTransposed,          /**< Whether we're dealing with the transpose. */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing the submatrix. */
  CMR_TU_STATS* stats,        /**< Statistics. */
  double timeLimit            /**< Time limit to impose. */
)
{
  assert(cmr);
  assert(matrix);
  assert(pisTotallyUnimodular);

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "tuEulerian called for a %zux%zu matrix\n", matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
  fflush(stdout);
#endif /* CMR_DEBUG */

  CMRassertStackConsistency(cmr);

  CMR_ERROR error = CMR_OKAY;

  /* Transpose if more rows than columns. */
  if (matrix->numRows > matrix->numColumns)
  {
    assert(!isTransposed);

    CMR_CHRMAT* transpose = NULL;
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMRdbgMsg(6, "Transposing matrix to have fewer rows than columns.\n");

    error = tuEulerian(cmr, transpose, true, pisTotallyUnimodular, psubmatrix, stats, timeLimit);

    /* Transpose the violator. */
    if (psubmatrix && *psubmatrix)
      CMR_CALL( CMRsubmatTranspose(*psubmatrix) );

    CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    return error;
  }

  CMR_TU_ENUMERATION enumeration;
  enumeration.cmr = cmr;
  enumeration.matrix = matrix;
  enumeration.pisTotallyUnimodular = pisTotallyUnimodular;
  enumeration.psubmatrix = psubmatrix;
  enumeration.stats = stats;
  enumeration.timeLimit = timeLimit;
  enumeration.isTransposed = isTransposed;
  enumeration.startClock = clock();
  enumeration.subsetRows = NULL;
  enumeration.usableColumns = NULL;
  enumeration.subsetUsable = NULL;
  enumeration.rowsNumNonzeros = NULL;
  enumeration.columnsNumNonzeros = NULL;
  enumeration.sumEntries = 0;
  *pisTotallyUnimodular = true;

  CMR_CALL( CMRallocStackArray(cmr, &enumeration.subsetRows, matrix->numRows) );
  CMR_CALL( CMRallocStackArray(cmr, &enumeration.usableColumns, matrix->numColumns) );
  CMR_CALL( CMRallocStackArray(cmr, &enumeration.subsetUsable, matrix->numColumns) );
  CMR_CALL( CMRallocStackArray(cmr, &enumeration.rowsNumNonzeros, matrix->numRows) );
  CMR_CALL( CMRallocStackArray(cmr, &enumeration.columnsNumNonzeros, matrix->numColumns) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    enumeration.rowsNumNonzeros[row] = 0;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    enumeration.columnsNumNonzeros[column] = 0;

  CMRdbgMsg(6, "Starting enumeration algorithm with a time limit of %g.\n", timeLimit);

  CMRassertStackConsistency(cmr);

  for (enumeration.cardinality = 2; enumeration.cardinality <= matrix->numRows; enumeration.cardinality++)
  {
    CMRassertStackConsistency(cmr);
    CMRdbgMsg(8, "Considering submatrices of size %zux%zu.\n", enumeration.cardinality, enumeration.cardinality);
    CMRassertStackConsistency(cmr);
    CMR_CALL( tuEulerianRows(&enumeration, 0) );
    if (stats)
      stats->enumerationTime += (clock() - enumeration.startClock) * 1.0 / CLOCKS_PER_SEC;
    if (!*pisTotallyUnimodular || enumeration.timeLimit <= 0)
      break;
  }

  CMRassertStackConsistency(cmr);

  CMR_CALL( CMRfreeStackArray(cmr, &enumeration.columnsNumNonzeros) );
  CMR_CALL( CMRfreeStackArray(cmr, &enumeration.rowsNumNonzeros) );
  CMR_CALL( CMRfreeStackArray(cmr, &enumeration.subsetUsable) );
  CMR_CALL( CMRfreeStackArray(cmr, &enumeration.usableColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &enumeration.subsetRows) );

  CMRassertStackConsistency(cmr);

  if (enumeration.timeLimit <= 0)
    return CMR_ERROR_TIMEOUT;

  return CMR_OKAY;
}

/**
 * \brief Recursively assigns +1 or -1 to each row that is part of the subset to test Ghouila-Houri.
 *
 * \return whether a feasible assignment was found.
 */

static
bool tuPartitionSearch(
  CMR* cmr,           /**< \ref CMR environment */
  CMR_CHRMAT* matrix, /**< Matrix \f$ M \f$. */
  int8_t* selection,  /**< Array with selection. */
  size_t current,     /**< Index to decide for selection. */
  int* columnSum      /**< Array for computing column sums. */
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
    bool found = tuPartitionSearch(cmr, matrix, selection, current + 1, columnSum);
    if (found)
      return true;

    /* Recurse by making current row a -1. */
    size_t first = matrix->rowSlice[current];
    size_t beyond = matrix->rowSlice[current + 1];

    selection[current] = -1;
    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] -= 2 * matrix->entryValues[i];

    found = tuPartitionSearch(cmr, matrix, selection, current + 1, columnSum);

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
int tuPartitionSubset(
  CMR* cmr,             /**< \ref CMR environment */
  CMR_CHRMAT* matrix,   /**< Matrix \f$ M \f$. */
  bool transposed,      /**< Whether we're dealing with the transpose. */
  int8_t* selection,    /**< Array with selection. */
  size_t current,       /**< Index to decide for selection. */
  int* columnSum,       /**< Array for computing column sums. */
  CMR_TU_STATS* stats,  /**< Statistics. */
  clock_t startClock,   /**< Clock for start for computation. */
  double timeLimit      /**< Time limit for computation. */
)
{
  assert(cmr);
  assert(matrix);
  assert(selection);

  if (current < matrix->numRows)
  {
    /* Recurse by not selecting a column. */
    selection[current] = 0;
    int result = tuPartitionSubset(cmr, matrix, transposed, selection, current + 1, columnSum, stats, startClock, timeLimit);
    if (result <= 0)
      return result;

    /* Recurse by selecting a column unless we need to abort. */
    selection[current] = 1;
    size_t first = matrix->rowSlice[current];
    size_t beyond = matrix->rowSlice[current + 1];
    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] += matrix->entryValues[i];

    result = tuPartitionSubset(cmr, matrix, transposed, selection, current + 1, columnSum, stats, startClock, timeLimit);

    for (size_t i = first; i < beyond; ++i)
      columnSum[matrix->entryColumns[i]] -= matrix->entryValues[i];

    return result;
  }
  else
  {
    if ((clock() - startClock) * 1.0 / CLOCKS_PER_SEC > timeLimit)
      return -1;

#if defined(CMR_DEBUG)
    size_t count = 0;
    for (size_t row = 0; row < matrix->numRows; ++row)
    {
      if (selection[row])
        ++count;
    }
    CMRdbgMsg(4, "-> Selected a %s subset of cardinality %zu.\n", transposed ? "column" : "row", count);
#endif /* CMR_DEBUG */

    if (stats)
    {
      if (transposed)
        stats->partitionColumnSubsets++;
      else
        stats->partitionRowSubsets++;
    }

    bool foundPartition = tuPartitionSearch(cmr, matrix, selection, 0, columnSum);
    return foundPartition ? 1 : 0;
  }
}

/**
 * \brief Partition test based on Ghouila-Houri.
 */

static
CMR_ERROR tuPartition(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool transposed,            /**< Whether we're dealing with the transposed .*/
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  CMR_TU_STATS* stats,        /**< Statistics. */
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
    CMR_CALL( tuPartition(cmr, transpose, true, pisTotallyUnimodular, stats, timeLimit) );
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
  int result = tuPartitionSubset(cmr, matrix, transposed, selection, 0, columnSum, stats, startClock, timeLimit);
  if (result < 0)
    error = CMR_ERROR_TIMEOUT;
  else
    *pisTotallyUnimodular = (result > 0);

  CMR_CALL( CMRfreeStackArray(cmr, &columnSum) );
  CMR_CALL( CMRfreeStackArray(cmr, &selection) );

  if (stats)
    stats->partitionTime += (clock() - startClock) * 1.0 / CLOCKS_PER_SEC;

  return error;
}

CMR_ERROR CMRtuTest(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_SEYMOUR_NODE** proot,
  CMR_SUBMAT** psubmatrix, CMR_TU_PARAMS* params, CMR_TU_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  CMRdbgConsistencyAssert( CMRchrmatConsistency(matrix) );

  CMR_TU_PARAMS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRtuParamsInit(&defaultParams) );
    params = &defaultParams;
  }

  clock_t totalClock = clock();

  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
  {
    if (pisTotallyUnimodular)
      *pisTotallyUnimodular = false;
    return CMR_OKAY;
  }

  double remainingTime = timeLimit - ((clock() - totalClock) * 1.0 / CLOCKS_PER_SEC);

  CMRdbgMsg(0, "CMRtuTest called with algorithm = %d.\n", params->algorithm);

  if (params->algorithm == CMR_TU_ALGORITHM_DECOMPOSITION)
  {
    if (!params->ternary && params->camionFirst)
    {
      CMRdbgMsg(2, "Testing Camion signs before constructing a Seymour decomposition.\n");
      CMR_CALL( CMRcamionTestSigns(cmr, matrix, pisTotallyUnimodular, psubmatrix,
        stats ? &stats->camion : NULL, remainingTime) );

      if (!*pisTotallyUnimodular)
        return CMR_OKAY;
    }

    double remainingTime = timeLimit - ((clock() - totalClock) * 1.0 / CLOCKS_PER_SEC);

    CMR_SEYMOUR_NODE* root = NULL;
    CMR_ERROR error = CMRseymourDecompose(cmr, matrix, params->ternary, &root, &(params->seymour),
      stats ? &stats->seymour : NULL, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
    {
      assert( root == NULL);
      return error;
    }
    CMR_CALL( error );

    int8_t regularity = CMRseymourRegularity(root);
    if (regularity != 0)
      *pisTotallyUnimodular = regularity > 0;
    if (proot)
      *proot = root;
    else
      CMR_CALL( CMRseymourRelease(cmr, &root) );

    if (regularity < 0 && psubmatrix)
    {
      assert(!*psubmatrix);
      remainingTime = timeLimit - (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
      if (params->naiveSubmatrix)
        CMR_CALL( CMRtestHereditaryPropertyNaive(cmr, matrix, tuDecomposition, stats, psubmatrix, remainingTime) );
      else
        CMR_CALL( CMRtestHereditaryPropertyGreedy(cmr, matrix, tuDecomposition, stats, psubmatrix, remainingTime) );

      return CMR_OKAY;
    }

    if (regularity > 0 && !params->ternary && !params->camionFirst)
    {
      CMRdbgMsg(2, "Testing Camion signs afterward constructing a Seymour decomposition.\n");
      CMR_CALL( CMRcamionTestSigns(cmr, matrix, pisTotallyUnimodular, psubmatrix,
        stats ? &stats->camion : NULL, remainingTime) );
    }

  }
  else if (params->algorithm == CMR_TU_ALGORITHM_EULERIAN)
  {
    CMR_ERROR error = tuEulerian(cmr, matrix, false, pisTotallyUnimodular, psubmatrix, stats, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
      return error;
    CMR_CALL( error );
  }
  else if (params->algorithm == CMR_TU_ALGORITHM_PARTITION)
  {
    CMR_ERROR error = tuPartition(cmr, matrix, false, pisTotallyUnimodular, stats, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
      return error;
    CMR_CALL( error );
  }
  else
  {
    return CMR_ERROR_INVALID;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtuCompleteDecomposition(CMR* cmr, CMR_SEYMOUR_NODE* dec, CMR_TU_PARAMS* params, CMR_TU_STATS* stats,
  double timeLimit)
{
  assert(cmr);
  assert(dec);
  assert(timeLimit > 0);

  CMR_TU_PARAMS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRtuParamsInit(&defaultParams) );
    params = &defaultParams;
  }

  if (params->algorithm != CMR_TU_ALGORITHM_DECOMPOSITION)
    return CMR_ERROR_INPUT;

  if (!CMRseymourIsTernary(dec))
    return CMR_ERROR_INPUT;

  CMR_CALL( CMRregularityCompleteDecomposition(cmr, dec, &params->seymour, stats ? &stats->seymour : NULL,
    timeLimit) );

  return CMR_OKAY;
}

