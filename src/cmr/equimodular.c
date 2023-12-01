#define CMR_DEBUG /* Uncomment to debug. */

#include <cmr/equimodular.h>
#include <cmr/tu.h>

#include <assert.h>
#include <float.h>
#include <time.h>

#include "env_internal.h"
#include "linalg.h"

CMR_ERROR CMRparamsEquimodularityInit(CMR_EQUIMODULAR_PARAMETERS* params)
{
  assert(params);

  CMR_CALL( CMRparamsTotalUnimodularityInit(&params->tu) );

  return CMR_OKAY;
}

CMR_ERROR CMRstatsEquimodularityInit(CMR_EQUIMODULAR_STATISTICS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  stats->linalgTime = 0.0;
  CMR_CALL( CMRstatsTotalUnimodularityInit(&stats->tu) );

  return CMR_OKAY;
}

CMR_ERROR CMRstatsEquimodularityPrint(FILE* stream, CMR_EQUIMODULAR_STATISTICS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Equimodularity recognition:\n");
    prefix = "  ";
  }

  fprintf(stream, "%slinear algebra: %f seconds\n", prefix, stats->linalgTime);

  char subPrefix[256];
  snprintf(subPrefix, 256, "%stu ", prefix);
  CMR_CALL( CMRstatsTotalUnimodularityPrint(stream, &stats->tu, subPrefix) );

  fprintf(stream, "%stotal: %ld in %f seconds\n", prefix, stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}

CMR_ERROR CMRtestEquimodularity(CMR* cmr, CMR_INTMAT* matrix, bool* pisEquimodular, int64_t* pgcdDet,
  CMR_EQUIMODULAR_PARAMETERS* params, CMR_EQUIMODULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(pisEquimodular);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );

  CMR_EQUIMODULAR_PARAMETERS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRparamsEquimodularityInit(&defaultParams) );
    params = &defaultParams;
  }
  if (stats)
    stats->totalCount++;

  CMR_ERROR result = CMR_OKAY;

  clock_t totalClock = clock();

  /* Transform matrix to upper-diagonal matrix with diagonally dominant columns. */
  size_t rank;
  CMR_SUBMAT* basisPermutation = NULL;
  CMR_INTMAT* transformed_matrix = NULL;
  CMR_INTMAT* transformed_transpose = NULL;
  CMR_ERROR error = CMRintmatComputeUpperDiagonal(cmr, matrix, true, &rank, &basisPermutation, &transformed_matrix,
    &transformed_transpose);
  if (error == CMR_ERROR_OVERFLOW)
    return CMR_ERROR_OVERFLOW;
  CMR_CALL(error);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "Transformed matrix has rank %ld.\n", rank);
  CMRintmatPrintDense(cmr, transformed_matrix, stdout, '0', false);
  CMRintmatPrintDense(cmr, transformed_transpose, stdout, '0', false);
  fflush(stdout);
#endif /* CMR_DEBUG */

  assert(rank <= matrix->numRows);
  int64_t gcdDet = 1;
  for (size_t row = 0; row < rank; ++row)
  {
    assert(transformed_matrix->rowSlice[row] < transformed_matrix->rowSlice[row + 1]);
    assert(transformed_matrix->entryColumns[transformed_matrix->rowSlice[row]] == row);

    int diagonalElement = transformed_matrix->entryValues[transformed_matrix->rowSlice[row]];
    assert(diagonalElement > 0);
    int64_t old = gcdDet;
    gcdDet *= diagonalElement;
    if (gcdDet / diagonalElement != old)
    {
      /* We caught an overflow. */
      result = CMR_ERROR_OVERFLOW;
      goto cleanup;
    }
  }

  /* Test for a particular determinant gcd if requested. */
  CMRdbgMsg(2, "The determinant gcd is %ld.\n", gcdDet);
  if (pgcdDet && *pgcdDet && *pgcdDet != gcdDet)
  {
    *pisEquimodular = false;
    *pgcdDet = gcdDet;
    goto cleanup;
  }

  /* Construct the transpose (since we create it column-wises) of the pseudo-inverse that we feed into the TU test. */
  int64_t* denseColumn = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumn, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    denseColumn[column] = 0;
  size_t numDenseColumnNonzeros;
  size_t* denseColumnNonzeros = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &denseColumnNonzeros, matrix->numColumns) );
  CMR_CHRMAT* transposed_pseudo_inverse = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &transposed_pseudo_inverse, matrix->numColumns, matrix->numRows,
    matrix->numRows + 2 * transformed_matrix->numNonzeros) );

  size_t numNonzeros = 0;
  transposed_pseudo_inverse->rowSlice[0] = 0;
  for (size_t column = 0; column < matrix->numColumns; ++column)
  {
    /* We first copy the entries of the column into denseColumn. */
    size_t begin = transformed_transpose->rowSlice[column];
    size_t beyond = transformed_transpose->rowSlice[column + 1];
    for (size_t i = begin; i < beyond; ++i)
      denseColumn[transformed_transpose->entryColumns[i]] = transformed_transpose->entryValues[i];
    numDenseColumnNonzeros = 0;

    /* We the go through all pivots in the basis matrix. */
    for (size_t i = rank; i > 0; --i)
    {
      size_t pivotRow = i - 1;
      size_t begin = transformed_matrix->rowSlice[pivotRow] + 1;
      size_t beyond = transformed_matrix->rowSlice[pivotRow + 1];
      int64_t numerator = denseColumn[pivotRow];
      for (size_t j = begin + 1; j < beyond; ++j)
      {
        size_t basisColumn = transformed_matrix->entryColumns[j];
        if (basisColumn >= rank)
          break;

        numerator -= transformed_matrix->entryValues[j] * denseColumn[pivotRow];
        /* TODO: check for overflow here! */
      }

      /* Entry in column should be -1, 0 or +1. */
      int64_t denominator = transformed_matrix->entryValues[transformed_matrix->rowSlice[pivotRow]];
      if (numerator != 0 && numerator != denominator && numerator != -denominator)
      {
        /* Resulting entry is not in {-1, 0, +1}. */
        *pisEquimodular = false;
        CMR_CALL( CMRfreeStackArray(cmr, &denseColumnNonzeros) );
        CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );
        CMRchrmatFree(cmr, &transposed_pseudo_inverse);
        goto cleanup;
      }

      char x = numerator / denominator;
      denseColumn[pivotRow] = x;
      if (x)
        denseColumnNonzeros[numDenseColumnNonzeros++] = pivotRow;
    }

    /* Ensure that we have enough memory for the column and add it to the pseudo-inverse. */
    numNonzeros += numDenseColumnNonzeros;
    if (numNonzeros > transposed_pseudo_inverse->numNonzeros)
      CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, transposed_pseudo_inverse, 2 * numNonzeros) );
    begin = transposed_pseudo_inverse->rowSlice[column];
    transposed_pseudo_inverse->rowSlice[column + 1] = numNonzeros;
    for (size_t j = 0; j < numDenseColumnNonzeros; ++j)
    {
      transposed_pseudo_inverse->entryColumns[numNonzeros - j - 1] = denseColumnNonzeros[j];
      transposed_pseudo_inverse->entryValues[numNonzeros - j - 1] = denseColumn[denseColumnNonzeros[j]];
    }

    for (size_t i = 0; i < numDenseColumnNonzeros; ++i)
      denseColumn[denseColumnNonzeros[i]] = 0;
  }
  transposed_pseudo_inverse->numNonzeros = numNonzeros;


#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "Transposed pseudo-inverse is:\n");
  CMRchrmatPrintDense(cmr, transposed_pseudo_inverse, stdout, '0', false);
  fflush(stdout);
#endif /* CMR_DEBUG */

  CMR_CALL( CMRfreeStackArray(cmr, &denseColumnNonzeros) );
  CMR_CALL( CMRfreeStackArray(cmr, &denseColumn) );

  double linalgTime = (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  if (stats)
    stats->linalgTime += linalgTime;
  double remainingTime = timeLimit - linalgTime;
  if (remainingTime <= 0)
  {
    CMRchrmatFree(cmr, &transposed_pseudo_inverse);
    result = CMR_ERROR_TIMEOUT;
    goto cleanup;
  }

  CMR_CALL( CMRtestTotalUnimodularity(cmr, transposed_pseudo_inverse, pisEquimodular, NULL, NULL, &params->tu,
    stats ? &stats->tu : NULL, remainingTime) );

  if (pgcdDet)
    *pgcdDet = (*pisEquimodular) ? gcdDet : 0;

  CMRchrmatFree(cmr, &transposed_pseudo_inverse);

cleanup:

  CMRintmatFree(cmr, &transformed_transpose);
  CMRintmatFree(cmr, &transformed_matrix);
  CMRsubmatFree(cmr, &basisPermutation);

  if (stats)
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;

  return result;
  
}

CMR_ERROR CMRtestStrongEquimodularity(CMR* cmr, CMR_INTMAT* matrix, bool* pisStronglyEquimodular, int64_t* pgcdDet,
  CMR_EQUIMODULAR_PARAMETERS* params, CMR_EQUIMODULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(pisStronglyEquimodular);
  assert(params);
  assert(stats);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );

  /* If not supplied, we maintain a local gcd determinant variable. */
  int64_t gcdDet = 0;
  if (!pgcdDet)
    pgcdDet = &gcdDet;

  clock_t startClock = clock();
  CMR_CALL( CMRtestEquimodularity(cmr, matrix, pisStronglyEquimodular, pgcdDet, params, stats, timeLimit) );
  double remainingTime = timeLimit - ((clock() - startClock) * 1.0 / CLOCKS_PER_SEC);
  if (remainingTime <= 0)
    return CMR_ERROR_TIMEOUT;

  /* If it is equimodular, then we also test the transpose. */
  if (*pisStronglyEquimodular)
  {
    CMR_INTMAT* transpose = NULL;
    CMR_CALL( CMRintmatTranspose(cmr, matrix, &transpose) );
    CMR_CALL( CMRtestEquimodularity(cmr, transpose, pisStronglyEquimodular, pgcdDet, params, stats, remainingTime) );
    CMR_CALL( CMRintmatFree(cmr, &transpose) );
  }

  return CMR_OKAY;
}


CMR_EXPORT
CMR_ERROR CMRtestUnimodularity(CMR* cmr, CMR_INTMAT* matrix, bool* pisUnimodular, CMR_EQUIMODULAR_PARAMETERS* params,
  CMR_EQUIMODULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(pisUnimodular);
  assert(params);
  assert(stats);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );

  int64_t gcdDet = 1;
  CMR_CALL( CMRtestEquimodularity(cmr, matrix, pisUnimodular, &gcdDet, params, stats, timeLimit) );

  return CMR_OKAY;
}

CMR_EXPORT
CMR_ERROR CMRtestStrongUnimodularity(CMR* cmr, CMR_INTMAT* matrix, bool* pisStronglyUnimodular,
  CMR_EQUIMODULAR_PARAMETERS* params, CMR_EQUIMODULAR_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(pisStronglyUnimodular);
  assert(params);
  assert(stats);
  CMRconsistencyAssert( CMRintmatConsistency(matrix) );

  int64_t gcdDet = 1;
  CMR_CALL( CMRtestStrongEquimodularity(cmr, matrix, pisStronglyUnimodular, &gcdDet, params, stats, timeLimit) );

  return CMR_OKAY;
}
