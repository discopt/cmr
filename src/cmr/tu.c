#include <cmr/tu.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "sign_internal.h"

#include <stdlib.h>
#include <assert.h>


/**
 * \brief Tests the 1-sum of char matrices for total unimodularity.
 *
 * Returns \c true if and only if the 1-sum of the given \p components is TU.
 *
 * If \p pdecomposition is not \c NULL and the algorithm has to test regularity of the support matrix, then
 * \c *pdecomposition will point to a decomposition tree for which the caller must use \ref CMRdecFree to free memory.
 * It is set to \c NULL in case regularity of the support matrix does not need to be determined.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute determinant larger than
 * 1 will be searched, which may cause extra computational effort. In this case, *\p psubmatrix will point to this
 * submatrix for which the caller must use \ref CMRsubmatFree to free memory. It is set to \c NULL otherwise.
 */

static
CMR_ERROR testTotalUnimodularityOneSum(
  CMR* cmr,                         /**< \ref CMR environment. */
  int numComponents,                /**< Number of 1-connected components. */
  CMR_ONESUM_COMPONENT* components, /**< 1-sum decomposition of matrix to be tested. */
  bool* pisTU,                      /**< Pointer for storing whether matrix is TU.*/
  CMR_TU_DEC** pdec,                /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix           /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
)
{
  assert(numComponents >= 0);
  assert(components);

  /* Sort components by number of nonzeros. */

  /* Check each component for base cases first. If enabled, ask for forbidden minors. */

  /* Check non-base cases by decomposing, using information about minors. If enabled, ask for
     rows/columns that cause irregularity. */

  /* If not regular and submatrix is required, start search in smallest non-regular component. */

  assert("TU test not implemented, yet." == 0);

  for (int c = 0; c < numComponents; ++c)
  {
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
    CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestTotalUnimodularityDbl(CMR* cmr, CMR_DBLMAT* matrix, double epsilon, bool* pisTU, CMR_TU_DEC** pdec,
  CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);

  int numComponents;
  CMR_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!CMRisTernaryDbl(cmr, matrix, epsilon, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(double), sizeof(char), &numComponents,
    &components, NULL, NULL, NULL, NULL);

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMR_SUBMAT* compSubmatrix;
    char modification;
    CMR_CALL( signSequentiallyConnected(cmr, (CMR_CHRMAT*) &components[comp].matrix,
      (CMR_CHRMAT*) &components[comp].transpose, false, &modification, psubmatrix ? &compSubmatrix : NULL) );

    if (modification)
    {
      if (psubmatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        CMRsortSubmatrix(compSubmatrix);
        *psubmatrix = compSubmatrix;
      }

      for (int c = 0; c < numComponents; ++c)
      {
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
        CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
        CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
      }

      return false;
    }
  }

  return testTotalUnimodularityOneSum(cmr, numComponents, components, pisTU, pdec, psubmatrix);
}

CMR_ERROR CMRtestTotalUnimodularityInt(CMR* cmr, CMR_INTMAT* matrix, bool* pisTU, CMR_TU_DEC** pdec, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);

  int numComponents;
  CMR_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!CMRisTernaryInt(cmr, matrix, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(int), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMR_SUBMAT* compSubmatrix;
    char modified;
    CMR_CALL( signSequentiallyConnected(cmr, (CMR_CHRMAT*) &components[comp].matrix,
      (CMR_CHRMAT*) &components[comp].transpose, false, &modified, psubmatrix ? &compSubmatrix : NULL) );

    if (modified)
    {
      if (psubmatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        CMRsortSubmatrix(compSubmatrix);
        *psubmatrix = compSubmatrix;
      }

      for (int c = 0; c < numComponents; ++c)
      {
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
        CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
        CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
      }

      return false;
    }
  }

  return testTotalUnimodularityOneSum(cmr, numComponents, components, pisTU, pdec, psubmatrix);
}

CMR_ERROR CMRtestTotalUnimodularityChr(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTU, CMR_TU_DEC** pdec, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);

  int numComponents;
  CMR_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!CMRisTernaryChr(cmr, matrix, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMR_SUBMAT* compSubmatrix;
    char modified;
    signSequentiallyConnected(cmr, (CMR_CHRMAT*) &components[comp].matrix,
      (CMR_CHRMAT*) &components[comp].transpose, false, &modified, psubmatrix ? &compSubmatrix : NULL);

    if (modified)
    {
      if (psubmatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        CMRsortSubmatrix(compSubmatrix);
        *psubmatrix = compSubmatrix;
      }

      for (int c = 0; c < numComponents; ++c)
      {
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
        CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
        CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
        CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
      }

      return false;
    }
  }

  return testTotalUnimodularityOneSum(cmr, numComponents, components, pisTU, pdec, psubmatrix);
}
