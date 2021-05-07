#include <tu/tu.h>

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
 * \c *pdecomposition will point to a decomposition tree for which the caller must use \ref TUdecFree to free memory.
 * It is set to \c NULL in case regularity of the support matrix does not need to be determined.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute determinant larger than
 * 1 will be searched, which may cause extra computational effort. In this case, *\p psubmatrix will point to this
 * submatrix for which the caller must use \ref TUsubmatFree to free memory. It is set to \c NULL otherwise.
 */

static
TU_ERROR testTotalUnimodularityOneSum(
  TU* tu,                           /**< \ref TU environment. */
  int numComponents,                /**< Number of 1-connected components. */
  TU_ONESUM_COMPONENT* components,  /**< 1-sum decomposition of matrix to be tested. */
  bool* pisTU,                      /**< Pointer for storing whether matrix is TU.*/
  TU_DEC** pdec,          /**< Pointer for storing the decomposition tree (may be \c NULL). */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
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
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }

  return TU_OKAY;
}

TU_ERROR TUtestTotalUnimodularityDbl(TU* tu, TU_DBLMAT* matrix, double epsilon, bool* pisTU, TU_DEC** pdec,
  TU_SUBMAT** psubmatrix)
{
  assert(tu);
  assert(matrix);

  int numComponents;
  TU_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!TUisTernaryDbl(tu, matrix, epsilon, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(double), sizeof(char), &numComponents,
    &components, NULL, NULL, NULL, NULL);

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMAT* compSubmatrix;
    char modification;
    TU_CALL( signSequentiallyConnected(tu, (TU_CHRMAT*) &components[comp].matrix,
      (TU_CHRMAT*) &components[comp].transpose, false, &modification, psubmatrix ? &compSubmatrix : NULL) );

    if (modification)
    {
      if (psubmatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        TUsortSubmatrix(compSubmatrix);
        *psubmatrix = compSubmatrix;
      }

      for (int c = 0; c < numComponents; ++c)
      {
        TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
        TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
        TUfreeBlockArray(tu, &components[c].rowsToOriginal);
        TUfreeBlockArray(tu, &components[c].columnsToOriginal);
      }

      return false;
    }
  }

  return testTotalUnimodularityOneSum(tu, numComponents, components, pisTU, pdec, psubmatrix);
}

TU_ERROR TUtestTotalUnimodularityInt(TU* tu, TU_INTMAT* matrix, bool* pisTU, TU_DEC** pdec, TU_SUBMAT** psubmatrix)
{
  assert(tu);
  assert(matrix);

  int numComponents;
  TU_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!TUisTernaryInt(tu, matrix, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(int), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMAT* compSubmatrix;
    char modified;
    TU_CALL( signSequentiallyConnected(tu, (TU_CHRMAT*) &components[comp].matrix,
      (TU_CHRMAT*) &components[comp].transpose, false, &modified, psubmatrix ? &compSubmatrix : NULL) );

    if (modified)
    {
      if (psubmatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        TUsortSubmatrix(compSubmatrix);
        *psubmatrix = compSubmatrix;
      }

      for (int c = 0; c < numComponents; ++c)
      {
        TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
        TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
        TUfreeBlockArray(tu, &components[c].rowsToOriginal);
        TUfreeBlockArray(tu, &components[c].columnsToOriginal);
      }

      return false;
    }
  }

  return testTotalUnimodularityOneSum(tu, numComponents, components, pisTU, pdec, psubmatrix);
}

TU_ERROR TUtestTotalUnimodularityChr(TU* tu, TU_CHRMAT* matrix, bool* pisTU, TU_DEC** pdec, TU_SUBMAT** psubmatrix)
{
  assert(tu);
  assert(matrix);

  int numComponents;
  TU_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!TUisTernaryChr(tu, matrix, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMAT* compSubmatrix;
    char modified;
    signSequentiallyConnected(tu, (TU_CHRMAT*) &components[comp].matrix,
      (TU_CHRMAT*) &components[comp].transpose, false, &modified, psubmatrix ? &compSubmatrix : NULL);

    if (modified)
    {
      if (psubmatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        TUsortSubmatrix(compSubmatrix);
        *psubmatrix = compSubmatrix;
      }

      for (int c = 0; c < numComponents; ++c)
      {
        TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
        TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
        TUfreeBlockArray(tu, &components[c].rowsToOriginal);
        TUfreeBlockArray(tu, &components[c].columnsToOriginal);
      }

      return false;
    }
  }

  return testTotalUnimodularityOneSum(tu, numComponents, components, pisTU, pdec, psubmatrix);
}
