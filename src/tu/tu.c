#include <tu/tu.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "sign_internal.h"

#include <stdlib.h>
#include <assert.h>


/**
 * \brief Tests the 1-sum of char matrices for total unimodularity.
 *
 * Returns \c true if and only if \p matrix is TU.
 *
 * If \p decomposition is not \c NULL and the algorithm has to test regularity of the support
 * matrix, then \c *decomposition will point to a decomposition tree for which the caller must use
 * \ref TUfreeDec to free memory. It is set to \c NULL otherwise.
 *
 * If \p submatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute
 * determinant larger than 1 will be searched, which may cause extra computational effort. In this
 * case, \c *submatrix will point to this submatrix for which the caller must use
 * \ref TUfreeSubmatrix to free memory. It is set to \c NULL otherwise.
 */

static bool testTotalUnimodularityOneSum(TU* tu, int numComponents,
  TU_ONESUM_COMPONENT* components, TU_DEC** decomposition, TU_SUBMATRIX** submatrix)
{
  assert(numComponents >= 0);
  assert(components);

  /* Sort components by number of nonzeros. */

  /* Check each component for base cases first. If enabled, ask for forbidden minors. */

  /* Check non-base cases by decomposing, using information about minors. If enabled, ask for
     rows/columns that cause irregularity. */

  /* If not regular and submatrix is required, start search in smallest non-regular component. */

  assert(false);

  for (int c = 0; c < numComponents; ++c)
  {
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }

  return true;
}

bool TUtestTotalUnimodularityDouble(TU* tu, TU_DBLMAT* matrix, double epsilon,
  TU_DEC** decomposition, TU_SUBMATRIX** submatrix)
{
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(tu);
  assert(matrix);

  /* Check entries. */

  if (!TUisTernaryDouble(tu, matrix, epsilon, submatrix))
    return false;

  /* Perform 1-sum decomposition. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(double), sizeof(char), &numComponents,
    &components, NULL, NULL, NULL, NULL);

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* compSubmatrix;
    char signFailed = signSequentiallyConnected(tu, (TU_CHRMAT*) &components[comp].matrix,
      (TU_CHRMAT*) &components[comp].transpose, false, submatrix ? &compSubmatrix : NULL);

    if (signFailed)
    {
      if (submatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        TUsortSubmatrix(compSubmatrix);
        *submatrix = compSubmatrix;
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

  return testTotalUnimodularityOneSum(tu, numComponents, components, decomposition, submatrix);
}

bool TUtestTotalUnimodularityInt(TU* tu, TU_INTMAT* matrix, TU_DEC** decomposition,
  TU_SUBMATRIX** submatrix)
{
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(tu);
  assert(matrix);

  /* Check entries. */

  if (!TUisTernaryInt(tu, matrix, submatrix))
    return false;

  /* Perform 1-sum decomposition. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(int), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* compSubmatrix;
    char signFailed = signSequentiallyConnected(tu, (TU_CHRMAT*) &components[comp].matrix,
      (TU_CHRMAT*) &components[comp].transpose, false, submatrix ? &compSubmatrix : NULL);

    if (signFailed)
    {
      if (submatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        TUsortSubmatrix(compSubmatrix);
        *submatrix = compSubmatrix;
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

  return testTotalUnimodularityOneSum(tu, numComponents, components, decomposition, submatrix);
}

bool TUtestTotalUnimodularityChar(TU* tu, TU_CHRMAT* matrix, TU_DEC** decomposition,
  TU_SUBMATRIX** submatrix)
{
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(tu);
  assert(matrix);

  /* Check entries. */

  if (!TUisTernaryChar(tu, matrix, submatrix))
    return false;

  /* Perform 1-sum decomposition. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* compSubmatrix;
    char signFailed = signSequentiallyConnected(tu, (TU_CHRMAT*) &components[comp].matrix,
      (TU_CHRMAT*) &components[comp].transpose, false, submatrix ? &compSubmatrix : NULL);

    if (signFailed)
    {
      if (submatrix)
      {
        /* Translate component indices to indices of whole matrix and sort them again. */
        for (int r = 0; r < compSubmatrix->numRows; ++r)
          compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
        for (int c = 0; c < compSubmatrix->numColumns; ++c)
          compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
        TUsortSubmatrix(compSubmatrix);
        *submatrix = compSubmatrix;
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

  return testTotalUnimodularityOneSum(tu, numComponents, components, decomposition, submatrix);
}
