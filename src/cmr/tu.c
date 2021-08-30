#include <cmr/tu.h>

#include "matrix_internal.h"
#include "one_sum.h"
#include "camion_internal.h"
#include "regular_internal.h"

#include <stdlib.h>
#include <assert.h>

#include "interface.h"

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
  CMR_DEC** pdec,                /**< Pointer for storing the decomposition tree (may be \c NULL). */
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

CMR_ERROR CMRtestTotalUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisTotallyUnimodular, CMR_DEC** pdec,
  CMR_SUBMAT** psubmatrix, bool checkPlanarity, bool completeTree)
{
  assert(cmr);
  assert(matrix);

  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  CMR_CALL( CMRinterfaceTU(cmr, matrix, pisTotallyUnimodular, pdec, psubmatrix) );
  
  CMR_MINOR* minor = NULL;
//   CMR_CALL( CMRtestRegular(cmr, matrix, true, pisTotallyUnimodular, pdec, psubmatrix ? &minor : NULL, checkPlanarity,
//     completeTree) );
  if (minor)
  {
    assert(minor->numPivots == 0);
    *psubmatrix = minor->remainingSubmatrix;
    minor->remainingSubmatrix = NULL;
  }
  
  return CMR_OKAY;


  size_t numComponents;
  CMR_ONESUM_COMPONENT* components;

  /* Check entries. */

  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return false;

  /* Perform 1-sum decomposition. */

  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  /* Check correct signing for each component. */

  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMR_SUBMAT* compSubmatrix;
    char modified;
    CMR_CALL( CMRcomputeCamionSignSequentiallyConnected(cmr, (CMR_CHRMAT*) &components[comp].matrix,
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
        CMRsortSubmatrix(cmr, compSubmatrix);
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

  return testTotalUnimodularityOneSum(cmr, numComponents, components, pisTotallyUnimodular, pdec, psubmatrix);
}
