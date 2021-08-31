#define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "dec_internal.h"
#include "regular_internal.h"

CMR_ERROR CMRregularInitParameters(CMR_REGULAR_PARAMETERS* params)
{
  assert(params);

  params->planarityCheck = false;
  params->seriesParallel = 2;
  params->completeTree = false;
  params->matrices = CMR_DEC_CONSTRUCT_NONE;
  params->transposes = CMR_DEC_CONSTRUCT_NONE;
  params->graphs = CMR_DEC_CONSTRUCT_NONE;

  return CMR_OKAY;
}

static 
CMR_ERROR testRegularOneConnected(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC* dec,                   /**< Decomposition node. */
  bool ternary,                   /**< Whether signs matter. */
  bool *pisRegular,               /**< Pointer for storing whether \p matrix is regular. */
  CMR_MINOR** pminor,             /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMETERS* params  /**< Parameters for the computation. */
)
{
  assert(cmr);
  assert(dec);
  assert(dec->matrix);

  CMRdbgMsg(2, "Testing binary %dx%d 1-connected matrix for regularity.\n", dec->matrix->numRows,
    dec->matrix->numColumns);

  if (params->seriesParallel)
  {
    CMRdbgMsg(4, "Splitting off series-parallel elements.\n");
    CMR_SUBMAT* submatrix = NULL;
    CMR_CALL( CMRregularDecomposeSeriesParallel(cmr, &dec, ternary, &submatrix, params) );
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestRegular(CMR* cmr, CMR_CHRMAT* matrix, bool ternary, bool *pisRegular, CMR_DEC** pdec,
  CMR_MINOR** pminor, CMR_REGULAR_PARAMETERS* params)
{
  assert(cmr);
  assert(matrix);
  assert(params);

  CMRdbgMsg(0, "Testing %s %dx%d matrix for regularity.\n", ternary ? "ternary" : "binary", matrix->numRows,
    matrix->numColumns);

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
  if (dec->numChildren)
  {
    for (size_t c = 0; c < dec->numChildren; ++c)
    {
      if (isRegular || params->completeTree)
      {
        bool childIsRegular = true;
        CMR_CALL( testRegularOneConnected(cmr, dec->children[c], ternary, &childIsRegular, pminor, params) );
        if (*pminor)
          CMR_CALL( CMRdecTranslateMinorToParent(dec->children[c], *pminor) );

        isRegular = isRegular && childIsRegular;
      }
    }
  }
  else
  {
    CMR_CALL( testRegularOneConnected(cmr, dec, ternary, &isRegular, pminor, params) );
  }

  if (pisRegular)
    *pisRegular = isRegular;

  CMR_CALL( CMRdecComputeRegularity(dec) );

  /* If the root must have a matrix, we eventually copy the one passed by the user. */
  if (params->matrices == CMR_DEC_CONSTRUCT_ALL
    || (params->matrices == CMR_DEC_CONSTRUCT_LEAVES && dec->numChildren == 0))
  {
    CMR_CALL( CMRchrmatCopy(cmr, dec->matrix, &dec->matrix) );
  }
  else
    dec->matrix = NULL;

  /* Either store or free the decomposition. */
  if (pdec)
    *pdec = dec;
  else
    CMR_CALL( CMRdecFree(cmr, &dec) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestBinaryRegularConnected(CMR* cmr, CMR_DEC* dec, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose,
  bool checkPlanarity, bool certify, bool constructDecomposition, bool* pisRegular)
{
  assert(cmr);
  assert(dec);
  assert(matrix);

  CMRdbgMsg(2, "Testing binary %dx%d 1-connected matrix for regularity.\n", matrix->numRows, matrix->numColumns);


  return CMR_OKAY;
}

CMR_ERROR CMRtestBinaryRegular(CMR* cmr, CMR_CHRMAT* matrix, bool *pisRegular, CMR_DEC** pdec, CMR_MINOR** pminor,
  CMR_REGULAR_PARAMETERS* params)
{
  assert(cmr);
  assert(matrix);

  CMR_REGULAR_PARAMETERS defaultParams;
  if (!params)
  {
    CMR_CALL( CMRregularInitParameters(&defaultParams) );
    params = &defaultParams;
  }

  CMR_SUBMAT* submatrix = NULL;
  if (!CMRchrmatIsBinary(cmr, matrix, &submatrix))
  {
    CMR_CALL( CMRminorCreate(cmr, pminor, 0, submatrix) );
    return CMR_OKAY;
  }

  CMR_CALL( CMRtestRegular(cmr, matrix, false, pisRegular, pdec, pminor, params) );

  return CMR_OKAY;
}
