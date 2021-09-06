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

  CMRdbgMsg(4, "Checking for graphicness...");
  CMR_SUBMAT* submatrix = NULL;
  bool isGraphic;
  CMR_CALL( CMRregularTestGraphic(cmr, &dec->matrix, &dec->transpose, ternary, &isGraphic, &dec->graph,
    &dec->graphForest, &dec->graphCoforest, &dec->graphArcsReversed, &submatrix) );
  if (isGraphic)
  {
    CMRdbgMsg(0, " graphic.\n");
    dec->type = CMR_DEC_GRAPHIC;
    if (!params->planarityCheck)
      return CMR_OKAY;
  }
  CMRdbgMsg(0, " NOT graphic.\n");

  CMRdbgMsg(4, "Checking for cographicness...");
  bool isCographic;
  CMR_CALL( CMRregularTestGraphic(cmr, &dec->transpose, &dec->matrix, ternary, &isCographic, &dec->cograph,
    &dec->cographForest, &dec->cographCoforest, &dec->cographArcsReversed, &submatrix) );
  if (isCographic)
  {
    CMRdbgMsg(0, " cographic.\n");
    dec->type = (dec->type == CMR_DEC_GRAPHIC) ? CMR_DEC_PLANAR : CMR_DEC_COGRAPHIC;
    return CMR_OKAY;
  }
  CMRdbgMsg(0, " NOT cographic.\n");

  if (submatrix)
    CMR_CALL( CMRsubmatTranspose(submatrix) );
  
  CMRdbgMsg(4, "Splitting off series-parallel elements...");
  CMR_CALL( CMRregularDecomposeSeriesParallel(cmr, &dec, ternary, &submatrix, params) );

  if (dec->type == CMR_DEC_IRREGULAR)
  {
    CMRdbgMsg(0, " NOT regular.\n");
    return CMR_OKAY;
  }

  if (dec->type == CMR_DEC_TWO_SUM)
  {
    CMRdbgMsg(0, " Encountered a 2-separation.\n");
    assert(dec->numChildren == 2);
    CMR_CALL( testRegularOneConnected(cmr, dec->children[0], ternary, pisRegular, pminor, params) );

    if (params->completeTree || *pisRegular)
      CMR_CALL( testRegularOneConnected(cmr, dec->children[1], ternary, pisRegular, pminor, params) );

    return CMR_OKAY;
  }

  CMRdbgMsg(0, " Found a W_k minor.\n");

  /* No 2-sum found, so we have a wheel submatrix. */

  // TODO: Extract a W_k minor.
  // TODO: Try out extracting a W_3 minor via pivots to then have a smaller non-(co)graphic minor, reducing search
  //       effort for 3-separations.

  CMR_CALL( CMRsubmatFree(cmr, &submatrix) );

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
        if (pminor && *pminor)
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
