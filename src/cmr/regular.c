#define CMR_DEBUG /** Uncomment to debug the regularity check. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "dec_internal.h"
#include "regular_internal.h"

CMR_ERROR CMRtestRegularOneConnected(CMR* cmr, CMR_DEC* dec, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose, bool ternary,
  bool *pisRegular, CMR_MINOR** pminor, bool checkPlanarity, bool completeTree)
{
  assert(cmr);
  assert(dec);
  assert(matrix);

  CMRdbgMsg(2, "Testing binary %dx%d 1-connected matrix for regularity.\n", matrix->numRows, matrix->numColumns);

  CMRdbgMsg(4, "Splitting off series-parallel elements.\n");
  CMR_CALL( CMRregularDecomposeSeriesParallel(cmr, &dec, &matrix, ternary) );

  assert(false);

  return CMR_OKAY;
}

CMR_ERROR CMRtestRegular(CMR* cmr, CMR_CHRMAT* matrix, bool ternary, bool *pisRegular, CMR_DEC** pdec,
  CMR_MINOR** pminor, bool checkPlanarity, bool completeTree)
{
  assert(cmr);
  assert(matrix);

  CMRdbgMsg(0, "Testing %s %dx%d matrix for regularity.\n", ternary ? "ternary" : "binary", matrix->numRows,
    matrix->numColumns);

  CMR_DEC* dec = NULL;
  CMR_CALL( CMRdecCreate(cmr, NULL, matrix->numRows, NULL, matrix->numColumns, NULL, &dec) );
  assert(dec);

  CMR_CALL( CMRregularDecomposeOneSum(cmr, dec, matrix) );

  CMRdbgMsg(2, "1-sum decomposition yields %d components.\n", dec->numChildren == 0 ? 1 : dec->numChildren);
#if defined(CMR_DEBUG)
  CMR_CALL( CMRdecPrint(stdout, dec, 2) );
#endif /* CMR_DEBUG */

  bool isRegular = true;
  if (dec->numChildren)
  {
    for (size_t c = 0; c < dec->numChildren; ++c)
    {
      if (isRegular || completeTree)
      {
        bool childIsRegular = true;
        CMR_CALL( CMRtestRegularOneConnected(cmr, dec, dec->children[c]->matrix, dec->children[c]->transpose, ternary,
          &childIsRegular, pminor, checkPlanarity, completeTree) );
        if (*pminor)
          CMR_CALL( CMRdecTranslateMinorToParent(dec->children[c], *pminor) );

        isRegular = isRegular && childIsRegular;
      }
    }
  }
  else
  {
    CMR_CALL( CMRtestRegularOneConnected(cmr, dec, matrix, NULL, ternary, &isRegular, pminor, checkPlanarity,
      completeTree) );
  }

  if (pisRegular)
    *pisRegular = isRegular;

  CMR_CALL( CMRdecComputeRegularity(dec) );

  *pdec = dec;

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
  bool checkPlanarity, bool completeTree)
{
  assert(cmr);
  assert(matrix);

  CMR_SUBMAT* submatrix = NULL;
  if (!CMRchrmatIsBinary(cmr, matrix, &submatrix))
  {
    CMR_CALL( CMRminorCreate(cmr, pminor, 0, submatrix) );
    return CMR_OKAY;
  }

  CMR_CALL( CMRtestRegular(cmr, matrix, false, pisRegular, pdec, pminor, checkPlanarity, completeTree) );

  return CMR_OKAY;
}
