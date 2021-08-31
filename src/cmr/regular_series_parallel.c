// #define CMR_DEBUG /** Uncomment to debug this file. */

#include "regular_internal.h"

#include <cmr/series_parallel.h>

#include "env_internal.h"
#include "dec_internal.h"

CMR_ERROR CMRregularDecomposeSeriesParallel(CMR* cmr, CMR_DEC** pdec, bool ternary, CMR_SUBMAT** psubmatrix,
  CMR_REGULAR_PARAMETERS* params)
{
  assert(cmr);
  assert(pdec);
  assert(*pdec);
  assert(params);

  CMR_DEC* dec = *pdec;

  bool isSeriesParallel;
  CMR_SP_REDUCTION* reductions = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &reductions, dec->matrix->numRows + dec->matrix->numColumns) );
  size_t numReductions;

  CMR_SUBMAT* reducedSubmatrix = NULL;
  CMR_SEPA* separation = NULL;
  if (ternary)
  {
    CMR_CALL( CMRdecomposeTernarySeriesParallel(cmr, dec->matrix, &isSeriesParallel, reductions, &numReductions,
      &reducedSubmatrix, psubmatrix, &separation, NULL) );
  }
  else
  {
    CMR_CALL( CMRdecomposeBinarySeriesParallel(cmr, dec->matrix, &isSeriesParallel, reductions, &numReductions,
      &reducedSubmatrix, psubmatrix, &separation, NULL) );
  }

  /* Did we find a 2-by-2 submatrix? If yes, it has determinant -2 or +2? */
  if (psubmatrix && *psubmatrix && (*psubmatrix)->numRows == 2)
    dec->type = CMR_DEC_IRREGULAR;

  /* Modify the decomposition to reflect the SP reductions. */
  if (numReductions > 0 && dec->type != CMR_DEC_IRREGULAR)
  {
    dec->type = CMR_DEC_SERIES_PARALLEL;
    dec->numReductions = numReductions;
    CMR_CALL( CMRduplicateBlockArray(cmr, &dec->reductions, numReductions, reductions) );

    if (!isSeriesParallel)
    {
      CMR_CALL( CMRdecSetNumChildren(cmr, dec, 1) );
      CMR_CALL( CMRdecCreate(cmr, dec, reducedSubmatrix->numRows, reducedSubmatrix->rows, reducedSubmatrix->numColumns,
        reducedSubmatrix->columns, &dec->children[0]) );
      CMR_CALL( CMRchrmatFilterSubmat(cmr, dec->matrix, reducedSubmatrix, &dec->children[0]->matrix) );
      dec = dec->children[0];
      *pdec = dec;
    }
  }

  /* Modify the decomposition for the 2-separation. */
  if (separation)
  {
    dec->type = CMR_DEC_TWO_SUM;
    CMR_CALL(  CMRdecSetNumChildren(cmr, dec, 2) );

    
  }

  CMR_CALL( CMRfreeStackArray(cmr, &reductions) );
  CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );  

  return CMR_OKAY;
}
