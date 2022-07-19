// #define CMR_DEBUG /** Uncomment to debug this file. */

#include "regular_internal.h"

#include <cmr/series_parallel.h>

#include "env_internal.h"
#include "dec_internal.h"

CMR_ERROR CMRregularDecomposeSeriesParallel(CMR* cmr, CMR_DEC** pdec, bool ternary, CMR_SUBMAT** psubmatrix,
  CMR_REGULAR_PARAMETERS* params, CMR_REGULAR_STATISTICS* stats)
{
  assert(cmr);
  assert(pdec);
  assert(*pdec);
  assert(params);

  CMR_DEC* dec = *pdec;
  
  CMRdbgMsg(0, "\nCMRregularDecomposeSeriesParallel called for matrix\n");
#if defined(CMR_DEBUG)
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  bool isSeriesParallel = true;
  CMR_SP_REDUCTION* reductions = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &reductions, dec->matrix->numRows + dec->matrix->numColumns) );
  size_t numReductions;

  CMR_SUBMAT* reducedSubmatrix = NULL;
  CMR_SEPA* separation = NULL;
  if (ternary)
  {
    CMR_CALL( CMRdecomposeTernarySeriesParallel(cmr, dec->matrix, &isSeriesParallel, reductions, &numReductions,
      &reducedSubmatrix, psubmatrix, &separation, stats ? &stats->seriesParallel : NULL) );
  }
  else
  {
    CMR_CALL( CMRdecomposeBinarySeriesParallel(cmr, dec->matrix, &isSeriesParallel, reductions, &numReductions,
      &reducedSubmatrix, psubmatrix, &separation, stats ? &stats->seriesParallel : NULL) );
  }

  /* Did we find a 2-by-2 submatrix? If yes, is has determinant -2 or +2! */
  if (psubmatrix && *psubmatrix && (*psubmatrix)->numRows == 2)
    dec->type = CMR_DEC_IRREGULAR;

  CMRdbgMsg(0, "[CMRregularDecomposeSeriesParallel -> %ld SP reductions]", numReductions);

  /* Modify the decomposition to reflect the SP reductions. */
  if (numReductions > 0 && dec->type != CMR_DEC_IRREGULAR)
  {
    if (params->seriesParallel)
    {
      dec->type = CMR_DEC_SERIES_PARALLEL;
      dec->numReductions = numReductions;
      CMR_CALL( CMRduplicateBlockArray(cmr, &dec->reductions, numReductions, reductions) );

      if (!isSeriesParallel)
      {
        CMR_CALL( CMRdecSetNumChildren(cmr, dec, 1) );
        CMR_CALL( CMRdecCreate(cmr, dec, reducedSubmatrix->numRows, reducedSubmatrix->rows, reducedSubmatrix->numColumns,
          reducedSubmatrix->columns, &dec->children[0]) );
        CMR_CALL( CMRchrmatZoomSubmat(cmr, dec->matrix, reducedSubmatrix, &dec->children[0]->matrix) );
        if (psubmatrix && *psubmatrix)
        {
          /* Zoom into submatrix. */
          CMR_SUBMAT* zoomedSubmatrix = NULL;
          CMR_CALL( CMRsubmatZoomSubmat(cmr, reducedSubmatrix, *psubmatrix, &zoomedSubmatrix));
          CMR_CALL( CMRsubmatFree(cmr, psubmatrix) );
          *psubmatrix = zoomedSubmatrix;
        }
        dec = dec->children[0];
        *pdec = dec;
      }
    }
    else
    {
      /* We have to carry out each SP reduction as a 2-separation. */
      CMRdbgMsg(0, "\n");
      for (size_t r = 0; r < numReductions; ++r)
      {
        char buffer[16];
        CMRdbgMsg(0, "Reduction %ld is (%s,%s).\n", r, CMRelementString(reductions[r].element, 0),
          CMRelementString(reductions[r].mate, buffer));
      }

      assert(!separation);
      size_t parentNumRows = dec->matrix->numRows;
      size_t parentNumColumns = dec->matrix->numColumns;
      CMR_CALL( CMRsepaCreate(cmr, parentNumRows, parentNumColumns, &separation) );
      for (size_t r = 0; r < parentNumRows; ++r)
        separation->rowsToPart[r] = 0;
      for (size_t c = 0; c < parentNumColumns; ++c)
        separation->columnsToPart[c] = 0;
      if (CMRspIsRow(reductions[0]))
        separation->rowsToPart[CMRelementToRowIndex(reductions[0].element)] = 1;
      else
        separation->columnsToPart[CMRelementToColumnIndex(reductions[0].element)] = 1;
      if (CMRelementIsRow(reductions[0].mate))
        separation->rowsToPart[CMRelementToRowIndex(reductions[0].mate)] = 1;
      else if (CMRelementIsColumn(reductions[0].mate))
        separation->columnsToPart[CMRelementToColumnIndex(reductions[0].mate)] = 1;
      CMR_CALL( CMRsepaInitializeMatrix(cmr, separation, dec->matrix, 1) );
    }
  }

  /* Modify the decomposition for the 2-separation. */
  if (separation)
    CMR_CALL( CMRdecApplySeparation(cmr, dec, separation) );

  CMR_CALL( CMRfreeStackArray(cmr, &reductions) );
  CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );  

  return CMR_OKAY;
}
