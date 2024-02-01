#define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "matroid_internal.h"
#include "regularity_internal.h"

CMR_ERROR CMRregularityDecomposeThreeSum(
  CMR* cmr,
  DecompositionTask* task,
  DecompositionTask** punprocessed,
  CMR_SEPA* separation
)
{
  assert(cmr);
  assert(task);
  assert(punprocessed);
  assert(separation);

  CMR_MATROID_DEC* dec = task->dec;
  assert(dec);
  size_t* rowsToChild = NULL;
  size_t* columnsToChild = NULL;

#if defined(CMR_DEBUG)
  CMRdbgMsg(8, "Processing 3-separation to produce a 3-sum for the following matrix:\n");
  CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true);
  for (size_t row = 0; row < dec->numRows; ++row)
  {
    int part = ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? 0 : 1;
    int rank1Flag = (separation->rowsFlags[row] & CMR_SEPA_FLAG_RANK1) ? 1 : 0;
    int rank2Flag = (separation->rowsFlags[row] & CMR_SEPA_FLAG_RANK2) ? 1 : 0;
    CMRdbgMsg(10, "Row r%zu belongs to part %d; flags = %d/%d\n", row+1, part, rank1Flag, rank2Flag);
  }
  for (size_t column = 0; column < dec->numColumns; ++column)
  {
    int part = ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? 0 : 1;
    int rank1Flag = (separation->columnsFlags[column] & CMR_SEPA_FLAG_RANK1) ? 1 : 0;
    int rank2Flag = (separation->columnsFlags[column] & CMR_SEPA_FLAG_RANK2) ? 1 : 0;
    CMRdbgMsg(10, "Column c%zu belongs to part %d; flags = %d/%d\n", column+1, part, rank1Flag, rank2Flag);
  }
#endif /* CMR_DEBUG */

  size_t numPivots = 0;
  size_t* pivotRows = NULL;
  size_t* pivotColumns = NULL;
  size_t maxNumPivots = dec->numRows < dec->numColumns ? dec->numRows : dec->numColumns;
  CMR_CALL( CMRallocStackArray(cmr, &pivotRows, maxNumPivots) );
  CMR_CALL( CMRallocStackArray(cmr, &pivotColumns, maxNumPivots) );

  size_t pivotRow = SIZE_MAX;
  size_t pivotColumn = SIZE_MAX;

  size_t extraRows[2][3];
  size_t extraColumns[2][3];
  CMR_CALL( CMRsepaGetRepresentatives(cmr, separation, extraRows, extraColumns) );

  if (separation->type == CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS &&
    (task->params->threeSumStrategy & CMR_MATROID_DEC_THREESUM_FLAG_CONCENTRATED_RANK))
  {
    /* We have distributed ranks but requested concentrated rank; we find a top-right nonzero. */
    pivotRow = extraRows[1][0];
    assert(pivotRow != SIZE_MAX);

    pivotColumn = extraColumns[0][0];
    assert(pivotColumn != SIZE_MAX);
  }

  if (separation->type == CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK &&
    (task->params->threeSumStrategy & CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS))
  {
    /* We have concentrated rank but requested distributed ranks; we find a bottom-left nonzero. */
    pivotRow = extraRows[0][0];
    assert(pivotRow != SIZE_MAX);

    pivotColumn = extraColumns[1][0];
    assert(pivotColumn != SIZE_MAX);
  }

  CMR_CHRMAT* goodRankMatrix = NULL;
  CMR_CHRMAT* goodRankTranspose = NULL;
  if (pivotRow != SIZE_MAX)
  {
    CMRdbgMsg(10, "-> We pivot on r%zu,c%zu.\n", pivotRow+1, pivotColumn+1);

    if (pivotColumns)
    {
      pivotRows[0] = pivotRow;
      pivotColumns[0] = pivotColumn;
      ++numPivots;
    }

    if (dec->isTernary)
      CMR_CALL( CMRchrmatTernaryPivot(cmr, dec->matrix, pivotRow, pivotColumn, &goodRankMatrix) );
    else
      CMR_CALL( CMRchrmatBinaryPivot(cmr, dec->matrix, pivotRow, pivotColumn, &goodRankMatrix) );

    CMR_CALL( CMRchrmatTranspose(cmr, goodRankMatrix, &goodRankTranspose) );

    /* Swap association to first/second for pivot row and column. */
    if ((separation->rowsFlags[pivotRow] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      separation->rowsFlags[pivotRow] = CMR_SEPA_SECOND;
    else
      separation->rowsFlags[pivotRow] = CMR_SEPA_FIRST;
    if ((separation->columnsFlags[pivotColumn] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      separation->columnsFlags[pivotColumn] = CMR_SEPA_SECOND;
    else
      separation->columnsFlags[pivotColumn] = CMR_SEPA_FIRST;

    /* We recompute all representatives of the low-rank submatrices. */
    CMR_SUBMAT* violatorSubmatrix = NULL;
    CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, goodRankMatrix, goodRankTranspose, NULL,
      dec->isTernary ? &violatorSubmatrix : NULL) );

    if (violatorSubmatrix)
    {
      CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");

      CMR_CALL( CMRmatroiddecUpdateSubmatrix(cmr, dec, violatorSubmatrix, CMR_MATROID_DEC_TYPE_DETERMINANT) );
      assert(dec->type != CMR_MATROID_DEC_TYPE_DETERMINANT);

      CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );

      goto cleanup;
    }

    CMR_CALL( CMRsepaGetRepresentatives(cmr, separation, extraRows, extraColumns) );
  }
  else
  {
    CMRdbgMsg(10, "-> No need to pivot.\n");
    goodRankMatrix = dec->matrix;
    goodRankTranspose = dec->transpose;
  }

  /* Now the ranks are good for goodRankMatrix and goodRankTranspose, aligned with separation. */

  CMR_CHRMAT* pivotedMatrix = NULL;
  CMR_CHRMAT* pivotedTranspose = NULL;
  if (task->params->threeSumPivotChildren)
  {
    assert(false);
  }
  else
  {
    pivotedMatrix = goodRankMatrix;
    pivotedTranspose = goodRankTranspose;
  }

  if (numPivots > 0)
  {
    CMR_CALL( CMRmatroiddecUpdatePivots(cmr, dec, numPivots, pivotRows, pivotColumns, pivotedMatrix, pivotedTranspose) );
    dec = dec->children[0];
  }

  /* Initialize the 3-sum node. */
  CMR_CALL( CMRmatroiddecUpdateThreeSumInit(cmr, dec) );
  CMR_CALL( CMRallocStackArray(cmr, &rowsToChild, dec->numRows) );
  CMR_CALL( CMRallocStackArray(cmr, &columnsToChild, dec->numColumns) );
  size_t numChildBaseRows;
  size_t numChildBaseColumns;

  if (separation->type == CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS)
  {
    CMR_CALL( CMRsepaGetProjection(separation, 0, rowsToChild, columnsToChild, &numChildBaseRows,
      &numChildBaseColumns) );

    if (task->params->threeSumStrategy & CMR_MATROID_DEC_THREESUM_FLAG_FIRST_TALL)
    {
      /* First child is tall. */

      assert(false);
    }
    else
    {
      CMR_CALL( CMRmatroiddecUpdateThreeSumCreateWideFirstChild(cmr, dec, separation, rowsToChild, columnsToChild,
        numChildBaseRows, numChildBaseColumns, extraRows[0][0], extraColumns[0][0], extraColumns[0][0],
        dec->isTernary ? 0 : 1) );
    }

    CMR_CALL( CMRsepaGetProjection(separation, 1, rowsToChild, columnsToChild, &numChildBaseRows,
      &numChildBaseColumns) );

    if (task->params->threeSumStrategy & CMR_MATROID_DEC_THREESUM_FLAG_SECOND_TALL)
    {
      /* Second child is tall. */

      assert(false);
    }
    else
    {
      CMR_CALL( CMRmatroiddecUpdateThreeSumCreateWideSecondChild(cmr, dec, separation, rowsToChild, columnsToChild,
        numChildBaseRows, numChildBaseColumns, extraRows[1][0], extraColumns[1][0], extraColumns[1][0],
        dec->isTernary ? 0 : 1) );
    }
  }
  else
  {
    assert(separation->type == CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK);

    if (task->params->threeSumStrategy & CMR_MATROID_DEC_THREESUM_FLAG_FIRST_ALLREPR)
    {
      /* First child is all-repr`. */

      assert(false);
    }
    else
    {
      /* First child is mixed. */

      assert(false);
    }
  }

cleanup:

  if (rowsToChild)
  {
    CMR_CALL( CMRfreeStackArray(cmr, &columnsToChild) );
    CMR_CALL( CMRfreeStackArray(cmr, &rowsToChild) );
  }

  /* Free the good-rank matrices if they were created explicitly but are not needed due to further pivots. */

  if (pivotRow != SIZE_MAX && goodRankMatrix != pivotedMatrix)
  {
    CMRchrmatFree(cmr, &goodRankTranspose);
    CMRchrmatFree(cmr, &goodRankMatrix);
  }

  if (pivotColumns)
  {
    CMR_CALL( CMRfreeStackArray(cmr, &pivotColumns) );
    CMR_CALL( CMRfreeStackArray(cmr, &pivotRows) );
  }

  if (dec->type == CMR_MATROID_DEC_TYPE_THREE_SUM)
  {
    DecompositionTask* childTasks[2] = { task, NULL };
    CMR_CALL( CMRregularityTaskCreateRoot(cmr, dec->children[1], &childTasks[1], task->params, task->stats,
      task->startClock, task->timeLimit) );

    childTasks[0]->dec = dec->children[0];
    dec->children[0]->testedSeriesParallel = false;
    dec->children[1]->testedSeriesParallel = false;

    /* Add both child tasks to the list. */
    childTasks[0]->next = childTasks[1];
    childTasks[1]->next = *punprocessed;
    *punprocessed = childTasks[0];
  }
  else
  {
    assert(!"UNIMPLEMENTED"); /* TODO: Tasks? */
  }

  return CMR_OKAY;
}
