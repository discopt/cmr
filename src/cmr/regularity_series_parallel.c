// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/series_parallel.h>

#include "env_internal.h"
#include "seymour_internal.h"

#include <time.h>

CMR_ERROR CMRregularityDecomposeSeriesParallel(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

    CMR_SEYMOUR_NODE* dec = task->node;
    CMR_SEYMOUR_NODE* decReduced = dec;
  assert(dec);

#if defined(CMR_DEBUG)
  CMRdbgMsg(6, "Testing the following matrix for series-parallel reductions:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stndout, '0', true) );
#endif /* CMR_DEBUG */

  double remainingTime = task->timeLimit - (clock() - task->startClock) * 1.0 / CLOCKS_PER_SEC;

  bool isSeriesParallel = true;
  CMR_SP_REDUCTION* reductions = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &reductions, dec->matrix->numRows + dec->matrix->numColumns) );
  size_t numReductions;

  CMR_SUBMAT* violatorSubmatrix = NULL;
  CMR_SUBMAT* reducedSubmatrix = NULL;
  CMR_SEPA* separation = NULL;
  if (dec->isTernary)
  {
    CMR_CALL( CMRspDecomposeTernary(cmr, dec->matrix, &isSeriesParallel, reductions,
      task->params->seriesParallel ? SIZE_MAX : 1, &numReductions, &reducedSubmatrix, &violatorSubmatrix, &separation,
      task->stats ? &task->stats->seriesParallel : NULL, remainingTime) );

    assert(violatorSubmatrix || separation || (numReductions == dec->numRows + dec->numColumns));
  }
  else
  {
    CMR_CALL( CMRspDecomposeBinary(cmr, dec->matrix, &isSeriesParallel, reductions,
      task->params->seriesParallel ? SIZE_MAX : 1,  &numReductions, &reducedSubmatrix, &violatorSubmatrix, &separation,
      task->stats ? &task->stats->seriesParallel : NULL, remainingTime) );
  }

  /* Did we find a 2-by-2 submatrix? If yes, is has determinant -2 or +2! */
  if (violatorSubmatrix && violatorSubmatrix->numRows == 2)
  {
    CMRdbgMsg(8, "-> 2x2 violator.\n");
    assert(dec->isTernary);

    CMR_CALL( CMRseymourUpdateViolator(cmr, dec, violatorSubmatrix) );

    /* Task is done. */
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );

    CMR_CALL( CMRfreeStackArray(cmr, &reductions) );
    CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );
    queue->foundIrregularity = true;

    return CMR_OKAY;
  }

  /* Modify the decomposition to reflect the SP reductions. */
  if (numReductions > 0)
  {
    if (task->params->seriesParallel)
    {
      CMRdbgMsg(8, "-> %zu series-parallel reductions.\n", numReductions);

      dec->type = CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL;
      dec->numSeriesParallelReductions = numReductions;
      CMR_CALL( CMRduplicateBlockArray(cmr, &dec->seriesParallelReductions, numReductions, reductions) );

      if (!isSeriesParallel)
      {
        CMRdbgMsg(8, "-> replacing current node by submatrix child.\n", numReductions);

        CMR_CALL( CMRseymourUpdateSeriesParallel(cmr, dec, reducedSubmatrix) );
        assert(dec->numChildren == 1);
        decReduced = dec->children[0];
        task->node = decReduced;
        task->node->testedSeriesParallel = true;
      }
    }
    else
    {
      /* We carry out the first SP reduction as a 2-separation. */
      assert(numReductions == SIZE_MAX);

      char buffer[16];
      CMRdbgMsg(8, "-> applying reduction (%s,%s) as a 2-separation.\n", CMRelementString(reductions[0].element, NULL),
        CMRelementString(reductions[0].mate, buffer));

      assert(!separation);
      size_t parentNumRows = dec->matrix->numRows;
      size_t parentNumColumns = dec->matrix->numColumns;
      CMR_CALL( CMRsepaCreate(cmr, parentNumRows, parentNumColumns, &separation) );

      bool spIsRow = CMRspIsRow(reductions[0]);
      bool mateIsRow = CMRelementIsRow(reductions[0].mate);
      for (size_t r = 0; r < parentNumRows; ++r)
        separation->rowsFlags[r] = mateIsRow ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
      for (size_t c = 0; c < parentNumColumns; ++c)
        separation->columnsFlags[c] = mateIsRow ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;

      if (spIsRow && mateIsRow)
      {
        /* Duplicate row: big part is top left. */
        separation->rowsFlags[CMRelementToRowIndex(reductions[0].element)] = CMR_SEPA_SECOND;
        separation->rowsFlags[CMRelementToRowIndex(reductions[0].mate)] = CMR_SEPA_SECOND;
      }
      if (spIsRow && !mateIsRow)
      {
        /* Unit row: big part is bottom right. */
        separation->rowsFlags[CMRelementToRowIndex(reductions[0].element)] = CMR_SEPA_FIRST;
        separation->columnsFlags[CMRelementToColumnIndex(reductions[0].mate)] = CMR_SEPA_FIRST;
      }
      if (!spIsRow && mateIsRow)
      {
        /* Unit column: big part is top left. */
        separation->columnsFlags[CMRelementToColumnIndex(reductions[0].element)] = CMR_SEPA_SECOND;
        separation->rowsFlags[CMRelementToRowIndex(reductions[0].mate)] = CMR_SEPA_SECOND;
      }
      if (!spIsRow && !mateIsRow)
      {
        /* Duplicate column: big part is bottom right. */
        separation->columnsFlags[CMRelementToColumnIndex(reductions[0].element)] = CMR_SEPA_FIRST;
        separation->columnsFlags[CMRelementToColumnIndex(reductions[0].mate)] = CMR_SEPA_FIRST;
      }

      CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, dec->matrix, dec->transpose, NULL, NULL) );
      assert(separation->type == CMR_SEPA_TYPE_TWO);
    }
  }
  else
  {
    CMRdbgMsg(8, "-> no series-parallel reductions applicable.\n");
    task->node->testedSeriesParallel = true;
  }

  /* Note: decReduced is the node for the SP-reduced submatrix (equal to dec if there are no SP-reductions). */

  /* Modify the decomposition for the 2-separation. */
  if (separation)
  {
    CMRdbgMsg(8, "-> 2-separation found.\n");
    CMR_CALL( CMRseymourUpdateTwoSum(cmr, decReduced, separation) );
    CMR_CALL( CMRsepaFree(cmr, &separation) );

    DecompositionTask* childTasks[2] = { task, NULL };
    CMR_CALL( CMRregularityTaskCreateRoot(cmr, decReduced->children[1], &childTasks[1], task->params, task->stats,
      task->startClock, task->timeLimit) );

    childTasks[0]->node = decReduced->children[0];

    /* We mark the children as not SP-tested because we did not find a W_3 minor. */
    decReduced->children[0]->testedSeriesParallel = false;
    decReduced->children[1]->testedSeriesParallel = false;

    /* Add both child tasks to the list. */
    CMRregularityQueueAdd(queue, childTasks[0]);
    CMRregularityQueueAdd(queue, childTasks[1]);
  }
  else
  {
    /* We can use the violator submatrix to start the sequence of nested minors. */
    assert(violatorSubmatrix);
    CMRdbgMsg(8, "-> wheel-submatrix representing W_%zu was found.\n", violatorSubmatrix->numRows);

    if (numReductions > 0)
    {
      /* violatorSubmatrix refers to the given matrix, but we need it in terms of the reduced matrix. */
      size_t* originalToReduced = NULL;
      CMR_CALL( CMRallocStackArray(cmr, &originalToReduced,
        dec->numRows > dec->numColumns ? dec->numRows : dec->numColumns) );

      /* Rows */
#if !defined(NDEBUG)
      for (size_t row = 0; row < dec->numRows; ++row)
        originalToReduced[row] = SIZE_MAX;
#endif /* !NDEBUG */

      for (size_t reducedRow = 0; reducedRow < reducedSubmatrix->numRows; ++reducedRow)
        originalToReduced[reducedSubmatrix->rows[reducedRow]] = reducedRow;
      for (size_t violatorRow = 0; violatorRow < violatorSubmatrix->numRows; ++violatorRow)
      {
        assert(originalToReduced[violatorSubmatrix->rows[violatorRow]] != SIZE_MAX);
        violatorSubmatrix->rows[violatorRow] = originalToReduced[violatorSubmatrix->rows[violatorRow]];
      }

      /* Columns */
#if !defined(NDEBUG)
      for (size_t column = 0; column < dec->numColumns; ++column)
        originalToReduced[column] = SIZE_MAX;
#endif /* !NDEBUG */

      for (size_t reducedColumn = 0; reducedColumn < reducedSubmatrix->numColumns; ++reducedColumn)
        originalToReduced[reducedSubmatrix->columns[reducedColumn]] = reducedColumn;
      for (size_t violatorColumn = 0; violatorColumn < violatorSubmatrix->numColumns; ++violatorColumn)
      {
        assert(originalToReduced[violatorSubmatrix->columns[violatorColumn]] != SIZE_MAX);
        violatorSubmatrix->columns[violatorColumn] = originalToReduced[violatorSubmatrix->columns[violatorColumn]];
      }

      CMR_CALL( CMRfreeStackArray(cmr, &originalToReduced) );
    }

    remainingTime = task->timeLimit - (clock() - task->startClock) * 1.0 / CLOCKS_PER_SEC;
    task->node = decReduced;
    decReduced->testedTwoConnected = true;
    decReduced->graphicness = dec->graphicness;
    decReduced->cographicness = dec->cographicness;
    CMR_CALL( CMRregularityInitNestedMinorSequence(cmr, task, violatorSubmatrix) );

    /* Just mark the current task as SP-reduced and add it back to the list of unprocessed tasks. */
    CMRregularityQueueAdd(queue, task);
  }

  CMR_CALL( CMRfreeStackArray(cmr, &reductions) );
  CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );
  CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );

  return CMR_OKAY;
}
