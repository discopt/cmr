// #define CMR_DEBUG /* Uncomment to debug this file. */
// #define CMR_DEBUG_MATRICES /* Uncomment to print matrices. */

#include "env_internal.h"
#include "seymour_internal.h"

#include <time.h>

CMR_ERROR CMRregularitySimpleSearchThreeSeparation(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMRdbgMsg(6, "Starting search for simple 3-separations for %zux%zu-matrix.\n", task->node->matrix->numRows,
    task->node->matrix->numColumns);
#if defined(CMR_DEBUG_MATRICES)
  CMR_CALL( CMRchrmatPrintDense(cmr, task->node->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG_MATRICES */

  clock_t startTime = clock();
  task->stats->simpleThreeSeparationsCount++;

  /* We need the transpose matrix. */
  CMR_SEYMOUR_NODE* node = task->node;
  CMR_CHRMAT* matrix = node->matrix;
  if (!node->transpose)
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &node->transpose) );
  CMR_CHRMAT* transpose = node->transpose;

  for (size_t row1 = 0; row1 < matrix->numRows; ++row1)
  {
    size_t first = matrix->rowSlice[row1];
    size_t beyond = matrix->rowSlice[row1 + 1];

    if (beyond != first + 2)
      continue;

    size_t* rowColumns = &matrix->entryColumns[first];
    for (size_t c = 0; c < 2; ++c)
    {
      size_t column1 = rowColumns[c];
      size_t numRowColumns = transpose->rowSlice[column1 + 1] - transpose->rowSlice[column1];
      if (numRowColumns != 2)
        continue;

      /* We have found a nonzero at (r,c) such that row r has 2 nonzeros and column c has 2 nonzeros. */
      size_t column2 = rowColumns[1-c];
      size_t row2 = SIZE_MAX;
      for (size_t e = transpose->rowSlice[column1]; e < transpose->rowSlice[column1 + 1]; ++e)
      {
        size_t row = transpose->entryColumns[e];
        if (row != row1)
          row2 = row;
        else
          assert(row == row1);
      }

      CMRdbgMsg(8, "Found a simple 3-separation for rows {r%zu,r%zu} and columns {c%zu,c%zu}.\n", row1+1, row2+1,
        column1+1, column2+1);

      task->stats->simpleThreeSeparationsSuccess++;
      CMR_SEPA* separation = NULL;
      CMR_CALL( CMRsepaCreate(cmr, matrix->numRows, matrix->numColumns, &separation) );

      separation->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;
      for (size_t row = 0; row < matrix->numRows; ++row)
        separation->rowsFlags[row] = (row == row1 || row == row2) ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;
      for (size_t column = 0; column < matrix->numColumns; ++column)
        separation->columnsFlags[column] = (column == column1 || column == column2) ? CMR_SEPA_FIRST : CMR_SEPA_SECOND;

      CMR_SUBMAT* violatorSubmatrix = NULL;
      CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, matrix, transpose, NULL,
        node->isTernary ? &violatorSubmatrix : NULL) );

      if (violatorSubmatrix)
      {
        CMRdbgMsg(8, "-> 2x2 submatrix with bad determinant.\n");

        CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );

        CMR_CALL( CMRregularityTaskFree(cmr, &task) );
      }
      else
      {
        CMR_CALL( CMRregularityDecomposeThreeSum(cmr, task, queue, separation) );
      }

      CMR_CALL( CMRsepaFree(cmr, &separation) );

      goto cleanup;
    }
  }

  /* We have found nothing, so we note this and re-add the task to the queue. */
  task->node->testedSimpleThreeSeparations = true;
  CMRregularityQueueAdd(queue, task);


cleanup:

  task->stats->simpleThreeSeparationsTime += (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;

  return CMR_OKAY;
}
