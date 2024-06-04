// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "matroid_internal.h"
#include "regularity_internal.h"

typedef enum
{
  ELEMENT_TYPE_NONE,    /**< Row/column does not belong to the submatrix. */
  ELEMENT_TYPE_NORMAL,  /**< Row/column belongs to the submatrix. */
  ELEMENT_TYPE_SOURCE,  /**< Row/column is a source. */
  ELEMENT_TYPE_TARGET   /**< Row/column is a target. */
} ElementType;

typedef struct
{
  size_t predecessor; /**< Predecessor node (row or column index). */
  int8_t status;      /**< 0: unknown, 1: in queue, 2: processed. */
  int8_t edgeValue;   /**< Matrix entry of the edge to the predecessor. */
} GraphNode;

/**
 * \brief Finds a shortest path from any source to any target element.
 */

static
CMR_ERROR findSubmatrixCycle(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Matrix. */
  CMR_CHRMAT* transpose,    /**< Transpose of \p matrix. */
  ElementType* rowTypes,    /**< Array with rows' types. */
  ElementType* columnTypes, /**< Array with columns' types. */
  CMR_ELEMENT* ppathSource, /**< Pointer for storing the source row/column. */
  CMR_ELEMENT* ppathTarget, /**< Pointer for storing the target row/column. */
  int* pentrySum            /**< Pointer for storing the sum of the path's entries. */
)
{
  assert(cmr);
  assert(matrix);
  assert(transpose);
  assert(rowTypes);
  assert(columnTypes);
  assert(ppathSource);
  assert(ppathTarget);
  assert(pentrySum);

  size_t numRows = matrix->numRows;
  size_t numColumns = matrix->numColumns;

  CMRdbgMsg(12, "findSubmatrixCycle.\n");

  /* Initialize row / column arrays. */

  GraphNode* rowData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowData, numRows) );
  GraphNode* columnData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnData, numColumns) );

  /* Initialize queue. */
  CMR_ELEMENT* queue = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &queue, numRows + numColumns) );
  size_t firstQueued = 0;
  size_t beyondQueued = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    rowData[row].predecessor = SIZE_MAX;
    if (rowTypes[row] == ELEMENT_TYPE_SOURCE)
    {
      queue[beyondQueued++] = CMRrowToElement(row);
      rowData[row].status = 1;
    }
    else
      rowData[row].status = 0;
    CMRdbgMsg(14, "r%zu has status %d and type %d.\n", row+1, rowData[row].status, rowTypes[row]);
  }
  for (size_t column = 0; column < numColumns; ++column)
  {
    columnData[column].predecessor = SIZE_MAX;
    if (columnTypes[column] == ELEMENT_TYPE_SOURCE)
    {
      queue[beyondQueued++] = CMRcolumnToElement(column);
      columnData[column].status = 1;
    }
    else
      columnData[column].status = 0;
    CMRdbgMsg(14, "c%zu has status %d and type %d.\n", column+1, columnData[column].status, columnTypes[column]);
  }

  /* Start BFS. */
  bool done = false;
  while (!done && firstQueued < beyondQueued)
  {
    CMR_ELEMENT current = queue[firstQueued++];

    CMRdbgMsg(12, "findSubmatrixCycle processes %s.\n", CMRelementString(current, NULL));

    if (CMRelementIsRow(current))
    {
      size_t row = CMRelementToRowIndex(current);
      rowData[row].status = 2;

      size_t first = matrix->rowSlice[row];
      size_t beyond = matrix->rowSlice[row + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        size_t column = matrix->entryColumns[e];
        ElementType type = columnTypes[column];
        if (type == ELEMENT_TYPE_NONE || columnData[column].status)
          continue;

        columnData[column].predecessor = row;
        columnData[column].edgeValue = matrix->entryValues[e];

        if (type == ELEMENT_TYPE_TARGET)
        {
          done = true;
          *ppathTarget = CMRcolumnToElement(column);
          break;
        }
        columnData[column].status = 1;
        queue[beyondQueued++] = CMRcolumnToElement(column);
      }
    }
    else
    {
      size_t column = CMRelementToColumnIndex(current);
      columnData[column].status = 2;

      size_t first = transpose->rowSlice[column];
      size_t beyond = transpose->rowSlice[column + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        size_t row = transpose->entryColumns[e];
        ElementType type = rowTypes[row];
        if (type == ELEMENT_TYPE_NONE || rowData[row].status)
          continue;

        rowData[row].predecessor = column;
        rowData[row].edgeValue = transpose->entryValues[e];

        if (type == ELEMENT_TYPE_TARGET)
        {
          done = true;
          *ppathTarget = CMRrowToElement(row);
          break;
        }
        rowData[row].status = 1;
        queue[beyondQueued++] = CMRrowToElement(row);
      }
    }
  }

  /* Trace back the path. */
  assert(done);
  CMR_ELEMENT current = *ppathTarget;
  *pentrySum = 0;
  while (true)
  {
    if (CMRelementIsRow(current))
    {
      size_t row = CMRelementToRowIndex(current);
      size_t column = rowData[row].predecessor;
      CMRdbgMsg(12, "Predecessor of r%zu is c%zu.\n", row + 1, column + 1);
      if (column < SIZE_MAX)
      {
        (*pentrySum) += rowData[row].edgeValue;
        current = CMRcolumnToElement(column);
      }
      else
      {
        *ppathSource = current;
        break;
      }
    }
    else
    {
      size_t column = CMRelementToColumnIndex(current);
      size_t row = columnData[column].predecessor;
      CMRdbgMsg(12, "Predecessor of c%zu is r%zu.\n", column + 1, row + 1);
      if (row < SIZE_MAX)
      {
        (*pentrySum) += columnData[column].edgeValue;
        current = CMRrowToElement(row);
      }
      else
      {
        *ppathSource = current;
        break;
      }
    }
  }

  char buffer[16];
  CMRdbgMsg(10, "findSubmatrixCycle found a path from %s to %s with entry sum %d.\n",
    CMRelementString(*ppathSource, NULL), CMRelementString(*ppathTarget, buffer), *pentrySum);

  /* Cleanup */

  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnData) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowData) );

  return CMR_OKAY;
}

CMR_ERROR CMRregularityDecomposeThreeSum(
  CMR* cmr,
  DecompositionTask* task,
  DecompositionQueue* queue,
  CMR_SEPA* separation
)
{
  assert(cmr);
  assert(task);
  assert(queue);
  assert(separation);

    CMR_SEYMOUR_NODE* dec = task->dec;
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
    (task->params->threeSumStrategy & CMR_SEYMOUR_NODE_THREESUM_FLAG_CONCENTRATED_RANK))
  {
    /* We have distributed ranks but requested concentrated rank; we find a top-right nonzero. */
    pivotRow = extraRows[1][0];
    assert(pivotRow != SIZE_MAX);

    size_t first = dec->matrix->rowSlice[pivotRow];
    size_t beyond = dec->matrix->rowSlice[pivotRow + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = dec->matrix->entryColumns[e];
      int flags = separation->columnsFlags[column];
      if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
      {
        pivotColumn = column;
        break;
      }
    }
    assert(pivotColumn != SIZE_MAX);
  }

  if (separation->type == CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK &&
    (task->params->threeSumStrategy & CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS))
  {
    /* We have concentrated rank but requested distributed ranks; we find a bottom-left nonzero. */
    pivotRow = extraRows[0][0];
    assert(pivotRow != SIZE_MAX);

    size_t first = dec->matrix->rowSlice[pivotRow];
    size_t beyond = dec->matrix->rowSlice[pivotRow + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t column = dec->matrix->entryColumns[e];
      int flags = separation->columnsFlags[column];
      if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
      {
        pivotColumn = column;
        break;
      }
    }
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

      CMR_CALL( CMRmatroiddecUpdateSubmatrix(cmr, dec, violatorSubmatrix, CMR_SEYMOUR_NODE_TYPE_DETERMINANT) );
      assert(dec->type != CMR_SEYMOUR_NODE_TYPE_DETERMINANT);

      CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );

      /* TODO: Add as unittest. */
      queue->foundIrregularity = true;

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

  ElementType* rowTypes = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowTypes, dec->numRows) );
  ElementType* columnTypes = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnTypes, dec->numColumns) );

  if (separation->type == CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS)
  {
    CMR_CALL( CMRsepaGetProjection(separation, 0, rowsToChild, columnsToChild, &numChildBaseRows,
      &numChildBaseColumns) );

    if (task->params->threeSumStrategy & CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_TALL)
    {
      /* First child is tall. */

      assert(false);
    }
    else
    {
      /* Prepare shortest-path search in bottom-right submatrix. */
      CMR_ELEMENT sourceRowElement;
      CMR_ELEMENT targetColumnElement;
      int extraEntry;
      if (dec->isTernary)
      {
        for (size_t row = 0; row < dec->numRows; ++row)
        {
          CMR_SEPA_FLAGS flags = separation->rowsFlags[row];
          if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
            rowTypes[row] = ELEMENT_TYPE_NONE;
          else if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
            rowTypes[row] = ELEMENT_TYPE_SOURCE;
          else
            rowTypes[row] = ELEMENT_TYPE_NORMAL;
        }
        for (size_t column = 0; column < dec->numColumns; ++column)
        {
          CMR_SEPA_FLAGS flags = separation->columnsFlags[column];
          if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
            columnTypes[column] = ELEMENT_TYPE_NONE;
          else if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
            columnTypes[column] = ELEMENT_TYPE_TARGET;
          else
            columnTypes[column] = ELEMENT_TYPE_NORMAL;
        }
        int sumEntries;
        CMR_CALL( findSubmatrixCycle(cmr, dec->matrix, dec->transpose, rowTypes, columnTypes, &sourceRowElement,
          &targetColumnElement, &sumEntries) );
        sumEntries = ((sumEntries % 4) + 4) % 4;
        assert(sumEntries == 1 || sumEntries == 3);
        extraEntry = (sumEntries == 1) ? +1 : -1;
      }
      else
      {
        sourceRowElement = CMRrowToElement(extraRows[0][0]);
        targetColumnElement = CMRcolumnToElement(extraColumns[0][0]);
        extraEntry = 1;

        /* TODO: Is there always a path between these particular representatives that does not use any other
         * representative? */
      }

      assert( CMRelementIsRow(sourceRowElement) );
      assert( CMRelementIsColumn(targetColumnElement) );

      CMR_CALL( CMRmatroiddecUpdateThreeSumCreateWideFirstChild(cmr, dec, separation, rowsToChild, columnsToChild,
        numChildBaseRows, numChildBaseColumns, CMRelementToRowIndex(sourceRowElement),
        CMRelementToColumnIndex(targetColumnElement), CMRelementToColumnIndex(targetColumnElement), extraEntry) );
    }

    CMR_CALL( CMRsepaGetProjection(separation, 1, rowsToChild, columnsToChild, &numChildBaseRows,
      &numChildBaseColumns) );

    if (task->params->threeSumStrategy & CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_TALL)
    {
      /* Second child is tall. */

      assert(false);
    }
    else
    {
      /* Prepare shortest-path search in top-left submatrix. */
      CMR_ELEMENT sourceRowElement;
      CMR_ELEMENT targetColumnElement;
      int extraEntry;
      if (dec->isTernary)
      {
        for (size_t row = 0; row < dec->numRows; ++row)
        {
          CMR_SEPA_FLAGS flags = separation->rowsFlags[row];
          if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
            rowTypes[row] = ELEMENT_TYPE_NONE;
          else if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
            rowTypes[row] = ELEMENT_TYPE_SOURCE;
          else
            rowTypes[row] = ELEMENT_TYPE_NORMAL;
        }
        for (size_t column = 0; column < dec->numColumns; ++column)
        {
          CMR_SEPA_FLAGS flags = separation->columnsFlags[column];
          if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
            columnTypes[column] = ELEMENT_TYPE_NONE;
          else if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
            columnTypes[column] = ELEMENT_TYPE_TARGET;
          else
            columnTypes[column] = ELEMENT_TYPE_NORMAL;
        }
        int sumEntries;
        CMR_CALL( findSubmatrixCycle(cmr, dec->matrix, dec->transpose, rowTypes, columnTypes, &sourceRowElement,
          &targetColumnElement, &sumEntries) );
        sumEntries = ((sumEntries % 4) + 4) % 4;
        assert(sumEntries == 1 || sumEntries == 3);
        extraEntry = (sumEntries == 1) ? +1 : -1;
      }
      else
      {
        sourceRowElement = CMRrowToElement(extraRows[1][0]);
        targetColumnElement = CMRcolumnToElement(extraColumns[1][0]);
        extraEntry = 1;

        /* TODO: Is there always a path between these particular representatives that does not use any other
         * representative? */
      }

      assert( CMRelementIsRow(sourceRowElement) );
      assert( CMRelementIsColumn(targetColumnElement) );

      CMR_CALL( CMRmatroiddecUpdateThreeSumCreateWideSecondChild(cmr, dec, separation, rowsToChild, columnsToChild,
        numChildBaseRows, numChildBaseColumns, CMRelementToRowIndex(sourceRowElement),
        CMRelementToColumnIndex(targetColumnElement), CMRelementToColumnIndex(targetColumnElement), extraEntry) );
    }

    dec->threesumFlags = CMR_SEYMOUR_NODE_THREESUM_FLAG_DISTRIBUTED_RANKS;
  }
  else
  {
    assert(separation->type == CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK);

    CMR_CALL( CMRsepaGetProjection(separation, 0, rowsToChild, columnsToChild, &numChildBaseRows,
      &numChildBaseColumns) );

    if (task->params->threeSumStrategy & CMR_SEYMOUR_NODE_THREESUM_FLAG_FIRST_ALLREPR)
    {
      /* First child is all-repr`. */

      assert(false);
    }
    else
    {
      /* First child is mixed. */

      /* Prepare shortest-path search in bottom-right submatrix. */
      CMR_ELEMENT sourceRowElement;
      CMR_ELEMENT targetRowElement;
      int extraEntry;

      for (size_t row = 0; row < dec->numRows; ++row)
      {
        CMR_SEPA_FLAGS flags = separation->rowsFlags[row];
        if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          rowTypes[row] = ELEMENT_TYPE_NONE;
        else if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
          rowTypes[row] = ELEMENT_TYPE_SOURCE;
        else if ((flags & CMR_SEPA_MASK_EXTRA) & CMR_SEPA_FLAG_RANK2)
          rowTypes[row] = ELEMENT_TYPE_TARGET;
        else
          rowTypes[row] = ELEMENT_TYPE_NORMAL;
      }
      for (size_t column = 0; column < dec->numColumns; ++column)
      {
        CMR_SEPA_FLAGS flags = separation->columnsFlags[column];
        if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          columnTypes[column] = ELEMENT_TYPE_NONE;
        else
          columnTypes[column] = ELEMENT_TYPE_NORMAL;
      }
      int sumEntries;
      CMR_CALL( findSubmatrixCycle(cmr, dec->matrix, dec->transpose, rowTypes, columnTypes, &sourceRowElement,
        &targetRowElement, &sumEntries) );
      sumEntries = ((sumEntries % 4) + 4) % 4;
      assert(sumEntries == 0 || sumEntries == 2);
      extraEntry = (sumEntries == 2) ? +1 : -1;

      assert( CMRelementIsRow(sourceRowElement) );
      assert( CMRelementIsRow(targetRowElement) );

      CMR_CALL( CMRmatroiddecUpdateThreeSumCreateMixedFirstChild(cmr, dec, separation, rowsToChild, columnsToChild,
        numChildBaseRows, numChildBaseColumns, CMRelementToRowIndex(sourceRowElement),
        CMRelementToRowIndex(targetRowElement), extraEntry) );
    }

    CMR_CALL( CMRsepaGetProjection(separation, 1, rowsToChild, columnsToChild, &numChildBaseRows,
      &numChildBaseColumns) );

    if (task->params->threeSumStrategy & CMR_SEYMOUR_NODE_THREESUM_FLAG_SECOND_ALLREPR)
    {
      /* Second child is all-repr. */

      assert(false);
    }
    else
    {
      /* Prepare shortest-path search in top-left submatrix. */
      CMR_ELEMENT sourceColumnElement;
      CMR_ELEMENT targetColumnElement;
      int extraEntry;

      for (size_t row = 0; row < dec->numRows; ++row)
      {
        CMR_SEPA_FLAGS flags = separation->rowsFlags[row];
        if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          rowTypes[row] = ELEMENT_TYPE_NONE;
        else
          rowTypes[row] = ELEMENT_TYPE_NORMAL;

      }
      for (size_t column = 0; column < dec->numColumns; ++column)
      {
        CMR_SEPA_FLAGS flags = separation->columnsFlags[column];
        if ((flags & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          columnTypes[column] = ELEMENT_TYPE_NONE;
        else if ((flags & CMR_SEPA_MASK_EXTRA) == CMR_SEPA_FLAG_RANK1)
          columnTypes[column] = ELEMENT_TYPE_SOURCE;
        else if ((flags & CMR_SEPA_MASK_EXTRA) & CMR_SEPA_FLAG_RANK2)
          columnTypes[column] = ELEMENT_TYPE_TARGET;
        else
          columnTypes[column] = ELEMENT_TYPE_NORMAL;
      }
      int sumEntries;
      CMR_CALL( findSubmatrixCycle(cmr, dec->matrix, dec->transpose, rowTypes, columnTypes, &sourceColumnElement,
        &targetColumnElement, &sumEntries) );
      sumEntries = ((sumEntries % 4) + 4) % 4;
      assert(sumEntries == 0 || sumEntries == 2);
      extraEntry = (sumEntries == 2) ? +1 : -1;

      assert( CMRelementIsColumn(sourceColumnElement) );
      assert( CMRelementIsColumn(targetColumnElement) );

      CMR_CALL( CMRmatroiddecUpdateThreeSumCreateMixedSecondChild(cmr, dec, separation, rowsToChild, columnsToChild,
        numChildBaseRows, numChildBaseColumns, CMRelementToColumnIndex(sourceColumnElement),
        CMRelementToColumnIndex(targetColumnElement), extraEntry) );
    }

    dec->threesumFlags = CMR_SEYMOUR_NODE_THREESUM_FLAG_CONCENTRATED_RANK;
  }

cleanup:

  if (rowsToChild)
  {
    CMR_CALL( CMRfreeStackArray(cmr, &columnTypes) );
    CMR_CALL( CMRfreeStackArray(cmr, &rowTypes) );
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

  if (dec->type == CMR_SEYMOUR_NODE_TYPE_THREE_SUM)
  {
    DecompositionTask* childTasks[2] = { task, NULL };
    CMR_CALL( CMRregularityTaskCreateRoot(cmr, dec->children[1], &childTasks[1], task->params, task->stats,
      task->startClock, task->timeLimit) );

    childTasks[0]->dec = dec->children[0];
    dec->children[0]->testedSeriesParallel = false;
    dec->children[1]->testedSeriesParallel = false;

    /* Add both child tasks to the list. */
    CMRregularityQueueAdd(queue, childTasks[0]);
    CMRregularityQueueAdd(queue, childTasks[1]);
  }

  return CMR_OKAY;
}
