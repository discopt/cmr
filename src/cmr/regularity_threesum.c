// #define CMR_DEBUG /** Uncomment to debug this file. */

#include <cmr/regular.h>

#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"
#include "seymour_internal.h"

typedef enum
{
  ELEMENT_TYPE_NONE,  /**< Row/column does not belong to the submatrix. */
  ELEMENT_TYPE_1,     /**< Row/column is a source/target of type 1. */
  ELEMENT_TYPE_2,     /**< Row/column is a source/target of type 2. */
  ELEMENT_TYPE_3,     /**< Row/column is a source/target of type 3. */
  ELEMENT_TYPE_NORMAL /**< Row/column belongs to the submatrix. */
} ElementType;

typedef struct
{
  size_t predecessor; /**< Predecessor node (row or column index). */
  int8_t status;      /**< 0: unknown, 1: in queue, 2: processed. */
  int8_t edgeValue;   /**< Matrix entry of the edge to the predecessor. */
} GraphNode;

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

  CMR_SEYMOUR_NODE* node = task->node;
  assert(node);

#if defined(CMR_DEBUG)
  CMRdbgMsg(8, "Processing 3-separation to produce a 3-sum for the following matrix:\n");
  CMRchrmatPrintDense(cmr, node->matrix, stdout, '0', true);
  for (size_t row = 0; row < node->numRows; ++row)
  {
    int part = ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? 0 : 1;
    int rank1Flag = (separation->rowsFlags[row] & CMR_SEPA_FLAG_RANK1) ? 1 : 0;
    int rank2Flag = (separation->rowsFlags[row] & CMR_SEPA_FLAG_RANK2) ? 1 : 0;
    CMRdbgMsg(10, "Row r%zu belongs to part %d; flags = %d/%d\n", row+1, part, rank1Flag, rank2Flag);
  }
  for (size_t column = 0; column < node->numColumns; ++column)
  {
    int part = ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST) ? 0 : 1;
    int rank1Flag = (separation->columnsFlags[column] & CMR_SEPA_FLAG_RANK1) ? 1 : 0;
    int rank2Flag = (separation->columnsFlags[column] & CMR_SEPA_FLAG_RANK2) ? 1 : 0;
    CMRdbgMsg(10, "Column c%zu belongs to part %d; flags = %d/%d\n", column+1, part, rank1Flag, rank2Flag);
  }
#endif /* CMR_DEBUG */

  if (((separation->type == CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS)
    && ((task->params->decomposeStrategy & CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_MASK)
    == CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_PIVOT)) || ((separation->type == CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK)
    && ((task->params->decomposeStrategy & CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_MASK)
    == CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT)))
  {
    CMRdbgMsg(10, "Pivoting for rank distribution.\n");

    /* Find nonzero in low-rank matrix. */
    size_t pivotRow = SIZE_MAX;
    size_t pivotColumn = SIZE_MAX;
    bool pivotBottomLeft;
    for (size_t row = 0; (row < node->numRows) && pivotRow == SIZE_MAX; ++row)
    {
      /* Neither a representative for first nor for second. */
      if (!(separation->rowsFlags[row] & CMR_SEPA_MASK_EXTRA))
        continue;

      int rowChild = separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD;
      size_t beyond = node->matrix->rowSlice[row + 1];
      for (size_t e = node->matrix->rowSlice[row]; e < beyond; ++e)
      {
        size_t column = node->matrix->entryColumns[e];
        int columnChild = separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD;
        if (columnChild != rowChild)
        {
          pivotRow = row;
          pivotColumn = column;
          pivotBottomLeft = rowChild == 1;
          break;
        }
      }
    }

    assert(pivotRow != SIZE_MAX);
    assert(pivotColumn != SIZE_MAX);

    CMRdbgMsg(10, "Pivot entry is r%zu,c%zu belonging to %s part.\n", pivotRow+1, pivotColumn+1,
      pivotBottomLeft ? "bottom-left" : "top-right");

    CMR_CHRMAT* childMatrix = NULL;
    CMR_SUBMAT* violatorSubmatrix = NULL;
    if (node->isTernary)
      CMR_CALL( CMRchrmatRegularPivot(cmr, node->matrix, pivotRow, pivotColumn, &violatorSubmatrix, &childMatrix) );
    else
      CMR_CALL( CMRchrmatBinaryPivot(cmr, node->matrix, pivotRow, pivotColumn, &childMatrix) );

    if (violatorSubmatrix)
    {
      CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
      assert(node->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

      CMRdbgMsg(10, "Pivoting produced a non-ternary entry. Computed 2x2 violator.\n");

      /* Tested in ThreeSumTruemperPivotHighRank unittest. */
      queue->foundIrregularity = true;

      return CMR_OKAY;
    }

    CMR_CHRMAT* childTransposed = NULL;
    CMR_CALL( CMRchrmatTranspose(cmr, childMatrix, &childTransposed) );

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
    CMR_CALL( CMRsepaFindBinaryRepresentatives(cmr, separation, childMatrix, childTransposed, NULL,
      node->isTernary ? &violatorSubmatrix : NULL) );

    if (violatorSubmatrix)
    {
      CMRdbgMsg(10, "-> %zux%zu submatrix with bad determinant.\n", violatorSubmatrix->numRows,
        violatorSubmatrix->numColumns);

      CMR_CALL( CMRseymourUpdateViolator(cmr, node, violatorSubmatrix) );
      assert(node->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

      /* TODO: Add as unittest. */
      queue->foundIrregularity = true;

      return CMR_OKAY;
    }

    CMR_CALL( CMRseymourUpdatePivots(cmr, node, 1, &pivotRow, &pivotColumn, childMatrix, childTransposed) );

    /* Recurse on child. */
    task->node = node->children[0];
    CMR_CALL( CMRregularityDecomposeThreeSum(cmr, task, queue, separation) );

    // TODO: how about row/col indices from child back to parent?

    return CMR_OKAY;
  }

  if (separation->type == CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS)
  {
    int distributedStrategy = task->params->decomposeStrategy & CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_MASK;
    if (distributedStrategy == CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM)
    {
      CMRdbgMsg(10, "Carrying out Delta-sum for a distributed-rank 3-separation.\n");

      node->type = CMR_SEYMOUR_NODE_TYPE_DELTASUM;
      CMR_CALL( CMRseymourSetNumChildren(cmr, node, 2) );

      char epsilon = 0;
      if (node->isTernary)
        CMR_CALL( CMRdeltasumDecomposeEpsilon(cmr, node->matrix, node->transpose, separation, &epsilon) );
      else
        epsilon = 1;

      /* Temporary data. */
      size_t* rowsToChild = NULL;
      size_t* columnsToChild = NULL;
      size_t* childRowsToParent = NULL;
      size_t* childColumnsToParent = NULL;
      CMR_CALL( CMRallocStackArray(cmr, &rowsToChild, node->matrix->numRows) );
      CMR_CALL( CMRallocStackArray(cmr, &columnsToChild, node->matrix->numColumns) );
      CMR_CALL( CMRallocStackArray(cmr, &childRowsToParent, node->matrix->numRows) );
      CMR_CALL( CMRallocStackArray(cmr, &childColumnsToParent, node->matrix->numColumns) );

      /* First child. */
      CMR_CHRMAT* first = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows[0], 1) );
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns[0], 2) );
      CMR_CALL( CMRdeltasumDecomposeFirst(cmr, node->matrix, separation, epsilon, &first, childRowsToParent,
        childColumnsToParent, rowsToChild, columnsToChild, node->childSpecialRows[0], node->childSpecialColumns[0]) );

#if defined(CMR_DEBUG)
      CMRdbgMsg(12, "First child matrix:\n");
      CMRchrmatPrintDense(cmr, first, stdout, '0', false);
#endif /* CMR_DEBUG */

      /* Create first decomposition node. */
      CMR_CALL( CMRseymourCreate(cmr, &node->children[0], node->isTernary, first, false) );

      /* Mapping from first child rows to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[0], first->numRows) );
      for (size_t row = 0; row < first->numRows; ++row)
        node->childRowsToParent[0][row] = CMRrowToElement(childRowsToParent[row]);

      /* Mapping from first child columns to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[0], first->numColumns) );
      for (size_t column = 0; column < first->numColumns; ++column)
        node->childColumnsToParent[0][column] = CMRcolumnToElement(childColumnsToParent[column]);

      /* Mapping from parent rows to first child rows. */
      for (size_t row = 0; row < node->matrix->numRows; ++row)
      {
        if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          node->rowsToChild[row] = rowsToChild[row];
      }

      /* Mapping from parent columns to first child columns. */
      for (size_t column = 0; column < node->matrix->numColumns; ++column)
      {
        if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          node->columnsToChild[column] = columnsToChild[column];
      }

      /* Second child. */
      CMR_CHRMAT* second = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows[1], 1) );
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns[1], 2) );
      CMR_CALL( CMRdeltasumDecomposeSecond(cmr, node->matrix, separation, epsilon, &second, childRowsToParent,
        childColumnsToParent, rowsToChild, columnsToChild, node->childSpecialRows[1], node->childSpecialColumns[1]) );

#if defined(CMR_DEBUG)
      CMRdbgMsg(12, "Second child matrix:\n");
      CMRchrmatPrintDense(cmr, second, stdout, '0', false);
#endif /* CMR_DEBUG */

      /* Create second decomposition node. */
      CMR_CALL( CMRseymourCreate(cmr, &node->children[1], node->isTernary, second, false) );

      /* Mapping from second child rows to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[1], second->numRows) );
      for (size_t row = 0; row < second->numRows; ++row)
        node->childRowsToParent[1][row] = CMRrowToElement(childRowsToParent[row]);

      /* Mapping from second child columns to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[1], second->numColumns) );
      for (size_t column = 0; column < second->numColumns; ++column)
        node->childColumnsToParent[1][column] = CMRcolumnToElement(childColumnsToParent[column]);

      /* Mapping from parent rows to second child rows. */
      for (size_t row = 0; row < node->matrix->numRows; ++row)
      {
        if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          node->rowsToChild[row] = rowsToChild[row];
      }

      /* Mapping from parent columns to second child columns. */
      for (size_t column = 0; column < node->matrix->numColumns; ++column)
      {
        if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          node->columnsToChild[column] = columnsToChild[column];
      }

      CMR_CALL( CMRfreeStackArray(cmr, &childColumnsToParent) );
      CMR_CALL( CMRfreeStackArray(cmr, &childRowsToParent) );
      CMR_CALL( CMRfreeStackArray(cmr, &columnsToChild) );
      CMR_CALL( CMRfreeStackArray(cmr, &rowsToChild) );
    }
    else if (distributedStrategy == CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM)
    {
      CMRdbgMsg(10, "Carrying out Y-sum for a distributed-rank 3-separation.\n");

      node->type = CMR_SEYMOUR_NODE_TYPE_YSUM;
      CMR_CALL( CMRseymourSetNumChildren(cmr, node, 2) );

      char epsilon = 0;
      if (node->isTernary)
        CMR_CALL( CMRysumDecomposeEpsilon(cmr, node->matrix, node->transpose, separation, &epsilon) );
      else
        epsilon = 1;

      /* Temporary data. */
      size_t* rowsToChild = NULL;
      size_t* columnsToChild = NULL;
      size_t* childRowsToParent = NULL;
      size_t* childColumnsToParent = NULL;
      CMR_CALL( CMRallocStackArray(cmr, &rowsToChild, node->matrix->numRows) );
      CMR_CALL( CMRallocStackArray(cmr, &columnsToChild, node->matrix->numColumns) );
      CMR_CALL( CMRallocStackArray(cmr, &childRowsToParent, node->matrix->numRows) );
      CMR_CALL( CMRallocStackArray(cmr, &childColumnsToParent, node->matrix->numColumns) );

      /* First child. */
      CMR_CHRMAT* first = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows[0], 2) );
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns[0], 1) );
      CMR_CALL( CMRysumDecomposeFirst(cmr, node->matrix, separation, epsilon, &first, childRowsToParent,
        childColumnsToParent, rowsToChild, columnsToChild, node->childSpecialRows[0], node->childSpecialColumns[0]) );

#if defined(CMR_DEBUG)
      CMRdbgMsg(12, "First child matrix:\n");
      CMRchrmatPrintDense(cmr, first, stdout, '0', false);
#endif /* CMR_DEBUG */

      /* Create first decomposition node. */
      CMR_CALL( CMRseymourCreate(cmr, &node->children[0], node->isTernary, first, false) );

      /* Mapping from first child rows to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[0], first->numRows) );
      for (size_t row = 0; row < first->numRows; ++row)
        node->childRowsToParent[0][row] = CMRrowToElement(childRowsToParent[row]);

      /* Mapping from first child columns to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[0], first->numColumns) );
      for (size_t column = 0; column < first->numColumns; ++column)
        node->childColumnsToParent[0][column] = CMRcolumnToElement(childColumnsToParent[column]);

      /* Mapping from parent rows to first child rows. */
      for (size_t row = 0; row < node->matrix->numRows; ++row)
      {
        if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          node->rowsToChild[row] = rowsToChild[row];
      }

      /* Mapping from parent columns to first child columns. */
      for (size_t column = 0; column < node->matrix->numColumns; ++column)
      {
        if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          node->columnsToChild[column] = columnsToChild[column];
      }

      /* Second child. */
      CMR_CHRMAT* second = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows[1], 2) );
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns[1], 1) );
      CMR_CALL( CMRysumDecomposeSecond(cmr, node->matrix, separation, epsilon, &second, childRowsToParent,
        childColumnsToParent, rowsToChild, columnsToChild, node->childSpecialRows[1], node->childSpecialColumns[1]) );

#if defined(CMR_DEBUG)
      CMRdbgMsg(12, "Second child matrix:\n");
      CMRchrmatPrintDense(cmr, second, stdout, '0', false);
#endif /* CMR_DEBUG */

      /* Create second decomposition node. */
      CMR_CALL( CMRseymourCreate(cmr, &node->children[1], node->isTernary, second, false) );

      /* Mapping from second child rows to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[1], second->numRows) );
      for (size_t row = 0; row < second->numRows; ++row)
        node->childRowsToParent[1][row] = CMRrowToElement(childRowsToParent[row]);

      /* Mapping from second child columns to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[1], second->numColumns) );
      for (size_t column = 0; column < second->numColumns; ++column)
        node->childColumnsToParent[1][column] = CMRcolumnToElement(childColumnsToParent[column]);

      /* Mapping from parent rows to second child rows. */
      for (size_t row = 0; row < node->matrix->numRows; ++row)
      {
        if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          node->rowsToChild[row] = rowsToChild[row];
      }

      /* Mapping from parent columns to second child columns. */
      for (size_t column = 0; column < node->matrix->numColumns; ++column)
      {
        if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          node->columnsToChild[column] = columnsToChild[column];
      }

      CMR_CALL( CMRfreeStackArray(cmr, &childColumnsToParent) );
      CMR_CALL( CMRfreeStackArray(cmr, &childRowsToParent) );
      CMR_CALL( CMRfreeStackArray(cmr, &columnsToChild) );
      CMR_CALL( CMRfreeStackArray(cmr, &rowsToChild) );
    }
    else
    {
      assert(0 == "Invalid sum strategy for distributed ranks!");
      return CMR_ERROR_INVALID;
    }
  }
  else
  {
    assert(separation->type == CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK);
    int concentratedStrategy = task->params->decomposeStrategy & CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_MASK;
    if (concentratedStrategy == CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM)
    {
      CMRdbgMsg(10, "Carrying out 3-sum for a concentrated-rank 3-separation.\n");

      node->type = CMR_SEYMOUR_NODE_TYPE_THREESUM;
      CMR_CALL( CMRseymourSetNumChildren(cmr, node, 2) );

      size_t specialRows[3];
      size_t specialColumns[3];
      char gamma, beta;
      CMR_CALL( CMRthreesumDecomposeConnecting(cmr, node->matrix, node->transpose, separation, specialRows,
        specialColumns, &gamma, &beta) );

      /* Temporary data. */
      size_t* rowsToChild = NULL;
      size_t* columnsToChild = NULL;
      size_t* childRowsToParent = NULL;
      size_t* childColumnsToParent = NULL;
      CMR_CALL( CMRallocStackArray(cmr, &rowsToChild, node->matrix->numRows) );
      CMR_CALL( CMRallocStackArray(cmr, &columnsToChild, node->matrix->numColumns) );
      CMR_CALL( CMRallocStackArray(cmr, &childRowsToParent, node->matrix->numRows) );
      CMR_CALL( CMRallocStackArray(cmr, &childColumnsToParent, node->matrix->numColumns) );

      /* First child. */
      CMR_CHRMAT* first = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows[0], 3) );
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns[0], 3) );
      CMR_CALL( CMRthreesumDecomposeFirst(cmr, node->matrix, separation, specialRows, specialColumns, beta,
        &first, childRowsToParent, childColumnsToParent, rowsToChild, columnsToChild, node->childSpecialRows[0],
        node->childSpecialColumns[0]) );

      /* Make it binary if necessary. */
      if (!node->isTernary)
        CMR_CALL( CMRchrmatSupport(cmr, first, &first) );

#if defined(CMR_DEBUG)
      CMRdbgMsg(12, "First child matrix:\n");
      CMRchrmatPrintDense(cmr, first, stdout, '0', false);
#endif /* CMR_DEBUG */

      /* Create first decomposition node. */
      CMR_CALL( CMRseymourCreate(cmr, &node->children[0], node->isTernary, first, false) );

      /* Mapping from first child rows to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[0], first->numRows) );
      for (size_t row = 0; row < first->numRows; ++row)
        node->childRowsToParent[0][row] = CMRrowToElement(childRowsToParent[row]);

      /* Mapping from first child columns to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[0], first->numColumns) );
      for (size_t column = 0; column < first->numColumns; ++column)
        node->childColumnsToParent[0][column] = CMRcolumnToElement(childColumnsToParent[column]);

      /* Mapping from parent rows to first child rows. */
      for (size_t row = 0; row < node->matrix->numRows; ++row)
      {
        if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          node->rowsToChild[row] = rowsToChild[row];
      }

      /* Mapping from parent columns to first child columns. */
      for (size_t column = 0; column < node->matrix->numColumns; ++column)
      {
        if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_FIRST)
          node->columnsToChild[column] = columnsToChild[column];
      }

      /* Second child. */
      CMR_CHRMAT* second = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialRows[1], 3) );
      CMR_CALL( CMRallocBlockArray(cmr, &node->childSpecialColumns[1], 3) );
      CMR_CALL( CMRthreesumDecomposeSecond(cmr, node->matrix, separation, specialRows, specialColumns, gamma,
        &second, childRowsToParent, childColumnsToParent, rowsToChild, columnsToChild, node->childSpecialRows[1],
        node->childSpecialColumns[1]) );

      /* Make it binary if necessary. */
      if (!node->isTernary)
        CMR_CALL( CMRchrmatSupport(cmr, second, &second) );

#if defined(CMR_DEBUG)
      CMRdbgMsg(12, "Second child matrix:\n");
      CMRchrmatPrintDense(cmr, second, stdout, '0', false);
#endif /* CMR_DEBUG */

      /* Create second decomposition node. */
      CMR_CALL( CMRseymourCreate(cmr, &node->children[1], node->isTernary, second, false) );

      /* Mapping from second child rows to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childRowsToParent[1], second->numRows) );
      for (size_t row = 0; row < second->numRows; ++row)
        node->childRowsToParent[1][row] = CMRrowToElement(childRowsToParent[row]);

      /* Mapping from second child columns to parent elements. */
      CMR_CALL( CMRallocBlockArray(cmr, &node->childColumnsToParent[1], second->numColumns) );
      for (size_t column = 0; column < second->numColumns; ++column)
        node->childColumnsToParent[1][column] = CMRcolumnToElement(childColumnsToParent[column]);

      /* Mapping from parent rows to second child rows. */
      for (size_t row = 0; row < node->matrix->numRows; ++row)
      {
        if ((separation->rowsFlags[row] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          node->rowsToChild[row] = rowsToChild[row];
      }

      /* Mapping from parent columns to second child columns. */
      for (size_t column = 0; column < node->matrix->numColumns; ++column)
      {
        if ((separation->columnsFlags[column] & CMR_SEPA_MASK_CHILD) == CMR_SEPA_SECOND)
          node->columnsToChild[column] = columnsToChild[column];
      }

      CMR_CALL( CMRfreeStackArray(cmr, &childColumnsToParent) );
      CMR_CALL( CMRfreeStackArray(cmr, &childRowsToParent) );
      CMR_CALL( CMRfreeStackArray(cmr, &columnsToChild) );
      CMR_CALL( CMRfreeStackArray(cmr, &rowsToChild) );
    }
    else
    {
      assert(0 == "Invalid 3-sum strategy for concentrated rank!");
      return CMR_ERROR_INVALID;
    }
  }

  if (!node->isTernary)
  {
    assert( CMRchrmatIsBinary(cmr, node->children[0]->matrix, NULL) );
    assert( CMRchrmatIsBinary(cmr, node->children[1]->matrix, NULL) );
  }

  DecompositionTask* childTasks[2] = { task, NULL };
  CMR_CALL( CMRregularityTaskCreateRoot(cmr, node->children[1], &childTasks[1], task->params, task->stats,
    task->startClock, task->timeLimit) );

  /* TODO: Find case where this can be SP-reducible. */

  childTasks[0]->node = node->children[0];
  node->children[0]->testedSeriesParallel = false;
  node->children[1]->testedSeriesParallel = false;

  /* Add both child tasks to the list. */
  CMRregularityQueueAdd(queue, childTasks[0]);
  CMRregularityQueueAdd(queue, childTasks[1]);

  return CMR_OKAY;
}
