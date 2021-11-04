// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/camion.h>

#include "camion_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "env_internal.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

/**
 * \brief Graph node for BFS in signing algorithm.
 */

typedef struct
{
  int status;             /**< \brief 0: not visited, 1: in queue, 2: processed */
  int predecessorNode;    /**< \brief Node number of predecessor. */
  char predecessorValue;  /**< \brief Value of matrix entry of predecessor. */
  char targetValue;       /**< \brief Entry in current row if a target node, and 0 otherwise. */
} GRAPH_NODE;

CMR_ERROR CMRcomputeCamionSignSequentiallyConnected(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< The matrix to be signed. */
  CMR_CHRMAT* transpose,  /**< The transpose of \p matrix. */
  bool change,            /**< Whether to modify the matrix. */
  char* pmodification,    /**< Pointer for storing which matrix was modified.*/
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
)
{
  assert(cmr);
  assert(matrix);
  assert(transpose);
  assert(pmodification);

  bool isTranspose;
  CMR_CALL( CMRchrmatCheckTranspose(cmr, matrix, transpose, &isTranspose) );
  assert(isTranspose);
  assert(CMRchrmatIsTernary(cmr, matrix, NULL));

  /* If we have more rows than columns, we work with the transpose. */
  if (matrix->numRows > matrix->numColumns)
  {
    CMR_CALL( CMRcomputeCamionSignSequentiallyConnected(cmr, transpose, matrix, change, pmodification, psubmatrix) );
    assert(*pmodification == 0 || *pmodification == 'm');
    if (psubmatrix && *psubmatrix)
    {
      assert((*psubmatrix)->numRows == (*psubmatrix)->numColumns);
      size_t* tmp = (*psubmatrix)->rows;
      (*psubmatrix)->rows = (*psubmatrix)->columns;
      (*psubmatrix)->columns = tmp;
    }
    if (*pmodification == 'm')
      *pmodification = 't';
    return CMR_OKAY;
  }

  CMRdbgMsg(2, "signSequentiallyConnected.\n");

  *pmodification = 0;
  const int firstRowNode = matrix->numColumns;
  GRAPH_NODE* graphNodes = NULL;
  int* bfsQueue = NULL;
  int bfsQueueBegin = 0;
  int bfsQueueEnd = 0;

  CMRallocStackArray(cmr, &graphNodes, matrix->numColumns + matrix->numRows);
  CMRallocStackArray(cmr, &bfsQueue, matrix->numColumns + matrix->numRows);

  /* Main loop iterates over the rows. */
  for (int row = 1; row < matrix->numRows; ++row)
  {
    CMRdbgMsg(2, "Before processing row %d:\n", row);
#if defined(CMR_DEBUG)
    CMRchrmatPrintDense(cmr, matrix, stdout, ' ', true);
#endif

    for (int v = 0; v < matrix->numColumns + matrix->numRows; ++v)
    {
      graphNodes[v].targetValue = 0;
      graphNodes[v].status = 0;
      graphNodes[v].predecessorNode = -1;
    }

    bool rowChanged = false;
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    if (first == beyond)
    {
      CMRdbgMsg(2, "Empty row.\n");
      continue;
    }

    /* First nonzero in row determines start column node. */
    int startNode = matrix->entryColumns[first];
    /* All columns of the row's nonzeros are target column nodes. */
    for (int e = first; e < beyond; ++e)
      graphNodes[matrix->entryColumns[e]].targetValue = matrix->entryValues[e];
    bfsQueue[0] = startNode;
    graphNodes[startNode].status = 1;
    bfsQueueBegin = 0;
    bfsQueueEnd = 1;

    while (bfsQueueBegin < bfsQueueEnd)
    {
      int currentNode = bfsQueue[bfsQueueBegin];
      assert(graphNodes[currentNode].status == 1);
      graphNodes[currentNode].status = 2;
      ++bfsQueueBegin;

      if (currentNode >= firstRowNode)
      {
        int r = currentNode - firstRowNode;
        CMRdbgMsg(4, "Current node is %d (row r%d), queue length is %d\n", currentNode, r+1,
          bfsQueueEnd - bfsQueueBegin);

        /* Iterate over outgoing edges. */
        first = matrix->rowSlice[r];
        beyond = matrix->rowSlice[r + 1];
        for (int e = first; e < beyond; ++e)
        {
          int c = matrix->entryColumns[e];
          if (graphNodes[c].status == 0)
          {
            graphNodes[c].status = 1;
            graphNodes[c].predecessorNode = currentNode;
            graphNodes[c].predecessorValue = matrix->entryValues[e];
            bfsQueue[bfsQueueEnd++] = c;
            /* If we reach a target node for the first time, we trace back to the previous target
               node (which might be the starting node). */
            if (graphNodes[c].targetValue != 0)
            {
              int length = 2;
              int sum = graphNodes[c].targetValue;
              CMRdbgMsg(8, "sum = %d\n", sum);
              int pathNode = c;
              do
              {
                sum += graphNodes[pathNode].predecessorValue;
                CMRdbgMsg(8, "sum = %d\n", sum);
                pathNode = graphNodes[pathNode].predecessorNode;
                ++length;
              }
              while (graphNodes[pathNode].targetValue == 0);
              sum += graphNodes[pathNode].targetValue;
              CMRdbgMsg(6, "Found a chordless cycle between c%d and c%d with sum %d of length %d\n", c+1, pathNode+1,
                sum, length);

              if (sum % 4 != 0)
              {
                assert(sum % 4 == -2 || sum % 4 == 2);

                /* If we didn't find a submatrix yet: */
                if (psubmatrix && *psubmatrix == NULL)
                {
                  int i = 1;
                  int j = 1;
                  CMR_CALL( CMRsubmatCreate(cmr, length/2, length/2, psubmatrix) );
                  CMR_SUBMAT* submatrix = *psubmatrix;
                  pathNode = c;
                  submatrix->columns[0] = c;
                  submatrix->rows[0] = row;
                  do
                  {
                    pathNode = graphNodes[pathNode].predecessorNode;
                    if (pathNode >= firstRowNode)
                      submatrix->rows[i++] = pathNode - firstRowNode;
                    else
                      submatrix->columns[j++] = pathNode;
                  }
                  while (graphNodes[pathNode].targetValue == 0);
                  CMR_CALL( CMRsortSubmatrix(cmr, submatrix) );

                  CMRdbgMsg(6, "Submatrix filled with %d rows and %d columns.\n", i, j);
                }
                CMRdbgMsg(6, "Sign change required.\n");
                graphNodes[c].targetValue *= -1;
                *pmodification = 'm';
                if (change)
                  rowChanged = true;
                else
                {
                  CMRfreeStackArray(cmr, &bfsQueue);
                  CMRfreeStackArray(cmr, &graphNodes);
                  return CMR_OKAY;
                }
              }
            }
          }
        }
      }
      else
      {
        int c = currentNode;
        CMRdbgMsg(4, "Current node is %d (column c%d), queue length is %d\n", currentNode, c+1,
          bfsQueueEnd - bfsQueueBegin);

        /* Iterate over outgoing edges. */
        first = transpose->rowSlice[c];
        beyond = transpose->rowSlice[c + 1];
        for (size_t e = first; e < beyond; ++e)
        {
          size_t r = transpose->entryColumns[e];
          /* Only rows before current iteration row participate. */
          if (r >= row)
            break;
          if (graphNodes[firstRowNode + r].status == 0)
          {
            graphNodes[firstRowNode + r].status = 1;
            graphNodes[firstRowNode + r].predecessorNode = currentNode;
            graphNodes[firstRowNode + r].predecessorValue = transpose->entryValues[e];
            bfsQueue[bfsQueueEnd++] = firstRowNode + r;
          }
        }
      }
    }

#if defined(CMR_DEBUG)
    for (int v = 0; v < matrix->numColumns + row; ++v)
    {
      if (v == startNode)
        CMRdbgMsg(4, "Source node ");
      else if (graphNodes[v].targetValue != 0)
        CMRdbgMsg(4, "Target node ");
      else
        CMRdbgMsg(4, "Node ");
      CMRdbgMsg(0, "%d is %s%d and has predecessor %d.\n", v, v >= firstRowNode ? "row r": "column c",
        (v >= firstRowNode ? v-firstRowNode : v) + 1, graphNodes[v].predecessorNode);
    }
#endif

    if (rowChanged)
    {
      first = matrix->rowSlice[row];
      beyond = matrix->rowSlice[row + 1];
      for (int e = first; e < beyond; ++e)
      {
        int column = matrix->entryColumns[e];
        if (matrix->entryValues[e] != graphNodes[column].targetValue)
        {
          CMRdbgMsg(2, "Sign change at r%d,c%d.\n", row+1, column+1);
          matrix->entryValues[e] = graphNodes[column].targetValue;
          size_t entry = SIZE_MAX;
          CMR_CALL( CMRchrmatFindEntry(cmr, transpose, column, row, &entry) );
          assert(entry < SIZE_MAX);
          assert(transpose->entryValues[entry] == -graphNodes[column].targetValue);
          transpose->entryValues[entry] = graphNodes[column].targetValue;
        }
      }
    }
  }

#if defined(CMR_DEBUG)
  if (change)
  {
    CMRdbgMsg(2, "After signing:\n");
    CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, ' ', true) );
  }
#endif /* CMR_DEBUG */

  CMRfreeStackArray(cmr, &bfsQueue);
  CMRfreeStackArray(cmr, &graphNodes);

  return CMR_OKAY;
}

/**
 * \brief Signs a given matrix.
 */

static
CMR_ERROR sign(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix \f$ M \f$. */
  bool change,            /**< Whether the signs of \f$ M \f$ shall be modified. */
  bool* pisCamionSigned,  /**< Pointer for storing whether \f$ M \f$ was already [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-camion submatrix (may be \c NULL). */
)
{
  assert(cmr);
  assert(matrix);
  assert(!psubmatrix || !*psubmatrix);

  size_t numComponents;
  CMR_ONESUM_COMPONENT* components = NULL;

  assert(CMRchrmatIsTernary(cmr, matrix, NULL));

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "sign:\n");
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  /* Decompose into 1-connected components. */

  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  if (pisCamionSigned)
    *pisCamionSigned = true;
  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMR_SUBMAT* compSubmatrix = NULL;

    CMRdbgMsg(2, "-> Component %d of size %dx%d\n", comp, components[comp].matrix->numRows,
      components[comp].matrix->numColumns);

    char modified;
    CMR_CALL( CMRcomputeCamionSignSequentiallyConnected(cmr, (CMR_CHRMAT*) components[comp].matrix,
      (CMR_CHRMAT*) components[comp].transpose, change, &modified,
      (psubmatrix && !*psubmatrix) ? &compSubmatrix : NULL) );

    CMRdbgMsg(2, "-> Component %d yields: %c\n", comp, modified ? modified : '0');

    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    if (pisCamionSigned)
      *pisCamionSigned = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(psubmatrix && !*psubmatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (size_t r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (size_t c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      CMRsortSubmatrix(cmr, compSubmatrix);
      *psubmatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
    {
      CMRdbgMsg(2, "Aborting early because we don't change the matrix.\n");
      break;
    }

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    CMR_CHRMAT* sourceMatrix = copyTranspose ?
      (CMR_CHRMAT*) components[comp].transpose :
      (CMR_CHRMAT*) components[comp].matrix;

    /* We have to copy the changes back to the original matrix. */
    for (size_t sourceRow = 0; sourceRow < sourceMatrix->numRows; ++sourceRow)
    {
      size_t sourceFirst = sourceMatrix->rowSlice[sourceRow];
      size_t sourceBeyond = sourceMatrix->rowSlice[sourceRow + 1];
      for (size_t  sourceEntry = sourceFirst; sourceEntry < sourceBeyond; ++sourceEntry)
      {
        size_t sourceColumn = sourceMatrix->entryColumns[sourceEntry];
        size_t compRow = copyTranspose ? sourceColumn : sourceRow;
        size_t compColumn = copyTranspose ? sourceRow : sourceColumn;
        size_t row = components[comp].rowsToOriginal[compRow];
        size_t column = components[comp].columnsToOriginal[compColumn];

        CMRdbgMsg(4, "Searching entry for row %d and column %d.\n", row, column);

        /* Perform binary search in row of original matrix to find the column. */

        size_t lower = matrix->rowSlice[row];
        size_t upper = matrix->rowSlice[row + 1];
        while (lower < upper)
        {
          size_t entry = (lower + upper) / 2;
          size_t searchColumn = matrix->entryColumns[entry];
          if (column < searchColumn)
            upper = entry;
          else if (column > searchColumn)
            lower = entry + 1;
          else
          {
            CMRdbgMsg(4, "Original matrix entry %d (%d) is replaced by source matrix entry %d (%d).\n", entry,
              matrix->entryValues[entry], sourceEntry, sourceMatrix->entryValues[sourceEntry]);
            matrix->entryValues[entry] = sourceMatrix->entryValues[sourceEntry];
            break;
          }
        }
        assert(lower < upper);
      }
    }
  }

#if defined(CMR_DEBUG)
  if (pisCamionSigned && !*pisCamionSigned && change)
  {
    CMRdbgMsg(0, "Modified original matrix:\n");
    CMRchrmatPrintDense(cmr, matrix, stdout, ' ', true);
  }
#endif /* CMR_DEBUG */

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
    CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &components);

  return CMR_OKAY;
}

CMR_ERROR CMRtestCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bool* pisCamionSigned, CMR_SUBMAT** psubmatrix)
{
  return sign(cmr, matrix, false, pisCamionSigned, psubmatrix);
}

CMR_ERROR CMRcomputeCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bool* pwasCamionSigned, CMR_SUBMAT** psubmatrix)
{
  return sign(cmr, matrix, true, pwasCamionSigned, psubmatrix);
}
