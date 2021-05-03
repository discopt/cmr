#define TU_DEBUG /* Uncomment to debug graphic. */

#include <tu/sign.h>
#include "sign_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "env_internal.h"

#include <assert.h>
#include <stdlib.h>

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

TU_ERROR signSequentiallyConnected(
  TU* tu,                 /**< \brief \ref TU environment. */
  TU_CHRMAT* matrix,      /**< \brief The matrix to be signed. */
  TU_CHRMAT* transpose,   /**< \brief The transpose of \p matrix. */
  bool change,            /**< \brief Whether to modify the matrix. */
  char* pmodification,    /**< \brief Pointer for storing which matrix was modified.*/
  TU_SUBMAT** psubmatrix  /**< \brief If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  assert(tu);
  assert(matrix);
  assert(transpose);
  assert(pmodification);

  assert(TUchrmatCheckTranspose(matrix, transpose));
  assert(TUisTernaryChr(tu, matrix, NULL));

  /* If we have more rows than columns, we work with the transpose. */
  if (matrix->numRows > matrix->numColumns)
  {
    TU_CALL( signSequentiallyConnected(tu, transpose, matrix, change, pmodification, psubmatrix) );
    assert(*pmodification == 0 || *pmodification == 'm');
    if (psubmatrix && *psubmatrix)
    {
      assert((*psubmatrix)->numRows == (*psubmatrix)->numColumns);
      int* tmp = (*psubmatrix)->rows;
      (*psubmatrix)->rows = (*psubmatrix)->columns;
      (*psubmatrix)->columns = tmp;
    }
    if (*pmodification == 'm')
      *pmodification = 't';
    return TU_OKAY;
  }

  TUdbgMsg(2, "signSequentiallyConnected.\n");

  *pmodification = 0;
  const int firstRowNode = matrix->numColumns;
  GRAPH_NODE* graphNodes = NULL;
  int* bfsQueue = NULL;
  int bfsQueueBegin = 0;
  int bfsQueueEnd = 0;

  TUallocStackArray(tu, &graphNodes, matrix->numColumns + matrix->numRows);
  TUallocStackArray(tu, &bfsQueue, matrix->numColumns + matrix->numRows);

  /* Main loop iterates over the rows. */
  for (int row = 1; row < matrix->numRows; ++row)
  {
    TUdbgMsg(2, "Before processing row %d:\n", row);
#if defined(TU_DEBUG)
    TUchrmatPrintDense(stdout, matrix, ' ', true);
#endif

    for (int v = 0; v < matrix->numColumns + matrix->numRows; ++v)
    {
      graphNodes[v].targetValue = 0;
      graphNodes[v].status = 0;
      graphNodes[v].predecessorNode = -1;
    }

    bool rowChanged = false;
    int begin = matrix->rowStarts[row];
    int end = matrix->rowStarts[row+1];
    if (begin == end)
    {
      TUdbgMsg(2, "Empty row.\n");
      continue;
    }

    /* First nonzero in row determines start column node. */
    int startNode = matrix->entryColumns[begin];
    /* All columns of the row's nonzeros are target column nodes. */
    for (int e = begin; e < end; ++e)
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
        TUdbgMsg(4, "Current node is %d (row %d), queue length is %d\n", currentNode, r, bfsQueueEnd - bfsQueueBegin);

        /* Iterate over outgoing edges. */
        begin = matrix->rowStarts[r];
        end = matrix->rowStarts[r+1];
        for (int e = begin; e < end; ++e)
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
              int pathNode = c;
              do
              {
                sum += graphNodes[pathNode].predecessorValue;
                pathNode = graphNodes[pathNode].predecessorNode;
                ++length;
              }
              while (graphNodes[pathNode].targetValue == 0);
              sum += graphNodes[pathNode].targetValue;
              TUdbgMsg(6, "Found a chordless cycle between %d and %d with sum %d of length %d\n", c, pathNode, sum,
                length);

              if (sum % 4 != 0)
              {
                assert(sum % 4 == -2 || sum % 4 == 2);

                /* If we didn't find a submatrix yet: */
                if (psubmatrix && *psubmatrix == NULL)
                {
                  int i = 1;
                  int j = 1;
                  TUsubmatCreate(tu, psubmatrix, length/2, length/2);
                  pathNode = c;
                  (*psubmatrix)->columns[0] = c;
                  (*psubmatrix)->rows[0] = row;
                  do
                  {
                    pathNode = graphNodes[pathNode].predecessorNode;
                    if (pathNode >= firstRowNode)
                      (*psubmatrix)->rows[i++] = pathNode - firstRowNode;
                    else
                      (*psubmatrix)->columns[j++] = pathNode;
                  }
                  while (graphNodes[pathNode].targetValue == 0);
                  TUsortSubmatrix(*psubmatrix);

                  TUdbgMsg(6, "Submatrix filled with %d rows and %d columns.\n", i, j);
                }
                TUdbgMsg(6, "Sign change required.\n");
                graphNodes[c].targetValue *= -1;
                *pmodification = 'm';
                if (change)
                  rowChanged = true;
                else
                {
                  TUfreeStackArray(tu, &bfsQueue);
                  TUfreeStackArray(tu, &graphNodes);
                  return TU_OKAY;
                }
              }
            }
          }
        }
      }
      else
      {
        int c = currentNode;
        TUdbgMsg(4, "Current node is %d (column %d), queue length is %d\n", currentNode, c,
          bfsQueueEnd - bfsQueueBegin);

        /* Iterate over outgoing edges. */
        begin = transpose->rowStarts[c];
        end = transpose->rowStarts[c+1];
        for (int e = begin; e < end; ++e)
        {
          int r = transpose->entryColumns[e];
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

#if defined(TU_DEBUG)
    for (int v = 0; v < matrix->numColumns + row; ++v)
    {
      if (v == startNode)
        TUdbgMsg(4, "Source node ");
      else if (graphNodes[v].targetValue != 0)
        TUdbgMsg(4, "Target node ");
      else
        TUdbgMsg(4, "Node ");
      TUdbgMsg(0, "%d is %s%d and has predecessor %d.\n", v, v >= firstRowNode ? "row ": "column ",
        v >= firstRowNode ? v-firstRowNode : v, graphNodes[v].predecessorNode);
    }
#endif

    if (rowChanged)
    {
      begin = matrix->rowStarts[row];
      end = matrix->rowStarts[row+1];
      for (int e = begin; e < end; ++e)
      {
        int column = matrix->entryColumns[e];
        if (matrix->entryValues[e] != graphNodes[column].targetValue)
          TUdbgMsg(2, "Sign change at %d,%d.\n", row, column);
        matrix->entryValues[e] = graphNodes[column].targetValue;
      }
    }
  }

#if defined(TU_DEBUG)
  if (change)
  {
    TUdbgMsg(2, "After signing:\n");
    TU_CALL( TUchrmatPrintDense(stdout, matrix, ' ', true) );
  }
#endif /* TU_DEBUG */

  TUfreeStackArray(tu, &bfsQueue);
  TUfreeStackArray(tu, &graphNodes);

  return TU_OKAY;
}

/**
 * \brief Signs a given ternary double matrix.
 */

static
TU_ERROR signDbl(
  TU* tu,                 /**< \brief \ref TU environment. */
  TU_DBLMAT* matrix,      /**< \brief Sparse double matrix. */
  bool change,            /**< \brief Whether the signs of \p matrix shall be modified. */
  bool* palreadySigned,   /**< \brief Pointer for storing whether \p matrix was already signed correctly. */
  TU_SUBMAT** psubmatrix  /**< \brief If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  assert(tu);
  assert(matrix);
  assert(palreadySigned);

  int numComponents;
  TU_ONESUM_COMPONENT* components = NULL;

  assert(TUisTernaryDbl(tu, matrix, 1.0e-3, NULL));

#if defined(TU_DEBUG)
  TUdbgMsg(0, "sign:\n");
  TUdblmatPrintDense(stdout, matrix, ' ', true);
#endif /* TU_DEBUG */

  /* Decompose into 1-connected components. */

  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(double), sizeof(double), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  *palreadySigned = true;
  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMAT* compSubmatrix = NULL;

    TUdbgMsg(2, "-> Component %d of size %dx%d\n", comp, components[comp].matrix->numRows,
      components[comp].matrix->numColumns);

    char modified;
    TU_CALL( signSequentiallyConnected(tu, (TU_CHRMAT*) components[comp].matrix,
      (TU_CHRMAT*) components[comp].transpose, change, &modified,
      (psubmatrix && !*psubmatrix) ? &compSubmatrix : NULL) );

    TUdbgMsg(2, "-> Component %d yields: %c\n", comp, modified ? modified : '0');

    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    *palreadySigned = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(psubmatrix && !*psubmatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (int r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (int c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      TUsortSubmatrix(compSubmatrix);
      *psubmatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_DBLMAT* sourceMatrix = copyTranspose ?
      (TU_DBLMAT*) components[comp].transpose :
      (TU_DBLMAT*) components[comp].matrix;

    /* We have to copy the changes back to the original matrix. */
    for (int sourceRow = 0; sourceRow < sourceMatrix->numRows; ++sourceRow)
    {
      int sourceBegin = sourceMatrix->rowStarts[sourceRow];
      int sourceEnd = sourceMatrix->rowStarts[sourceRow + 1];
      for (int sourceEntry = sourceBegin; sourceEntry < sourceEnd; ++sourceEntry)
      {
        int sourceColumn = sourceMatrix->entryColumns[sourceEntry];
        int compRow = copyTranspose ? sourceColumn : sourceRow;
        int compColumn = copyTranspose ? sourceRow : sourceColumn;
        int row = components[comp].rowsToOriginal[compRow];
        int column = components[comp].columnsToOriginal[compColumn];

        /* Perform binary search in row of original matrix to find the column. */

        int lower = matrix->rowStarts[row];
        int upper = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
        while (lower < upper)
        {
          int entry = (lower + upper) / 2;
          int searchColumn = matrix->entryColumns[entry];
          if (column < searchColumn)
            upper = entry;
          else if (column > searchColumn)
            lower = entry + 1;
          else
          {
            matrix->entryValues[entry] = sourceMatrix->entryValues[sourceEntry];
            break;
          }
        }
        assert(lower < upper);
      }
    }
  }

#if defined(TU_DEBUG)
  if (!*palreadySigned && change)
  {
    TUdbgMsg(0, "Modified original matrix:\n");
    TUdblmatPrintDense(stdout, matrix, ' ', true);
  }
#endif /* TU_DEBUG */

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return TU_OKAY;
}

TU_ERROR TUtestSignDbl(TU* tu, TU_DBLMAT* matrix, bool* palreadySigned, TU_SUBMAT** psubmatrix)
{
  return signDbl(tu, matrix, false, palreadySigned, psubmatrix);
}

TU_ERROR TUcorrectSignDbl(TU* tu, TU_DBLMAT* matrix, bool* palreadySigned, TU_SUBMAT** psubmatrix)
{
  return signDbl(tu, matrix, true, palreadySigned, psubmatrix);
}

/**
 * \brief Signs a given ternary int matrix.
 */

static
TU_ERROR signInt(
  TU* tu,                 /**< \brief \ref TU environment. */
  TU_INTMAT* matrix,      /**< \brief Sparse int matrix. */
  bool change,            /**< \brief Whether the signs of \p matrix shall be modified. */
  bool* palreadySigned,   /**< \brief Pointer for storing whether \p matrix was already signed correctly. */
  TU_SUBMAT** psubmatrix  /**< \brief If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  assert(tu);
  assert(matrix);
  assert(palreadySigned);

  int numComponents;
  TU_ONESUM_COMPONENT* components = NULL;

  assert(TUisTernaryInt(tu, matrix, NULL));

#if defined(TU_DEBUG)
  TUdbgMsg(0, "sign:\n");
  TUintmatPrintDense(stdout, matrix, ' ', true);
#endif /* TU_DEBUG */

  /* Decompose into 1-connected components. */

  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(int), sizeof(int), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  *palreadySigned = true;
  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMAT* compSubmatrix = NULL;

    TUdbgMsg(2, "-> Component %d of size %dx%d\n", comp, components[comp].matrix->numRows,
      components[comp].matrix->numColumns);

    char modified;
    TU_CALL( signSequentiallyConnected(tu, (TU_CHRMAT*) components[comp].matrix,
      (TU_CHRMAT*) components[comp].transpose, change, &modified,
      (psubmatrix && !*psubmatrix) ? &compSubmatrix : NULL) );

    TUdbgMsg(2, "-> Component %d yields: %c\n", comp, modified ? modified : '0');

    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    *palreadySigned = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(psubmatrix && !*psubmatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (int r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (int c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      TUsortSubmatrix(compSubmatrix);
      *psubmatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_CHRMAT* sourceMatrix = copyTranspose ?
      (TU_CHRMAT*) components[comp].transpose :
      (TU_CHRMAT*) components[comp].matrix;

    /* We have to copy the changes back to the original matrix. */
    for (int sourceRow = 0; sourceRow < sourceMatrix->numRows; ++sourceRow)
    {
      int sourceBegin = sourceMatrix->rowStarts[sourceRow];
      int sourceEnd = sourceMatrix->rowStarts[sourceRow + 1];
      for (int sourceEntry = sourceBegin; sourceEntry < sourceEnd; ++sourceEntry)
      {
        int sourceColumn = sourceMatrix->entryColumns[sourceEntry];
        int compRow = copyTranspose ? sourceColumn : sourceRow;
        int compColumn = copyTranspose ? sourceRow : sourceColumn;
        int row = components[comp].rowsToOriginal[compRow];
        int column = components[comp].columnsToOriginal[compColumn];

        /* Perform binary search in row of original matrix to find the column. */

        int lower = matrix->rowStarts[row];
        int upper = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
        while (lower < upper)
        {
          int entry = (lower + upper) / 2;
          int searchColumn = matrix->entryColumns[entry];
          if (column < searchColumn)
            upper = entry;
          else if (column > searchColumn)
            lower = entry + 1;
          else
          {
            matrix->entryValues[entry] = sourceMatrix->entryValues[sourceEntry];
            break;
          }
        }
        assert(lower < upper);
      }
    }
  }

#if defined(TU_DEBUG)
  if (!*palreadySigned && change)
  {
    TUdbgMsg(0, "Modified original matrix:\n");
    TUintmatPrintDense(stdout, matrix, ' ', true);
  }
#endif /* TU_DEBUG */

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return TU_OKAY;
}

TU_ERROR TUtestSignInt(TU* tu, TU_INTMAT* matrix, bool* palreadySigned, TU_SUBMAT** psubmatrix)
{
  return signInt(tu, matrix, false, palreadySigned, psubmatrix);
}

TU_ERROR TUcorrectSignInt(TU* tu, TU_INTMAT* matrix, bool* palreadySigned, TU_SUBMAT** psubmatrix)
{
  return signInt(tu, matrix, true, palreadySigned, psubmatrix);
}


/**
 * \brief Signs a given ternary char matrix.
 */

static
TU_ERROR signChr(
  TU* tu,                 /**< \brief \ref TU environment. */
  TU_CHRMAT* matrix,      /**< \brief Sparse char matrix. */
  bool change,            /**< \brief Whether the signs of \p matrix shall be modified. */
  bool* palreadySigned,   /**< \brief Pointer for storing whether \p matrix was already signed correctly. */
  TU_SUBMAT** psubmatrix  /**< \brief If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  assert(tu);
  assert(matrix);
  assert(palreadySigned);

  int numComponents;
  TU_ONESUM_COMPONENT* components = NULL;

  assert(TUisTernaryChr(tu, matrix, NULL));

#if defined(TU_DEBUG)
  TUdbgMsg(0, "sign:\n");
  TUchrmatPrintDense(stdout, matrix, ' ', true);
#endif /* TU_DEBUG */

  /* Decompose into 1-connected components. */

  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL) );

  *palreadySigned = true;
  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMAT* compSubmatrix = NULL;

    TUdbgMsg(2, "-> Component %d of size %dx%d\n", comp, components[comp].matrix->numRows,
      components[comp].matrix->numColumns);

    char modified;
    TU_CALL( signSequentiallyConnected(tu, (TU_CHRMAT*) components[comp].matrix,
      (TU_CHRMAT*) components[comp].transpose, change, &modified,
      (psubmatrix && !*psubmatrix) ? &compSubmatrix : NULL) );

    TUdbgMsg(2, "-> Component %d yields: %c\n", comp, modified ? modified : '0');

    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    *palreadySigned = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(psubmatrix && !*psubmatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (int r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (int c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      TUsortSubmatrix(compSubmatrix);
      *psubmatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_CHRMAT* sourceMatrix = copyTranspose ?
      (TU_CHRMAT*) components[comp].transpose :
      (TU_CHRMAT*) components[comp].matrix;

    /* We have to copy the changes back to the original matrix. */
    for (int sourceRow = 0; sourceRow < sourceMatrix->numRows; ++sourceRow)
    {
      int sourceBegin = sourceMatrix->rowStarts[sourceRow];
      int sourceEnd = sourceMatrix->rowStarts[sourceRow + 1];
      for (int sourceEntry = sourceBegin; sourceEntry < sourceEnd; ++sourceEntry)
      {
        int sourceColumn = sourceMatrix->entryColumns[sourceEntry];
        int compRow = copyTranspose ? sourceColumn : sourceRow;
        int compColumn = copyTranspose ? sourceRow : sourceColumn;
        int row = components[comp].rowsToOriginal[compRow];
        int column = components[comp].columnsToOriginal[compColumn];

        /* Perform binary search in row of original matrix to find the column. */

        int lower = matrix->rowStarts[row];
        int upper = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;
        while (lower < upper)
        {
          int entry = (lower + upper) / 2;
          int searchColumn = matrix->entryColumns[entry];
          if (column < searchColumn)
            upper = entry;
          else if (column > searchColumn)
            lower = entry + 1;
          else
          {
            matrix->entryValues[entry] = sourceMatrix->entryValues[sourceEntry];
            break;
          }
        }
        assert(lower < upper);
      }
    }
  }

#if defined(TU_DEBUG)
  if (!*palreadySigned && change)
  {
    TUdbgMsg(0, "Modified original matrix:\n");
    TUchrmatPrintDense(stdout, matrix, ' ', true);
  }
#endif /* TU_DEBUG */

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return TU_OKAY;
}

TU_ERROR TUtestSignChr(TU* tu, TU_CHRMAT* matrix, bool* palreadySigned, TU_SUBMAT** psubmatrix)
{
  return signChr(tu, matrix, false, palreadySigned, psubmatrix);
}

TU_ERROR TUcorrectSignChr(TU* tu, TU_CHRMAT* matrix, bool* palreadySigned, TU_SUBMAT** psubmatrix)
{
  return signChr(tu, matrix, true, palreadySigned, psubmatrix);
}
