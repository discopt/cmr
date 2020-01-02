#include <tu/sign.h>
#include "sign_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"

#include <assert.h>
#include <stdlib.h>

/**
 * \brief Graph node for BFS in signing algorithm.
 */

typedef struct
{
  int status; /**< 0: not visited, 1: in queue, 2: processed */
  int predecessorNode; /**< Node number of predecessor. */
  char predecessorValue; /**< Value of matrix entry of predecessor. */
  char targetValue; /**< Entry in current row if a target node, and 0 otherwise. */
} GRAPH_NODE;

char signSequentiallyConnected(
  TU* tu,
  TU_MATRIX_CHAR* matrix,
  TU_MATRIX_CHAR* transpose,
  bool change,
  TU_SUBMATRIX** submatrix
  )
{
  bool matrixChanged = false;
  
  assert(TUcheckMatrixTransposeChar(matrix, transpose));
  assert(TUisTernaryChar(matrix, NULL));

  /* If we have more rows than columns, we work with the transpose. */
  if (matrix->numRows > matrix->numColumns)
  {
    char modified = signSequentiallyConnected(tu, transpose, matrix, change, submatrix);
    assert(modified == 0 || modified == 'm');
    if (submatrix && *submatrix)
    {
      assert((*submatrix)->numRows == (*submatrix)->numColumns);
      int* tmp = (*submatrix)->rows;
      (*submatrix)->rows = (*submatrix)->columns;
      (*submatrix)->columns = tmp;
    }
    return modified == 0 ? 0 : 't';
  }

#ifdef DEBUG_SIGN
  printf("signSequentiallyConnected.\n");
#endif

  const int firstRowNode = matrix->numColumns;

  GRAPH_NODE* graphNodes = (GRAPH_NODE*) malloc((matrix->numColumns + matrix->numRows) * sizeof(GRAPH_NODE));
  int* bfsQueue = (int*) malloc((matrix->numColumns + matrix->numRows) * sizeof(int));
  int bfsQueueBegin = 0;
  int bfsQueueEnd = 0;

  /* Main loop iterates over the rows. */
  for (int row = 1; row < matrix->numRows; ++row)
  {
#ifdef DEBUG_SIGN
    printf("  Before processing row %d:\n", row);
    TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
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
#ifdef DEBUG_SIGN
      printf("  Empty row.\n");
#endif
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
#ifdef DEBUG_SIGN
        printf("    Current node is %d (row %d), queue length is %d\n", currentNode, r,
          bfsQueueEnd - bfsQueueBegin);
#endif

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
#ifdef DEBUG_SIGN
              printf("      Found a chordless cycle between %d and %d with sum %d of length %d\n",
                c, pathNode, sum, length);
#endif
              if (sum % 4 != 0)
              {
                assert(sum % 4 == -2 || sum % 4 == 2);

                /* If we didn't find a submatrix yet: */
                if (submatrix && *submatrix == NULL)
                {
                  int i = 1;
                  int j = 1;
                  TUcreateSubmatrix(submatrix, length/2, length/2);
                  pathNode = c;
                  (*submatrix)->columns[0] = c;
                  (*submatrix)->rows[0] = row;
                  do
                  {
                    pathNode = graphNodes[pathNode].predecessorNode;
                    if (pathNode >= firstRowNode)
                      (*submatrix)->rows[i++] = pathNode - firstRowNode;
                    else
                      (*submatrix)->columns[j++] = pathNode;
                  }
                  while (graphNodes[pathNode].targetValue == 0);
                  TUsortSubmatrix(*submatrix);
                  
#ifdef DEBUG_SIGN
                  printf("      Submatrix filled with %d rows and %d columns.\n", i, j);
#endif
                }
#ifdef DEBUG_SIGN
                printf("      Sign change required.\n");
#endif
                graphNodes[c].targetValue *= -1;
                if (change)
                {
                  rowChanged = true;
                  matrixChanged = true;
                }
                else
                {
                  free(bfsQueue);
                  free(graphNodes);
                  return 'm';
                }
              }
            }
          }
        }
      }
      else
      {
        int c = currentNode;
#ifdef DEBUG_SIGN
        printf("    Current node is %d (column %d), queue length is %d\n", currentNode, c,
          bfsQueueEnd - bfsQueueBegin);
#endif

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

#ifdef DEBUG_SIGN
    for (int v = 0; v < matrix->numColumns + row; ++v)
    {
      if (v == startNode)
        printf("    Source node ");
      else if (graphNodes[v].targetValue != 0)
        printf("    Target node ");
      else
        printf("    Node ");
      printf("%d is %s%d and has predecessor %d.\n", v, v >= firstRowNode ? "row ": "column ",
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
#ifdef DEBUG_SIGN
        if (matrix->entryValues[e] != graphNodes[column].targetValue)
          printf("  Sign change at %d,%d.\n", row, column);
#endif
        matrix->entryValues[e] = graphNodes[column].targetValue;
      }
    }
  }

#ifdef DEBUG_SIGN
  if (change)
  {
    printf("  After signing:\n");
    TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
  }
#endif

  free(bfsQueue);
  free(graphNodes);

  return matrixChanged ? 'm' : 0;
}

/**
 * \brief Signs a given ternary double matrix.
 * 
 * Returns \c true if and only if the signs were already correct.
 */

bool signDouble(
  TU* tu,                   /**< TU environment. */
  TU_MATRIX_DOUBLE* matrix, /**< Sparse double matrix. */
  bool change,              /**< Whether the signs of \p matrix shall be modified. */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  bool wasCorrect = true;
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(TUisTernaryDouble(matrix, 0.1, NULL));

#ifdef DEBUG_SIGN
  printf("sign:\n");
  TUprintSparseAsDenseDouble(stdout, matrix, ' ', true);
#endif

  /* Decompose into 1-connected components. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(double), sizeof(char), &numComponents, &components, NULL, NULL,
    NULL, NULL);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* compSubmatrix = NULL;

#ifdef DEBUG_SIGN
    printf("-> Component %d of size %dx%d\n", comp, components[comp].matrix.numRows,
      components[comp].matrix.numColumns);
#endif
    char modified = signSequentiallyConnected(tu, (TU_MATRIX_CHAR*) &components[comp].matrix,
      (TU_MATRIX_CHAR*) &components[comp].transpose, change, (submatrix && !*submatrix) ? &compSubmatrix : NULL);
#ifdef DEBUG_SIGN
    printf("-> Component %d yields: %c\n", comp, modified ? modified : '0');
#endif
    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    wasCorrect = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(submatrix && !*submatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (int r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (int c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      TUsortSubmatrix(compSubmatrix);
      *submatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_MATRIX_CHAR* sourceMatrix = copyTranspose ?
      (TU_MATRIX_CHAR*) &components[comp].transpose :
      (TU_MATRIX_CHAR*) &components[comp].matrix;

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

#ifdef DEBUG_SIGN
  if (!wasCorrect && change)
  {
    printf("Modified original matrix:\n");
    TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
  }
#endif

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].matrix);
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

  return wasCorrect;
}

bool TUtestSignDouble(TU* tu, TU_MATRIX_DOUBLE* matrix, TU_SUBMATRIX** submatrix)
{
  return signDouble(tu, matrix, false, submatrix);
}

bool TUcorrectSignDouble(TU* tu, TU_MATRIX_DOUBLE* matrix, TU_SUBMATRIX** submatrix)
{
  return signDouble(tu, matrix, true, submatrix);
}

/**
 * \brief Signs a given ternary int matrix.
 * 
 * Returns \c true if and only if the signs were already correct.
 */

bool signInt(
  TU* tu,                   /**< TU environment. */
  TU_MATRIX_INT* matrix,    /**< Sparse int matrix. */
  bool change,              /**< Whether the signs of \p matrix shall be modified. */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  bool wasCorrect = true;
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(TUisTernaryInt(matrix, NULL));

#ifdef DEBUG_SIGN
  printf("sign:\n");
  TUprintSparseAsDenseInt(stdout, matrix, ' ', true);
#endif

  /* Decompose into 1-connected components. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(int), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* compSubmatrix = NULL;

#ifdef DEBUG_SIGN
    printf("-> Component %d of size %dx%d\n", comp, components[comp].matrix.numRows,
      components[comp].matrix.numColumns);
#endif
    char modified = signSequentiallyConnected(tu, (TU_MATRIX_CHAR*) &components[comp].matrix,
      (TU_MATRIX_CHAR*) &components[comp].transpose, change, (submatrix &&
      !*submatrix) ? &compSubmatrix : NULL);
#ifdef DEBUG_SIGN
    printf("-> Component %d yields: %c\n", comp, modified ? modified : '0');
#endif
    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    wasCorrect = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(submatrix && !*submatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (int r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (int c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      TUsortSubmatrix(compSubmatrix);
      *submatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_MATRIX_CHAR* sourceMatrix = copyTranspose ?
      (TU_MATRIX_CHAR*) &components[comp].transpose :
      (TU_MATRIX_CHAR*) &components[comp].matrix;

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

#ifdef DEBUG_SIGN
  if (!wasCorrect && change)
  {
    printf("Modified original matrix:\n");
    TUprintSparseAsDenseInt(stdout, matrix, ' ', true);
  }
#endif

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].matrix);
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

  return wasCorrect;
}

bool TUtestSignInt(TU* tu, TU_MATRIX_INT* matrix, TU_SUBMATRIX** submatrix)
{
  return signInt(tu, matrix, false, submatrix);
}

bool TUcorrectSignInt(TU* tu, TU_MATRIX_INT* matrix, TU_SUBMATRIX** submatrix)
{
  return signInt(tu, matrix, true, submatrix);
}

/**
 * \brief Signs a given ternary char matrix.
 * 
 * Returns \c true if and only if the signs were already correct.
 */

bool signChar(
  TU* tu,                   /**< TU environment. */
  TU_MATRIX_CHAR* matrix,   /**< Sparse char matrix. */
  bool change,              /**< Whether the signs of \p matrix shall be modified. */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
)
{
  bool wasCorrect = true;
  int numComponents;
  TU_ONESUM_COMPONENT* components;

  assert(TUisTernaryChar(matrix, NULL));

#ifdef DEBUG_SIGN
  printf("sign:\n");
  TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
#endif

  /* Decompose into 1-connected components. */

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_SUBMATRIX* compSubmatrix = NULL;

#ifdef DEBUG_SIGN
    printf("-> Component %d of size %dx%d\n", comp, components[comp].matrix.numRows,
      components[comp].matrix.numColumns);
#endif
    char modified = signSequentiallyConnected(tu, (TU_MATRIX_CHAR*) &components[comp].matrix,
      (TU_MATRIX_CHAR*) &components[comp].transpose, change, (submatrix && !*submatrix) ? &compSubmatrix : NULL);
#ifdef DEBUG_SIGN
    printf("-> Component %d yields: %c\n", comp, modified ? modified : '0');
#endif
    if (modified == 0)
    {
      assert(compSubmatrix == NULL);
      continue;
    }

    wasCorrect = false;

    /* If we found a submatrix for the first time: */
    if (compSubmatrix)
    {
      assert(submatrix && !*submatrix);
      /* Translate component indices to indices of whole matrix and sort them again. */
      for (int r = 0; r < compSubmatrix->numRows; ++r)
        compSubmatrix->rows[r] = components[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (int c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = components[comp].columnsToOriginal[compSubmatrix->columns[c]];
      TUsortSubmatrix(compSubmatrix);
      *submatrix = compSubmatrix;
    }

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_MATRIX_CHAR* sourceMatrix = copyTranspose ?
      (TU_MATRIX_CHAR*) &components[comp].transpose :
      (TU_MATRIX_CHAR*) &components[comp].matrix;

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

#ifdef DEBUG_SIGN
  if (!wasCorrect && change)
  {
    printf("Modified original matrix:\n");
    TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
  }
#endif

  /* Clean-up */

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].matrix);
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

  return wasCorrect;
}

bool TUtestSignChar(TU* tu, TU_MATRIX_CHAR* matrix, TU_SUBMATRIX** submatrix)
{
  return signChar(tu, matrix, false, submatrix);
}

bool TUcorrectSignChar(TU* tu, TU_MATRIX_CHAR* matrix, TU_SUBMATRIX** submatrix)
{
  return signChar(tu, matrix, true, submatrix);
}
