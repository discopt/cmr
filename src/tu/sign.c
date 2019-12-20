#include <tu/sign.h>
#include "sign_internal.h"
#include "one_sum.h"

#include <assert.h>
#include <stdlib.h>

// TODO: merge wasVisited and isTarget into one structure.

typedef struct
{
  int status; /**< 0: not visited, 1: in queue, 2: processed */
  int predecessorNode; /**< Node number of predecessor. */
  char predecessorValue; /**< Value of matrix entry of predecessor. */
  char targetValue; /**< Entry in current row if a target node, and 0 otherwise. */
} GRAPH_NODE;

char signSequentiallyConnected(
  TU* tu,
  TU_SPARSE_CHAR* matrix,
  TU_SPARSE_CHAR* transpose,
  bool change
)
{
  assert(TUcheckSparseTransposeChar(matrix, transpose));
  assert(TUisTernaryChar(matrix));

  /* If we have more rows than columns, we work with the transpose. */
  if (matrix->numRows > matrix->numColumns)
  {
    char modified = signSequentiallyConnected(tu, transpose, matrix, change);
    assert(modified == 0 || modified == 'm');
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
    printf("Before processing row %d:\n", row);
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
      printf("Empty row.\n");
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
        printf("Current node is %d (row %d), queue length is %d\n", currentNode, r,
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
              int sum = graphNodes[c].targetValue;
              int pathNode = c;
              do
              {
                sum += graphNodes[pathNode].predecessorValue;
                pathNode = graphNodes[pathNode].predecessorNode;
              }
              while (graphNodes[pathNode].targetValue == 0);
              sum += graphNodes[pathNode].targetValue;
#ifdef DEBUG_SIGN
              printf("Found a chordless cycle between %d and %d with sum %d\n", c, pathNode, sum);
#endif
              if (sum % 4 != 0)
              {
                assert(sum % 4 == -2 || sum % 4 == 2);
                graphNodes[c].targetValue *= -1;
                if (change)
                  rowChanged = true;
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
        printf("Current node is %d (column %d), queue length is %d\n", currentNode, c,
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
        printf("Source ");
      if (graphNodes[v].targetValue != 0)
        printf("Target ");
      printf("Node %d is %s%d and has predecessor %d.\n", v, v >= firstRowNode ? "row ": "column ",
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
          printf("Sign change at %d,%d.\n", row, column);
#endif
        matrix->entryValues[e] = graphNodes[column].targetValue;
      }
    }
  }

#ifdef DEBUG_SIGN
  if (change)
  {
    printf("After signing:\n");
    TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
  }
#endif

  free(bfsQueue);
  free(graphNodes);

  return 'm'; // Or 0?
}

/**
 * \brief Signs a given ternary matrix.
 */

bool sign(
  TU* tu,                 /**< TU environment. */
  TU_SPARSE_CHAR* matrix, /**< Sparse matrix. */
  bool change             /**< Whether the signs of \p matrix shall be modified. */
)
{
  bool wasCorrect = true;
  int numComponents;
  TU_ONESUM_COMPONENT_CHAR* components;

  assert(TUisTernaryChar(matrix));

#ifdef DEBUG_SIGN
  printf("sign:\n");
  TUprintSparseAsDenseChar(stdout, matrix, ' ', true);
#endif

  /* Decompose into 1-connected components. */

  decomposeOneSumCharToChar(tu, matrix, &numComponents, &components, NULL, NULL, NULL, NULL);

  for (int comp = 0; comp < numComponents; ++comp)
  {
#ifdef DEBUG_SIGN
    printf("-> Component %d of size %dx%d\n", comp, components[comp].matrix.numRows,
      components[comp].matrix.numColumns);
#endif
    char modified = signSequentiallyConnected(tu, &components[comp].matrix,
      &components[comp].transpose, change);
    if (modified == 0)
      continue;

    wasCorrect = false;

    /* As we don't modify, we can abort early. */
    if (!change)
      break;

    assert(modified == 'm' || modified == 't');
    bool copyTranspose = modified == 't';

    /* Either the matrix or its transposed was modified. */
    TU_SPARSE_CHAR* sourceMatrix = 
      copyTranspose ? &components[comp].transpose : &components[comp].matrix;

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
    TUclearSparseChar(&components[c].matrix);
    TUclearSparseChar(&components[c].transpose);
  }
  free(components);

  return wasCorrect;
}

bool TUtestSignChar(TU* tu, TU_SPARSE_CHAR* matrix)
{
  return sign(tu, matrix, false);
}

bool TUcorrectSignChar(TU* tu, TU_SPARSE_CHAR* matrix)
{
  return sign(tu, matrix, true);
}
