#include "one_sum.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

struct GraphNode
{
  int adjacencyStart; /**< Index of first outgoing arc. */
  int component; /**< Index of component of matrix. */
  int degree; /**< Used to count the degree. */
  int order; /**< Corresponding row/column in component. */
};
typedef struct GraphNode GRAPH_NODE;

void decomposeOneSum(TU* tu, TU_MATRIX* matrix, size_t matrixType, size_t targetType,
  int* numComponents, TU_ONESUM_COMPONENT** components, int* rowsToComponents,
  int* columnsToComponents, int* rowsToComponentRows, int* columnsToComponentColumns)
{
  GRAPH_NODE* graphNodes = NULL;
  int* graphAdjacencies = NULL;
  int* queue = NULL;
  int queueLength = 0;
  int numNodes = matrix->numRows + matrix->numColumns;
  int countComponents = 0;
  const int firstColumnNode = matrix->numRows;
  int i;

  assert(tu != NULL);
  assert(matrix != NULL);
  assert(numComponents != NULL);
  assert(components != NULL);

  graphNodes = (GRAPH_NODE*) malloc((numNodes + 1) * sizeof(GRAPH_NODE));
  graphAdjacencies = (int*) malloc(2 * matrix->numNonzeros * sizeof(int));
  queue = (int*) malloc(numNodes * sizeof(int));

  for (int node = 0; node < numNodes; ++node)
  {
    graphNodes[node].component = -1;
    graphNodes[node].degree = 0;
  }

  /* Count degrees */
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int start = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ?  matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int e = start; e < end; ++e)
    {
      bool nonzero;
      if (matrixType == sizeof(double))
        nonzero = round(((double*)matrix->entryValues)[e]) != 0;
      else if (matrixType == sizeof(int))
        nonzero = ((int*)matrix->entryValues)[e] != 0;
      else if (matrixType == sizeof(char))
        nonzero = ((char*)matrix->entryValues)[e] != 0;
      else
        assert("Invalid matrixType parameter." == 0);
      if (nonzero)
      {
        int column = matrix->entryColumns[e];
        graphNodes[row].degree++;
        graphNodes[firstColumnNode + column].degree++;
      }
    }
  }

  /* Compute ranges for adjacencies */
  i = 0;
  for (int node = 0; node < numNodes; ++node)
  {
    graphNodes[node].adjacencyStart = i;
    i += graphNodes[node].degree;
  }
  graphNodes[numNodes].adjacencyStart = i;

  /* Scan entries and create adjacencies. */
  for (int row = 0; row < matrix->numRows; ++row)
  {
    int start = matrix->rowStarts[row];
    int end = row + 1 < matrix->numRows ? matrix->rowStarts[row+1] : matrix->numNonzeros;
    for (int e = start; e < end; ++e)
    {
      bool nonzero;
      if (matrixType == sizeof(double))
        nonzero = round(((double*)matrix->entryValues)[e]) != 0;
      else if (matrixType == sizeof(int))
        nonzero = ((int*)matrix->entryValues)[e] != 0;
      else if (matrixType == sizeof(char))
        nonzero = ((char*)matrix->entryValues)[e] != 0;
      else
        assert("Invalid matrixType parameter." == 0);      
      if (nonzero)
      {
        int column = matrix->entryColumns[e];
        int columnNode = firstColumnNode + column;
        graphAdjacencies[graphNodes[row + 1].adjacencyStart - graphNodes[row].degree] = columnNode;
        graphNodes[row].degree--;
        graphAdjacencies[graphNodes[columnNode + 1].adjacencyStart - graphNodes[columnNode].degree]
          = row;
        graphNodes[columnNode].degree--;
      }
    }
  }

  /* 
   * We decremented the degree entries, so they should be 0. From now on we can query 
   * graphNodes[node+1].adjacencyStart - graphNodes[node].adjacencyStart if necessary.
   * 
   * TODO: make degree a union with something else to save memory and improve cache behavior.
   */

  for (int node = 0; node < numNodes; ++node)
    assert(graphNodes[node].degree == 0);

  /* Run DFS. */
  for (int startNode = 0; startNode < numNodes; ++startNode)
  {
    if (graphNodes[startNode].component < 0)
    {
      /* Start a new component. */
      int currentOrderRow = 0;
      int currentOrderColumn = 0;

      graphNodes[startNode].component = countComponents;
      graphNodes[startNode].order = 0;
      if (startNode < firstColumnNode)
        currentOrderRow++;
      else
        currentOrderColumn++;
      queueLength = 1;
      queue[0] = startNode;
      while (queueLength > 0)
      {
        int currentNode = queue[--queueLength];
        int start = graphNodes[currentNode].adjacencyStart;
        int end = graphNodes[currentNode + 1].adjacencyStart;

        for (int i = start; i < end; ++i)
        {
          int endNode = graphAdjacencies[i];
          if (graphNodes[endNode].component < 0)
          {
            graphNodes[endNode].component = countComponents;
            if (endNode < firstColumnNode)
              graphNodes[endNode].order = currentOrderRow++;
            else
              graphNodes[endNode].order = currentOrderColumn++;
            queue[queueLength] = endNode;
            ++queueLength;
          }
        }
      }

      ++countComponents;
    }
  }

  *numComponents = countComponents;

  /* Allocate component data. */
  *components = (TU_ONESUM_COMPONENT*) malloc(countComponents*sizeof(TU_ONESUM_COMPONENT));

  /* Compute sizes. */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_MATRIX* compMatrix = &(*components)[component].matrix;
    compMatrix->numRows = 0;
    compMatrix->numColumns = 0;
    compMatrix->numNonzeros = 0;
  }
  for (int node = 0; node < numNodes; ++node)
  {
    int component = graphNodes[node].component;
    int start = graphNodes[node].adjacencyStart;
    int end = graphNodes[node + 1].adjacencyStart;
    assert(component >= 0);
    if (node < firstColumnNode)
    {
      (*components)[component].matrix.numRows++;
      (*components)[component].matrix.numNonzeros += end - start;
    }
    else
      (*components)[component].matrix.numColumns++;
  }

  /* Allocate memory */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_MATRIX* compMatrix = &(*components)[component].matrix;
    TU_MATRIX* compTranspose = &(*components)[component].transpose;

    (*components)[component].rowsToOriginal = (int*) malloc(compMatrix->numRows*sizeof(int));
    (*components)[component].columnsToOriginal = (int*) malloc(compMatrix->numColumns*sizeof(int));

    compMatrix->rowStarts = (int*) malloc((compMatrix->numRows + 1) * sizeof(int));
    compMatrix->entryColumns = (int*) malloc(compMatrix->numNonzeros * sizeof(int));
    compMatrix->entryValues = (int*) malloc(compMatrix->numNonzeros * targetType);

    compTranspose->numRows = compMatrix->numColumns;
    compTranspose->numColumns = compMatrix->numRows;
    compTranspose->numNonzeros = compMatrix->numNonzeros;
    compTranspose->rowStarts = (int*) malloc((compTranspose->numRows + 1) * sizeof(int));
    compTranspose->entryColumns = (int*) malloc(compTranspose->numNonzeros * sizeof(int));
    compTranspose->entryValues = (int*) malloc(compTranspose->numNonzeros * targetType);
  }

  /* Fill mapping arrays. */
  for (int node = 0; node < numNodes; ++node)
  {
    int component = graphNodes[node].component;
    int order = graphNodes[node].order;
    if (node < firstColumnNode)
      (*components)[component].rowsToOriginal[order] = node;
    else
      (*components)[component].columnsToOriginal[order] = node - firstColumnNode;
  }

  /* We can now fill the matrices of each component. */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_MATRIX* compTranspose = &(*components)[component].transpose;

    /* Compute the slices in the transposed component matrix from the graph. */
    int countNonzeros = 0;
    for (int compColumn = 0; compColumn < compTranspose->numRows; ++compColumn)
    {
      int column = (*components)[component].columnsToOriginal[compColumn];
      int node = firstColumnNode + column;
      compTranspose->rowStarts[compColumn] = countNonzeros;
      countNonzeros += graphNodes[node+1].adjacencyStart - graphNodes[node].adjacencyStart;
    }

    /* Fill the slices. To ensure that it is sorted, we iterate row-wise. */
    for (int compRow = 0; compRow < compTranspose->numColumns; ++compRow)
    {
      int row = (*components)[component].rowsToOriginal[compRow];
      int start = matrix->rowStarts[row];
      int end = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;

      /* Iterate over all entries of that row. */
      for (int matrixEntry = start; matrixEntry < end; ++matrixEntry)
      {
        bool nonzero;
        if (matrixType == sizeof(double))
          nonzero = round(((double*)matrix->entryValues)[matrixEntry]) != 0;
        else if (matrixType == sizeof(int))
          nonzero = ((int*)matrix->entryValues)[matrixEntry] != 0;
        else if (matrixType == sizeof(char))
          nonzero = ((char*)matrix->entryValues)[matrixEntry] != 0;
        else
          assert("Invalid matrixType parameter." == 0);
        if (nonzero)
        {
          int column = matrix->entryColumns[matrixEntry];
          int compColumn = graphNodes[firstColumnNode + column].order;
          int compEntry = compTranspose->rowStarts[compColumn];
          compTranspose->entryColumns[compEntry] = compRow;
          compTranspose->rowStarts[compColumn]++;
          if (targetType == sizeof(double) && matrixType == sizeof(double))
          {
            ((double*)compTranspose->entryValues)[compEntry] =
              ((double*)matrix->entryValues)[matrixEntry];
          }
          else if (targetType == sizeof(int) && matrixType == sizeof(double))
          {
            ((int*)compTranspose->entryValues)[compEntry] =
              (int)(round(((double*)matrix->entryValues)[matrixEntry]) + 0.5);
          }
          else if (targetType == sizeof(char) && matrixType == sizeof(double))
          {
            ((char*)compTranspose->entryValues)[compEntry] =
              (char)(round(((double*)matrix->entryValues)[matrixEntry]) + 0.5);
          }
          else if (targetType == sizeof(double) && matrixType == sizeof(int))
          {
            ((double*)compTranspose->entryValues)[compEntry] =
              ((int*)matrix->entryValues)[matrixEntry];
          }
          else if (targetType == sizeof(int) && matrixType == sizeof(int))
          {
            ((int*)compTranspose->entryValues)[compEntry] =
              ((int*)matrix->entryValues)[matrixEntry];
          }
          else if (targetType == sizeof(char) && matrixType == sizeof(int))
          {
            ((char*)compTranspose->entryValues)[compEntry] =
              ((int*)matrix->entryValues)[matrixEntry];
          }
          else if (targetType == sizeof(double) && matrixType == sizeof(char))
          {
            ((double*)compTranspose->entryValues)[compEntry] =
              ((char*)matrix->entryValues)[matrixEntry];
          }
          else if (targetType == sizeof(int) && matrixType == sizeof(char))
          {
            ((int*)compTranspose->entryValues)[compEntry] =
              ((char*)matrix->entryValues)[matrixEntry];
          }
          else if (targetType == sizeof(char) && matrixType == sizeof(char))
          {
            ((char*)compTranspose->entryValues)[compEntry] =
              ((char*)matrix->entryValues)[matrixEntry];
          }
          else
            assert("Invalid targetType / matrixType parameter combination." == 0);
        }
      }
    }

    /* Since we incremented the rowStarts for each nonzero, the array is shifted by one entry.
     * We restore this now. */
    for (int compColumn = compTranspose->numRows; compColumn > 0; --compColumn)
      compTranspose->rowStarts[compColumn] = compTranspose->rowStarts[compColumn-1];
    compTranspose->rowStarts[0] = 0;
  }

  /* We now create the row-wise representation from the column-wise one. */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_MATRIX* compMatrix = &(*components)[component].matrix;
    TU_MATRIX* compTranspose = &(*components)[component].transpose;

    /* Compute the slices in the component matrix from the graph. */
    int countNonzeros = 0;
    for (int compRow = 0; compRow < compMatrix->numRows; ++compRow)
    {
      int row = (*components)[component].rowsToOriginal[compRow];
      int node = row;
      compMatrix->rowStarts[compRow] = countNonzeros;
      countNonzeros += graphNodes[node+1].adjacencyStart - graphNodes[node].adjacencyStart;
    }

    /* Fill the slices. To ensure that it is sorted, we iterate column-wise. */
    for (int compColumn = 0; compColumn < compMatrix->numColumns; ++compColumn)
    {
      int start = compTranspose->rowStarts[compColumn];
      int end = compTranspose->rowStarts[compColumn + 1];

      /* Iterate over all entries of that row. */
      for (int compTransposeEntry = start; compTransposeEntry < end; ++compTransposeEntry)
      {
        int compRow = compTranspose->entryColumns[compTransposeEntry];
        int compMatrixEntry = compMatrix->rowStarts[compRow];
        compMatrix->entryColumns[compMatrixEntry] = compColumn;
        compMatrix->rowStarts[compRow]++;
        if (targetType == sizeof(double))
          ((double*)compMatrix->entryValues)[compMatrixEntry] = ((double*)compTranspose->entryValues)[compTransposeEntry];
        else if (targetType == sizeof(int))
          ((int*)compMatrix->entryValues)[compMatrixEntry] = ((int*)compTranspose->entryValues)[compTransposeEntry];
        else if (targetType == sizeof(char))
          ((char*)compMatrix->entryValues)[compMatrixEntry] = ((char*)compTranspose->entryValues)[compTransposeEntry];
        else
          assert("Invalid targetType parameter." == 0);
      }
    }

    /* Since we incremented the rowStarts for each nonzero, the array is shifted by one entry.
     * We restore this now. */
    for (int compRow = compMatrix->numRows; compRow > 0; --compRow)
      compMatrix->rowStarts[compRow] = compMatrix->rowStarts[compRow-1];
    compMatrix->rowStarts[0] = 0;
  }

  /* Fill arrays for original matrix viewpoint. */
  if (rowsToComponents != NULL)
  {
    for (int row = 0; row < matrix->numRows; ++row)
      rowsToComponents[row] = graphNodes[row].component;
  }
  if (columnsToComponents != NULL)
  {
    for (int column = 0; column < matrix->numColumns; ++column)
      columnsToComponents[column] = graphNodes[firstColumnNode + column].component;
  }
  if (rowsToComponentRows != NULL)
  {
    for (int row = 0; row < matrix->numRows; ++row)
      rowsToComponentRows[row] = graphNodes[row].order;
  }
  if (columnsToComponentColumns != NULL)
  {
    for (int column = 0; column < matrix->numColumns; ++column)
      columnsToComponentColumns[column] = graphNodes[firstColumnNode + column].order;
  }

  free(queue);
  free(graphAdjacencies);
  free(graphNodes);
}
