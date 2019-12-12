#include "one_sum.h"

#include <assert.h>
#include <stdlib.h>

struct GraphNode
{
  int adjacencyStart; /**< Index of first outgoing arc. */
  int component; /**< Index of component of matrix. */
  int degree; /**< Used to count the degree. */
  int order; /**< Corresponding row/column in component. */
};
typedef struct GraphNode GRAPH_NODE;

void decomposeOneSumIntInt(TU* tu, TU_SPARSE_INT* matrix, int* numComponents,
  TU_SPARSE_INT** compMatrices, TU_SPARSE_INT** compTransposes, int*** rowMapping,
  int*** columnMapping)
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
  assert(compMatrices != NULL);
  assert(compTransposes != NULL);
  assert(rowMapping != NULL);
  assert(columnMapping != NULL);

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
      if (matrix->entryValues[e] != 0)
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
      if (matrix->entryValues[e] != 0)
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
  *rowMapping = (int**) malloc(countComponents * sizeof(int*));
  *columnMapping = (int**) malloc(countComponents * sizeof(int*));
  *compMatrices = (TU_SPARSE_INT*) malloc(countComponents * sizeof(TU_SPARSE_INT));
  *compTransposes = (TU_SPARSE_INT*) malloc(countComponents * sizeof(TU_SPARSE_INT));

  /* Compute sizes. */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_SPARSE_INT* compMatrix = &(*compMatrices)[component];
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
      (*compMatrices)[component].numRows++;
      (*compMatrices)[component].numNonzeros += end - start;
    }
    else
      (*compMatrices)[component].numColumns++;    
  }

  /* Allocate memory */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_SPARSE_INT* compMatrix = &(*compMatrices)[component];
    TU_SPARSE_INT* compTranspose = &(*compTransposes)[component];

    (*rowMapping)[component] = (int*) malloc(compMatrix->numRows * sizeof(int));
    (*columnMapping)[component] = (int*) malloc(compMatrix->numColumns * sizeof(int));

    compMatrix->rowStarts = (int*) malloc((compMatrix->numRows + 1) * sizeof(int));
    compMatrix->entryColumns = (int*) malloc(compMatrix->numNonzeros * sizeof(int));
    compMatrix->entryValues = (int*) malloc(compMatrix->numNonzeros * sizeof(int));

    compTranspose->numRows = compMatrix->numColumns;
    compTranspose->numColumns = compMatrix->numRows;
    compTranspose->numNonzeros = compMatrix->numNonzeros;
    compTranspose->rowStarts = (int*) malloc((compTranspose->numRows + 1) * sizeof(int));
    compTranspose->entryColumns = (int*) malloc(compTranspose->numNonzeros * sizeof(int));
    compTranspose->entryValues = (int*) malloc(compTranspose->numNonzeros * sizeof(int));
  }

  /* Fill mapping arrays. */
  for (int node = 0; node < numNodes; ++node)
  {
    if (node < firstColumnNode)
      (*rowMapping)[graphNodes[node].component][graphNodes[node].order] = node;
    else
      (*columnMapping)[graphNodes[node].component][graphNodes[node].order] = node - firstColumnNode;
  }

  /* We can now fill the matrices of each component. */
  for (int component = 0; component < countComponents; ++component)
  {
    TU_SPARSE_INT* compTranspose = &(*compTransposes)[component];

    /* Compute the slices in the transposed component matrix from the graph. */
    int countNonzeros = 0;
    for (int compColumn = 0; compColumn < compTranspose->numRows; ++compColumn)
    {
      int column = (*columnMapping)[component][compColumn];
      int node = firstColumnNode + column;
      compTranspose->rowStarts[compColumn] = countNonzeros;
      countNonzeros += graphNodes[node+1].adjacencyStart - graphNodes[node].adjacencyStart;
    }

    /* Fill the slices. To ensure that it is sorted, we iterate row-wise. */
    for (int compRow = 0; compRow < compTranspose->numColumns; ++compRow)
    {
      int row = (*rowMapping)[component][compRow];
      int start = matrix->rowStarts[row];
      int end = row + 1 < matrix->numRows ? matrix->rowStarts[row + 1] : matrix->numNonzeros;

      /* Iterate over all entries of that row. */
      for (int matrixEntry = start; matrixEntry < end; ++matrixEntry)
      {
        int value = matrix->entryValues[matrixEntry];
        if (value != 0)
        {
          int column = matrix->entryColumns[matrixEntry];
          int compColumn = graphNodes[firstColumnNode + column].order;
          int compEntry = compTranspose->rowStarts[compColumn];
          compTranspose->entryColumns[compEntry] = compRow;
          compTranspose->entryValues[compEntry] = value;
          compTranspose->rowStarts[compColumn]++;
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
    TU_SPARSE_INT* compMatrix = &(*compMatrices)[component];
    TU_SPARSE_INT* compTranspose = &(*compTransposes)[component];

    /* Compute the slices in the component matrix from the graph. */
    int countNonzeros = 0;
    for (int compRow = 0; compRow < compMatrix->numRows; ++compRow)
    {
      int row = (*rowMapping)[component][compRow];
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
        compMatrix->entryValues[compMatrixEntry] = compTranspose->entryValues[compTransposeEntry];
        compMatrix->rowStarts[compRow]++;
      }
    }

    /* Since we incremented the rowStarts for each nonzero, the array is shifted by one entry.
     * We restore this now. */
    for (int compRow = compMatrix->numRows; compRow > 0; --compRow)
      compMatrix->rowStarts[compRow] = compMatrix->rowStarts[compRow-1];
    compMatrix->rowStarts[0] = 0;
  }

  free(queue);
  free(graphAdjacencies);
  free(graphNodes);
}
