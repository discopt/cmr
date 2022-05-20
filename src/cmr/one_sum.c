// #define CMR_DEBUG /* Uncomment to debug the 1-sum decomposition. */

#include "one_sum.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "env_internal.h"

struct GraphNode
{
  int adjacencyStart; /**< \brief Index of first outgoing arc. */
  int component;      /**< \brief Index of component of matrix. */
  int degree;         /**< \brief Used to count the degree. */
  int order;          /**< \brief Corresponding row/column in component. */
};
typedef struct GraphNode GRAPH_NODE;

CMR_ERROR decomposeOneSum(CMR* cmr, CMR_MATRIX* matrix, size_t matrixType, size_t targetType,
  size_t* pnumComponents, CMR_ONESUM_COMPONENT** pcomponents, size_t* rowsToComponents,
  size_t* columnsToComponents, size_t* rowsToComponentRows, size_t* columnsToComponentColumns)
{
  GRAPH_NODE* graphNodes = NULL;
  int* graphAdjacencies = NULL;
  int* queue = NULL;
  int queueLength = 0;
  int numNodes = matrix->numRows + matrix->numColumns;
  int countComponents = 0;
  const int firstColumnNode = matrix->numRows;
  int i;

  assert(cmr);
  assert(matrix);
  assert(pnumComponents);
  assert(pcomponents);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "decomposeOneSum:\n");
  if (matrixType == sizeof(double))
    CMRdblmatPrintDense(stdout, (CMR_DBLMAT*) matrix, '0', true);
  else if (matrixType == sizeof(int))
    CMRintmatPrintDense(stdout, (CMR_INTMAT*) matrix, '0', true);
  else if (matrixType == sizeof(char))
    CMRchrmatPrintDense(stdout, (CMR_CHRMAT*) matrix, '0', true);
#endif

  CMR_CALL( CMRallocStackArray(cmr, &graphNodes, numNodes + 1) );
  CMR_CALL( CMRallocStackArray(cmr, &graphAdjacencies, 2 * matrix->numNonzeros) );
  CMR_CALL( CMRallocStackArray(cmr, &queue, numNodes) );

  for (int node = 0; node < numNodes; ++node)
  {
    graphNodes[node].component = -1;
    graphNodes[node].degree = 0;
  }

  /* Count degrees */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      bool isNonzero;
      if (matrixType == sizeof(double))
        isNonzero = ((double*)matrix->entryValues)[e] != 0.0;
      else if (matrixType == sizeof(int))
        isNonzero = ((int*)matrix->entryValues)[e] != 0;
      else
      {
        assert(matrixType == sizeof(char));
        isNonzero = ((char*)matrix->entryValues)[e] != 0;
      }
      if (isNonzero)
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
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row+1];
    for (size_t e = first; e < beyond; ++e)
    {
      bool isNonzero = false;
      if (matrixType == sizeof(double))
        isNonzero = ((double*)matrix->entryValues)[e] != 0.0;
      else if (matrixType == sizeof(int))
        isNonzero = ((int*)matrix->entryValues)[e] != 0;
      else if (matrixType == sizeof(char))
        isNonzero = ((char*)matrix->entryValues)[e] != 0;
      else
        assert("Invalid matrixType parameter." == 0);

      if (isNonzero)
      {
        int column = matrix->entryColumns[e];
        int columnNode = firstColumnNode + column;
        CMRdbgMsg(2, "Nonzero (%d,%d) with row node %d (%d neighbors missing) and column node %d (%d neighbors missing).\n",
          row, column, row, graphNodes[row].degree, columnNode, graphNodes[columnNode].degree);
        graphAdjacencies[graphNodes[row + 1].adjacencyStart - graphNodes[row].degree] = columnNode;
        graphNodes[row].degree--;
        graphAdjacencies[graphNodes[columnNode + 1].adjacencyStart - graphNodes[columnNode].degree] = row;
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

  *pnumComponents = countComponents;

#if defined(CMR_DEBUG)
  printf("DFS found %d components.\n", countComponents);
  for (int node = 0; node < numNodes; ++node)
  {
    printf("Node %d has component %d.\n", node, graphNodes[node].component);
  }
#endif

  /* Allocate component data. */
  CMR_CALL( CMRallocBlockArray(cmr, pcomponents, countComponents) );
  CMR_ONESUM_COMPONENT* components = *pcomponents;

  /* Compute sizes. */
  for (int comp = 0; comp < countComponents; ++comp)
  {
    components[comp].matrix = NULL;
    components[comp].transpose = NULL;
    components[comp].rowsToOriginal = NULL;
    components[comp].columnsToOriginal = NULL;
    CMR_CALL( CMRchrmatCreate(cmr, (CMR_CHRMAT**) &components[comp].matrix, 0, 0, 0) );
  }

  for (int node = 0; node < numNodes; ++node)
  {
    int comp = graphNodes[node].component;
    int start = graphNodes[node].adjacencyStart;
    int end = graphNodes[node + 1].adjacencyStart;
    assert(comp >= 0);
    if (node < firstColumnNode)
    {
      components[comp].matrix->numRows++;
      components[comp].matrix->numNonzeros += end - start;
    }
    else
      components[comp].matrix->numColumns++;
  }

  /* Allocate memory */
  for (int comp = 0; comp < countComponents; ++comp)
  {
    CMR_MATRIX* compMatrix = components[comp].matrix;

#if defined(CMR_DEBUG)
    printf("Component %d has %dx%d matrix with %d nonzeros.\n", comp, compMatrix->numRows,
      compMatrix->numColumns, compMatrix->numNonzeros);
#endif

    CMR_CALL( CMRallocBlockArray(cmr, &components[comp].rowsToOriginal, compMatrix->numRows) );
    CMR_CALL( CMRallocBlockArray(cmr, &components[comp].columnsToOriginal, compMatrix->numColumns) );
    CMR_CALL( CMRreallocBlockArray(cmr, &compMatrix->rowSlice, compMatrix->numRows + 1) );
    if (compMatrix->numNonzeros > 0)
      CMR_CALL( CMRallocBlockArray(cmr, &compMatrix->entryColumns, compMatrix->numNonzeros) );

    if (targetType == sizeof(char))
    {
      if (compMatrix->numNonzeros > 0)
        CMR_CALL( CMRallocBlockArray(cmr, (char**) &compMatrix->entryValues, compMatrix->numNonzeros) );
      CMR_CALL( CMRchrmatCreate(cmr, (CMR_CHRMAT**) &components[comp].transpose,
        compMatrix->numColumns, compMatrix->numRows, compMatrix->numNonzeros) );
    }
    else if (targetType == sizeof(int))
    {
      if (compMatrix->numNonzeros > 0)
        CMR_CALL( CMRallocBlockArray(cmr, (int**) &compMatrix->entryValues, compMatrix->numNonzeros) );
      CMR_CALL( CMRintmatCreate(cmr, (CMR_INTMAT**) &components[comp].transpose,
        compMatrix->numColumns, compMatrix->numRows, compMatrix->numNonzeros) );
    }
    else
    {
      assert(targetType == sizeof(double));
      if (compMatrix->numNonzeros > 0)
        CMR_CALL( CMRallocBlockArray(cmr, (double**) &compMatrix->entryValues, compMatrix->numNonzeros) );
      CMR_CALL( CMRdblmatCreate(cmr, (CMR_DBLMAT**) &components[comp].transpose,
        compMatrix->numColumns, compMatrix->numRows, compMatrix->numNonzeros) );
    }
  }

  /* Fill mapping arrays. */
  for (int node = 0; node < numNodes; ++node)
  {
    int comp = graphNodes[node].component;
    int order = graphNodes[node].order;
    if (node < firstColumnNode)
      components[comp].rowsToOriginal[order] = node;
    else
      components[comp].columnsToOriginal[order] = node - firstColumnNode;
  }

#if defined(CMR_DEBUG)
  for (int comp = 0; comp < countComponents; ++comp)
  {
    printf("Component %d's rows map to original rows:", comp);
    for (int row = 0; row < components[comp].matrix->numRows; ++row)
      printf(" %d", components[comp].rowsToOriginal[row]);
    printf("\n");
    printf("Component %d's columns map to original columns:", comp);
    for (int column = 0; column < components[comp].matrix->numColumns; ++column)
      printf(" %d", components[comp].columnsToOriginal[column]);
    printf("\n");
  }
#endif

  /* We can now fill the matrices of each component. */
  for (int comp = 0; comp < countComponents; ++comp)
  {
    CMR_MATRIX* compTranspose = components[comp].transpose;

    /* Compute the slices in the transposed component matrix from the graph. */
    int countNonzeros = 0;
    for (size_t compColumn = 0; compColumn < compTranspose->numRows; ++compColumn)
    {
      int column = components[comp].columnsToOriginal[compColumn];
      int node = firstColumnNode + column;
      compTranspose->rowSlice[compColumn] = countNonzeros;
#if defined(CMR_DEBUG)
      printf("Component %d's column %d (row of transposed) starts at component entry %d.\n", comp, compColumn,
        countNonzeros);
#endif
      countNonzeros += graphNodes[node+1].adjacencyStart - graphNodes[node].adjacencyStart;
    }

    /* Fill the slices. To ensure that it is sorted, we iterate row-wise. */
    for (size_t compRow = 0; compRow < compTranspose->numColumns; ++compRow)
    {
      size_t row = components[comp].rowsToOriginal[compRow];
      size_t start = matrix->rowSlice[row];
      size_t end = matrix->rowSlice[row + 1];

      /* Iterate over all entries of that row. */
      for (size_t matrixEntry = start; matrixEntry < end; ++matrixEntry)
      {
        bool isNonzero = false;
        if (matrixType == sizeof(double))
          isNonzero = round(((double*)matrix->entryValues)[matrixEntry]) != 0;
        else if (matrixType == sizeof(int))
          isNonzero = ((int*)matrix->entryValues)[matrixEntry] != 0;
        else if (matrixType == sizeof(char))
          isNonzero = ((char*)matrix->entryValues)[matrixEntry] != 0;
        else
          assert("Invalid matrixType parameter." == 0);

        if (isNonzero)
        {
          size_t column = matrix->entryColumns[matrixEntry];
          int compColumn = graphNodes[firstColumnNode + column].order;
          int compEntry = compTranspose->rowSlice[compColumn];
#if defined(CMR_DEBUG)
          printf("Component %d contains matrix entry %d in at %d,%d.",
            comp, matrixEntry, row, column);
          printf(" It will be component entry %d at %d,%d\n", compEntry, compRow, compColumn);
#endif
          compTranspose->entryColumns[compEntry] = compRow;
          compTranspose->rowSlice[compColumn]++;
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

    /* Since we incremented the rowSlice for each nonzero, the array is shifted by one entry.
     * We restore this now. */
    for (int compColumn = compTranspose->numRows; compColumn > 0; --compColumn)
      compTranspose->rowSlice[compColumn] = compTranspose->rowSlice[compColumn-1];
    compTranspose->rowSlice[0] = 0;

#if defined(CMR_DEBUG)
    printf("Component %d's transpose:\n", comp);
    if (targetType == sizeof(double))
      CMRdblmatPrintDense(stdout, (CMR_DBLMAT*) compTranspose, '0', true);
    else if (targetType == sizeof(int))
      CMRintmatPrintDense(stdout, (CMR_INTMAT*) compTranspose, '0', true);
    else if (targetType == sizeof(char))
      CMRchrmatPrintDense(stdout, (CMR_CHRMAT*) compTranspose, '0', true);
#endif
  }

  /* We now create the row-wise representation from the column-wise one. */
  for (int comp = 0; comp < countComponents; ++comp)
  {
    CMR_MATRIX* compMatrix = components[comp].matrix;
    CMR_MATRIX* compTranspose = components[comp].transpose;

    /* Compute the slices in the component matrix from the graph. */
    int countNonzeros = 0;
    for (size_t compRow = 0; compRow < compMatrix->numRows; ++compRow)
    {
      int row = components[comp].rowsToOriginal[compRow];
      int node = row;
      compMatrix->rowSlice[compRow] = countNonzeros;
      countNonzeros += graphNodes[node+1].adjacencyStart - graphNodes[node].adjacencyStart;
    }

    /* Fill the slices. To ensure that it is sorted, we iterate column-wise. */
    for (size_t compColumn = 0; compColumn < compMatrix->numColumns; ++compColumn)
    {
      int start = compTranspose->rowSlice[compColumn];
      int end = compTranspose->rowSlice[compColumn + 1];

      /* Iterate over all entries of that column. */
      for (int compTransposeEntry = start; compTransposeEntry < end; ++compTransposeEntry)
      {
        int compRow = compTranspose->entryColumns[compTransposeEntry];
        int compMatrixEntry = compMatrix->rowSlice[compRow];
        compMatrix->entryColumns[compMatrixEntry] = compColumn;
        compMatrix->rowSlice[compRow]++;
#if defined(CMR_DEBUG)
        printf("Component matrix entry %d (%d,%d in component) copied from component transpose entry %d.\n",
          compMatrixEntry, compRow, compColumn, compTransposeEntry);
#endif
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

    /* Since we incremented the rowSlice for each nonzero, the array is shifted by one entry.
     * We restore this now. */
    for (int compRow = compMatrix->numRows; compRow > 0; --compRow)
      compMatrix->rowSlice[compRow] = compMatrix->rowSlice[compRow-1];
    compMatrix->rowSlice[0] = 0;
    bool isTranspose;
    if (targetType == sizeof(double))
    {
      CMR_CALL( CMRdblmatCheckTranspose(cmr, (CMR_DBLMAT*) components[comp].matrix,
        (CMR_DBLMAT*) components[comp].transpose, &isTranspose) ); 
    }
    else if (targetType == sizeof(int))
    {
      CMR_CALL( CMRintmatCheckTranspose(cmr, (CMR_INTMAT*) components[comp].matrix,
        (CMR_INTMAT*) components[comp].transpose, &isTranspose) ); 
    }
    else if (targetType == sizeof(char))
    {
      CMR_CALL( CMRchrmatCheckTranspose(cmr, (CMR_CHRMAT*) components[comp].matrix,
        (CMR_CHRMAT*) components[comp].transpose, &isTranspose) );
    }
    else
    {
      isTranspose = false;
    }
    assert(isTranspose);
  }

  /* Fill arrays for original matrix viewpoint. */
  if (rowsToComponents)
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      rowsToComponents[row] = graphNodes[row].component;
  }
  if (columnsToComponents)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnsToComponents[column] = graphNodes[firstColumnNode + column].component;
  }
  if (rowsToComponentRows)
  {
    for (size_t row = 0; row < matrix->numRows; ++row)
      rowsToComponentRows[row] = graphNodes[row].order;
  }
  if (columnsToComponentColumns)
  {
    for (size_t column = 0; column < matrix->numColumns; ++column)
      columnsToComponentColumns[column] = graphNodes[firstColumnNode + column].order;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &graphAdjacencies) );
  CMR_CALL( CMRfreeStackArray(cmr, &graphNodes) );

  return CMR_OKAY;
}
