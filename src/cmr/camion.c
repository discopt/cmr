// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/camion.h>

#include "camion_internal.h"
#include "matrix_internal.h"
#include "block_decomposition.h"
#include "env_internal.h"

#if defined(CMR_DEBUG)
#include <cmr/graphic.h>
#endif /* CMR_DEBUG */

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

CMR_ERROR CMRcamionStatsInit(CMR_CAMION_STATISTICS* stats)
{
  assert(stats);

  stats->generalCount = 0;
  stats->generalTime = 0.0;
  stats->graphCount = 0;
  stats->graphTime = 0.0;
  stats->totalCount = 0;
  stats->totalTime = 0.0;

  return CMR_OKAY;
}

CMR_ERROR CMRcamionStatsPrint(FILE* stream, CMR_CAMION_STATISTICS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Camion signing:\n");
    prefix = "  ";
  }
  fprintf(stream, "%sgeneral: %lu in %f seconds\n", prefix, (unsigned long)stats->generalCount, stats->generalTime);
  fprintf(stream, "%sgraph: %lu in %f seconds\n", prefix, (unsigned long)stats->graphCount, stats->graphTime);
  fprintf(stream, "%stotal: %lu in %f seconds\n", prefix, (unsigned long)stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}

/**
 * \brief Graph node for BFS in signing algorithm.
 */

typedef struct
{
  int status;                   /**< \brief 0: not visited, 1: in queue, 2: processed */
  int predecessorNode;          /**< \brief Node number of predecessor. */
  signed char predecessorValue; /**< \brief Value of matrix entry of predecessor. */
  signed char targetValue;      /**< \brief Entry in current row if a target node, and 0 otherwise. */
} GRAPH_NODE;

CMR_ERROR CMRcamionComputeSignSequentiallyConnected(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< The matrix to be signed. */
  CMR_CHRMAT* transpose,    /**< The transpose of \p matrix. */
  bool change,              /**< Whether to modify the matrix. */
  char* pmodification,      /**< Pointer for storing which matrix was modified.*/
  CMR_SUBMAT** psubmatrix,  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
  double timeLimit          /**< Time limit to impose. */
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
    CMR_CALL( CMRcamionComputeSignSequentiallyConnected(cmr, transpose, matrix, change, pmodification, psubmatrix,
      timeLimit) );
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
  const size_t firstRowNode = matrix->numColumns;
  GRAPH_NODE* graphNodes = NULL;
  int* bfsQueue = NULL;
  int bfsQueueBegin = 0;
  int bfsQueueEnd = 0;
  clock_t time = clock();

  CMR_CALL(CMRallocStackArray(cmr, &graphNodes, matrix->numColumns + matrix->numRows));
  CMR_CALL(CMRallocStackArray(cmr, &bfsQueue, matrix->numColumns + matrix->numRows));

  /* Main loop iterates over the rows. */
  size_t clockRows = matrix->numRows / 100 + 1;
  for (size_t row = 1; row < matrix->numRows; ++row)
  {
    if ((row % clockRows) == 0 && ((clock() - time) * 1.0 / CLOCKS_PER_SEC > timeLimit))
    {
      CMRfreeStackArray(cmr, &bfsQueue);
      CMRfreeStackArray(cmr, &graphNodes);
      return CMR_ERROR_TIMEOUT;
    }
    
    CMRdbgMsg(2, "Before processing row %d:\n", row);
#if defined(CMR_DEBUG)
    CMRchrmatPrintDense(cmr, matrix, stdout, ' ', true);
#endif

    for (size_t v = 0; v < matrix->numColumns + matrix->numRows; ++v)
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
    size_t startNode = matrix->entryColumns[first];
    /* All columns of the row's nonzeros are target column nodes. */
    for (size_t e = first; e < beyond; ++e)
      graphNodes[matrix->entryColumns[e]].targetValue = matrix->entryValues[e];
    bfsQueue[0] = startNode;
    graphNodes[startNode].status = 1;
    bfsQueueBegin = 0;
    bfsQueueEnd = 1;

    while (bfsQueueBegin < bfsQueueEnd)
    {
      size_t currentNode = bfsQueue[bfsQueueBegin];
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
        for (size_t e = first; e < beyond; ++e)
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
              size_t pathNode = c;
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
    for (size_t v = 0; v < matrix->numColumns + row; ++v)
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
      for (size_t e = first; e < beyond; ++e)
      {
        int column = matrix->entryColumns[e];
        if (matrix->entryValues[e] != graphNodes[column].targetValue)
        {
          CMRdbgMsg(2, "Sign change at r%d,c%d.\n", row+1, column+1);
          matrix->entryValues[e] = graphNodes[column].targetValue;
          size_t entry = SIZE_MAX;
          CMR_CALL( CMRchrmatFindEntry(transpose, column, row, &entry) );
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
CMR_ERROR signCamion(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Matrix \f$ M \f$. */
  bool change,                  /**< Whether the signs of \f$ M \f$ shall be modified. */
  bool* pisCamionSigned,        /**< Pointer for storing whether \f$ M \f$ was already [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix,      /**< Pointer for storing a non-camion submatrix (may be \c NULL). */
  CMR_CAMION_STATISTICS* stats, /**< Statistics for the computation (may be \c NULL). */
  double timeLimit              /**< Time limit to impose. */
)
{
  assert(cmr);
  assert(matrix);
  assert(!psubmatrix || !*psubmatrix);

  clock_t totalClock = clock();

  size_t numBlocks;
  CMR_BLOCK* blocks = NULL;

  assert(CMRchrmatIsTernary(cmr, matrix, NULL));

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "signCamion:\n");
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  /* Decompose into 1-connected components. */

  CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) matrix, sizeof(signed char), sizeof(signed char), &numBlocks, &blocks, NULL,
    NULL, NULL, NULL) );

  if (pisCamionSigned)
    *pisCamionSigned = true;
  for (size_t comp = 0; comp < numBlocks; ++comp)
  {
    CMR_SUBMAT* compSubmatrix = NULL;

    CMRdbgMsg(2, "-> Block %d of size %dx%d\n", comp, blocks[comp].matrix->numRows,
      blocks[comp].matrix->numColumns);

    double remainingTime = timeLimit - ((clock() - totalClock) * 1.0 / CLOCKS_PER_SEC);
    char modified;
    CMR_CALL( CMRcamionComputeSignSequentiallyConnected(cmr, (CMR_CHRMAT*) blocks[comp].matrix,
      (CMR_CHRMAT*) blocks[comp].transpose, change, &modified,
      (psubmatrix && !*psubmatrix) ? &compSubmatrix : NULL, remainingTime) );

    CMRdbgMsg(2, "-> Block %d yields: %c\n", comp, modified ? modified : '0');

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
        compSubmatrix->rows[r] = blocks[comp].rowsToOriginal[compSubmatrix->rows[r]];
      for (size_t c = 0; c < compSubmatrix->numColumns; ++c)
        compSubmatrix->columns[c] = blocks[comp].columnsToOriginal[compSubmatrix->columns[c]];
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
      (CMR_CHRMAT*) blocks[comp].transpose :
      (CMR_CHRMAT*) blocks[comp].matrix;

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
        size_t row = blocks[comp].rowsToOriginal[compRow];
        size_t column = blocks[comp].columnsToOriginal[compColumn];

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

  for (size_t c = 0; c < numBlocks; ++c)
  {
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &blocks[c].matrix);
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &blocks[c].transpose);
    CMRfreeBlockArray(cmr, &blocks[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &blocks[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &blocks);

  if (stats)
  {
    double time = (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
    stats->generalCount++;
    stats->generalTime += time;
    stats->totalCount++;
    stats->totalTime += time;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRcamionTestSigns(CMR* cmr, CMR_CHRMAT* matrix, bool* pisCamionSigned, CMR_SUBMAT** psubmatrix,
  CMR_CAMION_STATISTICS* stats, double timeLimit)
{
  return signCamion(cmr, matrix, false, pisCamionSigned, psubmatrix, stats, timeLimit);
}

CMR_ERROR CMRcamionComputeSigns(CMR* cmr, CMR_CHRMAT* matrix, bool* pwasCamionSigned, CMR_SUBMAT** psubmatrix,
  CMR_CAMION_STATISTICS* stats, double timeLimit)
{
  return signCamion(cmr, matrix, true, pwasCamionSigned, psubmatrix, stats, timeLimit);
}

typedef struct
{
  CMR_ELEMENT element;        /**< Row or column of this edge. */
  CMR_GRAPH_EDGE causingEdge; /**< Another edge that triggered the fixation of this edge (-1 if not fixed). */
} OrientationSearchEdgeData;

typedef enum
{
  UNKNOWN = 0,    /**< \brief The node was not considered by the shortest-path, yet. */
  SEEN = 1,       /**< \brief Some path to the node is known. */
  COMPLETED = 2,  /**< \brief The shortest path to the node is known. */
  BASIC = 3,      /**< \brief The rootEdge of that node belongs to the spanning forest. */
} OrientationSearchStage;

typedef struct
{
  OrientationSearchStage stage; /**< Stage in BFS. */
  CMR_GRAPH_NODE predecessor;   /**< Predecessor node (\c SIZE_MAX if not present). */
  CMR_GRAPH_EDGE edge;          /**< Edge connecting to predecessor node. */
  size_t distance;              /**< Combinatorial distance to the BFS root. */
  int8_t sign;                  /**< Sign of this tree edge with respect to current column. */
  bool fixed;                   /**< Whether the orientation of this edge is already fixed. */
} OrientationSearchNodeData;

static
CMR_ERROR constructNonCamionSubmatrix(
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_GRAPH* cograph,                   /**< Cograph we consider. */
  OrientationSearchEdgeData* edgeData,  /**< Edge data array. */
  CMR_GRAPH_EDGE conflictEdge1,         /**< First edge of conflict. */
  CMR_GRAPH_EDGE conflictEdge2,         /**< Second edge of conflict. */
  CMR_SUBMAT** psubmatrix               /**< Pointer for storing the submatrix. */
)
{
  assert(cmr);
  assert(edgeData);
  assert(psubmatrix);

  CMRdbgMsg(6, "Extracting a non-Camion submatrix from an orientation conflict.\n");

  /* Trace back from first edge. */

  CMR_GRAPH_EDGE* trace1 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &trace1, CMRgraphMemEdges(cograph)) );
  size_t sizeTrace1 = 0;
  CMR_GRAPH_EDGE reason = conflictEdge1;
  while (reason != -1)
  {
    CMRdbgMsg(8, "1: Reason edge for edge %d is edge %d.\n", reason, edgeData[reason].causingEdge);
    trace1[sizeTrace1++] = reason;
    reason = edgeData[reason].causingEdge;
  }

  /* Trace back from second edge. */

  CMR_GRAPH_EDGE* trace2 = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &trace2, CMRgraphMemEdges(cograph)) );
  size_t sizeTrace2 = 0;
  reason = conflictEdge2;
  while (reason != -1)
  {
    CMRdbgMsg(8, "2: Reason edge for edge %d is edge %d.\n", reason, edgeData[reason].causingEdge);
    trace2[sizeTrace2++] = reason;
    reason = edgeData[reason].causingEdge;
  }

  /* Remove common end of traces. */
  while (sizeTrace1 > 0 && sizeTrace2 > 0 && trace1[sizeTrace1-1] == trace2[sizeTrace2-1])
  {
    --sizeTrace1;
    --sizeTrace2;
  }
  ++sizeTrace1; /* The last joint element must also participate (once). */

  assert((sizeTrace1 + sizeTrace2) % 2 == 0);

  size_t size = (sizeTrace1 + sizeTrace2) / 2;
  CMR_CALL( CMRsubmatCreate(cmr, size, size, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  submatrix->numRows = 0;
  submatrix->numColumns = 0;
  for (size_t i = 0; i < sizeTrace1; ++i)
  {
    CMR_ELEMENT element = edgeData[trace1[i]].element;
    if (CMRelementIsRow(element))
      submatrix->rows[submatrix->numRows++] = CMRelementToRowIndex(element);
    else
      submatrix->columns[submatrix->numColumns++] = CMRelementToColumnIndex(element);
  }
  for (size_t i = 0; i < sizeTrace2; ++i)
  {
    CMR_ELEMENT element = edgeData[trace2[i]].element;
    if (CMRelementIsRow(element))
      submatrix->rows[submatrix->numRows++] = CMRelementToRowIndex(element);
    else
      submatrix->columns[submatrix->numColumns++] = CMRelementToColumnIndex(element);
  }

  CMR_CALL( CMRfreeStackArray(cmr, &trace2) );
  CMR_CALL( CMRfreeStackArray(cmr, &trace1) );

  return CMR_OKAY;
}

CMR_ERROR CMRcamionCographicOrient(CMR* cmr, CMR_CHRMAT* matrix, CMR_GRAPH* cograph, CMR_GRAPH_EDGE* forestEdges,
  CMR_GRAPH_EDGE* coforestEdges, bool* arcsReversed, bool* pisCamionSigned, CMR_SUBMAT** psubmatrix,
  CMR_CAMION_STATISTICS* stats)
{
  assert(cmr);
  assert(matrix);
  assert(cograph);
  assert(forestEdges);
  assert(coforestEdges);
  assert(arcsReversed);
  assert(pisCamionSigned);
  assert(!psubmatrix || !*psubmatrix);

  clock_t totalClock = clock();

#if defined(CMR_DEBUG)

  CMRdbgMsg(2, "CMRcamionCographicOrient for matrix:\n");
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
  CMRdbgMsg(4, "Given cograph:\n");
  CMRgraphPrint(cograph, stdout);
  CMRdbgMsg(4, "Forest edges (in order):");
  for (size_t i = 0; i < CMRgraphNumNodes(cograph)-1; ++i)
    CMRdbgMsg(0, " e%d", forestEdges[i]);
  CMRdbgMsg(0, "\n");
  CMRdbgMsg(4, "Coforest edges (in order):");
  for (size_t i = 0; i < CMRgraphNumEdges(cograph) - CMRgraphNumNodes(cograph) + 1; ++i)
    CMRdbgMsg(0, " e%d", coforestEdges[i]);
  CMRdbgMsg(0, "\n");

  CMR_CHRMAT* graph_transpose = NULL;
  bool isCorrectForest;
  CMR_CALL( CMRgraphicComputeMatrix(cmr, cograph, NULL, &graph_transpose, CMRgraphNumNodes(cograph)-1,
    forestEdges, CMRgraphNumEdges(cograph) - CMRgraphNumNodes(cograph) + 1, coforestEdges, &isCorrectForest) );

  bool same = matrix->numNonzeros == graph_transpose->numNonzeros;
  if (same)
  {
    for (size_t row = 0; row <= matrix->numRows; ++row)
    {
      if (matrix->rowSlice[row] != graph_transpose->rowSlice[row])
      {
        same = false;
        break;
      }
    }
  }
  if (same)
  {
    for (size_t e = 0; e < matrix->numNonzeros; ++e)
    {
      if (matrix->entryColumns[e] != graph_transpose->entryColumns[e])
        same = false;
    }
  }
  if (!same)
  {
    CMRdbgMsg(4, "Error: the given cograph represents the following *different* matrix:\n");
    CMR_CALL( CMRchrmatPrintDense(cmr, graph_transpose, stdout, '0', true) );
    assert(same);
  }
  else
  {
    CMRdbgMsg(4, "-> Confirmed that the given graph/tree pair has the same representation matrix up to signs.\n");
  }

  CMR_CALL( CMRchrmatFree(cmr, &graph_transpose) );

#endif /* CMR_DEBUG */

  /* Decompose into blocks. */
  size_t numBlocks;
  CMR_BLOCK* blocks = NULL;
  CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) matrix, sizeof(signed char), sizeof(signed char), &numBlocks, &blocks, NULL,
    NULL, NULL, NULL) );

  /* Allocate and initialize auxiliary data for nodes. */
  OrientationSearchNodeData* nodeData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodeData, CMRgraphMemNodes(cograph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(cograph); CMRgraphNodesValid(cograph, v);
    v = CMRgraphNodesNext(cograph, v))
  {
    nodeData[v].stage = UNKNOWN;
    nodeData[v].fixed = false;
    nodeData[v].predecessor = -1;
    nodeData[v].distance = 0;
    nodeData[v].sign = 0;
    nodeData[v].edge = -1;
  }

  /* Allocate and initialize auxiliary data for edges. */
  OrientationSearchEdgeData* edgeData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &edgeData, CMRgraphMemEdges(cograph)) );
  CMRassertStackConsistency(cmr);
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(cograph); CMRgraphEdgesValid(cograph, i);
    i = CMRgraphEdgesNext(cograph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(cograph, i);
    edgeData[e].causingEdge = -1;
    arcsReversed[e] = false;
  }
  for (size_t row = 0; row < matrix->numRows; ++row)
    edgeData[coforestEdges[row]].element = CMRrowToElement(row);
  for (size_t column = 0; column < matrix->numColumns; ++column)
    edgeData[forestEdges[column]].element = CMRcolumnToElement(column);

  /* Allocate and initialize a queue for BFS. */
  CMR_GRAPH_NODE* queue = NULL;
  size_t queueFirst;
  size_t queueBeyond;
  CMR_CALL(CMRallocStackArray(cmr, &queue, matrix->numColumns + matrix->numRows));
  CMRassertStackConsistency(cmr);

  /* Process each block separately. */
  for (size_t b = 0; b < numBlocks; ++b)
  {
    CMR_CHRMAT* blockMatrix = (CMR_CHRMAT*) blocks[b].matrix;

#if defined(CMR_DEBUG)
    CMRdbgMsg(2, "Processing block #%zu of %zu.\n", b, numBlocks);
    for (size_t row = 0; row < blockMatrix->numRows; ++row)
      CMRdbgMsg(4, "Component row %zu corresponds to original row %zu.\n", row, blocks[b].rowsToOriginal[row]);
    for (size_t column = 0; column < blockMatrix->numColumns; ++column)
    {
      CMRdbgMsg(4, "Component column %zu corresponds to original column %zu.\n", column,
        blocks[b].columnsToOriginal[column]);
    }
    CMR_CALL( CMRchrmatPrintDense(cmr, blockMatrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

    /* If there are no nonzeros then also nothing must be signed. */
    if (blockMatrix->numNonzeros == 0)
      continue;

    assert(blockMatrix->numRows > 0);
    assert(blockMatrix->numColumns > 0);

    /* Run BFS on the component of the graph induced by this matrix block.
     * We use some node from one of the columns as a starting node. */
    size_t componentColumn = blocks[b].columnsToOriginal[0];
    CMR_GRAPH_EDGE e = forestEdges[componentColumn];
    CMR_GRAPH_NODE start = CMRgraphEdgeU(cograph, e);
    CMRdbgMsg(4, "Starting BFS at node %d.\n", start);
    queue[0] = start;
    queueFirst = 0;
    queueBeyond = 1;
    assert(nodeData[start].stage == UNKNOWN);
    nodeData[start].stage = SEEN;

    /* Process BFS queue until it is empty. */
    while (queueFirst < queueBeyond)
    {
      CMR_GRAPH_NODE v = queue[queueFirst];
      ++queueFirst;
      CMRdbgMsg(6, "Processing node %d.\n", v);
      nodeData[v].stage = COMPLETED;
      for (CMR_GRAPH_ITER i = CMRgraphIncFirst(cograph, v); CMRgraphIncValid(cograph, i);
        i = CMRgraphIncNext(cograph, i))
      {
        assert(CMRgraphIncSource(cograph, i) == v);
        CMR_GRAPH_NODE w = CMRgraphIncTarget(cograph, i);

        /* Skip if already completed. */
        if (nodeData[w].stage == COMPLETED)
          continue;

        CMR_GRAPH_EDGE e = CMRgraphIncEdge(cograph, i);

        /* We skip cotree edges. */
        if (CMRelementIsRow(edgeData[e].element))
          continue;

        if (nodeData[w].stage == UNKNOWN)
        {
          CMRdbgMsg(6, "Found new node via tree arc (%d,%d).\n", v, w);
          nodeData[w].stage = SEEN;
          nodeData[w].predecessor = v;
          nodeData[w].distance = nodeData[v].distance + 1;
          nodeData[w].edge = e;
          queue[queueBeyond] = w;
          ++queueBeyond;
        }
      }
    }

    /* We now go through the rows of the matrix and inspect the signs. */
    for (size_t componentRow = 0; componentRow < blockMatrix->numRows; ++componentRow)
    {
      size_t row = blocks[b].rowsToOriginal[componentRow];

      CMR_GRAPH_EDGE rowEdge = coforestEdges[row];
      CMR_GRAPH_NODE s = CMRgraphEdgeU(cograph, rowEdge);
      CMR_GRAPH_NODE t = CMRgraphEdgeV(cograph, rowEdge);

      CMRdbgMsg(4, "Inspecting signs of row r%zu corresponding to %d={%d,%d}.\n", row+1, rowEdge, s, t);

      size_t first = matrix->rowSlice[row];
      size_t beyond = matrix->rowSlice[row + 1];
      size_t minDistance = SIZE_MAX; /* The depth in the BFS tree that the s-r and t-r paths have in common. */
      for (size_t entry = first; entry < beyond; ++entry)
      {
        CMRdbgMsg(6, "Entry #%zu is in column c%zu with value %d.\n", entry, matrix->entryColumns[entry]+1,
          matrix->entryValues[entry]);

        CMR_GRAPH_EDGE rowEdge = forestEdges[matrix->entryColumns[entry]];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(cograph, rowEdge);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(cograph, rowEdge);
        if (nodeData[v].predecessor == u)
        {
          /* (u,v) */
          if (nodeData[u].distance < minDistance)
            minDistance = nodeData[u].distance;
          nodeData[v].sign = matrix->entryValues[entry];
        }
        else
        {
          /* (v,u) */
          assert(nodeData[u].predecessor == v);
          if (nodeData[v].distance < minDistance)
            minDistance = nodeData[v].distance;
          nodeData[u].sign = matrix->entryValues[entry];
        }
      }

      CMRdbgMsg(6, "Minimum distance is %d.\n", minDistance);

      /* Follow s-r path up to minDistance. If we encounter a fixed edge, then we decide whether we have to revert the
       * row edge. */
      CMR_GRAPH_NODE v = s;
      bool foundFixed = false;
      bool reversedRowEdge = false;
      while (nodeData[v].distance > minDistance)
      {
        if (nodeData[v].fixed)
        {
          int8_t currentSign = CMRgraphEdgeU(cograph, nodeData[v].edge) == v ? 1 : -1;
          if (arcsReversed[nodeData[v].edge])
            currentSign *= -1;
          foundFixed = true;
          reversedRowEdge = currentSign != nodeData[v].sign;
          edgeData[rowEdge].causingEdge = nodeData[v].edge;
          break;
        }
        v = nodeData[v].predecessor;
      }

      if (!foundFixed)
      {
        /* Since we were not successful with the s-r path, we now follow the t-r path up to minDistance. Again, if we
         * encounter a fixed edge, then we decide whether we have to revert the column edge. */
        v = t;
        while (nodeData[v].distance > minDistance)
        {
          if (nodeData[v].fixed)
          {
            int8_t currentSign = CMRgraphEdgeU(cograph, nodeData[v].edge) == v ? -1 : 1;
            if (arcsReversed[nodeData[v].edge])
              currentSign *= -1;
            foundFixed = true;
            reversedRowEdge = currentSign != nodeData[v].sign;
            edgeData[rowEdge].causingEdge = nodeData[v].edge;
            break;
          }
          v = nodeData[v].predecessor;
        }
      }

      /* Store whether we reversed the row edge. */
      arcsReversed[rowEdge] = reversedRowEdge;
      CMRdbgMsg(6, "Found a fixed tree edge: %s. Row edge reversed = %s\n", foundFixed ? "yes" : "no",
        reversedRowEdge ? "yes" : "no");

      /* Again we follow the s-r path up to minDistance to reorder the tree edges. */
      v = s;
      while (nodeData[v].distance > minDistance)
      {
        signed char currentSign = CMRgraphEdgeU(cograph, nodeData[v].edge) == v ? 1 : -1;

        if (reversedRowEdge)
          currentSign *= -1;

        bool shouldBeReversed = currentSign != nodeData[v].sign;
        if (nodeData[v].fixed)
        {
          if (arcsReversed[nodeData[v].edge] != shouldBeReversed)
          {
            CMRdbgMsg(6, "Found a contradiction in the orientation, i.e., the matrix is not a network matrix.\n");

            *pisCamionSigned = false;
            if (psubmatrix)
              CMR_CALL( constructNonCamionSubmatrix(cmr, cograph, edgeData, rowEdge, nodeData[v].edge, psubmatrix) );

            goto cleanup;
          }
        }
        else
        {
          arcsReversed[nodeData[v].edge] = shouldBeReversed;
          CMRdbgMsg(6, "Path from %d towards root: tree edge (%d,%d) is edge {%d,%d}; graph imposed sign"
            " (with row edge reverting) is %d; matrix sign is %d; reversed = %s\n", s, nodeData[v].predecessor, v,
            CMRgraphEdgeU(cograph, nodeData[v].edge), CMRgraphEdgeV(cograph, nodeData[v].edge),
            currentSign, nodeData[v].sign, shouldBeReversed ? "yes" : "no");
          nodeData[v].fixed = true;
          edgeData[nodeData[v].edge].causingEdge = rowEdge;
        }

#if !defined(NDEBUG)
        nodeData[v].sign = 0; /* For debugging we make all signs 0 again. */
#endif /* !NDEBUG */

        v = nodeData[v].predecessor;
      }

      /* Finally, we follow the t-r path up to minDistance to reorder the tree edges. */
      v = t;
      while (nodeData[v].distance > minDistance)
      {
        signed char currentSign = CMRgraphEdgeU(cograph, nodeData[v].edge) == v ? -1 : 1;
        if (reversedRowEdge)
          currentSign *= -1;

        bool shouldBeReversed = currentSign != nodeData[v].sign;
        if (nodeData[v].fixed)
        {
          if (arcsReversed[nodeData[v].edge] != shouldBeReversed)
          {
            CMRdbgMsg(6, "Found a contradiction in the orientation, i.e., the matrix is not a network matrix.\n");

            *pisCamionSigned = false;
            if (psubmatrix)
              CMR_CALL( constructNonCamionSubmatrix(cmr, cograph, edgeData, rowEdge, nodeData[v].edge, psubmatrix) );

            goto cleanup;
          }
        }
        else
        {
          arcsReversed[nodeData[v].edge] = shouldBeReversed;
          CMRdbgMsg(6, "Path from %d towards root: tree edge (%d,%d) is edge {%d,%d}; graph imposed sign"
            " (with row edge reverting) is %d; matrix sign is %d; reversed = %s\n", t, nodeData[v].predecessor, v,
            CMRgraphEdgeU(cograph, nodeData[v].edge), CMRgraphEdgeV(cograph, nodeData[v].edge),
            currentSign, nodeData[v].sign, shouldBeReversed ? "yes" : "no");
          nodeData[v].fixed = true;
          edgeData[nodeData[v].edge].causingEdge = rowEdge;
        }

#if !defined(NDEBUG)
          nodeData[v].sign = 0; /* For debugging we make all signs 0 again. */
#endif /* !NDEBUG */

        v = nodeData[v].predecessor;
      }
    }
  }

#if defined(CMR_DEBUG)
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(cograph); CMRgraphEdgesValid(cograph, i);
    i = CMRgraphEdgesNext(cograph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(cograph, i);
    CMRdbgMsg(2, "Edge %d={%d,%d} reversed = %s\n", e, CMRgraphEdgeU(cograph, e), CMRgraphEdgeV(cograph, e),
      arcsReversed[e] ? "yes" : "no");
  }
#endif /* CMR_DEBUG */

  *pisCamionSigned = true;

cleanup:

  /* Free search data. */
  CMRassertStackConsistency(cmr);
  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &edgeData) );
  CMR_CALL( CMRfreeStackArray(cmr, &nodeData) );
  CMRassertStackConsistency(cmr);

  /* Free memory of block decomposition. */
  for (size_t c = 0; c < numBlocks; ++c)
  {
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &blocks[c].matrix);
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &blocks[c].transpose);
    CMRfreeBlockArray(cmr, &blocks[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &blocks[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &blocks);

  if (stats)
  {
    double time = (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
    stats->graphCount++;
    stats->graphTime += time;
    stats->totalCount++;
    stats->totalTime += time;
  }

  return CMR_OKAY;
}

