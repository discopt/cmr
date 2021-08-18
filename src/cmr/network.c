// #define CMR_DEBUG /* Uncomment to debug network. */

#include <cmr/graphic.h>
#include <cmr/network.h>
#include <cmr/camion.h>

#include "graphic_internal.h"
#include "env_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "heap.h"
#include "sort.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

typedef enum
{
  UNKNOWN = 0,    /**< \brief The node was not considered by the shortest-path, yet. */
  SEEN = 1,       /**< \brief Some path to the node is known. */
  COMPLETED = 2,  /**< \brief The shortest path to the node is known. */
  BASIC = 3,      /**< \brief The rootEdge of that node belongs to the spanning forest. */
} DIJKSTRA_STAGE;

CMR_ERROR CMRcomputeNetworkMatrix(CMR* cmr, CMR_GRAPH* digraph, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose,
  bool* arcsReversed, int numForestArcs, CMR_GRAPH_EDGE* forestArcs, int numCoforestArcs,
  CMR_GRAPH_EDGE* coforestArcs, bool* pisCorrectForest)
{
  assert(cmr);
  assert(digraph);
  assert(pmatrix || ptranspose);
  assert(!pmatrix || !*pmatrix);
  assert(!ptranspose || !*ptranspose);

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRcomputeRepresentationMatrix(cmr, digraph, true, &transpose, arcsReversed, numForestArcs, forestArcs,
    numCoforestArcs, coforestArcs, pisCorrectForest) );

  if (pmatrix)
    CMR_CALL( CMRchrmatTranspose(cmr, transpose, pmatrix) );

  /* Return or free the transpose matrix. */
  if (ptranspose)
    *ptranspose = transpose;
  else
    CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}


typedef struct
{
  int forestIndex;
} NetworkEdgeData;

typedef struct
{
  DIJKSTRA_STAGE stage; /* Stage in BFS. */
  int predecessor;      /* Predecessor node. */
  CMR_GRAPH_EDGE edge;   /* Edge connecting to predecessor node. */
  int distance;         /* Combinatorial distance to the BFS root. */
  char sign;            /* Sign of this tree edge with respect to current column. */
  bool fixed;           /* Whether the orientation of this edge is already fixed. */
} NetworkNodeData;



CMR_ERROR CMRtestConetworkMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisConetwork, CMR_GRAPH** pdigraph,
  CMR_GRAPH_EDGE** pforestArcs, CMR_GRAPH_EDGE** pcoforestArcs, bool** parcsReversed, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pforestArcs || pdigraph);
  assert(!pcoforestArcs || pdigraph);
  assert(!parcsReversed || pdigraph);
  assert(pisConetwork);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRtestConetworkMatrix called for a %dx%d matrix \n", matrix->numRows,
    matrix->numColumns);
  CMR_CALL( CMRchrmatPrintDense(cmr, stdout, matrix, '0', true) );
#endif /* CMR_DEBUG */

  bool isCamionSigned;
  CMR_CALL( CMRtestCamionSigned(cmr, matrix, &isCamionSigned, psubmatrix) );
  if (!isCamionSigned)
  {
    *pisConetwork = false;
    return CMR_OKAY;
  }

  CMR_GRAPH_EDGE* forestEdges = NULL;
  CMR_GRAPH_EDGE* coforestEdges = NULL;
  CMR_CALL( CMRtestCographicMatrix(cmr, matrix, pisConetwork, pdigraph, &forestEdges, &coforestEdges, psubmatrix) );
  if (pforestArcs)
    *pforestArcs = forestEdges;
  if (pcoforestArcs)
    *pcoforestArcs = coforestEdges;
  if (!*pisConetwork || !pdigraph || !parcsReversed)
  {
    /* We have to free (co)forest information if the caller didn't ask for it. */
    if (!pforestArcs)
      CMR_CALL( CMRfreeBlockArray(cmr, &forestEdges) );
    if (!pcoforestArcs)
      CMR_CALL( CMRfreeBlockArray(cmr, &coforestEdges) );
    return CMR_OKAY;
  }

  /* We have to find out which edges are reversed. */
  CMR_GRAPH* graph = *pdigraph;
  CMRdbgMsg(0, "Matrix is graphic. Computing reversed edges.");
  CMR_CALL( CMRallocBlockArray(cmr, parcsReversed, CMRgraphMemEdges(graph)) );

#if defined(CMR_DEBUG)
  CMRgraphPrint(stdout, *pdigraph);
  for (int b = 0; b < matrix->numColumns; ++b)
    CMRdbgMsg(2, "Forest #%d is %d.\n", b, (*pforestArcs)[b]);
  for (int b = 0; b < matrix->numRows; ++b)
    CMRdbgMsg(2, "Coforest #%d is %d.\n", b, (*pcoforestArcs)[b]);
#endif /* CMR_DEBUG */

  /* Decompose into 1-connected components. */
  size_t numComponents;
  CMR_ONESUM_COMPONENT* components = NULL;
  CMR_CALL( decomposeOneSum(cmr, (CMR_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components, NULL,
    NULL, NULL, NULL) );

  /* Allocate and initialize auxiliary data for nodes. */
  NetworkNodeData* nodeData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodeData, CMRgraphMemNodes(graph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
  {
    nodeData[v].stage = UNKNOWN;
    nodeData[v].fixed = false;
    nodeData[v].predecessor = -1;
    nodeData[v].distance = 0;
    nodeData[v].sign = 0;
    nodeData[v].edge = -1;
  }

  /* Allocate and initialize auxiliary data for edges. */
  NetworkEdgeData* edgeData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &edgeData, CMRgraphMemEdges(graph)) );
  CMRassertStackConsistency(cmr);
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, i); i = CMRgraphEdgesNext(graph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, i);
    edgeData[e].forestIndex = -1;
    (*parcsReversed)[e] = false;
  }
  for (int b = 0; b < matrix->numColumns; ++b)
    edgeData[forestEdges[b]].forestIndex = b;

  /* Allocate and initialize a queue for BFS. */
  int* queue = NULL;
  int queueFirst;
  int queueBeyond;
  CMRallocStackArray(cmr, &queue, matrix->numColumns + matrix->numRows);
  CMRassertStackConsistency(cmr);

  /* Process 1-connected components of the (transposed) matrix. */
  for (int comp = 0; comp < numComponents; ++comp)
  {
    CMR_CHRMAT* componentMatrix = (CMR_CHRMAT*) components[comp].transpose;

#if defined(CMR_DEBUG)
    CMRdbgMsg(2, "Processing component #%d of %d.\n", comp, numComponents);
    for (int row = 0; row < componentMatrix->numRows; ++row)
      CMRdbgMsg(4, "Component row %d corresponds to original row %d.\n", row, components[comp].columnsToOriginal[row]);
    for (int column = 0; column < componentMatrix->numColumns; ++column)
      CMRdbgMsg(4, "Component column %d corresponds to original column %d.\n", column,
        components[comp].rowsToOriginal[column]);
    CMR_CALL( CMRchrmatPrintDense(cmr, stdout, componentMatrix, '0', true) );
#endif /* CMR_DEBUG */

    /* If there are no nonzeros then also no signs can be wrong. */
    if (componentMatrix->numNonzeros == 0)
      continue;

    assert(componentMatrix->numRows > 0);
    assert(componentMatrix->numColumns > 0);

    /* Run BFS on the component of the graph induced by this 1-connected matrix component.
     * We use some node from one of the rows as a starting node. */
    int componentRow = components[comp].columnsToOriginal[0];
    CMR_GRAPH_EDGE e = forestEdges[componentRow];
    CMR_GRAPH_NODE start = CMRgraphEdgeU(graph, e);
    CMRdbgMsg(4, "Starting BFS at node %d.\n", start);
    queue[0] = start;
    queueFirst = 0;
    queueBeyond = 1;
    assert(nodeData[start].stage == UNKNOWN);
    nodeData[start].stage = SEEN;
    while (queueFirst < queueBeyond)
    {
      CMR_GRAPH_NODE v = queue[queueFirst];
      ++queueFirst;
      CMRdbgMsg(6, "Processing node %d.\n", v);
      nodeData[v].stage = COMPLETED;
      for (CMR_GRAPH_ITER i = CMRgraphIncFirst(graph, v); CMRgraphIncValid(graph, i); i = CMRgraphIncNext(graph, i))
      {
        assert(CMRgraphIncSource(graph, i) == v);
        CMR_GRAPH_NODE w = CMRgraphIncTarget(graph, i);

        /* Skip if already completed. */
        if (nodeData[w].stage == COMPLETED)
          continue;

        CMR_GRAPH_EDGE e = CMRgraphIncEdge(graph, i);
        if (edgeData[e].forestIndex < 0)
          continue;

        if (nodeData[w].stage == UNKNOWN)
        {
          CMRdbgMsg(6, "Found new node via arc (%d,%d).\n", v, w);
          nodeData[w].stage = SEEN;
          nodeData[w].predecessor = v;
          nodeData[w].distance = nodeData[v].distance + 1;
          nodeData[w].edge = e;
          queue[queueBeyond] = w;
          ++queueBeyond;
        }
      }
    }

    /* We now go through the columns of the matrix and inspect the signs. */
    for (int componentColumn = 0; componentColumn < componentMatrix->numColumns; ++componentColumn)
    {
      int column = components[comp].rowsToOriginal[componentColumn];

      CMR_GRAPH_EDGE columnEdge = coforestEdges[column];
      CMR_GRAPH_NODE s = CMRgraphEdgeU(graph, columnEdge);
      CMR_GRAPH_NODE t = CMRgraphEdgeV(graph, columnEdge);

      CMRdbgMsg(4, "Inspecting signs of column %d corresponding to %d={%d,%d}.\n", column, columnEdge, s, t);

      int first = matrix->rowStarts[column];
      int beyond = column == matrix->numRows ? matrix->numNonzeros : matrix->rowStarts[column+1];
      int minDistance = INT_MAX; /* The depth in the BFS tree that the s-r and t-r paths have in common. */
      for (int entry = first; entry < beyond; ++entry)
      {
        CMRdbgMsg(6, "Entry %d is in row %d with value %d.\n", entry, matrix->entryColumns[entry],
          matrix->entryValues[entry]);

        CMR_GRAPH_EDGE rowEdge = forestEdges[matrix->entryColumns[entry]];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, rowEdge);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, rowEdge);
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
       * column edge. */
      CMR_GRAPH_NODE v = s;
      bool foundFixed = false;
      bool reversedColumnEdge = false;
      while (nodeData[v].distance > minDistance)
      {
        if (nodeData[v].fixed)
        {
          char currentSign = CMRgraphEdgeU(graph, nodeData[v].edge) == v ? 1 : -1;
          if ((*parcsReversed)[nodeData[v].edge])
            currentSign *= -1;
          foundFixed = true;
          reversedColumnEdge = currentSign != nodeData[v].sign;
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
            char currentSign = CMRgraphEdgeU(graph, nodeData[v].edge) == v ? -1 : 1;
            if ((*parcsReversed)[nodeData[v].edge])
              currentSign *= -1;
            foundFixed = true;
            reversedColumnEdge = currentSign != nodeData[v].sign;
            break;
          }
          v = nodeData[v].predecessor;
        }
      }
      (*parcsReversed)[columnEdge] = reversedColumnEdge;
      CMRdbgMsg(6, "Found a fixed tree edge: %s. Column edge reversed = %s\n", foundFixed ? "yes" : "no",
        reversedColumnEdge ? "yes" : "no");

      /* Again we follow the s-r path up to minDistance to reorder the tree edges. */
      v = s;
      while (nodeData[v].distance > minDistance)
      {
        char currentSign = CMRgraphEdgeU(graph, nodeData[v].edge) == v ? 1 : -1;

        if (reversedColumnEdge)
          currentSign *= -1;
        assert(!nodeData[v].fixed || (*parcsReversed)[nodeData[v].edge] == (currentSign != nodeData[v].sign));
        (*parcsReversed)[nodeData[v].edge] = currentSign != nodeData[v].sign;
        CMRdbgMsg(6, "Path from %d towards root: tree edge (%d,%d) is edge {%d,%d}; graph imposed sign (with column edge reverting) is %d; matrix sign is %d; reversed = %s\n",
          s, nodeData[v].predecessor, v, CMRgraphEdgeU(graph, nodeData[v].edge), CMRgraphEdgeV(graph, nodeData[v].edge),
          currentSign, nodeData[v].sign, (*parcsReversed)[nodeData[v].edge] ? "yes" : "no");
        nodeData[v].fixed = true;
#if !defined(NDEBUG)
        nodeData[v].sign = 0; /* For debugging we make all signs 0 again. */
#endif /* !NDEBUG */
        v = nodeData[v].predecessor;
      }
      /* Finally, we follow the t-r path up to minDistance to reorder the tree edges. */
      v = t;
      while (nodeData[v].distance > minDistance)
      {
        char currentSign = CMRgraphEdgeU(graph, nodeData[v].edge) == v ? -1 : 1;
        if (reversedColumnEdge)
          currentSign *= -1;
        assert(!nodeData[v].fixed || (*parcsReversed)[nodeData[v].edge] == (currentSign != nodeData[v].sign));
        (*parcsReversed)[nodeData[v].edge] = currentSign != nodeData[v].sign;
        CMRdbgMsg(6, "Path from %d towards root: tree edge (%d,%d) is edge {%d,%d}; graph imposed sign (with column edge reverting) is %d; matrix sign is %d; reversed = %s\n",
          t, nodeData[v].predecessor, v, CMRgraphEdgeU(graph, nodeData[v].edge), CMRgraphEdgeV(graph, nodeData[v].edge),
          currentSign, nodeData[v].sign, (*parcsReversed)[nodeData[v].edge] ? "yes" : "no");
        nodeData[v].fixed = true;
#if !defined(NDEBUG)
        nodeData[v].sign = 0; /* For debugging we make all signs 0 again. */
#endif /* !NDEBUG */
        v = nodeData[v].predecessor;
      }
    }
  }

#if defined(CMR_DEBUG)
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, i); i = CMRgraphEdgesNext(graph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, i);
    CMRdbgMsg(2, "Edge %d={%d,%d} reversed = %s\n", e, CMRgraphEdgeU(graph, e), CMRgraphEdgeV(graph, e),
      (*parcsReversed)[e] ? "yes" : "no");
  }
#endif /* CMR_DEBUG */

  CMRassertStackConsistency(cmr);
  CMR_CALL( CMRfreeStackArray(cmr, &queue) );
  CMR_CALL( CMRfreeStackArray(cmr, &edgeData) );
  CMR_CALL( CMRfreeStackArray(cmr, &nodeData) );
  CMRassertStackConsistency(cmr);

  /* Free memory of 1-sum decomposition. */
  for (int c = 0; c < numComponents; ++c)
  {
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
    CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &components);

  /* We have to free (co)forest information if the caller didn't ask for it. */
  if (!pforestArcs)
    CMR_CALL( CMRfreeBlockArray(cmr, &forestEdges) );
  if (!pcoforestArcs)
    CMR_CALL( CMRfreeBlockArray(cmr, &coforestEdges) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestNetworkMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisNetwork, CMR_GRAPH** pdigraph,
  CMR_GRAPH_EDGE** pforestArcs, CMR_GRAPH_EDGE** pcoforestArcs, bool** parcsReversed, CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(matrix);
  assert(pisNetwork);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pforestArcs || pdigraph);
  assert(!pcoforestArcs || pdigraph);
  assert(!parcsReversed || pdigraph);

  /* Create transpose of matrix. */
  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

  CMR_CALL( CMRtestConetworkMatrix(cmr, transpose, pisNetwork, pdigraph, pforestArcs, pcoforestArcs, parcsReversed,
    psubmatrix) );

  /* Transpose minimal non-conetwork matrix to become a minimal non-network matrix. */
  if (psubmatrix && *psubmatrix)
  {
    size_t* temp = (*psubmatrix)->rows;
    (*psubmatrix)->rows = (*psubmatrix)->columns;
    (*psubmatrix)->columns = temp;

    size_t n = (*psubmatrix)->numRows;
    (*psubmatrix)->numRows = (*psubmatrix)->numColumns;
    (*psubmatrix)->numColumns = n;
  }

  CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}

/**@}*/
