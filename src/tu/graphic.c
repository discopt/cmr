// #define TU_DEBUG_GRAPHIC /* Uncomment to debug graphic. */

#include <tu/graphic.h>

#include "env_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "heap.h"
#include <tu/tdec.h>

#include <assert.h>
#include <limits.h>

/**
 * \brief Tests a 1-sum decomposition of a matrix for graphicness.
 */

static TU_ERROR testGraphicness(
  TU* tu,                           /**< \ref TU environment. */
  int numComponents,                /**< Number of 1-connected components. */
  TU_ONESUM_COMPONENT* components,  /**< Array of 1-connected components. */
  bool* isGraphic,                  /**< Set to \c true iff all components are graphic. */
  TU_GRAPH** pgraph,                /**< Either \c NULL or a pointer for storing a graph represented by the matrix. */
  TU_GRAPH_EDGE** pbasis,           /**< Either \c NULL or a pointer for storing an array of (basis) edges corresponding to the rows. */
  TU_GRAPH_EDGE** pcobasis,         /**< Either \c NULL or a pointer for storing an array of (cobasis) edges corresponding to the columns . */
  TU_SUBMAT** psubmatrix,           /**< Either \c NULL or a pointer for storing a minimally non-graphic submatrix. */
  int numRows,                      /**< Number of rows of the decomposed matrix. */
  int numColumns                    /**< Number of columns of the decomposed matrix. */
)
{
  TU_GRAPH* graph = NULL;
  TU_GRAPH* componentGraph = NULL;
  if (pgraph)
  {
    TU_CALL( TUgraphCreateEmpty(tu, pgraph, numRows + numComponents, numRows + numColumns) );
    TU_CALL( TUgraphCreateEmpty(tu, &componentGraph, numRows + numComponents, numRows + numColumns) );
    graph = *pgraph;
  }
  TU_GRAPH_EDGE* componentBasis = NULL;
  if (pbasis)
  {
    TU_CALL( TUallocBlockArray(tu, pbasis, numRows) );
    TU_CALL( TUallocStackArray(tu, &componentBasis, numRows) );

#ifndef NDEBUG /* We initialize with something far off to find bugs. */
    for (int r = 0; r < numRows; ++r)
      componentBasis[r] = INT_MIN;
#endif /* !NDEBUG */
  }
  TU_GRAPH_EDGE* componentCobasis = NULL;
  if (pcobasis)
  {
    TU_CALL( TUallocBlockArray(tu, pcobasis, numColumns) );
    TU_CALL( TUallocStackArray(tu, &componentCobasis, numColumns) );

#ifndef NDEBUG /* We initialize with something far off to find bugs. */
    for (int c = 0; c < numColumns; ++c)
      componentCobasis[c] = INT_MIN;
#endif /* !NDEBUG */
  }

  *isGraphic = true;
  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_CALL( testGraphicnessTDecomposition(tu, (TU_CHRMAT*)components[comp].matrix,
      (TU_CHRMAT*)components[comp].transpose, isGraphic, componentGraph, componentBasis,
      componentCobasis, psubmatrix) );

    if (!*isGraphic)
      break;

    if (componentGraph)
    {
#if defined(TU_DEBUG_GRAPHIC)
      printf("testGraphicnessTDecomposition returned the following graph:\n");
      TUgraphPrint(stdout, componentGraph);
#endif /* TU_DEBUG_GRAPHIC */

      assert(graph);
      TU_GRAPH_NODE* componentNodesToNodes = NULL;
      TU_CALL( TUallocStackArray(tu, &componentNodesToNodes, TUgraphMemNodes(componentGraph)) );
      TU_GRAPH_EDGE* componentEdgesToEdges = NULL;
      TU_CALL( TUallocStackArray(tu, &componentEdgesToEdges, TUgraphMemEdges(componentGraph)) );

#if defined(TU_DEBUG_GRAPHIC)
      printf("Copying %d nodes.\n", TUgraphNumNodes(componentGraph));
#endif /* TU_DEBUG_GRAPHIC */
      for (TU_GRAPH_NODE v = TUgraphNodesFirst(componentGraph);
        TUgraphNodesValid(componentGraph, v); v = TUgraphNodesNext(componentGraph, v))
      {
        assert(v >= 0 && v < TUgraphMemNodes(componentGraph));
        TU_CALL( TUgraphAddNode(tu, graph, &componentNodesToNodes[v]) );
#if defined(TU_DEBUG_GRAPHIC)
        printf("component node %d is mapped to node %d.\n", v, componentNodesToNodes[v]);
#endif /* TU_DEBUG_GRAPHIC */
      }

      for (TU_GRAPH_ITER i = TUgraphEdgesFirst(componentGraph);
        TUgraphEdgesValid(componentGraph, i); i = TUgraphEdgesNext(componentGraph, i))
      {
        TU_GRAPH_EDGE e = TUgraphEdgesEdge(componentGraph, i);
        assert(e >= 0 && e < TUgraphMemEdges(componentGraph));
        TU_CALL( TUgraphAddEdge(tu, graph, componentNodesToNodes[TUgraphEdgeU(componentGraph, e)],
          componentNodesToNodes[TUgraphEdgeV(componentGraph, e)], &componentEdgesToEdges[e]) );
      }
      if (componentBasis)
      {
        for (int r = 0; r < components[comp].matrix->numRows; ++r)
        {
          assert(componentBasis[r] >= 0);
          assert(componentEdgesToEdges[componentBasis[r]] >= 0);
          (*pbasis)[components[comp].rowsToOriginal[r]] = componentEdgesToEdges[componentBasis[r]];
        }
      }
      if (componentCobasis)
      {
        for (int c = 0; c < components[comp].matrix->numColumns; ++c)
        {
          assert(componentCobasis[c] >= 0);
          assert(componentEdgesToEdges[componentCobasis[c]] >= 0);
          (*pcobasis)[components[comp].columnsToOriginal[c]] =
            componentEdgesToEdges[componentCobasis[c]];
        }
      }

      TU_CALL( TUfreeStackArray(tu, &componentEdgesToEdges) );
      TU_CALL( TUfreeStackArray(tu, &componentNodesToNodes) );

      TU_CALL( TUgraphClear(tu, componentGraph) );
    }
  }

  if (pcobasis)
    TU_CALL( TUfreeStackArray(tu, &componentCobasis) );
  if (pbasis)
    TU_CALL( TUfreeStackArray(tu, &componentBasis) );
  if (componentGraph)
    TU_CALL( TUgraphFree(tu, &componentGraph) );
  if (!isGraphic && graph)
  {
    TU_CALL( TUgraphFree(tu, &graph) );
    *pgraph = NULL;
  }

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TU_CALL( TUchrmatFree(tu, (TU_CHRMAT**) &components[comp].matrix) );
    TU_CALL( TUchrmatFree(tu, (TU_CHRMAT**) &components[comp].transpose) );
    TU_CALL( TUfreeBlockArray(tu, &components[comp].rowsToOriginal) );
    TU_CALL( TUfreeBlockArray(tu, &components[comp].columnsToOriginal) );
  }
  TU_CALL( TUfreeBlockArray(tu, &components) );

  return TU_OKAY;
}

TU_ERROR TUtestGraphicnessChr(TU* tu, TU_CHRMAT* matrix, bool* isGraphic, TU_GRAPH** pgraph,
  TU_GRAPH_EDGE** pbasis, TU_GRAPH_EDGE** pcobasis, TU_SUBMAT** psubmatrix)
{
  int numComponents;
  TU_ONESUM_COMPONENT* components = NULL;

  assert(tu);
  assert(matrix);
  assert(!pgraph || !*pgraph);
  assert(!pbasis || !*pbasis);
  assert(!pcobasis || !*pcobasis);
  assert(!psubmatrix || !*psubmatrix);

  /* Perform 1-sum decomposition. */

  TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents,
    &components, NULL, NULL, NULL, NULL) );

  /* Process all components. */

  TU_CALL( testGraphicness(tu, numComponents, components, isGraphic, pgraph, pbasis, pcobasis,
    psubmatrix, matrix->numRows, matrix->numColumns) );

  return TU_OKAY;
}

typedef enum
{
  UNKNOWN = 0,    /**< The node was not considered by the shortest-path, yet. */
  SEEN = 1,       /**< Some path to the node is known. */
  COMPLETED = 2,  /**< The shortest path to the node is known. */
  BASIC = 3,      /**< The rootEdge of that node belongs to the spanning forest. */
} STAGE;

/**
 * \brief Node information for shortest-path computation in \ref TUconvertGraphToBinaryMatrix.
 */

typedef struct
{
  STAGE stage;            /**< At which stage of the algorithm is this node? */
  int predecessor;        /**< Predecessor node in shortest-path branching, or -1 for a root. */
  TU_GRAPH_EDGE rootEdge; /**< The actual edge towards the predecessor, or -1 for a root./ */
} NodeData;

/**
 * \brief Comparator for sorting ints increasingly.
 */

int compareInt(const void* A, const void* B)
{
  int* a = (int*) A;
  int* b = (int*) B;
  return *a - *b;
}

TU_ERROR TUconvertGraphToBinaryMatrix(TU* tu, TU_GRAPH* graph, TU_CHRMAT** pmatrix,
  int numBasisEdges, TU_GRAPH_EDGE* basisEdges, int numCobasisEdges, TU_GRAPH_EDGE* cobasisEdges)
{
  assert(tu);
  assert(graph);
  assert(pmatrix);
  assert(!*pmatrix);
  assert(numBasisEdges == 0 || basisEdges);
  assert(numCobasisEdges == 0 || cobasisEdges);

  NodeData* nodeData = NULL;
  TU_CALL( TUallocStackArray(tu, &nodeData, TUgraphMemNodes(graph)) );
  TU_INTHEAP heap;
  TU_CALL( TUintheapInitStack(tu, &heap, TUgraphMemNodes(graph)) );
  int* lengths = NULL;
  TU_CALL( TUallocStackArray(tu, &lengths, TUgraphMemEdges(graph)) );
  for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
    v = TUgraphNodesNext(graph, v))
  {
    nodeData[v].stage = UNKNOWN;
  }
  for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i);
    i = TUgraphEdgesNext(graph, i))
  {
    TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
    lengths[e] = 1;
  }
  for (int b = 0; b < numBasisEdges; ++b)
  {
    printf("basis[%d] = %d\n", b, basisEdges[b]);
    fflush(stdout);
    lengths[basisEdges[b]] = 0;
  }

  /* Start Dijkstra's algorithm at each node. */
  int countComponents = 0;
  for (TU_GRAPH_NODE s = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, s);
    s = TUgraphNodesNext(graph, s))
  {
    if (nodeData[s].stage != UNKNOWN)
      continue;

    nodeData[s].predecessor = -1;
    nodeData[s].rootEdge = -1;
    ++countComponents;
    TUintheapInsert(&heap, s, 0);
    while (!TUintheapEmpty(&heap))
    {
      int distance = TUintheapMinimumValue(&heap);
      TU_GRAPH_NODE v = TUintheapExtractMinimum(&heap);
      nodeData[v].stage = COMPLETED;
      for (TU_GRAPH_ITER i = TUgraphIncFirst(graph, v); TUgraphIncValid(graph, i);
        i = TUgraphIncNext(graph, i))
      {
        assert(TUgraphIncSource(graph, i) == v);
        TU_GRAPH_NODE w = TUgraphIncTarget(graph, i);

        /* Skip if already completed. */
        if (nodeData[w].stage == COMPLETED)
          continue;

        TU_GRAPH_EDGE e = TUgraphIncEdge(graph, i);
        int newDistance = distance + lengths[e];
        if (newDistance < TUintheapGetValueInfinity(&heap, w))
        {
          nodeData[w].stage = SEEN;
          nodeData[w].predecessor = v;
          nodeData[w].rootEdge = e;
          TUintheapDecreaseInsert(&heap, w, newDistance);
        }
      }
    }
  }

  TU_CALL( TUfreeStackArray(tu, &lengths) );
  TU_CALL( TUintheapClearStack(tu, &heap) );

  /* Now nodeData[.].predecessor is an arborescence for each connected component. */

  TU_GRAPH_NODE* nodesRows = NULL; /* Non-root node v is mapped to row of edge {v,predecessor(v)}. */
  TU_CALL( TUallocStackArray(tu, &nodesRows, TUgraphMemNodes(graph)) );
  for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
    v = TUgraphNodesNext(graph, v))
  {
    nodesRows[v] = -1;
  }
  int numRows = 0;
  for (int i = 0; i < numBasisEdges; ++i)
  {
    TU_GRAPH_NODE u = TUgraphEdgeU(graph, basisEdges[i]);
    TU_GRAPH_NODE v = TUgraphEdgeV(graph, basisEdges[i]);
    if (nodeData[u].predecessor == v)
    {
      nodesRows[u] = numRows;
      ++numRows;
      nodeData[u].stage = BASIC;
    }
    else if (nodeData[v].predecessor == u)
    {
      nodesRows[v] = numRows;
      ++numRows;
      nodeData[v].stage = BASIC;
    }
  }
  if (numRows < TUgraphNumNodes(graph) - countComponents)
  {
    /* Some basis edges are missing. */
    for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
      v = TUgraphNodesNext(graph, v))
    {
      if (nodeData[v].predecessor >= 0 && nodeData[v].stage != BASIC)
      {
        nodesRows[v] = numRows;
        ++numRows;
        nodeData[v].stage = BASIC;
      }
    }
  }

  TU_CHRMAT* transposed = NULL;
  TU_CALL( TUchrmatCreate(tu, &transposed, TUgraphNumEdges(graph) - numRows, numRows,
    16 * numRows) );
  int numNonzeros = 0; /* Current number of nonzeros. transpose->numNonzeros is the memory. */
  int numColumns = 0;
  TU_GRAPH_EDGE* edgeColumns = NULL;
  TU_CALL( TUallocStackArray(tu, &edgeColumns, TUgraphMemEdges(graph)) );
  for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i);
    i = TUgraphEdgesNext(graph, i))
  {
    TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
    TU_GRAPH_NODE u = TUgraphEdgeU(graph, e);
    TU_GRAPH_NODE v = TUgraphEdgeV(graph, e);
    edgeColumns[TUgraphEdgesEdge(graph, i)] =
      (nodeData[u].rootEdge == e || nodeData[v].rootEdge == e) ? -1 : -2;
  }
  TU_GRAPH_NODE* uPath = NULL;
  TU_CALL( TUallocStackArray(tu, &uPath, numRows) );
  TU_GRAPH_NODE* vPath = NULL;
  TU_CALL( TUallocStackArray(tu, &vPath, numRows) );
  TU_GRAPH_ITER iter = TUgraphEdgesFirst(graph);
  int cobasicIndex = 0;
  while (TUgraphEdgesValid(graph, iter))
  {
    TU_GRAPH_EDGE e;
    if (cobasicIndex < numCobasisEdges)
    {
      e = cobasisEdges[cobasicIndex];
      ++cobasicIndex;
    }
    else
    {
      e = TUgraphEdgesEdge(graph, iter);
      iter = TUgraphEdgesNext(graph, iter);
    }

    if (edgeColumns[e] >= -1)
      continue;

    TU_GRAPH_NODE u = TUgraphEdgeU(graph, e);
    TU_GRAPH_NODE v = TUgraphEdgeV(graph, e);

    transposed->rowStarts[numColumns] = numNonzeros;
    edgeColumns[e] = numColumns;

    /* Enlarge space for nonzeros if necessary. */
    if (numNonzeros + numRows > transposed->numNonzeros)
      TU_CALL( TUchrmatChangeNumNonzeros(tu, transposed, 2 * transposed->numNonzeros) );

    /* Compute u-root path. */
    int uPathLength = 0;
    TU_GRAPH_NODE w = u;
    while (nodeData[w].predecessor != -1)
    {
      uPath[uPathLength] = w;
      ++uPathLength;
      w = nodeData[w].predecessor;
    }

    /* Compute v-root path. */
    int vPathLength = 0;
    w = v;
    while (nodeData[w].predecessor != -1)
    {
      vPath[vPathLength] = w;
      ++vPathLength;
      w = nodeData[w].predecessor;
    }

    /* Remove common part of u-root path and v-root path. */
    while (uPathLength > 0 && vPathLength > 0 && uPath[uPathLength-1] == vPath[vPathLength-1])
    {
      --uPathLength;
      --vPathLength;
    }

    for (int j = 0; j < uPathLength; ++j)
    {
      assert(nodesRows[uPath[j]] >= 0);
      transposed->entryColumns[numNonzeros] = nodesRows[uPath[j]];
      transposed->entryValues[numNonzeros] = 1;
      ++numNonzeros;
    }
    for (int j = 0; j < vPathLength; ++j)
    {
      assert(nodesRows[vPath[j]] >= 0);
      transposed->entryColumns[numNonzeros] = nodesRows[vPath[j]];
      transposed->entryValues[numNonzeros] = 1;
      ++numNonzeros;
    }
    qsort(&transposed->entryColumns[transposed->rowStarts[numColumns]], uPathLength + vPathLength,
      sizeof(int), compareInt);

    ++numColumns;
  }

  TU_CALL( TUfreeStackArray(tu, &vPath) );
  TU_CALL( TUfreeStackArray(tu, &uPath) );
  TU_CALL( TUfreeStackArray(tu, &edgeColumns) );

  transposed->rowStarts[numColumns] = numNonzeros;
  if (numNonzeros == 0 && transposed->numNonzeros > 0)
  {
    TU_CALL( TUfreeBlockArray(tu, &transposed->entryColumns) );
    TU_CALL( TUfreeBlockArray(tu, &transposed->entryValues) );
  }
  transposed->numNonzeros = numNonzeros;

  TU_CALL( TUchrmatTranspose(tu, transposed, pmatrix) );
  TU_CALL( TUchrmatFree(tu, &transposed) );

  /* We now process the nonbasic edges. */

  TU_CALL( TUfreeStackArray(tu, &nodesRows) );
  TU_CALL( TUfreeStackArray(tu, &nodeData) );

  return TU_OKAY;
}
