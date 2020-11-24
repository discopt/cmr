#include <tu/graphic.h>

#include "env_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "bixby_wagner.h"

#include <assert.h>
#include <limits.h>

static bool testGraphicness(TU* tu, int numComponents, TU_ONESUM_COMPONENT* components,
  TU_GRAPH** pgraph, TU_GRAPH_EDGE** pbasis, TU_GRAPH_EDGE** pcobasis,
  TU_SUBMAT** psubmatrix, int numRows, int numColumns)
{
  TU_GRAPH* graph = NULL;
  TU_GRAPH* componentGraph = NULL;
  if (pgraph)
  {
    TUgraphCreateEmpty(tu, pgraph, numRows + numComponents, numRows + numColumns);
    TUgraphCreateEmpty(tu, &componentGraph, numRows + numComponents, numRows + numColumns);
    graph = *pgraph;
  }
  TU_GRAPH_EDGE* componentBasis = NULL;
  if (pbasis)
  {
    TUallocBlockArray(tu, pbasis, numRows);
    TUallocBlockArray(tu, &componentBasis, numRows);

#ifndef NDEBUG /* We initialize with something far off to find bugs. */
    for (int r = 0; r < numRows; ++r)
      componentBasis[r] = INT_MIN;
#endif /* !NDEBUG */
  }
  TU_GRAPH_EDGE* componentCobasis = NULL;
  if (pcobasis)
  {
    TUallocBlockArray(tu, pcobasis, numColumns);
    TUallocBlockArray(tu, &componentCobasis, numColumns);

#ifndef NDEBUG /* We initialize with something far off to find bugs. */
    for (int c = 0; c < numColumns; ++c)
      componentCobasis[c] = INT_MIN;
#endif /* !NDEBUG */
  }

  bool isGraphic = true;
  for (int comp = 0; comp < numComponents; ++comp)
  {
    isGraphic = testGraphicnessBixbyWagner(tu, (TU_CHRMAT*)components[comp].matrix,
      (TU_CHRMAT*)components[comp].transpose, componentGraph, componentBasis, componentCobasis,
      psubmatrix);

    if (!isGraphic)
    {
      break;
    }

    if (componentGraph)
    {
      printf("testGraphicnessBixbyWagner returned the following graph:\n");
      TUgraphPrint(stdout, componentGraph);

      assert(graph);
      TU_GRAPH_NODE* componentNodesToNodes = NULL;
      TUallocBlockArray(tu, &componentNodesToNodes, TUgraphMemNodes(componentGraph));
      TU_GRAPH_EDGE* componentEdgesToEdges = NULL;
      TUallocBlockArray(tu, &componentEdgesToEdges, TUgraphMemEdges(componentGraph));

      printf("Copying %d nodes.\n", TUgraphNumNodes(componentGraph));
      for (TU_GRAPH_NODE v = TUgraphNodesFirst(componentGraph);
        TUgraphNodesValid(componentGraph, v); v = TUgraphNodesNext(componentGraph, v))
      {
        assert(v >= 0 && v < TUgraphMemNodes(componentGraph));
        componentNodesToNodes[v] = TUgraphAddNode(tu, graph);
        printf("component node %d is mapped to node %d.\n", v, componentNodesToNodes[v]);
      }

      for (TU_GRAPH_ITER i = TUgraphEdgesFirst(componentGraph);
        TUgraphEdgesValid(componentGraph, i); i = TUgraphEdgesNext(componentGraph, i))
      {
        TU_GRAPH_EDGE e = TUgraphEdgesEdge(componentGraph, i);
        assert(e >= 0 && e < TUgraphMemEdges(componentGraph));
        printf("Edge %d is {%d,%d} and mapped to {%d,%d}\n", e, TUgraphEdgeU(componentGraph, e),
          TUgraphEdgeV(componentGraph, e), componentNodesToNodes[TUgraphEdgeU(componentGraph, e)], componentNodesToNodes[TUgraphEdgeV(componentGraph, e)] );
        componentEdgesToEdges[e] = TUgraphAddEdge(tu, graph, componentNodesToNodes[
          TUgraphEdgeU(componentGraph, e)], componentNodesToNodes[TUgraphEdgeV(componentGraph, e)]);
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

      TUfreeBlockArray(tu, &componentEdgesToEdges);
      TUfreeBlockArray(tu, &componentNodesToNodes);

      TUgraphClear(tu, componentGraph);
    }
  }

  if (pcobasis)
    TUfreeBlockArray(tu, &componentCobasis);
  if (pbasis)
    TUfreeBlockArray(tu, &componentBasis);
  if (componentGraph)
    TUgraphFree(tu, &componentGraph);

  for (int comp = 0; comp < numComponents; ++comp)
  {
    TUchrmatFree(tu, (TU_CHRMAT**) &components[comp].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[comp].transpose);
    TUfreeBlockArray(tu, &components[comp].rowsToOriginal);
    TUfreeBlockArray(tu, &components[comp].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  return isGraphic;
}

bool TUtestGraphicnessChr(TU* tu, TU_CHRMAT* matrix, TU_GRAPH** pgraph,
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

  decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    NULL, NULL, NULL, NULL);

  /* Process all components. */

  return testGraphicness(tu, numComponents, components, pgraph, pbasis, pcobasis, psubmatrix,
    matrix->numRows, matrix->numColumns);
}
