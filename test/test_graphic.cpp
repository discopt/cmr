#include <gtest/gtest.h>

#include "common.h"
#include <tu/graphic.h>

TEST(Graphic, ToBinaryMatrix)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);
  
  TU_GRAPH* graph = NULL;
  TUgraphCreateEmpty(tu, &graph, 1, 1);

  TU_GRAPH_NODE v1 = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE v2 = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE v3 = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE v4 = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE v5 = TUgraphAddNode(tu, graph);

  TU_GRAPH_EDGE e1 = TUgraphAddEdge(tu, graph, v1, v1);
  TU_GRAPH_EDGE e2 = TUgraphAddEdge(tu, graph, v1, v2);
  TU_GRAPH_EDGE e3 = TUgraphAddEdge(tu, graph, v2, v3);
  TU_GRAPH_EDGE e4 = TUgraphAddEdge(tu, graph, v3, v4);
  TU_GRAPH_EDGE e5 = TUgraphAddEdge(tu, graph, v4, v5);
  TU_GRAPH_EDGE e6 = TUgraphAddEdge(tu, graph, v1, v4);
  TU_GRAPH_EDGE e7 = TUgraphAddEdge(tu, graph, v1, v5);
  TU_GRAPH_EDGE e8 = TUgraphAddEdge(tu, graph, v2, v4);
  
  TU_GRAPH_EDGE basis[4] = { e2, e3, e4, e5 };
  TU_GRAPH_EDGE cobasis[4] = { e1, e7, e6, e8 };

  TU_CHRMAT* matrix = NULL;
  TUconvertGraphToBinaryMatrix(tu, graph, &matrix, 4, basis, 4, cobasis);

  TUchrmatPrintDense(stdout, matrix, ' ', true);
  TUchrmatFree(tu, &matrix);

  TUgraphFree(tu, &graph);
  
  TUfreeEnvironment(&tu);
}
