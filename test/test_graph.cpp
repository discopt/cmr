#include <gtest/gtest.h>

#include "common.h"
#include <tu/graph.h>

TEST(Graph, Modifications)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);
  
  TU_GRAPH* graph = NULL;
  TUgraphCreateEmpty(tu, &graph, 1, 1);

  TU_GRAPH_NODE a = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE b = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE c = TUgraphAddNode(tu, graph);
  TU_GRAPH_NODE d = TUgraphAddNode(tu, graph);

  TU_GRAPH_EDGE ab = TUgraphAddEdge(tu, graph, a, b);
  TU_GRAPH_EDGE ac = TUgraphAddEdge(tu, graph, a, c);
  TU_GRAPH_EDGE ad = TUgraphAddEdge(tu, graph, a, d);
  TU_GRAPH_EDGE bc = TUgraphAddEdge(tu, graph, b, c);
  TU_GRAPH_EDGE bd = TUgraphAddEdge(tu, graph, b, d);
  TU_GRAPH_EDGE cd = TUgraphAddEdge(tu, graph, c, d);

  ASSERT_EQ(TUgraphNumNodes(graph), 4);
  ASSERT_EQ(TUgraphNumEdges(graph), 6);

  int countNodes = 0;
  for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
    v = TUgraphNodesNext(graph, v))
  {
    ++countNodes;
  }
  ASSERT_EQ(countNodes, TUgraphNumNodes(graph));

  int countIncidentEdges = 0;
  for (TU_GRAPH_ITER i = TUgraphIncFirst(graph, b);
    TUgraphIncValid(graph, i); i = TUgraphIncNext(graph, i))
  {
    TU_GRAPH_EDGE e = TUgraphIncEdge(graph, i);
    ASSERT_GE(e, 0);
    ASSERT_LT(e, graph->memEdges);
    ++countIncidentEdges;
  }
  ASSERT_EQ(countIncidentEdges, 3);

  TUgraphDeleteEdge(tu, graph, bc);

  TUgraphDeleteNode(tu, graph, a);

  ASSERT_EQ(TUgraphNumNodes(graph), 3);
  ASSERT_EQ(TUgraphNumEdges(graph), 2);

  TUgraphFree(tu, &graph);
  
  TUfreeEnvironment(&tu);
}
