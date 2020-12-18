#include <gtest/gtest.h>

#include "common.h"
#include <tu/graph.h>

TEST(Graph, Modifications)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);
  
  TU_GRAPH* graph = NULL;
  TUgraphCreateEmpty(tu, &graph, 1, 1);

  TU_GRAPH_NODE a,b,c,d;
  TUgraphAddNode(tu, graph, &a);
  TUgraphAddNode(tu, graph, &b);
  TUgraphAddNode(tu, graph, &c);
  TUgraphAddNode(tu, graph, &d);

  TU_GRAPH_EDGE ab, ac, ad, bc, bd, cd;
  TUgraphAddEdge(tu, graph, a, b, &ab);
  TUgraphAddEdge(tu, graph, a, c, &ac);
  TUgraphAddEdge(tu, graph, a, d, &ad);
  TUgraphAddEdge(tu, graph, b, c, &bc);
  TUgraphAddEdge(tu, graph, b, d, &bd);
  TUgraphAddEdge(tu, graph, c, d, &cd);

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
