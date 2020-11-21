#include <gtest/gtest.h>

#include "common.h"
#include <tu/graph.h>

TEST(Graph, Modifications)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_LISTGRAPH* graph = NULL;
  TUlistgraphCreateEmpty(tu, &graph, 1, 1);

  TU_LISTGRAPH_NODE a = TUlistgraphAddNode(tu, graph);
  TU_LISTGRAPH_NODE b = TUlistgraphAddNode(tu, graph);
  TU_LISTGRAPH_NODE c = TUlistgraphAddNode(tu, graph);
  TU_LISTGRAPH_NODE d = TUlistgraphAddNode(tu, graph);

  TU_LISTGRAPH_EDGE ab = TUlistgraphAddEdge(tu, graph, a, b);
  TU_LISTGRAPH_EDGE ac = TUlistgraphAddEdge(tu, graph, a, c);
  TU_LISTGRAPH_EDGE ad = TUlistgraphAddEdge(tu, graph, a, d);
  TU_LISTGRAPH_EDGE bc = TUlistgraphAddEdge(tu, graph, b, c);
  TU_LISTGRAPH_EDGE bd = TUlistgraphAddEdge(tu, graph, b, d);
  TU_LISTGRAPH_EDGE cd = TUlistgraphAddEdge(tu, graph, c, d);

  ASSERT_EQ(TUlistgraphNumNodes(graph), 4);
  ASSERT_EQ(TUlistgraphNumEdges(graph), 6);

  int countNodes = 0;
  for (TU_LISTGRAPH_NODE v = TUlistgraphNodesFirst(graph); TUlistgraphNodesValid(graph, v);
    TUlistgraphNodesNext(graph, v))
  {
    ++countNodes;
  }
  ASSERT_EQ(countNodes, TUlistgraphNumNodes(graph));

  int countIncidentEdges = 0;
  for (TU_LISTGRAPH_INCIDENT i= TUlistgraphIncidentFirst(graph, b);
    TUlistgraphIncidentValid(graph, i); TUlistgraphIncidentNext(graph, i))
  {
    TU_LISTGRAPH_EDGE e = TUlistgraphIncidentEdge(graph, i);
    ++countIncidentEdges;
  }
  ASSERT_EQ(countIncidentEdges, 3);

  TUlistgraphDeleteEdge(tu, graph, bc);

  TUlistgraphDeleteNode(tu, graph, a);

  ASSERT_EQ(TUlistgraphNumNodes(graph), 3);
  ASSERT_EQ(TUlistgraphNumEdges(graph), 2);

  TUlistgraphFree(tu, &graph);
  
  TUfreeEnvironment(&tu);
}
