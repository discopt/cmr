#include <gtest/gtest.h>

#include "common.h"
#include <cmr/graph.h>

TEST(Graph, Modifications)
{
  CMR* cmr = NULL;
  CMRcreateEnvironment(&cmr);
  
  CMR_GRAPH* graph = NULL;
  CMRgraphCreateEmpty(cmr, &graph, 1, 1);

  CMR_GRAPH_NODE a,b,c,d;
  CMRgraphAddNode(cmr, graph, &a);
  CMRgraphAddNode(cmr, graph, &b);
  CMRgraphAddNode(cmr, graph, &c);
  CMRgraphAddNode(cmr, graph, &d);

  CMR_GRAPH_EDGE ab, ac, ad, bc, bd, cd;
  CMRgraphAddEdge(cmr, graph, a, b, &ab);
  CMRgraphAddEdge(cmr, graph, a, c, &ac);
  CMRgraphAddEdge(cmr, graph, a, d, &ad);
  CMRgraphAddEdge(cmr, graph, b, c, &bc);
  CMRgraphAddEdge(cmr, graph, b, d, &bd);
  CMRgraphAddEdge(cmr, graph, c, d, &cd);

  ASSERT_EQ(CMRgraphNumNodes(graph), 4UL);
  ASSERT_EQ(CMRgraphNumEdges(graph), 6UL);

  size_t countNodes = 0;
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v);
    v = CMRgraphNodesNext(graph, v))
  {
    ++countNodes;
  }
  ASSERT_EQ(countNodes, CMRgraphNumNodes(graph));

  int countIncidentEdges = 0;
  for (CMR_GRAPH_ITER i = CMRgraphIncFirst(graph, b);
    CMRgraphIncValid(graph, i); i = CMRgraphIncNext(graph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphIncEdge(graph, i);
    ASSERT_GE(e, 0);
    ASSERT_LT(e, (int) graph->memEdges);
    ++countIncidentEdges;
  }
  ASSERT_EQ(countIncidentEdges, 3);

  CMRgraphDeleteEdge(cmr, graph, bc);

  CMRgraphDeleteNode(cmr, graph, a);

  ASSERT_EQ(CMRgraphNumNodes(graph), 3UL);
  ASSERT_EQ(CMRgraphNumEdges(graph), 2UL);

  CMRgraphFree(cmr, &graph);
  
  CMRfreeEnvironment(&cmr);
}
