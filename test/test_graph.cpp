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

TEST(Graph, EdgeList)
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

  FILE* graphFile = NULL;
  graphFile = fopen("graph.edgelist", "w");
  const char* nodeLabels[4] = { "a", "b", "c", "d" };
  CMR_ELEMENT edgeElements[6] = { -1, -2, -3, 1, 2, 3 };
  ASSERT_CMR_CALL( CMRgraphWriteEdgeList(cmr, graph, edgeElements, nodeLabels, graphFile) );
  fclose(graphFile);

  graphFile = fopen("graph.edgelist", "r");
  ASSERT_TRUE( graphFile );
  CMR_GRAPH* graph2 = NULL;
  CMR_ELEMENT* edgeElements2 = NULL;
  char** nodeLabels2 = NULL;
  ASSERT_CMR_CALL( CMRgraphCreateFromEdgeList(cmr, &graph2, &edgeElements2, &nodeLabels2, graphFile) );
  fclose(graphFile);

  // Write edge lists as a vector of strings and sort it for comparison.

  std::vector<std::string> edges;
  for (auto iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter); iter = CMRgraphEdgesNext(graph, iter))
  {
    auto edge = CMRgraphEdgesEdge(graph, iter);
    auto u = CMRgraphEdgeU(graph, edge);
    auto v = CMRgraphEdgeV(graph, edge);
    edges.push_back(std::string(nodeLabels[u]) + "-" + std::string(nodeLabels[v]) + "="
      + std::string(CMRelementString(edgeElements[edge], 0)));
  }
  std::sort(edges.begin(), edges.end());

  std::vector<std::string> edges2;
  for (auto iter = CMRgraphEdgesFirst(graph2); CMRgraphEdgesValid(graph2, iter); iter = CMRgraphEdgesNext(graph2, iter))
  {
    auto edge = CMRgraphEdgesEdge(graph2, iter);
    auto u = CMRgraphEdgeU(graph2, edge);
    auto v = CMRgraphEdgeV(graph2, edge);
    edges2.push_back(std::string(nodeLabels2[u]) + "-" + std::string(nodeLabels2[v]) + "="
      + std::string(CMRelementString(edgeElements2[edge], 0)));
  }
  std::sort(edges2.begin(), edges2.end());

  ASSERT_EQ(edges, edges2);

  free(edgeElements2);
  free(nodeLabels2);
  CMRgraphFree(cmr, &graph);
  CMRgraphFree(cmr, &graph2);

  CMRfreeEnvironment(&cmr);
}
