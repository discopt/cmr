#ifndef TU_GRAPH_H
#define TU_GRAPH_H

#include <tu/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int TU_LISTGRAPH_NODE;
typedef int TU_LISTGRAPH_EDGE;
typedef int TU_LISTGRAPH_INCIDENT;

typedef struct
{
  int prev;
  int next;
  int firstOut;
} _TU_LISTGRAPH_NODE;

typedef struct
{
  int target;
  int prev;
  int next;
} _TU_LISTGRAPH_ARC;

typedef struct
{
  int numNodes;
  int memNodes;
  _TU_LISTGRAPH_NODE* nodes;
  int firstNode;
  int freeNode;

  int numEdges;
  int memEdges;
  _TU_LISTGRAPH_ARC* arcs;
  int freeEdge;
} TU_LISTGRAPH;

#define TUlistgraphNumNodes(graph) \
  ((graph)->numNodes)

#define TUlistgraphMemNodes(graph) \
  ((graph)->memNodes)

#define TUlistgraphNumEdges(graph) \
  ((graph)->numEdges)

#define TUlistgraphMemEdges(graph) \
  ((graph)->memEdges)

#define TUlistgraphEdgeU(graph, e) \
  ((graph)->arcs[2*e + 1].target)
#define TUlistgraphEdgeV(graph, e) \
  ((graph)->arcs[2*e].target)
  
void TUlistgraphCreateEmpty(
  TU* tu,                 /**< TU environment. */
  TU_LISTGRAPH** pgraph,  /**< Pointer to graph structure. */
  int memNodes,           /**< Allocate memory for this number of nodes. */
  int memEdges            /**< Allocate memory for this number of edges. */
);

void TUlistgraphFree(
  TU* tu,               /**< TU environment. */
  TU_LISTGRAPH** pgraph /**< Pointer to graph structure. */
);

/**
 * \brief Adds a node to a listgraph.
 * 
 * Adds a node to a listgraph.
 * 
 * \return Node structure of new node.
 */

TU_LISTGRAPH_NODE TUlistgraphAddNode(
  TU* tu,             /**< TU environment. */
  TU_LISTGRAPH* graph /**< Graph structure. */
);

/**
 * \brief Adds an edge to a listgraph.
 * 
 * Adds an edge to a listgraph.
 * 
 * \return Edge structure of new edge.
 */

TU_LISTGRAPH_EDGE TUlistgraphAddEdge(
  TU* tu,              /**< TU environment. */
  TU_LISTGRAPH* graph, /**< Graph structure. */
  TU_LISTGRAPH_NODE u, /**< One node of the edge. */
  TU_LISTGRAPH_NODE v  /**< Other node of the edge. */
);

void TUlistgraphDeleteNode(
  TU* tu,               /**< TU environment. */
  TU_LISTGRAPH* graph,  /**< Graph structure. */
  TU_LISTGRAPH_NODE v   /**< Node to be deleted. */
);

void TUlistgraphDeleteEdge(
  TU* tu,               /**< TU environment. */
  TU_LISTGRAPH* graph,  /**< Graph structure. */
  TU_LISTGRAPH_EDGE e   /**< Edge to be deleted. */
);

#define TUlistgraphNodesFirst(graph) \
  (graph->firstNode)
#define TUlistgraphNodesValid(graph, v) \
  (v >= 0)
#define TUlistgraphNodesNext(graph, v) \
  v = (graph)->nodes[v].next

#define TUlistgraphIncidentFirst(graph, v) \
  (graph)->nodes[v].firstOut
#define TUlistgraphIncidentValid(graph, i) \
  (i >= 0)
#define TUlistgraphIncidentNext(graph, i) \
  i = (graph)->arcs[i].next
#define TUlistgraphIncidentEdge(graph, i) \
  ((i)/2)
#define TUlistgraphIncidentSource(graph, i) \
  ((graph)->arcs[i ^ 1].target)
#define TUlistgraphIncidentTarget(graph, i) \
  ((graph)->arcs[i].target)

void TUlistgraphPrint(
  TU_LISTGRAPH* graph  /**< Graph structure. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPH_H */
