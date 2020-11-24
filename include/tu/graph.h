#ifndef TU_GRAPH_H
#define TU_GRAPH_H

#include <tu/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int TU_GRAPH_NODE;
typedef int TU_GRAPH_EDGE;
typedef int TU_GRAPH_ITER;

typedef struct
{
  int prev;
  int next;
  int firstOut;
} _TU_GRAPH_NODE;

typedef struct
{
  int target;
  int prev;
  int next;
} _TU_GRAPH_ARC;

typedef struct
{
  int numNodes;
  int memNodes;
  _TU_GRAPH_NODE* nodes;
  int firstNode;
  int freeNode;

  int numEdges;
  int memEdges;
  _TU_GRAPH_ARC* arcs;
  int freeEdge;
} TU_GRAPH;

#define TUgraphNumNodes(graph) \
  ((graph)->numNodes)

#define TUgraphMemNodes(graph) \
  ((graph)->memNodes)

#define TUgraphNumEdges(graph) \
  ((graph)->numEdges)

#define TUgraphMemEdges(graph) \
  ((graph)->memEdges)

#define TUgraphEdgeU(graph, e) \
  ((graph)->arcs[2*e + 1].target)
#define TUgraphEdgeV(graph, e) \
  ((graph)->arcs[2*e].target)
  
void TUgraphCreateEmpty(
  TU* tu,                 /**< TU environment. */
  TU_GRAPH** pgraph,  /**< Pointer to graph structure. */
  int memNodes,           /**< Allocate memory for this number of nodes. */
  int memEdges            /**< Allocate memory for this number of edges. */
);

void TUgraphFree(
  TU* tu,               /**< TU environment. */
  TU_GRAPH** pgraph /**< Pointer to graph structure. */
);

/**
 * \brief Removes all nodes and columns, keeping the memory.
 */

void TUgraphClear(
  TU* tu,             /**< TU environment. */
  TU_GRAPH* graph /**< Graph structure. */
);

/**
 * \brief Adds a node to a listgraph.
 * 
 * Adds a node to a listgraph.
 * 
 * \return Node structure of new node.
 */

TU_GRAPH_NODE TUgraphAddNode(
  TU* tu,             /**< TU environment. */
  TU_GRAPH* graph /**< Graph structure. */
);

/**
 * \brief Adds an edge to a listgraph.
 * 
 * Adds an edge to a listgraph.
 * 
 * \return Edge structure of new edge.
 */

TU_GRAPH_EDGE TUgraphAddEdge(
  TU* tu,              /**< TU environment. */
  TU_GRAPH* graph, /**< Graph structure. */
  TU_GRAPH_NODE u, /**< One node of the edge. */
  TU_GRAPH_NODE v  /**< Other node of the edge. */
);

void TUgraphDeleteNode(
  TU* tu,               /**< TU environment. */
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_NODE v   /**< Node to be deleted. */
);

void TUgraphDeleteEdge(
  TU* tu,               /**< TU environment. */
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_EDGE e   /**< Edge to be deleted. */
);

#define TUgraphNodesFirst(graph) \
  ((graph)->firstNode)
#define TUgraphNodesValid(graph, v) \
  (v >= 0)
#define TUgraphNodesNext(graph, v) \
  (graph)->nodes[v].next

TU_GRAPH_ITER TUgraphIncFirst(TU_GRAPH* graph, TU_GRAPH_NODE v);

#define TUgraphIncValid(graph, i) \
  (i >= 0)

TU_GRAPH_ITER TUgraphIncNext(TU_GRAPH* graph, TU_GRAPH_ITER e);

#define TUgraphIncEdge(graph, i) \
  ((i)/2)
#define TUgraphIncSource(graph, i) \
  ((graph)->arcs[i ^ 1].target)
#define TUgraphIncTarget(graph, i) \
  ((graph)->arcs[i].target)


#define TUgraphEdgesValid(graph, e) \
  (e >= 0)

TU_GRAPH_ITER TUgraphEdgesFirst(
  TU_GRAPH* graph
);

TU_GRAPH_ITER TUgraphEdgesNext(
  TU_GRAPH* graph,
  TU_GRAPH_ITER e
);

#define TUgraphEdgesEdge(graph, i) \
  ((i)/2)

void TUgraphPrint(
  FILE* stream,        /**< Stream. */
  TU_GRAPH* graph  /**< Graph structure. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPH_H */
