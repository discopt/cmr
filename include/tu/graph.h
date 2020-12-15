#ifndef TU_GRAPH_H
#define TU_GRAPH_H

#include <tu/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int TU_GRAPH_NODE; /**< Reference to a node of \ref TU_GRAPH. */
typedef int TU_GRAPH_EDGE; /**< Reference to an edge of \ref TU_GRAPH. */
typedef int TU_GRAPH_ITER; /**< Reference to an edge iterator of \ref TU_GRAPH. */

typedef struct
{
  int prev;     /*< Next node in node list. */
  int next;     /*< Previous node in node list. */
  int firstOut; /*< First out-arc of this node. */
} TU_GRAPH_NODE_DATA;

typedef struct
{
  int target; /*< Target node of this arc. */
  int prev;   /*< Next arc in out-arc list of source node. */
  int next;   /*< Previous arc in out-arc list of source node. */
} TU_GRAPH_ARC_DATA;

typedef struct
{
  int numNodes;               /**< Number of nodes. */
  int memNodes;               /**< Number of nodes for which memory is allocated. */
  TU_GRAPH_NODE_DATA* nodes;  /**< Array containing node data. */
  int firstNode;              /**< Index of first node. */
  int freeNode;               /**< Beginning of free-list of nodes. */

  int numEdges;               /**< Number of edges. */
  int memEdges;               /**< Number of edges for which memory is allocated. */
  TU_GRAPH_ARC_DATA* arcs;    /**< Array containing arc data. */
  int freeEdge;               /**< Beginning of free-list of arc. */
} TU_GRAPH;

static inline
size_t TUgraphMemNodes(TU_GRAPH* graph)
{
  return graph->memNodes;
}

static inline
int TUgraphNumNodes(TU_GRAPH* graph)
{
  return graph->numNodes;
}

static inline
size_t TUgraphMemEdges(TU_GRAPH* graph)
{
  return graph->memEdges;
}

static inline
int TUgraphNumEdges(TU_GRAPH* graph)
{
  return graph->numEdges;
}

static inline
TU_GRAPH_NODE TUgraphEdgeU(TU_GRAPH* graph, TU_GRAPH_EDGE e)
{
  return graph->arcs[2*e+1].target;
}

static inline
TU_GRAPH_NODE TUgraphEdgeV(TU_GRAPH* graph, TU_GRAPH_EDGE e)
{
  return graph->arcs[2*e].target;
}

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

static inline
TU_GRAPH_NODE TUgraphNodesFirst(
  TU_GRAPH* graph /**< Graph structure. */
)
{
  return graph->firstNode; 
}

static inline
bool TUgraphNodesValid(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_NODE v   /**< Node. */
)
{
  return v >= 0;
}

static inline
TU_GRAPH_NODE TUgraphNodesNext(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_NODE v   /**< Node. */
)
{
  return graph->nodes[v].next;
}

static inline
TU_GRAPH_ITER TUgraphIncFirst(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_NODE v   /**< Node. */
)
{
  TU_GRAPH_ITER i = graph->nodes[v].firstOut;
  while (true)
  {
    if (i < 0)
      return -1;
    if ((graph->arcs[i].target != (graph)->arcs[i ^ 1].target) || !(i & 0x1))
      return i;
    i = graph->arcs[i].next;
  }
}

static inline
bool TUgraphIncValid(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return i >= 0;
}

static inline
TU_GRAPH_ITER TUgraphIncNext(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  while (true)
  {
    i = graph->arcs[i].next;
    if (i < 0)
      return -1;
    if (((graph)->arcs[i].target != (graph)->arcs[i ^ 1].target) || !(i & 0x1))
      return i;
  }
}

static inline
TU_GRAPH_EDGE TUgraphIncEdge(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return i/2;
}

static inline
TU_GRAPH_NODE TUgraphIncSource(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return graph->arcs[i^1].target;
}

static inline
TU_GRAPH_NODE TUgraphIncTarget(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return graph->arcs[i].target;
}

/**
 * \brief Returns iterator of next edge in list of all edges.
 */

static inline
TU_GRAPH_ITER TUgraphEdgesNext(
  TU_GRAPH* graph,  /*< Graph. */
  TU_GRAPH_ITER i   /*< Current edge iterator. */
)
{
  while (true)
  {
    TU_GRAPH_ITER j = graph->arcs[i].next;
    while (j >= 0 && (j & 0x1))
      j = graph->arcs[j].next;
    if (j >= 0)
      return j;

    TU_GRAPH_NODE source = graph->arcs[i ^ 1].target;
    source = graph->nodes[source].next;
    while (true)
    {
      if (source < 0)
        return -1;
      i = graph->nodes[source].firstOut;
      if (i >= 0)
      {
        if (!(i & 0x1))
          return i;
        else
          break;
      }
      source = graph->nodes[source].next;
    }
  } 
}

static inline
TU_GRAPH_ITER TUgraphEdgesFirst(
  TU_GRAPH* graph /**< Graph structure. */
)
{
  if (graph->firstNode < 0)
    return -1;

  TU_GRAPH_ITER i = graph->nodes[graph->firstNode].firstOut;
  if (i & 0x1)
    i = TUgraphEdgesNext(graph, i);
  return i;
}

static inline
bool TUgraphEdgesValid(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return i >= 0;
}


static inline
TU_GRAPH_EDGE TUgraphEdgesEdge(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges. */
)
{
  return i/2;
}

void TUgraphPrint(
  FILE* stream,   /*< Stream. */
  TU_GRAPH* graph /*< Graph structure. */
);

/**
 * \brief Merges two nodes \p u and \p v.
 */

void TUgraphMergeNodes(
  TU* tu,           /*< TU environment. */
  TU_GRAPH* graph,  /*< Graph. */
  TU_GRAPH_NODE u,  /*< First node. */
  TU_GRAPH_NODE v   /*< Second node. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPH_H */
