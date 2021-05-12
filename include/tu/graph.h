#ifndef TU_GRAPH_H
#define TU_GRAPH_H

#include <tu/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup Graph Graph
 *
 * Undirected graphs.
 *
 * @{
 */

typedef int TU_GRAPH_NODE; /**< \brief Reference to a node of \ref TU_GRAPH. */
typedef int TU_GRAPH_EDGE; /**< \brief Reference to an edge of \ref TU_GRAPH. */
typedef int TU_GRAPH_ITER; /**< \brief Reference to an edge iterator of \ref TU_GRAPH. */

typedef struct
{
  int prev;     /**< \brief Next node in node list. */
  int next;     /**< \brief Previous node in node list. */
  int firstOut; /**< \brief First out-arc of this node. */
} TU_GRAPH_NODE_DATA;

typedef struct
{
  int target; /**< \brief Target node of this arc. */
  int prev;   /**< \brief Next arc in out-arc list of source node. */
  int next;   /**< \brief Previous arc in out-arc list of source node. */
} TU_GRAPH_ARC_DATA;

typedef struct
{
  int numNodes;               /**< \brief Number of nodes. */
  int memNodes;               /**< \brief Number of nodes for which memory is allocated. */
  TU_GRAPH_NODE_DATA* nodes;  /**< \brief Array containing node data. */
  int firstNode;              /**< \brief Index of first node. */
  int freeNode;               /**< \brief Beginning of free-list of nodes. */

  int numEdges;               /**< \brief Number of edges. */
  int memEdges;               /**< \brief Number of edges for which memory is allocated. */
  TU_GRAPH_ARC_DATA* arcs;    /**< \brief Array containing arc data. */
  int freeEdge;               /**< \brief Beginning of free-list of arc. */
} TU_GRAPH;

/**
 * \brief Returns number of nodes for which memory is allocated.
 */

static inline
size_t TUgraphMemNodes(
  TU_GRAPH* graph /**< Graph. */
)
{
  return graph->memNodes;
}

/**
 * \brief Returns number of nodes.
 */

static inline
int TUgraphNumNodes(
  TU_GRAPH* graph /**< Graph. */
)
{
  return graph->numNodes;
}

/**
 * \brief Returns number of edges for which memory is allocated.
 */

static inline
size_t TUgraphMemEdges(
  TU_GRAPH* graph /**< Graph. */
)
{
  return graph->memEdges;
}

/**
 * \brief Returns number of edges.
 */

static inline
int TUgraphNumEdges(
  TU_GRAPH* graph /**< Graph. */
)
{
  return graph->numEdges;
}

/**
 * \brief Returns node u of edge \p e.
 */

static inline
TU_GRAPH_NODE TUgraphEdgeU(
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_EDGE e   /**< Edge of \p graph. */
)
{
  return graph->arcs[2*e+1].target;
}

/**
 * \brief Returns node v of edge \p e.
 */

static inline
TU_GRAPH_NODE TUgraphEdgeV(
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_EDGE e   /**< Edge of \p graph. */
)
{
  return graph->arcs[2*e].target;
}

/**
 * \brief Creates an empty graph.
 */

TU_ERROR TUgraphCreateEmpty(
  TU* tu,             /**< \ref TU environment. */
  TU_GRAPH** pgraph,  /**< Pointer for storing the graph. */
  int memNodes,       /**< Allocate memory for this number of nodes. */
  int memEdges        /**< Allocate memory for this number of edges. */
);

/**
 * \brief Frees a graph.
 */

TU_ERROR TUgraphFree(
  TU* tu,           /**< \ref TU environment. */
  TU_GRAPH** pgraph /**< Pointer to graph. */
);

/**
 * \brief Removes all nodes and columns, keeping the memory.
 */

TU_ERROR TUgraphClear(
  TU* tu,         /**< \ref TU environment. */
  TU_GRAPH* graph /**< Graph structure. */
);

/**
 * \brief Adds a node to a graph.
 *
 * Adds a node to a graph.
 *
 * \return Node structure of new node.
 */

TU_ERROR TUgraphAddNode(
  TU* tu,               /**< \ref TU environment. */
  TU_GRAPH* graph,      /**< Graph. */
  TU_GRAPH_NODE* pnode  /**< Pointer for storing the new node, or \c NULL. */
);

/**
 * \brief Adds an edge to a graph.
 *
 * Adds an edge to a graph.
 *
 * \return Edge of new edge.
 */

TU_ERROR TUgraphAddEdge(
  TU* tu,               /**< \ref TU environment. */
  TU_GRAPH* graph,      /**< Graph. */
  TU_GRAPH_NODE u,      /**< One node of the edge. */
  TU_GRAPH_NODE v,      /**< Other node of the edge. */
  TU_GRAPH_EDGE* pedge  /**< Pointer for storinge the new edge, or \c NULL.*/
);

/**
 * \brief Removes node \p v and all its incident edges from \p graph.
 */

TU_ERROR TUgraphDeleteNode(
  TU* tu,           /**< \ref TU environment. */
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_NODE v   /**< Node to be deleted. */
);

/**
 * \brief Removes edge \p e from \p graph.
 */

TU_ERROR TUgraphDeleteEdge(
  TU* tu,           /**< \ref TU environment. */
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_EDGE e   /**< Edge to be deleted. */
);

/**
 * \brief Returns a node iterator for iterating over all nodes.
 */

static inline
TU_GRAPH_NODE TUgraphNodesFirst(
  TU_GRAPH* graph /**< Graph. */
)
{
  return graph->firstNode;
}

/**
 * \brief Return \c true if and only if this \p v is not the last node in node iteration.
 */

static inline
bool TUgraphNodesValid(
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_NODE v   /**< Node. */
)
{
  return v >= 0;
}

/**
 * \brief Returns the next node after \p v.
 */

static inline
TU_GRAPH_NODE TUgraphNodesNext(
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_NODE v   /**< Node. */
)
{
  return graph->nodes[v].next;
}

/**
 * \brief Returns an iterator for all edges incident to node \p v.
 */

static inline
TU_GRAPH_ITER TUgraphIncFirst(
  TU_GRAPH* graph,  /**< Graph. */
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

/**
 * \brief Returns \c true if iterator \p i for all incident edges of some node is valid.
 */

static inline
bool TUgraphIncValid(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return i >= 0;
}

/**
 * \brief Returns the iterator following iterator \p i for all edges incident to some node.
 */

static inline
TU_GRAPH_ITER TUgraphIncNext(
  TU_GRAPH* graph,  /**< Graph. */
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

/**
 * \brief Converts an iterator for edges incident to a node to the actual edge.
 */

static inline
TU_GRAPH_EDGE TUgraphIncEdge(
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return i/2;
}

/**
 * \brief Returns the node of which iterator \p i traverses through incident edges.
 */

static inline
TU_GRAPH_NODE TUgraphIncSource(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return graph->arcs[i^1].target;
}

/**
 * \brief Returns the end node of the edge corresponding to this iterator \p i.
 */

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
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_ITER i   /**< Current edge iterator. */
)
{
  TU_GRAPH_NODE source = graph->arcs[i ^ 1].target;
  TU_GRAPH_ITER j = graph->arcs[i].next;
  while (true)
  {
    while (j >= 0)
    {
      if (!(j & 0x1))
        return j;
      j = graph->arcs[j].next;
    }
    source = graph->nodes[source].next;
    if (source < 0)
      return -1;

    j = graph->nodes[source].firstOut;
  }
}

/**
 * \brief Returns iterator for all edges of \p graph.
 */

static inline
TU_GRAPH_ITER TUgraphEdgesFirst(
  TU_GRAPH* graph /**< Graph. */
)
{
  for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
    v = TUgraphNodesNext(graph, v))
  {
    TU_GRAPH_ITER i = graph->nodes[v].firstOut;
    if (i >= 0)
    {
      /* We found some node with an incident edge. */
      if (i & 0x1)
        i = TUgraphEdgesNext(graph, i);
      return i;
    }
  }
  return -1;
}

/**
 * \brief Returns \c true if and only if iterator \p i for all edges is valid.
 */

static inline
bool TUgraphEdgesValid(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return i >= 0;
}

/**
 * \brief Returns the actual edge, iterator \p i represents.
 */

static inline
TU_GRAPH_EDGE TUgraphEdgesEdge(
  TU_GRAPH* graph,  /**< Graph structure. */
  TU_GRAPH_ITER i   /**< Iterator for edges. */
)
{
  return i/2;
}

/**
 * \brief Prints the \p graph, writing to \p stream.
 */

TU_ERROR TUgraphPrint(
  FILE* stream,   /**< Stream. */
  TU_GRAPH* graph /**< Graph structure. */
);

/**
 * \brief Merges two nodes \p u and \p v of \p graph.
 */

TU_ERROR TUgraphMergeNodes(
  TU* tu,           /**< \ref TU environment. */
  TU_GRAPH* graph,  /**< Graph. */
  TU_GRAPH_NODE u,  /**< First node. */
  TU_GRAPH_NODE v   /**< Second node. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPH_H */
