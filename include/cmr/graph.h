#ifndef CMR_GRAPH_H
#define CMR_GRAPH_H

/**
 * \file graph.h
 *
 * \author Matthias Walter
 *
 * \brief Functionality for graphs.
 */


#include <cmr/env.h>
#include <cmr/element.h>

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

typedef int CMR_GRAPH_NODE; /**< \brief Reference to a node of \ref CMR_GRAPH. */
typedef int CMR_GRAPH_EDGE; /**< \brief Reference to an edge of \ref CMR_GRAPH. */
typedef int CMR_GRAPH_ITER; /**< \brief Reference to an edge iterator of \ref CMR_GRAPH. */

typedef struct
{
  int prev;     /**< \brief Next node in node list. */
  int next;     /**< \brief Previous node in node list. */
  int firstOut; /**< \brief First out-arc of this node. */
} CMR_GRAPH_NODE_DATA;

typedef struct
{
  int target; /**< \brief Target node of this arc. */
  int prev;   /**< \brief Next arc in out-arc list of source node. */
  int next;   /**< \brief Previous arc in out-arc list of source node. */
} CMR_GRAPH_ARC_DATA;

typedef struct
{
  size_t numNodes;            /**< \brief Number of nodes. */
  size_t memNodes;            /**< \brief Number of nodes for which memory is allocated. */
  CMR_GRAPH_NODE_DATA* nodes; /**< \brief Array containing node data. */
  int firstNode;              /**< \brief Index of first node. */
  int freeNode;               /**< \brief Beginning of free-list of nodes. */

  size_t numEdges;            /**< \brief Number of edges. */
  size_t memEdges;            /**< \brief Number of edges for which memory is allocated. */
  CMR_GRAPH_ARC_DATA* arcs;   /**< \brief Array containing arc data. */
  int freeEdge;               /**< \brief Beginning of free-list of arc. */
} CMR_GRAPH;

/**
 * \brief Returns number of nodes for which memory is allocated.
 */

static inline
size_t CMRgraphMemNodes(
  CMR_GRAPH* graph /**< Graph. */
)
{
  assert(graph);

  return graph->memNodes;
}

/**
 * \brief Returns number of nodes.
 */

static inline
size_t CMRgraphNumNodes(
  CMR_GRAPH* graph /**< Graph. */
)
{
  assert(graph);

  return graph->numNodes;
}

/**
 * \brief Returns number of edges for which memory is allocated.
 */

static inline
size_t CMRgraphMemEdges(
  CMR_GRAPH* graph /**< Graph. */
)
{
  assert(graph);

  return graph->memEdges;
}

/**
 * \brief Returns number of edges.
 */

static inline
size_t CMRgraphNumEdges(
  CMR_GRAPH* graph /**< Graph. */
)
{
  assert(graph);

  return graph->numEdges;
}

/**
 * \brief Returns node u of edge \p e.
 */

static inline
CMR_GRAPH_NODE CMRgraphEdgeU(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_EDGE e   /**< Edge of \p graph. */
)
{
  assert(graph);

  return graph->arcs[2*e+1].target;
}

/**
 * \brief Returns node v of edge \p e.
 */

static inline
CMR_GRAPH_NODE CMRgraphEdgeV(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_EDGE e   /**< Edge of \p graph. */
)
{
  assert(graph);

  return graph->arcs[2*e].target;
}

/**
 * \brief Creates an empty graph.
 */

CMR_EXPORT
CMR_ERROR CMRgraphCreateEmpty(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_GRAPH** pgraph,  /**< Pointer for storing the graph. */
  int memNodes,       /**< Allocate memory for this number of nodes. */
  int memEdges        /**< Allocate memory for this number of edges. */
);

/**
 * \brief Frees a graph.
 */

CMR_EXPORT
CMR_ERROR CMRgraphFree(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_GRAPH** pgraph /**< Pointer to graph. */
);

/**
 * \brief Removes all nodes and columns, keeping the memory.
 */

CMR_EXPORT
CMR_ERROR CMRgraphClear(
  CMR* cmr,         /**< \ref CMR environment. */
  CMR_GRAPH* graph /**< Graph structure. */
);

/**
 * \brief Adds a node to a graph.
 *
 * Adds a node to a graph.
 *
 * \return Node structure of new node.
 */

CMR_EXPORT
CMR_ERROR CMRgraphAddNode(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_GRAPH* graph,      /**< Graph. */
  CMR_GRAPH_NODE* pnode  /**< Pointer for storing the new node, or \c NULL. */
);

/**
 * \brief Adds an edge to a graph.
 *
 * Adds an edge to a graph.
 *
 * \return Edge of new edge.
 */

CMR_EXPORT
CMR_ERROR CMRgraphAddEdge(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_GRAPH* graph,      /**< Graph. */
  CMR_GRAPH_NODE u,      /**< One node of the edge. */
  CMR_GRAPH_NODE v,      /**< Other node of the edge. */
  CMR_GRAPH_EDGE* pedge  /**< Pointer for storinge the new edge, or \c NULL.*/
);

/**
 * \brief Removes node \p v and all its incident edges from \p graph.
 */

CMR_EXPORT
CMR_ERROR CMRgraphDeleteNode(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_NODE v   /**< Node to be deleted. */
);

/**
 * \brief Removes edge \p e from \p graph.
 */

CMR_EXPORT
CMR_ERROR CMRgraphDeleteEdge(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_EDGE e   /**< Edge to be deleted. */
);

/**
 * \brief Returns a node iterator for iterating over all nodes.
 */

static inline
CMR_GRAPH_NODE CMRgraphNodesFirst(
  CMR_GRAPH* graph /**< Graph. */
)
{
  assert(graph);

  return graph->firstNode;
}

/**
 * \brief Return \c true if and only if this \p v is not the last node in node iteration.
 */

static inline
bool CMRgraphNodesValid(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_NODE v   /**< Node. */
)
{
  CMR_UNUSED(graph);

  assert(graph);

  return v >= 0;
}

/**
 * \brief Returns the next node after \p v.
 */

static inline
CMR_GRAPH_NODE CMRgraphNodesNext(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_NODE v   /**< Node. */
)
{
  assert(graph);

  return graph->nodes[v].next;
}

/**
 * \brief Returns an iterator for all edges incident to node \p v.
 */

static inline
CMR_GRAPH_ITER CMRgraphIncFirst(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_NODE v   /**< Node. */
)
{
  assert(graph);

  CMR_GRAPH_ITER i = graph->nodes[v].firstOut;
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
bool CMRgraphIncValid(
  CMR_GRAPH* graph,  /**< Graph structure. */
  CMR_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  CMR_UNUSED(graph);

  assert(graph);

  return i >= 0;
}

/**
 * \brief Returns the iterator following iterator \p i for all edges incident to some node.
 */

static inline
CMR_GRAPH_ITER CMRgraphIncNext(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  assert(graph);

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
CMR_GRAPH_EDGE CMRgraphIncEdge(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  CMR_UNUSED(graph);

  assert(graph);

  return i/2;
}

/**
 * \brief Returns the node of which iterator \p i traverses through incident edges.
 */

static inline
CMR_GRAPH_NODE CMRgraphIncSource(
  CMR_GRAPH* graph,  /**< Graph structure. */
  CMR_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return graph->arcs[i^1].target;
}

/**
 * \brief Returns the end node of the edge corresponding to this iterator \p i.
 */

static inline
CMR_GRAPH_NODE CMRgraphIncTarget(
  CMR_GRAPH* graph,  /**< Graph structure. */
  CMR_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  return graph->arcs[i].target;
}

/**
 * \brief Returns iterator of next edge in list of all edges.
 */

static inline
CMR_GRAPH_ITER CMRgraphEdgesNext(
  CMR_GRAPH* graph,  /**< Graph. */
  CMR_GRAPH_ITER i   /**< Current edge iterator. */
)
{
  CMR_GRAPH_NODE source = graph->arcs[i ^ 1].target;
  CMR_GRAPH_ITER j = graph->arcs[i].next;
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
CMR_GRAPH_ITER CMRgraphEdgesFirst(
  CMR_GRAPH* graph /**< Graph. */
)
{
  assert(graph);

  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v);
    v = CMRgraphNodesNext(graph, v))
  {
    CMR_GRAPH_ITER i = graph->nodes[v].firstOut;
    if (i >= 0)
    {
      /* We found some node with an incident edge. */
      if (i & 0x1)
        i = CMRgraphEdgesNext(graph, i);
      return i;
    }
  }
  return -1;
}

/**
 * \brief Returns \c true if and only if iterator \p i for all edges is valid.
 */

static inline
bool CMRgraphEdgesValid(
  CMR_GRAPH* graph,  /**< Graph structure. */
  CMR_GRAPH_ITER i   /**< Iterator for edges incident to a node. */
)
{
  CMR_UNUSED(graph);

  assert(graph);

  return i >= 0;
}

/**
 * \brief Returns the actual edge, iterator \p i represents.
 */

static inline
CMR_GRAPH_EDGE CMRgraphEdgesEdge(
  CMR_GRAPH* graph,  /**< Graph structure. */
  CMR_GRAPH_ITER i   /**< Iterator for edges. */
)
{
  CMR_UNUSED(graph);

  assert(graph);

  return i/2;
}

/**
 * \brief Prints the \p graph, writing to \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRgraphPrint(
  CMR_GRAPH* graph, /**< Graph structure. */
  FILE* stream      /**< Stream. */
);

/**
 * \brief Merges two nodes \p u and \p v of \p graph.
 */

CMR_EXPORT
CMR_ERROR CMRgraphMergeNodes(
  CMR* cmr,         /**< \ref CMR environment. */
  CMR_GRAPH* graph, /**< Graph. */
  CMR_GRAPH_NODE u, /**< First node. */
  CMR_GRAPH_NODE v  /**< Second node. */
);


CMR_EXPORT
CMR_ERROR CMRgraphCreateFromEdgeList(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH** pgraph,           /**< Pointer for storing the graph. */
  CMR_ELEMENT** pedgeElements,  /**< Pointer for storing element of each edge (may be \c NULL). */
  char*** pnodeLabels,          /**< Pointer for storing string node labels (may be \c NULL). */
  FILE* stream                  /**< File stream to read from. */
);

CMR_EXPORT
CMR_ERROR CMRgraphCopy(
  CMR* cmr,         /**< \ref CMR environment. */
  CMR_GRAPH* graph, /**< Graph structure. */
  CMR_GRAPH** pcopy /**< Pointer for storing the copied graph. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_GRAPH_H */
