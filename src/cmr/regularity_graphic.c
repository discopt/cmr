// #define CMR_DEBUG /* Uncomment to debug this file. */

#include <cmr/graphic.h>
#include <cmr/network.h>

#include "seymour_internal.h"
#include "env_internal.h"
#include "hashtable.h"

#include <stdint.h>
#include <time.h>

/**
 * \brief Recursive DFS for finding all articulation points of a graph.
 */

static
int dfsArticulationPoint(
  CMR_GRAPH* graph,               /**< Graph. */
  bool* edgesEnabled,             /**< Edge array indicating whether an edge is enabled. */
  CMR_GRAPH_NODE node,            /**< Current node. */
  bool* nodesVisited,             /**< Node array indicating whether a node was already visited. */
  int* nodesDiscoveryTime,        /**< Node array indicating at which time a node was visited. */
  int* ptime,                     /**< Pointer to current time. */
  CMR_GRAPH_NODE parentNode,      /**< Parent node in DFS arborescence. */
  size_t* nodesArticulationPoint  /**< Node array indicating whether a node is an articulation point. */
)
{
  assert(graph);
  assert(nodesVisited);
  assert(nodesDiscoveryTime);
  assert(ptime);
  assert(nodesArticulationPoint);

  size_t numChildren = 0;
  nodesVisited[node] = true;
  ++(*ptime);
  nodesDiscoveryTime[node] = *ptime;
  int earliestReachableTime = *ptime;

  for (CMR_GRAPH_ITER iter = CMRgraphIncFirst(graph, node); CMRgraphIncValid(graph, iter);
    iter = CMRgraphIncNext(graph, iter))
  {
    assert(CMRgraphIncSource(graph, iter) == node);
    if (!edgesEnabled[CMRgraphIncEdge(graph, iter)])
      continue;

    CMR_GRAPH_NODE v = CMRgraphIncTarget(graph, iter);
    if (!nodesVisited[v])
    {
      ++numChildren;
      int childEarliestReachableTime = dfsArticulationPoint(graph, edgesEnabled, v, nodesVisited, nodesDiscoveryTime,
        ptime, node, nodesArticulationPoint);
      if (childEarliestReachableTime < earliestReachableTime)
        earliestReachableTime = childEarliestReachableTime;
      if (parentNode >= 0 && childEarliestReachableTime >= nodesDiscoveryTime[node])
        nodesArticulationPoint[node] = 1;
    }
    else if (v != parentNode && nodesDiscoveryTime[v] < earliestReachableTime)
      earliestReachableTime = nodesDiscoveryTime[v];
  }

  if (parentNode < 0)
  {
    if (numChildren > 1)
      nodesArticulationPoint[node] = 1;
  } 

  return earliestReachableTime;
}

static
CMR_ERROR findArticulationPoints(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_GRAPH* graph,               /**< Graph. */
  CMR_GRAPH_EDGE* columnEdges,    /**< Array with with map from columns to edges. */
  size_t* nodesArticulationPoint, /**< Node array indicating whether a node is an articulation point. */
  size_t* nonzeroColumns,         /**< Array with columns containing a nonzero. */
  size_t numNonzeroColumns        /**< Length of \p nonzeroColumns. */
)
{
  assert(graph);
  assert(columnEdges);
  assert(nodesArticulationPoint);
  assert(nonzeroColumns);

  bool* nodesVisited = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodesVisited, CMRgraphMemNodes(graph)) );
  int* nodesDiscoveryTime = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodesDiscoveryTime, CMRgraphMemNodes(graph)) );
  bool* edgesEnabled = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &edgesEnabled, CMRgraphMemEdges(graph)) );

  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
  {
    nodesArticulationPoint[v] = 0;
    nodesVisited[v] = false;
    nodesDiscoveryTime[v] = 0;
  }
  for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
    iter = CMRgraphEdgesNext(graph, iter))
  {
    edgesEnabled[CMRgraphEdgesEdge(graph, iter)] = true;
  }

  for (size_t i = 0; i < numNonzeroColumns; ++i)
    edgesEnabled[columnEdges[nonzeroColumns[i]]] = false;  

  int time = 0;
  dfsArticulationPoint(graph, edgesEnabled, CMRgraphNodesFirst(graph), nodesVisited, nodesDiscoveryTime, &time, -1,
    nodesArticulationPoint);

  CMR_CALL( CMRfreeStackArray(cmr, &edgesEnabled) );
  CMR_CALL( CMRfreeStackArray(cmr, &nodesDiscoveryTime) );
  CMR_CALL( CMRfreeStackArray(cmr, &nodesVisited) );

  return CMR_OKAY;
}

static
void dfsTree(
  CMR_GRAPH* graph,             /**< Graph. */
  bool* edgesTree,              /**< Edge array indicating whether an edge is enabled. */
  bool* nodesVisited,           /**< Node array indicating whether a node was already visited. */
  CMR_GRAPH_NODE* nodesParent,  /**< Node array indicating the parent node of each node. */
  CMR_GRAPH_NODE node           /**< Current node. */
)
{
  assert(graph);
  assert(edgesTree);
  assert(nodesVisited);
  assert(nodesParent);
  assert(node >= 0);

  nodesVisited[node] = true;
  for (CMR_GRAPH_ITER iter = CMRgraphIncFirst(graph, node); CMRgraphIncValid(graph, iter);
    iter = CMRgraphIncNext(graph, iter))
  {
    assert(CMRgraphIncSource(graph, iter) == node);
    if (edgesTree[CMRgraphIncEdge(graph, iter)])
    {
      CMR_GRAPH_NODE v = CMRgraphIncTarget(graph, iter);
      if (!nodesVisited[v])
      {
        nodesParent[v] = node;
        dfsTree(graph, edgesTree, nodesVisited, nodesParent, v);
      }
    }
  }
}

static
CMR_ERROR findTreeParents(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_GRAPH* graph,           /**< Graph. */
  CMR_GRAPH_EDGE* rowEdges,   /**< Array with with map from rows to edges. */
  size_t numRows,             /**< Number of rows. */
  CMR_GRAPH_NODE* nodesParent /**< Array to be filled with map from nodes to parent nodes. */
)
{
  assert(graph);
  assert(rowEdges);
  assert(nodesParent);

  bool* nodesVisited = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodesVisited, CMRgraphMemNodes(graph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
    nodesVisited[v] = false;

  bool* edgesTree = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &edgesTree, CMRgraphMemEdges(graph)) );
  for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
    iter = CMRgraphEdgesNext(graph, iter))
  {
    edgesTree[ CMRgraphEdgesEdge(graph, iter) ] = false;
  }

  /* Enable tree edge. */
  for (size_t row = 0; row < numRows; ++row)
    edgesTree[ rowEdges[row] ] = true;

  CMR_GRAPH_NODE root = CMRgraphNodesFirst(graph);
  nodesParent[root] = -1;
  dfsTree(graph, edgesTree, nodesVisited, nodesParent, root);

  CMR_CALL( CMRfreeStackArray(cmr, &edgesTree) );
  CMR_CALL( CMRfreeStackArray(cmr, &nodesVisited) );

  return CMR_OKAY;
}

static
void dfsComponents(
  CMR_GRAPH* graph,
  bool* edgesEnabled,
  size_t* nodesComponent,
  CMR_GRAPH_NODE node,
  size_t component
)
{
  assert(graph);
  assert(edgesEnabled);
  assert(nodesComponent);
  assert(node >= 0);
  assert(component < SIZE_MAX);
  
  nodesComponent[node] = component;
  for (CMR_GRAPH_ITER iter = CMRgraphIncFirst(graph, node); CMRgraphIncValid(graph, iter);
    iter = CMRgraphIncNext(graph, iter))
  {
    assert(CMRgraphIncSource(graph, iter) == node);
    if (edgesEnabled[CMRgraphIncEdge(graph, iter)])
    {
      CMR_GRAPH_NODE v = CMRgraphIncTarget(graph, iter);
      if (nodesComponent[v] == SIZE_MAX)
        dfsComponents(graph, edgesEnabled, nodesComponent, v, component);
    }
  }
}

static
CMR_ERROR findComponents(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* graph,             /**< Graph. */
  CMR_GRAPH_EDGE* columnEdges,  /**< Array with with map from columns to edges. */
  CMR_GRAPH_NODE removedNode,   /**< Node that shall be considered as removed. */
  size_t* nodesComponent,       /**< Node array indicating the components. */
  size_t* pnumComponents,       /**< Pointer for storing the number of connected components. */
  size_t* nonzeroColumns,       /**< Array with columns containing a nonzero. */
  size_t numNonzeroColumns      /**< Length of \p nonzeroColumns. */
)
{
  assert(graph);
  assert(columnEdges);
  assert(nonzeroColumns);

  bool* edgesEnabled = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &edgesEnabled, CMRgraphMemEdges(graph)) );
  for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
    iter = CMRgraphEdgesNext(graph, iter))
  {
    edgesEnabled[ CMRgraphEdgesEdge(graph, iter) ] = true;
  }

  /* Disable edges around special node. */
  for (CMR_GRAPH_ITER iter = CMRgraphIncFirst(graph, removedNode); CMRgraphIncValid(graph, iter);
    iter = CMRgraphIncNext(graph, iter))
  {
    edgesEnabled[ CMRgraphIncEdge(graph, iter) ] = false;
  }

  /* Disable 1-edges. */
  for (size_t i = 0; i < numNonzeroColumns; ++i)
    edgesEnabled[ columnEdges[ nonzeroColumns[i] ] ] = false;

  /* Initialize components. */
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
    nodesComponent[v] = SIZE_MAX;

  size_t component = 0;
  for (CMR_GRAPH_NODE source = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, source);
    source = CMRgraphNodesNext(graph, source))
  {
    if (nodesComponent[source] == SIZE_MAX && source != removedNode)
    {
      dfsComponents(graph, edgesEnabled, nodesComponent, source, component);
      ++component;
    }
  }
  *pnumComponents = component;

  CMR_CALL( CMRfreeStackArray(cmr, &edgesEnabled) );

  return CMR_OKAY;
}

/**
 * \brief DFS for searching for a bipartition.
 */

static
bool dfsBipartite(
  CMR_GRAPH* graph,   /**< Graph. */
  bool* nodesVisited, /**< Node array indicating whether a node was visited already. */
  int* bipartition,   /**< Node array for storing the bipartition. */
  CMR_GRAPH_NODE node /**< Current node. */
)
{
  nodesVisited[node] = true;
  for (CMR_GRAPH_ITER iter = CMRgraphIncFirst(graph, node); CMRgraphIncValid(graph, iter);
    iter = CMRgraphIncNext(graph, iter))
  {
    assert(CMRgraphIncSource(graph, iter) == node);
    CMR_GRAPH_NODE v = CMRgraphIncTarget(graph, iter);
    if (nodesVisited[v])
    {
      if (bipartition[v] == bipartition[node])
        return false;
    }
    else
    {
      bipartition[v] = 1 - bipartition[node];
      bool isBipartite = dfsBipartite(graph, nodesVisited, bipartition, v);
      if (!isBipartite)
        return false;
    }
  }
  return true;
}

/**
 * \brief Finds a bipartition of a graph.
 */

static
CMR_ERROR findBipartition(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_GRAPH* graph,   /**< Graph. */
  int* bipartition,   /**< Node array indicating color class. */
  bool* pisBipartite  /**< Pointer for storing whether \p is bipartite. */
)
{
  assert(cmr);
  assert(graph);
  assert(bipartition);

  bool* nodesVisited = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodesVisited, CMRgraphMemNodes(graph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
    nodesVisited[v] = false;

  bool isBipartite = true;
  for (CMR_GRAPH_NODE source = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, source );
    source = CMRgraphNodesNext(graph, source))
  {
    if (!nodesVisited[source])
    {
      bipartition[source] = 0;
      if (!dfsBipartite(graph, nodesVisited, bipartition, source))
      {
        isBipartite = false;
        break;
      }
    }
  }

  if (pisBipartite)
    *pisBipartite = isBipartite;

  CMR_CALL( CMRfreeStackArray(cmr, &nodesVisited) );

  return CMR_OKAY;
}


/**
 * \brief Extends \p graph for a submatrix augmented by 1 row.
 */

static
CMR_ERROR addToGraph1Row(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* graph,             /**< Empty graph to be filled. */
  CMR_GRAPH_EDGE* rowEdges,     /**< Array to be filled with map from rows to edges. */
  CMR_GRAPH_EDGE* columnEdges,  /**< Array to be filled with map from columns to edges. */
  size_t baseNumRows,           /**< Number of rows already processed. */
  size_t* nonzeroColumns,       /**< Array with columns containing a nonzero. */
  size_t numNonzeroColumns,     /**< Length of \p nonzeroColumns. */
  bool* pisGraphic              /**< Pointer for storing whether this extension was graphic. */
)
{
  assert(cmr);
  assert(graph);
  assert(baseNumRows >= 3);
  assert(rowEdges);
  assert(columnEdges);

  
  /* We first check whether the edges of columns with a nonzero form a star. */
  CMR_GRAPH_NODE starNode = -1;
  CMR_GRAPH_NODE u1 = -1;
  CMR_GRAPH_NODE v1 = -1;
  for (size_t i = 0; i < numNonzeroColumns; ++i)
  {
    size_t column = nonzeroColumns[i];
    CMR_GRAPH_EDGE e = columnEdges[column];
    if (i == 0)
    {
      u1 = CMRgraphEdgeU(graph, e);
      v1 = CMRgraphEdgeV(graph, e);
    }
    else if (i == 1)
    {
      starNode = CMRgraphEdgeU(graph, e);
      if (starNode == u1 || starNode == v1)
        continue;
      starNode = CMRgraphEdgeV(graph, e);
      if (starNode != u1 && starNode != v1)
      {
        starNode = -1;
        break;
      }
    }
    else
    {
      u1 = CMRgraphEdgeU(graph, e);
      v1 = CMRgraphEdgeV(graph, e);
      if (u1 != starNode && v1 != starNode)
      {
        starNode = -1;
        break;
      }
    }
  }

  /* The edges form a star. */
  if (starNode >= 0)
  {
    *pisGraphic = true;
    CMR_GRAPH_NODE newNode;
    CMR_CALL( CMRgraphAddNode(cmr, graph, &newNode) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, starNode, newNode, &rowEdges[baseNumRows]) );
    for (size_t i = 0; i < numNonzeroColumns; ++i)
    {
      CMR_GRAPH_EDGE edge = columnEdges[nonzeroColumns[i]];
      CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, edge);
      CMR_GRAPH_EDGE newEdge;
      if (u == starNode)
        u = CMRgraphEdgeV(graph, edge);
      CMR_CALL( CMRgraphDeleteEdge(cmr, graph, edge) );
      CMR_CALL( CMRgraphAddEdge(cmr, graph, u, newNode, &newEdge) );
      assert(newEdge == edge);
    }

    return CMR_OKAY;
  }

  size_t* nodesCandidate = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodesCandidate, CMRgraphMemNodes(graph)) );

  CMR_CALL( findArticulationPoints(cmr, graph, columnEdges, nodesCandidate, nonzeroColumns, numNonzeroColumns) );
  
  size_t countCandidates = 0;
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
  {
    if (nodesCandidate[v])
      ++countCandidates;
  }

  CMRdbgMsg(12, "Found %ld articulation points.\n", countCandidates);

  *pisGraphic = false;
  if (countCandidates > 0)
  {
    /* We need a rooted arborescence along the row (tree) edges. */
    CMR_GRAPH_NODE* nodesParent = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &nodesParent, CMRgraphMemNodes(graph)) );
    CMR_CALL( findTreeParents(cmr, graph, rowEdges, baseNumRows, nodesParent) );

    /* Ensure that the fundamental cycles induced by the column-edges with a 1-entry go through the articular points. */
    CMR_GRAPH_NODE* nodeStacks[2] = { NULL, NULL };
    CMR_CALL( CMRallocStackArray(cmr, &nodeStacks[0], baseNumRows+1) );
    CMR_CALL( CMRallocStackArray(cmr, &nodeStacks[1], baseNumRows+1) );
    size_t nodeStackSizes[2];
    for (size_t i = 0; i < numNonzeroColumns; ++i)
    {
      CMR_GRAPH_EDGE columnEdge = columnEdges[nonzeroColumns[i]];
      CMR_GRAPH_NODE nodes[2] = { CMRgraphEdgeU(graph, columnEdge), CMRgraphEdgeV(graph, columnEdge) };
      for (int j = 0; j < 2; ++j)
      {
        nodeStackSizes[j] = 0;
        for (CMR_GRAPH_NODE v = nodes[j]; v >= 0; v = nodesParent[v])
          nodeStacks[j][nodeStackSizes[j]++] = v;
      }
      CMRdbgMsg(12, "For nonzero c%ld, paths to root %ld have lengths %ld and %ld.\n", nonzeroColumns[i]+1,
        nodeStacks[0][nodeStackSizes[0]-1], nodeStackSizes[0], nodeStackSizes[1]);

      while (nodeStackSizes[0] > 0 && nodeStackSizes[1] > 0
        && nodeStacks[0][nodeStackSizes[0]-1] == nodeStacks[1][nodeStackSizes[1]-1])
      {
        nodeStackSizes[0]--;
        nodeStackSizes[1]--;
      }
      nodeStackSizes[0]++;

      CMRdbgMsg(12, "For nonzero c%ld, pruned paths have lengths %ld and %ld.\n", nonzeroColumns[i]+1,
        nodeStackSizes[0], nodeStackSizes[1]);

      countCandidates = 0;
      for (int j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < nodeStackSizes[j]; ++k)
        {
          CMR_GRAPH_NODE v = nodeStacks[j][k];
          if (nodesCandidate[v] == i+1)
          {
            nodesCandidate[v]++;
            countCandidates++;
          }
        }
      }

      CMRdbgMsg(12, "Number of candidate points is %ld.\n", countCandidates);
      if (countCandidates == 0)
        break;
    }

    CMR_CALL( CMRfreeStackArray(cmr, &nodeStacks[1]) );
    CMR_CALL( CMRfreeStackArray(cmr, &nodeStacks[0]) );
    CMR_CALL( CMRfreeStackArray(cmr, &nodesParent) );

    if (countCandidates == 1 || countCandidates == 2)
    {
      for (CMR_GRAPH_NODE splitNode = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, splitNode);
        splitNode = CMRgraphNodesNext(graph, splitNode))
      {
        if (nodesCandidate[splitNode] != numNonzeroColumns+1)
          continue;
      
        CMRdbgMsg(12, "Candidate node is %ld.\n", splitNode);

#if defined(CMR_DEBUG)
        CMRgraphPrint(graph, stdout);
        fflush(stdout);
#endif /* CMR_DEBUG */

        for (size_t i = 0; i < numNonzeroColumns; ++i)
        {
          CMRdbgMsg(14, "1-edge {%ld,%ld}\n", CMRgraphEdgeU(graph, columnEdges[nonzeroColumns[i]]),
            CMRgraphEdgeV(graph, columnEdges[nonzeroColumns[i]]));
        }

        size_t numComponents = SIZE_MAX;
        size_t* nodesComponent = NULL;
        CMR_CALL( CMRallocStackArray(cmr, &nodesComponent, CMRgraphNumNodes(graph)) );

        CMR_CALL( findComponents(cmr, graph, columnEdges, splitNode, nodesComponent, &numComponents, nonzeroColumns,
          numNonzeroColumns) );
        for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
        {
          CMRdbgMsg(14, "Node %ld belongs to component %ld of %ld.\n", v, nodesComponent[v], numComponents);
        }
        assert(numComponents >= 2);

        CMR_GRAPH* auxiliaryGraph = NULL;
        CMR_CALL( CMRgraphCreateEmpty(cmr, &auxiliaryGraph, numComponents, numNonzeroColumns) );
        CMR_GRAPH_NODE* componentAuxiliaryNodes = NULL;
        CMR_CALL( CMRallocStackArray(cmr, &componentAuxiliaryNodes, numComponents) );
        for (size_t c = 0; c < numComponents; ++c)
          CMR_CALL( CMRgraphAddNode(cmr, auxiliaryGraph, &componentAuxiliaryNodes[c]) );

        for (size_t i = 0; i < numNonzeroColumns; ++i)
        {
          CMR_GRAPH_EDGE e = columnEdges[nonzeroColumns[i]];
          size_t components[2] = { nodesComponent[CMRgraphEdgeU(graph, e)], nodesComponent[CMRgraphEdgeV(graph, e)] };
          if (components[0] < SIZE_MAX && components[1] < SIZE_MAX)
          {
            CMR_CALL( CMRgraphAddEdge(cmr, auxiliaryGraph, componentAuxiliaryNodes[components[0]],
              componentAuxiliaryNodes[components[1]], NULL) );
          }
        }

  #if defined(CMR_DEBUG)
        CMRdbgMsg(14, "Constructed auxiliary graph.\n");
        CMRgraphPrint(auxiliaryGraph, stdout);
        fflush(stdout);
  #endif /* CMR_DEBUG */

        bool isBipartite = false;
        int* bipartition = NULL;
        CMR_CALL( CMRallocStackArray(cmr, &bipartition, CMRgraphMemNodes(auxiliaryGraph)) );

        CMR_CALL( findBipartition(cmr, auxiliaryGraph, bipartition, &isBipartite) );
        if (isBipartite)
        {
          *pisGraphic = true;

          for (size_t c = 0; c < numComponents; ++c)
          {
            CMRdbgMsg(16, "Component %ld belongs to bipartition %d.\n", c, bipartition[componentAuxiliaryNodes[c]]);
          }

          /* we carry out the re-assignment. */

          CMR_GRAPH_NODE sisterNode;
          CMR_CALL( CMRgraphAddNode(cmr, graph, &sisterNode) );

          /* We mark the 1-edges. */
          bool* edges1 = NULL;
          CMR_CALL( CMRallocStackArray(cmr, &edges1, CMRgraphMemEdges(graph)) );
          for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
            iter = CMRgraphEdgesNext(graph, iter))
          {
            edges1[CMRgraphEdgesEdge(graph, iter)] = false;
          }
          for (size_t i = 0; i < numNonzeroColumns; ++i)
            edges1[columnEdges[nonzeroColumns[i]]] = true;

          /* We store the incident edges since we change that list. */
          size_t numIncidentEdges = 0;
          CMR_GRAPH_EDGE* incidentEdges = NULL;
          CMR_CALL( CMRallocStackArray(cmr, &incidentEdges, CMRgraphNumNodes(graph)) );
          for (CMR_GRAPH_ITER iter = CMRgraphIncFirst(graph, splitNode); CMRgraphIncValid(graph, iter);
            iter = CMRgraphIncNext(graph, iter))
          {
            incidentEdges[numIncidentEdges++] = CMRgraphIncEdge(graph, iter);
          }

          for (size_t i = 0; i < numIncidentEdges; ++i)
          {
            CMR_GRAPH_EDGE edge = incidentEdges[i];
            CMR_GRAPH_NODE v = CMRgraphEdgeU(graph, edge);
            if (v == splitNode)
              v = CMRgraphEdgeV(graph, edge);
            int side = bipartition[componentAuxiliaryNodes[nodesComponent[v]]];
            CMRdbgMsg(16, "Node %ld of edge {%ld,%ld} belongs to bipartition side %d.\n", v, CMRgraphEdgeU(graph, edge),
              CMRgraphEdgeV(graph, edge), side);

            /* Complement decision for 1-edges. */
            if (edges1[edge])
              side = 1-side;

            if (side)
            {
              /* Reconnect the edge to the sister node. */
              CMR_CALL( CMRgraphDeleteEdge(cmr, graph, edge) );
              CMR_GRAPH_EDGE modifiedEdge;
              CMR_CALL( CMRgraphAddEdge(cmr, graph, v, sisterNode, &modifiedEdge) );
              assert(modifiedEdge == edge);
            }
          }
          
          /* Finally, connect the split node and the sister node. */
          CMR_CALL( CMRgraphAddEdge(cmr, graph, splitNode, sisterNode, &rowEdges[baseNumRows]) );
    
          CMR_CALL( CMRfreeStackArray(cmr, &incidentEdges) );
          CMR_CALL( CMRfreeStackArray(cmr, &edges1) );
        }
      
        CMR_CALL( CMRfreeStackArray(cmr, &bipartition) );

        CMR_CALL( CMRfreeStackArray(cmr, &componentAuxiliaryNodes) );
        CMR_CALL( CMRfreeStackArray(cmr, &nodesComponent) );
        CMR_CALL( CMRgraphFree(cmr, &auxiliaryGraph) );
      }
    }
  }

  CMR_CALL( CMRfreeStackArray(cmr, &nodesCandidate) );

  return CMR_OKAY;
}

/**
 * \brief Find element in submatrix parallel to vector.
 */

static
CMR_ELEMENT findParallel(
  CMR_CHRMAT* matrix,
  size_t row,
  size_t numRows,
  size_t numColumns,
  long long* rowHashValues,
  long long* hashVector
)
{
  assert(matrix);
  assert(rowHashValues);

  long long hashValue = 0;
  size_t first = matrix->rowSlice[row];
  size_t beyond = matrix->rowSlice[row + 1];
  size_t countNonzeros = 0;
  for (size_t e = first; e < beyond; ++e)
  {
    size_t column = matrix->entryColumns[e];
    if (column < numColumns)
    {
      hashValue = projectSignedHash(hashValue + hashVector[column]);
      ++countNonzeros;
    }
    else
      break;
  }

  assert(countNonzeros >= 1);
  if (countNonzeros == 1)
    return CMRcolumnToElement(matrix->entryColumns[first]);

  for (size_t row2 = 0; row2 < numRows; ++row2)
  {
    if (rowHashValues[row2] != hashValue)
      continue;

    size_t first2 = matrix->rowSlice[row2];
    size_t beyond2 = matrix->rowSlice[row2 + 1];
    size_t e = first;
    size_t e2 = first2;
    bool isParallel = true;
    while (e < beyond && e2 < beyond2)
    {
      size_t column = matrix->entryColumns[e];
      size_t column2 = matrix->entryColumns[e2];
      if (column >= numColumns && column2 >= numColumns)
        break;
      if (column != column2)
      {
        isParallel = false;
        break;
      }
      ++e;
      ++e2;
    }
    if (isParallel)
      return CMRrowToElement(row2);
  }

  return 0;
}


/**
 * \brief Creates a hash vector to speed-up recognition of parallel vectors.
 */

static
CMR_ERROR createHashVector(
  CMR* cmr,                 /**< \ref CMR environment. */
  long long** phashVector,  /**< Pointer for storing the hash vector. */
  size_t size               /**< Size of hash vector. */
)
{
  assert(cmr);

  CMR_CALL( CMRallocStackArray(cmr, phashVector, size) );
  long long* hashVector = *phashVector;
  size_t h = 1;
  for (size_t e = 0; e < size; ++e)
  {
    hashVector[e] = h;
    h = projectSignedHash(3 * h);
  }

  return CMR_OKAY;
}

/**
 * \brief Update the hash values of rows/column of submatrix that is grown by a number of rows.
 */

static
CMR_ERROR updateHashValues(
  CMR_CHRMAT* matrix,         /**< Matrix. */
  long long* majorHashValues, /**< Map for hash values of major indices. */
  long long* minorHashValues, /**< Map for hash values of minor indices. */
  long long* hashVector,      /**< Hash vector. */
  size_t majorFirst,          /**< First new major index. */
  size_t majorBeyond,         /**< Last new major index plus 1. */
  size_t minorSize            /**< Number of minor indices in submatrix. */
)
{
  assert(matrix);
  assert(hashVector);

  for (size_t major = majorFirst; major < majorBeyond; ++major)
  {
    size_t first = matrix->rowSlice[major];
    size_t beyond = matrix->rowSlice[major + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      size_t minor = matrix->entryColumns[e];
      if (minor < minorSize)
      {
        majorHashValues[major] = projectSignedHash(majorHashValues[major] + hashVector[minor]);
        minorHashValues[minor] = projectSignedHash(minorHashValues[minor] + hashVector[major]);
      }
      else
        break;
    }
  }

  return CMR_OKAY;
}


/**
 * \brief Returns \c true if two edges \p e and \p f are adjacent.
 */

static
bool checkEdgesAdjacent(
  CMR_GRAPH* graph,         /**< Graph. */
  CMR_GRAPH_EDGE e,         /**< First edge. */
  CMR_GRAPH_EDGE f,         /**< Second edge. */
  CMR_GRAPH_NODE* pcommon,  /**< Pointer for storing the common endnode. */
  CMR_GRAPH_NODE* peOther,  /**< Pointer for storing the endnode of \p e that is not common. */
  CMR_GRAPH_NODE* pfOther   /**< Pointer for storing the endnode of \p f that is not common. */
)
{
  assert(graph);
  assert(pcommon);
  assert(peOther);
  assert(pfOther);

  CMR_GRAPH_NODE eNodes[2] = { CMRgraphEdgeU(graph, e), CMRgraphEdgeV(graph, e) };
  CMR_GRAPH_NODE fNodes[2] = { CMRgraphEdgeU(graph, f), CMRgraphEdgeV(graph, f) };
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      if (eNodes[i] == fNodes[j])
      {
        *pcommon = eNodes[i];
        *peOther = eNodes[1-i];
        *pfOther = fNodes[1-j];
        return true;
      }
    }
  }
  return false;
}

/**
 * \brief Extends \p graph for a submatrix augmented by 1 row and 1 column.
 */

static
CMR_ERROR addToGraph1Row1Column(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* graph,             /**< Empty graph to be filled. */
  CMR_GRAPH_EDGE* rowEdges,     /**< Array to be filled with map from rows to edges. */
  CMR_GRAPH_EDGE* columnEdges,  /**< Array to be filled with map from columns to edges. */
  size_t baseNumRows,           /**< Number of rows already processed. */
  size_t baseNumColumns,        /**< Number of columns already processed. */
  CMR_ELEMENT rowParallel,      /**< Element to which the row is parallel. */
  CMR_ELEMENT columnParallel,   /**< Element to which the column is parallel. */
  bool* pisGraphic              /**< Pointer for storing whether this extension was graphic. */
)
{
  assert(cmr);
  assert(graph);
  assert(baseNumRows >= 3);
  assert(baseNumColumns >= 3);
  assert(rowEdges);
  assert(columnEdges);
  assert(CMRelementIsValid(rowParallel));
  assert(CMRelementIsValid(columnParallel));

  CMR_GRAPH_EDGE rowEdge, columnEdge;
  if (CMRelementIsRow(rowParallel))
    rowEdge = rowEdges[CMRelementToRowIndex(rowParallel)];
  else
    rowEdge = columnEdges[CMRelementToColumnIndex(rowParallel)];
  CMRdbgMsg(12, "Existing row edge is {%ld,%ld}.\n", CMRgraphEdgeU(graph, rowEdge), CMRgraphEdgeV(graph, rowEdge));
  if (CMRelementIsRow(columnParallel))
    columnEdge = rowEdges[CMRelementToRowIndex(columnParallel)];
  else
    columnEdge = columnEdges[CMRelementToColumnIndex(columnParallel)];
  CMRdbgMsg(12, "Existing column edge is {%ld,%ld}.\n", CMRgraphEdgeU(graph, columnEdge),
    CMRgraphEdgeV(graph, columnEdge));

  CMR_GRAPH_NODE common, rowOther, columnOther;
  if (checkEdgesAdjacent(graph, rowEdge, columnEdge, &common, &rowOther, &columnOther))
  {
    *pisGraphic = true;
    CMR_GRAPH_NODE rowSplit;
    CMR_CALL( CMRgraphAddNode(cmr, graph, &rowSplit) );
    CMR_CALL( CMRgraphDeleteEdge(cmr, graph, rowEdge) );
    CMR_GRAPH_EDGE modifiedRowEdge, newRowEdge, newColumnEdge;
    CMR_CALL( CMRgraphAddEdge(cmr, graph, rowOther, rowSplit, &modifiedRowEdge) );
    assert(modifiedRowEdge == rowEdge);
    CMR_CALL( CMRgraphAddEdge(cmr, graph, rowSplit, common, &newRowEdge) );
    rowEdges[baseNumRows] = newRowEdge;
    CMRdbgMsg(12, "Edge corresponding to row is {%ld,%ld}.\n", CMRgraphEdgeU(graph, newRowEdge),
      CMRgraphEdgeV(graph, newRowEdge));
    CMR_CALL( CMRgraphAddEdge(cmr, graph, rowSplit, columnOther, &newColumnEdge) );
    columnEdges[baseNumColumns] = newColumnEdge;
    CMRdbgMsg(12, "Edge corresponding to column is {%ld,%ld}.\n", CMRgraphEdgeU(graph, newColumnEdge),
      CMRgraphEdgeV(graph, newColumnEdge));
  }
  else
    *pisGraphic = false;

  return CMR_OKAY;
}

/**
 * \brief Extends \p graph for a submatrix augmented by 2 rows and 1 column.
 */

static
CMR_ERROR addToGraph2Rows1Column(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* graph,             /**< Empty graph to be filled. */
  CMR_GRAPH_EDGE* rowEdges,     /**< Array to be filled with map from rows to edges. */
  CMR_GRAPH_EDGE* columnEdges,  /**< Array to be filled with map from columns to edges. */
  size_t baseNumRows,           /**< Number of rows already processed. */
  size_t baseNumColumns,        /**< Number of columns already processed. */    
  CMR_ELEMENT row1Parallel,     /**< Element to which row1 is parallel. */
  CMR_ELEMENT row2Parallel,     /**< Element to which row2 is parallel. */
  bool* pisGraphic              /**< Pointer for storing whether this extension was graphic. */
)
{
  assert(cmr);
  assert(graph);
  assert(baseNumRows >= 3);
  assert(baseNumColumns >= 3);
  assert(rowEdges);
  assert(columnEdges);
  assert(CMRelementIsValid(row1Parallel));
  assert(CMRelementIsValid(row2Parallel));

  CMR_GRAPH_EDGE row1Edge, row2Edge;
  if (CMRelementIsRow(row1Parallel))
    row1Edge = rowEdges[CMRelementToRowIndex(row1Parallel)];
  else
    row1Edge = columnEdges[CMRelementToColumnIndex(row1Parallel)];
  CMRdbgMsg(12, "Row1's edge is {%ld,%ld}.\n", CMRgraphEdgeU(graph, row1Edge), CMRgraphEdgeV(graph, row1Edge));
  if (CMRelementIsRow(row2Parallel))
    row2Edge = rowEdges[CMRelementToRowIndex(row2Parallel)];
  else
    row2Edge = columnEdges[CMRelementToColumnIndex(row2Parallel)];
  CMRdbgMsg(12, "Row2's edge is {%ld,%ld}.\n", CMRgraphEdgeU(graph, row2Edge), CMRgraphEdgeV(graph, row2Edge));

  CMR_GRAPH_NODE common, other1, other2;
  if (checkEdgesAdjacent(graph, row1Edge, row2Edge, &common, &other1, &other2))
  {
    *pisGraphic = true;

    CMR_GRAPH_NODE row1Split;
    CMR_CALL( CMRgraphAddNode(cmr, graph, &row1Split) );
    CMR_CALL( CMRgraphDeleteEdge(cmr, graph, row1Edge) );
    CMR_GRAPH_EDGE modifiedRow1Edge;
    CMR_CALL( CMRgraphAddEdge(cmr, graph, other1, row1Split, &modifiedRow1Edge) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, row1Split, common, &rowEdges[baseNumRows]) );
    assert(modifiedRow1Edge == row1Edge);

    CMRdbgMsg(12, "Row1's edge {%ld,%ld} is subdivided with new node %ld.\n", other1, common, row1Split);

    CMR_GRAPH_NODE row2Split;
    CMR_CALL( CMRgraphAddNode(cmr, graph, &row2Split) );
    CMR_CALL( CMRgraphDeleteEdge(cmr, graph, row2Edge) );
    CMR_GRAPH_EDGE modifiedRow2Edge;
    CMR_CALL( CMRgraphAddEdge(cmr, graph, other2, row2Split, &modifiedRow2Edge) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, row2Split, common, &rowEdges[baseNumRows+1]) );
    assert(modifiedRow2Edge == row2Edge);

    CMRdbgMsg(12, "Row2's edge {%ld,%ld} is subdivided with new node %ld.\n", other2, common, row2Split);

    CMR_CALL( CMRgraphAddEdge(cmr, graph, row1Split, row2Split, &columnEdges[baseNumColumns]) );
  }
  else
    *pisGraphic = false;

  return CMR_OKAY;
}

/**
 * \brief Extends \p graph for a submatrix augmented by 1 row and 2 columns.
 */

static
CMR_ERROR addToGraph1Row2Columns(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* graph,             /**< Empty graph to be filled. */
  CMR_GRAPH_EDGE* rowEdges,     /**< Array to be filled with map from rows to edges. */
  CMR_GRAPH_EDGE* columnEdges,  /**< Array to be filled with map from columns to edges. */
  size_t baseNumRows,           /**< Number of rows already processed. */
  size_t baseNumColumns,        /**< Number of columns already processed. */    
  CMR_ELEMENT column1Parallel,  /**< Element to which column1 is parallel. */
  CMR_ELEMENT column2Parallel,  /**< Element to which column2 is parallel. */
  bool* pisGraphic              /**< Pointer for storing whether this extension was graphic. */
)
{
  assert(cmr);
  assert(graph);
  assert(baseNumRows >= 3);
  assert(baseNumColumns >= 3);
  assert(rowEdges);
  assert(columnEdges);
  assert(CMRelementIsValid(column1Parallel));
  assert(CMRelementIsValid(column2Parallel));

  CMR_GRAPH_EDGE column1Edge, column2Edge;
  if (CMRelementIsRow(column1Parallel))
    column1Edge = rowEdges[CMRelementToRowIndex(column1Parallel)];
  else
    column1Edge = columnEdges[CMRelementToColumnIndex(column1Parallel)];
  CMRdbgMsg(12, "Column1's edge is {%ld,%ld}.\n", CMRgraphEdgeU(graph, column1Edge), CMRgraphEdgeV(graph, column1Edge));
  if (CMRelementIsRow(column2Parallel))
    column2Edge = rowEdges[CMRelementToRowIndex(column2Parallel)];
  else
    column2Edge = columnEdges[CMRelementToColumnIndex(column2Parallel)];
  CMRdbgMsg(12, "Column2's edge is {%ld,%ld}.\n", CMRgraphEdgeU(graph, column2Edge), CMRgraphEdgeV(graph, column2Edge));

  CMR_GRAPH_NODE common, other1, other2;
  if (checkEdgesAdjacent(graph, column1Edge, column2Edge, &common, &other1, &other2))
  {
    *pisGraphic = true;
    
    CMR_GRAPH_NODE newNode;
    CMR_CALL( CMRgraphAddNode(cmr, graph, &newNode) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, other1, newNode, &columnEdges[baseNumColumns]) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, other2, newNode, &columnEdges[baseNumColumns + 1]) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, common, newNode, &rowEdges[baseNumRows]) );
  }
  else
    *pisGraphic = false;

  return CMR_OKAY;
}

/**
 * \brief Extends \p graph for a submatrix augmented by 1 column.
 */

static
CMR_ERROR addToGraph1Column(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* graph,             /**< Empty graph to be filled. */
  CMR_GRAPH_EDGE* rowEdges,     /**< Array to be filled with map from rows to edges. */
  CMR_GRAPH_EDGE* columnEdges,  /**< Array to be filled with map from columns to edges. */
  size_t baseNumColumns,        /**< Number of columns already processed. */
  size_t* nonzeroRows,          /**< Array with rows containing a nonzero. */
  size_t numNonzeroRows,        /**< Length of \p nonzeroRows. */
  bool* pisGraphic              /**< Pointer for storing whether this extension was graphic. */
)
{
  assert(cmr);
  assert(graph);
  assert(baseNumColumns >= 3);
  assert(rowEdges);
  assert(columnEdges);

  size_t* nodeDegrees = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodeDegrees, CMRgraphMemNodes(graph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
    nodeDegrees[v] = 0;

  size_t countLeaves = 0;
  for (size_t i = 0; i < numNonzeroRows; ++i)
  {
    CMR_GRAPH_EDGE e = rowEdges[nonzeroRows[i]];
    size_t deg = ++nodeDegrees[CMRgraphEdgeU(graph, e)];
    if (deg == 1)
      ++countLeaves;
    else if (deg == 2)
      --countLeaves;
    deg = ++nodeDegrees[CMRgraphEdgeV(graph, e)];
    if (deg == 1)
      ++countLeaves;
    else if (deg == 2)
      --countLeaves;
  }

  *pisGraphic = countLeaves == 2;

  if (*pisGraphic)
  {
    CMR_GRAPH_NODE nodes[2] = { -1, -1 };
    for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
    {
      if (nodeDegrees[v] == 1)
      {
        if (nodes[0] < 0)
          nodes[0] = v;
        else
          nodes[1] = v;
      }
    }
    CMR_CALL( CMRgraphAddEdge(cmr, graph, nodes[0], nodes[1], &columnEdges[baseNumColumns]) );
  }

  CMR_CALL( CMRfreeStackArray(cmr, &nodeDegrees) );

  return CMR_OKAY;
}

/**
 * \brief Creates the wheel graph for a wheel submatrix.
 */

static
CMR_ERROR createWheel(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_GRAPH* graph,           /**< Empty graph to be filled. */
  CMR_CHRMAT* matrix,         /**< Matrix. */
  CMR_CHRMAT* transpose,      /**< Transpose of \p matrix. */
  size_t wheelSize,           /**< Size of wheel. */
  CMR_GRAPH_EDGE* rowEdges,   /**< Array to be filled with map from rows to edges. */
  CMR_GRAPH_EDGE* columnEdges /**< Array to be filled with map from columns to edges. */
)
{
  assert(graph);
  assert(matrix);
  assert(transpose);
  assert(wheelSize <= matrix->numRows);
  assert(wheelSize <= matrix->numColumns);
  assert(rowEdges);
  assert(columnEdges);

  CMRdbgMsg(8, "Creating wheel graph W_%ld for first minor.\n", wheelSize);

  /* Check which row contains 3 indices (if any). */
  size_t rowWithThree = SIZE_MAX;
  for (size_t row = 0; row < wheelSize; ++row)
  {
    size_t count = 0;
    size_t first = matrix->rowSlice[row];
    size_t beyond = matrix->rowSlice[row + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      if (matrix->entryColumns[e] < wheelSize)
        ++count;
      else
        break;
    }
    assert(count == 2 || count == 3);
    if (count == 3)
    {
      assert(rowWithThree == SIZE_MAX);
      rowWithThree = row;
    }
  }

  /* Check which column contains 3 indices (if any). */
  size_t columnWithThree = SIZE_MAX;
  for (size_t column = 0; column < wheelSize; ++column)
  {
    size_t count = 0;
    size_t first = transpose->rowSlice[column];
    size_t beyond = transpose->rowSlice[column + 1];
    for (size_t e = first; e < beyond; ++e)
    {
      if (transpose->entryColumns[e] < wheelSize)
        ++count;
      else
        break;
    }
    assert(count == 2 || count == 3);
    if (count == 3)
    {
      assert(columnWithThree == SIZE_MAX);
      columnWithThree = column;
    }
  }

  assert((rowWithThree == SIZE_MAX && columnWithThree == SIZE_MAX)
    || (rowWithThree < SIZE_MAX && columnWithThree < SIZE_MAX));

  CMR_GRAPH_NODE centerNode, firstRimNode;
  CMR_CALL( CMRgraphAddNode(cmr, graph, &centerNode) );
  CMR_CALL( CMRgraphAddNode(cmr, graph, &firstRimNode) );
  CMR_GRAPH_NODE lastRimNode = firstRimNode;

  size_t lastRow = 0;
  size_t lastColumn = matrix->entryColumns[matrix->rowSlice[0]];
  size_t nextRow = SIZE_MAX;
  size_t nextColumn = SIZE_MAX;
  while (nextRow)
  {
    size_t e = matrix->rowSlice[lastRow];

    /* Find next column. */
    if (lastRow == rowWithThree)
    {
      nextColumn = matrix->entryColumns[e];
      if (nextColumn == lastColumn || nextColumn == columnWithThree)
        nextColumn = matrix->entryColumns[e + 1];
      if (nextColumn == lastColumn || nextColumn == columnWithThree)
        nextColumn = matrix->entryColumns[e + 2];
    }
    else
    {
      nextColumn = matrix->entryColumns[e];
      if (nextColumn == lastColumn)
        nextColumn = matrix->entryColumns[e + 1];
    }

    e = transpose->rowSlice[nextColumn];

    /* Find next row. */
    if (nextColumn == columnWithThree)
    {
      nextRow = transpose->entryColumns[e];
      if (nextRow == lastRow || nextRow == rowWithThree)
        nextRow = transpose->entryColumns[e + 1];
      if (nextRow == lastRow || nextRow == rowWithThree)
        nextRow = transpose->entryColumns[e + 2];
    }
    else
    {
      nextRow = transpose->entryColumns[e];
      if (nextRow == lastRow)
        nextRow = transpose->entryColumns[e + 1];
    }

    CMRdbgMsg(10, "Last column is y%ld, last row is x%ld, new column is y%ld and new row is x%ld\n", lastColumn+1,
      lastRow+1, nextColumn+1, nextRow+1);

    CMR_GRAPH_NODE nextRimNode;
    CMR_GRAPH_EDGE rimEdge;
    CMR_GRAPH_EDGE spokeEdge;
    if (nextRow == 0)
      nextRimNode = firstRimNode;
    else
      CMR_CALL( CMRgraphAddNode(cmr, graph, &nextRimNode) );
    CMR_CALL( CMRgraphAddEdge(cmr, graph, lastRimNode, nextRimNode, &rimEdge) );

    CMRdbgMsg(10, "Added rim {%ld,%ld} for column y%ld.\n", lastRimNode, nextRimNode, lastColumn+1);

    CMR_CALL( CMRgraphAddEdge(cmr, graph, centerNode, nextRimNode, &spokeEdge) );

    CMRdbgMsg(10, "Added spoke {%ld,%ld} for row x%ld.\n", centerNode, nextRimNode, lastRow+1);

    if (rowWithThree < SIZE_MAX && lastRow != rowWithThree && nextRow != rowWithThree)
    {
      columnEdges[lastColumn] = spokeEdge;
      rowEdges[lastRow] = rimEdge;
      CMRdbgMsg(10, "Spoke is assigned to y%ld and rim to x%ld.\n", lastColumn+1, lastRow+1);
    }
    else
    {
      columnEdges[lastColumn] = rimEdge;
      rowEdges[lastRow] = spokeEdge;
      CMRdbgMsg(10, "Spoke is assigned to x%ld and rim to y%ld.\n", lastRow+1, lastColumn+1);
    }

    lastRimNode = nextRimNode;
    lastRow = nextRow;
    lastColumn = nextColumn;
  }

  return CMR_OKAY;
}


CMR_ERROR CMRregularSequenceGraphic(CMR* cmr, CMR_CHRMAT* matrix, CMR_CHRMAT* transpose, size_t lengthSequence,
  size_t* sequenceNumRows, size_t* sequenceNumColumns, size_t* plastGraphicMinor, CMR_GRAPH** pgraph,
  CMR_ELEMENT** pedgeElements, CMR_SEYMOUR_STATS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(transpose);
  assert(sequenceNumRows);
  assert(sequenceNumColumns);
  assert(plastGraphicMinor);
  assert(pgraph);
  assert(!*pgraph);
  assert(pedgeElements);
  assert(!*pedgeElements);

  CMRdbgMsg(8, "Testing sequence for (co)graphicness.\n");

  clock_t time = clock();

  CMR_CALL( CMRgraphCreateEmpty(cmr, pgraph, matrix->numRows, matrix->numRows + matrix->numColumns) );
  CMR_GRAPH* graph = *pgraph;

  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );
  CMR_GRAPH_EDGE* rowEdges = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowEdges, matrix->numRows) );
  CMR_GRAPH_EDGE* columnEdges = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnEdges, matrix->numColumns) );
  long long* rowHashValues = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowHashValues, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    rowHashValues[row] = 0;
  long long* columnHashValues = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnHashValues, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnHashValues[column] = 0;

  /* Create graph for first minor. */

  assert(sequenceNumRows[0] == sequenceNumColumns[0]);
  CMR_CALL( createWheel(cmr, graph, matrix, transpose, sequenceNumRows[0], rowEdges, columnEdges) );
  *plastGraphicMinor = 0;

  CMR_CALL( updateHashValues(matrix, rowHashValues, columnHashValues, hashVector, 0, sequenceNumRows[0],
    sequenceNumColumns[0]) );

  size_t extensionTimeFactor = lengthSequence / 100 + 1;
  for (size_t extension = 1; extension < lengthSequence; ++extension)
  { 
    if ((extension % extensionTimeFactor == 0) && (clock() - time) * 1.0 / CLOCKS_PER_SEC > timeLimit)
    {
      CMR_CALL( CMRfreeStackArray(cmr, &columnHashValues) );
      CMR_CALL( CMRfreeStackArray(cmr, &rowHashValues) );
      CMR_CALL( CMRfreeStackArray(cmr, &columnEdges) );
      CMR_CALL( CMRfreeStackArray(cmr, &rowEdges) );
      CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );
      return CMR_ERROR_TIMEOUT;
    }
    
    size_t newRows = sequenceNumRows[extension] - sequenceNumRows[extension-1];
    size_t newColumns = sequenceNumColumns[extension] - sequenceNumColumns[extension-1];

    CMRdbgMsg(10, "Processing extension step %ld with %ld new rows and %ld new columns.\n", extension, newRows,
      newColumns);

    bool isGraphic;
    if (newRows == 1 && newColumns == 1)
    {
      CMR_ELEMENT rowParallel = findParallel(matrix, sequenceNumRows[extension-1], sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowHashValues, hashVector);
      CMR_ELEMENT columnParallel = CMRelementTranspose(findParallel(transpose, sequenceNumColumns[extension-1],
        sequenceNumColumns[extension-1], sequenceNumRows[extension-1], columnHashValues, hashVector));

      CMRdbgMsg(10, "The new row is parallel to %c%ld", CMRelementIsRow(rowParallel) ? 'x' : 'y',
        CMRelementIsRow(rowParallel) ? CMRelementToRowIndex(rowParallel) + 1 : CMRelementToColumnIndex(rowParallel) + 1);
      CMRdbgMsg(0, " and the new column is parallel to %c%ld.\n", CMRelementIsRow(columnParallel) ? 'x' : 'y',
        CMRelementIsRow(columnParallel) ? CMRelementToRowIndex(columnParallel) + 1
        : CMRelementToColumnIndex(columnParallel) + 1);

      CMR_CALL( addToGraph1Row1Column(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowParallel, columnParallel, &isGraphic) );
    }
    else if (newRows == 2 && newColumns == 1)
    {
      CMR_ELEMENT row1Parallel = findParallel(matrix, sequenceNumRows[extension-1], sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowHashValues, hashVector);
      CMR_ELEMENT row2Parallel = findParallel(matrix, sequenceNumRows[extension-1] + 1, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowHashValues, hashVector);

      CMRdbgMsg(10, "Row 1 is parallel to %s", CMRelementString(row1Parallel, 0));
      CMRdbgMsg(0, " and row 2 is parallel to %s.\n", CMRelementString(row2Parallel, 0));

      CMR_CALL( addToGraph2Rows1Column(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], row1Parallel, row2Parallel, &isGraphic) );
    }
    else if (newRows == 1 && newColumns == 2)
    {
      CMR_ELEMENT column1Parallel = CMRelementTranspose(findParallel(transpose, sequenceNumColumns[extension-1],
        sequenceNumColumns[extension-1], sequenceNumRows[extension-1], columnHashValues, hashVector));
      CMR_ELEMENT column2Parallel = CMRelementTranspose(findParallel(transpose, sequenceNumColumns[extension-1] + 1,
        sequenceNumColumns[extension-1], sequenceNumRows[extension-1], columnHashValues, hashVector));

      CMRdbgMsg(10, "Column 1 is parallel to %s", CMRelementString(column1Parallel, 0));
      CMRdbgMsg(0, " and column 2 is parallel to %s.\n", CMRelementString(column2Parallel, 0));

      CMR_CALL( addToGraph1Row2Columns(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], column1Parallel, column2Parallel, &isGraphic) );
    }
    else if (newRows == 0 && newColumns == 1)
    {
      size_t first = transpose->rowSlice[sequenceNumColumns[extension-1]];
      size_t beyond = transpose->rowSlice[sequenceNumColumns[extension-1] + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        if (transpose->entryColumns[e] >= sequenceNumRows[extension-1])
          beyond = e;
      }
      CMR_CALL( addToGraph1Column(cmr, graph, rowEdges, columnEdges, sequenceNumColumns[extension-1],
        &transpose->entryColumns[first], beyond-first, &isGraphic) );
    }
    else
    {
      assert(newRows == 1 && newColumns == 0);

      size_t first = matrix->rowSlice[sequenceNumRows[extension-1]];
      size_t beyond = matrix->rowSlice[sequenceNumRows[extension-1] + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        if (matrix->entryColumns[e] >= sequenceNumColumns[extension-1])
          beyond = e;
      }
      CMR_CALL( addToGraph1Row(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        &matrix->entryColumns[first], beyond-first, &isGraphic) );
    }

    if (isGraphic)
      *plastGraphicMinor = extension;
    else
      break;

    CMR_CALL( updateHashValues(matrix, rowHashValues, columnHashValues, hashVector, sequenceNumRows[extension-1],
      sequenceNumRows[extension], sequenceNumColumns[extension - 1]) );
    CMR_CALL( updateHashValues(transpose, columnHashValues, rowHashValues, hashVector, sequenceNumColumns[extension-1],
      sequenceNumColumns[extension], sequenceNumRows[extension]) );
  }

  if (*plastGraphicMinor == lengthSequence - 1)
  {
    CMR_CALL( CMRallocBlockArray(cmr, pedgeElements, matrix->numRows + matrix->numColumns) );
    CMR_ELEMENT* edgeElements = *pedgeElements;
    for (size_t e = 0; e < matrix->numRows + matrix->numColumns; ++e)
      edgeElements[e] = 0;

    for (size_t row = 0; row < matrix->numRows; ++row)
      edgeElements[rowEdges[row]] = CMRrowToElement(row);
    for (size_t column = 0; column < matrix->numColumns; ++column)
      edgeElements[columnEdges[column]] = CMRcolumnToElement(column);
  }
  else
  {
    CMR_CALL( CMRgraphFree(cmr, pgraph) );
  }

  if (stats)
  {
    stats->sequenceGraphicCount++;
    stats->sequenceGraphicTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &columnHashValues) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowHashValues) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnEdges) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowEdges) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );

  return CMR_OKAY;
}

/**
 * \brief Tests each minor of the sequence of nested 3-connected minors for graphicness.
 */

static
CMR_ERROR sequenceGraphicness(
  CMR* cmr,                     /**< \ref CMR environment. */
  DecompositionTask* task,      /**< Task to be processed; already removed from the list of unprocessed tasks. */
  CMR_CHRMAT* matrix,           /**< Matrix that displays the nested minor sequences. */
  CMR_CHRMAT* transpose,        /**< Transpose of \p matrix. */
  size_t length,                /**< Length of sequence of nested minors. */
  size_t* sequenceNumRows,      /**< Number of rows of sequence of nested minors. */
  size_t* sequenceNumColumns,   /**< Number of columns of sequence of nested minors. */
  bool cographicness,           /**< Whether we actually received the transpose as input. */
  CMR_GRAPH** pgraph,           /**< Pointer for storing the constructed graph. */
  CMR_ELEMENT** pedgeElements,  /**< Pointer for storing the mapping from edges to elements. */
  size_t* plastGraphicMinor     /**< Pointer for storing the last graphic minor. */
)
{
  assert(cmr);
  assert(task);
  assert(matrix);
  assert(transpose);
  assert(sequenceNumRows);
  assert(sequenceNumColumns);
  assert(pgraph);
  assert(pedgeElements);
  assert(plastGraphicMinor);

  CMRdbgMsg(6, "Testing sequence of nested minors of length %zu for %sgraphicness.\n", length,
    cographicness ? "co" : "");

  CMR_CALL( CMRgraphCreateEmpty(cmr, pgraph, matrix->numRows, matrix->numRows + matrix->numColumns) );
  CMR_GRAPH* graph = *pgraph;

  long long* hashVector = NULL;
  CMR_CALL( createHashVector(cmr, &hashVector,
    matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns) );
  CMR_GRAPH_EDGE* rowEdges = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowEdges, matrix->numRows) );
  CMR_GRAPH_EDGE* columnEdges = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnEdges, matrix->numColumns) );
  long long* rowHashValues = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &rowHashValues, matrix->numRows) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    rowHashValues[row] = 0;
  long long* columnHashValues = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &columnHashValues, matrix->numColumns) );
  for (size_t column = 0; column < matrix->numColumns; ++column)
    columnHashValues[column] = 0;

  /* Create graph for first minor. */

  assert(sequenceNumRows[0] == sequenceNumColumns[0]);
  CMR_CALL( createWheel(cmr, graph, matrix, transpose, sequenceNumRows[0], rowEdges, columnEdges) );
  *plastGraphicMinor = 0;

  CMR_CALL( updateHashValues(matrix, rowHashValues, columnHashValues, hashVector, 0, sequenceNumRows[0],
    sequenceNumColumns[0]) );

//   size_t extensionTimeFactor = length / 100 + 1;
  for (size_t extension = 1; extension < length; ++extension)
  {
//     if ((extension % extensionTimeFactor == 0) && (clock() - time) * 1.0 / CLOCKS_PER_SEC > timeLimit)
//     {
//       CMR_CALL( CMRfreeStackArray(cmr, &columnHashValues) );
//       CMR_CALL( CMRfreeStackArray(cmr, &rowHashValues) );
//       CMR_CALL( CMRfreeStackArray(cmr, &columnEdges) );
//       CMR_CALL( CMRfreeStackArray(cmr, &rowEdges) );
//       CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );
//       return CMR_ERROR_TIMEOUT;
//     }

    size_t newRows = sequenceNumRows[extension] - sequenceNumRows[extension-1];
    size_t newColumns = sequenceNumColumns[extension] - sequenceNumColumns[extension-1];

    CMRdbgMsg(8, "Processing extension step %zu with %zu new rows and %zu new columns.\n", extension, newRows,
      newColumns);

    bool isGraphic;
    if (newRows == 1 && newColumns == 1)
    {
      CMR_ELEMENT rowParallel = findParallel(matrix, sequenceNumRows[extension-1], sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowHashValues, hashVector);
      CMR_ELEMENT columnParallel = CMRelementTranspose(findParallel(transpose, sequenceNumColumns[extension-1],
        sequenceNumColumns[extension-1], sequenceNumRows[extension-1], columnHashValues, hashVector));

      CMRdbgMsg(10, "The new row is parallel to %c%zu", CMRelementIsRow(rowParallel) ? 'x' : 'y',
        CMRelementIsRow(rowParallel) ? CMRelementToRowIndex(rowParallel) + 1 : CMRelementToColumnIndex(rowParallel) + 1);
      CMRdbgMsg(0, " and the new column is parallel to %c%zu.\n", CMRelementIsRow(columnParallel) ? 'x' : 'y',
        CMRelementIsRow(columnParallel) ? CMRelementToRowIndex(columnParallel) + 1
        : CMRelementToColumnIndex(columnParallel) + 1);

      CMR_CALL( addToGraph1Row1Column(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowParallel, columnParallel, &isGraphic) );
    }
    else if (newRows == 2 && newColumns == 1)
    {
      CMR_ELEMENT row1Parallel = findParallel(matrix, sequenceNumRows[extension-1], sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowHashValues, hashVector);
      CMR_ELEMENT row2Parallel = findParallel(matrix, sequenceNumRows[extension-1] + 1, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], rowHashValues, hashVector);

      CMRdbgMsg(10, "Row 1 is parallel to %s", CMRelementString(row1Parallel, 0));
      CMRdbgMsg(0, " and row 2 is parallel to %s.\n", CMRelementString(row2Parallel, 0));

      CMR_CALL( addToGraph2Rows1Column(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], row1Parallel, row2Parallel, &isGraphic) );
    }
    else if (newRows == 1 && newColumns == 2)
    {
      CMR_ELEMENT column1Parallel = CMRelementTranspose(findParallel(transpose, sequenceNumColumns[extension-1],
        sequenceNumColumns[extension-1], sequenceNumRows[extension-1], columnHashValues, hashVector));
      CMR_ELEMENT column2Parallel = CMRelementTranspose(findParallel(transpose, sequenceNumColumns[extension-1] + 1,
        sequenceNumColumns[extension-1], sequenceNumRows[extension-1], columnHashValues, hashVector));

      CMRdbgMsg(10, "Column 1 is parallel to %s", CMRelementString(column1Parallel, 0));
      CMRdbgMsg(0, " and column 2 is parallel to %s.\n", CMRelementString(column2Parallel, 0));

      CMR_CALL( addToGraph1Row2Columns(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        sequenceNumColumns[extension-1], column1Parallel, column2Parallel, &isGraphic) );
    }
    else if (newRows == 0 && newColumns == 1)
    {
      size_t first = transpose->rowSlice[sequenceNumColumns[extension-1]];
      size_t beyond = transpose->rowSlice[sequenceNumColumns[extension-1] + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        if (transpose->entryColumns[e] >= sequenceNumRows[extension-1])
          beyond = e;
      }
      CMR_CALL( addToGraph1Column(cmr, graph, rowEdges, columnEdges, sequenceNumColumns[extension-1],
        &transpose->entryColumns[first], beyond-first, &isGraphic) );
    }
    else
    {
      assert(newRows == 1 && newColumns == 0);

      size_t first = matrix->rowSlice[sequenceNumRows[extension-1]];
      size_t beyond = matrix->rowSlice[sequenceNumRows[extension-1] + 1];
      for (size_t e = first; e < beyond; ++e)
      {
        if (matrix->entryColumns[e] >= sequenceNumColumns[extension-1])
          beyond = e;
      }
      CMR_CALL( addToGraph1Row(cmr, graph, rowEdges, columnEdges, sequenceNumRows[extension-1],
        &matrix->entryColumns[first], beyond-first, &isGraphic) );
    }

    if (isGraphic)
      *plastGraphicMinor = extension;
    else
      break;

    CMR_CALL( updateHashValues(matrix, rowHashValues, columnHashValues, hashVector, sequenceNumRows[extension-1],
      sequenceNumRows[extension], sequenceNumColumns[extension - 1]) );
    CMR_CALL( updateHashValues(transpose, columnHashValues, rowHashValues, hashVector, sequenceNumColumns[extension-1],
      sequenceNumColumns[extension], sequenceNumRows[extension]) );
  }

  if (*plastGraphicMinor == length - 1)
  {
    CMR_CALL( CMRallocBlockArray(cmr, pedgeElements, matrix->numRows + matrix->numColumns) );
    CMR_ELEMENT* edgeElements = *pedgeElements;
    for (size_t e = 0; e < matrix->numRows + matrix->numColumns; ++e)
      edgeElements[e] = 0;

    if (cographicness)
    {
      for (size_t row = 0; row < matrix->numRows; ++row)
        edgeElements[rowEdges[row]] = task->node->nestedMinorsColumnsOriginal[row];
      for (size_t column = 0; column < matrix->numColumns; ++column)
        edgeElements[columnEdges[column]] = task->node->nestedMinorsRowsOriginal[column];
    }
    else
    {
      for (size_t row = 0; row < matrix->numRows; ++row)
        edgeElements[rowEdges[row]] = task->node->nestedMinorsRowsOriginal[row];
      for (size_t column = 0; column < matrix->numColumns; ++column)
        edgeElements[columnEdges[column]] = task->node->nestedMinorsColumnsOriginal[column];
    }
  }
  else
  {
    CMR_CALL( CMRgraphFree(cmr, pgraph) );
  }

  if (task->stats)
  {
    task->stats->sequenceGraphicCount++;
//     task->stats->sequenceGraphicTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  }

  CMR_CALL( CMRfreeStackArray(cmr, &columnHashValues) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowHashValues) );
  CMR_CALL( CMRfreeStackArray(cmr, &columnEdges) );
  CMR_CALL( CMRfreeStackArray(cmr, &rowEdges) );
  CMR_CALL( CMRfreeStackArray(cmr, &hashVector) );

  CMRdbgMsg(8, "Sequence of nested minors is %sgraphic up to and including minor 0 <= %zu <= %zu.\n",
    cographicness ? "co" : "", *plastGraphicMinor, length - 1);

  return CMR_OKAY;
}

CMR_ERROR CMRregularityNestedMinorSequenceGraphicness(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMR_SEYMOUR_NODE* dec = task->node;
  assert(dec);

  CMR_GRAPH* graph = NULL;
  CMR_ELEMENT* edgeElements = NULL;

  CMR_CALL( sequenceGraphicness(cmr, task, dec->nestedMinorsMatrix, dec->nestedMinorsTranspose, dec->nestedMinorsLength,
    dec->nestedMinorsSequenceNumRows, dec->nestedMinorsSequenceNumColumns, false, &graph, &edgeElements,
    &dec->nestedMinorsLastGraphic) );

  if (dec->nestedMinorsLastGraphic + 1 == dec->nestedMinorsLength)
  {
    dec->type = (dec->type == CMR_SEYMOUR_NODE_TYPE_COGRAPH) ? CMR_SEYMOUR_NODE_TYPE_PLANAR : CMR_SEYMOUR_NODE_TYPE_GRAPH;
    CMRdbgMsg(8, "Whole sequence is %s.\n", dec->type == CMR_SEYMOUR_NODE_TYPE_PLANAR ? "planar" : "graphic");
    dec->graph = graph;
    dec->graphForest = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &dec->graphForest, dec->matrix->numRows) );
    dec->graphCoforest = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &dec->graphCoforest, dec->matrix->numColumns) );
    dec->graphArcsReversed = NULL;
    if (dec->isTernary)
      CMR_CALL( CMRallocBlockArray(cmr, &dec->graphArcsReversed, CMRgraphMemEdges(dec->graph)) );

    assert(CMRgraphNumEdges(graph) == dec->matrix->numRows + dec->matrix->numColumns);
    for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
      iter = CMRgraphEdgesNext(graph, iter))
    {
      CMR_GRAPH_EDGE edge = CMRgraphEdgesEdge(graph, iter);
      CMR_ELEMENT element = edgeElements[edge];
      if (CMRelementIsRow(element))
        dec->graphForest[CMRelementToRowIndex(element)] = edge;
      else
        dec->graphCoforest[CMRelementToColumnIndex(element)] = edge;
    }

    CMR_CALL( CMRfreeBlockArray(cmr, &edgeElements) );

    if (dec->isTernary)
    {
      /* Check Camion signs. */

      if (dec->transpose == NULL)
        CMR_CALL( CMRchrmatTranspose(cmr, dec->matrix, &dec->transpose) );

      bool isCamionSigned;
      CMR_SUBMAT* violatorSubmatrix = NULL;
      CMR_CALL( CMRcamionCographicOrient(cmr, dec->transpose, dec->graph, dec->graphForest,
        dec->graphCoforest, dec->graphArcsReversed, &isCamionSigned, &violatorSubmatrix, NULL) );

      if (violatorSubmatrix)
      {
        CMR_CALL( CMRsubmatTranspose(violatorSubmatrix) );

        CMRdbgMsg(8, "-> %zux%zu submatrix with bad determinant.\n", violatorSubmatrix->numRows,
          violatorSubmatrix->numColumns);

        CMR_CALL( CMRseymourUpdateViolator(cmr, dec, violatorSubmatrix) );
        assert(dec->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

        CMR_CALL( CMRgraphFree(cmr, &dec->graph) );
        CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphForest) );
        CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphCoforest) );
        CMR_CALL( CMRfreeBlockArray(cmr, &dec->graphArcsReversed) );
        queue->foundIrregularity = true;
      }
    }

    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }
  else
  {
    dec->graphicness = -1;
    queue->foundNongraphicness = true;

    /* Add task back to list of unprocessed tasks to test cographicness. */
    CMRregularityQueueAdd(queue, task);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRregularityNestedMinorSequenceCographicness(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMR_SEYMOUR_NODE* dec = task->node;
  assert(dec);

  CMR_GRAPH* cograph = NULL;
  CMR_ELEMENT* edgeElements = NULL;

  CMR_CALL( sequenceGraphicness(cmr, task, dec->nestedMinorsTranspose, dec->nestedMinorsMatrix, dec->nestedMinorsLength,
    dec->nestedMinorsSequenceNumColumns, dec->nestedMinorsSequenceNumRows, true, &cograph, &edgeElements,
    &dec->nestedMinorsLastCographic) );

  if (dec->nestedMinorsLastCographic + 1 == dec->nestedMinorsLength)
  {
    dec->type = (dec->type == CMR_SEYMOUR_NODE_TYPE_GRAPH) ? CMR_SEYMOUR_NODE_TYPE_PLANAR : CMR_SEYMOUR_NODE_TYPE_COGRAPH;
    CMRdbgMsg(8, "Whole sequence is %s.\n", dec->type == CMR_SEYMOUR_NODE_TYPE_PLANAR ? "planar" : "cographic");
    dec->cograph = cograph;
    dec->cographForest = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &dec->cographForest, dec->matrix->numColumns) );
    dec->cographCoforest = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &dec->cographCoforest, dec->matrix->numRows) );
    dec->graphArcsReversed = NULL;
    if (dec->isTernary)
      CMR_CALL( CMRallocBlockArray(cmr, &dec->cographArcsReversed, CMRgraphMemEdges(dec->cograph)) );

    assert(CMRgraphNumEdges(cograph) == dec->matrix->numRows + dec->matrix->numColumns);
    for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(cograph); CMRgraphEdgesValid(cograph, iter);
      iter = CMRgraphEdgesNext(cograph, iter))
    {
      CMR_GRAPH_EDGE edge = CMRgraphEdgesEdge(cograph, iter);
      CMR_ELEMENT element = edgeElements[edge];
      if (CMRelementIsColumn(element))
        dec->cographForest[CMRelementToColumnIndex(element)] = edge;
      else
        dec->cographCoforest[CMRelementToRowIndex(element)] = edge;
    }

    CMR_CALL( CMRfreeBlockArray(cmr, &edgeElements) );

    if (dec->isTernary)
    {
      /* Check Camion signs. */

      bool isCamionSigned;
      CMR_SUBMAT* violatorSubmatrix = NULL;
      CMR_CALL( CMRcamionCographicOrient(cmr, dec->matrix, dec->cograph, dec->cographForest,
        dec->cographCoforest, dec->cographArcsReversed, &isCamionSigned, &violatorSubmatrix, NULL) );

      if (violatorSubmatrix)
      {
        CMRdbgMsg(8, "-> %zux%zu submatrix with bad determinant.\n", violatorSubmatrix->numRows,
          violatorSubmatrix->numColumns);

        CMR_CALL( CMRseymourUpdateViolator(cmr, dec, violatorSubmatrix) );

        CMR_CALL( CMRgraphFree(cmr, &dec->cograph) );
        CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographForest) );
        CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographCoforest) );
        CMR_CALL( CMRfreeBlockArray(cmr, &dec->cographArcsReversed) );
        queue->foundIrregularity = true;
      }
    }

    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }
  else
  {
    dec->cographicness = -1;
    queue->foundNoncographicness = true;

    /* Add task back to list of unprocessed tasks to search for 3-separations. */
    CMRregularityQueueAdd(queue, task);
  }

  return CMR_OKAY;
}


CMR_ERROR CMRregularityTestGraphicness(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMR_SEYMOUR_NODE* dec = task->node;
  assert(dec);

#if defined(CMR_DEBUG)
  CMRdbgMsg(6, "Testing the following matrix for %s:\n", dec->isTernary ? "being network" : "graphicness");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  if (!dec->transpose)
  {
    assert(dec->matrix);
    CMR_CALL( CMRchrmatTranspose(cmr, dec->matrix, &dec->transpose) );
  }

  double remainingTime = task->timeLimit - (clock() - task->startClock) * 1.0 / CLOCKS_PER_SEC;
  bool isGraphic;
  if (dec->isTernary)
  {
    CMR_SUBMAT* violatorSubmatrix = NULL;
    bool supportGraphic;
    CMR_ERROR error = CMRnetworkTestTranspose(cmr, dec->transpose, &isGraphic, &supportGraphic, &dec->graph,
      &dec->graphForest, &dec->graphCoforest, &dec->graphArcsReversed, &violatorSubmatrix,
      task->stats ? &task->stats->network : NULL, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
      return error;
    CMR_CALL( error );

    if (violatorSubmatrix)
    {
      if (supportGraphic)
      {
        CMR_CALL( CMRsubmatTranspose(violatorSubmatrix) );

        CMRdbgMsg(8, "-> %zux%zu submatrix with bad determinant.\n", violatorSubmatrix->numRows,
          violatorSubmatrix->numColumns);

        CMR_CALL( CMRseymourUpdateViolator(cmr, dec, violatorSubmatrix) );
        assert(dec->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

        CMR_CALL( CMRregularityTaskFree(cmr, &task) );
        queue->foundIrregularity = true;
        return CMR_OKAY;
      }

      CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
    }
  }
  else
  {
    CMR_ERROR error = CMRgraphicTestTranspose(cmr, dec->transpose, &isGraphic, &dec->graph, &dec->graphForest,
      &dec->graphCoforest, NULL, task->stats ? &task->stats->graphic : NULL, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
      return error;
    CMR_CALL( error );
  }

  CMRdbgMsg(8, "-> %s%s\n", isGraphic ? "" : "NOT ", dec->isTernary ? "network" : "graphic");

  dec->graphicness = isGraphic ? 1 : -1;
  if (isGraphic)
    dec->type = (dec->type == CMR_SEYMOUR_NODE_TYPE_COGRAPH) ? CMR_SEYMOUR_NODE_TYPE_PLANAR : CMR_SEYMOUR_NODE_TYPE_GRAPH;
  else
    queue->foundNongraphicness = true;

  if ((isGraphic && (!task->params->planarityCheck || dec->cographicness)) || dec->cographicness > 0)
  {
    /* Task is done. */
    CMRdbgMsg(8, "Marking task as complete.\n");
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }
  else
  {
    /* Re-insert task. */
    CMRregularityQueueAdd(queue, task);
  }

  return CMR_OKAY;
}

CMR_ERROR CMRregularityTestCographicness(CMR* cmr, DecompositionTask* task, DecompositionQueue* queue)
{
  assert(cmr);
  assert(task);
  assert(queue);

  CMR_SEYMOUR_NODE* dec = task->node;
  assert(dec);

#if defined(CMR_DEBUG)
  CMRdbgMsg(6, "Testing the following matrix for %s:\n", dec->isTernary ? "being conetwork" : "cographicness");
  CMR_CALL( CMRchrmatPrintDense(cmr, dec->matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  if (!dec->matrix)
  {
    assert(dec->transpose);
    CMR_CALL( CMRchrmatTranspose(cmr, dec->transpose, &dec->matrix) );
  }

  double remainingTime = task->timeLimit - (clock() - task->startClock) * 1.0 / CLOCKS_PER_SEC;
  bool isCographic;
  if (dec->isTernary)
  {
    CMR_SUBMAT* violatorSubmatrix = NULL;
    bool supportCographic;

    CMR_ERROR error = CMRnetworkTestTranspose(cmr, dec->matrix, &isCographic, &supportCographic, &dec->cograph,
      &dec->cographForest, &dec->cographCoforest, &dec->cographArcsReversed, &violatorSubmatrix,
      task->stats ? &task->stats->network : NULL, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
      return error;
    CMR_CALL( error );

    if (violatorSubmatrix)
    {
      if (supportCographic)
      {
        CMRdbgMsg(8, "-> %zux%zu submatrix with bad determinant.\n", violatorSubmatrix->numRows,
          violatorSubmatrix->numColumns);

        CMR_CALL( CMRseymourUpdateViolator(cmr, dec, violatorSubmatrix) );
        assert(dec->type == CMR_SEYMOUR_NODE_TYPE_IRREGULAR);

        CMR_CALL( CMRregularityTaskFree(cmr, &task) );
        queue->foundIrregularity = true;
        return CMR_OKAY;
      }

      CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
    }
  }
  else
  {
    CMR_ERROR error = CMRgraphicTestTranspose(cmr, dec->matrix, &isCographic, &dec->cograph, &dec->cographForest,
      &dec->cographCoforest, NULL, task->stats ? &task->stats->graphic : NULL, remainingTime);
    if (error == CMR_ERROR_TIMEOUT)
      return error;
    CMR_CALL( error );
  }

  CMRdbgMsg(8, "-> %s%s\n", isCographic ? "" : "NOT ", dec->isTernary ? "conetwork" : "cographic");

  dec->cographicness = isCographic ? 1 : -1;
  if (isCographic)
    dec->type = (dec->type == CMR_SEYMOUR_NODE_TYPE_GRAPH) ? CMR_SEYMOUR_NODE_TYPE_PLANAR : CMR_SEYMOUR_NODE_TYPE_COGRAPH;
  else
    queue->foundNoncographicness = true;

  if ((isCographic && (!task->params->planarityCheck || dec->graphicness)) || dec->graphicness > 0)
  {
    CMRdbgMsg(8, "Marking task as complete.\n");

    /* Task is done. */
    CMR_CALL( CMRregularityTaskFree(cmr, &task) );
  }
  else
  {
    /* Re-insert task. */
    CMRregularityQueueAdd(queue, task);
  }

  return CMR_OKAY;
}
