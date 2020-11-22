// #define DEBUG_GRAPH

#include <tu/graph.h>

#include "env_internal.h"

#include <assert.h>

#define isValid(nodeOrArc) \
  ((nodeOrArc) >= 0)

void TUlistgraphEnsureConsistent(TU* tu, TU_LISTGRAPH* graph)
{
  assert(tu);
  assert(graph);

#ifdef DEBUG_GRAPH
  printf("Ensuring consistency of listgraph with %d nodes and %d edges.\n",
    TUlistgraphNumNodes(graph), TUlistgraphNumEdges(graph));
#endif /* DEBUG_GRAPH */

  /* Count nodes and check prev/next linked lists. */

  int countNodes = 0;
  TU_LISTGRAPH_NODE u = -1;
  for (TU_LISTGRAPH_NODE v = TUlistgraphNodesFirst(graph); TUlistgraphNodesValid(graph, v);
    TUlistgraphNodesNext(graph, v))
  {
    assert(graph->nodes[v].prev == u);
    u = v;
    ++countNodes;
  }
  assert(graph->numNodes == countNodes);

  /* Count free nodes. */

  int countFree = 0;
  TU_LISTGRAPH_NODE v = graph->freeNode;
  while (isValid(v))
  {
    v = graph->nodes[v].next;
    ++countFree;
  }
  assert(countFree + countNodes == graph->memNodes);

  /* Check outgoing and incoming arcs. */

  int countIncident = 0;
  for (TU_LISTGRAPH_NODE v = TUlistgraphNodesFirst(graph); TUlistgraphNodesValid(graph, v);
    TUlistgraphNodesNext(graph, v))
  {
#ifdef DEBUG_GRAPH
    printf("First out-arc of node %d is %d\n", v, graph->nodes[v].firstOut);
#endif /* DEBUG_GRAPH */

    /* Check lists for outgoing arcs. */
    int j = -1;
    for (TU_LISTGRAPH_INCIDENT i = TUlistgraphIncidentFirst(graph, v);
      TUlistgraphIncidentValid(graph, i); TUlistgraphIncidentNext(graph, i))
    {
#ifdef DEBUG_GRAPH
      printf("Current arc is %d = (%d,%d), other is %d, prev is %d, next is %d.\n", i,
        graph->arcs[i ^ 1].target, graph->arcs[i].target, i^1, graph->arcs[i].prev,
        graph->arcs[i].next);
#endif /* DEBUG_GRAPH */

      assert(graph->arcs[i ^ 1].target == v);
      assert(graph->arcs[i].prev == j);
      j = i;
      ++countIncident;
    }
  }
  assert(countIncident == 2 * graph->numEdges);

  countFree = 0;
  int e = graph->freeEdge;

#ifdef DEBUG_GRAPH
  printf("freeEdge = %d\n", e);
#endif /* DEBUG_GRAPH */

  while (isValid(e))
  {
    e = graph->arcs[2*e].next;
    ++countFree;
  }

  assert(countIncident + 2 * countFree == 2 * graph->memEdges);
}

void TUlistgraphCreateEmpty(TU* tu, TU_LISTGRAPH** pgraph, int memNodes, int memEdges)
{
  assert(tu);
  assert(pgraph);
  assert(*pgraph == NULL);

  TUallocBlock(tu, pgraph);
  TU_LISTGRAPH* graph = *pgraph;
  graph->numNodes = 0;
  graph->memNodes = memNodes;
  graph->nodes = NULL;
  TUallocBlockArray(tu, &graph->nodes, memNodes);
  graph->numEdges = 0;
  graph->memEdges = memEdges;
  graph->arcs = NULL;
  TUallocBlockArray(tu, &graph->arcs, 2 * memEdges);
  graph->firstNode = -1;
  graph->freeNode = (memNodes > 0) ? 0 : -1;
  for (int v = 0; v < graph->memNodes - 1; ++v)
    graph->nodes[v].next = v+1;
  graph->nodes[graph->memNodes-1].next = -1;
  graph->freeEdge = (memEdges > 0) ? 0 : -1;
  for (int e = 0; e < graph->memEdges - 1; ++e)
    graph->arcs[2*e].next = e+1;
  graph->arcs[2*graph->memEdges-2].next = -1;

  TUlistgraphEnsureConsistent(tu, graph);
}

void TUlistgraphFree(TU* tu, TU_LISTGRAPH** pgraph)
{
  assert(pgraph);

#ifdef DEBUG_GRAPH
  printf("TUlistgraphFree(|V|=%d, |E|=%d)\n", TUlistgraphNumNodes(*pgraph),
    TUlistgraphNumEdges(*pgraph));
#endif /* DEBUG_GRAPH */

  TUlistgraphEnsureConsistent(tu, *pgraph);

  TU_LISTGRAPH* graph = *pgraph;

  TUfreeBlockArray(tu, &graph->nodes);
  TUfreeBlockArray(tu, &graph->arcs);

  TUfreeBlock(tu, pgraph);
  *pgraph = NULL;
}

TU_LISTGRAPH_NODE TUlistgraphAddNode(TU* tu, TU_LISTGRAPH *graph)
{
#ifdef DEBUG_GRAPH
  printf("TUlistgraphAddNode\n");
#endif /* DEBUG_GRAPH */

  TUlistgraphEnsureConsistent(tu, graph);

  /* If the free list is empty, we have reallocate. */
  if (!isValid(graph->freeNode))
  {
    int mem = (graph->memNodes < 256 ? 0 : 256) + 2 * graph->memNodes;
    TUreallocBlockArray(tu, &graph->nodes, mem);
    assert(graph->nodes);
    for (int v = graph->memNodes; v < mem-1; ++v)
      graph->nodes[v].next = v+1;
    graph->nodes[mem-1].next = -1;
    graph->freeNode = graph->memNodes;
    graph->memNodes = mem;
  }

  /* Add to list. */

  int result = graph->freeNode;
  graph->freeNode = graph->nodes[result].next;
  graph->numNodes++;
  graph->nodes[result].firstOut = -1;
  graph->nodes[result].prev = -1;
  graph->nodes[result].next = graph->firstNode;
  if (isValid(graph->firstNode))
    graph->nodes[graph->firstNode].prev = result;
  graph->firstNode = result;

  TUlistgraphEnsureConsistent(tu, graph);

  return result;
}

TU_LISTGRAPH_EDGE TUlistgraphAddEdge(TU* tu, TU_LISTGRAPH* graph, TU_LISTGRAPH_NODE u,
  TU_LISTGRAPH_NODE v)
{
#ifdef DEBUG_GRAPH
  printf("TUlistgraphAddEdge\n");
#endif /* DEBUG_GRAPH */

  TUlistgraphEnsureConsistent(tu, graph);
  assert(u >= 0);
  assert(u < graph->numNodes);
  assert(v >= 0);
  assert(v < graph->numNodes);

  /* If the free list is empty, we have reallocate. */

  if (!isValid(graph->freeEdge))
  {
    int newMemEdges = (graph->memEdges < 1024 ? 0 : 1024) + 2 * graph->memEdges;
    TUreallocBlockArray(tu, &graph->arcs, 2*newMemEdges);
    assert(graph->arcs);
    for (int e = graph->memEdges; e < newMemEdges-1; ++e)
      graph->arcs[2*e].next = (e+1);
    graph->arcs[2*newMemEdges-2].next = -1;
    graph->freeEdge = graph->memEdges;
    graph->memEdges = newMemEdges;
  }
  
  /* Add to list. */

  int result = graph->freeEdge;
  graph->freeEdge = graph->arcs[2*result].next;
  graph->numEdges++;
  int arc = 2*result;

  /* Add arc to list of outgoing arcs from u. */

  graph->arcs[arc].target = v;
  int firstOut = graph->nodes[u].firstOut;
  graph->arcs[arc].prev = -1;
  graph->arcs[arc].next = firstOut;
  if (isValid(firstOut))
    graph->arcs[firstOut].prev = arc;
  graph->nodes[u].firstOut = arc;

  ++arc;

  /* Add arc to list of incoming arcs to v. */
  graph->arcs[arc].target = u;
  firstOut = graph->nodes[v].firstOut;
  graph->arcs[arc].prev = -1;
  graph->arcs[arc].next = firstOut;
  if (isValid(firstOut))
    graph->arcs[firstOut].prev = arc;
  graph->nodes[v].firstOut = arc;

  TUlistgraphEnsureConsistent(tu, graph);

  return result;
}

void TUlistgraphDeleteNode(TU* tu, TU_LISTGRAPH* graph, TU_LISTGRAPH_NODE v)
{
#ifdef DEBUG_GRAPH
  printf("TUlistgraphDeleteNode(|V|=%d, |E|=%d, v=%d)\n", TUlistgraphNumNodes(graph),
    TUlistgraphNumEdges(graph), v);
#endif /* DEBUG_GRAPH */

  TUlistgraphEnsureConsistent(tu, graph);

  /* Remove incident edges of which v is the source. */
  while (isValid(graph->nodes[v].firstOut))
    TUlistgraphDeleteEdge(tu, graph, graph->nodes[v].firstOut/2);

  TUlistgraphEnsureConsistent(tu, graph);

  /* Remove node from node list. */
  TU_LISTGRAPH_NODE prev = graph->nodes[v].prev;
  TU_LISTGRAPH_NODE next = graph->nodes[v].next;
  if (isValid(prev))
    graph->nodes[prev].next = next;
  else
    graph->firstNode = next;
  if (isValid(next))
    graph->nodes[next].prev = prev;

  graph->nodes[v].next = graph->freeNode;
  graph->freeNode = v;
  graph->numNodes--;

  TUlistgraphEnsureConsistent(tu, graph);
}

void TUlistgraphDeleteEdge(TU* tu, TU_LISTGRAPH* graph, TU_LISTGRAPH_EDGE e)
{
#ifdef DEBUG_GRAPH
  printf("TUlistgraphDeleteEdge(|V|=%d, |E|=%d, %d", TUlistgraphNumNodes(graph),
    TUlistgraphNumEdges(graph), e);
  fflush(stdout);
#endif /* DEBUG_GRAPH */

  assert(isValid(e));

  int arc = 2*e;
  TU_LISTGRAPH_NODE u = graph->arcs[arc+1].target;
  TU_LISTGRAPH_NODE v = graph->arcs[arc].target;

#ifdef DEBUG_GRAPH
  printf(" = {%d,%d})\n", u, v);
#endif /* DEBUG_GRAPH */

  TUlistgraphEnsureConsistent(tu, graph);

  /* Remove from u's list of outgoing arcs. */  
  TU_LISTGRAPH_EDGE prev = graph->arcs[arc].prev;
  TU_LISTGRAPH_EDGE next = graph->arcs[arc].next;
  if (isValid(prev))
    graph->arcs[prev].next = next;
  else
    graph->nodes[u].firstOut = next;
  if (isValid(next))
    graph->arcs[next].prev = prev;

  /* Add to free list. */
  graph->arcs[arc].next = graph->freeEdge;
  graph->freeEdge = e;
  graph->numEdges--;

  /* Remove from v's list of outgoing arcs. */
  ++arc;
  prev = graph->arcs[arc].prev;
  next = graph->arcs[arc].next;
  if (isValid(prev))
    graph->arcs[prev].next = next;
  else
    graph->nodes[v].firstOut = next;
  if (isValid(next))
    graph->arcs[next].prev = prev;

  TUlistgraphEnsureConsistent(tu, graph);
}


void TUlistgraphPrint(FILE* stream, TU_LISTGRAPH* graph)
{
  for (TU_LISTGRAPH_NODE v = TUlistgraphNodesFirst(graph); TUlistgraphNodesValid(graph, v);
    TUlistgraphNodesNext(graph, v))
  {
    fprintf(stream, "Node %d:\n", v);
    for (TU_LISTGRAPH_EDGE e = TUlistgraphIncidentFirst(graph, v);
      TUlistgraphIncidentValid(graph, e); TUlistgraphIncidentNext(graph, e))
    {
      fprintf(stream, "  Edge %d: {%d,%d}\n", e, TUlistgraphIncidentSource(graph, e),
        TUlistgraphIncidentTarget(graph, e));
    }
  }
}
