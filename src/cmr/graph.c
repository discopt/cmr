// #define CMR_DEBUG /* Uncomment to debug graph operations. */
// #define CMR_DEBUG_CONSISTENCY /* Uncomment to check consistency of t-decompositions. */

#include <cmr/graph.h>

#include "env_internal.h"
#include "hashtable.h"
#include "io_internal.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#define isValid(nodeOrArc) \
  ((nodeOrArc) >= 0)

void CMRgraphEnsureConsistent(CMR* cmr, CMR_GRAPH* graph)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(graph);

  CMRdbgMsg(0, "Ensuring consistency of listgraph with %d nodes and %d edges.\n", CMRgraphNumNodes(graph),
    CMRgraphNumEdges(graph));

  /* Count nodes and check prev/next linked lists. */

  size_t countNodes = 0;
#if !defined(NDEBUG)
  CMR_GRAPH_NODE u = -1;
#endif /* !NDEBUG */
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
  {
    assert(graph->nodes[v].prev == u);
#if !defined(NDEBUG)
    u = v;
#endif /* !NDEBUG */
    ++countNodes;
  }
  assert(graph->numNodes == countNodes);

  /* Count free nodes. */

  int countFree = 0;
  CMR_GRAPH_NODE v = graph->freeNode;
  while (isValid(v))
  {
    v = graph->nodes[v].next;
    ++countFree;
  }
  assert(countFree + countNodes == graph->memNodes);

  /* Check outgoing and incoming arcs. */

  size_t countIncident = 0;
  size_t countLoops = 0;
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
  {
    CMRdbgMsg(0, "First out-arc of node %d is %d\n", v, graph->nodes[v].firstOut);

    /* Check lists for outgoing arcs. */
    for (CMR_GRAPH_ITER i = CMRgraphIncFirst(graph, v); CMRgraphIncValid(graph, i); i = CMRgraphIncNext(graph, i))
    {
      CMRdbgMsg(0, "Current arc is %d = (%d,%d), opposite arc is %d, prev is %d, next is %d.\n", i,
        graph->arcs[i ^ 1].target, graph->arcs[i].target, i^1, graph->arcs[i].prev, graph->arcs[i].next);

      if (CMRgraphIncSource(graph, i) == CMRgraphIncTarget(graph, i))
        ++countLoops;
      assert(graph->arcs[i ^ 1].target == v);
      ++countIncident;
    }
  }
  assert(countIncident + countLoops == 2 * graph->numEdges);

  countFree = 0;
  int e = graph->freeEdge;

  CMRdbgMsg(0, "freeEdge = %d\n", e);

  while (isValid(e))
  {
    e = graph->arcs[2*e].next;
    ++countFree;
  }

  assert(countIncident + countLoops + 2 * countFree == 2 * graph->memEdges);

  CMRdbgMsg(0, "Consistency checked.\n");
}

CMR_ERROR CMRgraphCreateEmpty(CMR* cmr, CMR_GRAPH** pgraph, int memNodes, int memEdges)
{
  assert(cmr);
  assert(pgraph);
  assert(*pgraph == NULL);

  CMR_CALL( CMRallocBlock(cmr, pgraph) );
  CMR_GRAPH* graph = *pgraph;
  graph->numNodes = 0;
  if (memNodes <= 0)
    memNodes = 1;
  graph->memNodes = memNodes;
  graph->nodes = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &graph->nodes, memNodes) );
  graph->numEdges = 0;
  if (memEdges <= 0)
    memEdges = 1;
  graph->memEdges = memEdges;
  graph->arcs = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &graph->arcs, 2 * memEdges) );
  graph->firstNode = -1;
  graph->freeNode = (memNodes > 0) ? 0 : -1;
  for (size_t v = 0; v < graph->memNodes - 1; ++v)
    graph->nodes[v].next = v+1;
  graph->nodes[graph->memNodes-1].next = -1;
  graph->freeEdge = (memEdges > 0) ? 0 : -1;
  for (size_t e = 0; e < graph->memEdges - 1; ++e)
    graph->arcs[2*e].next = e+1;
  graph->arcs[2*graph->memEdges-2].next = -1;

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

CMR_ERROR CMRgraphFree(CMR* cmr, CMR_GRAPH** pgraph)
{
  CMR_UNUSED(cmr);

  assert(pgraph);

  CMR_GRAPH* graph = *pgraph;
  if (!graph)
    return CMR_OKAY;

  CMRdbgMsg(0, "CMRgraphFree(|V|=%d, |E|=%d)\n", CMRgraphNumNodes(graph), CMRgraphNumEdges(graph));

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  CMR_CALL( CMRfreeBlockArray(cmr, &graph->nodes) );
  CMR_CALL( CMRfreeBlockArray(cmr, &graph->arcs) );

  CMR_CALL( CMRfreeBlock(cmr, pgraph) );
  *pgraph = NULL;

  return CMR_OKAY;
}

CMR_ERROR CMRgraphClear(CMR* cmr, CMR_GRAPH* graph)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(graph);

  graph->numNodes = 0;
  graph->numEdges = 0;
  graph->firstNode = -1;
  graph->freeNode = (graph->memNodes > 0) ? 0 : -1;
  for (size_t v = 0; v < graph->memNodes - 1; ++v)
    graph->nodes[v].next = v+1;
  graph->nodes[graph->memNodes-1].next = -1;
  graph->freeEdge = (graph->memEdges > 0) ? 0 : -1;
  for (size_t e = 0; e < graph->memEdges - 1; ++e)
    graph->arcs[2*e].next = e+1;
  graph->arcs[2*graph->memEdges-2].next = -1;

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

CMR_ERROR CMRgraphAddNode(CMR* cmr, CMR_GRAPH *graph, CMR_GRAPH_NODE* pnode)
{
  assert(cmr);
  assert(graph);

  CMRdbgMsg(0, "CMRgraphAddNode().\n");

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  /* If the free list is empty, we have reallocate. */
  if (!isValid(graph->freeNode))
  {
    int mem = (graph->memNodes < 256 ? 0 : 256) + 2 * graph->memNodes;
    CMR_CALL( CMRreallocBlockArray(cmr, &graph->nodes, mem) );
    assert(graph->nodes);
    for (int v = graph->memNodes; v < mem-1; ++v)
      graph->nodes[v].next = v+1;
    graph->nodes[mem-1].next = -1;
    graph->freeNode = graph->memNodes;
    graph->memNodes = mem;
  }

  /* Add to list. */

  CMR_GRAPH_NODE node = graph->freeNode;
  graph->freeNode = graph->nodes[node].next;
  graph->numNodes++;
  graph->nodes[node].firstOut = -1;
  graph->nodes[node].prev = -1;
  graph->nodes[node].next = graph->firstNode;
  if (isValid(graph->firstNode))
    graph->nodes[graph->firstNode].prev = node;
  graph->firstNode = node;

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  if (pnode)
    *pnode = node;

  return CMR_OKAY;
}

CMR_ERROR CMRgraphAddEdge(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE u, CMR_GRAPH_NODE v,
  CMR_GRAPH_EDGE* pedge)
{
  CMRdbgMsg(0, "CMRgraphAddEdge(%d,%d).\n", u, v);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */
  assert(u >= 0);
  assert(u < (long) graph->numNodes);
  assert(v >= 0);
  assert(v < (long) graph->numNodes);

  /* If the free list is empty, we have reallocate. */

  if (!isValid(graph->freeEdge))
  {
    int newMemEdges = (graph->memEdges < 1024 ? 0 : 1024) + 2 * graph->memEdges;
    CMR_CALL( CMRreallocBlockArray(cmr, &graph->arcs, 2*newMemEdges) );
    assert(graph->arcs);
    for (int e = graph->memEdges; e < newMemEdges-1; ++e)
      graph->arcs[2*e].next = (e+1);
    graph->arcs[2*newMemEdges-2].next = -1;
    graph->freeEdge = graph->memEdges;
    graph->memEdges = newMemEdges;
  }

  /* Add to list. */

  CMR_GRAPH_EDGE edge = graph->freeEdge;
  int arc = 2*edge;
  graph->freeEdge = graph->arcs[arc].next;
  graph->numEdges++;

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

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  if (pedge)
    *pedge = edge;

  return CMR_OKAY;
}

CMR_ERROR CMRgraphDeleteNode(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE v)
{
  CMRdbgMsg(0, "CMRgraphDeleteNode(|V|=%d, |E|=%d, v=%d)\n", CMRgraphNumNodes(graph), CMRgraphNumEdges(graph), v);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  /* Remove incident edges of which v is the source. */
  while (isValid(graph->nodes[v].firstOut))
    CMRgraphDeleteEdge(cmr, graph, graph->nodes[v].firstOut/2);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  /* Remove node from node list. */
  CMR_GRAPH_NODE prev = graph->nodes[v].prev;
  CMR_GRAPH_NODE next = graph->nodes[v].next;
  if (isValid(prev))
    graph->nodes[prev].next = next;
  else
    graph->firstNode = next;
  if (isValid(next))
    graph->nodes[next].prev = prev;

  graph->nodes[v].next = graph->freeNode;
  graph->freeNode = v;
  graph->numNodes--;

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

CMR_ERROR CMRgraphDeleteEdge(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_EDGE e)
{
  CMR_UNUSED(cmr);

  CMRdbgMsg(0, "CMRgraphDeleteEdge(|V|=%d, |E|=%d, %d", CMRgraphNumNodes(graph), CMRgraphNumEdges(graph), e);

  assert(isValid(e));

  int arc = 2*e;
  CMR_GRAPH_NODE u = graph->arcs[arc+1].target;
  CMR_GRAPH_NODE v = graph->arcs[arc].target;

  CMRdbgMsg(0, " = {%d,%d})\n", u, v);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  /* Remove from u's list of outgoing arcs. */
  CMR_GRAPH_EDGE prev = graph->arcs[arc].prev;
  CMR_GRAPH_EDGE next = graph->arcs[arc].next;
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

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}


CMR_ERROR CMRgraphPrint(CMR_GRAPH* graph, FILE* stream)
{
  assert(stream);
  assert(graph);

  fprintf(stream, "Graph with %zu nodes and %zu edges.\n", CMRgraphNumNodes(graph), CMRgraphNumEdges(graph));
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
  {
    fprintf(stream, "Node %d:\n", v);
    for (CMR_GRAPH_ITER i = CMRgraphIncFirst(graph, v); CMRgraphIncValid(graph, i); i = CMRgraphIncNext(graph, i))
    {
      fprintf(stream, "  Edge %d: {%d,%d} {arc = %d}\n", CMRgraphIncEdge(graph, i), CMRgraphIncSource(graph, i),
        CMRgraphIncTarget(graph, i), i);
    }
  }

  return CMR_OKAY;
}

CMR_ERROR CMRgraphMergeNodes(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH_NODE u, CMR_GRAPH_NODE v)
{
  CMR_UNUSED(cmr);

  assert(graph);
  assert(u >= 0);
  assert(u < (long) graph->memNodes);
  assert(v >= 0);
  assert(v < (long) graph->memNodes);
  assert(u != v);

  int a;
  while ((a = graph->nodes[v].firstOut) >= 0)
  {
    graph->nodes[v].firstOut = graph->arcs[a].next;
    graph->arcs[a ^ 1].target = u;
    graph->arcs[a].next = graph->nodes[u].firstOut;
    if (graph->nodes[u].firstOut >= 0)
      graph->arcs[graph->nodes[u].firstOut].prev = a;
    graph->arcs[a].prev = -1;
    graph->nodes[u].firstOut = a;
  }

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRgraphEnsureConsistent(cmr, graph);
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

CMR_ERROR CMRgraphCreateFromEdgeList(CMR* cmr, CMR_GRAPH** pgraph, CMR_ELEMENT** pedgeElements, char*** pnodeLabels,
  FILE* stream)
{
  assert(cmr);
  assert(pgraph);
  assert(!*pgraph);
  assert(!pedgeElements || !*pedgeElements);
  assert(!pnodeLabels || !*pnodeLabels);
  assert(stream);

  char* line = NULL;
  size_t length = 0;
  ssize_t numRead;

  char* uToken = NULL;
  char* vToken = NULL;
  char* elementToken = NULL;

  CMR_CALL( CMRgraphCreateEmpty(cmr, pgraph, 256, 1024) );
  CMR_GRAPH* graph = *pgraph;
  size_t memNodeLabels = 256;
  if (pnodeLabels)
    CMR_CALL( CMRallocBlockArray(cmr, pnodeLabels, memNodeLabels) );
  size_t memEdgeElements = 256;
  if (pedgeElements)
    CMR_CALL( CMRallocBlockArray(cmr, pedgeElements, memEdgeElements) );

  CMR_LINEARHASHTABLE_ARRAY* nodeNames = NULL;
  CMR_CALL( CMRlinearhashtableArrayCreate(cmr, &nodeNames, 8, 1024) );
  
  while ((numRead = getline(&line, &length, stream)) != -1)
  {
    char* s = line;

    /* Scan for whitespace */
    while (*s && isspace(*s))
      ++s;
    if (!*s)
      break;

    /* Scan name of node u. */
    uToken = s;
    while (*s && !isspace(*s))
      ++s;
    *s = '\0';
    ++s;

    /* Scan for whitespace */
    while (*s && isspace(*s))
      ++s;
    if (!*s)
      break;

    /* Scan name of node v. */
    vToken = s;
    while (*s && !isspace(*s))
      ++s;
    *s = '\0';
    ++s;

    /* Scan for whitespace */
    while (*s && isspace(*s))
      ++s;
    if (*s)
    {
      /* Scan element. */
      elementToken = s;
      while (*s && !isspace(*s))
        ++s;
      *s = '\0';
      ++s;
    }
    else
    {
      elementToken = NULL;
    }

    /* Figure out node u if it exists. */

    CMR_LINEARHASHTABLE_BUCKET bucket;
    CMR_LINEARHASHTABLE_HASH hash;
    CMR_GRAPH_NODE uNode; 
    if (CMRlinearhashtableArrayFind(nodeNames, uToken, strlen(uToken), &bucket, &hash))
      uNode = (CMR_GRAPH_NODE) (size_t) CMRlinearhashtableArrayValue(nodeNames, bucket);
    else
    {
      CMR_CALL( CMRgraphAddNode(cmr, graph, &uNode) );
      CMR_CALL( CMRlinearhashtableArrayInsertBucketHash(cmr, nodeNames, uToken, strlen(uToken), bucket, hash,
        (void*) (size_t) uNode) );

      /* Add node label. */
      if (pnodeLabels)
      {
        if (uNode >= (int)memNodeLabels)
        {
          memNodeLabels *= 2;
          CMR_CALL( CMRreallocBlockArray(cmr, pnodeLabels, memNodeLabels) );
        }

        (*pnodeLabels)[uNode] = strdup(uToken);
      }
    }

    /* Figure out node v if it exists. */

    CMR_GRAPH_NODE vNode; 
    if (CMRlinearhashtableArrayFind(nodeNames, vToken, strlen(vToken), &bucket, &hash))
      vNode = (CMR_GRAPH_NODE) (size_t) CMRlinearhashtableArrayValue(nodeNames, bucket);
    else
    {
      CMR_CALL( CMRgraphAddNode(cmr, graph, &vNode) );
      CMR_CALL( CMRlinearhashtableArrayInsertBucketHash(cmr, nodeNames, vToken, strlen(vToken), bucket, hash,
        (void*) (size_t) vNode) );
      
      /* Add node label. */
      if (pnodeLabels)
      {
        if (vNode >= (int)memNodeLabels)
        {
          memNodeLabels *= 2;
          CMR_CALL( CMRreallocBlockArray(cmr, pnodeLabels, memNodeLabels) );
        }

        (*pnodeLabels)[vNode] = strdup(vToken);
      }
    }

    /* Extract element. */

    CMR_ELEMENT element = 0;
    if (elementToken)
    {
      if (strchr("rRtT-", elementToken[0]))
      {
        sscanf(&elementToken[1], "%d", &element);
        element = -element;
      }
      else if (strchr("cC", elementToken[0]))
        sscanf(&elementToken[1], "%d", &element);
      else
        sscanf(elementToken, "%d", &element);
    }

    CMR_GRAPH_EDGE edge;
    CMR_CALL( CMRgraphAddEdge(cmr, graph, uNode, vNode, &edge) );

    if (pedgeElements)
    {
      if (edge >= (int)memEdgeElements)
      {
        do
        {
          memEdgeElements *= 2;
        }
        while (edge >= (int)memEdgeElements);
        CMR_CALL( CMRreallocBlockArray(cmr, pedgeElements, memEdgeElements) );
      }

      (*pedgeElements)[edge] = element;
    }
  }

  CMR_CALL( CMRlinearhashtableArrayFree(cmr, &nodeNames) );
  free(line);

  return CMR_OKAY;
}

CMR_ERROR CMRgraphWriteEdgeList(CMR* cmr, CMR_GRAPH* graph, CMR_ELEMENT* edgeElements, const char** nodeLabels, FILE* stream)
{
  assert(cmr);
  assert(graph);

  for (CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, iter);
    iter = CMRgraphEdgesNext(graph, iter))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, iter);
    CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
    CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
    if (nodeLabels)
      fprintf(stream, "%s", nodeLabels[u]);
    else
      fprintf(stream, "v%d", u);
    if (nodeLabels)
      fprintf(stream, " %s", nodeLabels[v]);
    else
      fprintf(stream, " v%d", v);
    if (edgeElements)
      fprintf(stream, " %s", CMRelementString(edgeElements[e], 0));
    else
      fprintf(stream, " e%d", e);
    fprintf(stream, "\n");
  }

  return CMR_OKAY;
}


CMR_ERROR CMRgraphCopy(CMR* cmr, CMR_GRAPH* graph, CMR_GRAPH** pcopy)
{
  assert(cmr);
  assert(graph);
  assert(pcopy);

  size_t n = CMRgraphMemNodes(graph);
  size_t m = CMRgraphMemEdges(graph);
  CMR_CALL( CMRgraphCreateEmpty(cmr, pcopy, n, m) );
  CMR_GRAPH* copy = *pcopy;
  assert(copy);

  copy->firstNode = graph->firstNode;
  copy->freeEdge = graph->freeEdge;
  copy->freeNode = graph->freeNode;
  for (size_t v = 0; v < n; ++v)
    copy->nodes[v] = graph->nodes[v];
  for (size_t a = 0; a < 2 * m; ++m)
    copy->arcs[a] = graph->arcs[a];

  return CMR_OKAY;
}

