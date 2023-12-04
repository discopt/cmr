// #define CMR_DEBUG /* Uncomment to debug graphic. */
// #define CMR_DEBUG_DOT /* Uncomment to write dot files of t-decompositions. */
// #define CMR_DEBUG_CONSISTENCY /* Uncomment to check consistency of t-decompositions. */

#include <cmr/graphic.h>

#include "env_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "heap.h"
#include "sort.h"
#include "hereditary_property.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define SWAP_INTS(a, b) \
  do \
  { \
    int tmp = a; \
    a = b; \
    b = tmp; \
  } \
  while (false)

CMR_ERROR CMRstatsGraphicInit(CMR_GRAPHIC_STATISTICS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  stats->checkCount = 0;
  stats->checkTime = 0.0;
  stats->applyCount = 0;
  stats->applyTime = 0.0;
  stats->transposeCount = 0;
  stats->transposeTime = 0.0;

  return CMR_OKAY;
}

CMR_ERROR CMRstatsGraphicPrint(FILE* stream, CMR_GRAPHIC_STATISTICS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Graphicness recognition:\n");
    prefix = "  ";
  }
  fprintf(stream, "%stranspositions: %lu in %f seconds\n", prefix, (unsigned long)stats->transposeCount,
    stats->transposeTime);
  fprintf(stream, "%scolumn checks: %lu in %f seconds\n", prefix, (unsigned long)stats->checkCount, stats->checkTime);
  fprintf(stream, "%scolumn additions: %lu in %f seconds\n", prefix, (unsigned long)stats->applyCount,
    stats->applyTime);
  fprintf(stream, "%stotal: %lu in %f seconds\n", prefix, (unsigned long)stats->totalCount, stats->totalTime);

  return CMR_OKAY;
}

typedef enum
{
  UNKNOWN = 0,    /**< \brief The node was not considered by the shortest-path, yet. */
  SEEN = 1,       /**< \brief Some path to the node is known. */
  COMPLETED = 2,  /**< \brief The shortest path to the node is known. */
  BASIC = 3,      /**< \brief The rootEdge of that node belongs to the spanning forest. */
} DIJKSTRA_STAGE;

/**
 * \brief Node information for shortest-path computation in \ref CMRcomputeRepresentationMatrix().
 */

typedef struct
{
  DIJKSTRA_STAGE stage;     /**< \brief At which stage of the algorithm is this node? */
  int predecessor;          /**< \brief Predecessor node in shortest-path branching, or -1 for a root. */
  CMR_GRAPH_EDGE rootEdge;  /**< \brief The actual edge towards the predecessor, or -1 for a root./ */
  bool reversed;            /**< \brief Whether the edge towards the predecessor is reversed. */
} DijkstraNodeData;

CMR_ERROR CMRcomputeRepresentationMatrix(CMR* cmr, CMR_GRAPH* digraph, bool ternary, CMR_CHRMAT** ptranspose,
  bool* arcsReversed, int numForestArcs, CMR_GRAPH_EDGE* forestArcs, int numCoforestArcs, CMR_GRAPH_EDGE* coforestArcs,
  bool* pisCorrectForest)
{
  assert(cmr);
  assert(digraph);
  assert(ptranspose && !*ptranspose);
  assert(numForestArcs == 0 || forestArcs);
  assert(numForestArcs == 0 || coforestArcs);
  CMRassertStackConsistency(cmr);

  CMRdbgMsg(0, "Computing %s representation matrix.\n", ternary ? "ternary" : "binary");

  DijkstraNodeData* nodeData = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &nodeData, CMRgraphMemNodes(digraph)) );
  CMR_INTHEAP heap;
  CMR_CALL( CMRintheapInitStack(cmr, &heap, CMRgraphMemNodes(digraph)) );
  int* lengths = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &lengths, CMRgraphMemEdges(digraph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(digraph); CMRgraphNodesValid(digraph, v);
    v = CMRgraphNodesNext(digraph, v))
  {
    nodeData[v].stage = UNKNOWN;
  }
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(digraph); CMRgraphEdgesValid(digraph, i);
    i = CMRgraphEdgesNext(digraph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(digraph, i);
    lengths[e] = 1;
  }
  for (int b = 0; b < numForestArcs; ++b)
  {
    if (forestArcs[b] >= 0)
    {
      CMRdbgMsg(0, "forest element %d is edge %d = {%d,%d}\n", b, forestArcs[b], CMRgraphEdgeU(digraph, forestArcs[b]),
        CMRgraphEdgeV(digraph, forestArcs[b]));
      lengths[forestArcs[b]] = 0;
    }
  }

  CMRassertStackConsistency(cmr);

  /* Start Dijkstra's algorithm at each node. */
  int countComponents = 0;
  for (CMR_GRAPH_NODE s = CMRgraphNodesFirst(digraph); CMRgraphNodesValid(digraph, s);
    s = CMRgraphNodesNext(digraph, s))
  {
    if (nodeData[s].stage != UNKNOWN)
      continue;

    CMRdbgMsg(2, "Executing Dijkstra at starting node %d.\n", s);
    nodeData[s].predecessor = -1;
    nodeData[s].rootEdge = -1;
    ++countComponents;
    CMRintheapInsert(&heap, s, 0);
    while (!CMRintheapEmpty(&heap))
    {
      int distance = CMRintheapMinimumValue(&heap);
      CMR_GRAPH_NODE v = CMRintheapExtractMinimum(&heap);
      CMRdbgMsg(4, "Processing node %d at distance %d.\n", v, distance);
      nodeData[v].stage = COMPLETED;
      for (CMR_GRAPH_ITER i = CMRgraphIncFirst(digraph, v); CMRgraphIncValid(digraph, i);
        i = CMRgraphIncNext(digraph, i))
      {
        assert(CMRgraphIncSource(digraph, i) == v);
        CMR_GRAPH_NODE w = CMRgraphIncTarget(digraph, i);

        /* Skip if already completed. */
        if (nodeData[w].stage == COMPLETED)
          continue;

        CMR_GRAPH_EDGE e = CMRgraphIncEdge(digraph, i);
        int newDistance = distance + lengths[e];
        if (newDistance < CMRintheapGetValueInfinity(&heap, w))
        {
          CMRdbgMsg(6, "Updating distance of (%d,%d) from %d to %d.\n", v, w, CMRintheapGetValueInfinity(&heap, w),
            newDistance);
          nodeData[w].stage = SEEN;
          nodeData[w].predecessor = v;
          nodeData[w].rootEdge = e;
          nodeData[w].reversed = arcsReversed ? arcsReversed[e] : false;
          if (w == CMRgraphEdgeU(digraph, e))
            nodeData[w].reversed = !nodeData[w].reversed;
          CMRintheapDecreaseInsert(&heap, w, newDistance);
        }
      }
    }
  }

  CMRassertStackConsistency(cmr);
  CMR_CALL( CMRfreeStackArray(cmr, &lengths) );
  CMR_CALL( CMRintheapClearStack(cmr, &heap) );

  /* Now nodeData[.].predecessor is an arborescence for each connected component. */

  CMR_GRAPH_NODE* nodesRows = NULL; /* Non-root node v is mapped to row of edge {v,predecessor(v)}. */
  CMR_CALL( CMRallocStackArray(cmr, &nodesRows, CMRgraphMemNodes(digraph)) );
  char* nodesReversed = NULL; /* Non-root node v is mapped to +1 or -1 depending on the direction of {v,predecessor(v)}. */
  CMR_CALL( CMRallocStackArray(cmr, &nodesReversed, CMRgraphMemNodes(digraph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(digraph); CMRgraphNodesValid(digraph, v); v = CMRgraphNodesNext(digraph, v))
  {
    nodesRows[v] = -1;
    nodesReversed[v] = 1;
  }

  size_t numRows = 0;
  if (pisCorrectForest)
    *pisCorrectForest = true;
  for (int i = 0; i < numForestArcs; ++i)
  {
    CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, forestArcs[i]);
    CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, forestArcs[i]);
    CMRdbgMsg(2, "Forest edge %d = {%d,%d}.\n", forestArcs[i], u, v);
    if (nodeData[u].predecessor == v)
    {
      nodesRows[u] = numRows;
      nodesReversed[u] = nodeData[u].reversed ? -1 : 1;
      ++numRows;
      nodeData[u].stage = BASIC;
      CMRdbgMsg(2, "Basic edge (%d,%d): %d is predecessor of %d; node %d is row %d; reversed = %s.\n", v, u, v, u, u,
        nodesRows[u], nodeData[u].reversed ? "true" : "false");
    }
    else if (nodeData[v].predecessor == u)
    {
      nodesRows[v] = numRows;
      nodesReversed[v] = nodeData[v].reversed ? -1 : 1;
      ++numRows;
      nodeData[v].stage = BASIC;
      CMRdbgMsg(2, "Basic edge (%d,%d): %d is predecessor of %d; node %d is row %d; reversed = %s\n", u, v, u, v, v,
        nodesRows[v], nodeData[v].reversed ? "true" : "false");
    }
    else
    {
      /* A provided basic edge is not part of the spanning forest. */
      if (pisCorrectForest)
        *pisCorrectForest = false;
    }
  }
  if (numRows + countComponents < CMRgraphNumNodes(digraph))
  {
    /* Some edge from the spanning forest is not a forest edge. */
    if (pisCorrectForest)
      *pisCorrectForest = false;

    for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(digraph); CMRgraphNodesValid(digraph, v); v = CMRgraphNodesNext(digraph, v))
    {
      if (nodeData[v].predecessor >= 0 && nodeData[v].stage != BASIC)
      {
        nodesRows[v] = numRows;
        nodesReversed[v] = nodeData[v].reversed ? -1 : 1;
        ++numRows;
        nodeData[v].stage = BASIC;
        CMRdbgMsg(2, "Predecessor edge {%d,%d} not basic; node %d is row %d.\n", nodeData[v].predecessor, v, v,
          nodesRows[v]);
      }
    }
  }

  CMRassertStackConsistency(cmr);

  CMR_CALL( CMRchrmatCreate(cmr, ptranspose, CMRgraphNumEdges(digraph) - numRows, numRows,
    16 * numRows) );
  CMR_CHRMAT* transpose = *ptranspose;
  size_t numNonzeros = 0; /* Current number of nonzeros. transpose->numNonzeros is the memory. */
  size_t numColumns = 0;
  CMR_GRAPH_EDGE* edgeColumns = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &edgeColumns, CMRgraphMemEdges(digraph)) );
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(digraph); CMRgraphEdgesValid(digraph, i);
    i = CMRgraphEdgesNext(digraph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(digraph, i);
    CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
    CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
    edgeColumns[CMRgraphEdgesEdge(digraph, i)] =
      (nodeData[u].rootEdge == e || nodeData[v].rootEdge == e) ? -1 : -2;
  }
  CMR_GRAPH_NODE* uPath = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &uPath, numRows) );
  CMR_GRAPH_NODE* vPath = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &vPath, numRows) );
  CMR_GRAPH_ITER iter = CMRgraphEdgesFirst(digraph);
  int cobasicIndex = 0;
  while (CMRgraphEdgesValid(digraph, iter))
  {
    CMR_GRAPH_EDGE e = -1;
    while (cobasicIndex < numCoforestArcs && e < 0)
    {
      e = coforestArcs[cobasicIndex];
      ++cobasicIndex;
    }
    if (cobasicIndex >= numCoforestArcs)
    {
      e = CMRgraphEdgesEdge(digraph, iter);
      iter = CMRgraphEdgesNext(digraph, iter);
    }

    if (edgeColumns[e] >= -1)
      continue;

    CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
    CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
    if (arcsReversed && arcsReversed[e])
      SWAP_INTS(u, v);
    CMRdbgMsg(4, "Edge %d = (%d,%d) (including reverse)\n", e, u, v);

    transpose->rowSlice[numColumns] = numNonzeros;
    edgeColumns[e] = numColumns;

    /* Enlarge space for nonzeros if necessary. */
    if (numNonzeros + numRows > transpose->numNonzeros)
      CMR_CALL( CMRchrmatChangeNumNonzeros(cmr, transpose, 2 * transpose->numNonzeros) );

    /* Compute u-root path. */
    int uPathLength = 0;
    CMR_GRAPH_NODE w = u;
    while (nodeData[w].predecessor != -1)
    {
      uPath[uPathLength] = w;
      ++uPathLength;
      w = nodeData[w].predecessor;
    }

    /* Compute v-root path. */
    int vPathLength = 0;
    w = v;
    while (nodeData[w].predecessor != -1)
    {
      vPath[vPathLength] = w;
      ++vPathLength;
      w = nodeData[w].predecessor;
    }

    /* Remove common part of u-root path and v-root path. */
    while (uPathLength > 0 && vPathLength > 0 && uPath[uPathLength-1] == vPath[vPathLength-1])
    {
      --uPathLength;
      --vPathLength;
    }

    /* Create nonzeros. */
    for (int j = 0; j < uPathLength; ++j)
    {
      assert(nodesRows[uPath[j]] >= 0);
      transpose->entryColumns[numNonzeros] = nodesRows[uPath[j]];
      transpose->entryValues[numNonzeros] = ternary ? -nodesReversed[uPath[j]] : 1;
      CMRdbgMsg(4, "u: row %d. nodesReversed = %d, result = %d\n", nodesRows[uPath[j]],
        nodesReversed[uPath[j]], transpose->entryValues[numNonzeros]);
      ++numNonzeros;
    }
    for (int j = 0; j < vPathLength; ++j)
    {
      assert(nodesRows[vPath[j]] >= 0);
      transpose->entryColumns[numNonzeros] = nodesRows[vPath[j]];
      transpose->entryValues[numNonzeros] = ternary ? nodesReversed[vPath[j]] : 1;
      CMRdbgMsg(4, "v: row %d. nodesReversed = %d, result = %d\n", nodesRows[vPath[j]],
        nodesReversed[vPath[j]], transpose->entryValues[numNonzeros]);
      ++numNonzeros;
    }

    ++numColumns;
  }
  transpose->rowSlice[numColumns] = numNonzeros;
  transpose->numNonzeros = numNonzeros;

  CMR_CALL( CMRchrmatSortNonzeros(cmr, transpose) );

  CMRassertStackConsistency(cmr);

  CMR_CALL( CMRfreeStackArray(cmr, &vPath) );
  CMR_CALL( CMRfreeStackArray(cmr, &uPath) );
  CMR_CALL( CMRfreeStackArray(cmr, &edgeColumns) );

  if (numNonzeros == 0 && transpose->numNonzeros > 0)
  {
    CMR_CALL( CMRfreeBlockArray(cmr, &transpose->entryColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &transpose->entryValues) );
  }
  transpose->numNonzeros = numNonzeros;

  /* We now process the nonbasic edges. */

  CMRassertStackConsistency(cmr);
  CMR_CALL( CMRfreeStackArray(cmr, &nodesReversed) );
  CMRassertStackConsistency(cmr);
  CMR_CALL( CMRfreeStackArray(cmr, &nodesRows) );
  CMRassertStackConsistency(cmr);
  CMR_CALL( CMRfreeStackArray(cmr, &nodeData) );

  CMRassertStackConsistency(cmr);

  return CMR_OKAY;
}

CMR_ERROR CMRcomputeGraphicMatrix(CMR* cmr, CMR_GRAPH* graph, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose,
  int numForestEdges, CMR_GRAPH_EDGE* forestEdges, int numCoforestEdges, CMR_GRAPH_EDGE* coforestEdges,
  bool* pisCorrectForest)
{
  assert(cmr);
  assert(graph);
  assert(pmatrix || ptranspose);
  assert(!pmatrix || !*pmatrix);
  assert(!ptranspose || !*ptranspose);

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRcomputeRepresentationMatrix(cmr, graph, false, &transpose, NULL, numForestEdges, forestEdges,
    numCoforestEdges, coforestEdges, pisCorrectForest) );

  if (pmatrix)
    CMR_CALL( CMRchrmatTranspose(cmr, transpose, pmatrix) );

  /* Return or free the transpose matrix. */
  if (ptranspose)
    *ptranspose = transpose;
  else
    CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}

typedef enum
{
  DEC_MEMBER_TYPE_INVALID = 0,
  DEC_MEMBER_TYPE_PARALLEL = 1, /**< \brief Indicate a parallel member, also known as a bond. */
  DEC_MEMBER_TYPE_SERIES = 2,   /**< \brief Indicate a series member, also known as a polygon. */
  DEC_MEMBER_TYPE_RIGID = 3,    /**< \brief Indicate a parallel member, also known as prime. */
  DEC_MEMBER_TYPE_LOOP = 4      /**< \brief Indicate a loop. */
} DEC_MEMBER_TYPE;

static inline
const char* memberTypeString(
  DEC_MEMBER_TYPE type
)
{
  switch(type)
  {
  case DEC_MEMBER_TYPE_PARALLEL:
    return "parallel";
  case DEC_MEMBER_TYPE_SERIES:
    return "series";
  case DEC_MEMBER_TYPE_RIGID:
    return "rigid";
  case DEC_MEMBER_TYPE_LOOP:
    return "loop";
  default:
    return "invalid";
  }
}

typedef size_t DEC_EDGE;   /**< \brief Type for refering to a decomposition edge. */
typedef size_t DEC_NODE;   /**< \brief Type for refering to a decomposition node. */
typedef size_t DEC_MEMBER; /**< \brief Type for refering to a decomposition member. */

typedef struct
{
  DEC_NODE representativeNode;  /**< \brief Next representative of same node towards root, or -1 if root. */
} DecNodeData;

typedef struct
{
  CMR_ELEMENT element;    /**< \brief Element corresponding to this edge.
                           *
                           * 1, 2, ..., m indicate rows, -1,-2, ..., -n indicate columns,
                           * and for (small) k >= 0, MAX_INT-k and -MAX_INT+k indicate
                           * markers of the parent and to the parent, respectively. */
  DEC_MEMBER member;      /**< \brief Member this edge belongs to or -1 if in free list. */
  DEC_NODE head;          /**< \brief Head node of this edge for a rigid member, -1 otherwise. */
  DEC_NODE tail;          /**< \brief Tail node of this edge for a rigid member, -1 otherwise. */
  DEC_EDGE prev;          /**< \brief Next edge of this member. */
  DEC_EDGE next;          /**< \brief Previous edge of this member. */
  DEC_MEMBER childMember; /**< \brief Child member linked to this edge, or -1. */
} DecEdgeData;

typedef struct
{
  DEC_MEMBER_TYPE type;                 /**< \brief Type of member. Only valid for representative member. */
  DEC_MEMBER representativeMember;      /**< \brief Representative of member, or -1 if this is a representative member. */
  DEC_MEMBER parentMember;              /**< \brief Parent member of this member or -1 for a root. Only valid for representative member. */
  int numEdges;                         /**< \brief Number of edges. Only valid for representative member. */
  DEC_EDGE markerToParent;              /**< \brief Parent marker edge. Only valid for representative member. */
  DEC_EDGE markerOfParent;              /**< \brief Child marker of parent to which this member is linked. Only valid if root representative. */
  DEC_EDGE firstEdge;                   /**< \brief First edge in doubly-linked edge list of this member. */
  size_t lastParallelParentChildVisit;  /**< \brief Visit counter for \ref parallelParentChildCheckReducedMembers. */
} DEC_MEMBER_DATA;

typedef struct
{
  DEC_EDGE edge;  /**< \brief Edge corresponding to this row or -1. */
} DecRowData;

typedef struct
{
  DEC_EDGE edge;  /**< \brief Edge corresponding to this column or -1. */
} DecColumnData;

typedef struct
{
  CMR* cmr;                         /**< \brief \ref CMR environment. */

  size_t memMembers;                /**< \brief Allocated memory for members. */
  size_t numMembers;                /**< \brief Number of members. */
  DEC_MEMBER_DATA* members;         /**< \brief Array of members. */

  size_t memEdges;                  /**< \brief Allocated memory for edges. */
  size_t numEdges;                  /**< \brief Number of used edges. */
  DecEdgeData* edges;               /**< \brief Array of edges. */
  DEC_EDGE firstFreeEdge;           /**< \brief First edge in free list or -1. */

  size_t memNodes;                  /**< \brief Allocated memory for nodes. */
  size_t numNodes;                  /**< \brief Number of nodes. */
  DecNodeData* nodes;               /**< \brief Array of nodes. */
  DEC_NODE firstFreeNode;           /**< \brief First node in free list or -1. */

  size_t memRows;                   /**< \brief Allocated memory for \c rowEdges. */
  size_t numRows;                   /**< \brief Number of rows. */
  DecRowData* rowEdges;             /**< \brief Maps each row to its edge. */

  size_t memColumns;                /**< \brief Allocated memory for \c columnEdges. */
  size_t numColumns;                /**< \brief Number of columns. */
  DecColumnData* columnEdges;       /**< \brief Maps each column to its edge. */

  size_t numMarkerPairs;            /**< \brief Number of marker edge pairs in t-decomposition. */
  size_t parallelParentChildVisit;  /**< \brief Visit counter for \ref parallelParentChildCheckReducedMembers. */
} Dec;

/**
 * \brief Returns \c true if and only \p member is a representative member.
 */

static inline
bool isRepresentativeMember(
  Dec* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< Member of \p dec. */
)
{
  assert(dec);

  return dec->members[member].representativeMember == SIZE_MAX;
}

/**
 * \brief Returns the representative member of \p member.
 */

static
DEC_MEMBER findMember(
  Dec* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< Member for which the representative shall be returned. */
)
{
  DEC_MEMBER current = member;
  DEC_MEMBER next;
  while ((next = dec->members[current].representativeMember) != SIZE_MAX)
    current = next;
  DEC_MEMBER root = current;
  current = member;
  while ((next = dec->members[current].representativeMember) != SIZE_MAX)
  {
    if (next != root)
      dec->members[current].representativeMember = root;
    current = next;
  }
  return root;
}

/**
 * \brief Returns the representative of the parent member of \p member.
 *
 * Assumes that \p member is a representative member.
 */

static inline
DEC_MEMBER findMemberParent(
  Dec* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< Member whose parent shall be returned. */
)
{
  assert(isRepresentativeMember(dec, member));

  DEC_MEMBER someParent = dec->members[member].parentMember;
  if (someParent != SIZE_MAX)
    return findMember(dec, someParent);
  else
    return -1;
}

/**
 * \brief Returns the representative member of \p edge.
 */

static inline
DEC_MEMBER findEdgeMember(
  Dec* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge. */
)
{
  return findMember(dec, dec->edges[edge].member);
}

#if defined(CMR_DEBUG_CONSISTENCY)

/**
 * \brief Checks whether \p dec has consistent edge data.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyEdges(
  Dec* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    DEC_EDGE edge = dec->members[member].firstEdge;
    int countEdges = 0;
    if (edge != SIZE_MAX)
    {
      do
      {
        if (edge == SIZE_MAX || edge >= dec->memEdges)
          return CMRconsistencyMessage("edge %d of member %d out of range.", member, edge);
        if (dec->edges[edge].next == SIZE_MAX || dec->edges[edge].next > dec->memEdges)
          return CMRconsistencyMessage("edge %d of member %d has next out of range", member, edge);
        if (dec->edges[dec->edges[edge].next].prev != edge)
          return CMRconsistencyMessage("member %d has inconsistent edge list", member);
        if (findEdgeMember(dec, edge) != member)
          return CMRconsistencyMessage("edge %d belongs to member %d but is in member %d's edge list.", edge,
            findEdgeMember(dec, edge), member);
        edge = dec->edges[edge].next;
        countEdges++;
      }
      while (edge != dec->members[member].firstEdge);
    }
    if (countEdges != dec->members[member].numEdges)
    {
      return CMRconsistencyMessage("member %d has %d edges, but numEdges %d", member, countEdges,
        dec->members[member].numEdges);
    }
  }

  return NULL;
}

/**
 * \brief Checks whether \p dec has consistent member data.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyMembers(
  Dec* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (dec->members[member].type != DEC_MEMBER_TYPE_PARALLEL
      && dec->members[member].type != DEC_MEMBER_TYPE_RIGID
      && dec->members[member].type != DEC_MEMBER_TYPE_SERIES
      && dec->members[member].type != DEC_MEMBER_TYPE_LOOP)
    {
      return CMRconsistencyMessage("member %d has invalid type", member);
    }
  }

  return NULL;
}

/**
 * \brief Checks whether \p dec has consistent node data.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyNodes(
  Dec* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    bool isRigid = dec->members[member].type == DEC_MEMBER_TYPE_RIGID;
    DEC_EDGE edge = dec->members[member].firstEdge;
    if (edge == SIZE_MAX)
      continue;
    do
    {
      DEC_NODE head = dec->edges[edge].head;
      DEC_NODE tail = dec->edges[edge].tail;
      if (isRigid)
      {
        if (head == SIZE_MAX)
          return CMRconsistencyMessage("edge %d of rigid member %d has invalid head node", edge, member);
        if (tail == SIZE_MAX)
          return CMRconsistencyMessage("edge %d of rigid member %d has invalid tail node", edge, member);
        if (head >= dec->memNodes)
          return CMRconsistencyMessage("edge %d of rigid member %d has head node out of range", edge, member);
        if (tail >= dec->memNodes)
          return CMRconsistencyMessage("edge %d of rigid member %d has tail node out of range", edge, member);
      }
      edge = dec->edges[edge].next;
    }
    while (edge != dec->members[member].firstEdge);
  }

  return NULL;
}

/**
 * \brief Checks whether the members of \p dec form a forest.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyTree(
  Dec* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    int length = 0;
    DEC_MEMBER current;
    for (current = dec->members[member].parentMember; current != SIZE_MAX; current = dec->members[current].parentMember)
    {
      ++length;
      if (length > dec->numMembers)
        return "infinite member parent loop";
    }
  }

  return NULL;
}

/**
 * \brief Checks whether \p dec has consistent parent/child structure of members.
 *
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyParentChild(
  Dec* dec  /**< Decomposition. */
)
{
  assert(dec);

  if (dec->memMembers < dec->numMembers)
    return CMRconsistencyMessage("member count and memory are inconsistent");

  int* countChildren = NULL;
  if (CMRallocStackArray(dec->cmr, &countChildren, dec->memMembers) != CMR_OKAY)
    return CMRconsistencyMessage("stack allocation in consistencyParentChild() failed");
  for (int m = 0; m < dec->memMembers; ++m)
    countChildren[m] = 0;

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    if (dec->members[member].parentMember >= dec->memMembers)
    {
      CMRfreeStackArray(dec->cmr, &countChildren);
      return CMRconsistencyMessage("parent member of %d is out of range", member);
    }
    if (dec->members[member].parentMember != SIZE_MAX)
      countChildren[dec->members[member].parentMember]++;
  }

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    DEC_EDGE edge = dec->members[member].firstEdge;
    if (edge == SIZE_MAX)
      continue;
    do
    {
      if (dec->edges[edge].childMember != SIZE_MAX)
      {
        countChildren[member]--;

        if (findMember(dec, dec->members[findMember(dec, dec->edges[edge].childMember)].parentMember) != findMember(dec, member))
        {
          CMRfreeStackArray(dec->cmr, &countChildren);
          return CMRconsistencyMessage("member %d has child edge %d for child %d whose parent member is %d",
            member, edge, findMember(dec, dec->edges[edge].childMember),
            findMember(dec, dec->members[findMember(dec, dec->edges[edge].childMember)].parentMember));
        }
        if (dec->members[findMember(dec, dec->edges[edge].childMember)].markerOfParent != edge)
        {
          CMRfreeStackArray(dec->cmr, &countChildren);
          return CMRconsistencyMessage("member %d has child edge %d for child %d whose parent's markerOfParent is %d",
            member, edge, findMember(dec, dec->edges[edge].childMember),
            dec->members[findMember(dec, dec->edges[edge].childMember)].markerOfParent);
        }
        DEC_EDGE markerChild = dec->members[findMember(dec, dec->edges[edge].childMember)].markerToParent;
        if (dec->edges[markerChild].element != -dec->edges[edge].element)
        {
          CMRfreeStackArray(dec->cmr, &countChildren);
          return CMRconsistencyMessage("marker edges %d and %d of members %d (parent) and %d (child) have names %d and %d.",
            edge, markerChild, member, findEdgeMember(dec, markerChild), dec->edges[edge].element,
            dec->edges[markerChild].element);
        }
      }
      edge = dec->edges[edge].next;
    }
    while (edge != dec->members[member].firstEdge);
  }

  if (CMRfreeStackArray(dec->cmr, &countChildren) != CMR_OKAY)
    return "stack deallocation in consistencyParentChild() failed";

  return NULL;
}

/**
 * \brief Checks whether \p dec is consistent.
 *
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* decConsistency(
  Dec* dec  /**< Decomposition. */
)
{
  char* message = NULL;
  if ((message = consistencyMembers(dec)))
    return message;
  if ((message = consistencyEdges(dec)))
    return message;
  if ((message = consistencyNodes(dec)))
    return message;
  if ((message = consistencyParentChild(dec)))
    return message;
  if ((message = consistencyTree(dec)))
    return message;

  return NULL;
}

#endif /* CMR_DEBUG_CONSISTENCY */

typedef enum
{
  TYPE_CYCLE_CHILD = 1,   /**< Parent marker edge plus path edges form a cycle.
                           *   If graphic, this means we can replace the child by its child marker, adding the latter to
                           *   the list of path edges. */
  TYPE_SINGLE_CHILD = 2,  /**< One node of the parent marker edge is a path end.
                           *   If graphic, this means that one terminal node is in this member or its children. */
  TYPE_DOUBLE_CHILD = 3,  /**< Path edges form two cycles and adding the parent marker edge yields one.
                           *   If graphic, this means that both terminal nodes are in this member or its children. */
  TYPE_ROOT = 4           /**< Root member of reduced decomposition. */
} Type;

/**
 * \brief Additional information specific to a path edge.
 */

typedef struct _PathEdge
{
  DEC_EDGE edge;                  /**< \brief The actual edge in the decomposition. */
  struct _PathEdge* nextSibling;  /**< \brief Next edge of this reduced member, or \c NULL. */
  struct _PathEdge* nextOverall;  /**< \brief Next path edge in general. */
} PathEdge;

/**
 * \brief Additional member information specfic to a given path.
 *
 * TODO: Maybe add parent reduced member as well, so we don't have to go via the membersToReducedMembers array.
 */

typedef struct _ReducedMember
{
  DEC_MEMBER member;                /**< \brief The member from the decomposition. */
  DEC_MEMBER rootMember;            /**< \brief The root member of this component of the decomposition. */
  int depth;                        /**< \brief Depth of this member in the reduced decomposition. */
  Type type;                        /**< \brief Type of this member. */
  struct _ReducedMember* parent;    /**< \brief Parent in the reduced decomposition. */
  int numChildren;                  /**< \brief Number of children in the reduced decomposition. */
  struct _ReducedMember** children; /**< \brief Children in the reduced decomposition. */
  PathEdge* firstPathEdge;          /**< \brief First edge in linked list of path edges of \p member. */
  DEC_NODE rigidEndNodes[4];        /**< \brief For rigid members, the end nodes of the paths inside the member (or -1). */
} ReducedMember;

/**
 * \brief A component of the reduced decomposition.
 */

typedef struct _ReducedComponent
{
  int rootDepth;                /**< \brief Depth of reduced root member. */
  ReducedMember* root;          /**< \brief Reduced root member. */
  DEC_NODE terminalNode[2];     /**< \brief Terminal nodes of path. */
  DEC_MEMBER terminalMember[2]; /**< \brief Terminal members of path. */
  int numTerminals;
} ReducedComponent;

/**
 * \brief Additional information for each member (specific to a new column).
 */

typedef struct
{
  ReducedMember* reducedMember;       /**< \brief Reduced member corresponding to this member. */
  ReducedMember* rootDepthMinimizer;  /**< \brief Reduced root of this reduced component this member belongs to. */
} MemberInfo;

/**
 * \brief Information for adding a new column.
 */

typedef struct
{
  bool remainsGraphic;                  /**< \brief Indicator whether adding this column maintains graphicness. */
  size_t memReducedMembers;             /**< \brief Allocated memory for \c reducedMembers. */
  size_t numReducedMembers;             /**< \brief Number of members in \c reducedMembers. */
  ReducedMember* reducedMembers;        /**< \brief Array of reduced members, sorted by increasing depth. */
  MemberInfo* memberInfo;               /**< \brief Additional information for each member. */

  ReducedComponent* reducedComponents;  /**< \brief Array with reduced root members. */
  size_t memReducedComponents;          /**< \brief Allocated memory for \c reducedComponents. */
  size_t numReducedComponents;          /**< \brief Number of reduced root members. */

  PathEdge* pathEdges;                  /**< \brief Storage for edge lists of path edges. */
  size_t memPathEdges;                  /**< \brief Allocated memory for \c pathEdges. */
  size_t numPathEdges;                  /**< \brief Number of stored edges in \c pathEdges. */

  ReducedMember** childrenStorage;      /**< \brief Storage for reduced members' arrays of children. */
  size_t usedChildrenStorage;           /**< \brief Number of stored children in \c childrenStorage. */
  size_t memChildrenStorage;            /**< \brief Allocated memory for \c childrenStorage. */

  int* nodesDegree;                     /**< \brief Map from nodes to degree w.r.t. path edges. */
  size_t memNodesDegree;                /**< \brief Allocated memory for \c nodesDegree. */

  bool* edgesInPath;                    /**< \brief Map from edges to indicator for being in the path. */
  size_t memEdgesInPath;                /**< \brief Allocated memory for \p edgesInPath. */
  PathEdge* firstPathEdge;              /**< \brief Root of singly-linked list of all path edges. */
} DEC_NEWCOLUMN;

/**
 * \brief Returns the representative node of \p node.
 *
 * Assumes \p node is indeed a node, i.e., not -1.
 */

static
DEC_NODE findNode(
  Dec* dec,
  DEC_NODE node
)
{
  DEC_NODE current = node;
  DEC_NODE next;
  while ((next = dec->nodes[current].representativeNode) != SIZE_MAX)
    current = next;
  DEC_NODE root = current;
  current = node;
  while ((next = dec->nodes[current].representativeNode) != SIZE_MAX)
  {
    if (next != root)
      dec->nodes[current].representativeNode = root;
    current = next;
  }
  return root;
}

/**
 * \brief Returns the representative node of the tail of \p edge.
 */

static inline
DEC_NODE findEdgeTail(
  Dec* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge. */
)
{
  assert(dec);
  assert(edge < dec->memEdges);
  assert(dec->edges[edge].tail < dec->memNodes);

  return findNode(dec, dec->edges[edge].tail);
}

/**
 * \brief Returns the representative node of the head of \p edge.
 */

static inline
DEC_NODE findEdgeHead(
  Dec* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge. */
)
{
  assert(dec);
  assert(edge < dec->memEdges);
  assert(dec->edges[edge].head < dec->memNodes);

  return findNode(dec, dec->edges[edge].head);
}

/**
 * \brief Creates a node for some rigid member of the decomposition \p dec.
 */

static
CMR_ERROR createNode(
  Dec* dec,       /**< Decomposition. */
  DEC_NODE* pnode /**< Pointer for storing new node. */
)
{
  assert(dec);
  assert(pnode);

  DEC_NODE node = dec->firstFreeNode;
  if (node != SIZE_MAX)
  {
    CMRdbgMsg(10, "createNode returns free node %d.\n", node);
    dec->firstFreeNode = dec->nodes[node].representativeNode;
  }
  else
  {
    /* No node in free list, so we enlarge the array. */

    int newSize = 2 * dec->memNodes + 16;
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &dec->nodes, newSize) );
    for (int v = dec->memNodes + 1; v < newSize; ++v)
      dec->nodes[v].representativeNode = v+1;
    dec->nodes[newSize-1].representativeNode = SIZE_MAX;
    dec->firstFreeNode = dec->memNodes + 1;
    node = dec->memNodes;
    dec->memNodes = newSize;
    CMRdbgMsg(12, "createNode enlarges node array to %d and returns node %d.\n", newSize, node);
  }
  dec->nodes[node].representativeNode = SIZE_MAX;
  dec->numNodes++;

  *pnode = node;

  return CMR_OKAY;
}

/**
 * \brief Adds \p edge to the edge list of its member.
 */

static
CMR_ERROR addEdgeToMembersEdgeList(
  Dec* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge to be added. */
)
{
  assert(dec);
  assert(edge != SIZE_MAX);

  DEC_MEMBER member = findEdgeMember(dec, edge);
  DEC_EDGE first = dec->members[member].firstEdge;
  if (first != SIZE_MAX)
  {
    assert(dec->members[member].numEdges > 0);
    DEC_EDGE last = dec->edges[first].prev;
    dec->edges[edge].next = first;
    dec->edges[edge].prev = last;
    dec->edges[first].prev = edge;
    dec->edges[last].next =  edge;
  }
  else
  {
    assert(dec->members[member].numEdges == 0);
    dec->edges[edge].next = edge;
    dec->edges[edge].prev = edge;
  }
  dec->members[member].firstEdge = edge;
  dec->members[member].numEdges++;

  return CMR_OKAY;
}

/**
 * \brief Removes \p edge from the edge list of its member.
 */

static
CMR_ERROR removeEdgeFromMembersEdgeList(
  Dec* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge to be added. */
)
{
  assert(dec);

  DEC_MEMBER member = findEdgeMember(dec, edge);
  if (dec->members[member].numEdges == 1)
    dec->members[member].firstEdge = -1;
  else
  {
    if (dec->members[member].firstEdge == edge)
      dec->members[member].firstEdge = dec->edges[edge].next;

    assert(dec->members[member].firstEdge != edge);

    dec->edges[dec->edges[edge].prev].next = dec->edges[edge].next;
    dec->edges[dec->edges[edge].next].prev = dec->edges[edge].prev;
  }

  dec->members[member].numEdges--;

  return CMR_OKAY;
}


/**
 * \brief Replaced \p oldEdge in its member by \p newEdge.
 *
 * Assumes that \p newEdge has the same member but is not in its edge list.
 */

static
CMR_ERROR replaceEdgeInMembersEdgeList(
  Dec* dec,         /**< Decomposition. */
  DEC_EDGE oldEdge, /**< Edge to be removed. */
  DEC_EDGE newEdge  /**< Edge to be added. */
)
{
  assert(dec);

  DEC_MEMBER member = findEdgeMember(dec, oldEdge);
  assert(findEdgeMember(dec, newEdge) == member);

  dec->edges[newEdge].tail = dec->edges[oldEdge].tail;
  dec->edges[newEdge].head = dec->edges[oldEdge].head;
  dec->edges[newEdge].next = dec->edges[oldEdge].next;
  dec->edges[newEdge].prev = dec->edges[oldEdge].prev;
  dec->edges[dec->edges[oldEdge].next].prev = newEdge;
  dec->edges[dec->edges[oldEdge].prev].next = newEdge;
  if (dec->members[member].firstEdge == oldEdge)
    dec->members[member].firstEdge = newEdge;

  return CMR_OKAY;
}

/**
 * \brief Creates a new edge in \p member.
 */

static
CMR_ERROR createEdge(
  Dec* dec,           /**< Decomposition. */
  DEC_MEMBER member,  /**< Member this edge belongs to. */
  DEC_EDGE* pedge     /**< Pointer for storing the new edge. */
)
{
  assert(dec);
  assert(pedge);
  assert(member == SIZE_MAX || isRepresentativeMember(dec, member));

  DEC_EDGE edge = dec->firstFreeEdge;
  if (edge != SIZE_MAX)
  {
    CMRdbgMsg(12, "Creating edge %d by using a free edge.\n", edge);
    dec->firstFreeEdge = dec->edges[edge].next;
  }
  else /* No edge in free list, so we enlarge the array. */
  {
    int newSize = 2 * dec->memEdges + 16;
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &dec->edges, newSize) );
    for (int e = dec->memEdges + 1; e < newSize; ++e)
    {
      dec->edges[e].next = e+1;
      dec->edges[e].member = -1;
    }
    dec->edges[newSize-1].next = -1;
    dec->firstFreeEdge = dec->memEdges + 1;
    edge = dec->memEdges;
    dec->memEdges = newSize;
    CMRdbgMsg(12, "Creating edge %d and reallocating edge to array %d elements.\n", edge, newSize);
  }

  dec->edges[edge].tail = -1;
  dec->edges[edge].head = -1;
  dec->edges[edge].element = 0;
  dec->edges[edge].member = member;
  dec->numEdges++;

  *pedge = edge;

  return CMR_OKAY;
}

/**
 * \brief Creates a pair of marker edges for linking parent \p parentMember to child \p childMember.
 */

static
CMR_ERROR createMarkerEdgePair(
  Dec* dec,                     /**< Decomposition. */
  DEC_MEMBER parentMember,      /**< Parent member. */
  DEC_EDGE* pMarkerOfParent,    /**< Pointer for storing the new child marker in the parent. */
  DEC_NODE markerOfParentTail,  /**< Tail node of this *\p pChildMarker. */
  DEC_NODE markerOfParentHead,  /**< Head node of this *\p pChildMarker. */
  DEC_MEMBER childMember,       /**< Child member. */
  DEC_EDGE* pMarkerToParent,    /**< Pointer for storing the new parent marker in the child. */
  DEC_NODE markerToParentTail,  /**< Tail node of this *\p pParentMarker. */
  DEC_NODE markerToParentHead   /**< Head node of this *\p pParentMarker. */
)
{
  assert(dec);
  assert(pMarkerOfParent);
  assert(pMarkerToParent);
  assert(isRepresentativeMember(dec, parentMember));
  assert(isRepresentativeMember(dec, childMember));

  /* Create the child marker edge of the parent member. */

  CMR_CALL( createEdge(dec, parentMember, pMarkerOfParent) );
  DecEdgeData* data = &dec->edges[*pMarkerOfParent];
  data->tail = markerOfParentTail;
  data->head = markerOfParentHead;
  data->childMember = childMember;
  data->element = -INT_MAX + dec->numMarkerPairs;
  CMRdbgMsg(12, "Created child marker edge {%d,%d} <%s> of parent member %d.\n", markerOfParentTail, markerOfParentHead,
    CMRelementString(data->element, NULL), parentMember);

  /* Create the parent marker edge of the child member. */

  CMR_CALL( createEdge(dec, childMember, pMarkerToParent) );
  data = &dec->edges[*pMarkerToParent];
  data->tail = markerToParentTail;
  data->head = markerToParentHead;
  data->childMember = -1;
  dec->members[childMember].parentMember = parentMember;
  dec->members[childMember].markerOfParent = *pMarkerOfParent;
  dec->members[childMember].markerToParent = *pMarkerToParent;
  data->element = INT_MAX - dec->numMarkerPairs;
  CMRdbgMsg(12, "Created child marker edge {%d,%d} <%s> of child member %d.\n", markerToParentTail, markerToParentHead,
    CMRelementString(data->element, NULL), childMember);

  /* Increase counter of used marker pairs. */
  dec->numMarkerPairs++;

  return CMR_OKAY;
}

/**
 * \brief Creates a new member.
 */

static
CMR_ERROR createMember(
  Dec* dec,             /**< Decomposition. */
  DEC_MEMBER_TYPE type, /**< Type of member. */
  DEC_MEMBER* pmember   /**< Pointer for storing the new member. */
)
{
  assert(dec);

  if (dec->numMembers == dec->memMembers)
  {
    dec->memMembers = 16 + 2 * dec->memMembers;
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &dec->members, dec->memMembers) );
  }

  DEC_MEMBER_DATA* data = &dec->members[dec->numMembers];
  data->markerOfParent = -1;
  data->markerToParent = -1;
  data->firstEdge = -1;
  data->representativeMember = -1;
  data->numEdges = 0;
  data->parentMember = -1;
  data->type = type;
  data->lastParallelParentChildVisit = 0;
  *pmember = dec->numMembers;
  dec->numMembers++;

  CMRdbgMsg(10, "Creating %s member %d.\n", memberTypeString(type), *pmember);

  return CMR_OKAY;
}

/**
 * \brief Creates an empty decomposition.
 */

CMR_ERROR decCreate(
  CMR* cmr,           /**< \ref CMR environment. */
  Dec** pdec,       /**< Pointer to new decomposition. .*/
  size_t memEdges,     /**< Initial memory for edges of the decomposition. */
  size_t memNodes,     /**< Initial memory for nodes of the decomposition. */
  size_t memMembers,   /**< Initial memory for members of the decomposition. */
  size_t memRows,      /**< Initial memory for rows. */
  size_t memColumns    /**< Initial memory for columns. */
)
{
  assert(cmr);
  assert(pdec);
  assert(!*pdec);

  CMR_CALL( CMRallocBlock(cmr, pdec) );
  Dec* dec = *pdec;
  dec->cmr = cmr;
  dec->memMembers = memMembers;
  dec->numMembers = 0;
  dec->members = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->members, dec->memMembers) );

  if (memNodes < 1)
    memNodes = 1;
  dec->memNodes = memNodes;
  dec->nodes = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->nodes, memNodes) );
  dec->numNodes = 0;
  for (size_t v = 0; v < memNodes; ++v)
    dec->nodes[v].representativeNode = v+1;
  dec->nodes[memNodes-1].representativeNode = SIZE_MAX;
  dec->firstFreeNode = 0;

  if (memEdges < 1)
    memEdges = 1;
  dec->memEdges = memEdges;
  dec->edges = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->edges, memEdges) );
  dec->numEdges = 0;
  dec->numMarkerPairs = 0;
  dec->parallelParentChildVisit = 0;

  /* Initialize free list with unused edges. */
  if (memEdges > dec->numEdges)
  {
    for (size_t e = dec->numEdges; e < memEdges; ++e)
    {
      dec->edges[e].next = e+1;
      dec->edges[e].member = -1;
    }
    dec->edges[memEdges-1].next = -1;
    dec->firstFreeEdge = dec->numEdges;
  }
  else
    dec->firstFreeEdge = -1;

  dec->numRows = 0;
  dec->memRows = memRows;
  dec->rowEdges = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->rowEdges, dec->memRows) );
  for (size_t r = 0; r < dec->numRows; ++r)
    dec->rowEdges[r].edge = -1;

  dec->numColumns = 0;
  dec->memColumns = memColumns;
  dec->columnEdges = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &dec->columnEdges, dec->memColumns) );
  for (size_t c = 0; c < dec->numColumns; ++c)
    dec->columnEdges[c].edge = -1;

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

/**
 * \brief Frees the decomposition \p *pdec.
 */

CMR_ERROR decFree(
  Dec** pdec /**< Pointer to decomposition. */
)
{
  assert(pdec);
  assert(*pdec);

  Dec* dec = *pdec;
  CMR_CALL( CMRfreeBlockArray(dec->cmr, &dec->members) );
  CMR_CALL( CMRfreeBlockArray(dec->cmr, &dec->edges) );
  CMR_CALL( CMRfreeBlockArray(dec->cmr, &dec->nodes) );
  CMR_CALL( CMRfreeBlockArray(dec->cmr, &dec->rowEdges) );
  CMR_CALL( CMRfreeBlockArray(dec->cmr, &dec->columnEdges) );
  CMR_CALL( CMRfreeBlock(dec->cmr, pdec) );

  return CMR_OKAY;
}

/**
 * \brief Creates a graph represented by given decomposition.
 */

static
CMR_ERROR decToGraph(
  Dec* dec,                       /**< Decomposition. */
  CMR_GRAPH* graph,               /**< Graph to be filled. */
  bool merge,                     /**< Merge and remove corresponding parent and child markers. */
  CMR_GRAPH_EDGE* forestEdges,    /**< If not \c NULL, the edges of a spanning tree are stored here. */
  CMR_GRAPH_EDGE* coforestEdges,  /**< If not \c NULL, the non-basis edges are stored here. */
  CMR_ELEMENT* edgeElements       /**< If not \c NULL, the elements for each edge are stored here. */
)
{
  assert(dec);
  assert(graph);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */


  CMRdbgMsg(0, "CMRdecToGraph for t-decomposition.\n");

  CMR_CALL( CMRgraphClear(dec->cmr, graph) );

  CMR_ELEMENT* localEdgeElements = NULL;
  if (edgeElements)
    localEdgeElements = edgeElements;
  else if (forestEdges || coforestEdges)
    CMR_CALL( CMRallocStackArray(dec->cmr, &localEdgeElements, dec->memEdges) );
  CMR_GRAPH_NODE* decNodesToGraphNodes = NULL;
  CMR_CALL( CMRallocStackArray(dec->cmr, &decNodesToGraphNodes, dec->memNodes) );
  CMR_GRAPH_EDGE* decEdgesToGraphEdges = NULL;
  CMR_CALL( CMRallocStackArray(dec->cmr, &decEdgesToGraphEdges, dec->memEdges) );

  for (size_t v = 0; v < dec->memNodes; ++v)
  {
    if (dec->nodes[v].representativeNode == SIZE_MAX)
    {
      CMR_CALL( CMRgraphAddNode(dec->cmr, graph, &decNodesToGraphNodes[v]) );
    }
    else
      decNodesToGraphNodes[v] = -1;
  }

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    DEC_MEMBER_TYPE type = dec->members[member].type;
    CMRdbgMsg(2, "Member %d is %s with %d edges.\n", member, memberTypeString(type), dec->members[member].numEdges);

    CMR_GRAPH_EDGE graphEdge;
    DEC_EDGE edge = dec->members[member].firstEdge;
    if (type == DEC_MEMBER_TYPE_RIGID)
    {
      do
      {
        DEC_NODE head = findEdgeHead(dec, edge);
        DEC_NODE tail = findEdgeTail(dec, edge);
        CMR_CALL( CMRgraphAddEdge(dec->cmr, graph, decNodesToGraphNodes[head], decNodesToGraphNodes[tail],
          &graphEdge) );
        decEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = dec->edges[edge].element;
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (type == DEC_MEMBER_TYPE_PARALLEL)
    {
      CMR_GRAPH_NODE graphHead, graphTail;
      CMR_CALL( CMRgraphAddNode(dec->cmr, graph, &graphHead) );
      CMR_CALL( CMRgraphAddNode(dec->cmr, graph, &graphTail) );
      do
      {
        CMR_CALL( CMRgraphAddEdge(dec->cmr, graph, graphHead, graphTail, &graphEdge) );
        decEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = dec->edges[edge].element;
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (type == DEC_MEMBER_TYPE_SERIES)
    {
      CMR_GRAPH_NODE firstNode, v;
      CMR_CALL( CMRgraphAddNode(dec->cmr, graph, &firstNode) );
      v = firstNode;
      edge = dec->edges[edge].next;
      while (edge != dec->members[member].firstEdge)
      {
        CMR_GRAPH_NODE w;
        CMR_CALL( CMRgraphAddNode(dec->cmr, graph, &w) );
        CMR_CALL( CMRgraphAddEdge(dec->cmr, graph, v, w, &graphEdge) );
        decEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = dec->edges[edge].element;

        edge = dec->edges[edge].next;
        v = w;
      }
      CMR_CALL( CMRgraphAddEdge(dec->cmr, graph, v, firstNode, &graphEdge) );
      decEdgesToGraphEdges[edge] = graphEdge;
      if (localEdgeElements)
        localEdgeElements[graphEdge] = dec->edges[edge].element;
    }
    else
    {
      assert(type == DEC_MEMBER_TYPE_LOOP);

      CMR_GRAPH_NODE v;
      CMR_CALL( CMRgraphAddNode(dec->cmr, graph, &v) );
      CMR_CALL( CMRgraphAddEdge(dec->cmr, graph, v, v, &graphEdge) );
      decEdgesToGraphEdges[edge] = graphEdge;
      if (localEdgeElements)
        localEdgeElements[graphEdge] = dec->edges[edge].element;
    }
  }

  /* Merge respective parent and child edges. */

  if (merge)
  {
    CMRdbgMsg(2, "Before merging, the graph has %d nodes and %d edges.\n", CMRgraphNumNodes(graph),
      CMRgraphNumEdges(graph));

    for (size_t m = 0; m < dec->numMembers; ++m)
    {
      if (!isRepresentativeMember(dec, m) || dec->members[m].parentMember == SIZE_MAX)
        continue;

      CMR_GRAPH_EDGE parent = decEdgesToGraphEdges[dec->members[m].markerOfParent];
      CMR_GRAPH_EDGE child = decEdgesToGraphEdges[dec->members[m].markerToParent];
      CMR_GRAPH_NODE parentU = CMRgraphEdgeU(graph, parent);
      CMR_GRAPH_NODE parentV = CMRgraphEdgeV(graph, parent);
      CMR_GRAPH_NODE childU = CMRgraphEdgeU(graph, child);
      CMR_GRAPH_NODE childV = CMRgraphEdgeV(graph, child);

      CMRdbgMsg(2, "Merging edges %d = {%d,%d} <%s>", parent, parentU, parentV,
        CMRelementString(dec->edges[dec->members[m].markerOfParent].element, NULL));
      CMRdbgMsg(0, " and %d = {%d,%d} <%s>.\n", child, childU, childV,
        CMRelementString(dec->edges[dec->members[m].markerToParent].element, NULL));

      CMR_CALL( CMRgraphMergeNodes(dec->cmr, graph, parentU, childU) );
      CMR_CALL( CMRgraphDeleteNode(dec->cmr, graph, childU) );
      CMR_CALL( CMRgraphMergeNodes(dec->cmr, graph, parentV, childV) );
      CMR_CALL( CMRgraphDeleteNode(dec->cmr, graph, childV) );

      CMR_CALL( CMRgraphDeleteEdge(dec->cmr, graph, parent) );
      CMR_CALL( CMRgraphDeleteEdge(dec->cmr, graph, child) );
    }
  }

  // TODO: Remove nodes with degree 0 or 1?!
  bool* nodesUsed = NULL;
  CMR_CALL( CMRallocStackArray(dec->cmr, &nodesUsed, CMRgraphMemNodes(graph)) );
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); v = CMRgraphNodesNext(graph, v))
    nodesUsed[v] = false;
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, i); i = CMRgraphEdgesNext(graph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, i);
    nodesUsed[CMRgraphEdgeU(graph, e)] = true;
    nodesUsed[CMRgraphEdgeV(graph, e)] = true;
  }
  for (CMR_GRAPH_NODE v = CMRgraphNodesFirst(graph); CMRgraphNodesValid(graph, v); )
  {
    CMR_GRAPH_NODE next = CMRgraphNodesNext(graph, v);
    if (!nodesUsed[v])
    {
      CMRdbgMsg(2, "Removing degree-0 node %d.\n", v);
      CMRgraphDeleteNode(dec->cmr, graph, v);
    }
    v = next;
  }
  CMR_CALL( CMRfreeStackArray(dec->cmr, &nodesUsed) );

  /* Construct (co)forest. */

  if (forestEdges || coforestEdges)
  {
#if !defined(NDEBUG)
    /* This is only relevant if a 1-separation exists. */
    for (size_t r = 0; r < dec->numRows; ++r)
      forestEdges[r] = INT_MIN;
    for (size_t c = 0; c < dec->numColumns; ++c)
      coforestEdges[c] = INT_MIN;
#endif /* !NDEBUG */

    for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, i);
      i = CMRgraphEdgesNext(graph, i))
    {
      CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, i);

      CMRdbgMsg(2, "Graph edge %d = {%d,%d} <%s>\n", e, CMRgraphEdgeU(graph, e), CMRgraphEdgeV(graph, e),
        CMRelementString(localEdgeElements[e], NULL));

      CMR_ELEMENT element = localEdgeElements[e];
      if (CMRelementIsRow(element) && forestEdges)
      {
        CMRdbgMsg(0, "Edge corresponds to element %d = row %d.\n", element, CMRelementToRowIndex(element));
        forestEdges[CMRelementToRowIndex(element)] = e;
      }
      else if (CMRelementIsColumn(element) && coforestEdges)
      {
        CMRdbgMsg(0, "Edge corresponds to element %d = column %d.\n", element, CMRelementToColumnIndex(element));
        coforestEdges[CMRelementToColumnIndex(element)] = e;
      }
    }

#if !defined(NDEBUG)
    /* These assertions indicate a 1-separable input matrix. */
    for (size_t r = 0; r < dec->numRows; ++r)
      assert(forestEdges[r] >= 0);
    for (size_t c = 0; c < dec->numColumns; ++c)
      assert(coforestEdges[c] >= 0);
#endif /* !NDEBUG */
  }

  CMR_CALL( CMRfreeStackArray(dec->cmr, &decEdgesToGraphEdges) );
  CMR_CALL( CMRfreeStackArray(dec->cmr, &decNodesToGraphNodes) );
  if (localEdgeElements != edgeElements)
    CMR_CALL( CMRfreeStackArray(dec->cmr, &localEdgeElements) );

  return CMR_OKAY;
}

/**
 * Prints an edge in \c dot format to \p stream.
 */

static
void edgeToDot(
  FILE* stream,       /**< File stream. */
  Dec* dec,           /**< Decomposition. */
  DEC_MEMBER member,  /**< Member this edge belongs to. */
  DEC_EDGE edge,      /**< Edge. */
  int u,              /**< First node. */
  int v,              /**< Second node. */
  bool red            /**< Whether to color it red. */
)
{
  assert(stream);

  member = findMember(dec, member);

  char type = (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL) ?
    'P' : (dec->members[member].type == DEC_MEMBER_TYPE_SERIES ? 'S' : 'R');
  const char* redStyle = red ? ",color=red" : "";
  if (dec->members[member].markerToParent == edge)
  {
    fprintf(stream, "    %c_%zu_%d -- %c_p_%zu [label=\"%zu\",style=dashed%s];\n", type, member, u, type, member, edge, redStyle);
    fprintf(stream, "    %c_p_%zu -- %c_%zu_%d [label=\"%zu\",style=dashed%s];\n", type, member, type, member, v, edge, redStyle);
    fprintf(stream, "    %c_%zu_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%zu_%d [shape=box];\n", type, member, v);
    fprintf(stream, "    %c_p_%zu [style=dashed];\n", type, member);
  }
  else if (dec->edges[edge].childMember != SIZE_MAX)
  {
    DEC_MEMBER child = findMember(dec, dec->edges[edge].childMember);
    char childType = (dec->members[child].type == DEC_MEMBER_TYPE_PARALLEL) ?
      'P' : (dec->members[child].type == DEC_MEMBER_TYPE_SERIES ? 'S' : 'R');
    fprintf(stream, "    %c_%zu_%d -- %c_c_%zu [label=\"%zu\",style=dotted%s];\n", type, member, u, type, child, edge, redStyle);
    fprintf(stream, "    %c_c_%zu -- %c_%zu_%d [label=\"%zu\",style=dotted%s];\n", type, child, type, member, v, edge, redStyle);
    fprintf(stream, "    %c_%zu_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%zu_%d [shape=box];\n", type, member, v);
    fprintf(stream, "    %c_c_%zu [style=dotted];\n", type, child);

    fprintf(stream, "    %c_p_%zu -- %c_c_%zu [style=dashed,dir=forward];\n", childType, child, type, child);
  }
  else
  {
    fflush(stdout);
    fprintf(stream, "    %c_%zu_%d -- %c_%zu_%d [label=\"%zu <%s>\",style=bold%s];\n", type, member, u, type, member, v,
      edge, CMRelementString(dec->edges[edge].element, NULL), redStyle);
    fprintf(stream, "    %c_%zu_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%zu_%d [shape=box];\n", type, member, v);
  }
}

/**
 * \brief Visualizes a decomposition in \c dot format.
 */

CMR_ERROR CMRdecToDot(
  Dec* dec,               /**< Decomposition. */
  FILE* stream,           /**< Stream to write to. */
  bool* edgesHighlighted  /**< Indicator for edges to be highlighted. */
)
{
  assert(dec);
  assert(stream);

  fprintf(stream, "// t-decomposition\n");
  fprintf(stream, "graph dec {\n");
  fprintf(stream, "  compound = true;\n");
  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    fprintf(stream, "  subgraph member%zu {\n", member);
    if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
    {
      DEC_EDGE edge = dec->members[member].firstEdge;
      do
      {
        edgeToDot(stream, dec, member, edge, 0, 1, edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
    {
      DEC_EDGE edge = dec->members[member].firstEdge;
      do
      {
        DEC_NODE u = findEdgeHead(dec, edge);
        DEC_NODE v = findEdgeTail(dec, edge);
        edgeToDot(stream, dec, member, edge, u, v,
          edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (dec->members[member].type == DEC_MEMBER_TYPE_SERIES)
    {
      DEC_EDGE edge = dec->members[member].firstEdge;
      int i = 0;
      do
      {
        edgeToDot(stream, dec, member, edge, i, (i+1) % dec->members[member].numEdges,
          edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = dec->edges[edge].next;
        i++;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else
    {
      assert(dec->members[member].type == DEC_MEMBER_TYPE_LOOP);

      edgeToDot(stream, dec, member, dec->members[member].firstEdge, 0, 0, false);
    }
    fprintf(stream, "  }\n");
  }
  fprintf(stream, "}\n");

  return CMR_OKAY;
}

/**<
 * \brief Writes decomposition to a dot file in the current directory.
 *
 * The file name is dec-NUMBER.dot where NUMBER is increased for each call and reset for a new decomposition.
 */

static inline
void debugDot(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
  CMR_UNUSED(dec);
  CMR_UNUSED(newcolumn);

  assert(dec);
  assert(newcolumn || !newcolumn);

#if defined(CMR_DEBUG_DOT)
  static int dotFileCounter = 1;
  char name[256];
  snprintf(name, 256, "dec-%03d.dot", dotFileCounter);
  CMRdbgMsg(0, "Writing <%s>...", name);
  FILE* dotFile = fopen(name, "w");
  CMRdecToDot(dec, dotFile, newcolumn ? newcolumn->edgesInPath : NULL);
  fclose(dotFile);
  CMRdbgMsg(0, " done.\n");

  dotFileCounter++;
#endif /* CMR_DEBUG_DOT */
}

/**
 * \brief Creates a \ref DEC_NEWCOLUMN structure.
 */

CMR_ERROR newcolumnCreate(
  CMR* cmr,                     /**< \ref CMR environment. */
  DEC_NEWCOLUMN** pnewcolumn  /**< Pointer for storing the newcolumn structure. */
)
{
  assert(cmr);

  CMR_CALL( CMRallocBlock(cmr, pnewcolumn) );
  DEC_NEWCOLUMN* newcolumn = *pnewcolumn;
  newcolumn->remainsGraphic = true;
  newcolumn->memReducedMembers = 0;
  newcolumn->numReducedMembers = 0;
  newcolumn->reducedMembers = NULL;
  newcolumn->memberInfo = NULL;

  newcolumn->numReducedComponents = 0;
  newcolumn->memReducedComponents = 0;
  newcolumn->reducedComponents = NULL;

  newcolumn->pathEdges = NULL;
  newcolumn->memPathEdges = 0;
  newcolumn->numPathEdges = 0;

  newcolumn->memChildrenStorage = 0;
  newcolumn->usedChildrenStorage = 0;
  newcolumn->childrenStorage = NULL;

  newcolumn->nodesDegree = NULL;
  newcolumn->memNodesDegree = 0;

  newcolumn->edgesInPath = NULL;
  newcolumn->memEdgesInPath = 0;
  newcolumn->firstPathEdge = NULL;

  return CMR_OKAY;
}

/**
 * \brief Frees a \ref DEC_NEWCOLUMN structure.
 */

CMR_ERROR newcolumnFree(
  CMR* cmr,                     /**< \ref CMR environment. */
  DEC_NEWCOLUMN** pnewcolumn  /**< Pointer to newcolumn structure. */
)
{
  assert(cmr);
  assert(*pnewcolumn);

  DEC_NEWCOLUMN* newcolumn = *pnewcolumn;

  if (newcolumn->edgesInPath)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->edgesInPath) );

  if (newcolumn->nodesDegree)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->nodesDegree) );

  if (newcolumn->reducedComponents)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->reducedComponents) );
  if (newcolumn->reducedMembers)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->reducedMembers) );

  if (newcolumn->memberInfo)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->memberInfo) );
  if (newcolumn->pathEdges)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->pathEdges) );
  if (newcolumn->childrenStorage)
    CMR_CALL( CMRfreeBlockArray(cmr, &newcolumn->childrenStorage) );

  CMR_CALL( CMRfreeBlock(cmr, pnewcolumn) );

  return CMR_OKAY;
}


/**
 * \brief Removes all path edges.
 */

static
CMR_ERROR removeAllPathEdges(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
  assert(dec);
  assert(newcolumn);

  for (PathEdge* pathEdge = newcolumn->firstPathEdge; pathEdge; pathEdge = pathEdge->nextOverall)
  {
    DEC_EDGE edge = pathEdge->edge;
    DEC_NODE tail = dec->edges[edge].tail;
    if (tail != SIZE_MAX)
      tail = findNode(dec, tail);
    DEC_NODE head = dec->edges[edge].head;
    if (head != SIZE_MAX)
      head = findNode(dec, head);
    CMRdbgMsg(6, "Removing path edge %d = {%d,%d}.\n", pathEdge->edge, tail, head);
    newcolumn->edgesInPath[edge] = false;

    DEC_MEMBER member = findEdgeMember(dec, edge);
    if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
    {
      DEC_NODE tail = findEdgeTail(dec, edge);
      if (tail == SIZE_MAX)
        continue;
      DEC_NODE head = findEdgeHead(dec, edge);
      if (head == SIZE_MAX)
        continue;
      newcolumn->nodesDegree[tail] = 0;
      newcolumn->nodesDegree[head] = 0;
      CMRdbgMsg(8, "Set nodesDegree of nodes of edge %d = {%d,%d} to 0.\n", edge, tail, head);
    }
  }
  newcolumn->firstPathEdge = NULL;
  newcolumn->numPathEdges = 0;

  return CMR_OKAY;
}

/**
 * \brief Initializes a \ref DEC_NEWCOLUMN structure in order to check for a column.
 */

static
CMR_ERROR initializeNewColumn(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn  /**< newcolumn. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(newcolumn->numReducedMembers == 0);

  newcolumn->remainsGraphic = true;

#if defined(CMR_DEBUG)
  for (size_t i = 0; i < newcolumn->memNodesDegree; ++i)
  {
    if (dec->nodes[i].representativeNode != SIZE_MAX)
      continue;

    if (newcolumn->nodesDegree[i] != 0)
      CMRdbgMsg(6, "Node %d still has degree %d after clean-up.\n", i, newcolumn->nodesDegree[i]);
    assert(newcolumn->nodesDegree[i] == 0);
  }
#endif /* CMR_DEBUG || !NDEBUG */

  /* Enlarge nodesDegree array and initialize new portion to zero. */
  if (newcolumn->memNodesDegree < dec->memNodes + dec->memEdges)
  {
    size_t newSize = 16 + 2 * newcolumn->memNodesDegree;
    while (newSize < dec->memNodes + dec->memEdges)
      newSize *= 2;
    CMRdbgMsg(6, "Enlarging nodesDegree array to %ld.\n", newSize);
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->nodesDegree, newSize) );
    for (size_t i = newcolumn->memNodesDegree; i < newSize; ++i)
      newcolumn->nodesDegree[i] = 0;
    newcolumn->memNodesDegree = newSize;
  }

  return CMR_OKAY;
}

/**
 * \brief Ensures that the child marker in \p member for \p childMember is not parallel to the parent marker of
 * \p member.
 *
 * Note that this can only happen when a component is reordered (receiving a new root) because reduced components shall
 * be joined.
 */

static
CMR_ERROR parallelParentChildCheckMember(
  Dec* dec,              /**< Decomposition. */
  DEC_MEMBER member,     /**< Parent of \p childMember. */
  DEC_MEMBER childMember /**< Child member of \p member. */
)
{
  assert(dec);
  assert(childMember != SIZE_MAX);
  assert(member == findMemberParent(dec, childMember));

  /* Stop if childMember was already processed at some point. */
  if (dec->members[childMember].lastParallelParentChildVisit == dec->parallelParentChildVisit)
    return CMR_OKAY;

  dec->members[childMember].lastParallelParentChildVisit = dec->parallelParentChildVisit;
  DEC_MEMBER parentMember = findMemberParent(dec, member);
  CMRdbgMsg(10, "Consider child member %d of %d with parent %d.\n", childMember, member, parentMember);
  if (parentMember == SIZE_MAX)
    return CMR_OKAY;

  if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
  {
    CMRdbgMsg(10, "Checking if the child marker of %d for child member %d is not parallel to %d's parent marker.\n",
      member, childMember, member);

    DEC_EDGE childMarkerEdge = dec->members[childMember].markerOfParent;
    DEC_NODE nodes[4] = {
      findEdgeTail(dec, childMarkerEdge),
      findEdgeHead(dec, childMarkerEdge),
      findEdgeTail(dec, dec->members[member].markerToParent),
      findEdgeHead(dec, dec->members[member].markerToParent),
    };

    if ((nodes[0] == nodes[2] && nodes[1] == nodes[3]) || (nodes[0] == nodes[3] && nodes[1] == nodes[2]))
    {
      CMRdbgMsg(12, "Child marker edge %d is parallel to parent marker edge %d.\n",
        dec->members[childMember].markerOfParent, dec->members[member].markerToParent);

      if (dec->members[parentMember].type != DEC_MEMBER_TYPE_PARALLEL)
      {
        CMRdbgMsg(12, "Parent member is not a parallel member, so we create one for edges %d and %d.\n",
          dec->members[member].markerToParent, dec->members[member].markerOfParent);

        DEC_MEMBER newParallel;
        CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &newParallel) );
        DEC_EDGE markerOfParent = dec->members[member].markerOfParent;
        DEC_EDGE newMarkerOfParent, newMarkerToParent;
        CMR_CALL( createMarkerEdgePair(dec, parentMember, &newMarkerOfParent, dec->edges[markerOfParent].tail,
          dec->edges[markerOfParent].head, newParallel, &newMarkerToParent, -1, -1) );

        CMR_CALL( replaceEdgeInMembersEdgeList(dec, markerOfParent, newMarkerOfParent) );
        dec->edges[markerOfParent].childMember = member;
        dec->edges[markerOfParent].member = newParallel;
        dec->edges[markerOfParent].tail = -1;
        dec->edges[markerOfParent].head = -1;
        CMR_CALL( addEdgeToMembersEdgeList(dec, markerOfParent) );
        CMR_CALL( addEdgeToMembersEdgeList(dec, newMarkerToParent) );
        dec->members[member].parentMember = newParallel;
        parentMember = newParallel;
      }

      assert(dec->members[parentMember].type == DEC_MEMBER_TYPE_PARALLEL);

      CMRdbgMsg(12, "Removing child marker edge %d from member %d and adding it to the new parallel member %d.\n",
        childMarkerEdge, member, parentMember);

      CMR_CALL( removeEdgeFromMembersEdgeList(dec, childMarkerEdge) );
      dec->edges[childMarkerEdge].member = parentMember;
      CMR_CALL( addEdgeToMembersEdgeList(dec, childMarkerEdge) );
      dec->edges[childMarkerEdge].tail = -1;
      dec->edges[childMarkerEdge].head = -1;
      dec->members[childMember].parentMember = parentMember;

      debugDot(dec, NULL);
    }
  }

  CMR_CALL( parallelParentChildCheckMember(dec, parentMember, member) );

  return CMR_OKAY;
}

/**
 * \brief Ensures that no child marker in the reduced decomposition is parallel to the parent marker.
 */

static
CMR_ERROR parallelParentChildCheckReducedMembers(
  Dec* dec,           /**< Decomposition. */
  size_t* entryRows,  /**< Array of rows of new column's enries. */
  size_t numEntries   /**< Length of \p entryRows. */
)
{
  assert(dec);
  assert(entryRows);

  dec->parallelParentChildVisit++;

  /* Start processing at each member containing a row element. */
  for (size_t p = 0; p < numEntries; ++p)
  {
    size_t row = entryRows[p];
    if (row >= dec->numRows)
      continue;
    DEC_EDGE edge = dec->rowEdges[row].edge;
    if (edge == SIZE_MAX)
      continue;
    DEC_MEMBER member = findEdgeMember(dec, edge);
    CMRdbgMsg(6, "Entry %d is row %d of %d and corresponds to edge %d of member %d.\n", p, row, dec->numRows, edge,
      member);
    DEC_MEMBER parentMember = findMemberParent(dec, member);
    if (parentMember != SIZE_MAX)
      CMR_CALL( parallelParentChildCheckMember(dec, parentMember, member) );
  }

  return CMR_OKAY;
}

/**
 * \brief Creates, if necessary, the reduced member for \p member and calls itself for the parent.
 */

static
CMR_ERROR createReducedMembers(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  DEC_MEMBER member,                  /**< Member to create reduced member for. */
  ReducedMember** preducedMember      /**< Pointer for storing the created reduced member. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(member != SIZE_MAX);
  assert(preducedMember);

  CMRdbgMsg(8, "Attempting to create reduced member %d.\n", member);

  assert(member < newcolumn->memReducedMembers);
  ReducedMember* reducedMember = newcolumn->memberInfo[member].reducedMember;
  if (reducedMember)
  {
    /* This member is a known reduced member. If we meet an existing path of low depth, we remember
     * that. */

    CMRdbgMsg(10, "Reduced member exists.\n");

    if (!newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer ||
      reducedMember->depth < newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer->depth)
    {
      CMRdbgMsg(8, "Updating depth to %d.\n", reducedMember->depth);
      newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer = reducedMember;
    }
  }
  else
  {
    reducedMember = &newcolumn->reducedMembers[newcolumn->numReducedMembers];
    newcolumn->numReducedMembers++;
    newcolumn->memberInfo[member].reducedMember = reducedMember;
    assert(isRepresentativeMember(dec, member));
    reducedMember->member = member;
    reducedMember->numChildren = 0;
    reducedMember->rigidEndNodes[0] = -1;
    reducedMember->rigidEndNodes[1] = -1;
    reducedMember->rigidEndNodes[2] = -1;
    reducedMember->rigidEndNodes[3] = -1;

    CMRdbgMsg(10, "Reduced member is new.\n");

    DEC_MEMBER parentMember = findMemberParent(dec, member);
    if (parentMember != SIZE_MAX)
    {
      ReducedMember* parentReducedMember;
      CMR_CALL( createReducedMembers(dec, newcolumn, parentMember, &parentReducedMember) );
      reducedMember->parent = parentReducedMember;
      reducedMember->depth = parentReducedMember->depth + 1;
      reducedMember->rootMember = parentReducedMember->rootMember;
      reducedMember->type = 0;
      parentReducedMember->numChildren++;
    }
    else
    {
      reducedMember->parent = NULL;
      reducedMember->depth = 0;
      reducedMember->rootMember = member;
      reducedMember->type = 0;

      /* We found a new component. We temporarily store the root member and later compute the actual reduced root. */
      if (newcolumn->memReducedComponents == newcolumn->numReducedComponents)
      {
        newcolumn->memReducedComponents = 2 * newcolumn->memReducedComponents + 16;
        CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->reducedComponents, newcolumn->memReducedComponents) );
      }
      CMRdbgMsg(8, "Initializing the new reduced component %d.\n", newcolumn->numReducedComponents);

      newcolumn->reducedComponents[newcolumn->numReducedComponents].root = reducedMember;
      assert(isRepresentativeMember(dec, reducedMember->member));
      newcolumn->numReducedComponents++;
    }

    CMRdbgMsg(8, "The root member of %d is %d.\n", reducedMember->member, reducedMember->rootMember);
  }

  *preducedMember = reducedMember;

  return CMR_OKAY;
}

/**
 * \brief Computes the reduced decomposition.
 */

static
CMR_ERROR computeReducedDecomposition(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< newcolumn. */
  size_t* entryRows,        /**< Array of rows of new column's enries. */
  size_t numEntries         /**< Length of \p entryRows. */
)
{
  /* Identify all members on the path. For the induced sub-arborescence we also compute the
   * depths. After the computation, its root has depth pathRootDepth. */
  CMRdbgMsg(4, "Computing reduced t-decomposition.\n");

  /* Enlarge members array. */
  size_t maxRow = 0;
  for (size_t p = 0; p < numEntries; ++p)
  {
    if (entryRows[p] > maxRow)
      maxRow = entryRows[p];
  }
  size_t requiredNumReducedMembers = dec->numMembers + maxRow + 1;
  CMRdbgMsg(6, "Requiring storage for at most %ld reduced members.\n", requiredNumReducedMembers);
  if (requiredNumReducedMembers >= newcolumn->memReducedMembers)
  {
    size_t newSize = 16 + 2 * newcolumn->memReducedMembers;
    while (newSize < requiredNumReducedMembers)
      newSize *= 2;
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->reducedMembers, newSize) );
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->memberInfo, newSize) );
    for (size_t m = newcolumn->memReducedMembers; m < newSize; ++m)
    {
      newcolumn->memberInfo[m].reducedMember = NULL;
      newcolumn->memberInfo[m].rootDepthMinimizer = NULL;
    }
    newcolumn->memReducedMembers = newSize;
  }

  newcolumn->numReducedMembers = 0;
  for (size_t p = 0; p < numEntries; ++p)
  {
    size_t row = entryRows[p];
    DEC_EDGE edge = (row < dec->numRows) ? dec->rowEdges[row].edge : SIZE_MAX;
    CMRdbgMsg(6, "Entry %d is row %d of %d and corresponds to edge %d.\n", p, row, dec->numRows, edge);
    if (edge != SIZE_MAX)
    {
      DEC_MEMBER member = findEdgeMember(dec, edge);
      CMRdbgMsg(8, "Edge %d exists and belongs to member %d.\n", edge, member);
      ReducedMember* reducedMember;
      CMR_CALL( createReducedMembers(dec, newcolumn, member, &reducedMember) );

      /* For the first edge of this member, we set the depth minimizer to the new reduced member. */
      if (!newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer)
      {
        newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer = reducedMember;
      }
    }
  }

  /* We set the reduced roots according to the minimizers. */
  for (size_t i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    CMRdbgMsg(6, "Considering reduced component %d with initial root member %d.\n", i,
      newcolumn->reducedComponents[i].root->member);
    CMRdbgMsg(8, "The minimizer is %d.\n",
      newcolumn->memberInfo[newcolumn->reducedComponents[i].root->member].rootDepthMinimizer->member);

    newcolumn->reducedComponents[i].rootDepth =
      newcolumn->memberInfo[newcolumn->reducedComponents[i].root->member].rootDepthMinimizer->depth;
    newcolumn->reducedComponents[i].root =
      newcolumn->memberInfo[newcolumn->reducedComponents[i].root->member].rootDepthMinimizer;
    assert(isRepresentativeMember(dec, newcolumn->reducedComponents[i].root->member));
    newcolumn->reducedComponents[i].numTerminals = 0;
    CMRdbgMsg(8, "Member %d is the new root of the reduced decomposition of this component.\n",
      newcolumn->reducedComponents[i].root->member);
  }

  /* Allocate memory for children. */
  if (newcolumn->memChildrenStorage < newcolumn->numReducedMembers)
  {
    newcolumn->memChildrenStorage = 2*newcolumn->numReducedMembers;
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->childrenStorage, newcolumn->memChildrenStorage) );
  }
  newcolumn->usedChildrenStorage = 0;

  /* Set children memory pointer of each reduced member. */
  for (size_t m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[m];
    if (reducedMember->depth >= newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer->depth)
    {
      CMRdbgMsg(6, "Member %d's depth is greater than that of the root, and it has %d children.\n",
        reducedMember->member, reducedMember->numChildren);
      reducedMember->children = &newcolumn->childrenStorage[newcolumn->usedChildrenStorage];
      newcolumn->usedChildrenStorage += reducedMember->numChildren;
      reducedMember->numChildren = 0;
    }
    else
    {
      CMRdbgMsg(6, "Member %d's depth is smaller than that of its new root.\n", reducedMember->member);
      continue;
    }
  }

  CMRdbgMsg(4, "Total number of children is %d with space for %d.\n", newcolumn->usedChildrenStorage, newcolumn->memChildrenStorage);

  /* Set children of each reduced member. */
  for (size_t m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[m];
    if (reducedMember->depth <= newcolumn->memberInfo[reducedMember->rootMember].rootDepthMinimizer->depth)
    {
      CMRdbgMsg(6, "Member %d's depth is smaller than or equal to that of its reduced root.\n", reducedMember->member);
      continue;
    }

    DEC_MEMBER parentMember = findMemberParent(dec, newcolumn->reducedMembers[m].member);
    ReducedMember* parentReducedMember = parentMember != SIZE_MAX ? newcolumn->memberInfo[parentMember].reducedMember : NULL;
    CMRdbgMsg(6, "Member %d's depth is greater than that of its reduced root. Its parent is %d, and reduced parent %p.\n",
      reducedMember->member, parentMember, parentReducedMember);

    if (parentReducedMember)
    {
      CMRdbgMsg(6, "Reduced member %ld (= member %d) has %d (= member %d) as child %d.\n",
        (parentReducedMember - newcolumn->reducedMembers),
        parentReducedMember->member, m, newcolumn->reducedMembers[m].member, parentReducedMember->numChildren);
      parentReducedMember->children[parentReducedMember->numChildren] = &newcolumn->reducedMembers[m];
      parentReducedMember->numChildren++;
    }
  }

  /* Clean-up rootDepthMinimizer entries. */
  for (size_t i = 0; i < newcolumn->numReducedMembers; ++i)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[i];
    assert(reducedMember);
    DEC_MEMBER rootMember = reducedMember->rootMember;
    assert(rootMember < dec->memMembers);
    newcolumn->memberInfo[rootMember].rootDepthMinimizer = NULL;
  }

  return CMR_OKAY;
}

/**
 * \brief Creates a path edge structure.
 *
 * Adds it to the singly-linked path edge list of \p reducedMember and to the overall singly-linked path edge list.
 * The created \ref PathEdge can be accessed via \ref ReducedMember::firstPathEdge.
 */

static
CMR_ERROR createPathEdge(
  Dec* dec,                     /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,     /**< new column. */
  DEC_EDGE edge,                /**< The edge it refers to. */
  ReducedMember* reducedMember  /**< Reduced member this path edge belongs to. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(edge < newcolumn->memEdgesInPath);

  assert(newcolumn->numPathEdges < newcolumn->memPathEdges);
  PathEdge* pathEdge = &newcolumn->pathEdges[newcolumn->numPathEdges];
  pathEdge->edge = edge;
  pathEdge->nextSibling = reducedMember->firstPathEdge;
  reducedMember->firstPathEdge = pathEdge;
  pathEdge->nextOverall = newcolumn->firstPathEdge;
  newcolumn->firstPathEdge = pathEdge;
  newcolumn->numPathEdges++;

  newcolumn->edgesInPath[edge] = true;

  /* Increase node degrees of nodes in a rigid parent. */
  if (dec->members[reducedMember->member].type == DEC_MEMBER_TYPE_RIGID)
  {
    CMRdbgMsg(14, "Creating path edge for %d = {%d,%d} in rigid member %d.\n", edge, findEdgeTail(dec, edge),
      findEdgeHead(dec, edge), reducedMember->member);
    DEC_EDGE tail = findEdgeTail(dec, edge);
    assert(tail < newcolumn->memNodesDegree);
    newcolumn->nodesDegree[tail]++;
    CMRdbgMsg(16, "Increasing node degree of %d to %d.\n", tail, newcolumn->nodesDegree[tail]);
    DEC_EDGE head = findEdgeHead(dec, edge);
    assert(head < newcolumn->memNodesDegree);
    newcolumn->nodesDegree[head]++;
    CMRdbgMsg(16, "Increasing node degree of %d to %d.\n", head, newcolumn->nodesDegree[head]);
  }
  else
    CMRdbgMsg(14, "Creating path edge for %d in non-rigid member %d.\n", edge, reducedMember->member);

  return CMR_OKAY;
}

/**
 * \brief Creates members and reduced members of new edges.
 */

static
CMR_ERROR completeReducedDecomposition(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< newcolumn. */
  size_t* rows,             /**< Array of rows (of new column's entries). */
  size_t numRows            /**< Length of \p rows. */
)
{
  assert(dec);
  assert(newcolumn);

  /* Check if we need new rows. */

  size_t newNumRows = dec->numRows;
  for (size_t p = 0; p < numRows; ++p)
  {
    size_t row = rows[p];
    DEC_EDGE edge = (row < dec->numRows) ? dec->rowEdges[row].edge : SIZE_MAX;
    if (edge == SIZE_MAX)
    {
      if (row >= newNumRows)
        newNumRows = row + 1;
    }
  }

  CMRdbgMsg(4, "Completing reduced decomposition: increasing #rows from %d to %d.\n", dec->numRows, newNumRows);

  /* Create single-edge parallel members for all new rows. */

  if (newNumRows > dec->numRows)
  {
    if (newNumRows > dec->memRows)
    {
      dec->memRows = 2*newNumRows;
      CMR_CALL( CMRreallocBlockArray(dec->cmr, &dec->rowEdges, dec->memRows) );
    }

    for (size_t r = dec->numRows; r < newNumRows; ++r)
    {
      DEC_MEMBER member;
      CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &member) );

      DEC_EDGE edge;
      CMR_CALL( createEdge(dec, member, &edge) );
      CMR_CALL( addEdgeToMembersEdgeList(dec, edge) );
      dec->edges[edge].element = CMRrowToElement(r);
      dec->edges[edge].head = -1;
      dec->edges[edge].tail = -1;
      dec->edges[edge].childMember = -1;

      CMRdbgMsg(8, "New row %d is edge %d of member %d.\n", r, edge, member);

      dec->rowEdges[r].edge = edge;
    }
  }

  /* Create reduced members. */
  for (size_t p = 0; p < numRows; ++p)
  {
    size_t row = rows[p];
    if (row >= dec->numRows)
    {
      /* Edge and member for this row were created. We create the reduced member. */
      DEC_EDGE edge = dec->rowEdges[row].edge;
      DEC_MEMBER member = findEdgeMember(dec, edge);

      CMRdbgMsg(4, "Creating reduced member %ld for edge %d of member %d.\n", newcolumn->numReducedMembers, edge,
        member);

      assert(newcolumn->numReducedMembers < newcolumn->memReducedMembers);
      ReducedMember* reducedMember = &newcolumn->reducedMembers[newcolumn->numReducedMembers];
      newcolumn->numReducedMembers++;
      reducedMember->numChildren = 0;
      assert(isRepresentativeMember(dec, member));
      reducedMember->member = member;
      reducedMember->depth = 0;
      reducedMember->rootMember = -1;
      reducedMember->type = TYPE_ROOT;
      reducedMember->firstPathEdge = NULL;
      CMR_CALL( createPathEdge(dec, newcolumn, edge, reducedMember) );

      if (newcolumn->numReducedComponents == newcolumn->memReducedComponents)
      {
        newcolumn->memReducedComponents = 2 * newcolumn->memReducedComponents + 16;
        CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->reducedComponents,
          newcolumn->memReducedComponents) );
      }

      newcolumn->reducedComponents[newcolumn->numReducedComponents].root = reducedMember;
      assert(isRepresentativeMember(dec, reducedMember->member));
      newcolumn->reducedComponents[newcolumn->numReducedComponents].rootDepth = 0;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].numTerminals = 0;
#if !defined(NDEBUG)
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalMember[0] = INT_MIN;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalMember[1] = INT_MIN;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalNode[0] = INT_MIN;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalNode[1] = INT_MIN;
#endif /* !NDEBUG */
      newcolumn->numReducedComponents++;
    }
  }

  dec->numRows = newNumRows;

  return CMR_OKAY;
}

/**
 * \brief Creates path edges in reduced decomposition.
 */

static
CMR_ERROR createReducedDecompositionPathEdges(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< newcolumn. */
  size_t* rows,             /**< Array of rows (of new column's enries). */
  size_t numRows            /**< Length of \p rows. */
)
{
  assert(dec);
  assert(newcolumn);

  CMRdbgMsg(4, "Initializing edge lists for members of reduced decomposition.\n");

  assert(newcolumn->numPathEdges == 0);
  size_t maxNumPathEdges = numRows + dec->numMembers;
  if (newcolumn->memPathEdges < maxNumPathEdges)
  {
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->pathEdges, maxNumPathEdges) );
    newcolumn->memPathEdges = maxNumPathEdges;
  }

  /* Start with empty lists. */
  for (size_t i = 0; i < newcolumn->numReducedMembers; ++i)
    newcolumn->reducedMembers[i].firstPathEdge = NULL;

  /* Fill edge lists. */
  for (size_t p = 0; p < numRows; ++p)
  {
    size_t row = rows[p];
    DEC_EDGE edge = (row < dec->numRows) ? dec->rowEdges[row].edge : SIZE_MAX;
    if (edge != SIZE_MAX)
    {
      DEC_MEMBER member = findEdgeMember(dec, edge);
      assert(member != SIZE_MAX);
      ReducedMember* reducedMember = newcolumn->memberInfo[member].reducedMember;

      assert(reducedMember);
      CMR_CALL( createPathEdge(dec, newcolumn, edge, reducedMember) );

      CMRdbgMsg(6, "Edge %d <%s> belongs to reduced member %ld which is %s member %d.\n", edge,
        CMRelementString(dec->edges[edge].element, NULL), (reducedMember - newcolumn->reducedMembers),
        memberTypeString(reducedMember->member), reducedMember->member);
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Count the number of children of a reduced member having certain types.
 */

static
CMR_ERROR countChildrenTypes(
  Dec* dec,                     /**< Decomposition. */
  ReducedMember* reducedMember, /**< Reduced member. */
  int* pNumOneEnd,              /**< Number of children that (recursively) must contain one path end. */
  int* pNumTwoEnds,             /**< Number of children that (recursively) must contain two path ends. */
  DEC_EDGE childMarkerEdges[2]  /**< Array for storing a child marker edges containing one/two path ends. */
)
{
  assert(dec);
  assert(reducedMember);
  assert(pNumOneEnd);
  assert(pNumTwoEnds);
  assert(childMarkerEdges);

  *pNumOneEnd = 0;
  *pNumTwoEnds = 0;
  childMarkerEdges[0] = -1;
  childMarkerEdges[1] = -1;
  int nextChildMarker = 0;

  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    assert(child);

    if (child->type == TYPE_SINGLE_CHILD)
    {
      if (nextChildMarker < 2)
      {
        childMarkerEdges[nextChildMarker] = dec->members[findMember(dec, child->member)].markerOfParent;
        nextChildMarker++;
      }
      (*pNumOneEnd)++;
    }
    else if (child->type == TYPE_DOUBLE_CHILD)
    {
      if (nextChildMarker < 2)
      {
        childMarkerEdges[nextChildMarker] = dec->members[findMember(dec, child->member)].markerOfParent;
        nextChildMarker++;
      }
      (*pNumTwoEnds)++;
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Determines the type of a parallel member.
 */

static
CMR_ERROR determineTypeParallel(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  int depth                           /**< Depth of member in reduced decomposition. */
)
{
  CMR_UNUSED(dec);

  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(dec->members[findMember(dec, reducedMember->member)].type == DEC_MEMBER_TYPE_PARALLEL);

  if (depth == 0)
  {
    /* A parallel root always works. */
    reducedMember->type = TYPE_ROOT;
    return CMR_OKAY;
  }

  /* No children, but a path edge. */
  if (2*numTwoEnds + numOneEnd == 0 && reducedMember->firstPathEdge)
    reducedMember->type = TYPE_CYCLE_CHILD;
  else if (numOneEnd == 1)
    reducedMember->type = TYPE_SINGLE_CHILD;
  else if (numOneEnd + 2 * numTwoEnds == 2)
  {
    if (reducedMember->firstPathEdge)
    {
      /* Tested in TypingInnerParallelWithEdge */
      newcolumn->remainsGraphic = false;
    }
    else
      reducedMember->type = TYPE_DOUBLE_CHILD;
  }
  else
  {
    /* Since no child contains path edges, the parallel must be a leaf of the reduced decomposition and contains one. */
    assert(reducedMember->firstPathEdge);
    reducedMember->type = TYPE_CYCLE_CHILD;
  }

  return CMR_OKAY;
}

/**
 * \brief Determines the type of a series member.
 */

static
CMR_ERROR determineTypeSeries(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(dec->members[findMember(dec, reducedMember->member)].type == DEC_MEMBER_TYPE_SERIES);

  DEC_MEMBER member = findMember(dec, reducedMember->member);

  int countPathEdges = 0;
  for (PathEdge* edge = reducedMember->firstPathEdge; edge != NULL; edge = edge->nextSibling)
    ++countPathEdges;
  int numEdges = dec->members[member].numEdges;

  CMRdbgMsg(8+2*depth,
    "Determining type of series with %d edges, %d path edges, %d 1-end children and %d 2-end children.\n", numEdges,
    countPathEdges, numOneEnd, numTwoEnds);

  if (depth == 0)
  {
    /* We assume that we are not the root of the whole decomposition. */
    assert(dec->members[member].parentMember != SIZE_MAX);

    /* Tested in TypingRootSeriesDoubleChild */
    newcolumn->remainsGraphic = (numTwoEnds == 0);
    reducedMember->type = (countPathEdges == numEdges - 1) ? TYPE_CYCLE_CHILD : TYPE_ROOT;
    return CMR_OKAY;
  }

  if (countPathEdges == numEdges - 1)
  {
    reducedMember->type = TYPE_CYCLE_CHILD;
  }
  else if (countPathEdges + numTwoEnds == numEdges - 1)
  {
    assert(numTwoEnds == 1);
    reducedMember->type = TYPE_DOUBLE_CHILD;
  }
  else if (numTwoEnds == 1)
  {
    /* Tested in TypingInnerSeries2Child. */
    newcolumn->remainsGraphic = false;
  }
  else if (numOneEnd == 1)
  {
    reducedMember->type = TYPE_SINGLE_CHILD;
  }
  else if (numOneEnd == 2)
  {
    reducedMember->type = TYPE_DOUBLE_CHILD;
  }
  else
  {
    assert(numOneEnd == 0);
    assert(numTwoEnds == 0);
    reducedMember->type = TYPE_SINGLE_CHILD;
  }

  return CMR_OKAY;
}

/**
 * \brief Determines the type of a rigid member.
 */

static
CMR_ERROR determineTypeRigid(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  DEC_EDGE childMarkerEdges[2],       /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(dec->members[findMember(dec, reducedMember->member)].type == DEC_MEMBER_TYPE_RIGID);

  DEC_MEMBER member = findMember(dec, reducedMember->member);

  /* Collect nodes of parent marker and of child markers containing one/two path end nodes. */
  DEC_NODE parentMarkerNodes[2] = {
    depth == 0 ? SIZE_MAX : findEdgeTail(dec, dec->members[member].markerToParent),
    depth == 0 ? SIZE_MAX : findEdgeHead(dec, dec->members[member].markerToParent)
  };
  DEC_NODE childMarkerNodes[4] = {
    childMarkerEdges[0] == SIZE_MAX ? SIZE_MAX : findEdgeTail(dec, childMarkerEdges[0]),
    childMarkerEdges[0] == SIZE_MAX ? SIZE_MAX : findEdgeHead(dec, childMarkerEdges[0]),
    childMarkerEdges[1] == SIZE_MAX ? SIZE_MAX : findEdgeTail(dec, childMarkerEdges[1]),
    childMarkerEdges[1] == SIZE_MAX ? SIZE_MAX : findEdgeHead(dec, childMarkerEdges[1])
  };

  DEC_NODE* pathEndNodes = reducedMember->rigidEndNodes;
  int numPathEndNodes = 0;

  /* Check the node degrees (with respect to path edges) in this component. */
  for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->nextSibling)
  {
    DEC_NODE nodes[2] = { findEdgeHead(dec, reducedEdge->edge), findEdgeTail(dec, reducedEdge->edge) };
    CMRdbgMsg(6 + 2*depth, "Reduced edge {%d,%d}\n", nodes[0], nodes[1]);
    for (int i = 0; i < 2; ++i)
    {
      DEC_NODE v = nodes[i];
      CMRdbgMsg(6 + 2*depth, "Node %d of rigid member %d has path-degree %d.\n", v, reducedMember->member, newcolumn->nodesDegree[v]);
      if (newcolumn->nodesDegree[v] >= 3)
      {
        CMRdbgMsg(6 + 2*depth, "Node %d of rigid member %d has path-degree at least 3.\n", v, reducedMember->member);

        /* Tested in TypingAnyRigidDegree3. */
        newcolumn->remainsGraphic = false;
        return CMR_OKAY;
      }

      if (newcolumn->nodesDegree[v] == 1)
      {
        if (numPathEndNodes == 4)
        {
          CMRdbgMsg(6 + 2*depth, "Rigid member %d has at least five path end nodes: %d, %d, %d, %d and %d.\n",
            reducedMember->member, pathEndNodes[0], pathEndNodes[1], pathEndNodes[2], pathEndNodes[3], v);

          /* Tested in TypingAnyRigidManyPaths.*/
          newcolumn->remainsGraphic = false;
          return CMR_OKAY;
        }

        pathEndNodes[numPathEndNodes] = v;
        ++numPathEndNodes;
      }
    }
  }

  /* Exchange end nodes array using these rules:
   * If there are 4 edges, then one path connects 0 with 1 and the other 2 with 3.
   * For each path, parent marker nodes should come first (i.e., index 0 and 2).
   */

  if (numPathEndNodes == 4)
  {
    DEC_EDGE* nodeEdges = NULL;
    CMR_CALL(CMRallocStackArray(dec->cmr, &nodeEdges, 2*dec->memNodes));

    /* Initialize relevant entries to -1. */
    for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->nextSibling)
    {
      DEC_NODE nodes[2] = { findEdgeHead(dec, reducedEdge->edge), findEdgeTail(dec, reducedEdge->edge) };
      for (int i = 0; i < 2; ++i)
      {
        DEC_NODE v = nodes[i];
        nodeEdges[2*v] = -1;
        nodeEdges[2*v + 1] = -1;
      }
    }

    /* Store incident edges for every node. */
    for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->nextSibling)
    {
      DEC_NODE nodes[2] = { findEdgeHead(dec, reducedEdge->edge), findEdgeTail(dec, reducedEdge->edge) };
      for (int i = 0; i < 2; ++i)
      {
        DEC_NODE v = nodes[i];
        nodeEdges[2*v + (nodeEdges[2*v] == SIZE_MAX ? 0 : 1)] = reducedEdge->edge;
      }
    }

    /* Start at end node 0 and see where we end. */
    DEC_EDGE previousEdge = -1;
    DEC_NODE currentNode = pathEndNodes[0];
    while (true)
    {
      DEC_EDGE edge = nodeEdges[2*currentNode];
      if (edge == previousEdge)
        edge = nodeEdges[2*currentNode+1];
      if (edge == SIZE_MAX)
        break;
      previousEdge = edge;
      DEC_NODE v = findEdgeHead(dec, edge);
      currentNode = (v != currentNode) ? v : findEdgeTail(dec, edge);
    }
    CMRfreeStackArray(dec->cmr, &nodeEdges);

    /* Exchange such that we end nodes 0 and 1 are end nodes of the same path. */
    if (currentNode == pathEndNodes[2])
    {
      pathEndNodes[2] = pathEndNodes[1];
      pathEndNodes[1] = currentNode;
    }
    else if (currentNode == pathEndNodes[3])
    {
      pathEndNodes[3] = pathEndNodes[1];
      pathEndNodes[1] = currentNode;
    }

    /* Exchange such that end node 2 is at the parent marker. */
    if (pathEndNodes[2] != parentMarkerNodes[0] && pathEndNodes[2] != parentMarkerNodes[1])
      SWAP_INTS(pathEndNodes[2], pathEndNodes[3]);
  }

  /* Exchange such that end node 0 is at the parent marker. */
  if (numPathEndNodes >= 2 && pathEndNodes[0] != parentMarkerNodes[0] && pathEndNodes[0] != parentMarkerNodes[1])
  {
    SWAP_INTS(pathEndNodes[0], pathEndNodes[1]);
  }

  CMRdbgMsg(6 + 2*depth, "Rigid %smember %d has %d path end nodes:", depth == 0 ? "root " : "", member, numPathEndNodes);
  if (numPathEndNodes >= 2)
    CMRdbgMsg(0, " a path from %d to %d.", pathEndNodes[0], pathEndNodes[1]);
  if (numPathEndNodes == 4)
    CMRdbgMsg(0, " a path from %d to %d.", pathEndNodes[2], pathEndNodes[3]);
  CMRdbgMsg(0, "\n");

  if (depth == 0)
  {
    if (numPathEndNodes == 0)
    {
      CMRdbgMsg(6 + 2*depth, "Root without paths. Child markers are {%d,%d} and {%d,%d}.\n", childMarkerNodes[0], childMarkerNodes[1],
        childMarkerNodes[2], childMarkerNodes[3]);
      /* No path edges, so there should be two adjacent child marker edges. */
      if (numOneEnd == 2 && (childMarkerNodes[0] == childMarkerNodes[2] || childMarkerNodes[0] == childMarkerNodes[3]
          || childMarkerNodes[1] == childMarkerNodes[2] || childMarkerNodes[1] == childMarkerNodes[3]))
      {
        reducedMember->type = TYPE_ROOT;
      }
      else
      {
        /* Tested in TypingRootRigidNoPathsDisjointSingleChildren. */
        newcolumn->remainsGraphic = false;
      }
    }
    else if (numPathEndNodes == 2)
    {
      if (numOneEnd == 1)
      {
        bool pathAdjacentToChildMarker = false;
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 2; ++j)
          {
            if (pathEndNodes[i] == childMarkerNodes[j])
            {
              pathAdjacentToChildMarker = true;
              break;
            }
          }
        }
        if (pathAdjacentToChildMarker)
          reducedMember->type = TYPE_ROOT;
        else
        {
          /* Tested in TypingRootRigidOnePathSingleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numOneEnd == 2)
      {
        /* 1 path and 2 child markers, each containing one end. */
        bool childMarkerNodesMatched[2] = { false, false };
        bool endNodesMatched[2] = { false, false };
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 4; ++j)
          {
            if (reducedMember->rigidEndNodes[i] == childMarkerNodes[j])
            {
              endNodesMatched[i] = true;
              childMarkerNodesMatched[j/2] = true;
            }
          }
        }

        if (childMarkerNodesMatched[0] && childMarkerNodesMatched[1] && endNodesMatched[0] && endNodesMatched[1])
          reducedMember->type = TYPE_ROOT;
        else
        {
          /* Tested in TypingRootRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        assert(numOneEnd == 0);
        reducedMember->type = TYPE_ROOT;
      }
      else
      {
        assert(numOneEnd == 0);
        assert(numTwoEnds == 1);

        if ((childMarkerNodes[0] == pathEndNodes[0] && childMarkerNodes[1] == pathEndNodes[1])
          || (childMarkerNodes[0] == pathEndNodes[1] && childMarkerNodes[1] == pathEndNodes[0]))
        {
          reducedMember->type = TYPE_ROOT;
        }
        else
        {
          /* Tested in TypingRootRigidOnePathDoubleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else
    {
      assert(numPathEndNodes == 4);

      /* Tested in TypingRootRigidTwoPaths. */
      newcolumn->remainsGraphic = false;
    }
  }
  else
  {
    /* Non-root rigid member. */

    int parentMarkerDegrees[2] = {
      newcolumn->nodesDegree[parentMarkerNodes[0]],
      newcolumn->nodesDegree[parentMarkerNodes[1]]
    };

    if (numPathEndNodes == 0)
    {
      /* We have no path edges, so there must be at least one child containing one/two path ends. */
      assert(numOneEnd + numTwoEnds > 0);
      /* We should not have a child marker edge parallel to the parent marker edge! */
      assert(!(parentMarkerNodes[0] == childMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1])
        && !(parentMarkerNodes[0] == childMarkerNodes[1] && parentMarkerNodes[1] == childMarkerNodes[0]));

      if (numOneEnd == 0)
      {
        /* Even if parent and child marker (type 4) are adjacent, this is non-graphic. */

        /* Tested in TypingInnerRigidNoPathNoSingleChild. */
        newcolumn->remainsGraphic = false;
      }
      else if (numOneEnd == 1)
      {
        if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
          || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
        {
          reducedMember->type = TYPE_SINGLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidNoPathOneSingleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
      else
      {
        /* We're not the root but have two proper children. */
        int childMarkerParentNode[2] = { -1, -1 };
        bool isParallel = false;
        for (int i = 0; i < 4; ++i)
        {
          for (int j = 0; j < 2; ++j)
          {
            if (childMarkerNodes[i] == parentMarkerNodes[j])
            {
              if (childMarkerParentNode[i/2] >= 0)
                isParallel = true;
              childMarkerParentNode[i/2] = j;
            }
          }
        }

        CMRdbgMsg(6 + 2*depth, "No paths, %s, child[0] incident to %d and child[1] incident to %d.\n",
          isParallel ? "parallel" : "not parallel", childMarkerParentNode[0], childMarkerParentNode[1]);

        if (!isParallel && childMarkerParentNode[0] != -1 && childMarkerParentNode[1] != -1
          && childMarkerParentNode[0] != childMarkerParentNode[1])
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidNoPathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else if (numPathEndNodes == 2)
    {
      /* Exchange such that end node 0 is at the parent marker (or none of the end nodes is). */
      if (pathEndNodes[0] != parentMarkerNodes[0] && pathEndNodes[0] != parentMarkerNodes[1])
        SWAP_INTS(pathEndNodes[0], pathEndNodes[0]);

      if (numOneEnd == 1)
      {
        CMRdbgMsg(6 + 2*depth, "%d-%d-path, parent {%d,%d} and child {%d,%d}\n", pathEndNodes[0], pathEndNodes[1],
          parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0], childMarkerNodes[1]);

        if (parentMarkerNodes[0] != pathEndNodes[0])
        {
          SWAP_INTS(parentMarkerNodes[0], parentMarkerNodes[1]);
          SWAP_INTS(parentMarkerDegrees[0], parentMarkerDegrees[1]);
        }
        if (parentMarkerNodes[0] != pathEndNodes[0])
        {
          /* Tested in TypingInnerRigidOnePathOneSingleChild. */
          newcolumn->remainsGraphic = false;
          return CMR_OKAY;
        }

        if (parentMarkerNodes[1] == pathEndNodes[1])
        {
          /* Path closes a cycle with parent marker edge. */
          if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
            || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
          {
            reducedMember->type = TYPE_SINGLE_CHILD;
          }
          else
          {
            /* Tested in TypingInnerRigidOnePathOneSingleChild. */
            newcolumn->remainsGraphic = false;
          }
        }
        else
        {
          /* Path end is not incident to parent marker edge. */
          if (childMarkerNodes[0] == pathEndNodes[1] || childMarkerNodes[1] == pathEndNodes[1])
          {
            reducedMember->type = TYPE_SINGLE_CHILD;
          }
          else if (childMarkerNodes[0] == parentMarkerNodes[1] || childMarkerNodes[1] == parentMarkerNodes[1])
          {
            reducedMember->type = TYPE_DOUBLE_CHILD;
          }
          else
          {
            /* Tested in TypingInnerRigidOnePathOneSingleChild. */
            newcolumn->remainsGraphic = false;
          }
        }
      }
      else if (numOneEnd == 2)
      {
        CMRdbgMsg(6 + 2*depth, "%d-%d-path, parent marker {%d,%d} and child markers {%d,%d} and {%d,%d}\n",
          pathEndNodes[0], pathEndNodes[1], parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0],
          childMarkerNodes[1], childMarkerNodes[2], childMarkerNodes[3]);

        /* We know that pathEndNodes[0] is incident to the parent marker. */
        DEC_NODE otherParentNode;
        if (pathEndNodes[0] == parentMarkerNodes[0])
          otherParentNode = parentMarkerNodes[1];
        else if (pathEndNodes[0] == parentMarkerNodes[1])
          otherParentNode = parentMarkerNodes[0];
        else
        {
          /* Tested in TypingInnerRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
          return CMR_OKAY;
        }

        /* Two 1-ends and a path that closes a cycle with the parent marker is only allowed for root members. */
        if (pathEndNodes[1] == otherParentNode)
        {
          /* Tested in TypingInnerRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
          return CMR_OKAY;
        }

        bool childMatched[2] = { false, false };
        bool pathEndMatched = false;
        bool otherParentMatched = false;
        for (int i = 0; i < 4; ++i)
        {
          if (childMarkerNodes[i] == pathEndNodes[1])
          {
            childMatched[i/2] = true;
            pathEndMatched = true;
          }
          if (childMarkerNodes[i] == otherParentNode)
          {
            childMatched[i/2] = true;
            otherParentMatched = true;
          }
        }

        if (childMatched[0] && childMatched[1] && pathEndMatched && otherParentMatched)
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        if (parentMarkerDegrees[0] % 2 == 0 && parentMarkerDegrees[1] == 1)
        {
          reducedMember->type = TYPE_SINGLE_CHILD;
        }
        else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] % 2 == 0)
        {
          reducedMember->type = TYPE_SINGLE_CHILD;
        }
        else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] == 1)
        {
          reducedMember->type = TYPE_CYCLE_CHILD;
        }
        else
        {
          /* Both have degree 0 or 2. */
          /* Tested in TypingInnerRigidOnePathNoChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else
      {
        assert(numTwoEnds == 1);

        if ((pathEndNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[0]
          && childMarkerNodes[1] == pathEndNodes[1])
          || (pathEndNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1]
          && childMarkerNodes[0] == pathEndNodes[1])
          || (pathEndNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[0]
          && childMarkerNodes[1] == pathEndNodes[1])
          || (pathEndNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[1]
          && childMarkerNodes[0] == pathEndNodes[1]))
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidOnePathDoubleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else if (numPathEndNodes == 4)
    {
      if (reducedMember->rigidEndNodes[0] != parentMarkerNodes[0] && reducedMember->rigidEndNodes[0] != parentMarkerNodes[1])
      {
        CMRdbgMsg(6 + 2*depth, "First path does not start at parent marker edge.\n");
        /* Tested in TypingInnerRigidTwoPathsNonadjacentParent. */
        newcolumn->remainsGraphic = false;
        return CMR_OKAY;
      }
      if (reducedMember->rigidEndNodes[2] != parentMarkerNodes[0] && reducedMember->rigidEndNodes[2] != parentMarkerNodes[1])
      {
        CMRdbgMsg(6 + 2*depth, "Second path does not start at parent marker edge.\n");
        /* Tested in TypingInnerRigidTwoPathsNonadjacentParent. */
        newcolumn->remainsGraphic = false;
        return CMR_OKAY;
      }

      if (numOneEnd == 1)
      {
        bool pathConnects[2] = {
          reducedMember->rigidEndNodes[1] == childMarkerNodes[0]
            || reducedMember->rigidEndNodes[1] == childMarkerNodes[1],
          reducedMember->rigidEndNodes[3] == childMarkerNodes[0]
            || reducedMember->rigidEndNodes[3] == childMarkerNodes[1]
        };

        if (pathConnects[0] && pathConnects[1])
        {
          CMRdbgMsg(6 + 2*depth, "The both paths end at the child marker edge.\n");
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else if (pathConnects[0])
        {
          CMRdbgMsg(6 + 2*depth, "Only the first path ends at the child marker edge.\n");
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else if (pathConnects[1])
        {
          CMRdbgMsg(6 + 2*depth, "Only the second path ends at the child marker edge.\n");
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          CMRdbgMsg(6 + 2*depth, "No path ends at the child marker edge.\n");
          /* Tested in TypingInnerRigidTwoPathOneSingleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numOneEnd == 2)
      {
        bool pathConnected[2] = { false, false };
        bool childConnected[2] = { false, false };
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 4; ++j)
          {
            if (reducedMember->rigidEndNodes[1 + 2*i] == childMarkerNodes[j])
            {
              pathConnected[i] = true;
              childConnected[j/2] = true;
            }
          }
        }
        if (pathConnected[0] && pathConnected[1] && childConnected[0] && childConnected[1])
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          CMRdbgMsg(6 + 2*depth, "No pairing of paths to nodes of child marker edges possible.\n");
          /* Tested in TypingInnerRigidTwoPathsTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        reducedMember->type = TYPE_DOUBLE_CHILD;
      }
      else
      {
        /* One child marker with both ends. We already checked that path end nodes 0 and 2 are parent marker nodes. */
        assert(numTwoEnds == 1);

        if ((pathEndNodes[1] == childMarkerNodes[0] && pathEndNodes[3] == childMarkerNodes[1])
          || (pathEndNodes[1] == childMarkerNodes[1] && pathEndNodes[3] == childMarkerNodes[0]))
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          CMRdbgMsg(6 + 2*depth, "Paths do not end at child marker.");
          /* Tested in TypingInnerRigidTwoPathsDoubleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Determines the type of a member and all its children in the reduced decomposition.
 */

static
CMR_ERROR determineTypes(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);

  CMRdbgMsg(4 + 2*depth, "determineTypes(member %d = reduced member %ld)\n", reducedMember->member,
    reducedMember - &newcolumn->reducedMembers[0]);

  /* First handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* reducedChild = reducedMember->children[c];
    assert(reducedChild);
    CMR_CALL( determineTypes(dec, newcolumn, reducedComponent, reducedChild, depth + 1) );

    if (newcolumn->remainsGraphic)
    {
      CMRdbgMsg(6 + 2*depth, "Child member %d of %d has type %d\n", reducedMember->children[c]->member, reducedMember->member,
        reducedMember->children[c]->type);
    }
    else
    {
      CMRdbgMsg(6 + 2*depth, "Child prohibits graphicness.\n");
      return CMR_OKAY;
    }
  }

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  CMR_CALL( countChildrenTypes(dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

#if defined(CMR_DEBUG)
  int countReducedEdges = 0;
  for (PathEdge* e = reducedMember->firstPathEdge; e; e = e->nextSibling)
    ++countReducedEdges;
  CMRdbgMsg(6 + 2*depth, "Member %d has %d children with one end and %d with two ends and %d path edges.\n",
    reducedMember->member, numOneEnd, numTwoEnds, countReducedEdges);
#endif /* CMR_DEBUG */

  if (2*numTwoEnds + numOneEnd > 2)
  {
    /* Tested in Graphic.TypingManyChildrenTerminals */
    newcolumn->remainsGraphic = false;
    return CMR_OKAY;
  }

  bool isRoot = reducedMember == reducedComponent->root;
  DEC_MEMBER member = findMember(dec, reducedMember->member);
  if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
  {
    CMR_CALL( determineTypeParallel(dec, newcolumn, reducedMember, numOneEnd, numTwoEnds, depth) );
  }
  else if (dec->members[member].type == DEC_MEMBER_TYPE_SERIES)
  {
    CMR_CALL( determineTypeSeries(dec, newcolumn, reducedMember, numOneEnd, numTwoEnds, depth) );
  }
  else
  {
    assert(dec->members[member].type == DEC_MEMBER_TYPE_RIGID);
    CMR_CALL( determineTypeRigid(dec, newcolumn, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }

  CMRdbgMsg(6 + 2*depth, "Determined type %d (%s).\n", reducedMember->type, newcolumn->remainsGraphic ? "graphic" : "non-graphic");

  /* Parent marker edge closes cycle, so we propagate information to parent. */

  if (newcolumn->remainsGraphic && !isRoot && reducedMember->type == TYPE_CYCLE_CHILD)
  {
    DEC_MEMBER parentMember = findMemberParent(dec, reducedMember->member);
    ReducedMember* reducedParent = newcolumn->memberInfo[parentMember].reducedMember;
    DEC_EDGE markerOfParent = dec->members[member].markerOfParent;

    CMRdbgMsg(6 + 2*depth, "Marker edge closes cycle.\n");
    CMRdbgMsg(6 + 2*depth, "Parent member %d is reduced member %ld.\n", parentMember,
      (reducedParent - newcolumn->reducedMembers));

    /* Add marker edge of parent to reduced parent's path edges. */
    CMR_CALL( createPathEdge(dec, newcolumn, markerOfParent, reducedParent) );

    CMRdbgMsg(6 + 2*depth, "Added marker edge of parent to list of path edges.\n");
  }

  return CMR_OKAY;
}

/**
 * \brief Checks if a new column can be added to the decomposition.
 *
 * Information necessary for carrying out the addition is stored in \p newcolumn.
 */

CMR_ERROR addColumnCheck(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< newcolumn. */
  size_t* rows,             /**< Array of rows (with 1-entry in this column). */
  size_t numRows            /**< Length of \p rows. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(rows);

  CMRdbgMsg(0, "\n  Checking whether we can add a column with %d 1's.\n", numRows);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */

  /* Reset \ref nodesDegree to 0 and edgesInPath to false by inspecting path edges from previous iteration. */
  CMR_CALL( removeAllPathEdges(dec, newcolumn) );

#if defined(CMR_DEBUG)
  for (size_t e = 0; e < newcolumn->memEdgesInPath; ++e)
  {
    if (newcolumn->edgesInPath[e])
      CMRdbgMsg(6, "Edge %d is still marked as path edge after clean-up.\n", e);
    assert(!newcolumn->edgesInPath[e]);
  }
#endif /* CMR_DEBUG */

  /* Check for the (not yet computed) reduced decomposition whether there is a pair of parallel child/parent marker
   * edges in a non-parallel. */
  CMR_CALL( parallelParentChildCheckReducedMembers(dec, rows, numRows) );

  CMR_CALL( initializeNewColumn(dec, newcolumn) );

  /* Ensure that the edgesInPath array is long enough.
   * memEdges does not suffice since new edges can be created by splitting and we need those for new rows.
   * Each splitting introduces 4 new edges, and we might apply this twice for each series member.
   */
  size_t maxRow = 0;
  for (size_t i = 0; i < numRows; ++i)
    maxRow = rows[i] > maxRow ? rows[i] : maxRow;

  size_t requiredNumEdgesInPath = dec->memEdges + maxRow + 8*dec->numMembers;
  if (requiredNumEdgesInPath > newcolumn->memEdgesInPath)
  {
    size_t newSize = 16 + 2*newcolumn->memEdgesInPath;
    while (newSize < requiredNumEdgesInPath)
      newSize *= 2;
    CMR_CALL( CMRreallocBlockArray(dec->cmr, &newcolumn->edgesInPath, newSize) );
    for (size_t i = newcolumn->memEdgesInPath; i < newSize; ++i)
      newcolumn->edgesInPath[i] = false;
    newcolumn->memEdgesInPath = newSize;
  }

#if defined(CMR_DEBUG)
  for (size_t m = 0; m < newcolumn->memReducedMembers; ++m)
  {
    if (newcolumn->memberInfo[m].reducedMember)
      CMRdbgMsg(6, "Member %d is still mapped to reduced member before new reduced decomposition is computed.\n", m);
    assert(!newcolumn->memberInfo[m].reducedMember);
  }
#endif /* CMR_DEBUG */

  CMR_CALL( computeReducedDecomposition(dec, newcolumn, rows, numRows) );
  CMR_CALL( createReducedDecompositionPathEdges(dec, newcolumn, rows, numRows) );

  debugDot(dec, newcolumn);

  for (size_t i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    CMR_CALL( determineTypes(dec, newcolumn, &newcolumn->reducedComponents[i], newcolumn->reducedComponents[i].root,
      0) );

    CMRdbgMsg(6, "After inspecting reduced component %d, graphic = %d.\n", i, newcolumn->remainsGraphic);
  }

  /* Clean-up membersToReducedMembers. */
  for (size_t i = 0; i < newcolumn->numReducedMembers; ++i)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[i];
    newcolumn->memberInfo[reducedMember->member].reducedMember = NULL;
  }

  if (newcolumn->remainsGraphic)
    CMRdbgMsg(4, "Adding the column would maintain graphicness.\n");

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

/**
 * \brief Add a terminal of a reduced component.
 */

static
CMR_ERROR addTerminal(
  Dec* dec,                           /**< Decomposition. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  DEC_MEMBER member,                  /**< New terminal member. */
  DEC_NODE node                       /**< New terminal node. */
)
{
  CMR_UNUSED(dec);

  assert(reducedComponent);
  assert(member != SIZE_MAX);
  assert(isRepresentativeMember(dec, member));
  assert(dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL || node != SIZE_MAX);

  /* For parallels we don't need to provide a node. */
  assert(node != SIZE_MAX || dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL);
  assert(reducedComponent->numTerminals != 1 || node != SIZE_MAX || member == reducedComponent->terminalMember[0]);

  CMRdbgMsg(12, "Setting terminal node %d of member %d.\n", node, member);

  if (reducedComponent->numTerminals == 2)
  {
    CMRdbgMsg(8, "Attempted to add terminal but already 2 known. We should detect non-graphicness soon.\n");
  }
  else
  {
    reducedComponent->terminalMember[reducedComponent->numTerminals] = member;
    reducedComponent->terminalNode[reducedComponent->numTerminals] = node;
    reducedComponent->numTerminals++;
  }

  return CMR_OKAY;
}

/**
 * \brief Moves the reduced root upwards as long as it has a TYPE_DOUBLE_CHILD child.
 */

static
CMR_ERROR moveReducedRoot(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedComponent* reducedComponent  /**< Reduced component. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);

  ReducedMember* reducedMember = reducedComponent->root;
  DEC_MEMBER member = findMember(dec, reducedMember->member);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  CMR_CALL( countChildrenTypes(dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  bool cycleWithUniqueEndChild;
  if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
  {
    cycleWithUniqueEndChild = (numTwoEnds == 1 || numOneEnd == 1);
  }
  else if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
  {
    if (numTwoEnds == 1 || numOneEnd == 1)
    {
      /* For root rigid members, we have to check whether the path edges form a cycle with a two-end child. */
      DEC_NODE childMarkerNodes[2] = {
        findEdgeTail(dec, childMarkerEdges[0]),
        findEdgeHead(dec, childMarkerEdges[0])
      };
      cycleWithUniqueEndChild = (reducedMember->rigidEndNodes[2] == SIZE_MAX
        && ((reducedMember->rigidEndNodes[0] == childMarkerNodes[0]
        && reducedMember->rigidEndNodes[1] == childMarkerNodes[1])
        || (reducedMember->rigidEndNodes[0] == childMarkerNodes[1]
        && reducedMember->rigidEndNodes[1] == childMarkerNodes[0])));
    }
    else
      cycleWithUniqueEndChild = false;
  }
  else
  {
    /* For root series, the parent marker is not a path edge, so we cannot close a cycle with a child marker edge. */
    cycleWithUniqueEndChild = false;
  }

  if (!cycleWithUniqueEndChild)
  {
    CMRdbgMsg(6, "Reduced root member does not close a cycle with a 1- or 2-end child marker edge.\n");
    return CMR_OKAY;
  }

  while (cycleWithUniqueEndChild)
  {
    CMRdbgMsg(6, "Reduced %s member %d closes a cycle with a 1- or 2-end child marker edge.\n",
      memberTypeString(dec->members[member].type), member);

    DEC_MEMBER childMember = findMember(dec, dec->edges[childMarkerEdges[0]].childMember);

    /* Find the unique child member in order to process that. */
    for (int c = 0; c < reducedMember->numChildren; ++c)
    {
      Type type = reducedMember->children[c]->type;
      if (type == TYPE_SINGLE_CHILD || type == TYPE_DOUBLE_CHILD)
      {
        reducedMember = reducedMember->children[c];
        break;
      }
    }

    /* Now reduced member corresponds to the child, so we can mark its parent marker as a path edge. */
    assert(reducedMember->member == childMember);
//     newcolumn->edgesInPath[dec->members[childMember].markerToParent] = true;
    CMR_CALL( createPathEdge(dec, newcolumn, dec->members[childMember].markerToParent, reducedMember) );

    member = reducedMember->member;
    assert(isRepresentativeMember(dec, member));
    CMR_CALL( countChildrenTypes(dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

    if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
    {
      assert(reducedMember->firstPathEdge && !reducedMember->firstPathEdge->nextSibling);
      cycleWithUniqueEndChild = (numOneEnd == 1 || numTwoEnds == 1);
    }
    else if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
    {
      if (numTwoEnds == 1 || numOneEnd == 1)
      {
        /* For non-root rigid members we have to check whether the path edges together with the parent marker edge form
         * a cycle with a two-end child. */
        DEC_NODE parentMarkerNodes[2] = {
          findEdgeTail(dec, dec->members[member].markerToParent),
          findEdgeHead(dec, dec->members[member].markerToParent)
        };
        DEC_NODE childMarkerNodes[2] = {
          findEdgeTail(dec, childMarkerEdges[0]),
          findEdgeHead(dec, childMarkerEdges[0])
        };

        int numEndNodes = reducedMember->rigidEndNodes[0] == SIZE_MAX ? 0 : (reducedMember->rigidEndNodes[2] == SIZE_MAX ? 2 : 4);
        if (numEndNodes == 0)
        {
          /* Without path edges the child marker would have to be parallel to the parent marker, which is detected
           * during typing. */
          cycleWithUniqueEndChild = false;
          break;
        }

        /* Determine the end nodes of the path including the parent marker edge. */
        DEC_NODE endNodes[2];
        if (numEndNodes == 4)
        {
          endNodes[0] = reducedMember->rigidEndNodes[1];
          endNodes[1] = reducedMember->rigidEndNodes[3];
        }
        else if (reducedMember->rigidEndNodes[0] == parentMarkerNodes[0])
        {
          endNodes[0] = parentMarkerNodes[1];
          endNodes[1] = reducedMember->rigidEndNodes[1];
        }
        else
        {
          assert(reducedMember->rigidEndNodes[0] == parentMarkerNodes[1]);
          endNodes[0] = parentMarkerNodes[0];
          endNodes[1] = reducedMember->rigidEndNodes[1];
        }

        cycleWithUniqueEndChild = (endNodes[0] == childMarkerNodes[0] && endNodes[1] == childMarkerNodes[1])
          || (endNodes[0] == childMarkerNodes[1] && endNodes[1] == childMarkerNodes[0]);
      }
      else
        cycleWithUniqueEndChild = false;
    }
    else
    {
      assert(dec->members[member].type == DEC_MEMBER_TYPE_SERIES);
      if (numOneEnd == 1 || numTwoEnds == 1)
      {
        int countPathEdges = 0;
        for (PathEdge* pathEdge = reducedMember->firstPathEdge; pathEdge != NULL; pathEdge = pathEdge->nextSibling)
          ++countPathEdges;
        cycleWithUniqueEndChild = countPathEdges == dec->members[member].numEdges - 1;
      }
      else
        cycleWithUniqueEndChild = false;
    }
  }

  CMRdbgMsg(6, "Updated reduced root member is %d.\n", reducedMember->member);

  reducedComponent->root = reducedMember;

  return CMR_OKAY;
}

/**
 * \brief Sets the tail/head nodes of a \p edge to \p tail and \p head, respectively.
 */

static
CMR_ERROR setEdgeNodes(
  Dec* dec,       /**< Decomposition. */
  DEC_EDGE edge,  /**< Reduced component. */
  DEC_NODE tail,  /**< New tail node. */
  DEC_NODE head   /**< New head node. */
)
{
  assert(dec);
  assert(edge < dec->memEdges);

  dec->edges[edge].tail = tail;
  dec->edges[edge].head = head;

  return CMR_OKAY;
}

/**
 * \brief Merges \p member into its parent.
 *
 * If \p headToHead is \c true, then the two head nodes (and the two tail nodes) are identified with each other.
 * Otherwise, a head is identified with a tail.
 */

static
CMR_ERROR mergeMemberIntoParent(
  Dec* dec,          /**< Decomposition. */
  DEC_MEMBER member, /**< Reduced component. */
  bool headToHead    /**< Whether the heads of the edges shall be joined. */
)
{
  assert(dec);
  assert(member != SIZE_MAX);

  member = findMember(dec, member);
  DEC_MEMBER parentMember = findMemberParent(dec, member);
  assert(parentMember != SIZE_MAX);

  CMRdbgMsg(10, "Merging child member %d into its parent member %d.\n", member, parentMember);

#if defined(CMR_DEBUG)
  DEC_EDGE edge = dec->members[member].firstEdge;
  do
  {
    if (dec->edges[edge].head == SIZE_MAX || dec->edges[edge].tail == SIZE_MAX)
      CMRdbgMsg(10, "Edge %d of merge member %d does not have nodes.\n", edge, member);
    assert(dec->edges[edge].tail != SIZE_MAX);
    assert(dec->edges[edge].head != SIZE_MAX);
    edge = dec->edges[edge].next;
  }
  while (edge != dec->members[member].firstEdge);
#endif /* CMR_DEBUG */

  DEC_EDGE parentEdge = dec->members[member].markerOfParent;
  assert(parentEdge != SIZE_MAX);
  DEC_EDGE childEdge = dec->members[member].markerToParent;
  assert(childEdge != SIZE_MAX);

  DEC_NODE parentEdgeNodes[2] = { findEdgeTail(dec, parentEdge), findEdgeHead(dec, parentEdge) };
  DEC_NODE childEdgeNodes[2] = { findEdgeTail(dec, childEdge), findEdgeHead(dec, childEdge) };
  CMRdbgMsg(10, "Merging member %d into %d, identifying %d = {%d,%d} with %d = {%d,%d}.\n", member, parentMember,
    childEdge, childEdgeNodes[0], childEdgeNodes[1], parentEdge, parentEdgeNodes[0], parentEdgeNodes[1]);

  /* Identify nodes. */

  dec->nodes[childEdgeNodes[0]].representativeNode = parentEdgeNodes[headToHead ? 0 : 1];
  dec->nodes[childEdgeNodes[1]].representativeNode = parentEdgeNodes[headToHead ? 1 : 0];

  /* Identify members. */

  dec->members[member].representativeMember = parentMember;

  /* We add the member's edges to the parent's edge list and thereby remove the two marker edges. */
  if (dec->members[parentMember].firstEdge == parentEdge)
    dec->members[parentMember].firstEdge = dec->edges[parentEdge].next;

  dec->edges[dec->edges[parentEdge].next].prev = dec->edges[childEdge].prev;
  dec->edges[dec->edges[parentEdge].prev].next = dec->edges[childEdge].next;
  dec->edges[dec->edges[childEdge].next].prev = dec->edges[parentEdge].prev;
  dec->edges[dec->edges[childEdge].prev].next = dec->edges[parentEdge].next;
  dec->members[parentMember].numEdges += dec->members[member].numEdges - 2;
  dec->numEdges -= 2;
  dec->edges[parentEdge].next = dec->firstFreeEdge;
  dec->edges[childEdge].next = parentEdge;
  dec->firstFreeEdge = childEdge;
  dec->members[parentMember].type = DEC_MEMBER_TYPE_RIGID;

  return CMR_OKAY;
}

/**
 * \brief Create nodes of a parallel member.
 */

static
CMR_ERROR createParallelNodes(
  Dec* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< A parallel member. */
)
{
  assert(dec);
  assert(member != SIZE_MAX);
  member = findMember(dec, member);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL);

  DEC_EDGE edge = dec->members[member].firstEdge;
  assert(dec->edges[edge].head == SIZE_MAX);
  assert(dec->edges[edge].tail == SIZE_MAX);
  if (dec->edges[edge].head != SIZE_MAX)
  {
    assert(dec->edges[edge].tail != SIZE_MAX);
    return CMR_OKAY;
  }

  DEC_NODE tail, head;
  CMR_CALL( createNode(dec, &tail) );
  CMR_CALL( createNode(dec, &head) );

  do
  {
    assert(dec->edges[edge].tail == SIZE_MAX);
    assert(dec->edges[edge].head == SIZE_MAX);

    dec->edges[edge].tail = tail;
    dec->edges[edge].head = head;
    edge = dec->edges[edge].next;
  }
  while (edge != dec->members[member].firstEdge);

  return CMR_OKAY;
}

/**
 * \brief Splits a parallel member into two parallel members joined via a series member.
 */

static
CMR_ERROR splitParallel(
  Dec* dec,                   /**< Decomposition. */
  DEC_MEMBER parallel,        /**< A parallel member. */
  DEC_EDGE edge1,             /**< First edge to be moved into the child parallel. */
  DEC_EDGE edge2,             /**< Second edge to be moved into the child parallel. */
  DEC_MEMBER* pChildParallel  /**< Pointer to storing the newly created child parallel. */
)
{
  assert(dec);
  assert(parallel < dec->memMembers);
  assert(edge1 != SIZE_MAX);
  assert(edge1 < dec->memEdges);
  assert(edge2 != SIZE_MAX);
  assert(edge2 < dec->memEdges);

  DEC_MEMBER childParallel;
  CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &childParallel) );
  DEC_EDGE markerOfParentParallel, markerOfChildParallel;
  CMR_CALL( createMarkerEdgePair(dec, parallel, &markerOfParentParallel, -1, -1, childParallel, &markerOfChildParallel, -1, -1) );
  CMR_CALL( addEdgeToMembersEdgeList(dec, markerOfParentParallel) );
  CMR_CALL( addEdgeToMembersEdgeList(dec, markerOfChildParallel) );

  CMR_CALL( removeEdgeFromMembersEdgeList(dec, edge1) );
  dec->edges[edge1].member = childParallel;
  CMR_CALL( addEdgeToMembersEdgeList(dec, edge1) );
  dec->members[findMember(dec, dec->edges[edge1].childMember)].parentMember = childParallel;

  CMR_CALL( removeEdgeFromMembersEdgeList(dec, edge2) );
  dec->edges[edge2].member = childParallel;
  CMR_CALL( addEdgeToMembersEdgeList(dec, edge2) );
  dec->members[findMember(dec, dec->edges[edge2].childMember)].parentMember = childParallel;

  if (pChildParallel)
    *pChildParallel = childParallel;

  return CMR_OKAY;
}

/**
 * \brief Processes a parallel member when adding a column.
 */

static
CMR_ERROR addColumnProcessParallel(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  DEC_MEMBER member = findMember(dec, reducedMember->member);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  CMR_CALL( countChildrenTypes(dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  CMRdbgMsg(6 + 2*depth, "addColumnProcessParallel for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", member, (reducedMember - newcolumn->reducedMembers), numOneEnd, numTwoEnds, reducedMember->type);

  if (depth == 0)
  {
    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      /* Tested in UpdateRootParallelNoChildren. */

      assert(reducedMember->firstPathEdge);
      CMR_CALL( addTerminal(dec, reducedComponent, member, -1) );
      CMR_CALL( addTerminal(dec, reducedComponent, member, -1) );
      return CMR_OKAY;
    }
    else
    {
      /* Tested in UpdateRootParallelTwoSingleChildren. */

      /* If we have one single-child or a double-child then we should have moved the reduced root. */
      assert(numOneEnd == 2);
      assert(reducedComponent->numTerminals == 2);

      /* If the paralle contains more than 3 edges, we split the two child edges off, which creates a parallel with a
       * parent marker and the two children.
       * Test: Graphic.RootParallelTwoOneEnds */

      if (dec->members[member].numEdges > 3)
      {
        /* Tested in UpdateRootParallelNoChildrenSplit. */
        CMR_CALL( splitParallel(dec, member, childMarkerEdges[0], childMarkerEdges[1], &member) );
        assert(isRepresentativeMember(dec, member));
        reducedMember->member = member;
      }

      debugDot(dec, newcolumn);

      assert(dec->members[member].numEdges == 3);
      CMR_CALL( createParallelNodes(dec, member) );
      CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[1]].childMember,
        reducedMember->firstPathEdge == NULL && reducedMember->type != TYPE_DOUBLE_CHILD) );

      debugDot(dec, newcolumn);

      return CMR_OKAY;
    }
  }
  else
  {
    /* Tested in UpdateInnerParallelOneSingleChild. */

    /* Cannot be a leaf since then it would contain a path edge, i.e., it closes a cycle. */
    assert(numOneEnd == 1);
    assert(reducedComponent->numTerminals >= 1);

    DEC_NODE tail, head;
    CMR_CALL( createNode(dec, &tail) );
    CMR_CALL( createNode(dec, &head) );
    CMRdbgMsg(8 + 2*depth, "Parallel's tail node is %d and head node is %d.\n", tail, head);
    DEC_EDGE edge = dec->members[member].firstEdge;
    do
    {
      CMR_CALL( setEdgeNodes(dec, edge, tail, head) );
      edge = dec->edges[edge].next;
    }
    while (edge != dec->members[member].firstEdge);

    CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember,
      !reducedMember->firstPathEdge) );

    debugDot(dec, newcolumn);
  }

  return CMR_OKAY;
}

static inline
void flipEdge(
  Dec* dec,    /**< t-decomposition. */
  DEC_EDGE edge /**< edge. */
)
{
  assert(dec);
  assert(edge < dec->memEdges);

  SWAP_INTS(dec->edges[edge].tail, dec->edges[edge].head);
}

/**
 * \brief Processes a rigid member when adding a column.
 */

static
CMR_ERROR addColumnProcessRigid(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  DEC_MEMBER member = findMember(dec, reducedMember->member);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_RIGID);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  CMR_CALL( countChildrenTypes(dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  CMRdbgMsg(6 + 2*depth, "addColumnProcessRigid for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", reducedMember->member, (reducedMember - newcolumn->reducedMembers), numOneEnd,
    numTwoEnds, reducedMember->type);

  DEC_NODE parentMarkerNodes[2] = {
    dec->members[member].markerToParent == SIZE_MAX ? SIZE_MAX : findEdgeTail(dec, dec->members[member].markerToParent),
    dec->members[member].markerToParent == SIZE_MAX ? SIZE_MAX : findEdgeHead(dec, dec->members[member].markerToParent)
  };
  DEC_NODE childMarkerNodes[4] = {
    childMarkerEdges[0] == SIZE_MAX ? SIZE_MAX : findEdgeTail(dec, childMarkerEdges[0]),
    childMarkerEdges[0] == SIZE_MAX ? SIZE_MAX : findEdgeHead(dec, childMarkerEdges[0]),
    childMarkerEdges[1] == SIZE_MAX ? SIZE_MAX : findEdgeTail(dec, childMarkerEdges[1]),
    childMarkerEdges[1] == SIZE_MAX ? SIZE_MAX : findEdgeHead(dec, childMarkerEdges[1])
  };

  DEC_NODE* pathEndNodes = reducedMember->rigidEndNodes;
  int numPathEndNodes = pathEndNodes[0] == SIZE_MAX ? 0 : (pathEndNodes[2] == SIZE_MAX ? 2 : 4 );

  if (depth == 0)
  {
    /* Root rigid member. */

    if (dec->members[member].markerToParent != SIZE_MAX && newcolumn->edgesInPath[dec->members[member].markerToParent])
    {
      /* The parent marker is a path edge, so we modify the endNodes array. */
      if (numPathEndNodes == 0)
      {
        /* Tested in UpdateRootRigidParentOnly. */

        pathEndNodes[0] = parentMarkerNodes[0];
        pathEndNodes[1] = parentMarkerNodes[1];
        numPathEndNodes = 1;
      }
      else if (numPathEndNodes == 2)
      {
        if (pathEndNodes[0] == parentMarkerNodes[0])
          pathEndNodes[0] = parentMarkerNodes[1];
        else if (pathEndNodes[0] == parentMarkerNodes[1])
          pathEndNodes[0] = parentMarkerNodes[0];
      }
      else
      {
        /* Tested in UpdateRootRigidParentJoinsPaths. */

        pathEndNodes[0] = pathEndNodes[3];
        pathEndNodes[2] = -1;
        pathEndNodes[3] = -1;
        numPathEndNodes = 2;
      }
    }

    assert(numPathEndNodes <= 2);

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      /* Tested in UpdateRootRigidNoChildren. */

      assert(numPathEndNodes >= 2);
      CMR_CALL( addTerminal(dec, reducedComponent, member, pathEndNodes[0]) );
      CMR_CALL( addTerminal(dec, reducedComponent, member, pathEndNodes[1]) );
    }
    else if (numOneEnd == 1)
    {
      /* Tested in UpdateRootRigidOneSingleChild. */

      CMR_CALL( addTerminal(dec, reducedComponent, member,
        (pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[0] == childMarkerNodes[1])
        ? pathEndNodes[1] : pathEndNodes[0] ) );

      DEC_MEMBER childMember = findMember(dec, dec->edges[childMarkerEdges[0]].childMember);
      bool headToHead = pathEndNodes[0] == childMarkerNodes[1] || pathEndNodes[1] == childMarkerNodes[1];
      assert(headToHead || pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[0]);

      CMR_CALL( mergeMemberIntoParent(dec, childMember, headToHead) );
    }
    else
    {
      /* Tested in UpdateRootRigidTwoSingleChildren. */

      assert(numOneEnd == 2);
      assert(reducedComponent->numTerminals == 2);

      DEC_MEMBER childMember[2] = {
        findMember(dec, dec->edges[childMarkerEdges[0]].childMember),
        findMember(dec, dec->edges[childMarkerEdges[1]].childMember)
      };

      /* Count to how many path end nodes each child marker is incident. */
      int numIncidentPathNodes[2] = { 0, 0 };
      for (int c = 0; c < 2; ++c)
      {
        for (int i = 0; i < numPathEndNodes; ++i)
        {
          for (int j = 0; j < 2; ++j)
          {
            if (pathEndNodes[i] == childMarkerNodes[2*c + j])
              numIncidentPathNodes[c]++;
          }
        }
      }
      CMRdbgMsg(8 + 2*depth,
        "Child marker %d = {%d,%d} has %d incident path end nodes and child marker %d = {%d,%d} has %d.\n",
        childMarkerEdges[0], childMarkerNodes[0], childMarkerNodes[1], numIncidentPathNodes[0], childMarkerEdges[1],
        childMarkerNodes[2], childMarkerNodes[3], numIncidentPathNodes[1]);

      /* If a child marker is incident to both path ends, then we ensure it is the second one. */
      if (numIncidentPathNodes[0] == 2)
      {
        SWAP_INTS(childMember[0], childMember[1]);
        SWAP_INTS(childMarkerNodes[0], childMarkerNodes[2]);
        SWAP_INTS(childMarkerNodes[1], childMarkerNodes[3]);
      }

      /* Check if the two child markers are parallel. We then create a parallel with the two. */
      if ((childMarkerNodes[0] == childMarkerNodes[2] && childMarkerNodes[1] == childMarkerNodes[3])
        || (childMarkerNodes[0] == childMarkerNodes[3] && childMarkerNodes[1] == childMarkerNodes[2]))
      {
        /* Tested in UpdateRootRigidTwoParallelSingleChildren. */

        CMRdbgMsg(8, "Moving child marker edges %d = {%d,%d} and %d = {%d,%d} to a new parallel.\n", childMarkerEdges[0],
          childMarkerNodes[0], childMarkerNodes[1], childMarkerEdges[1], childMarkerNodes[2], childMarkerNodes[3]);

        DEC_MEMBER newParallel = -1;
        CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &newParallel) );
        dec->members[newParallel].parentMember = member;
        dec->members[childMember[0]].parentMember = newParallel;
        dec->members[childMember[1]].parentMember = newParallel;

        DEC_EDGE markerOfParent, markerToParent;
        CMR_CALL( createMarkerEdgePair(dec, member, &markerOfParent, childMarkerNodes[0], childMarkerNodes[1],
          newParallel, &markerToParent, -1, -1) );

        CMR_CALL( replaceEdgeInMembersEdgeList(dec, childMarkerEdges[0], markerOfParent) );
        CMR_CALL( addEdgeToMembersEdgeList(dec, markerToParent) );
        dec->edges[childMarkerEdges[0]].member = newParallel;
        CMR_CALL( addEdgeToMembersEdgeList(dec, childMarkerEdges[0]) );
        CMR_CALL( removeEdgeFromMembersEdgeList(dec, childMarkerEdges[1]) );
        dec->edges[childMarkerEdges[1]].member = newParallel;
        CMR_CALL( addEdgeToMembersEdgeList(dec, childMarkerEdges[1]) );

        /* The parallel has no nodes, so we have to get rid of these. */
        dec->edges[childMarkerEdges[0]].tail = -1;
        dec->edges[childMarkerEdges[0]].head = -1;
        dec->edges[childMarkerEdges[1]].tail = -1;
        dec->edges[childMarkerEdges[1]].head = -1;

        debugDot(dec, newcolumn);

        /* We have to merge in the parallel. */
        CMR_CALL( createParallelNodes(dec, newParallel) );
        CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember, true) );
        /* If there is no path, then we merge heads to heads. Otherwise, one head is mapped to tail, such that the path
         * is used. */
        CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[1]].childMember, numPathEndNodes == 0) );

        debugDot(dec, newcolumn);

        return CMR_OKAY;
      }

      if (numPathEndNodes == 0)
      {
        /* We fake existence of a path by setting both path end nodes to the node common to both child markers. */
        if (childMarkerNodes[0] == childMarkerNodes[2] || childMarkerNodes[0] == childMarkerNodes[3])
        {
          pathEndNodes[0] = childMarkerNodes[0];
          pathEndNodes[1] = childMarkerNodes[0];
        }
        else
        {
          assert(childMarkerNodes[1] == childMarkerNodes[2] || childMarkerNodes[1] == childMarkerNodes[3]);
          pathEndNodes[0] = childMarkerNodes[1];
          pathEndNodes[1] = childMarkerNodes[1];
        }
      }

      if (pathEndNodes[0] != childMarkerNodes[0] && pathEndNodes[0] != childMarkerNodes[1])
        SWAP_INTS(pathEndNodes[0], pathEndNodes[1]);

      CMRdbgMsg(8 + 2*depth, "After swapping, we have a %d-%d-path as well as child markers {%d,%d} and {%d,%d}.\n",
        pathEndNodes[0], pathEndNodes[1], childMarkerNodes[0], childMarkerNodes[1], childMarkerNodes[2], childMarkerNodes[3]);

      assert(pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[0] == childMarkerNodes[1]);
      CMR_CALL( mergeMemberIntoParent(dec, childMember[0], pathEndNodes[0] == childMarkerNodes[1]) );
      debugDot(dec, newcolumn);

      CMR_CALL( mergeMemberIntoParent(dec, childMember[1], pathEndNodes[1] == childMarkerNodes[3]) );
      debugDot(dec, newcolumn);
    }
  }
  else
  {
    /* Non-root rigid member. */

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      /* Tested in UpdateLeafRigid. */

      assert(reducedMember->firstPathEdge);
      assert(pathEndNodes[0] != SIZE_MAX);

      CMR_CALL( addTerminal(dec, reducedComponent, member, pathEndNodes[1]) );
      if (parentMarkerNodes[0] == pathEndNodes[0])
        flipEdge(dec, dec->members[member].markerToParent);
    }
    else
    {
      assert(numOneEnd == 1);

      if (numPathEndNodes >= 2)
      {
        /* Tested in UpdateInnerRigidOnePath. */

        CMRdbgMsg(8 + 2*depth, "%d-%d-path with parent marker {%d,%d} and child marker {%d,%d}.\n", pathEndNodes[0],
          pathEndNodes[1], parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0], childMarkerNodes[1]);

        /* Ensure that the child marker is incident to the path end node 1. */
        if (pathEndNodes[1] != childMarkerNodes[0] && pathEndNodes[1] != childMarkerNodes[1])
          SWAP_INTS(pathEndNodes[0], pathEndNodes[1]);
        assert(pathEndNodes[1] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[1]);

        /* Flip parent if necessary. */
        if (pathEndNodes[0] == parentMarkerNodes[0])
        {
          flipEdge(dec, dec->members[member].markerToParent);
          SWAP_INTS(parentMarkerNodes[0], parentMarkerNodes[1]);
        }

        assert(pathEndNodes[0] == parentMarkerNodes[1]);

        /* Merge child. */
        CMR_CALL( mergeMemberIntoParent(dec, findMember(dec, dec->edges[childMarkerEdges[0]].childMember),
          pathEndNodes[1] == childMarkerNodes[1]) );
        debugDot(dec, newcolumn);
      }
      else
      {
        /* Tested in UpdateInnerRigidNoPath. */

        /* Parent marker and child marker must be next to each other. */
        if (parentMarkerNodes[0] == childMarkerNodes[0] || parentMarkerNodes[0] == childMarkerNodes[1])
          flipEdge(dec, dec->members[member].markerToParent);

        CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember,
          parentMarkerNodes[0] == childMarkerNodes[1] || parentMarkerNodes[1] == childMarkerNodes[1]) );

        debugDot(dec, newcolumn);
      }
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Replaces an edge by a parallel containing it.
 *
 * The given member should have at least two edges.
 */

static
CMR_ERROR createEdgeParallel(
  Dec* dec,                 /**< Decomposition. */
  DEC_EDGE edge,            /**< Edge to be replaced by a parallel containing that edge. */
  DEC_MEMBER* pNewParallel  /**< Pointer for storing the new parallel. */
)
{
  assert(dec);

  CMRdbgMsg(8, "Creating parallel for edge %d.\n", edge);

  DEC_MEMBER parentMember = findEdgeMember(dec, edge);
  DEC_MEMBER newParallel = -1;
  CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &newParallel) );
  dec->members[newParallel].parentMember = parentMember;

  DEC_EDGE markerOfParent, markerToParent;
  CMR_CALL( createMarkerEdgePair(dec, parentMember, &markerOfParent, dec->edges[edge].tail, dec->edges[edge].head,
    newParallel, &markerToParent, -1, -1) );
  dec->edges[markerOfParent].next = dec->edges[edge].next;
  dec->edges[markerOfParent].prev = dec->edges[edge].prev;
  assert(dec->edges[markerOfParent].next != markerOfParent);
  dec->edges[dec->edges[markerOfParent].next].prev = markerOfParent;
  dec->edges[dec->edges[markerOfParent].prev].next = markerOfParent;
  if (dec->members[parentMember].firstEdge == edge)
    dec->members[parentMember].firstEdge = markerOfParent;

  CMR_CALL( addEdgeToMembersEdgeList(dec, markerToParent) );

  dec->edges[edge].member = newParallel;
  CMR_CALL( addEdgeToMembersEdgeList(dec, edge) );

  if (pNewParallel)
    *pNewParallel = newParallel;

  return CMR_OKAY;
}

/**
 * \brief Splits series member into two, connecting them via a parallel.
 *
 * Takes all edges of the series \p member for which \p edgesPredicate is the same as
 * \p predicateValue.
 */

static
CMR_ERROR splitSeries(
  Dec* dec,                       /**< Decomposition. */
  DEC_MEMBER member,              /**< Series member to be split. */
  bool* edgesPredicate,           /**< Map from edges to predicate. */
  bool predicateValue,            /**< Value of predicate. */
  DEC_EDGE* pRepresentativeEdge,  /**< Pointer for storing the child marker edge that links to new parallel. */
  DEC_MEMBER* pNewParallel,       /**< Pointer for storing the new parallel. */
  DEC_MEMBER* pNewSeries          /**< Pointer for storing the new parallel. */
)
{
  assert(dec);
  assert(member < dec->memMembers);
  assert(edgesPredicate);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_SERIES);

#if defined(CMR_DEBUG)
  CMRdbgMsg(8, "Checking series member %d for splitting... ", member);
#endif /* CMR_DEBUG */

  DEC_EDGE edge = dec->members[member].firstEdge;
  DEC_EDGE someSatisfyingEdge = -1;
  int numSatisfyingEdges = 0;
  do
  {
    bool value = edgesPredicate[edge];
    if ((value && predicateValue) || (!value && !predicateValue))
    {
      someSatisfyingEdge = edge;
      ++numSatisfyingEdges;
    }
    edge = dec->edges[edge].next;
  } while (edge != dec->members[member].firstEdge);

  CMRdbgMsg(0, "%d of %d edges satisfy the predicate.\n", numSatisfyingEdges, dec->members[member].numEdges);

  if (numSatisfyingEdges == 0)
  {
    if (pRepresentativeEdge)
      *pRepresentativeEdge = -1;
    if (pNewParallel)
      *pNewParallel = -1;
    if (pNewSeries)
      *pNewSeries = -1;
  }
  else if (numSatisfyingEdges == 1)
  {
    if (pNewSeries)
      *pNewSeries = -1;
    if (pRepresentativeEdge)
      *pRepresentativeEdge = someSatisfyingEdge;
    if (pNewParallel)
      *pNewParallel = -1;
  }
  else
  {
    /* Initialize new series member. */
    DEC_MEMBER series;
    CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_SERIES, &series) );

    /* Initialize new parallel. */
    DEC_MEMBER parallel;
    CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &parallel) );

    DEC_EDGE seriesParentMarker, parallelChildMarker;
    CMR_CALL( createMarkerEdgePair(dec, parallel, &parallelChildMarker, -1, -1, series, &seriesParentMarker, -1, -1) );
    CMR_CALL( addEdgeToMembersEdgeList(dec, seriesParentMarker) );
    CMR_CALL( addEdgeToMembersEdgeList(dec, parallelChildMarker) );

    /* Go through old series member. */
    DEC_EDGE firstEdge = dec->members[member].firstEdge;
    edge = firstEdge;
    bool encounteredStayingEdge = false;
    do
    {
#if defined(CMR_DEBUG_SPLITTING)
      CMRdbgMsg(8, "Edge %d <%d>", edge, dec->edges[edge].name);
      if (dec->edges[edge].childMember != SIZE_MAX)
        CMRdbgMsg(0, " (with child %d)", dec->edges[edge].childMember);
      if (edge == dec->members[member].markerToParent)
        CMRdbgMsg(0, " (with parent %d)", dec->members[member].parentMember);
      CMRdbgMsg(0, " (prev = %d, next = %d)", dec->edges[edge].prev, dec->edges[edge].next);
#endif /* CMR_DEBUG_SPLITTING*/

      /* Evaluate predicate. */
      bool value = edgesPredicate[edge];
      if ((value && !predicateValue) || (!value && predicateValue))
      {
#if defined(CMR_DEBUG_SPLITTING)
        CMRdbgMsg(" does not satisfy the predicate.\n");
#endif /* CMR_DEBUG_SPLITTING */
        edge = dec->edges[edge].next;
        encounteredStayingEdge = true;
        continue;
      }

#if defined(CMR_DEBUG_SPLITTING)
      CMRdbgMsg(0, " satisfies the predicate.\n");
#endif /* CMR_DEBUG_SPLITTING */

      assert(edge != dec->members[member].markerToParent);

      /* Remove edge from old edge list. */
      DEC_EDGE oldPrev = dec->edges[edge].prev;
      DEC_EDGE oldNext = dec->edges[edge].next;
      dec->edges[oldPrev].next = oldNext;
      dec->edges[oldNext].prev = oldPrev;
      dec->members[member].numEdges--;

      /* Add edge to new edge list. */
      DEC_EDGE newPrev = dec->edges[seriesParentMarker].prev;
      dec->edges[newPrev].next = edge;
      dec->edges[seriesParentMarker].prev = edge;
      dec->edges[edge].prev = newPrev;
      dec->edges[edge].next = seriesParentMarker;
      dec->edges[edge].member = series;
      if (dec->edges[edge].childMember != SIZE_MAX)
      {
        assert( dec->members[dec->edges[edge].childMember].parentMember == member);
        dec->members[dec->edges[edge].childMember].parentMember = series;
      }
      dec->members[series].numEdges++;

      /* Did we move the first edge of this member? */
      if (edge == firstEdge)
      {
        dec->members[member].firstEdge = oldNext;
        firstEdge = oldNext;
        edge = oldNext;
        continue;
      }

      edge = oldNext;
    }
    while (edge != firstEdge || !encounteredStayingEdge);

    DEC_EDGE memberChildMarker, parallelParentMarker;
    CMR_CALL( createMarkerEdgePair(dec, member, &memberChildMarker, -1, -1, parallel, &parallelParentMarker, -1, -1) );
    CMR_CALL( addEdgeToMembersEdgeList(dec, parallelParentMarker) );
    DEC_EDGE oldPrev = dec->edges[firstEdge].prev;
    dec->edges[memberChildMarker].next = firstEdge;
    dec->edges[memberChildMarker].prev = oldPrev;
    dec->edges[oldPrev].next = memberChildMarker;
    dec->edges[firstEdge].prev = memberChildMarker;
    dec->members[member].numEdges++;

#if defined(CMR_DEBUG_SPLITTING)
    CMRdbgMsg(8, "Updated old series member with these edges:");
    edge = firstEdge;
    do
    {
      CMRdbgMsg(0, " %d <%d>", edge, dec->edges[edge].name);
      edge = dec->edges[edge].next;
    }
    while (edge != firstEdge);
    CMRdbgMsg(0, ".\n");
    CMRdbgMsg(8, "New series member has these edges:");
    edge = seriesParentMarker;
    do
    {
      CMRdbgMsg(0, " %d <%d>", edge, dec->edges[edge].name);
      edge = dec->edges[edge].next;
    }
    while (edge != seriesParentMarker);
    CMRdbgMsg(0, ".\n");
#endif /* CMR_DEBUG_SPLITTING */

    CMRdbgMsg(8, "Connecting parallel is member %d and squeezed series member is member %d.\n", parallel, series);

    if (pRepresentativeEdge)
      *pRepresentativeEdge = memberChildMarker;
    if (pNewParallel)
      *pNewParallel = parallel;
    if (pNewSeries)
      *pNewSeries = series;
  }
  return CMR_OKAY;
}

/**
 * \brief Processes a series member when adding a column.
 */

static
CMR_ERROR addColumnProcessSeries(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< newcolumn. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);

  DEC_MEMBER member = findMember(dec, reducedMember->member);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  CMR_CALL( countChildrenTypes(dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  CMRdbgMsg(6 + 2*depth, "addColumnProcessSeries for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", member, (reducedMember - newcolumn->reducedMembers), numOneEnd, numTwoEnds, reducedMember->type);

  if (depth == 0)
  {
    /* If we have a child containing both ends then we should have moved the reduced root there. */
    assert(numTwoEnds == 0);

    if (numOneEnd == 0)
    {
      /* Root series member containing both ends. */
      assert(reducedMember->firstPathEdge);
      DEC_EDGE representativeEdge;
      if (newcolumn->edgesInPath[dec->members[member].markerToParent])
      {
        /* Tested in UpdateRootSeriesNoChildrenParent. */

        CMRdbgMsg(8 + 2*depth, "Series contains both terminal nodes and the parent marker edge is a path edge.\n");

        /* Squeeze off all non-path edges by moving them to a new series member and creating a parallel member to
         * connect it to the remaining series member. */
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, &representativeEdge, NULL, NULL) );
      }
      else if (reducedMember->type == TYPE_CYCLE_CHILD)
      {
        /* Tested in UpdateRootSeriesNoChildrenHamiltonianPath. */

        CMRdbgMsg(8 + 2*depth, "Series contains both terminal nodes which are the parent marker edge nodes.\n");

        DEC_MEMBER parentMember = findMemberParent(dec, member);
        DEC_EDGE markerOfParent = dec->members[member].markerOfParent;
        assert(parentMember != SIZE_MAX); /* A series member can only close a cycle if it has a parent. */
        if (dec->members[parentMember].type == DEC_MEMBER_TYPE_PARALLEL)
        {
          CMR_CALL( addTerminal(dec, reducedComponent, findEdgeMember(dec, markerOfParent), -1 ) );
          CMR_CALL( addTerminal(dec, reducedComponent, findEdgeMember(dec, markerOfParent), -1 ) );
        }
        else
        {
          assert(dec->members[parentMember].type == DEC_MEMBER_TYPE_RIGID);
          CMR_CALL( addTerminal(dec, reducedComponent, findEdgeMember(dec, markerOfParent),
            findEdgeTail(dec, markerOfParent) ) );
          CMR_CALL( addTerminal(dec, reducedComponent, findEdgeMember(dec, markerOfParent),
            findEdgeHead(dec, markerOfParent) ) );
        }

        return CMR_OKAY;
      }
      else
      {
        /* Tested in UpdateRootSeriesNoChildren. */
        CMRdbgMsg(8 + 2*depth, "Series member contains both terminal nodes and a non-path parent marker edge.\n");

        /* Squeeze off all path edges by moving them to a new series member and creating a parallel member to connect it
         * to the remaining series member. */
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, &representativeEdge, NULL, NULL) );
      }

      debugDot(dec, newcolumn);

      DEC_EDGE childMember = dec->edges[representativeEdge].childMember;
      DEC_NODE tail = SIZE_MAX;
      DEC_NODE head = SIZE_MAX;
      if (childMember == SIZE_MAX)
      {
        CMR_CALL( createEdgeParallel(dec, representativeEdge, &childMember) );
        debugDot(dec, newcolumn);
      }
      else if (dec->members[childMember].type == DEC_MEMBER_TYPE_RIGID)
      {
        tail = findEdgeTail(dec, dec->members[childMember].markerToParent);
        head = findEdgeHead(dec, dec->members[childMember].markerToParent);
      }

      assert(reducedComponent->numTerminals == 0);
      CMR_CALL( addTerminal(dec, reducedComponent, childMember, tail) );
      CMR_CALL( addTerminal(dec, reducedComponent, childMember, head) );
      debugDot(dec, newcolumn);
    }
    else if (numOneEnd == 1)
    {
      /* Root series member containing one end. */

      if (newcolumn->edgesInPath[dec->members[member].markerToParent])
      {
        /* Tested in UpdateRootSeriesOneSingleChildParent. */

        /* There is more than 1 path edge, so we split off all non-path edges and work in the new series member. */
        assert(reducedMember->firstPathEdge);
        if (reducedMember->firstPathEdge->nextSibling)
        {
          CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, NULL, NULL, &member) );
          assert(isRepresentativeMember(dec, member));
          reducedMember->member = member;
          CMR_CALL( createPathEdge(dec, newcolumn, dec->members[member].markerToParent, reducedMember) );
        }
        newcolumn->edgesInPath[childMarkerEdges[0]] = true;
        DEC_EDGE nonPathEdge;
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[childMarkerEdges[0]] = false;

        DEC_NODE a, b, c;
        CMR_CALL( createNode(dec, &a) );
        CMR_CALL( createNode(dec, &b) );
        if (dec->members[member].numEdges == 3)
        {
          CMR_CALL( createNode(dec, &c) );
          CMR_CALL( setEdgeNodes(dec, nonPathEdge, a, c) );
        }
        else
          c = a;
        CMR_CALL( setEdgeNodes(dec, dec->members[member].markerToParent, a, b) );
        CMR_CALL( setEdgeNodes(dec, childMarkerEdges[0], c, b) );
        CMR_CALL( addTerminal(dec, reducedComponent, member, a) );
        CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember, true) );
        dec->members[member].type = DEC_MEMBER_TYPE_RIGID;
      }
      else
      {
        /* Tested in UpdateRootSeriesOneSingleChild. */

        /* Parent marker edge is not a path edge. */
        assert(reducedMember->firstPathEdge);

        /* If there is more than 1 path edge, we squeeze off by moving them to a new series member and creating a
         * parallel member to connect it to the remaining series member. */
        DEC_EDGE pathEdge;
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
        if (pathEdge != SIZE_MAX)
        {
          CMR_CALL( createPathEdge(dec, newcolumn, pathEdge, reducedMember) );
        }
        debugDot(dec, newcolumn);

        /* Unless the series member consists of only the parent marker, the child marker (containing a path end) and a
         * representative edge, we squeeze off the representative edge and the child marker. */
        if (dec->members[member].numEdges > 3)
        {
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, NULL, NULL, &member) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          debugDot(dec, newcolumn);
        }

        assert(dec->members[member].numEdges == 3);

        DEC_NODE a, b, c;
        CMR_CALL( createNode(dec, &a) );
        CMR_CALL( createNode(dec, &b) );
        CMR_CALL( createNode(dec, &c) );
        CMR_CALL( setEdgeNodes(dec, dec->members[member].markerToParent, b, c) );
        CMR_CALL( setEdgeNodes(dec, pathEdge, a, b) );
        CMR_CALL( setEdgeNodes(dec, childMarkerEdges[0], c, a) );
        CMR_CALL( addTerminal(dec, reducedComponent, member, b) );
        CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      }
      debugDot(dec, newcolumn);
    }
    else
    {
      assert(numOneEnd == 2);

      /* If there is more than 1 path edge, we split off by moving them to a new series member and creating a parallel
       * member to connect it to the remaining series member. */
      DEC_EDGE pathEdge = -1;
      DEC_EDGE nonPathEdge = -1;
      if (reducedMember->type != TYPE_DOUBLE_CHILD)
      {
        /* Tested in UpdateRootSeriesTwoSingleChildren. */

        /* Parent marker is a non-path edge. We split off path edges if more than one. */
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );

        /* If there is another non-path edge, we split off the series member consisting of the two relevant child marker
         * edges and possibly the path edge. We then replace the current series member by the new one. */

        if (dec->members[member].numEdges > (pathEdge != SIZE_MAX ? 4 : 3))
        {
          if (pathEdge != SIZE_MAX)
          {
            CMR_CALL( createPathEdge(dec, newcolumn, pathEdge, reducedMember) );
          }
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          newcolumn->edgesInPath[childMarkerEdges[1]] = true;
          CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, NULL, NULL, &member) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          newcolumn->edgesInPath[childMarkerEdges[1]] = false;
          assert(isRepresentativeMember(dec, member));
          reducedMember->member = member;
        }

        nonPathEdge = dec->members[member].markerToParent;
      }
      else
      {
        /* Tested in UpdateRootSeriesTwoSingleChildrenParent. */

        assert(newcolumn->edgesInPath[dec->members[member].markerToParent]);
        assert(reducedMember->firstPathEdge);

        if (reducedMember->firstPathEdge->nextSibling)
        {
          /* Parent marker is one of several path edges. */
          CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, NULL, NULL, &member) );
          assert(isRepresentativeMember(dec, member));
          reducedMember->member = member;
          CMR_CALL( createPathEdge(dec, newcolumn, dec->members[member].markerToParent, reducedMember) );
        }
        pathEdge = dec->members[member].markerToParent;
        debugDot(dec, newcolumn);

        if (dec->members[member].numEdges > 3)
        {
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          newcolumn->edgesInPath[childMarkerEdges[1]] = true;
          CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          newcolumn->edgesInPath[childMarkerEdges[1]] = false;
        }
        else
          nonPathEdge = -1;
      }

      CMRdbgMsg(8 + 2*depth,
        "After splitting off, the (potential) path edge is %d and the (potential) non-path edge is %d.\n",
        pathEdge, nonPathEdge);
      debugDot(dec, newcolumn);

      assert(pathEdge != SIZE_MAX || nonPathEdge != SIZE_MAX);

      /* a <----- b ---- c -----> d -------- a
       *   child0   path   child1   non-path
       *
       * or
       *
       * a <----- b ---- c -----> d=a
       *   child0   path   child1
       *
       * or
       *
       * a <----- b=c -----> d -------- a
       *   child0     child1   non-path */
      DEC_NODE a, b, c, d;
      CMR_CALL( createNode(dec, &a) );
      CMR_CALL( createNode(dec, &b) );
      if (pathEdge != SIZE_MAX)
        CMR_CALL( createNode(dec, &c) );
      else
        c = b;
      if (nonPathEdge != SIZE_MAX)
        CMR_CALL( createNode(dec, &d) );
      else
        d = a;
      CMR_CALL( setEdgeNodes(dec, childMarkerEdges[0], a, b) );
      CMRdbgMsg(8 + 2*depth, "First child edge %d = {%d, %d}\n", childMarkerEdges[0], a, b);
      CMR_CALL( setEdgeNodes(dec, childMarkerEdges[1], d, c) );
      CMRdbgMsg(8 + 2*depth, "Second child edge %d = {%d, %d}\n", childMarkerEdges[1], d, c);
      if (pathEdge != SIZE_MAX)
      {
        CMR_CALL( setEdgeNodes(dec, pathEdge, b, c) );
        CMRdbgMsg(8 + 2*depth, "Path edge %d = {%d, %d}\n", pathEdge, b, c);
      }
      if (nonPathEdge != SIZE_MAX)
      {
        CMR_CALL( setEdgeNodes(dec, nonPathEdge, d, a) );
        CMRdbgMsg(8 + 2*depth, "Non-path edge %d = {%d, %d}\n", nonPathEdge, d, a);
      }

      CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[1]].childMember, true) );
      dec->members[member].type = DEC_MEMBER_TYPE_RIGID;
    }
  }
  else
  {
    /* Non-root series member. */
    assert(reducedMember->type != TYPE_CYCLE_CHILD); /* addColumnProcess should never consider such a member. */
    assert(reducedMember->type != TYPE_DOUBLE_CHILD); /* This should only happen at the root. */
    assert(reducedMember->type != TYPE_ROOT); /* We are not a root. */
    assert(reducedMember->type == TYPE_SINGLE_CHILD); /* Only remaining case. */

    if (numOneEnd == 0)
    {
      /* Tested in UpdateLeafSeries. */

      assert(numOneEnd == 0);
      assert(numTwoEnds == 0);
      assert(reducedComponent->numTerminals < 2);
      assert(reducedMember->firstPathEdge);

      CMRdbgMsg(6 + 2*depth, "Non-root series member %d without one-end children.\n", member);

      /* Squeeze off all path edges by moving them to a new series member and creating a parallel member to connect it
       * to the remaining series member. */
      DEC_EDGE pathEdge = -1;
      CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
      assert(pathEdge != SIZE_MAX);
      CMR_CALL( createPathEdge(dec, newcolumn, pathEdge, reducedMember) );

      /* If necessary, we squeeze off the non-path edges as well. */
      assert(dec->members[member].numEdges >= 3);
      DEC_EDGE nonPathEdge;
      if (dec->members[member].numEdges == 3)
      {
        nonPathEdge = dec->edges[dec->members[member].markerToParent].next;
        if (nonPathEdge == pathEdge)
          nonPathEdge = dec->edges[nonPathEdge].next;
      }
      else
      {
        /* We temporarily mark the parent edge to belong to the path. */
        DEC_EDGE markerToParent = dec->members[member].markerToParent;
        newcolumn->edgesInPath[markerToParent] = true;
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[markerToParent] = false;
      }
      assert(dec->members[member].numEdges == 3);
      debugDot(dec, newcolumn);

      /* We now create the nodes of the triangle so that the path leaves it via the parent marker edge's head node. */

      DEC_NODE a, b, c;
      CMR_CALL( createNode(dec, &a) );
      CMR_CALL( createNode(dec, &b) );
      CMR_CALL( createNode(dec, &c) );
      CMR_CALL( setEdgeNodes(dec, dec->members[reducedMember->member].markerToParent, a, b) );
      CMR_CALL( setEdgeNodes(dec, pathEdge, b, c) );
      CMR_CALL( setEdgeNodes(dec, nonPathEdge, c, a) );
      CMR_CALL( addTerminal(dec, reducedComponent, reducedMember->member, c) );

      return CMR_OKAY;
    }
    else
    {
      /* Tested in UpdateInnerSeriesOneSingleChild. */
      assert(numOneEnd == 1);

      /* Squeeze off all path edges by moving them to a new series member and creating a parallel member to connect it
       * to the remaining series member. */
      DEC_EDGE pathEdge = -1;
      CMRdbgMsg(8 + 2*depth, "Splitting of path edges.\n");
      CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
      if (pathEdge != SIZE_MAX)
      {
        CMR_CALL( createPathEdge(dec, newcolumn, pathEdge, reducedMember) );
      }
      debugDot(dec, newcolumn);

      /* If necessary, we squeeze off the non-path edges as well. */
      assert(dec->members[member].numEdges >= 3);
      int numNonPathEdges = dec->members[member].numEdges - 2 - (pathEdge != SIZE_MAX ? 1 : 0);
      CMRdbgMsg(8 + 2*depth, "Number of non-path edges is %d.\n", numNonPathEdges);
      DEC_EDGE nonPathEdge;
      if (numNonPathEdges == 0)
        nonPathEdge = -1;
      else if (numNonPathEdges == 1)
      {
        nonPathEdge = dec->edges[dec->members[member].markerToParent].next;
        while (nonPathEdge == childMarkerEdges[0] || nonPathEdge == pathEdge)
          nonPathEdge = dec->edges[nonPathEdge].next;
      }
      else
      {
        newcolumn->edgesInPath[dec->members[member].markerToParent] = true;
        newcolumn->edgesInPath[childMarkerEdges[0]] = true;
        CMR_CALL( splitSeries(dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[dec->members[member].markerToParent] = false;
        newcolumn->edgesInPath[childMarkerEdges[0]] = false;
        debugDot(dec, newcolumn);
      }

      CMRdbgMsg(8 + 2*depth, "After (potential) splitting: path edge is %d and non-path edge is %d.\n",
        pathEdge, nonPathEdge);
      assert(dec->members[member].numEdges <= 4);

      /* We now create the nodes of the triangle so that the path leaves it via the parent marker edge's head node. */

      DEC_NODE a, b, c, d;
      CMR_CALL( createNode(dec, &a) );
      CMR_CALL( createNode(dec, &b) );
      CMR_CALL( createNode(dec, &c) );
      CMR_CALL( setEdgeNodes(dec, dec->members[reducedMember->member].markerToParent, a, b) );
      if (dec->members[member].numEdges == 4)
      {
        CMR_CALL( createNode(dec, &d) );
        CMR_CALL( setEdgeNodes(dec, pathEdge, b, c) );
        CMR_CALL( setEdgeNodes(dec, childMarkerEdges[0], d, c) );
        CMR_CALL( setEdgeNodes(dec, nonPathEdge, d, a) );
      }
      else if (nonPathEdge == SIZE_MAX)
      {
        CMR_CALL( setEdgeNodes(dec, pathEdge, b, c) );
        CMR_CALL( setEdgeNodes(dec, childMarkerEdges[0], a, c) );
      }
      else
      {
        assert(pathEdge == SIZE_MAX);
        CMR_CALL( setEdgeNodes(dec, childMarkerEdges[0], c, b) );
        CMR_CALL( setEdgeNodes(dec, nonPathEdge, a, c) );
      }
      debugDot(dec, newcolumn);

      CMR_CALL( mergeMemberIntoParent(dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      debugDot(dec, newcolumn);

      return CMR_OKAY;
    }
  }

  return CMR_OKAY;
}

/**
 * \brief Processes a component of a reduced decomposition before the actual modification.
 *
 * Processes the reduced members in depth-first search manner and does the following:
 * - Series members are squeezed.
 * - Terminal nodes and (reduced) members are detected.
 * - Marker edges along unique path between terminal nodes are merged.
 */

static
CMR_ERROR addColumnProcessComponent(
  Dec* dec,                           /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);

  CMRdbgMsg(6 + 2*depth, "addColumnProcessComponent(member %d = reduced member %ld)\n", reducedMember->member,
    (reducedMember - &newcolumn->reducedMembers[0]));

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */


  /* If we are non-root type 1, then we don't need to do anything. */
  if (reducedMember->type == TYPE_CYCLE_CHILD && depth > 0)
  {
    return CMR_OKAY;
  }

  /* Handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    if (child->type != TYPE_CYCLE_CHILD)
    {
      CMR_CALL( addColumnProcessComponent(dec, newcolumn, reducedComponent, reducedMember->children[c], depth+1) );
    }
    else
    {
      CMRdbgMsg(8 + 2*depth, "Member %d is implicitly replaced by a path edge.\n", findMember(dec, child->member));
    }
  }

  /* Different behavior for parallel members, series members and rigid components. */
  if (dec->members[reducedMember->member].type == DEC_MEMBER_TYPE_PARALLEL)
    CMR_CALL( addColumnProcessParallel(dec, newcolumn, reducedComponent, reducedMember, depth) );
  else if (dec->members[reducedMember->member].type == DEC_MEMBER_TYPE_SERIES)
    CMR_CALL( addColumnProcessSeries(dec, newcolumn, reducedComponent, reducedMember, depth) );
  else
    CMR_CALL( addColumnProcessRigid(dec, newcolumn, reducedComponent, reducedMember, depth) );

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */


  return CMR_OKAY;
}

/**
 * \brief Recursive method to reorder a member by making \p newParent the new parent.
 */

static
CMR_ERROR doReorderComponent(
  Dec* dec,                   /**< Decomposition. */
  DEC_MEMBER member,          /**< Member to be processed. */
  DEC_MEMBER newParent,       /**< New parent of \p member. */
  DEC_EDGE newMarkerToParent, /**< New marker edge linked to new parent of \p member. */
  DEC_EDGE markerOfNewParent  /**< Counterpart to \p newMarkerToParent. */
)
{
  assert(dec);
  assert(member != SIZE_MAX);
  assert(newParent != SIZE_MAX);

  DEC_MEMBER oldParent = findMemberParent(dec, member);
  DEC_EDGE oldMarkerToParent = dec->members[member].markerToParent;
  DEC_EDGE oldMarkerOfParent = dec->members[member].markerOfParent;

  CMRdbgMsg(8, "Flipping parenting of member %d with old parent %d and new parent %d.\n", member, oldParent, newParent);

  dec->members[member].markerToParent = newMarkerToParent;
  dec->members[member].markerOfParent = markerOfNewParent;
  dec->members[member].parentMember = newParent;
  dec->edges[markerOfNewParent].childMember = member;
  dec->edges[newMarkerToParent].childMember = -1;

  if (oldMarkerToParent != SIZE_MAX)
    CMR_CALL( doReorderComponent(dec, oldParent, member, oldMarkerOfParent, oldMarkerToParent) );

  return CMR_OKAY;
}

/**
 * \brief Reorders a component, making \p newRoot the new root.
 */

static
CMR_ERROR reorderComponent(
  Dec* dec,           /**< Decomposition. */
  DEC_MEMBER newRoot  /**< The new root of the component. */
)
{
  assert(dec);
  assert(newRoot < dec->memMembers);
  assert(isRepresentativeMember(dec, newRoot));

  CMRdbgMsg(4, "Making member %d the new root of its component.\n", newRoot);

  if (dec->members[newRoot].parentMember != SIZE_MAX)
  {
    CMR_CALL( doReorderComponent(dec, findMemberParent(dec, newRoot), newRoot,
      dec->members[newRoot].markerOfParent, dec->members[newRoot].markerToParent) );
  }

  return CMR_OKAY;
}

/**
 * \brief Actually adds a column with ones in \p rows, using information from \p newcolumn.
 */

CMR_ERROR addColumnApply(
  Dec* dec,                 /**< Decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< newcolumn. */
  size_t column,            /**< Index of new column to be added. */
  size_t* rows,             /**< Array of rows with 1-entry in this column. */
  size_t numRows            /**< Length of \p rows. */
)
{
  assert(dec);
  assert(newcolumn);
  assert(newcolumn->remainsGraphic);

  CMRdbgMsg(0, "\n  Adding a column with %d 1's.\n", numRows);

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */


  /* Create reduced components for new edges. */
  CMR_CALL( completeReducedDecomposition(dec, newcolumn, rows, numRows) );

  /* Process all reduced components individually. */

  DEC_EDGE* componentNewEdges = NULL;
  CMR_CALL( CMRallocStackArray(dec->cmr, &componentNewEdges, newcolumn->numReducedComponents) );

  size_t maxDepthComponent = SIZE_MAX;;
  CMRdbgMsg(4, "Processing %d reduced components.\n", newcolumn->numReducedComponents);
  for (size_t i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    ReducedComponent* reducedComponent = &newcolumn->reducedComponents[i];

    CMRdbgMsg(4, "Moving root of reduced component %d.\n", i);

    CMR_CALL( moveReducedRoot(dec, newcolumn, reducedComponent) );
    debugDot(dec, newcolumn);

    CMRdbgMsg(4, "Processing reduced component %d of depth %d.\n", i, reducedComponent->rootDepth);

    if (maxDepthComponent == SIZE_MAX || reducedComponent->rootDepth
      > newcolumn->reducedComponents[maxDepthComponent].rootDepth)
    {
      maxDepthComponent = i;
    }

    CMR_CALL( addColumnProcessComponent(dec, newcolumn, reducedComponent, reducedComponent->root, 0) );

    assert(reducedComponent->numTerminals == 2);
    CMRdbgMsg(6, "Terminal members are %d and %d.\n", reducedComponent->terminalMember[0],
      reducedComponent->terminalMember[1]);
    CMRdbgMsg(6, "Terminal nodes are %d and %d.\n", reducedComponent->terminalNode[0],
      reducedComponent->terminalNode[1]);
    assert(findMember(dec, reducedComponent->terminalMember[0])
      == findMember(dec, reducedComponent->terminalMember[1]));

    /* Create new edge for this component. If there is one component, this is a column edge, and otherwise it is a
     * marker edge that will be linked to a new series member consisting of all these marker edges and the column edge.
     */
    DEC_EDGE newEdge;
    CMR_CALL( createEdge(dec, -1, &newEdge) );
    componentNewEdges[i] = newEdge;
    dec->edges[newEdge].childMember = -1;
    dec->edges[newEdge].member = findMember(dec, reducedComponent->terminalMember[0]);
    dec->edges[newEdge].head = reducedComponent->terminalNode[0];
    dec->edges[newEdge].tail = reducedComponent->terminalNode[1];
    dec->edges[newEdge].element = 0;
    CMR_CALL( addEdgeToMembersEdgeList(dec, newEdge) );
  }

  for (size_t c = 0; c < newcolumn->numReducedComponents; ++c)
  {
    if (c == maxDepthComponent)
      CMRdbgMsg(6, "Reduced component %d has maximum depth and will remain a root.\n", c);
    else
      CMR_CALL( reorderComponent(dec, findMember(dec, dec->edges[componentNewEdges[c]].member)) );
  }

  if (newcolumn->numReducedComponents == 0)
  {
    DEC_MEMBER loopMember;
    CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_LOOP, &loopMember) );
    DEC_MEMBER loopEdge;
    CMR_CALL( createEdge(dec, loopMember, &loopEdge) );
    CMR_CALL( addEdgeToMembersEdgeList(dec, loopEdge) );
    dec->edges[loopEdge].element = CMRcolumnToElement(column);
    dec->edges[loopEdge].childMember = -1;
  }
  else if (newcolumn->numReducedComponents == 1)
  {
    DEC_EDGE columnEdge = componentNewEdges[0];
    dec->edges[columnEdge].element = CMRcolumnToElement(column);
    dec->edges[columnEdge].childMember = -1;
  }
  else
  {
    /* We create another edge for the column as well as a series member containing all new edges and this one. */

    DEC_MEMBER series;
    CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_SERIES, &series) );

    DEC_EDGE columnEdge;
    CMR_CALL( createEdge(dec, series, &columnEdge) );
    dec->edges[columnEdge].childMember = -1;
    dec->edges[columnEdge].head = -1;
    dec->edges[columnEdge].tail = -1;
    dec->edges[columnEdge].element = CMRcolumnToElement(column);
    CMR_CALL( addEdgeToMembersEdgeList(dec, columnEdge) );

    for (size_t i = 0; i < newcolumn->numReducedComponents; ++i)
    {
      DEC_EDGE newEdge = componentNewEdges[i];
      DEC_EDGE markerEdge;
      CMR_CALL( createEdge(dec, series, &markerEdge) );
      CMR_CALL( addEdgeToMembersEdgeList(dec, markerEdge) );
      dec->edges[markerEdge].head = -1;
      dec->edges[markerEdge].tail = -1;

      DEC_MEMBER partnerMember = findEdgeMember(dec, newEdge);

      if (i == maxDepthComponent)
      {
        dec->edges[markerEdge].childMember = -1;
        dec->edges[markerEdge].element = INT_MAX - dec->numMarkerPairs;
        dec->members[series].parentMember = partnerMember;
        dec->members[series].markerToParent = markerEdge;
        dec->members[series].markerOfParent = newEdge;
        dec->edges[newEdge].element = -INT_MAX + dec->numMarkerPairs;
        dec->edges[newEdge].childMember = series;
      }
      else
      {
        dec->edges[markerEdge].childMember = partnerMember;
        dec->members[partnerMember].markerOfParent = markerEdge;
        dec->members[partnerMember].markerToParent = newEdge;
        dec->members[partnerMember].parentMember = series;
        dec->edges[markerEdge].element = INT_MAX - dec->numMarkerPairs;
        dec->edges[newEdge].element = -INT_MAX + dec->numMarkerPairs;
      }

      dec->numMarkerPairs++;
    }
  }

  debugDot(dec, newcolumn);

  CMR_CALL( CMRfreeStackArray(dec->cmr, &componentNewEdges) );

  newcolumn->numReducedMembers = 0;
  newcolumn->numReducedComponents = 0;

#if defined(CMR_DEBUG_CONSISTENCY)
  CMRconsistencyAssert( decConsistency(dec) );
#endif /* CMR_DEBUG_CONSISTENCY */

  return CMR_OKAY;
}

static
CMR_ERROR cographicnessTest(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,       /**< Some matrix to be tested for cographicness. */
  void* data,               /**< Additional data (must be \c NULL). */
  bool* pisCographic,       /**< Pointer for storing whether \p matrix is cographic. */
  CMR_SUBMAT** psubmatrix,  /**< Pointer for storing a proper non-cographic submatrix of \p matrix. */
  double timeLimit          /**< Time limit to impose. */
)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(matrix);
  assert(!data);
  assert(pisCographic);
  assert(!psubmatrix || !*psubmatrix);

  /* TODO: Algorithm for finding a bad minor. */
  CMR_UNUSED(psubmatrix);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "cographicnessTest called for a %dx%d matrix\n", matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', false);
#endif /* CMR_DEBUG */

  *pisCographic = true;
  Dec* dec = NULL;
  clock_t time = clock();
  if (matrix->numNonzeros > 0)
  {
    CMR_CALL( decCreate(cmr, &dec, 4096, 1024, 256, 256, 256) );

    /* Process each column. */
    DEC_NEWCOLUMN* newcolumn = NULL;
    CMR_CALL( newcolumnCreate(cmr, &newcolumn) );
    size_t columnTimeFactor = matrix->numRows / 100 + 1;
    for (size_t column = 0; column < matrix->numRows && *pisCographic; ++column)
    {
      if ((column % columnTimeFactor == 0) && (clock() - time) * 1.0 / CLOCKS_PER_SEC > timeLimit)
      {
        CMR_CALL( newcolumnFree(cmr, &newcolumn) );
        if (dec)
          CMR_CALL( decFree(&dec) );
        return CMR_ERROR_TIMEOUT;
      }
      
      CMR_CALL( addColumnCheck(dec, newcolumn, &matrix->entryColumns[matrix->rowSlice[column]],
        matrix->rowSlice[column+1] - matrix->rowSlice[column]) );

      if (newcolumn->remainsGraphic)
      {
        CMR_CALL( addColumnApply(dec, newcolumn, column, &matrix->entryColumns[matrix->rowSlice[column]],
          matrix->rowSlice[column+1] - matrix->rowSlice[column]) );
      }
      else
        *pisCographic = false;
    }

    CMR_CALL( newcolumnFree(cmr, &newcolumn) );
  }

  if (dec)
    CMR_CALL( decFree(&dec) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestCographicMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisCographic, CMR_GRAPH** pgraph,
  CMR_GRAPH_EDGE** pforestEdges, CMR_GRAPH_EDGE** pcoforestEdges, CMR_SUBMAT** psubmatrix,
  CMR_GRAPHIC_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pforestEdges || pgraph);
  assert(!pcoforestEdges || pgraph);
  assert(pisCographic);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRtestCographicMatrix called for a %dx%d matrix\n", matrix->numRows, matrix->numColumns);
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  clock_t time = clock();
  *pisCographic = true;

  Dec* dec = NULL;
  if (matrix->numNonzeros > 0)
  {
    CMR_CALL( decCreate(cmr, &dec, 4096, 1024, 256, 256, 256) );

    /* Process each column. */
    DEC_NEWCOLUMN* newcolumn = NULL;
    CMR_CALL( newcolumnCreate(cmr, &newcolumn) );
    for (size_t column = 0; column < matrix->numRows && *pisCographic; ++column)
    {
      clock_t checkClock = clock();
      double remainingTime = timeLimit - (checkClock - time) * 1.0 / CLOCKS_PER_SEC;
      if (remainingTime < 0)
      {
        CMR_CALL( newcolumnFree(cmr, &newcolumn) );
        if (dec)
          CMR_CALL( decFree(&dec) );
        return CMR_ERROR_TIMEOUT;
      }
      CMR_CALL( addColumnCheck(dec, newcolumn, &matrix->entryColumns[matrix->rowSlice[column]],
        matrix->rowSlice[column+1] - matrix->rowSlice[column]) );
      if (stats)
      {
        stats->checkCount++;
        stats->checkTime += (clock() - checkClock) * 1.0 / CLOCKS_PER_SEC;
      }

      debugDot(dec, newcolumn);

      if (newcolumn->remainsGraphic)
      {
        clock_t applyClock = (stats ? clock() : 0);

        CMR_CALL( addColumnApply(dec, newcolumn, column, &matrix->entryColumns[matrix->rowSlice[column]],
          matrix->rowSlice[column+1] - matrix->rowSlice[column]) );

        if (stats)
        {
          stats->applyCount++;
          stats->applyTime += (clock() - applyClock) * 1.0 / CLOCKS_PER_SEC;
        }
      }
      else
        *pisCographic = false;
    }

    CMR_CALL( newcolumnFree(cmr, &newcolumn) );
  }

  if (*pisCographic)
  {
    /* Allocate memory for graph, forest and coforest. */

    CMR_GRAPH* graph = NULL;
    if (pgraph)
    {
      if (!*pgraph)
      {
        CMR_CALL( CMRgraphCreateEmpty(cmr, pgraph, matrix->numColumns + 2 * matrix->numRows,
          matrix->numColumns + 3 * matrix->numRows) );
      }
      graph = *pgraph;
    }

    int* forest = NULL;
    if (pforestEdges)
    {
      if (!*pforestEdges)
        CMR_CALL( CMRallocBlockArray(cmr, pforestEdges, matrix->numColumns) );
      forest = *pforestEdges;
    }
    int* coforest = NULL;
    if (pcoforestEdges)
    {
      if (!*pcoforestEdges)
        CMR_CALL( CMRallocBlockArray(cmr, pcoforestEdges, matrix->numRows) );
      coforest = *pcoforestEdges;
    }

    if (graph)
    {
      if (matrix->numNonzeros > 0)
      {
        /* Add members and edges for empty rows. */
        if (dec->numRows < matrix->numColumns)
        {
          /* Reallocate if necessary. */
          if (dec->memRows < matrix->numColumns)
          {
            CMRreallocBlockArray(cmr, &dec->rowEdges, matrix->numColumns);
            dec->memRows = matrix->numColumns;
          }

          /* Add single-edge parallel for each missing row. */
          for (size_t r = dec->numRows; r < matrix->numColumns; ++r)
          {
            DEC_MEMBER member;
            CMR_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &member) );

            DEC_EDGE edge;
            CMR_CALL( createEdge(dec, member, &edge) );
            CMR_CALL( addEdgeToMembersEdgeList(dec, edge) );
            dec->edges[edge].element = CMRrowToElement(r);
            dec->edges[edge].head = -1;
            dec->edges[edge].tail = -1;
            dec->edges[edge].childMember = -1;

            CMRdbgMsg(8, "New empty row %d is edge %d of member %d.\n", r, edge, member);

            dec->rowEdges[r].edge = edge;
          }

          dec->numRows = matrix->numColumns;
        }

        CMR_CALL( decToGraph(dec, graph, true, forest, coforest, NULL) );
      }
      else
      {
        CMR_CALL( CMRgraphClear(cmr, graph) );
        /* Construct a path with numRows edges and with numColumns loops at 0. */

        CMR_GRAPH_NODE s;
        CMR_CALL( CMRgraphAddNode(cmr, graph, &s) );
        for (size_t c = 0; c < matrix->numRows; ++c)
        {
          CMR_GRAPH_EDGE e;
          CMR_CALL( CMRgraphAddEdge(cmr, graph, s, s, &e) );
          if (coforest)
            *coforest++ = e;
        }
        for (size_t r = 0; r < matrix->numColumns; ++r)
        {
          CMR_GRAPH_NODE t;
          CMR_CALL( CMRgraphAddNode(cmr, graph, &t) );
          CMR_GRAPH_EDGE e;
          CMR_CALL( CMRgraphAddEdge(cmr, graph, s, t, &e) );
          if (forest)
            *forest++ = e;
          s = t;
        }

        CMRdbgMsg(0, "Constructed graph with %d nodes and %d edges.\n", CMRgraphNumNodes(graph), CMRgraphNumEdges(graph));
      }
    }
  }

  if (dec)
    CMR_CALL( decFree(&dec) );

  if (!*pisCographic && psubmatrix)
  {
    /* Find submatrix. */
    double remainingTime = timeLimit - (clock() - time) * 1.0 / CLOCKS_PER_SEC;
    CMR_CALL( CMRtestHereditaryPropertySimple(cmr, matrix, cographicnessTest, NULL, psubmatrix, remainingTime) );
  }

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - time) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRtestBinaryGraphicColumnSubmatrixGreedy(CMR* cmr, CMR_CHRMAT* transpose, size_t* orderedColumns,
  CMR_SUBMAT** psubmatrix)
{
  assert(cmr);
  assert(transpose);
  assert(psubmatrix);

  size_t numRows = transpose->numColumns;
  size_t numColumns = transpose->numRows;

  CMRdbgMsg(0, "CMRtestBinaryGraphicColumnSubmatrixGreedy for %dx%d matrix with transpose\n", numRows, numColumns);
#if defined(CMR_DEBUG)
  CMR_CALL( CMRchrmatPrintDense(cmr, transpose, stdout, '0', true) );
#endif /* CMR_DEBUG */

  CMR_CALL( CMRsubmatCreate(cmr, numRows, numColumns, psubmatrix) );
  CMR_SUBMAT* submatrix = *psubmatrix;
  submatrix->numRows = 0;
  submatrix->numColumns = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    submatrix->rows[submatrix->numRows++] = row;
  }

  /* Try to add each column. */
  Dec* dec = NULL;
  CMR_CALL( decCreate(cmr, &dec, 4096, 1024, 256, 256, 256) );

  /* Process each column. */
  DEC_NEWCOLUMN* newcolumn = NULL;
  CMR_CALL( newcolumnCreate(cmr, &newcolumn) );
  for (size_t c = 0; c < numColumns; ++c)
  {
    size_t column = orderedColumns ? orderedColumns[c] : c;

    CMRdbgMsg(0, "!!! Trying to append column %d.\n", column);

    int lengthColumn = transpose->rowSlice[column+1] - transpose->rowSlice[column];
    CMR_CALL( addColumnCheck(dec, newcolumn, &transpose->entryColumns[transpose->rowSlice[column]], lengthColumn) );

    CMRdbgMsg(0, "!!! Appending column %s graphicness.\n", newcolumn->remainsGraphic ? "maintains" : " would destroy");

    if (newcolumn->remainsGraphic)
    {
      CMR_CALL( addColumnApply(dec, newcolumn, column, &transpose->entryColumns[transpose->rowSlice[column]],
        lengthColumn) );
      submatrix->columns[submatrix->numColumns] = column;
      submatrix->numColumns++;
    }
  }

  CMR_CALL( newcolumnFree(cmr, &newcolumn) );
  CMR_CALL( decFree(&dec) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestGraphicMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisGraphic, CMR_GRAPH** pgraph,
  CMR_GRAPH_EDGE** pforestEdges, CMR_GRAPH_EDGE** pcoforestEdges, CMR_SUBMAT** psubmatrix,
  CMR_GRAPHIC_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pforestEdges || pgraph);
  assert(!pcoforestEdges || pgraph);
  assert(pisGraphic);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRtestGraphicMatrix called for a %dx%d matrix \n", matrix->numRows,
    matrix->numColumns);
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  /* Create transpose of matrix. */
  clock_t transposeClock;
  double transposeTime;
  if (stats)
    transposeClock = clock();
  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

  if (stats)
  {
    stats->transposeCount++;
    transposeTime = (clock() - transposeClock) * 1.0 / CLOCKS_PER_SEC;
    stats->transposeTime += transposeTime;
  }

  CMR_CALL( CMRtestCographicMatrix(cmr, transpose, pisGraphic, pgraph, pforestEdges, pcoforestEdges, psubmatrix,
    stats, timeLimit) );

  if (stats)
    stats->totalTime += transposeTime;

  /* Transpose minimal non-cographic matrix to become a minimal non-graphic matrix. */
  if (psubmatrix && *psubmatrix)
    CMR_CALL( CMRsubmatTranspose(*psubmatrix) );

  CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}

/**@}*/
