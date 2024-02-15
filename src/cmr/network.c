// #define CMR_DEBUG /* Uncomment to debug network. */

#include <cmr/graphic.h>
#include <cmr/network.h>
#include <cmr/camion.h>

#include "graphic_internal.h"
#include "env_internal.h"
#include "matrix_internal.h"
#include "block_decomposition.h"
#include "heap.h"
#include "sort.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

CMR_ERROR CMRnetworkStatsInit(CMR_NETWORK_STATISTICS* stats)
{
  assert(stats);

  stats->totalCount = 0;
  stats->totalTime = 0.0;
  CMR_CALL( CMRcamionStatsInit(&stats->camion) );
  CMR_CALL( CMRgraphicStatsInit(&stats->graphic) );

  return CMR_OKAY;
}

CMR_ERROR CMRnetworkStatsPrint(FILE* stream, CMR_NETWORK_STATISTICS* stats, const char* prefix)
{
  assert(stream);
  assert(stats);

  if (!prefix)
  {
    fprintf(stream, "Network matrix recognition:\n");
    prefix = "  ";
  }

  char subPrefix[256];
  snprintf(subPrefix, 256, "%scamion ", prefix);
  CMR_CALL( CMRcamionStatsPrint(stream, &stats->camion, subPrefix) );
  snprintf(subPrefix, 256, "%sgraphic ", prefix);
  CMR_CALL( CMRgraphicStatsPrint(stream, &stats->graphic, subPrefix) );

  fprintf(stream, "%stotal: %lu in %f seconds\n", prefix, (unsigned long)stats->totalCount,
    stats->totalTime);

  return CMR_OKAY;
}

typedef enum
{
  UNKNOWN = 0,    /**< \brief The node was not considered by the shortest-path, yet. */
  SEEN = 1,       /**< \brief Some path to the node is known. */
  COMPLETED = 2,  /**< \brief The shortest path to the node is known. */
  BASIC = 3,      /**< \brief The rootEdge of that node belongs to the spanning forest. */
} DIJKSTRA_STAGE;

CMR_ERROR CMRnetworkComputeMatrix(CMR* cmr, CMR_GRAPH* digraph, CMR_CHRMAT** pmatrix, CMR_CHRMAT** ptranspose,
  bool* arcsReversed, int numForestArcs, CMR_GRAPH_EDGE* forestArcs, int numCoforestArcs,
  CMR_GRAPH_EDGE* coforestArcs, bool* pisCorrectForest)
{
  assert(cmr);
  assert(digraph);
  assert(pmatrix || ptranspose);
  assert(!pmatrix || !*pmatrix);
  assert(!ptranspose || !*ptranspose);

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRcomputeRepresentationMatrix(cmr, digraph, true, &transpose, arcsReversed, numForestArcs, forestArcs,
    numCoforestArcs, coforestArcs, pisCorrectForest) );

  CMRconsistencyAssert( CMRchrmatConsistency(transpose) );

  if (pmatrix)
    CMR_CALL( CMRchrmatTranspose(cmr, transpose, pmatrix) );

  /* Return or free the transpose matrix. */
  if (ptranspose)
    *ptranspose = transpose;
  else
    CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}


typedef struct
{
  int forestIndex;
} NetworkEdgeData;

typedef struct
{
  DIJKSTRA_STAGE stage; /* Stage in BFS. */
  int predecessor;      /* Predecessor node. */
  CMR_GRAPH_EDGE edge;   /* Edge connecting to predecessor node. */
  int distance;         /* Combinatorial distance to the BFS root. */
  char sign;            /* Sign of this tree edge with respect to current column. */
  bool fixed;           /* Whether the orientation of this edge is already fixed. */
} NetworkNodeData;

CMR_ERROR CMRnetworkTestTranspose(CMR* cmr, CMR_CHRMAT* matrix, bool* pisConetwork, CMR_GRAPH** pdigraph,
  CMR_GRAPH_EDGE** pforestArcs, CMR_GRAPH_EDGE** pcoforestArcs, bool** parcsReversed, CMR_SUBMAT** psubmatrix,
  CMR_NETWORK_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pforestArcs || pdigraph);
  assert(!pcoforestArcs || pdigraph);
  assert(!parcsReversed || pdigraph);
  assert(pisConetwork);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRnetworkTestTranspose called for a %dx%d matrix \n", matrix->numRows,
    matrix->numColumns);
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  clock_t totalClock = clock();

  double remainingTime = timeLimit - (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  CMR_GRAPH_EDGE* forestEdges = NULL;
  CMR_GRAPH_EDGE* coforestEdges = NULL;
  CMR_CALL( CMRgraphicTestTranspose(cmr, matrix, pisConetwork, pdigraph, pforestArcs ? &forestEdges : NULL,
    pcoforestArcs ? &coforestEdges : NULL, psubmatrix, stats ? &stats->graphic : NULL, remainingTime) );
#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "CMRtestCographicMatrix() returned %s.\n", (*pisConetwork) ? "TRUE": "FALSE");
#endif /* CMR_DEBUG */

  if (pforestArcs)
    *pforestArcs = forestEdges;
  if (pcoforestArcs)
    *pcoforestArcs = coforestEdges;
  if (!*pisConetwork || !pdigraph || !parcsReversed)
  {
    /* We have to free (co)forest information if the caller didn't ask for it. */
    if (!pforestArcs)
      CMR_CALL( CMRfreeBlockArray(cmr, &forestEdges) );
    if (!pcoforestArcs)
      CMR_CALL( CMRfreeBlockArray(cmr, &coforestEdges) );
    return CMR_OKAY;
  }

  /* We have to find out which edges are reversed. */
  CMR_GRAPH* graph = *pdigraph;
  CMRdbgMsg(0, "Matrix is graphic. Trying to compute reversed edges.");
  CMR_CALL( CMRallocBlockArray(cmr, parcsReversed, CMRgraphMemEdges(graph)) );

#if defined(CMR_DEBUG)
  CMRgraphPrint(stdout, *pdigraph);
  for (size_t b = 0; b < matrix->numColumns; ++b)
    CMRdbgMsg(2, "Forest #%zu is %d.\n", b, (*pforestArcs)[b]);
  for (size_t b = 0; b < matrix->numRows; ++b)
    CMRdbgMsg(2, "Coforest #%zu is %d.\n", b, (*pcoforestArcs)[b]);
#endif /* CMR_DEBUG */

  CMR_CALL( CMRcamionCographicOrient(cmr, matrix, *pdigraph, forestEdges, coforestEdges, *parcsReversed, pisConetwork,
    psubmatrix, stats ? &stats->camion : NULL) );

  /* We have to free (co)forest information if the caller didn't ask for it. */
  if (!pforestArcs)
    CMR_CALL( CMRfreeBlockArray(cmr, &forestEdges) );
  if (!pcoforestArcs)
    CMR_CALL( CMRfreeBlockArray(cmr, &coforestEdges) );

  if (stats)
  {
    stats->totalCount++;
    stats->totalTime += (clock() - totalClock) * 1.0 / CLOCKS_PER_SEC;
  }

  return CMR_OKAY;
}

CMR_ERROR CMRnetworkTestMatrix(CMR* cmr, CMR_CHRMAT* matrix, bool* pisNetwork, CMR_GRAPH** pdigraph,
  CMR_GRAPH_EDGE** pforestArcs, CMR_GRAPH_EDGE** pcoforestArcs, bool** parcsReversed, CMR_SUBMAT** psubmatrix,
  CMR_NETWORK_STATISTICS* stats, double timeLimit)
{
  assert(cmr);
  assert(matrix);
  assert(pisNetwork);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pforestArcs || pdigraph);
  assert(!pcoforestArcs || pdigraph);
  assert(!parcsReversed || pdigraph);

  /* Create transpose of matrix. */
  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
  
#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRnetworkTestMatrix called for a %dx%d matrix \n", matrix->numRows,
    matrix->numColumns);
  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
#endif /* CMR_DEBUG */

  CMR_CALL( CMRnetworkTestTranspose(cmr, transpose, pisNetwork, pdigraph, pforestArcs, pcoforestArcs, parcsReversed,
    psubmatrix, stats, timeLimit) );

  /* Transpose minimal non-conetwork matrix to become a minimal non-network matrix. */
  if (psubmatrix && *psubmatrix)
  {
    size_t* temp = (*psubmatrix)->rows;
    (*psubmatrix)->rows = (*psubmatrix)->columns;
    (*psubmatrix)->columns = temp;

    size_t n = (*psubmatrix)->numRows;
    (*psubmatrix)->numRows = (*psubmatrix)->numColumns;
    (*psubmatrix)->numColumns = n;
  }

  CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}

/**@}*/
