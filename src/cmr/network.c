// #define CMR_DEBUG /* Uncomment to debug network. */

#include <cmr/graphic.h>
#include <cmr/network.h>
#include <cmr/camion.h>

#if defined(CMR_DEBUG)
#include <cmr/linear_algebra.h>
#endif /* CMR_DEBUG */

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
  CMR_GRAPH* graph = NULL;
  bool isConetwork;
  CMR_CALL( CMRgraphicTestTranspose(cmr, matrix, &isConetwork, &graph, &forestEdges, &coforestEdges,
    psubmatrix, stats ? &stats->graphic : NULL, remainingTime) );

#if defined(CMR_DEBUG)
  CMRdbgMsg(2, "CMRtestCographicMatrix() returned %s.\n", isConetwork ? "TRUE": "FALSE");
#endif /* CMR_DEBUG */

  bool* arcsReversed = NULL;
  if (isConetwork)
  {
    /* We have to find out which edges are reversed. */
    CMRdbgMsg(0, "Matrix is graphic. Trying to compute reversed edges.");
    CMR_CALL( CMRallocBlockArray(cmr, &arcsReversed, CMRgraphMemEdges(graph)) );

#if defined(CMR_DEBUG)
    CMRgraphPrint(graph, stdout);
    for (size_t b = 0; b < matrix->numColumns; ++b)
      CMRdbgMsg(2, "Forest #%zu is %d.\n", b, forestEdges[b]);
    for (size_t b = 0; b < matrix->numRows; ++b)
      CMRdbgMsg(2, "Coforest #%zu is %d.\n", b, coforestEdges[b]);
#endif /* CMR_DEBUG */

    CMR_CALL( CMRcamionCographicOrient(cmr, matrix, graph, forestEdges, coforestEdges, arcsReversed, &isConetwork,
      psubmatrix, stats ? &stats->camion : NULL) );
  }

#if defined(CMR_DEBUG)
  if (psubmatrix && *psubmatrix)
  {
    CMR_CHRMAT* submatrix = NULL;
    CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, *psubmatrix, &submatrix) );
    CMR_CALL( CMRchrmatPrintDense(cmr, submatrix, stdout, '0', true) );

    int64_t determinant;
    CMR_CALL( CMRchrmatDeterminant(cmr, submatrix, &determinant) );
    CMRdbgMsg(2, "-> Returned submatrix has determinant %ld.\n", determinant);

    CMR_CALL( CMRchrmatFree(cmr, &submatrix) );
  }
#endif /* CMR_DEBUG */

  if (pisConetwork)
    *pisConetwork = isConetwork;
  if (isConetwork && pdigraph)
    *pdigraph = graph;
  else
    CMR_CALL( CMRgraphFree(cmr, &graph) );
  if (isConetwork && pforestArcs)
    *pforestArcs = forestEdges;
  else
    CMR_CALL( CMRfreeBlockArray(cmr, &forestEdges) );
  if (isConetwork && pcoforestArcs)
    *pcoforestArcs = coforestEdges;
  else
    CMR_CALL( CMRfreeBlockArray(cmr, &coforestEdges) );
  if (isConetwork && parcsReversed)
    *parcsReversed = arcsReversed;
  else
    CMR_CALL( CMRfreeBlockArray(cmr, &arcsReversed) );

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
