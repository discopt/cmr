#ifndef TU_GRAPHIC_H
#define TU_GRAPHIC_H

#include <tu/env.h>
#include <tu/matrix.h>
#include <tu/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Tests a char matrix for graphicness.
 *
 * Returns \c true if and only if \p matrix is graphic.
 *
 * If \p pgraph is not \c NULL and if \p matrix is graphic, then \c *pgraph will point to a graph
 * representing \p matrix, and it is set to \c NULL otherwise.
 * If in addition to \p pgraph also \pbasis is not \c NULL, then this will be an array of edges
 * forming the basis, i.e., indexing the rows of the matrix.
 *
 * If \p submatrix is not \c NULL and the matrix is not graphic, then a minimal nongraphic submatrix
 * will be searched, which may cause extra computational effort. In this case, \c *submatrix will
 * point to this submatrix for which the caller must use \ref TUsubmatrixFree to free memory.
 * It is set to \c NULL otherwise.
 */

TU_EXPORT
bool TUtestGraphicnessChr(
  TU* tu,                   /**< TU environment */
  TU_CHRMAT* matrix,        /**< Char matrix */
  TU_GRAPH** pgraph,        /**< If not \c NULL, the graph is stored. */
  TU_GRAPH_EDGE** pbasis,   /**< If not \c NULL, contains basis edges. */
  TU_GRAPH_EDGE** pcobasis, /**< If not \c NULL, contains cobasis edges. */
  TU_SUBMAT** psubmatrix    /**< If not \c NULL, containes a minimal nongraphic submatrix. */
);

/**
 * \brief Computes a binary matrix representing a given graph.
 */

TU_EXPORT
TU_ERROR TUconvertGraphToBinaryMatrix(
  TU* tu,                     /**< TU environment. */
  TU_GRAPH* graph,            /**< Graph. */
  TU_CHRMAT** matrix,         /**< Pointer for storing the binary representation matrix. */
  int numBasisEdges,          /**< Length of \p basisEdges (0 if \c basisEdges is NULL). */
  TU_GRAPH_EDGE* basisEdges,  /**< If not \c NULL, tries to use these edges for the basis. */
  int numCobasisEdges,        /**< Length of \p cobasisEdges (0 if \c cobasisEdges is NULL). */
  TU_GRAPH_EDGE* cobasisEdges /**< If not \c NULL, tries to order columns as specified. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPHIC_H */
