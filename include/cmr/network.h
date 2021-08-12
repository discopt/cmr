#ifndef CMR_NETWORK_H
#define CMR_NETWORK_H

#include <cmr/env.h>
#include <cmr/element.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Computes the network matrix for a given digraph.
 *
 * Let \f$ D = (V,A) \f$ be a directed graph with nodes \f$ V \f$ and arc \f$ A \f$ and let
 * \f$ T \subseteq E \f$ be an (arbitrarily) directed spanning forest of \f$ D \f$.
 * The **network matrix** \f$ M := M(D,T) \f$ is a matrix
 * \f$ M \in \{-1,0,+1\}^{T \times (A \setminus T)} \f$ with \f$ M_{a,\{v,w\}} = +1 \f$ (resp.\ \f$ M_{a,\{v,w\}} \f$ if
 * \f$ a \f$ is a forward edge (resp.\ backward edge) on the unique \f$ v \f$-\f$ w \f$-path in \f$ T \f$, and
 * \f$ M_{a,\{v,w\}} = 0 \f$ otherwise.
 *
 * Computes \f$ M(D,T) \f$ for given \f$ D \f$ and spanning forest \f$ T \f$ given by \p forestEdges.
 * The direction of the edges is that of \p digraph, but may be modified by specifying \p edgesReversed.
 * If \p forestEdges is \c NULL, an arbitrary spanning forest \f$ T \f$ of \f$ G \f$ is computed.
 * The ordering of the columns can be specified via \p coforestEdges.
 *
 * \note The function computes a network matrix of \f$ D \f$ regardless of whether \p forestEdges is a correct
 * (directed) spanning tree. This is indicated via *\p pisCorrectForest.
 */

CMR_EXPORT
CMR_ERROR CMRcomputeNetworkMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_GRAPH* digraph,             /**< \ref Graph \f$ G \f$. */
  CMR_CHRMAT** pmatrix,           /**< Pointer for storing \f$ M \f$ (may be \c NULL). */
  CMR_CHRMAT** ptranspose,        /**< Pointer for storing \f$ M^{\mathsf{T}} \f$ (may be \c NULL). */
  bool* edgesReversed,            /**< Indicates, for each edge \f$ \{u, v\}\f$, whether we consider \f$ (u, v)\f$  (if \c false)
                                   **  or \f$ (v,u)\f$  (if \c true). */
  int numForestEdges,             /**< Length of \p forestEdges (0 if \c forestEdges is \c NULL). */
  CMR_GRAPH_EDGE* forestEdges,    /**< If not \c NULL, spanning forest edges as rows in this order. */
  int numCoforestEdges,           /**< Length of \p coforestEdges (0 if \c coforestEdges is \c NULL). */
  CMR_GRAPH_EDGE* coforestEdges,  /**< If not \c NULL, complement of forest edges as columns in this order. */
  bool* pisCorrectForest          /**< Pointer for storing whether \c forestEdges is a spanning forest of \f$ G \f$ (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ A \f$ for being a network matrix.
 *
 * Tests if \f$ A = M(D,T) \f$ for some digraph \f$ D \f$ and some (directed) spanning forest \f$ T \f$ of D and sets
 * \p pisGraphic accordingly.
 * The matrix \f$ A \f$ is given by \f$ A^{\textsf{T}} := \f$ \p transpose.
 *
 * If \f$ A \f$ is such a representation matrix and \p pgraph != \c NULL, then one possible digraph \f$ D \f$ is
 * computed and stored in *\p pgraph.
 * The caller must release the memory via \ref CMRgraphFree.
 * If in addition to \p pgraph also \p pforestEdges (resp. \p pcoforestEdges) != \c NULL, then a corresponding
 * (directed) spanning forest \f$ T \f$ (resp.\ its complement \f$ A \setminus T \f$ is stored in \p pforestEdges
 * (resp. \p pcoforestEdges).
 * The caller must release the memory via \ref CMRfreeBlockArray.
 */

CMR_EXPORT
CMR_ERROR CMRtestNetworkMatrix(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* transpose,            /**< \f$ M^{\mathsf{T}} \f$ */
  bool* pisGraphic,                 /**< Returns true if and only if the matrix is a network matrix. */
  CMR_GRAPH** pdigraph,             /**< Pointer for storing the digraph \f$ D \f$ (if network). */
  CMR_GRAPH_EDGE** pforestEdges,    /**< Pointer for storing \f$ T \f$ (if network).  */
  CMR_GRAPH_EDGE** pcoforestEdges,  /**< Pointer for storing \f$ A \setminus T \f$ (if network). */
  bool** pedgesReversed,            /**< Pointer for storing indicators which edges are reversed for the correct sign. */
  CMR_SUBMAT** psubmatrix           /**< Pointer for storing a minimal nongraphic submatrix (if not network);
                                     **  This is not implemented, yet. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_NETWORK_H */
