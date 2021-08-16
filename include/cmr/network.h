#ifndef CMR_NETWORK_H
#define CMR_NETWORK_H

/**
 * \file network.h
 * 
 * \author Matthias Walter
 * 
 * \brief Computation and recognition of [network matrices](\ref network) and [conetwork matrices](\ref network).
 * 
 * The following notation is used throughout:
 *  - \f$ D = (V,A) \f$ for the digraph with nodes \f$ V \f$ and arcs \f$ A \f$.
 *  - \f$ T \subseteq A \f$ for a directed spanning forest of \f$ D \f$'s underlying undirected graph.
 *  - \f$ M \f$ for the (potential) [(co)network matrix](\ref network).
 */

#include <cmr/env.h>
#include <cmr/element.h>
#include <cmr/matrix.h>
#include <cmr/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Computes the network matrix of a given digraph \f$ D = (V,A) \f$.
 *
 * Computes the [network matrix](\ref network) \f$ M := M(D,T) \f$ for given \f$ D \f$ and optionally given (directed)
 * spanning forest \f$ T \subseteq A \f$.
 * If \f$ T \f$ is not given, an arbitrary (directed) spanning forest of \f$ D \f$ is used.
 * The direction of the edges is that of \p digraph, but may be flipped by specifying \p arcsReversed.
 * If \p forestArcs is \c NULL, an arbitrary (directed) spanning forest \f$ T \f$ of \f$ D \f$ is computed.
 * The ordering of the columns can be specified via \p coforestArcs.
 *
 * \note The function computes a network matrix of \f$ D \f$ (and \f$ T \f$) regardless of whether \p forestArcs is
 * a correct (directed) spanning forest. Whether this was the case is indicated via \p *pisCorrectForest.
 */

CMR_EXPORT
CMR_ERROR CMRcomputeNetworkMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_GRAPH* digraph,             /**< Digraph \f$ D = (V,A) \f$. */
  CMR_CHRMAT** pmatrix,           /**< Pointer for storing \f$ M \f$ (may be \c NULL). */
  CMR_CHRMAT** ptranspose,        /**< Pointer for storing \f$ M^{\mathsf{T}} \f$ (may be \c NULL). */
  bool* arcsReversed,             /**< Indicates, for each edge \f$ \{u, v\}\f$, whether we consider \f$ (u, v)\f$
                                   **  (if \c false) or \f$ (v,u)\f$  (if \c true). */
  int numForestArcs,              /**< \f$ |T| \f$ (0 if \c forestArcs is \c NULL). */
  CMR_GRAPH_EDGE* forestArcs,     /**< \f$ T \f$, ordered by the rows of \f$ M \f$ (may be \c NULL). */
  int numCoforestArcs,            /**< \f$ |A \setminus T| \f$ (0 if \c coforestArcs is \c NULL). */
  CMR_GRAPH_EDGE* coforestArcs,   /**< \f$ A \setminus T \f$, ordered by the columns of \f$ M \f$ (may be
                                   **  \c NULL). */
  bool* pisCorrectForest          /**< Pointer for storing whether \c forestArcs is a (directed) spanning forest of
                                   **  \f$ D \f$'s underlying undirected graph (may be \c NULL). */
);

/**
 * \brief Computes the conetwork matrix of a given digraph \f$ D = (V,A) \f$.
 *
 * Computes the [conetwork matrix](\ref network) \f$ M := M(D,T)^{\mathsf{T}} \f$ for given \f$ D \f$ and optionally
 * given (directed) spanning forest \f$ T \subseteq A \f$.
 * If \f$ T \f$ is not given, an arbitrary (directed) spanning forest of \f$ D \f$ is used.
 * The direction of the edges is that of \p digraph, but may be flipped by specifying \p arcsReversed.
 * If \p forestArcs is \c NULL, an arbitrary (directed) spanning forest \f$ T \f$ of \f$ D \f$ is computed.
 * The ordering of the rows can be specified via \p coforestArcs.
 *
 * \note The function computes a conetwork matrix of \f$ D \f$ (and \f$ T \f$) regardless of whether \p forestArcs is
 * a correct (directed) spanning forest. Whether this was the case is indicated via \p *pisCorrectForest.
 */

CMR_EXPORT
CMR_ERROR CMRcomputeConetworkMatrix(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_GRAPH* digraph,           /**< Digraph \f$ D = (V,A) \f$. */
  CMR_CHRMAT** pmatrix,         /**< Pointer for storing \f$ M \f$ (may be \c NULL). */
  CMR_CHRMAT** ptranspose,      /**< Pointer for storing \f$ M^{\mathsf{T}} \f$ (may be \c NULL). */
  bool* arcsReversed,           /**< Indicates, for each edge \f$ \{u, v\}\f$, whether we consider \f$ (u, v)\f$
                                 **  (if \c false) or \f$ (v,u)\f$  (if \c true). */
  int numForestArcs,            /**< \f$ |T| \f$ (0 if \c forestArcs is \c NULL). */
  CMR_GRAPH_EDGE* forestArcs,   /**< \f$ T \f$, ordered by the rows of \f$ M \f$ (may be \c NULL). */
  int numCoforestArcs,          /**< \f$ |A \setminus T| \f$ (0 if \c coforestArcs is \c NULL). */
  CMR_GRAPH_EDGE* coforestArcs, /**< \f$ A \setminus T \f$, ordered by the columns of \f$ M \f$ (may be
                                 **  \c NULL). */
  bool* pisCorrectForest        /**< Pointer for storing whether \c forestArcs is a (directed) spanning forest of
                                 **  \f$ D \f$'s underlying undirected graph (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being a [network matrix](\ref network).
 *
 * Tests if \f$ M = M(D,T) \f$ for some digraph \f$ D = (V,A) \f$ and some (directed) spanning forest
 * \f$ T \subseteq A \f$ of \f$ D \f$ and sets \p *pisNetwork accordingly.
 *
 * \note If a column-wise representation of \f$ M \f$ is available, it is recommended to call
 *       \ref CMRtestConetworkMatrix() for that. In fact, the implementation explicitly constructs
 *       \f$ M^{\mathsf{T}} \f$ before calling this function.
 *
 * If \f$ M \f$ is a network matrix and \p pdigraph != \c NULL, then one possible digraph \f$ D \f$ is computed and
 * stored in \p *pdigraph. The caller must release its memory via \ref CMRgraphFree.
 * If in addition to \p pdigraph also \p pforestArcs != \c NULL (resp. \p pcoforestArcs != \c NULL), then a
 * corresponding (directed) spanning forest \f$ T \f$ (resp.\ its complement \f$ A \setminus T \f$) is stored in
 * \p *pforestArcs (resp. \p *pcoforestArcs). The caller must release this memory via \ref CMRfreeBlockArray.
 * 
 * \note Retrieval of minimal non-network submatrices via \p *psubmatrix is not implemented, yet.
 */

CMR_EXPORT
CMR_ERROR CMRtestNetworkMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix \f$ M \f$. */
  bool* pisNetwork,               /**< Pointer for storing \c true if and only if \f$ M \f$ is network. */
  CMR_GRAPH** pdigraph,           /**< Pointer for storing the digraph \f$ D \f$ (if \f$ M \f$ is network). */
  CMR_GRAPH_EDGE** pforestArcs,   /**< Pointer for storing \f$ T \f$, indexed by the rows of \f$ M \f$ (if \f$ M \f$
                                   **  is network).  */
  CMR_GRAPH_EDGE** pcoforestArcs, /**< Pointer for storing \f$ A \setminus T \f$, indexed by the columns of \f$ M \f$
                                   **  (if \f$ M \f$ is network). */
  bool** parcsReversed,           /**< Pointer for storing indicators which arcs are reversed for the correct sign (if
                                   **  \f$ M \f$ is network). */
  CMR_SUBMAT** psubmatrix         /**< Pointer for storing a minimal non-network submatrix (if \f$ M \f$ is not
                                   **  network). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being a [conetwork matrix](\ref network).
 *
 * Tests if \f$ M = M(D,T)^{\mathsf{T}} \f$ for some digraph \f$ D = (V,A) \f$ and some (directed) spanning forest
 * \f$ T \subseteq A \f$ of \f$ D \f$ and sets \p *pisConetwork accordingly.
 *
 * If \f$ M \f$ is a conetwork matrix and \p pdigraph != \c NULL, then one possible digraph \f$ D \f$ is computed and
 * stored in \p *pdigraph. The caller must release its memory via \ref CMRgraphFree.
 * If in addition to \p pdigraph also \p pforestArcs != \c NULL (resp. \p pcoforestArcs != \c NULL), then a
 * corresponding (directed) spanning forest \f$ T \f$ (resp.\ its complement \f$ A \setminus T \f$) is stored in
 * \p *pforestArcs (resp. \p *pcoforestArcs). The caller must release this memory via \ref CMRfreeBlockArray.
 *
 * \note Retrieval of minimal non-conetwork submatrices via \p *psubmatrix is not implemented, yet.
 */

CMR_EXPORT
CMR_ERROR CMRtestConetworkMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix \f$ M \f$ */
  bool* pisConetwork,             /**< Returns true if and only if the matrix is a conetwork matrix. */
  CMR_GRAPH** pdigraph,           /**< Pointer for storing \c true if and only if \f$ M \f$ is conetwork. */
  CMR_GRAPH_EDGE** pforestArcs,   /**< Pointer for storing \f$ T \f$, indexed by the columns of \f$ M \f$ (if \f$ M \f$
                                   **  is network).  */
  CMR_GRAPH_EDGE** pcoforestArcs, /**< Pointer for storing \f$ A \setminus T \f$, indexed by the rows of \f$ M \f$
                                   **  (if \f$ M \f$ is conetwork). */
  bool** parcsReversed,           /**< Pointer for storing indicators which arcs are reversed for the correct sign (if
                                   **  \f$ M \f$ is conetwork). */
  CMR_SUBMAT** psubmatrix         /**< Pointer for storing a minimal non-conetwork submatrix (if \f$ M \f$ is not
                                   **  conetwork). */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_NETWORK_H */
