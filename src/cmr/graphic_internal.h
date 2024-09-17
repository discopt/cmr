#ifndef CMR_GRAPHIC_INTERNAL_H
#define CMR_GRAPHIC_INTERNAL_H

#include <cmr/graphic.h>

/**
 * \brief Computes the network or graphic matrix of a given (di)graph \f$ D = (V,A) \f$.
 *
 * Computes the [network matrix](\ref network) \f$ M := M(D,T) \f$ for given \f$ D \f$ and optionally given (directed)
 * spanning forest \f$ T \subseteq A \f$ or the support matrix of \f$ M(D,T) \f$.
 * If \f$ T \f$ is not given, an arbitrary (directed) spanning forest of \f$ D \f$ is used.
 * The direction of the edges is that of \p digraph, but may be flipped by specifying \p arcsReversed.
 * If \p forestArcs is \c NULL, an arbitrary (directed) spanning forest \f$ T \f$ of \f$ D \f$ is computed.
 * The ordering of the columns can be specified via \p coforestArcs.
 *
 * \note The function computes a network matrix of \f$ D \f$ (and \f$ T \f$) regardless of whether \p forestArcs is
 * a correct (directed) spanning forest. Whether this was the case is indicated via \p *pisCorrectForest.
 */

CMR_ERROR CMRcomputeRepresentationMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_GRAPH* digraph,             /**< Digraph \f$ D = (V,A) \f$. */
  bool ternary,                   /**< Whether we need to compute correct signs. */
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
 * \brief Tests the support matrix \f$ M \f$ of the input matrix for being a [cographic matrix](\ref graphic).
 *
 * Tests if \f$ M = M(G,T)^{\mathsf{T}} \f$ for some graph \f$ G = (V,E) \f$ and some spanning forest
 * \f$ T \subseteq E \f$ of \f$ G \f$ and sets \p *pisCographic accordingly.
 *
 * If \f$ M \f$ is a cographic matrix and \p pgraph != \c NULL, then one possible graph \f$ G \f$ is computed and
 * stored in \p *pgraph. The caller must release its memory via \ref CMRgraphFree.
 * If in addition to \p pgraph also \p pforestEdges != \c NULL (resp. \p pcoforestEdges != \c NULL), then a
 * corresponding spanning forest \f$ T \f$ (resp.\ its complement \f$ E \setminus T \f$) is stored in
 * \p *pforestEdges (resp. \p *pcoforestEdges). The caller must release this memory via \ref CMRfreeBlockArray.
 *
 * \note Retrieval of minimal non-cographic submatrices via \p *psubmatrix is not implemented, yet.
 */

CMR_ERROR CMRcographicTestSupport(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Matrix \f$ M \f$ */
  bool* pisCographic,               /**< Returns true if and only if \f$ M \f$ is a cographic matrix. */
  CMR_GRAPH** pgraph,               /**< Pointer for storing the graph \f$ G \f$ (if \f$ M \f$ is graphic). */
  CMR_GRAPH_EDGE** pforestEdges,    /**< Pointer for storing \f$ T \f$, indexed by the rows of \f$ M \f$ (if \f$ M \f$
                                     **  is graphic).  */
  CMR_GRAPH_EDGE** pcoforestEdges,  /**< Pointer for storing \f$ E \setminus T \f$, indexed by the columns of \f$ M \f$
                                     **  (if \f$ M \f$ is graphic). */
  CMR_SUBMAT** psubmatrix,          /**< Pointer for storing a minimal non-graphic submatrix (if \f$ M \f$ is not
                                     **  graphic). */
  CMR_GRAPHIC_STATISTICS* stats,    /**< Pointer to statistics (may be \c NULL). */
  double timeLimit                  /**< Time limit to impose. */
);

#endif /* CMR_GRAPHIC_INTERNAL_H */
