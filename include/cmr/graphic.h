#ifndef CMR_GRAPHIC_H
#define CMR_GRAPHIC_H

/**
 * \file graphic.h
 *
 * \author Matthias Walter
 *
 * \brief Computation and recognition of [graphic matrices](\ref graphic) and [cographic matrices](\ref graphic).
 *
 * The following notation is used throughout:
 *  - \f$ G = (V,E) \f$ for the graph with nodes \f$ V \f$ and edges \f$ E \f$.
 *  - \f$ T \subseteq E \f$ for a spanning forest of \f$ G \f$.
 *  - \f$ M \f$ for the (potential) [(co)graphic matrix](\ref graphic).
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
 * \brief Statistics for graphicness test.
 */

typedef struct
{
  size_t totalCount;      /**< Total number of invocations. */
  double totalTime;       /**< Total time of all invocations. */
  size_t checkCount;      /**< Number of calls to check algorithm. */
  double checkTime;       /**< Time of check algorithm calls. */
  size_t applyCount;      /**< Number of column additions. */
  double applyTime;       /**< Time of column additions. */
  size_t transposeCount;  /**< Number of matrix transpositions. */
  double transposeTime;   /**< Time for matrix transpositions. */
} CMR_GRAPHIC_STATISTICS;

/**
 * \brief Initializes all statistics for graphicness computations.
 */

CMR_EXPORT
CMR_ERROR CMRgraphicInitStatistics(
  CMR_GRAPHIC_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for graphicness computations.
 */

CMR_EXPORT
CMR_ERROR CMRgraphicPrintStatistics(
  FILE* stream,                 /**< File stream to print to. */
  CMR_GRAPHIC_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Computes the graphic matrix of a given graph \f$ G = (V,E) \f$.
 *
 * Computes the [graphic matrix](\ref graphic) \f$ M := M(G,T) \f$ for given \f$ G \f$ and optionally given spanning
 * forest \f$ T \subseteq E \f$.
 * If \f$ T \f$ is not given, an arbitrary spanning forest of \f$ G \f$ is used.
 * If \p forestEdges is \c NULL, an arbitrary spanning forest \f$ T \f$ of \f$ G \f$ is computed.
 * The ordering of the columns can be specified via \p coforestEdges.
 *
 * \note The function computes a graphic matrix of \f$ G \f$ (and \f$ T \f$) regardless of whether \p forestEdges is
 * a correct spanning forest. Whether this was the case is indicated via \p *pisCorrectForest.
 */

CMR_EXPORT
CMR_ERROR CMRcomputeGraphicMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_GRAPH* graph,               /**< Graph \f$ G = (V,E) \f$. */
  CMR_CHRMAT** pmatrix,           /**< Pointer for storing \f$ M \f$ (may be \c NULL). */
  CMR_CHRMAT** ptranspose,        /**< Pointer for storing \f$ M^{\mathsf{T}} \f$ (may be \c NULL). */
  int numForestEdges,             /**< \f$ |T| \f$ (0 if \c forestEdges is \c NULL). */
  CMR_GRAPH_EDGE* forestEdges,    /**< \f$ T \f$, ordered by the rows of \f$ M \f$ (may be \c NULL). */
  int numCoforestEdges,           /**< \f$ |E \setminus T| \f$ (0 if \c coforestEdges is \c NULL). */
  CMR_GRAPH_EDGE* coforestEdges,  /**< \f$ E \setminus T \f$, ordered by the columns of \f$ M \f$ (may be \c NULL). */
  bool* pisCorrectForest          /**< Pointer for storing whether \c forestEdges is a spanning forest of \f$ G \f$
                                   **  (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being a [graphic matrix](\ref graphic).
 *
 * Tests if \f$ M = M(G,T) \f$ for some graph \f$ G = (V,E) \f$ and some spanning forest \f$ T \subseteq E \f$ of
 * \f$ G \f$ and sets \p *pisGraphic accordingly.
 *
 * \note If a column-wise representation of \f$ M \f$ is available, it is recommended to call
 *       \ref CMRtestCographicMatrix() for that. In fact, the implementation explicitly constructs
 *       \f$ M^{\mathsf{T}} \f$ before calling this function.
 *
 * If \f$ M \f$ is a graphic matrix and \p pgraph != \c NULL, then one possible graph \f$ G \f$ is computed and
 * stored in \p *pgraph. The caller must release its memory via \ref CMRgraphFree.
 * If in addition to \p pgraph also \p pforestEdges != \c NULL (resp. \p pcoforestEdges != \c NULL), then a
 * corresponding spanning forest \f$ T \f$ (resp.\ its complement \f$ E \setminus T \f$) is stored in
 * \p *pforestEdges (resp. \p *pcoforestEdges). The caller must release this memory via \ref CMRfreeBlockArray.
 * 
 * \note Retrieval of minimal non-graphic submatrices via \p *psubmatrix is not implemented, yet.
 */

CMR_EXPORT
CMR_ERROR CMRtestGraphicMatrix(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Matrix \f$ M \f$. */
  bool* pisGraphic,                 /**< Pointer for storing \c true if and only if \f$ M \f$ is a graphic matrix. */
  CMR_GRAPH** pgraph,               /**< Pointer for storing the graph \f$ G \f$ (if \f$ M \f$ is graphic). */
  CMR_GRAPH_EDGE** pforestEdges,    /**< Pointer for storing \f$ T \f$, indexed by the rows of \f$ M \f$ (if \f$ M \f$
                                     **  is graphic).  */
  CMR_GRAPH_EDGE** pcoforestEdges,  /**< Pointer for storing \f$ E \setminus T \f$, indexed by the columns of \f$ M \f$
                                     **  (if \f$ M \f$ is graphic). */
  CMR_SUBMAT** psubmatrix,          /**< Pointer for storing a minimal non-graphic submatrix (if \f$ M \f$ is not
                                     **  graphic). */
  CMR_GRAPHIC_STATISTICS* stats     /**< Pointer to statistics (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being a [cographic matrix](\ref graphic).
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

CMR_EXPORT
CMR_ERROR CMRtestCographicMatrix(
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
  CMR_GRAPHIC_STATISTICS* stats     /**< Pointer to statistics (may be \c NULL). */
);

/**
 * \brief Finds an inclusion-wise maximal subset of columns that induces a graphic submatrix.
 *
 * Finds an inclusion-wise maximal subset \f$ J \f$ of columns of \f$ M \f$ such that \f$ M_{\star,J} \f$ is a graphic
 * representation matrix.
 * To achieve this, it tries to append columns in the order given by \p orderedColumns, maintaining graphicness.
 */

CMR_EXPORT
CMR_ERROR CMRtestGraphicColumnSubmatrixGreedy(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* transpose,   /**< \f$ M^{\mathsf{T}} \f$ */
  size_t* orderedColumns, /**< Permutation of column indices of \f$ M \f$. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing the submatrix. */
);

/**@}*/

#ifdef __cplusplus
}
#endif

#endif /* CMR_GRAPHIC_H */
