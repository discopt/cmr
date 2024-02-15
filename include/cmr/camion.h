#ifndef CMR_CAMION_H
#define CMR_CAMION_H

/**
 * \file camion.h
 *
 * \author Matthias Walter
 *
 * \brief Testing whether a matrix is [Camion-signed](\ref camion).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/matrix.h>
#include <cmr/graph.h>

#include <stdint.h>

/**
 * \brief Statistics for [Camion-signing](\ref camion) algorithm.
 */

typedef struct
{
  uint32_t generalCount;  /**< Number of invocations for general matrices. */
  double generalTime;     /**< Total time for all invocations for general matrices. */
  uint32_t graphCount;    /**< Number of invocations for graphic matrices. */
  double graphTime;       /**< Total time for all invocations for graphic matrices. */
  uint32_t totalCount;    /**< Total number of invocations. */
  double totalTime;       /**< Total time for all invocations. */
} CMR_CAMION_STATISTICS;

/**
 * \brief Initializes all statistics for [Camion-signing](\ref camion) algorithm.
 */

CMR_EXPORT
CMR_ERROR CMRcamionStatsInit(
  CMR_CAMION_STATISTICS* stats  /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for [Camion-signing](\ref camion) algorithm.
 */

CMR_EXPORT
CMR_ERROR CMRcamionStatsPrint(
  FILE* stream,                 /**< File stream to print to. */
  CMR_CAMION_STATISTICS* stats, /**< Pointer to statistics. */
  const char* prefix            /**< Prefix string to prepend to each printed line (may be \c NULL). */
);


/**
 * \brief Tests a matrix \f$ M \f$ for being a [Camion-signed](\ref camion).
 */

CMR_EXPORT
CMR_ERROR CMRcamionTestSigns(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Matrix \f$ M \f$. */
  bool* pisCamionSigned,        /**< Pointer for storing whether \f$ M \f$ is [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix,      /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
  CMR_CAMION_STATISTICS* stats, /**< Statistics for the computation (may be \c NULL). */
  double timeLimit              /**< Time limit to impose. */
);

/**
 * \brief Computes a [Camion-signed](\ref camion) version of a given ternary matrix \f$ M \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcamionComputeSigns(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Matrix \f$ M \f$ to be modified. */
  bool* pwasCamionSigned,       /**< Pointer for storing whether \f$ M \f$ was already [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix,      /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
  CMR_CAMION_STATISTICS* stats, /**< Statistics for the computation (may be \c NULL). */
  double timeLimit              /**< Time limit to impose. */
);


/**
 * \brief Orients the edges of the graph \p cograph such that the matrix \p matrix \f$ M \f$ is the corresponding
 *        network matrix, which implicitly tests if \p matrix is [Camion-signed](\ref camion).
 *
 * The cograph \f$ G = (V,E) \f$ has a spanning tree \f$ T \subseteq E \f$ indexed by the columns of \f$ M \f$.
 * Its complement \f$ E \setminus T \f$ is indexed by the rows of \f$ M \f$. The function assumes that
 * \f$ supp(M) = M(G,T)^\textsf{T} \f$ holds and attempts to compute an orientation \f$ A \f$ of the edges \f$ E \f$
 * (which is stored in \p arcsReversed) that corresponds to the signs of \f$ M \f$. \p *pisCamionSigned indicates
 * success. If successful, \f$ M \f$ is the network matrix of the digraph \f$ D = (V,A) \f$. Otherwise,
 * \p *psubmatrix will indicate a violating submatrix (if not \c NULL).
 */

CMR_EXPORT
CMR_ERROR CMRcamionCographicOrient(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix \f$ M \f$. */
  CMR_GRAPH* cograph,             /**< Cograph \f$ G = (V,E) \f$ claimed to correspond to \f$ M \f$. */
  CMR_GRAPH_EDGE* forestEdges,    /**< \f$ T \f$, ordered by the columns of \f$ M \f$. */
  CMR_GRAPH_EDGE* coforestEdges,  /**< \f$ E \setminus T \f$, ordered by the rows of \f$ M \f$. */
  bool* arcsReversed,             /**< Indicates, for each edge \f$ \{u, v\} \in E\f$, whether \f$ (u,v) \in A \f$
                                   **  (if \c false) or \f$ (v,u) \in A \f$  (if \c true). */
  bool* pisCamionSigned,          /**< Pointer for storing whether \f$ M \f$ is [Camion-signed](\ref camion). */
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a non-Camion submatrix (may be \c NULL). */
  CMR_CAMION_STATISTICS* stats    /**< Statistics for the computation (may be \c NULL). */
);


#ifdef __cplusplus
}
#endif

#endif /* CMR_CAMION_H */
