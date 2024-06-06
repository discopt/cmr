#ifndef CMR_REGULAR_H
#define CMR_REGULAR_H

/**
 * \file regular.h
 *
 * \author Matthias Walter and Klaus Truemper
 *
 * \brief Recognition of [regular matrices](\ref regular).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/matroid.h>
#include <cmr/seymour.h>
#include <cmr/graph.h>
#include <cmr/series_parallel.h>
#include <cmr/graphic.h>
#include <cmr/network.h>

typedef struct
{
  CMR_SEYMOUR_PARAMS seymour;
} CMR_REGULAR_PARAMS;

/**
 * \brief Initializes the default parameters for regularity testing.
 *
 * These are selected for minimum running time.
 */

CMR_EXPORT
CMR_ERROR CMRregularParamsInit(
  CMR_REGULAR_PARAMS* params  /**< Pointer to parameters. */
);

/**
 * \brief Statistics for regular matroid recognition algorithm.
 */

typedef struct
{
  CMR_SEYMOUR_STATS seymour;  /**< Statistics for Seymour decomposition computations. */
} CMR_REGULAR_STATS;


/**
 * \brief Initializes all statistics for regularity test computations.
 */

CMR_EXPORT
CMR_ERROR CMRregularStatsInit(
  CMR_REGULAR_STATS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for regularity test computations.
 */

CMR_EXPORT
CMR_ERROR CMRregularStatsPrint(
  FILE* stream,             /**< File stream to print to. */
  CMR_REGULAR_STATS* stats, /**< Pointer to statistics. */
  const char* prefix        /**< Prefix string to prepend to each printed line (may be \c NULL). */
);

/**
 * \brief Tests binary linear matroid for regularity.
 *
 * If \p pdec is not \c NULL, \c *pdec will be a (partial) decomposition tree.
 * If \p completeTree is \c true, then the decomposition tree is complete. Otherwise, it must only contain sufficient
 * information in order to determine regularity.
 *
 * If \p pminor is not \c NULL and \p matrix is not regular, then an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor is searched.
 * This causes additional computational effort!
 */

CMR_EXPORT
CMR_ERROR CMRregularTest(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Input matrix. */
  bool *pisRegular,           /**< Pointer for storing whether \p matrix is regular. */
  CMR_SEYMOUR_NODE** proot,   /**< Pointer for storing the Seymour decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,         /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMS* params, /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

/**
 * \brief Completes a subtree of an existing decomposition tree.
 *
 * Replace the node's subtree by a new one even if it exists. Note that different parameters may yield a different
 * subtree.
 */

CMR_EXPORT
CMR_ERROR CMRregularCompleteDecomposition(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* dec,      /**< Pointer to the decomposition node that is the root of the new subtree. */
  CMR_REGULAR_PARAMS* params, /**< Parameters for the computation (may be \c NULL). */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

/**
 * \brief Refines a list of decomposition nodes.
 *
 * Replace the nodes' subtrees by new ones even if they exist.
 */

CMR_EXPORT
CMR_ERROR CMRregularRefineDecomposition(
  CMR* cmr,                   /**< \ref CMR environment. */
  size_t numNodes,            /**< Number of nodes to refine. */
  CMR_SEYMOUR_NODE** nodes,   /**< Array of decomposition nodes to refine. */
  CMR_REGULAR_PARAMS* params, /**< Parameters for the computation (may be \c NULL). */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_REGULAR_H */

