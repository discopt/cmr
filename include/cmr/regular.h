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
#include <cmr/graph.h>
#include <cmr/series_parallel.h>
#include <cmr/graphic.h>
#include <cmr/network.h>

typedef enum
{
  CMR_DEC_CONSTRUCT_NONE = 0,
  CMR_DEC_CONSTRUCT_LEAVES = 1,
  CMR_DEC_CONSTRUCT_ALL = 2,
} CMR_DEC_CONSTRUCT;

/**
 * \brief Flags that control termination of the regular matroid decomposition algorithm.
 *
 * \see The desired types can be set by modifying \ref CMR_REGULAR_PARAMS.treeFlags.
 */

typedef enum
{
  CMR_REGULAR_TREE_FLAGS_RECURSE = 1,
    /**< Do not process child nodes. */
  CMR_REGULAR_TREE_FLAGS_STOP_IRREGULAR = 2,
    /**< Stop if a (grand-)child is irregular. */
  CMR_REGULAR_TREE_FLAGS_STOP_NONGRAPHIC = 4,
    /**< Stop if a (grand-)child is non-graphic. */
  CMR_REGULAR_TREE_FLAGS_STOP_NONCOGRAPHIC = 8,
    /**< Stop if a (grand-)child is non-cographic. */
  CMR_REGULAR_TREE_FLAGS_STOP_NONGRAPHIC_NONCOGRAPHIC = 16,
    /**< Stop if a (grand-)child is non-graphic **and** non-cographic. */

  CMR_REGULAR_TREE_FLAGS_DEFAULT = CMR_REGULAR_TREE_FLAGS_RECURSE | CMR_REGULAR_TREE_FLAGS_STOP_IRREGULAR

} CMR_REGULAR_TREE_FLAGS;

typedef struct
{
  bool directGraphicness;
  /**< \brief Whether to use fast graphicness routines; default: \c true */
  bool seriesParallel;
  /**< \brief Whether to allow series-parallel operations in the decomposition tree; default: \c true */
  bool planarityCheck;
  /**< \brief Whether minors identified as graphic should still be checked for cographicness; default: \c false. */

  int treeFlags;
  /**< \brief Flags controlling the decomposition algorithm. See \ref  */

  /**< \brief Whether to compute a complete decomposition tree (even if already non-regular; default: \c false. */
  bool threeSumPivotChildren;
  /**< \brief Whether pivots for 3-sums shall be applied such that the matrix contains both child matrices as
   **         submatrices, if possible. */
  int threeSumStrategy;
  /**< \brief Whether to perform pivots to change the rank distribution, and how to construct the children.
   **
   ** The value is a bit-wise or of three decisions. The first decision is that of the **rank distribution**:
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_NO_PIVOTS to not change the rank distribution (default), or
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS to enforce distributed ranks (1 + 1), or
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_CONCENTRATED_RANK to enforce concentrated ranks (2 + 0).
   **
   **  The second decision determines the layout of the **first child** matrix:
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_FIRST_WIDE for a wide first child (default) in case of distributed ranks,
   **     or
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_FIRST_TALL for a tall first child in that case.
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_FIRST_MIXED for a mixed first child (default) in case of concentrated
   **     ranks, or
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_FIRST_ALLREPR for a first child with all representing rows in that case.
   **
   **  Similarly, the third decision determines the layout of the **second child** matrix:
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_SECOND_WIDE for a wide second child (default) in case of distributed ranks,
   **     or
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_SECOND_TALL for a tall second child in that case.
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_SECOND_MIXED for a mixed second child (default) in case of concentrated
   **     ranks, or
   **   - \ref CMR_MATROID_DEC_THREESUM_FLAG_SECOND_ALLREPR for a first second with all representing rows in that case.
   **
   ** \see \ref matroid_decomposition for a description of these layouts.
   **
   ** A decomposition as described by Seymour can be selected via \ref CMR_MATROID_DEC_THREESUM_FLAG_SEYMOUR.
   ** A decomposition as used by Truemper can be selected via \ref CMR_MATROID_DEC_THREESUM_FLAG_TRUEMPER.
   ** The default is to not carry out any pivots and choose Seymour's or Truemper's definition depending on the rank
   ** distribution. */

  CMR_DEC_CONSTRUCT graphs;
  /**< \brief Which (co)graphs to construct; default: \ref CMR_DEC_CONSTRUCT_NONE. */
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
  uint32_t totalCount;                  /**< Total number of invocations. */
  double totalTime;                     /**< Total time of all invocations. */
  CMR_SP_STATISTICS seriesParallel;     /**< Statistics for series-parallel algorithm. */
  CMR_GRAPHIC_STATISTICS graphic;       /**< Statistics for direct (co)graphic checks. */
  CMR_NETWORK_STATISTICS network;       /**< Statistics for direct (co)network checks. */
  CMR_CAMION_STATISTICS camion;         /**< Statistics for Camion signing. */
  uint32_t sequenceExtensionCount;      /**< Number of extensions of sequences of nested minors. */
  double sequenceExtensionTime;         /**< Time of extensions of sequences of nested minors. */
  uint32_t sequenceGraphicCount;        /**< Number (co)graphicness tests applied to sequence of nested minors. */
  double sequenceGraphicTime;           /**< Time of (co)graphicness tests applied to sequence of nested minors. */
  uint32_t enumerationCount;            /**< Number of calls to enumeration algorithm for candidate 3-separations. */
  double enumerationTime;               /**< Time of enumeration of candidate 3-separations. */
  uint32_t enumerationCandidatesCount;  /**< Number of enumerated candidates for 3-separations. */
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
  CMR_MATROID_DEC** pdec,     /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,         /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMS* params, /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

/**
 * \brief Completes a subtree of an existing decomposition tree.
 *
 * Replace the node's subtree by a new one even if it exists. Note that different parameters may yield a different
 * subtree.s
 */

CMR_EXPORT
CMR_ERROR CMRregularCompleteDecomposition(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,       /**< Pointer to the decomposition node that is the root of the new subtree. */
  CMR_REGULAR_PARAMS* params, /**< Parameters for the computation. */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_REGULAR_H */

