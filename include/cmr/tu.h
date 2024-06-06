#ifndef CMR_TU_H
#define CMR_TU_H

/**
 * \file tu.h
 *
 * \author Matthias Walter and Klaus Truemper
 *
 * \brief Recognition of [totally unimodular matrices](\ref tu).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/regular.h>
#include <cmr/matrix.h>
#include <cmr/camion.h>

typedef enum
{
  CMR_TU_ALGORITHM_DECOMPOSITION = 0, /**< \brief Algorithm based on Seymour's decomposition of regular matroids. */
  CMR_TU_ALGORITHM_EULERIAN = 1,      /**< \brief Enumeration algorithm based on Eulerian submatrices. */
  CMR_TU_ALGORITHM_PARTITION = 2      /**< \brief Enumeration algorithm based on criterion of Ghouila-Houri. */
} CMR_TU_ALGORITHM;

typedef struct
{
  CMR_TU_ALGORITHM algorithm; /**< \brief Algorithm to use. */
  CMR_SEYMOUR_PARAMS seymour; /**< \brief Parameters for testing via Seymour decomposition. */
  bool ternary;               /**< \brief Whether to create a ternary Seymour decomposition tree (default: \c true). */
  bool camionFirst;           /**< \brief If \c ternary is \c false, then then whether to run the Camion test first. */
} CMR_TU_PARAMS;

/**
 * \brief Initializes the default parameters for recognition of [totally unimodular](\ref tu) matrices.
 *
 * These are selected for minimum running time.
 */

CMR_EXPORT
CMR_ERROR CMRtuParamsInit(
  CMR_TU_PARAMS* params  /**< Pointer to parameters. */
);

/**
 * \brief Statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

typedef struct
{
  CMR_SEYMOUR_STATS seymour;          /**< Statistics for Seymour decomposition computation. */
  CMR_CAMION_STATISTICS camion;       /**< Statistics for Camion signing. */

  uint32_t enumerationRowSubsets;     /**< Number of considered row subsets in enumeration algorithm. */
  uint32_t enumerationColumnSubsets;  /**< Number of considered column subsets in enumeration algorithm. */
  double enumerationTime;             /**< Total time of enumeration algorithm. */

  uint32_t partitionRowSubsets;       /**< Number of considered row subsets in partition algorithm. */
  uint32_t partitionColumnSubsets;    /**< Number of considered column subsets in partition algorithm. */
  double partitionTime;               /**< Total time of partition algorithm. */
} CMR_TU_STATS;

/**
 * \brief Initializes all statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRtuStatsInit(
  CMR_TU_STATS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRtuStatsPrint(
  FILE* stream,             /**< File stream to print to. */
  CMR_TU_STATS* stats, /**< Pointer to statistics. */
  const char* prefix        /**< Prefix string to prepend to each printed line (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [totally unimodular](\ref tu).
 *
 * Tests if matrix \f$ M \f$ is totally unimodular and sets \p *pisTotallyUnimodular accordingly.
 *
 * If \f$ M \f$ is totally unimodular and \p pdec != \c NULL, then \p *pdec will contain a decomposition tree of the
 * regular matroid. The caller must release it via \ref CMRdecFree().
 *
 * If \f$ M \f$ is not totally unimodular and \p psubmatrix != \c NULL, then \p *psubmatrix will indicate a submatrix
 * of \f$ M \f$ with determinant \f$ -2 \f$ or \f$ 2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRtuTest(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  CMR_SEYMOUR_NODE** proot,   /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a submatrix with non-ternary determinant (may be \c NULL). */
  CMR_TU_PARAMS* params,      /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_TU_STATS* stats,        /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

/**
 * \brief Completes a subtree of an existing decomposition tree.
 *
 * Replace the node's subtree by a new one even if it exists. Note that different parameters may yield a different
 * subtree.
 *
 * \note Requires \p params.algorithm to be \ref CMR_TU_ALGORITHM_DECOMPOSITION.
 */

CMR_EXPORT
CMR_ERROR CMRtuCompleteDecomposition(
  CMR* cmr,               /**< \ref CMR environment. */
    CMR_SEYMOUR_NODE* dec,   /**< Pointer to the decomposition node that is the root of the new subtree. */
  CMR_TU_PARAMS* params,  /**< Parameters for the computation. */
  CMR_TU_STATS* stats,    /**< Statistics for the computation (may be \c NULL). */
  double timeLimit        /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_TU_H */
