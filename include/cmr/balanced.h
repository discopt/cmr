#ifndef CMR_BALANCED_H
#define CMR_BALANCED_H

/**
 * \file balanced.h
 *
 * \author Henk Kraaij and Matthias Walter
 *
 * \brief Recognition of [balanced matrices](\ref balanced).
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/series_parallel.h>

#include <inttypes.h>

typedef enum
{
  CMR_BALANCED_ALGORITHM_AUTO = 0,      /**< \brief Automatically select a fast algorithm. */
  CMR_BALANCED_ALGORITHM_SUBMATRIX = 1, /**< \brief Exponential-time enumeration algorithm based on submatrices. */
  CMR_BALANCED_ALGORITHM_GRAPH = 2      /**< \brief Polynomial-time algorithm based on graphs. */
} CMR_BALANCED_ALGORITHM;

typedef struct
{
  CMR_BALANCED_ALGORITHM algorithm; /**< \brief Algorithm to use. */
  bool seriesParallel;              /**< \brief Whether to carry out series-parallel operations as preprocessing. */
} CMR_BALANCED_PARAMS;

/**
 * \brief Initializes the default parameters for recognition of [balanced](\ref balanced) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRbalancedParamsInit(
  CMR_BALANCED_PARAMS* params /**< Pointer to parameters. */
);

/**
 * \brief Statistics for recognition algorithm for [balanced](\ref balanced) matrices.
 */

typedef struct
{
  uint32_t totalCount;              /**< Total number of invocations. */
  double totalTime;                 /**< Total time of all invocations. */
  CMR_SP_STATISTICS seriesParallel; /**< Statistics for series-parallel algorithm. */
  size_t enumeratedRowSubsets;      /**< Number of enumerated row subsets. */
  size_t enumeratedColumnSubsets;   /**< Number of enumerated column subsets. */
} CMR_BALANCED_STATS;

/**
 * \brief Initializes all statistics for recognition algorithm for [balanced](\ref balanced) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRbalancedStatsInit(
  CMR_BALANCED_STATS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for recognition algorithm for [balanced](\ref balanced) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRbalancedStatsPrint(
  FILE* stream,               /**< File stream to print to. */
  CMR_BALANCED_STATS* stats,  /**< Pointer to statistics. */
  const char* prefix          /**< Prefix string to prepend to each printed line (may be \c NULL). */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [balanced](\ref balanced).
 *
 * Tests if matrix \f$ M \f$ is balanced and sets \p *pisBalanced accordingly.
 * Automatically decides which algorithm to use.
 *
 * If \f$ M \f$ is not balanced and \p psubmatrix != \c NULL, then \p *psubmatrix will indicate a submatrix
 * of \f$ M \f$ with exactly two nonzeros in each row and in each column and with determinant \f$ -2 \f$ or \f$ 2 \f$.
 */CMR_SP_STATISTICS seriesParallel;     /**< Statistics for series-parallel algorithm. */

CMR_EXPORT
CMR_ERROR CMRbalancedTest(
  CMR* cmr,                     /**< \ref CMR environment */
  CMR_CHRMAT* matrix,           /**< Matrix \f$ M \f$. */
  bool* pisBalanced,            /**< Pointer for storing whether \f$ M \f$ is balanced. */
  CMR_SUBMAT** psubmatrix,      /**< Pointer for storing a minimal nonbalanced submatrix (may be \c NULL). */
  CMR_BALANCED_PARAMS* params,  /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_BALANCED_STATS* stats,    /**< Statistics for the computation (may be \c NULL). */
  double timeLimit              /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_BALANCED_H */
