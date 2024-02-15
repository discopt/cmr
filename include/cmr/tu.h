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
  CMR_TU_ALGORITHM_SUBMATRIX = 1,     /**< \brief Enumeration algorithm based on submatrices. */
  CMR_TU_ALGORITHM_PARTITION = 2      /**< \brief Enumeration algorithm based on criterion of Ghouila-Houri. */
} CMR_TU_ALGORITHM;

typedef struct
{
  CMR_TU_ALGORITHM algorithm; /**< \brief Algorithm to use. */
  bool directCamion;          /**< \brief Whether to directly test signing of matrix (default: \c false). */
  CMR_REGULAR_PARAMS regular; /**< \brief Parameters for regularity test. */
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

typedef CMR_REGULAR_STATS CMR_TU_STATS;

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
  CMR_MATROID_DEC** pdec,     /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a submatrix with non-ternary determinant (may be \c NULL). */
  CMR_TU_PARAMS* params,      /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_TU_STATS* stats,        /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_TU_H */
