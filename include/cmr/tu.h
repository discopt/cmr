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
  CMR_TU_ALGORITHM algorithm;     /**< \brief Algorithm to use. */
  CMR_REGULAR_PARAMETERS regular; /**< \brief Parameters for regularity test. */
} CMR_TU_PARAMETERS;

/**
 * \brief Initializes the default parameters for recognition of [totally unimodular](\ref tu) matrices.
 *
 * These are selected for minimum running time.
 */

CMR_EXPORT
CMR_ERROR CMRparamsTotalUnimodularityInit(
  CMR_TU_PARAMETERS* params  /**< Pointer to parameters. */
);

/**
 * \brief Statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

typedef struct
{
  uint32_t totalCount;            /**< Total number of invocations. */
  double totalTime;               /**< Total time of all invocations. */
  CMR_CAMION_STATISTICS camion;   /**< Camion signing. */
  CMR_REGULAR_STATISTICS regular; /**< Regularity test. */
} CMR_TU_STATISTICS;

/**
 * \brief Initializes all statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRstatsTotalUnimodularityInit(
  CMR_TU_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRstatsTotalUnimodularityPrint(
  FILE* stream,             /**< File stream to print to. */
  CMR_TU_STATISTICS* stats, /**< Pointer to statistics. */
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
CMR_ERROR CMRtestTotalUnimodularity(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  CMR_DEC** pdec,             /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a submatrix with non-ternary determinant (may be \c NULL). */
  CMR_TU_PARAMETERS* params,  /**< Parameters for the computation (may be \c NULL for defaults). */
  CMR_TU_STATISTICS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);


#ifdef __cplusplus
}
#endif

#endif /* CMR_TU_H */
