#ifndef CMR_CTU_H
#define CMR_CTU_H

/**
 * \file ctu.h
 *
 * \author Matthias Walter and Klaus Truemper
 *
 * \brief Recognition of [complement totally unimodular matrices](\ref ctu).
 */

#ifdef __cplusplus
extern "C" {
#endif
  
#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/tu.h>

/**
 * \brief Statistics for recognition algorithm for [totally unimodular](\ref tu) matrices.
 */

typedef struct
{
  uint32_t totalCount;    /**< Total number of invocations. */
  double totalTime;     /**< Total time of all invocations. */
  CMR_TU_STATISTICS tu; /**< Total unimodularity test. */
} CMR_CTU_STATISTICS;

/**
 * \brief Initializes all statistics for recognition algorithm for [complement totally unimodular](\ref ctu) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRstatsComplementTotalUnimodularityInit(
  CMR_CTU_STATISTICS* stats /**< Pointer to statistics. */
);

/**
 * \brief Prints statistics for recognition algorithm for [complement totally unimodular](\ref ctu) matrices.
 */

CMR_EXPORT
CMR_ERROR CMRstatsComplementTotalUnimodularityPrint(
  FILE* stream,             /**< File stream to print to. */
  CMR_CTU_STATISTICS* stats, /**< Pointer to statistics. */
  const char* prefix        /**< Prefix string to prepend to each printed line (may be \c NULL). */
);

/**
 * \brief Carries out a row- and column-complement operations on the binary matrix.
 */

CMR_EXPORT
CMR_ERROR CMRcomplementRowColumn(
  CMR* cmr,                             /**< \ref CMR environment */
  CMR_CHRMAT* matrix,                   /**< Input matrix. */
  size_t complementRow,                 /**< Row to be complemented (\c SIZE_MAX for no row complement). */
  size_t complementColumn,              /**< Column to be complemented (\c SIZE_MAX for no column complement). */
  CMR_CHRMAT** presult                  /**< Resulting matrix. */
);

/**
 * \brief Tests a matrix \f$ M \f$ for being [complement totally unimodular](\ref ctu).
 *
 * Tests if matrix \f$ M \f$ is complement totally unimodular and sets \p *pisComplementTotallyUnimodular accordingly.
 *
 * If \f$ M \f$ is not complement totally unimodular and \p pcomplementRow != \c NULL and
 * \p pcomplementColumn != \c NULL, then \p *pcomplementRow and \p *pcomplementColumn will indicate the row and column
 * that need to be complemented for obtaining a matrix that is not [totally unimodular](\ref tu).
 * If no row/column needs to be complemented, then the respective variables are set to \c SIZE_MAX.
 */

CMR_EXPORT
CMR_ERROR CMRtestComplementTotalUnimodularity(
  CMR* cmr,                             /**< \ref CMR environment */
  CMR_CHRMAT* matrix,                   /**< Matrix \f$ M \f$. */
  bool* pisComplementTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is complement totally unimodular. */
  size_t* pcomplementRow,               /**< Pointer for storing the row to be complemented (may be \c NULL). */
  size_t* pcomplementColumn,            /**< Pointer for storing the column to be complemented (may be \c NULL). */
  CMR_CTU_STATISTICS* stats,            /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                      /**< Time limit to impose. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_CTU_H */
