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

/**
 * \brief Tests a matrix \f$ M \f$ for being [totally unimodular](\ref tu).
 *
 * Tests if matrix \f$ M \f$ is totally unimodular and sets \p *pisTotallyUnimodular accordingly.
 *
 * If \f$ M \f$ is totally unimodular and \p pdec != \c NULL, then \p *pdec will contain a decomposition tree of the
 * regular matroid. The caller must release it via \ref CMRtudecFree().
 *
 * If \f$ M \f$ is not totally unimodular and \p psubmatrix != \c NULL, then \p *psubmatrix will indicate a submatrix
 * of \f$ M \f$ with determinant \f$ -2 \f$ or \f$ 2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRtestTotalUnimodularity(
  CMR* cmr,                   /**< \ref CMR environment */
  CMR_CHRMAT* matrix,         /**< Matrix \f$ M \f$. */
  bool* pisTotallyUnimodular, /**< Pointer for storing whether \f$ M \f$ is totally unimodular. */
  CMR_TU_DEC** pdec,          /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix     /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_TU_H */
