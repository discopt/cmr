#ifndef CMR_TU_H
#define CMR_TU_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/env.h>
#include <cmr/regular.h>
#include <cmr/matrix.h>

/**
 * \brief Tests a double matrix for total unimodularity.
 *
 * Returns \c true if and only if \p matrix is TU, where we consider entries to be equal to -1, 0 or +1 if their
 * distance is at most \p epsilon.
 *
 * If \p pdecomposition is not \c NULL and the algorithm has to test regularity of the support matrix, then
 * \c *pdecomposition will point to a decomposition tree for which the caller must use \ref CMRdecFree to free memory.
 * It is set to \c NULL in case regularity of the support matrix does not need to be determined.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute determinant larger than
 * 1 will be searched, which may cause extra computational effort. In this case, *\p psubmatrix will point to this
 * submatrix for which the caller must use \ref CMRsubmatFree to free memory. It is set to \c NULL otherwise.
 */

CMR_EXPORT
CMR_ERROR CMRtestTotalUnimodularityDbl(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_DBLMAT* matrix,      /**< Double matrix */
  double epsilon,         /**< Absolute error tolerance */
  bool* pisTU,            /**< Pointer for storing whether \p matrix is TU.*/
  CMR_TU_DEC** pdec,          /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
);

/**
 * \brief Tests an int matrix for total unimodularity.
 *
 * Returns \c true if and only if \p matrix is TU.
 *
 * If \p pdecomposition is not \c NULL and the algorithm has to test regularity of the support matrix, then
 * \c *pdecomposition will point to a decomposition tree for which the caller must use \ref CMRdecFree to free memory.
 * It is set to \c NULL in case regularity of the support matrix does not need to be determined.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute determinant larger than
 * 1 will be searched, which may cause extra computational effort. In this case, *\p psubmatrix will point to this
 * submatrix for which the caller must use \ref CMRsubmatFree to free memory. It is set to \c NULL otherwise.
 */

CMR_EXPORT
CMR_ERROR CMRtestTotalUnimodularityInt(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_INTMAT* matrix,      /**< Int matrix */
  bool* pisTU,            /**< Pointer for storing whether \p matrix is TU.*/
  CMR_TU_DEC** pdec,          /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
);

/**
 * \brief Tests a char matrix for total unimodularity.
 *
 * Returns \c true if and only if \p matrix is TU.
 *
 * If \p pdecomposition is not \c NULL and the algorithm has to test regularity of the support matrix, then
 * \c *pdecomposition will point to a decomposition tree for which the caller must use \ref CMRdecFree to free memory.
 * It is set to \c NULL in case regularity of the support matrix does not need to be determined.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not TU, then a submatrix with an absolute determinant larger than
 * 1 will be searched, which may cause extra computational effort. In this case, *\p psubmatrix will point to this
 * submatrix for which the caller must use \ref CMRsubmatFree to free memory. It is set to \c NULL otherwise.
 */

CMR_EXPORT
CMR_ERROR CMRtestTotalUnimodularityChr(
  CMR* cmr,                /**< \ref CMR environment */
  CMR_CHRMAT* matrix,      /**< Char matrix to be tested. */
  bool* pisTU,            /**< Pointer for storing whether \p matrix is TU.*/
  CMR_TU_DEC** pdec,          /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a bad submatrix with a bad determinant (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_TU_H */
