#ifndef CMR_LINEAR_ALGEBRA_H
#define CMR_LINEAR_ALGEBRA_H

/**
 * \file linear_algebra.h
 *
 * \author Matthias Walter
 */

#include <cmr/env.h>
#include <cmr/matrix.h>

#include <stdio.h>
#include <stdint.h.>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Computes the determinant of an 8-bit integer matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatDeterminant(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  int64_t* pdeterminant   /**< Pointer for storing the determinant. */
);

/**
 * \brief Computes the determinant of an int matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatDeterminant(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,     /**< Matrix. */
  int64_t* pdeterminant   /**< Pointer for storing the determinant. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_LINEAR_ALGEBRA_H */
