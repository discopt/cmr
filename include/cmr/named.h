#ifndef CMR_NAMED_H
#define CMR_NAMED_H

/**
 * \file special.h
 *
 * \author Matthias Walter
 *
 * \brief Functionality for special matrix.
 */

#include <cmr/env.h>
#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Checks if the given \p matrix is an identity matrix.
 *
 * If \p matrix is not identity matrix then \p *porder indicates the order; Otherwise, it is set to \c SIZE_MAX.
 */

CMR_EXPORT
CMR_ERROR CMRisIdentityMatrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< Matrix. */
  size_t* porder      /**< Pointer for storing the order of the matrix. */
);

/**
 * \brief Constructs an identity matrix of given \p order.
 */

CMR_EXPORT
CMR_ERROR CMRcreateIdentityMatrix(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t order,         /**< Order of the matrix. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Checks if the given \p matrix represents the matroid \f$ R_{10} \f$.
 */

CMR_EXPORT
CMR_ERROR CMRisR10Matrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< Matrix. */
  bool* pisR10        /**< Pointer for storing whether \p matrix represents \f$ R_{10} \f$. */
);

/**
 * \brief Constructs a representation matrix for \f$ R_{10} \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcreateR10Matrix(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t index,         /**< Which of the two matrices to return; among {1,2}. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Constructs a representation matrix for \f$ R_{12} \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcreateR12Matrix(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t index,         /**< Which of the matrices to return; must be 1. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Constructs a representation matrix for \f$ M(K_5) \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcreateK5Matrix(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t index,         /**< Which of the matrices to return; must be 1. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Constructs a representation matrix for \f$ M(K_{3,3}) \f$.
 */

CMR_EXPORT
CMR_ERROR CMRcreateK33Matrix(
  CMR* cmr,             /**< \ref CMR environment. */
  size_t index,         /**< Which of the matrices to return; must be 1. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_NAMED_H */

