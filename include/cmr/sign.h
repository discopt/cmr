#ifndef CMR_SIGN_H
#define CMR_SIGN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cmr/matrix.h>

/**
 * \brief Tests if signs of double matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR CMRtestSignDbl(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,      /**< Sparse double matrix. */
  bool* pcorrectSign,     /**< Pointer for storing whether \p matrix can be signed. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

/**
 * \brief Modifies signs of double matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR CMRcorrectSignDbl(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_DBLMAT* matrix,      /**< Sparse double matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was modified. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

/**
 * \brief Tests if signs of int matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR CMRtestSignInt(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_INTMAT* matrix,      /**< Sparse int matrix */
  bool* pcorrectSign,     /**< Pointer for storing whether \p matrix was modified. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

/**
 * \brief Modifies signs of int matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR CMRcorrectSignInt(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_INTMAT* matrix,      /**< Sparse int matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was modified. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);


/**
 * \brief Tests if signs of char matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR CMRtestSignChr(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_CHRMAT* matrix,      /**< Sparse char matrix */
  bool* pcorrectSign,     /**< Pointer for storing whether \p matrix has the right sign. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);


/**
 * \brief Modifies signs of char matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR CMRcorrectSignChr(
  CMR* cmr,                 /**< \ref CMR environment */
  CMR_CHRMAT* matrix,      /**< Sparse char matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was modified. */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SIGN_H */
