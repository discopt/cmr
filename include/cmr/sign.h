#ifndef TU_SIGN_H
#define TU_SIGN_H

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
CMR_ERROR TUtestSignDbl(
  TU* tu,                 /**< \ref TU environment. */
  TU_DBLMAT* matrix,      /**< Sparse double matrix. */
  bool* pcorrectSign,     /**< Pointer for storing whether \p matrix can be signed. */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

/**
 * \brief Modifies signs of double matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR TUcorrectSignDbl(
  TU* tu,                 /**< \ref TU environment */
  TU_DBLMAT* matrix,      /**< Sparse double matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was modified. */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

/**
 * \brief Tests if signs of int matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR TUtestSignInt(
  TU* tu,                 /**< \ref TU environment */
  TU_INTMAT* matrix,      /**< Sparse int matrix */
  bool* pcorrectSign,     /**< Pointer for storing whether \p matrix was modified. */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

/**
 * \brief Modifies signs of int matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR TUcorrectSignInt(
  TU* tu,                 /**< \ref TU environment */
  TU_INTMAT* matrix,      /**< Sparse int matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was modified. */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);


/**
 * \brief Tests if signs of char matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR TUtestSignChr(
  TU* tu,                 /**< \ref TU environment */
  TU_CHRMAT* matrix,      /**< Sparse char matrix */
  bool* pcorrectSign,     /**< Pointer for storing whether \p matrix has the right sign. */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);


/**
 * \brief Modifies signs of char matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

CMR_EXPORT
CMR_ERROR TUcorrectSignChr(
  TU* tu,                 /**< \ref TU environment */
  TU_CHRMAT* matrix,      /**< Sparse char matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was modified. */
  TU_SUBMAT** psubmatrix  /**< Pointer for storing a submatrix with bad determinant (may be \c NULL). */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SIGN_H */
