#ifndef TU_SIGN_H
#define TU_SIGN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/matrix.h>

/**
 * \brief Tests if signs of double matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUtestSignDbl(
  TU* tu,                 /**< \ref TU environment. */
  TU_DBLMAT* matrix,      /**< Sparse double matrix. */
  bool* pisSignable,      /**< Pointer for storing whether \p matrix can be signed. */
  TU_SUBMAT** psubmatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

/**
 * \brief Modifies signs of double matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUcorrectSignDbl(
  TU* tu,               /**< \ref TU environment */
  TU_DBLMAT* matrix,    /**< Sparse double matrix */
  bool* pisSignable,    /**< Pointer for storing whether \p matrix was modified. */
  TU_SUBMAT** submatrix /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

/**
 * \brief Tests if signs of int matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUtestSignInt(
  TU* tu,                 /**< \ref TU environment */
  TU_INTMAT* matrix,      /**< Sparse int matrix */
  bool* pisSignable,      /**< Pointer for storing whether \p matrix can be signed. */
  TU_SUBMAT** psubmatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

/**
 * \brief Modifies signs of double matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUcorrectSignInt(
  TU* tu,                 /**< \ref TU environment */
  TU_INTMAT* matrix,      /**< Sparse int matrix */
  bool* pisSignable,      /**< Pointer for storing whether \p matrix can be signed. */
  TU_SUBMAT** psubmatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);


/**
 * \brief Tests if signs of char matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUtestSignChr(
  TU* tu,                 /**< \ref TU environment */
  TU_CHRMAT* matrix,      /**< Sparse char matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix has the right sign. */
  TU_SUBMAT** psubmatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);


/**
 * \brief Modifies signs of char matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 */

TU_EXPORT
TU_ERROR TUcorrectSignChr(
  TU* tu,                 /**< \ref TU environment */
  TU_CHRMAT* matrix,      /**< Sparse char matrix */
  bool* palreadySigned,   /**< Pointer for storing whether \p matrix was already signed correctly. */
  TU_SUBMAT** psubmatrix  /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SIGN_H */
