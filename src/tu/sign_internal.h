#ifndef TU_SIGN_INTERNAL_H
#define TU_SIGN_INTERNAL_H

#include <tu/env.h>
#include <tu/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Ensures that signs of sequentially connected char matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary. If sign changes are necessary, only \p matrix is
 * modified. In particular, \p transpose is remain unchanged. Expects \p matrix to have the
 * additional rowStart entry (\see TODO).
 *
 * The function returns 0 if the signs were already correct. Otherwise, it returns 'm' or 't'
 * depending on whether \p matrix or \p transpose was modified (or would have been modified).
 *
 * If \p submatrix is not \c NULL and sign changes are necessary, then a submatrix with determinant
 * -2 or +2 is stored in \c *submatrix and the caller must use \ref TUfreeSubmatrix to free its
 * memory. It is set to \c NULL if no sign changes are needed.
 */

char signSequentiallyConnected(
  TU* tu,                     /**< TU environment. */
  TU_SPARSE_CHAR* matrix,     /**< Sparse matrix. */
  TU_SPARSE_CHAR* transpose,  /**< Transpose of \p matrix. */
  bool change,                /**< Whether signs of \p matrix should be changed if necessary */
  TU_SUBMATRIX** submatrix    /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SIGN_INTERNAL_H */


