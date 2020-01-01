#ifndef TU_SIGN_H
#define TU_SIGN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/matrix.h>

/**
 * \brief Tests if signs of char matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 * 
 * Returns \c true if and only if the signs are correct.
 */

TU_EXPORT
bool TUtestSignChar(
  TU* tu, /**< TU environment */
  TU_SPARSE_CHAR* matrix, /**< Sparse char matrix */
  TU_SUBMATRIX** submatrix /**< If not \c NULL, a submatrix with bad determinant is stored. */
);


/**
 * \brief Modifies signs of char matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 *
 * Returns \c true if and only if the signs were already correct.
 */

TU_EXPORT
bool TUcorrectSignChar(
  TU* tu, /**< TU environment */
  TU_SPARSE_CHAR* matrix, /**< Sparse char matrix */
  TU_SUBMATRIX** submatrix /**< If not \c NULL, a submatrix with bad determinant is stored. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SIGN_H */
