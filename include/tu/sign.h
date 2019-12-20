#ifndef TU_SIGN_H
#define TU_SIGN_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/matrix.h>

/**
 * \brief Tests if signs of matrix nonzeros qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 * 
 * The function returns \c true if and only if the signs are correct.
 */

TU_EXPORT
bool TUtestSignChar(
  TU* tu, /**< TU environment */
  TU_SPARSE_CHAR* matrix /**< Sparse matrix */
);


/**
 * \brief Modifies signs of matrix nonzeros to qualify for being TU.
 *
 * The \p matrix is assumed to be ternary.
 *
 * The function returns \c true if and only if the signs were already correct.
 */

TU_EXPORT
bool TUcorrectSignChar(
  TU* tu, /**< TU environment */
  TU_SPARSE_CHAR* matrix /**< Sparse matrix */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_SIGN_H */
